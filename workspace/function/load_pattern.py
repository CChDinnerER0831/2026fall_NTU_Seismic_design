import openseespy.opensees as ops
import numpy as np
from abc import ABC, abstractmethod
from function.structure_model import StructureModel


class BaseLoadPattern(ABC):
    """載重模式抽象基底類別：定義初始化與每步迭代更新的統一介面"""

    def __init__(self, struct: StructureModel, pattern_tag: int = 2, ts_tag: int = 2):
        self.struct = struct
        self.pattern_tag = pattern_tag
        self.ts_tag = ts_tag
        self.current_loads: dict[int, tuple[float, float, float]] = {}

    @abstractmethod
    def compute_pattern(self, step: int) -> dict[int, tuple[float, float, float]]:
        """計算第 step 步的節點載重向量 {nodeID: (Fx, Fy, Mz)}"""
        pass

    def apply(self, step: int, is_update: bool = False):
        """將計算出的載重套用至 OpenSees（若為迭代更新則先鎖定前一步狀態並重置 Pattern）"""
        new_loads = self.compute_pattern(step)

        if is_update:
            ops.loadConst('-time', 0.0)
            ops.remove('loadPattern', self.pattern_tag)
            ops.remove('timeSeries', self.ts_tag)

        ops.timeSeries('Linear', self.ts_tag)
        ops.pattern('Plain', self.pattern_tag, self.ts_tag)
        for nid, load_vec in new_loads.items():
            ops.load(nid, *load_vec)

        self.current_loads = new_loads


class FixedIPTLoadPattern(BaseLoadPattern):
    """固定模式：直接使用 .ipt 中 *LATERAL_PATTERN 的分佈，不隨步數改變"""

    def compute_pattern(self, step: int) -> dict[int, tuple[float, float, float]]:
        return self.struct.lateral_pattern


class AdaptiveModalLoadPattern(BaseLoadPattern):
    """
    自適應振型側推 (Adaptive Modal Pushover)：
    每隔 update_interval 步重新執行特徵值分析 (ops.eigen)，
    依據結構進入塑性後的「當前第一振態振型 Φ1」與「質量 m」動態更新側力分佈 (F_i ∝ m_i * Φ_1i)。
    """

    def __init__(self, struct: StructureModel, update_interval: int = 10,
                 pattern_tag: int = 2, ts_tag: int = 2):
        super().__init__(struct, pattern_tag, ts_tag)
        self.update_interval = update_interval

    def should_update(self, step: int) -> bool:
        return step > 0 and (step % self.update_interval == 0)

    def compute_pattern(self, step: int) -> dict[int, tuple[float, float, float]]:
        # 求解當前切線勁度下的第一振態特徵值
        ops.eigen('-genBandArpack', 1)

        weights = {}
        for nid, (mx, _, _) in self.struct.masses.items():
            if nid in self.struct.fixities:
                continue
            phi_1 = abs(ops.nodeEigenvector(nid, 1, 1))
            weights[nid] = mx * phi_1

        total_w = sum(weights.values())
        if total_w == 0:
            return self.struct.lateral_pattern

        return {nid: (w / total_w, 0.0, 0.0) for nid, w in weights.items()}