import openseespy.opensees as ops
import math
import os
from pathlib import Path
from function.material_property import MaterialProperty
from function.section_property import SectionProperty, WSection, BoxSection
import matplotlib.pyplot as plt


class StructureModel:
    """
    讀取 FRAME26 格式的 .ipt 輸入檔並建立 OpenSees 2D 非線性構架模型
    """

    def __init__(self, mat: MaterialProperty, sec_prop: SectionProperty):
        self.mat = mat
        self.sec_prop = sec_prop

        self.control_node: int = -1
        self.nodes: dict[int, tuple[float, float]] = {}
        self.fixities: dict[int, tuple[int, int, int]] = {}
        self.elements: list[dict] = []
        self.masses: dict[int, tuple[float, float, float]] = {}
        self.gravity_loads: dict[int, tuple[float, float, float]] = {}
        self.lateral_pattern: dict[int, tuple[float, float, float]] = {}

        # 提供給 PushoverAnalyzer 使用
        self.base_nodes: list[int] = []
        self.total_height: float = 0.0

    def read_ipt(self, filepath: str | Path):
        """解析 .ipt 輸入檔"""
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        current_section = None
        mm = self.mat.mm
        kN = self.mat.kN
        MPa = self.mat.MPa
        GPa = self.mat.GPa

        with open(filepath, 'r', encoding='utf-8') as f:
            for line_num, raw_line in enumerate(f, start=1):
                # 去除註解 (# 或 !) 與前後空白
                line = raw_line.split('#')[0].split('!')[0].strip()
                if not line:
                    continue

                # 切換區段標籤
                if line.startswith('*'):
                    current_section = line.upper()
                    continue

                tokens = [t.strip() for t in line.split(',')]
                try:
                    if current_section == '*CONTROL':
                        self.control_node = int(tokens[0])

                    elif current_section == '*NODES':
                        nid = int(tokens[0])
                        x, y = float(tokens[1]), float(tokens[2])
                        self.nodes[nid] = (x, y)

                    elif current_section == '*FIXITIES':
                        nid = int(tokens[0])
                        dof = (int(tokens[1]), int(tokens[2]), int(tokens[3]))
                        self.fixities[nid] = dof
                        if dof[0] == 1 or dof[1] == 1:
                            self.base_nodes.append(nid)

                    elif current_section == '*MATERIALS':
                        mat_id = int(tokens[0])
                        mat_type = tokens[1].upper()
                        if mat_type == 'STEEL01':
                            fy = float(tokens[2]) * MPa
                            e = float(tokens[3]) * GPa
                            b = float(tokens[4])
                            self.mat.register(mat_id, 'Steel01', {'Fy': fy, 'E': e, 'b': b})
                        elif mat_type == 'STEEL02':
                            fy = float(tokens[2]) * MPa
                            e = float(tokens[3]) * GPa
                            b = float(tokens[4])
                            r0 = float(tokens[5]) if len(tokens) > 5 else 18.5
                            cr1 = float(tokens[6]) if len(tokens) > 6 else 0.925
                            cr2 = float(tokens[7]) if len(tokens) > 7 else 0.15
                            a1 = float(tokens[8]) if len(tokens) > 8 else 0.0
                            a2 = float(tokens[9]) if len(tokens) > 9 else 1.0
                            a3 = float(tokens[10]) if len(tokens) > 10 else 0.0
                            a4 = float(tokens[11]) if len(tokens) > 11 else 1.0
                            self.mat.register(mat_id, 'Steel02', {
                                'Fy': fy, 'E': e, 'b': b,
                                'R0': r0, 'cR1': cr1, 'cR2': cr2,
                                'a1': a1, 'a2': a2, 'a3': a3, 'a4': a4
                            })
                        elif mat_type == 'CONCRETE01':
                            fpc = float(tokens[2]) * MPa
                            epsc0 = float(tokens[3])
                            fpcu = float(tokens[4]) * MPa
                            epsu = float(tokens[5])
                            self.mat.register(mat_id, 'Concrete01', {
                                'fpc': fpc, 'epsc0': epsc0, 'fpcu': fpcu, 'epsU': epsu
                            })
                        else:
                            raise ValueError(f"Unknown material type '{mat_type}'")

                    elif current_section == '*SECTIONS':
                        name = tokens[0]
                        sec_type = tokens[1].upper()
                        sec_id = int(tokens[2])
                        mat_id = int(tokens[3])
                        if sec_type == 'BOX':
                            b, h, t = float(tokens[4]) * mm, float(tokens[5]) * mm, float(tokens[6]) * mm
                            self.sec_prop.register(name, BoxSection(sec_id, mat_id, b, h, t))
                        elif sec_type == 'W':
                            d, bf, tf, tw = (float(tokens[4]) * mm, float(tokens[5]) * mm,
                                             float(tokens[6]) * mm, float(tokens[7]) * mm)
                            self.sec_prop.register(name, WSection(sec_id, mat_id, d, bf, tf, tw))
                        else:
                            raise ValueError(f"Unknown section type '{sec_type}'")

                    elif current_section == '*ELEMENTS':
                        ele_id = int(tokens[0])
                        ni, nj = int(tokens[1]), int(tokens[2])
                        sec_name = tokens[3]
                        transf = tokens[4]
                        n_pts = int(tokens[5]) if len(tokens) > 5 else 5
                        self.elements.append({
                            'id': ele_id, 'ni': ni, 'nj': nj,
                            'sec_name': sec_name, 'transf': transf, 'n_pts': n_pts
                        })

                    elif current_section == '*MASSES':
                        nid = int(tokens[0])
                        self.masses[nid] = (float(tokens[1]), float(tokens[2]), float(tokens[3]))

                    elif current_section == '*GRAVITY_LOADS':
                        nid = int(tokens[0])
                        self.gravity_loads[nid] = (float(tokens[1]) * kN, float(tokens[2]) * kN, float(tokens[3]) * kN)

                    elif current_section == '*LATERAL_PATTERN':
                        nid = int(tokens[0])
                        self.lateral_pattern[nid] = (float(tokens[1]), float(tokens[2]), float(tokens[3]))

                except Exception as e:
                    raise ValueError(f"Error parsing {filepath.name} at line {line_num}: '{raw_line.strip()}' -> {e}")

        if self.nodes:
            y_vals = [coord[1] for coord in self.nodes.values()]
            self.total_height = max(y_vals) - min(y_vals)

    def build(self):
        """依據解析後的資料建立 OpenSees 模型"""
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        # 1. 材料與纖維斷面
        self.mat.build_materials()
        self.sec_prop.build_all()

        # 2. 節點
        for nid, (x, y) in self.nodes.items():
            ops.node(nid, x, y)

        # 3. 邊界條件
        for nid, dof in self.fixities.items():
            ops.fix(nid, *dof)

        # 4. 質量
        for nid, m_vec in self.masses.items():
            ops.mass(nid, *m_vec)

        # 5. 幾何轉換
        transf_map = {'Linear': 1, 'PDelta': 2, 'Corotational': 3}
        ops.geomTransf('Linear', 1)
        ops.geomTransf('PDelta', 2)
        ops.geomTransf('Corotational', 3)

        # 6. 非線性梁柱元素
        created_integ = set()
        for ele in self.elements:
            sec_id = self.sec_prop.get_id(ele['sec_name'])
            n_pts = ele['n_pts']
            integ_tag = sec_id * 100 + n_pts
            if integ_tag not in created_integ:
                ops.beamIntegration('Lobatto', integ_tag, sec_id, n_pts)
                created_integ.add(integ_tag)
            t_tag = transf_map.get(ele['transf'])
            if t_tag is None:
                raise ValueError(f"Unsupported geomTransf: {ele['transf']}")
            ops.element('forceBeamColumn', ele['id'], ele['ni'], ele['nj'], t_tag, integ_tag)

    def plot(self, show_node_id: bool = True, show_ele_id: bool = True,
             show_loads: bool = True, save_dir=None, show: bool = True,
             filename: str = "Structure_Model.png"):
        """繪製 2D 構架拓撲與側向力分佈"""
        if not self.nodes:
            raise ValueError("No nodes found. Please call read_ipt() first.")

        fig, ax = plt.subplots(figsize=(9, 6))

        sec_names = list({ele['sec_name'] for ele in self.elements})
        cmap = plt.get_cmap('tab10')
        color_map = {name: cmap(i % 10) for i, name in enumerate(sec_names)}
        plotted_labels = set()

        for ele in self.elements:
            ni, nj = ele['ni'], ele['nj']
            xi, yi = self.nodes[ni]
            xj, yj = self.nodes[nj]
            sec_name = ele['sec_name']
            label = sec_name if sec_name not in plotted_labels else None
            plotted_labels.add(sec_name)
            ax.plot([xi, xj], [yi, yj], color=color_map[sec_name], linewidth=2.8, label=label, zorder=2)

            if show_ele_id:
                xm, ym = (xi + xj) / 2.0, (yi + yj) / 2.0
                ax.text(xm, ym, f"E{ele['id']}", color='dimgray', fontsize=8,
                        ha='center', va='center',
                        bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='none', alpha=0.75))

        for nid, (x, y) in self.nodes.items():
            if nid in self.fixities:
                ax.scatter(x, y, marker='^', s=130, color='crimson', edgecolors='black', zorder=4,
                           label='Support' if 'Support' not in plotted_labels else None)
                plotted_labels.add('Support')
            elif nid == self.control_node:
                ax.scatter(x, y, marker='*', s=180, color='gold', edgecolors='black', zorder=5,
                           label='Control Node' if 'Control Node' not in plotted_labels else None)
                plotted_labels.add('Control Node')
            else:
                ax.scatter(x, y, marker='o', s=45, color='black', zorder=3)

            if show_node_id:
                ax.text(x + 0.12, y + 0.12, f"N{nid}", color='navy', fontsize=9, fontweight='bold')

        if show_loads and self.lateral_pattern:
            max_fx = max((abs(v[0]) for v in self.lateral_pattern.values()), default=1.0)
            x_span = max(c[0] for c in self.nodes.values()) - min(c[0] for c in self.nodes.values())
            arrow_scale = (x_span * 0.15) / max_fx if max_fx > 0 else 1.0
            for nid, (fx, fy, _) in self.lateral_pattern.items():
                if abs(fx) > 0:
                    x, y = self.nodes[nid]
                    dx = fx * arrow_scale
                    ax.annotate('', xy=(x, y), xytext=(x - dx, y),
                                arrowprops=dict(arrowstyle='->', color='darkorange', lw=2.0))
                    ax.text(x - dx, y + 0.1, f"{fx:.3f}", color='darkorange', fontsize=8, ha='right')

        ax.set_aspect('equal', adjustable='datalim')
        ax.margins(0.15)
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        ax.set_title("2D Frame Topology & Lateral Load Pattern")
        ax.legend(loc='upper right')
        fig.tight_layout()

        if save_dir is not None:
            save_dir = Path(save_dir)
            save_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_dir / filename, dpi=600, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig)
    def run_eigen(self, num_modes: int = 3, solver: str = "-genBandArpack") -> list[float]:
        """
        執行模態分析（特徵值分析），並直接回傳各階自然週期列表 (sec)。
        """
        eigenvalues = ops.eigen(solver, num_modes)

        if len(eigenvalues) < num_modes:
            raise RuntimeError(
                f"模態分析失敗：預期求解 {num_modes} 個模態，僅求得 {len(eigenvalues)} 個。"
                "請檢查節點質量 (ops.mass) 與邊界條件 (ops.fix) 是否設定正確。"
            )

        periods = []
        modal_details = []

        for mode_idx, lam in enumerate(eigenvalues, start=1):
            if lam <= 0:
                raise ValueError(f"第 {mode_idx} 階特徵值異常 (lambda = {lam})，請檢查剛度矩陣是否正定。")

            omega = math.sqrt(lam)             # 圓頻率 (rad/s)
            freq = omega / (2.0 * math.pi)     # 自然頻率 (Hz)
            period = 1.0 / freq                # 自然週期 (sec)

            periods.append(period)
            modal_details.append({
                "mode": mode_idx,
                "eigenvalue": lam,
                "omega": omega,
                "frequency": freq,
                "period": period
            })

        # 完整資訊仍保留在實例屬性中，供後續查閱或側推分析使用
        self.periods = periods
        self.modal_results = modal_details
        
        ops.record()
        
        # 直接回傳純浮點數列表，對齊 test.py 的 periods[0], periods[1], periods[2]
        return periods

    def _calculate_fema273_k(self, T1: float) -> float:
        """依據 FEMA 273 計算垂直分佈指數 k"""
        if T1 <= 0.5:
            return 1.0
        elif T1 >= 2.5:
            return 2.0
        return 1.0 + (T1 - 0.5) / 2.0

    def apply_fema273_lateral_pattern(
        self,
        lateral_load=None,
        T1: float = None,
        pattern_tag: int = 10,   # 改為 10，避開 test.py 的 tag 2
        ts_tag: int = 10,        # 改為 10，避開 test.py 的 tag 2
        distribute_to_all_floor_nodes: bool = True
    ) -> dict[int, float]:
        """
        依據 FEMA 273 規範與 .ipt 的 lateral_load 屬性建立樓層側向載重關係。

        參數:
            lateral_load: 從 .ipt 讀取的樓層地震重量/側載資料。
                          若未傳入，自動嘗試讀取 self.lateral_load 或 self.properties['lateral_load']。
            T1 (float): 第一模態自然週期 (sec)，若未傳入則讀取 self.periods[0]。
            pattern_tag (int): OpenSees 側推載重模式編號 (預設 2，接續在重力載重 pattern 1 之後)。
            ts_tag (int): OpenSees Linear TimeSeries 編號 (預設 2)。
            distribute_to_all_floor_nodes (bool): 
                True  -> 將該樓層計算出的側力 Fx 平均分配給該高程的所有節點（適合無剛性隔板設定時，避免局部梁軸壓）。
                False -> 僅將該樓層側力 Fx 集中施加在該層第一個節點（Grid 1，適合已設定 equalDOF 剛性樓板時）。
        """
        # 1. 取得第一模態週期 T1 與指數 k
        if T1 is None:
            if hasattr(self, "periods") and self.periods:
                T1 = self.periods[0]
            else:
                raise ValueError("請先執行 run_eigen() 取得第一模態週期，或手動傳入 T1。")
        k = self._calculate_fema273_k(T1)

        # 2. 自動從模型實例中尋找 .ipt 解析出的 lateral_load
        if lateral_load is None:
            if hasattr(self, "lateral_load"):
                lateral_load = self.lateral_load
            elif hasattr(self, "properties") and isinstance(self.properties, dict):
                lateral_load = self.properties.get("lateral_load")
            elif hasattr(self, "data") and isinstance(self.data, dict):
                lateral_load = self.data.get("lateral_load")

        if lateral_load is None:
            raise ValueError("找不到 lateral_load 資料，請確認 .ipt 讀取結果並傳入 apply_fema273_lateral_pattern。")

        # 3. 整理模型所有節點的高程 (Y座標)，建立 {高程: [node_tags...]} 映射表
        all_nodes = ops.getNodeTags()
        elev_to_nodes = {}
        for tag in all_nodes:
            y = round(ops.nodeCoord(tag)[1], 4)  # 取小數點後4位避免浮點數誤差
            elev_to_nodes.setdefault(y, []).append(tag)

        sorted_elevs = sorted(elev_to_nodes.keys())
        base_elev = sorted_elevs[0]  # 1F 基礎高程 (EL. +0.00)
        story_elevs = sorted_elevs[1:]  # 2F ~ RF 共 8 個樓層高程

        # 4. 解析 lateral_load 格式並配對樓層高度 h_x 與地震重量權重 w_x
        # 支援兩種常見的 .ipt 格式：
        # 格式 A: dict -> {高程或樓層索引或節點編號: w_x}
        # 格式 B: list -> [w_2F, w_3F, ..., w_RF] (長度為 8 的樓層重量列表)
        floor_data = []  # [(elev, h_x, w_x, target_nodes), ...]

        if isinstance(lateral_load, (list, tuple)):
            if len(lateral_load) != len(story_elevs):
                raise ValueError(f"lateral_load 長度 ({len(lateral_load)}) 與樓層數 ({len(story_elevs)}) 不符。")
            for elev, w_x in zip(story_elevs, lateral_load):
                h_x = elev - base_elev
                nodes = sorted(elev_to_nodes[elev])
                floor_data.append((elev, h_x, float(w_x), nodes))

        elif isinstance(lateral_load, dict):
            # 判斷 dict 的 key 是「節點編號」還是「樓層編號/高程」
            sample_key = next(iter(lateral_load.keys()))
            if sample_key in all_nodes:
                # Key 為節點編號：按節點所在高程加總該層權重
                elev_weight = {}
                for tag, w_x in lateral_load.items():
                    y = round(ops.nodeCoord(int(tag))[1], 4)
                    elev_weight[y] = elev_weight.get(y, 0.0) + float(w_x)
                for elev in sorted(elev_weight.keys()):
                    h_x = elev - base_elev
                    if h_x > 0:
                        floor_data.append((elev, h_x, elev_weight[elev], sorted(elev_to_nodes[elev])))
            else:
                # Key 為樓層名稱/編號 (如 '2F'~'RF' 或 1~8 或高程值)，依序對應 2F~RF
                if len(lateral_load) != len(story_elevs):
                    raise ValueError(f"lateral_load 樓層數 ({len(lateral_load)}) 與模型樓層數 ({len(story_elevs)}) 不符。")
                for elev, w_x in zip(story_elevs, lateral_load.values()):
                    h_x = elev - base_elev
                    nodes = sorted(elev_to_nodes[elev])
                    floor_data.append((elev, h_x, float(w_x), nodes))

        # 5. 計算各樓層 C_vx = (w_x * h_x^k) / sum(w_i * h_i^k)
        whk_list = [w_x * (h_x ** k) for (_, h_x, w_x, _) in floor_data]
        total_whk = sum(whk_list)

        # 6. 建立 OpenSees LoadPattern (總基底剪力標準化為 1.0，方便位移控制讀取基底剪力 V)
        try:
            ops.remove("loadPattern", pattern_tag)
            ops.remove("timeSeries", ts_tag)
        except Exception:
            pass

        ops.timeSeries("Linear", ts_tag)
        ops.pattern("Plain", pattern_tag, ts_tag)

        ndf = ops.getNDF()  # 2D 構架為 3 (Ux, Uy, Rz)
        applied_nodal_forces = {}
        self.story_lateral_forces = {}

        for (elev, h_x, w_x, nodes), whk in zip(floor_data, whk_list):
            c_vx = whk / total_whk
            story_fx = c_vx * 1.0  # 標準化總基底剪力 = 1.0
            self.story_lateral_forces[elev] = story_fx

            if distribute_to_all_floor_nodes:
                nodal_fx = story_fx / len(nodes)
                for n_tag in nodes:
                    load_vec = [0.0] * ndf
                    load_vec[0] = nodal_fx  # X 方向側向力
                    ops.load(n_tag, *load_vec)
                    applied_nodal_forces[n_tag] = nodal_fx
            else:
                master_node = nodes[0]
                load_vec = [0.0] * ndf
                load_vec[0] = story_fx
                ops.load(master_node, *load_vec)
                applied_nodal_forces[master_node] = story_fx

        self.fema273_k = k
        self.lateral_forces = applied_nodal_forces
        return applied_nodal_forces
    
    def plot_mode(
        self,
        modes: int | list[int] = [1, 2, 3],
        output_dir: str | Path = "output",
        scale: float = None,
        show: bool = False
    ):
        """
        繪製指定階數的模態振型圖並儲存至 output_dir。

        參數:
            modes (int | list[int]): 欲繪製的模態階數（例如 1 或 [1, 2, 3]）。
            output_dir (str | Path): 輸出資料夾路徑或完整檔案路徑（對齊 test.py 的第二個參數）。
            scale (float, optional): 變形放大倍率；若為 None 則自動依據建築高度計算。
            show (bool): 是否顯示圖表視窗。
        """
        if isinstance(modes, int):
            modes = [modes]

        if not hasattr(self, "periods") or len(self.periods) < max(modes):
            self.run_eigen(num_modes=max(modes))

        node_tags = ops.getNodeTags()
        ele_tags = ops.getEleTags()

        coords = {tag: ops.nodeCoord(tag) for tag in node_tags}
        xs = [c[0] for c in coords.values()]
        ys = [c[1] for c in coords.values()]
        span_x = max(xs) - min(xs) if xs else 1.0
        height_y = max(ys) - min(ys) if ys else 1.0

        fig, axes = plt.subplots(1, len(modes), figsize=(4.5 * len(modes), 6.5), squeeze=False)
        axes = axes.flatten()

        for idx, mode in enumerate(modes):
            ax = axes[idx]

            eig_vecs = {}
            max_disp = 0.0
            for tag in node_tags:
                ux = ops.nodeEigenvector(tag, mode, 1)
                uy = ops.nodeEigenvector(tag, mode, 2)
                eig_vecs[tag] = (ux, uy)
                max_disp = max(max_disp, abs(ux), abs(uy))

            mode_scale = scale if scale is not None else (0.12 * height_y / max_disp if max_disp > 0 else 1.0)

            # 1. 原始結構網格（淺灰色虛線）
            for ele in ele_tags:
                n1, n2 = ops.eleNodes(ele)
                x1, y1 = coords[n1][0], coords[n1][1]
                x2, y2 = coords[n2][0], coords[n2][1]
                ax.plot([x1, x2], [y1, y2], color="#B0B0B0", linestyle="--", linewidth=1.0, zorder=1)

            # 2. 變形後模態振型（藍色實線）
            for ele in ele_tags:
                n1, n2 = ops.eleNodes(ele)
                x1_def = coords[n1][0] + eig_vecs[n1][0] * mode_scale
                y1_def = coords[n1][1] + eig_vecs[n1][1] * mode_scale
                x2_def = coords[n2][0] + eig_vecs[n2][0] * mode_scale
                y2_def = coords[n2][1] + eig_vecs[n2][1] * mode_scale
                ax.plot([x1_def, x2_def], [y1_def, y2_def], color="#1f77b4", linewidth=2.0, zorder=2)

            # 3. 變形後節點
            def_xs = [coords[t][0] + eig_vecs[t][0] * mode_scale for t in node_tags]
            def_ys = [coords[t][1] + eig_vecs[t][1] * mode_scale for t in node_tags]
            ax.scatter(def_xs, def_ys, color="#d62728", s=18, zorder=3)

            T_val = self.periods[mode - 1]
            ax.set_title(f"Mode {mode}\n$T_{mode} = {T_val:.4f}$ s", fontsize=12, fontweight="bold")
            ax.set_xlabel("X", fontsize=10)
            if idx == 0:
                ax.set_ylabel("Elevation Y", fontsize=10)

            ax.set_aspect("equal", adjustable="box")
            ax.set_xlim(min(xs) - 0.25 * span_x, max(xs) + 0.25 * span_x)
            ax.set_ylim(min(ys) - 0.05 * height_y, max(ys) + 0.08 * height_y)
            ax.grid(True, linestyle=":", alpha=0.5)

        plt.tight_layout()

        # 處理 Path 或 str 輸出路徑：若傳入資料夾，自動命名為 Mode_1.png 或 Mode_Shapes.png
        if output_dir is not None:
            out_path = Path(output_dir)
            if out_path.suffix == "":
                os.makedirs(out_path, exist_ok=True)
                filename = f"Mode_{modes[0]}.png" if len(modes) == 1 else "Mode_Shapes.png"
                save_file = out_path / filename
            else:
                os.makedirs(out_path.parent, exist_ok=True)
                save_file = out_path
            plt.savefig(save_file, dpi=300, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig

    def apply_rigid_diaphragm_2d(self):
    #"""為 2D 構架的每個樓層（排除 1F 基礎）建立水平向剛性樓板 (equalDOF)"""
        all_nodes = ops.getNodeTags()
        elev_to_nodes = {}
        for tag in all_nodes:
            y = round(ops.nodeCoord(tag)[1], 4)
            elev_to_nodes.setdefault(y, []).append(tag)

        sorted_elevs = sorted(elev_to_nodes.keys())
    # 排除最底層 1F (index 0)，針對 2F ~ RF 建立剛性樓板
        for elev in sorted_elevs[1:]:
            nodes = sorted(elev_to_nodes[elev])
            if len(nodes) > 1:
                master = nodes[0]
                for slave in nodes[1:]:
                    ops.equalDOF(master, slave, 1)
                    
        return 0
    