import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import openseespy.opensees as ops
import opsvis as opsv

ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

from function.material_property import MaterialProperty
from function.section_property import SectionProperty, BoxSection, WSection
from function.structure_model import StructureModel
from function.pushover_analyze import PushoverAnalyzer


class Task1Evaluator:
    """FEMA 273 係數法評估與 Table 5-4 構件容量檢核器"""

    def __init__(self, struct: StructureModel, analyzer: PushoverAnalyzer):
        self.struct = struct
        self.analyzer = analyzer
        self.ele_history = []  # 儲存每一步的桿端內力與基本變形

    def run_pushover_with_history(self, target_drift: float = 0.05, d_disp: float = 0.001):
        """執行側推分析，並同步紀錄每一步的桿端內力與基本變形 (供 DBE/MCE 檢核)"""
        target_disp = target_drift * self.struct.total_height
        n_steps = int(abs(target_disp / d_disp))

        ops.timeSeries('Linear', 2)
        ops.pattern('Plain', 2, 2)
        for nid, pat_vec in self.struct.lateral_pattern.items():
            ops.load(nid, *pat_vec)

        ctrl_node = self.struct.control_node
        ops.integrator('DisplacementControl', ctrl_node, 1, d_disp)
        ops.analysis('Static')

        disp_list = [0.0]
        shear_list = [0.0]
        self.ele_history = [self._capture_ele_state()]

        for step in range(n_steps):
            ok = ops.analyze(1)
            if ok != 0:
                for algo in [('ModifiedNewton', '-initial'), ('KrylovNewton',), ('BFGS',)]:
                    ops.algorithm(*algo)
                    ok = ops.analyze(1)
                    if ok == 0:
                        ops.algorithm('Newton')
                        break
            if ok != 0:
                print(f"Pushover stopped early at step {step} (disp = {disp_list[-1]:.4f} m)")
                break

            ops.reactions()
            vb = -sum(ops.nodeReaction(n, 1) for n in self.struct.base_nodes)
            u_roof = ops.nodeDisp(ctrl_node, 1)

            disp_list.append(u_roof)
            shear_list.append(vb)
            self.ele_history.append(self._capture_ele_state())

        self.analyzer.roof_disp = np.array(disp_list)
        self.analyzer.base_shear = np.array(shear_list)
        return self.analyzer.roof_disp, self.analyzer.base_shear

    def _capture_ele_state(self) -> dict:
        """擷取當前分析步所有元素的內力、變形與節點位移"""
        state = {'forces': {}, 'deforms': {}, 'disp': {}}
        for nid in self.struct.nodes:
            state['disp'][nid] = ops.nodeDisp(nid)
        for ele in self.struct.elements:
            eid = ele['id']
            state['forces'][eid] = ops.eleForce(eid)
            # 2D forceBeamColumn basicDeformation: [eps_axial, theta_i, theta_j]
            state['deforms'][eid] = ops.eleResponse(eid, 'basicDeformation')
        return state

    def compute_target_displacement(self, T1: float, S_XS: float, T0: float, level_name: str = "DBE"):
        """
        依 FEMA 273 Coefficient Method 計算修正係數與目標位移 delta_t
        (台北二區: DBE -> S_DS=0.60, T0=1.30s; MCE -> S_MS=0.80, T0=1.30s)
        """
        disp = self.analyzer.roof_disp
        vb = self.analyzer.base_shear
        H = self.struct.total_height

        # 1. 總地震重量 W (全棟地震重量 = 2D 構架質量的 4 倍 * g)
        W_2d = sum(m[0] for m in self.struct.masses.values()) * self.struct.mat.acg

        # 2. 容量曲線雙線性化 (簡化等能量面積法)
        vb_max = np.max(vb)
        idx_max = np.argmax(vb)
        # 初始勁度 Ki (取 0.2 * Vb,max 處割線勁度)
        idx_02 = np.searchsorted(vb[:idx_max], 0.2 * vb_max)
        Ki = vb[idx_02] / disp[idx_02]

        # 降伏剪力 Vy 與有效勁度 Ke (取 0.6 * Vy 處割線勁度並調整使面積相等)
        Vy = 0.90 * vb_max
        idx_06 = np.searchsorted(vb[:idx_max], 0.6 * Vy)
        Ke = vb[idx_06] / disp[idx_06]
        dy = Vy / Ke
        Te = T1 * np.sqrt(Ki / Ke)

        # 3. 反應譜加速度 Sa(Te)
        if Te <= 0.2 * T0:
            Sa = S_XS * (0.4 + 3.0 * Te / T0)
        elif Te <= T0:
            Sa = S_XS
        else:
            Sa = (S_XS * T0) / Te

        # 4. 修正係數 C0, C1, C2, C3 (FEMA 273 Section 3.3.3.3)
        C0 = 1.35  # 8層樓抗彎構架典型第一振態參與係數 (亦可由模態向量精確計算)
        R = (Sa / (Vy / W_2d)) * C0
        if Te >= T0:
            C1 = 1.0
        else:
            C1 = min(1.5, (1.0 + (R - 1.0) * T0 / Te) / R)
        C1 = max(1.0, C1)

        # SMRF 韌性抗彎構架無明顯勁度/強度退化，取 C2 = 1.0
        C2 = 1.0

        # 考量 P-Delta 效應之 C3 (若降伏後勁度 alpha >= 0 則 C3 = 1.0)
        alpha = (vb[idx_max] - Vy) / (disp[idx_max] - dy) / Ke if disp[idx_max] > dy else 0.0
        if alpha >= 0:
            C3 = 1.0
        else:
            C3 = 1.0 + (abs(alpha) * (R - 1.0)**1.5) / Te

        # 5. 目標位移 delta_t
        delta_t = C0 * C1 * C2 * C3 * Sa * self.struct.mat.acg * (Te / (2.0 * np.pi))**2

        print(f"\n--- [{level_name}] FEMA 273 Target Displacement ---")
        print(f"  Ki = {Ki/1e6:.2f} kN/mm, Ke = {Ke/1e6:.2f} kN/mm, Vy = {Vy/1e3:.1f} kN")
        print(f"  T1 = {T1:.3f} s -> Te = {Te:.3f} s, Sa(Te) = {Sa:.3f} g")
        print(f"  C0 = {C0:.2f}, C1 = {C1:.3f}, C2 = {C2:.2f}, C3 = {C3:.3f}")
        print(f"  Target Displacement delta_t = {delta_t*100:.2f} cm (Drift = {delta_t/H*100:.2f} %)")
        return delta_t, {'Ke': Ke, 'Vy': Vy, 'dy': dy, 'Te': Te, 'C0': C0, 'C1': C1, 'C2': C2, 'C3': C3}

    def check_members_at_disp(self, target_disp: float, perf_level: str, save_dir: Path):
        """
        在目標位移處進行 FEMA 273 Table 5-4 塑鉸轉角檢核並繪製塑鉸分佈圖
        perf_level: 'LS' (for DBE) or 'CP' (for MCE)
        """
        idx = int(np.argmin(np.abs(self.analyzer.roof_disp - target_disp)))
        actual_disp = self.analyzer.roof_disp[idx]
        state = self.ele_history[idx]

        records = []
        hinge_beam, hinge_col = [], []

        for ele in self.struct.elements:
            eid = ele['id']
            ni, nj = ele['ni'], ele['nj']
            sec = self.struct.sec_prop.sections[ele['sec_name']]
            mat_info = self.struct.mat.materials.get(sec.matID)
            fy = mat_info['params']['Fy'] if mat_info else self.struct.mat.Fy
            E = mat_info['params']['E'] if mat_info else self.struct.mat.E

            xi, yi = self.struct.nodes[ni]
            xj, yj = self.struct.nodes[nj]
            dx, dy = xj - xi, yj - yi
            L = np.hypot(dx, dy)
            cos_a, sin_a = dx / L, dy / L

            forces = state['forces'][eid]
            deforms = state['deforms'][eid]  # [eps, theta_i, theta_j]
            p_axial = max(0.0, forces[0] * cos_a + forces[1] * sin_a)

            mp = sec.Zz * fy
            py = sec.area * fy
            is_col = isinstance(sec, BoxSection)

            if is_col:
                # 柱降伏轉角 Eq. (5-2) 與軸力折減 Mpc
                m_cap = min(mp, 1.18 * mp * (1.0 - p_axial / py))
                theta_y = (sec.Zz * fy * L) / (6.0 * E * sec.Iz) * (1.0 - p_axial / py)
                slenderness = sec.b / (2.0 * sec.t)  # 作業提示: box columns use b/(2t)
                limit_slend = 52.0 / np.sqrt(fy / self.struct.mat.MPa / 6.895)  # 52/sqrt(Fy_ksi)
                is_compact = slenderness <= limit_slend
                # FEMA 273 Table 5-4 柱塑性轉角倍數 (簡化對應緊密/非緊密斷面)
                mult = (6.0 if is_compact else 1.0) if perf_level == 'LS' else (8.0 if is_compact else 1.25)
            else:
                # 梁降伏轉角 Eq. (5-1)
                m_cap = mp
                theta_y = (sec.Zz * fy * L) / (6.0 * E * sec.Iz)
                slenderness = sec.bf / (2.0 * sec.tf)
                limit_slend = 52.0 / np.sqrt(fy / self.struct.mat.MPa / 6.895)
                is_compact = slenderness <= limit_slend
                # FEMA 273 Table 5-4 梁塑性轉角倍數
                mult = (6.0 if is_compact else 2.0) if perf_level == 'LS' else (8.0 if is_compact else 3.0)

            allow_rot = mult * theta_y

            # 檢查 i 端與 j 端 (basicDeformation 為總弦轉角，扣除 theta_y 得塑性轉角)
            for end_idx, (end_name, m_end, th_tot) in enumerate([
                ('I', abs(forces[2]), abs(deforms[1])),
                ('J', abs(forces[5]), abs(deforms[2]))
            ]):
                th_p = max(0.0, th_tot - theta_y)
                is_hinge = m_end >= 0.98 * m_cap
                dc_ratio = th_p / allow_rot if allow_rot > 0 else 0.0

                records.append({
                    'EleID': eid, 'Section': ele['sec_name'], 'Type': 'Col' if is_col else 'Beam',
                    'End': end_name, 'P_kN': round(p_axial / 1e3, 1),
                    'M_kNm': round(m_end / 1e3, 1), 'Mcap_kNm': round(m_cap / 1e3, 1),
                    'Hinge': is_hinge, 'theta_y': round(theta_y, 5),
                    'theta_p': round(th_p, 5), 'Allow_theta_p': round(allow_rot, 5),
                    'D/C': round(dc_ratio, 3), 'Status': 'PASS' if dc_ratio <= 1.0 else 'NG'
                })

                if is_hinge:
                    offset = 0.04 if end_name == 'I' else 0.96
                    hx = xi + offset * dx
                    hy = yi + offset * dy
                    (hinge_col if is_col else hinge_beam).append((hx, hy))

        # 輸出 CSV 檢核報表
        df = pd.DataFrame(records)
        csv_path = save_dir / f"Capacity_Check_{perf_level}.csv"
        df.to_csv(csv_path, index=False, encoding='utf-8-sig')
        print(f"  [{perf_level}] Capacity check saved to {csv_path} (Max D/C = {df['D/C'].max():.3f})")

        # 繪製塑鉸分佈圖 (結合 opsvis)
        fig, ax = plt.subplots(figsize=(9, 6))
        opsv.plot_model(node_labels=0, element_labels=0, ax=ax)
        if hinge_beam:
            bx, by = zip(*hinge_beam)
            ax.scatter(bx, by, s=85, color='dodgerblue', edgecolors='black', zorder=5, label='Beam Hinge ($M \geq M_p$)')
        if hinge_col:
            cx, cy = zip(*hinge_col)
            ax.scatter(cx, cy, s=105, color='crimson', edgecolors='black', zorder=6, label='Column Hinge ($M \geq M_{pc}$)')
        ax.set_title(f"Plastic Hinge Distribution - {perf_level} ($\\delta_t = {actual_disp*100:.1f}$ cm)")
        ax.set_xlabel("X (m)")
        ax.set_ylabel("Y (m)")
        if hinge_beam or hinge_col:
            ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(save_dir / f"Plastic_Hinges_{perf_level}.png", dpi=600, bbox_inches="tight")
        plt.close(fig)


def main():
    ipt_file = ROOT_DIR / 'program1' / "task_Model3.ipt"
    output_dir = ROOT_DIR / 'program1' / 'output_model3'
    output_dir.mkdir(parents=True, exist_ok=True)

    # 1. 初始化材料、斷面並建立 Model 3
    mat = MaterialProperty()
    sec_prop = SectionProperty()
    struct = StructureModel(mat, sec_prop)
    struct.read_ipt(ipt_file)
    struct.build()

    # 2. 模態分析求基本週期 T1 並依 FEMA 273 更新側力分佈
    periods = struct.run_eigen(num_modes=3)
    T1 = periods[0]
    print(f"Natural Periods: T1 = {periods[0]:.4f} s, T2 = {periods[1]:.4f} s, T3 = {periods[2]:.4f} s")
    #struct.apply_fema273_lateral_pattern(T1)

    # 繪製模型拓撲與第一振態圖 (使用 opsvis)
    struct.plot(save_dir=output_dir, show=False)
    struct.plot_mode(1,output_dir, show=False)
    struct.apply_rigid_diaphragm_2d()
    # 3. 執行靜載重與非線性側推分析 (紀錄桿端反應歷程)
    analyzer = PushoverAnalyzer(struct)
    analyzer.run_gravity(steps=10)

    
    evaluator = Task1Evaluator(struct, analyzer)
    disp, vb = evaluator.run_pushover_with_history(target_drift=0.04, d_disp=1.0 * mat.mm)
    print(f"Maximum Base Shear Vb,max = {vb.max() / mat.kN:.2f} kN")

    # 4. 計算 DBE 與 MCE 目標位移 (台北二區參數)
    dt_dbe, _ = evaluator.compute_target_displacement(T1, S_XS=0.60, T0=1.30, level_name="DBE")
    dt_mce, _ = evaluator.compute_target_displacement(T1, S_XS=0.80, T0=1.30, level_name="MCE")

    # 5. 繪製含 DBE / MCE 目標位移標記的側推曲線
    fig, ax = plt.subplots(figsize=(8, 5))
    
    ax.plot(disp, vb / mat.kN, color='firebrick', lw=2.0, label='Pushover Curve (Model 3)')
    ax.axvline(dt_dbe / struct.total_height * 100.0, color='royalblue', ls='--', lw=1.5, label=f'DBE $\\delta_t$ ({dt_dbe*100:.1f} cm)')
    ax.axvline(dt_mce / struct.total_height * 100.0, color='darkorange', ls='--', lw=1.5, label=f'MCE $\\delta_t$ ({dt_mce*100:.1f} cm)')
    ax.grid(True, ls='--', alpha=0.4)
    ax.set_xlabel("Roof Drift (m)")
    ax.set_ylabel("Base Shear (kN)")
    ax.set_title("Model 3 Pushover Curve & FEMA 273 Target Displacements")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "Pushover_DBE_MCE.png", dpi=600, bbox_inches="tight")
    plt.close(fig)

    # 6. 執行 DBE (LS) 與 MCE (CP) 塑鉸分佈圖與 Table 5-4 桿件容量檢核
    evaluator.check_members_at_disp(dt_dbe, perf_level="LS", save_dir=output_dir)
    evaluator.check_members_at_disp(dt_mce, perf_level="CP", save_dir=output_dir)
    print("\nAll Task 1 evaluations completed successfully!")

    # 7. 繪製所有斷面圖 (若斷面類別有 plot 方法)    

    sec_output_dir = output_dir / "sections"
    sec_output_dir.mkdir(parents=True, exist_ok=True)
    if hasattr(sec_prop, "plot_sections"):
        sec_prop.plot_sections(save_dir=sec_output_dir, show=False)
    else:
        for sec_name, sec in sec_prop.sections.items():
            if hasattr(sec, "plot"):
                sec.save_dir = sec_output_dir
                sec.show = False
                sec.plot()
            elif hasattr(sec, "fib_sec") and sec.fib_sec:
                opsv.plot_fiber_section(sec.fib_sec, fillflag=0)
                fig = plt.gcf()
                plt.axis("equal")
                plt.title(f"Fiber Section {sec.secID} ({sec_name})")
                fig.tight_layout()
                fig.savefig(sec_output_dir / f"Section_{sec.secID}_{sec_name}.png", dpi=600, bbox_inches="tight")
                plt.close(fig)

if __name__ == "__main__":
    main()