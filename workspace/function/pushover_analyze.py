import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from function.structure_model import StructureModel


class PushoverAnalyzer:
    """配合 StructureModel (.ipt) 執行重力與位移控制側推分析"""

    def __init__(self, struct: StructureModel):
        self.struct = struct
        self.roof_disp = np.array([])
        self.base_shear = np.array([])

    def run_gravity(self, steps: int = 10):
        """施加 *GRAVITY_LOADS 區塊定義之靜載重"""
        ops.timeSeries('Linear', 1)
        ops.pattern('Plain', 1, 1)

        for nid, load_vec in self.struct.gravity_loads.items():
            ops.load(nid, *load_vec)

        ops.constraints('Plain')
        ops.numberer('RCM')
        ops.system('BandGeneral')
        ops.test('NormDispIncr', 1.0e-8, 25)
        ops.algorithm('Newton')
        ops.integrator('LoadControl', 1.0 / steps)
        ops.analysis('Static')

        if ops.analyze(steps) != 0:
            raise RuntimeError("Gravity analysis failed to converge!")

        ops.loadConst('-time', 0.0)

    def run_pushover(self, target_drift: float = 0.05, d_disp: float = 0.0005):
        """依據 *LATERAL_PATTERN 區塊執行位移控制側推"""
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

        self.roof_disp = np.array(disp_list)
        self.base_shear = np.array(shear_list)
        return self.roof_disp, self.base_shear

    def plot_curve(self, save_dir=None, show=True, filename="Pushover_Curve.png"):
        kN = self.struct.mat.kN
        H_tot = self.struct.total_height

        drift_pct = (self.roof_disp / H_tot) * 100.0
        vb_kN = self.base_shear / kN

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(drift_pct, vb_kN, color='firebrick', linewidth=2.0, label='Pushover Curve')
        ax.grid(True, linestyle='--', alpha=0.4)
        ax.set_xlabel("Roof Drift Ratio (%)")
        ax.set_ylabel("Base Shear (kN)")
        ax.set_title("Nonlinear Static Pushover Curve")
        ax.legend()
        fig.tight_layout()

        if save_dir is not None:
            save_dir = Path(save_dir)
            save_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_dir / filename, dpi=600, bbox_inches="tight")

        if show:
            plt.show()
        else:
            plt.close(fig)