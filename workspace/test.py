import sys
from pathlib import Path

# 僅將專案根目錄 (c:\workspace) 加入 sys.path
ROOT_DIR = Path(__file__).resolve().parent.parent
if str(ROOT_DIR) not in sys.path:
    sys.path.append(str(ROOT_DIR))

# 全專案統一使用 function.xxx 前綴
from function.material_property import MaterialProperty
from function.section_property import SectionProperty
from function.structure_model import StructureModel
from function.pushover_analyze import PushoverAnalyzer

def main():
    ipt_file = ROOT_DIR /'workspace'/ "task" / "2.ipt"
    output_dir = ROOT_DIR / 'workspace'/ "output"

    # 1. Material & Section Property
    mat = MaterialProperty()
    sec_prop = SectionProperty()

    # 2. Structure Model (讀取 .ipt 並建立模型)
    struct = StructureModel(mat, sec_prop)
    struct.read_ipt(ipt_file)
    struct.build()
    struct.plot(save_dir=output_dir, show=True)
    sec_prop.plot_all(save_dir=output_dir, show=False)

    # 3. Pushover Analyze
    analyzer = PushoverAnalyzer(struct)
    analyzer.run_gravity(steps=10)
    disp, vb = analyzer.run_pushover(target_drift=0.05, d_disp=0.5*mat.mm)

    print(f"最大基底剪力 Vb,max = {vb.max() / mat.kN:.2f} kN")
    analyzer.plot_curve(save_dir=output_dir, show=True)


if __name__ == "__main__":
    main()