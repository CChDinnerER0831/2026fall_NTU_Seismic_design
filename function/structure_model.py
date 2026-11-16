import openseespy.opensees as ops
from pathlib import Path
from function.material_property import MaterialProperty
from function.section_property import SectionProperty, WSection, BoxSection
import matplotlib.pyplot as plt

class StructureModel:
    """
    透過讀取類似 FRAME26 的 .ipt 文字檔建立複雜 OpenSees 2D 結構模型
    """

    def __init__(self, mat: MaterialProperty, sec_prop: SectionProperty):
        self.mat = mat
        self.sec_prop = sec_prop

        # 解析後的原始資料結構
        self.control_node: int = -1
        self.nodes: dict[int, tuple[float, float]] = {}
        self.fixities: dict[int, tuple[int, int, int]] = {}
        self.elements: list[dict] = []
        self.masses: dict[int, tuple[float, float, float]] = {}
        self.gravity_loads: dict[int, tuple[float, float, float]] = {}
        self.lateral_pattern: dict[int, tuple[float, float, float]] = {}

        # 記錄支承點與總高度供 PushoverAnalyzer 使用
        self.base_nodes: list[int] = []
        self.total_height: float = 0.0

    def read_ipt(self, filepath: str | Path):
        """解析 .ipt 檔案並載入幾何、斷面、元素與載重設定"""
        filepath = Path(filepath)
        if not filepath.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        current_section = None
        mm = self.mat.mm
        kN = self.mat.kN

        with open(filepath, 'r', encoding='utf-8') as f:
            for line_num, raw_line in enumerate(f, start=1):
                # 移除行內註解 (# 或 !) 與前後空白
                line = raw_line.split('#')[0].split('!')[0].strip()
                if not line:
                    continue

                # 切換關鍵字區塊
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
        """依據解析完畢的 .ipt 資料建立 OpenSees 模型"""
        ops.wipe()
        ops.model('basic', '-ndm', 2, '-ndf', 3)

        # 1. 建立材料與斷面
        self.mat.build_materials()
        self.sec_prop.build_all()

        # 2. 建立節點
        for nid, (x, y) in self.nodes.items():
            ops.node(nid, x, y)

        # 3. 邊界條件
        for nid, dof in self.fixities.items():
            ops.fix(nid, *dof)

        # 4. 節點質量
        for nid, m_vec in self.masses.items():
            ops.mass(nid, *m_vec)

        # 5. 座標轉換定義
        transf_map = {'Linear': 1, 'PDelta': 2, 'Corotational': 3}
        ops.geomTransf('Linear', 1)
        ops.geomTransf('PDelta', 2)
        ops.geomTransf('Corotational', 3)

        # 6. 積分標籤與元素建立
        created_integ = set()
        for ele in self.elements:
            sec_id = self.sec_prop.get_id(ele['sec_name'])
            n_pts = ele['n_pts']
            integ_tag = sec_id * 100 + n_pts  # 避免相同斷面不同積分點數時標籤衝突

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
        """繪製 .ipt 解析後的 2D 結構拓撲、邊界條件、斷面分佈與側推載重方向"""
        if not self.nodes:
            raise ValueError("No nodes found. Please call read_ipt() first.")

        fig, ax = plt.subplots(figsize=(9, 6))

        # 1. 依斷面名稱自動分配不同顏色繪製構件
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

        # 2. 繪製節點與支承邊界條件
        for nid, (x, y) in self.nodes.items():
            if nid in self.fixities:
                # 固定端/支承以紅色三角形標示
                ax.scatter(x, y, marker='^', s=130, color='crimson', edgecolors='black', zorder=4,
                           label='Support' if 'Support' not in plotted_labels else None)
                plotted_labels.add('Support')
            elif nid == self.control_node:
                # 側推控制節點以金色星號標示
                ax.scatter(x, y, marker='*', s=180, color='gold', edgecolors='black', zorder=5,
                           label='Control Node' if 'Control Node' not in plotted_labels else None)
                plotted_labels.add('Control Node')
            else:
                ax.scatter(x, y, marker='o', s=45, color='black', zorder=3)

            if show_node_id:
                ax.text(x + 0.12, y + 0.12, f"N{nid}", color='navy', fontsize=9, fontweight='bold')

        # 3. 繪製側向推力箭頭 (*LATERAL_PATTERN)
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