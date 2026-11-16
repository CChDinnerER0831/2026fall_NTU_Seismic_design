import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from pathlib import Path
from abc import ABC, abstractmethod


class BaseSection(ABC):
    """抽象斷面基底類別"""

    def __init__(self, secID: int, matID: int):
        self.secID = secID
        self.matID = matID
        self.patches = []
        self.area = 0.0
        self.Iz = 0.0
        self.Iy = 0.0

    @abstractmethod
    def _validate(self):
        pass

    @abstractmethod
    def _compute_geometry(self):
        pass

    def build(self):
        ops.section('Fiber', self.secID)
        for p in self.patches:
            ops.patch(p[1], *p[2:])
        return self.area

    def plot(self, ax=None, save_dir=None, show=True, title=None):
        created_fig = False
        if ax is None:
            fig, ax = plt.subplots()
            created_fig = True
        else:
            fig = ax.figure

        for p in self.patches:
            _, _, _, sub_ij, sub_jk, yI, zI, yJ, zJ, yK, zK, yL, zL = p
            u = np.linspace(0.0, 1.0, sub_ij + 1)
            v = np.linspace(0.0, 1.0, sub_jk + 1)
            U, V = np.meshgrid(u, v, indexing='ij')

            Z = (1 - U) * (1 - V) * zI + U * (1 - V) * zJ + U * V * zK + (1 - U) * V * zL
            Y = (1 - U) * (1 - V) * yI + U * (1 - V) * yJ + U * V * yK + (1 - U) * V * yL

            for i in range(sub_ij + 1):
                ax.plot(Z[i, :], Y[i, :], color='steelblue', linewidth=0.5, alpha=0.7)
            for j in range(sub_jk + 1):
                ax.plot(Z[:, j], Y[:, j], color='steelblue', linewidth=0.5, alpha=0.7)

            outline = Polygon(
                [(zI, yI), (zJ, yJ), (zK, yK), (zL, yL)],
                closed=True, fill=True, facecolor='lightsteelblue',
                edgecolor='black', linewidth=1.2, alpha=0.35
            )
            ax.add_patch(outline)

        ax.autoscale_view()
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.set_title(title or f"Section {self.secID} ({self.__class__.__name__})")
        ax.set_xlabel("z (m)")
        ax.set_ylabel("y (m)")

        if created_fig:
            fig.tight_layout()
            if save_dir is not None:
                save_dir = Path(save_dir)
                save_dir.mkdir(parents=True, exist_ok=True)
                fig.savefig(save_dir / f"Section_{self.secID}.png", dpi=600, bbox_inches="tight")
            if show:
                plt.show()
            else:
                plt.close(fig)


class WSection(BaseSection):
    """W / H 型鋼斷面"""

    def __init__(self, secID: int, matID: int, d: float, bf: float, tf: float, tw: float,
                 nfdw: int = 20, nftw: int = 2, nfbf: int = 20, nftf: int = 2):
        super().__init__(secID, matID)
        self.d, self.bf, self.tf, self.tw = d, bf, tf, tw
        self.nfdw, self.nftw, self.nfbf, self.nftf = nfdw, nftw, nfbf, nftf
        self._validate()
        self._compute_geometry()

    def _validate(self):
        if min(self.d, self.bf, self.tf, self.tw) <= 0 or self.d <= 2.0 * self.tf or self.bf < self.tw:
            raise ValueError("Invalid H-section dimensions")

    def _compute_geometry(self):
        dw = self.d - 2.0 * self.tf
        self.area = self.d * self.bf - dw * (self.bf - self.tw)
        self.Iz = (self.bf * self.d**3 - (self.bf - self.tw) * dw**3) / 12.0
        self.Iy = (2.0 * self.tf * self.bf**3 + dw * self.tw**3) / 12.0

        y1, y2, y3, y4 = -self.d / 2.0, -dw / 2.0, dw / 2.0, self.d / 2.0
        z1, z2, z3, z4 = -self.bf / 2.0, -self.tw / 2.0, self.tw / 2.0, self.bf / 2.0

        self.patches = [
            ['patch', 'quad', self.matID, self.nfbf, self.nftf, y3, z4, y3, z1, y4, z1, y4, z4],
            ['patch', 'quad', self.matID, self.nftw, self.nfdw, y2, z3, y2, z2, y3, z2, y3, z3],
            ['patch', 'quad', self.matID, self.nfbf, self.nftf, y1, z4, y1, z1, y2, z1, y2, z4],
        ]


class BoxSection(BaseSection):
    """箱型 (RHS / SHS) 鋼斷面"""

    def __init__(self, secID: int, matID: int, b: float, h: float, t: float,
                 nfb: int = 20, nfh: int = 20, nft: int = 2):
        super().__init__(secID, matID)
        self.b, self.h, self.t = b, h, t
        self.nfb, self.nfh, self.nft = nfb, nfh, nft
        self._validate()
        self._compute_geometry()

    def _validate(self):
        if min(self.b, self.h, self.t) <= 0 or self.h <= 2.0 * self.t or self.b <= 2.0 * self.t:
            raise ValueError("Invalid Box-section dimensions")

    def _compute_geometry(self):
        bw = self.b - 2.0 * self.t
        hw = self.h - 2.0 * self.t
        self.area = self.b * self.h - bw * hw
        self.Iz = (self.b * self.h**3 - bw * hw**3) / 12.0
        self.Iy = (self.h * self.b**3 - hw * bw**3) / 12.0

        y1, y2, y3, y4 = -self.h / 2.0, -hw / 2.0, hw / 2.0, self.h / 2.0
        z1, z2, z3, z4 = -self.b / 2.0, -bw / 2.0, bw / 2.0, self.b / 2.0

        self.patches = [
            ['patch', 'quad', self.matID, self.nfb, self.nft, y3, z4, y3, z1, y4, z1, y4, z4],
            ['patch', 'quad', self.matID, self.nfb, self.nft, y1, z4, y1, z1, y2, z1, y2, z4],
            ['patch', 'quad', self.matID, self.nft, self.nfh, y2, z2, y2, z1, y3, z1, y3, z2],
            ['patch', 'quad', self.matID, self.nft, self.nfh, y2, z4, y2, z3, y3, z3, y3, z4],
        ]


class SectionProperty:
    """斷面庫管理器：負責註冊、建立與繪製全棟建築所需的斷面"""

    def __init__(self):
        self.sections: dict[str, BaseSection] = {}

    def register(self, name: str, section: BaseSection) -> BaseSection:
        if not isinstance(section, BaseSection):
            raise TypeError("Section must inherit from BaseSection")
        self.sections[name] = section
        return section

    def get_id(self, name: str) -> int:
        if name not in self.sections:
            raise KeyError(f"Section '{name}' not found in SectionProperty.")
        return self.sections[name].secID

    def build_all(self):
        for sec in self.sections.values():
            sec.build()

    def plot_all(self, save_dir=None, show=True, filename="All_Sections.png"):
        n = len(self.sections)
        if n == 0:
            return
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
        if n == 1:
            axes = [axes]

        for ax, (name, sec) in zip(axes, self.sections.items()):
            sec.plot(ax=ax, show=False, title=f"{name}\n({sec.__class__.__name__}, ID={sec.secID})")

        fig.tight_layout()
        if save_dir is not None:
            save_dir = Path(save_dir)
            save_dir.mkdir(parents=True, exist_ok=True)
            fig.savefig(save_dir / filename, dpi=600, bbox_inches="tight")
        if show:
            plt.show()
        else:
            plt.close(fig)