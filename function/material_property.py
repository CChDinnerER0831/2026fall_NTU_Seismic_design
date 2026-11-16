import openseespy.opensees as ops


class MaterialProperty:
    """管理單位系統與材料力學性質"""

    def __init__(self):
        self.m, self.N, self.kg, self.s = 1.0, 1.0, 1.0, 1.0
        self.acg = 9.81 * self.m / self.s**2
        self.Pa = self.N / self.m**2
        self.mm = self.m * 1.0e-3
        self.kN = self.N * 1.0e3
        self.MPa = self.Pa * 1.0e6
        self.GPa = self.Pa * 1.0e9

        # A36 鋼材預設參數
        self.steel_tag = 1
        self.E = 200.0 * self.GPa
        self.Fy = 250.0 * self.MPa
        self.b = 0.02

    def build_materials(self):
        ops.uniaxialMaterial('Steel01', self.steel_tag, self.Fy, self.E, self.b)