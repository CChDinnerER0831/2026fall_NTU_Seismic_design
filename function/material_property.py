import openseespy.opensees as ops


class MaterialProperty:
    """管理材料單位換算與 OpenSees 單軸材料註冊與建立"""

    def __init__(self):
        self.m, self.N, self.kg, self.s = 1.0, 1.0, 1.0, 1.0
        self.acg = 9.81 * self.m / self.s**2
        self.Pa = self.N / self.m**2
        self.mm = self.m * 1.0e-3
        self.kN = self.N * 1.0e3
        self.MPa = self.Pa * 1.0e6
        self.GPa = self.Pa * 1.0e9

        # 儲存從 .ipt 讀取的材料定義 {matID: {'type': str, 'params': dict}}
        self.materials: dict[int, dict] = {}

        # 預設 A36 參數 (向下相容無 *MATERIALS 區塊的舊 .ipt)
        self.steel_tag = 1
        self.E = 200.0 * self.GPa
        self.Fy = 250.0 * self.MPa
        self.b = 0.02

    def register(self, mat_id: int, mat_type: str, params: dict):
        """註冊單軸材料定義"""
        self.materials[mat_id] = {
            'type': mat_type,
            'params': params
        }

    def build_materials(self):
        """在 OpenSees 中建立所有已註冊的材料"""
        if not self.materials:
            ops.uniaxialMaterial('Steel01', self.steel_tag, self.Fy, self.E, self.b)
            return

        for mat_id, mat_info in self.materials.items():
            m_type = mat_info['type'].upper()
            p = mat_info['params']
            if m_type == 'STEEL01':
                ops.uniaxialMaterial('Steel01', mat_id, p['Fy'], p['E'], p['b'])
            elif m_type == 'STEEL02':
                ops.uniaxialMaterial(
                    'Steel02', mat_id,
                    p['Fy'], p['E'], p['b'],
                    p['R0'], p['cR1'], p['cR2'],
                    p['a1'], p['a2'], p['a3'], p['a4']
                )
            elif m_type == 'CONCRETE01':
                ops.uniaxialMaterial('Concrete01', mat_id, p['fpc'], p['epsc0'], p['fpcu'], p['epsU'])
            else:
                raise ValueError(f"Unsupported material type: {mat_info['type']}")