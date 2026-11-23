# 2026 Fall NTU Seismic Design — 非線性靜力側推分析與耐震性能評估

本專案為臺灣大學土木工程學系「耐震設計（Seismic Design）」課程專題程式庫，基於 **OpenSeesPy** 與 **opsvis** 建立 2D 鋼結構抗彎構架（SMRF）非線性纖維斷面模型，並依據 **FEMA 273 係數法（Coefficient Method）** 進行設計地震（DBE）與最大考量地震（MCE）下之目標位移計算與構件塑鉸容量檢核。

---

## 📁 專案目錄結構

```text
2026fall_NTU_Seismic_design/
├── README.md                     # 專案說明文件
├── task1.py                      # 主執行程式：整合側推分析、FEMA 273 目標位移與容量檢核
├── function/                     # 核心結構分析模組
│   ├── material_property.py      # 鋼材單軸非線性材料參數與單位轉換管理
│   ├── section_property.py       # 纖維斷面定義（W/H 型鋼梁、BOX 中空箱型柱）與網格繪製
│   ├── structure_model.py        # 讀取 .ipt 建模檔、建立 2D 構架、模態分析與剛性樓板設定
│   └── pushover_analyze.py       # 靜載重分析與位移控制非線性側推分析器
├── task/                         # 結構輸入檔目錄
│   └── 2.ipt                     # 建築幾何拓撲、節點質量、斷面配置與載重輸入檔
└── result/                       # 分析結果與視覺化圖表輸出目錄
    ├── All_Sections.png          # 全棟梁柱纖維斷面網格總覽圖
    ├── Structure_Model.png       # 2D 結構模型拓撲與邊界條件示意圖
    ├── Mode_1.png                # 第一自然振態（基本週期 T1）變形圖
    ├── Pushover_DBE_MCE.png      # 基底剪力－屋頂位移側推曲線（含 DBE/MCE 目標位移標記）
    ├── Plastic_Hinges_LS.png     # DBE（生命安全 LS 性能水準）塑鉸分佈圖
    ├── Plastic_Hinges_CP.png     # MCE（防止倒塌 CP 性能水準）塑鉸分佈圖
    ├── Capacity_Check_LS.csv     # DBE 目標位移下各桿件塑鉸轉角 D/C 檢核表
    └── Capacity_Check_CP.csv     # MCE 目標位移下各桿件塑鉸轉角 D/C 檢核表
```

---

## ⚙️ 環境需求與安裝

建議使用 **Python 3.10 – 3.12** 並建立獨立虛擬環境執行（Windows 環境請確保已安裝 Microsoft Visual C++ Redistributable 並將專案放置於純英文路徑下）：

```bash
# 建立並啟動虛擬環境
python -m venv venv
# Windows PowerShell:
.\venv\Scripts\Activate.ps1

# 安裝必要套件
pip install openseespy opsvis numpy pandas matplotlib
```

---

## 🚀 執行方式

直接執行主程式 `task1.py` 即可自動完成建模、模態分析、非線性側推分析、目標位移計算及所有圖表與報表輸出：

```bash
python task1.py
```

---

## 🔬 分析流程與理論依據

1. **纖維斷面與非線性桿件建模（`function/`）**
   * **梁構件（`WSection`）**：採用 3 塊 `quad` 纖維網格（上翼板、腹板、下翼板）建立 W/H 型鋼斷面。
   * **柱構件（`BoxSection`）**：採用 4 塊 `quad` 纖維網格（上下翼板、左右腹板）建立中空箱型（RHS/SHS）斷面，避免角隅纖維重疊。
   * 支援外部注入 `ops` 實例與自動化斷面幾何性質（$A, I_z, I_y, Z_z$）計算。

2. **模態分析與側向力分佈**
   * 求解結構特徵值取得基本振動週期 $T_1$，並設定 2D 剛性樓板（Rigid Diaphragm）束制同樓層水平自由度。

3. **非線性側推歷程追蹤（`Task1Evaluator.run_pushover_with_history`）**
   * 先施加重力靜載重，再以位移控制法（`DisplacementControl`）進行側推，具備自動切換收斂演算法（Newton / ModifiedNewton / KrylovNewton / BFGS）之機制。
   * 同步紀錄每一步的桿端內力（`eleForce`）與基本弦轉角變形（`basicDeformation`）。

4. **FEMA 273 目標位移與 Table 5-4 構件容量檢核**
   * **雙線性化與係數法**：由容量曲線求得有效勁度 $K_e$、降伏剪力 $V_y$ 與有效週期 $T_e$，計算修正係數 $C_0, C_1, C_2, C_3$ 與目標位移 $\delta_t$（涵蓋台北二區 DBE 與 MCE 反應譜參數）。
   * **塑鉸轉角檢核**：依據 FEMA 273 Eq. (5-1)、Eq. (5-2) 計算梁柱降伏轉角 $\theta_y$ 與軸力折減塑矩 $M_{pc}$，並依寬厚比（$b_f / 2t_f$ 與 $b / 2t$）判定緊密斷面與容許塑性轉角，輸出完整 `D/C` 檢核表（CSV）與塑鉸分佈圖。