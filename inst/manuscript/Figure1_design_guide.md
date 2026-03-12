# Figure 1: CellDeathAnalysis Package Workflow
# 
# 这个图需要用BioRender、Figma、AI或PPT手动绘制
# 以下是设计说明和元素

## 图表布局建议

```
┌─────────────────────────────────────────────────────────────────────────┐
│                     CellDeathAnalysis Workflow                          │
├─────────────────────────────────────────────────────────────────────────┤
│                                                                         │
│   ┌─────────┐     ┌─────────────┐     ┌────────────┐     ┌──────────┐  │
│   │  Input  │ ──► │  Gene Sets  │ ──► │  Scoring   │ ──► │ Analysis │  │
│   │  Data   │     │ (14 types)  │     │ (6 methods)│     │          │  │
│   └─────────┘     └─────────────┘     └────────────┘     └──────────┘  │
│        │                                                       │        │
│        │          ┌────────────────────────────────────────────┤        │
│        │          │                                            │        │
│        ▼          ▼                                            ▼        │
│   ┌─────────┐ ┌─────────┐ ┌─────────┐ ┌─────────┐     ┌──────────────┐ │
│   │RNA-seq  │ │Ferropt- │ │ssGSEA   │ │Survival │     │ Publication- │ │
│   │(bulk/sc)│ │osis     │ │GSVA     │ │Analysis │     │ Ready Plots  │ │
│   │         │ │Cupropt- │ │AUCell   │ │         │     │              │ │
│   │CSV/RDS  │ │osis     │ │Z-score  │ │Enrichmt │     │ CSV/PDF/HTML │ │
│   │Seurat   │ │Pyropt-  │ │Mean     │ │         │     │              │ │
│   │         │ │osis...  │ │Median   │ │ML Model │     │              │ │
│   └─────────┘ └─────────┘ └─────────┘ └─────────┘     └──────────────┘ │
│                                                                         │
│   ┌─────────────────────────────────────────────────────────────────┐   │
│   │                     Interactive Shiny App                        │   │
│   │   [Upload] [Score] [Visualize] [Analyze] [Export]               │   │
│   └─────────────────────────────────────────────────────────────────┘   │
│                                                                         │
└─────────────────────────────────────────────────────────────────────────┘
```

## 设计元素

### 1. 顶部标题
- "CellDeathAnalysis: Comprehensive Cell Death Pathway Analysis"
- 字体: Arial Bold, 24pt
- 颜色: 深蓝色 #2C3E50

### 2. 主要模块 (方框)
使用不同颜色区分:

| 模块 | 颜色 | 内容 |
|------|------|------|
| Input Data | 浅蓝 #E3F2FD | RNA-seq, Microarray, scRNA-seq |
| Gene Sets | 浅绿 #E8F5E9 | 14 death types, 500+ genes |
| Scoring | 浅黄 #FFF3E0 | ssGSEA, GSVA, AUCell, Z-score |
| Analysis | 浅紫 #F3E5F5 | Survival, Enrichment, ML |
| Output | 浅红 #FFEBEE | Figures, Tables, Reports |

### 3. 箭头
- 粗箭头连接主要流程
- 颜色: 深灰 #607D8B
- 宽度: 3pt

### 4. 图标建议 (可从Flaticon获取)
- Input: 📊 数据文件图标
- Gene Sets: 🧬 DNA螺旋图标
- Scoring: 🔢 计算器图标
- Analysis: 📈 图表图标
- Output: 📄 文档图标

### 5. 底部Shiny模块
- 背景: 渐变蓝色
- 图标: 桌面应用图标
- 文字: "Interactive Web Application - No Coding Required"

## 推荐工具

1. **BioRender** (推荐) - https://biorender.com
   - 专业生物医学图表
   - 有现成的模板
   - 导出高分辨率

2. **Figma** (免费) - https://figma.com
   - 专业设计工具
   - 可协作编辑

3. **PowerPoint/Keynote**
   - 简单快捷
   - 使用SmartArt

4. **draw.io** (免费) - https://draw.io
   - 在线流程图工具

## 颜色代码 (HEX)

```
蓝色系:
  深蓝: #2C3E50, #34495E
  中蓝: #3498DB, #2980B9
  浅蓝: #E3F2FD, #BBDEFB

绿色系:
  深绿: #27AE60, #229954
  浅绿: #E8F5E9, #C8E6C9

红色系:
  深红: #E74C3C, #C0392B
  浅红: #FFEBEE, #FFCDD2

黄色系:
  深黄: #F39C12, #E67E22
  浅黄: #FFF3E0, #FFE0B2

紫色系:
  深紫: #9B59B6, #8E44AD
  浅紫: #F3E5F5, #E1BEE7

灰色系:
  深灰: #607D8B, #455A64
  浅灰: #ECEFF1, #CFD8DC
```

## 文字内容

### 模块1: Input Data
- Bulk RNA-seq
- Single-cell RNA-seq  
- Microarray
- Formats: CSV, TSV, RDS
- Seurat/SCE objects

### 模块2: Gene Sets
- 14 Cell Death Types
- 500+ Curated Genes
- Updated to 2024
- Sources: FerrDb, MSigDB, KEGG

### 模块3: Scoring Methods
- ssGSEA (rank-based)
- GSVA (non-parametric)
- AUCell (single-cell)
- Z-score (fast)
- Mean/Median

### 模块4: Analysis
- Survival Analysis
  - Kaplan-Meier
  - Cox Regression
  - Time-ROC
- Enrichment
  - ORA
  - GSEA
- Machine Learning
  - Random Forest
  - LASSO/XGBoost

### 模块5: Output
- Publication-ready figures
- Statistical tables
- Interactive reports
- Downloadable results

## 示例图 (文字版)

```
╔════════════════════════════════════════════════════════════════════════════╗
║                                                                            ║
║     ╭──────────╮    ╭──────────╮    ╭──────────╮    ╭──────────╮          ║
║     │  INPUT   │───▶│   GENE   │───▶│ SCORING  │───▶│ ANALYSIS │          ║
║     │   DATA   │    │   SETS   │    │          │    │          │          ║
║     ╰──────────╯    ╰──────────╯    ╰──────────╯    ╰──────────╯          ║
║          │               │               │               │                 ║
║          ▼               ▼               ▼               ▼                 ║
║     ┌────────┐     ┌────────┐     ┌────────┐     ┌────────┐              ║
║     │RNA-seq │     │14 Types│     │ssGSEA  │     │Survival│     ╭──────╮ ║
║     │scRNA   │     │500+    │     │GSVA    │     │Enrich- │────▶│OUTPUT│ ║
║     │Seurat  │     │Genes   │     │AUCell  │     │ment    │     ╰──────╯ ║
║     └────────┘     └────────┘     │Z-score │     │ML      │              ║
║                                   └────────┘     └────────┘              ║
║                                                                            ║
║     ╔══════════════════════════════════════════════════════════════════╗  ║
║     ║              🖥️ Interactive Shiny Application                    ║  ║
║     ║         No coding required • Point-and-click interface           ║  ║
║     ╚══════════════════════════════════════════════════════════════════╝  ║
║                                                                            ║
╚════════════════════════════════════════════════════════════════════════════╝
```

## 导出设置

- 格式: PDF (矢量) 或 TIFF (300 dpi)
- 尺寸: 宽 18 cm × 高 12 cm (约 7" × 5")
- 字体: Arial或Helvetica
- 线宽: 最小1pt

---

如需帮助绘制此图,请联系: s1214844197@163.com
