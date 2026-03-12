# CellDeathAnalysis 论文投稿策略

## 两种投稿方案

根据您是否有时间运行真实TCGA数据分析，有两种投稿策略：

---

## 方案A：软件工具论文（推荐，快速发表）

### 适用情况
- 没有时间运行完整TCGA分析
- 想快速发表
- 重点展示软件功能

### 论文文件
`manuscript_software_paper.md`

### 推荐期刊

| 期刊 | 类型 | 影响因子 | 审稿周期 | 备注 |
|------|------|----------|----------|------|
| **Bioinformatics** | Application Note | ~5.8 | 4-8周 | 最经典，2页限制 |
| **BMC Bioinformatics** | Software | ~3.0 | 6-10周 | 开放获取，接受率高 |
| **JOSS** | Software | N/A | 2-4周 | 开源软件期刊，快速 |
| **F1000Research** | Software Tool | ~3.0 | 2-4周 | 开放同行评审 |
| **PeerJ** | Software | ~2.5 | 4-8周 | 快速审稿 |

### 论文结构（已完成）
```
1. Introduction (~400词)
2. Implementation (~600词)
3. Features and Usage (~400词)
4. Demonstration (~300词)
5. Comparison (~100词)
6. Conclusion (~100词)
总计: ~2000词 + 表格/图
```

### 所需图表
- Figure 1: 工作流程图（BioRender绘制）
- Figure 2: 基因集统计 + 评分演示
- Figure 3: 生存分析演示
- Figure 4: 相关性分析
- Figure 5: 富集分析演示
- Figure 6: Shiny界面截图

### 优点
✅ 无需真实数据分析  
✅ 1-2周可完成  
✅ 审稿周期短  
✅ 重点突出软件贡献  

---

## 方案B：完整研究论文（影响力更大）

### 适用情况
- 有时间运行TCGA分析（需要1-2天下载+运行）
- 想发更高影响因子期刊
- 想报告生物学发现

### 论文文件
`manuscript_draft.md`（原来的版本，需要填入真实数据）

### 推荐期刊

| 期刊 | 类型 | 影响因子 | 审稿周期 | 备注 |
|------|------|----------|----------|------|
| **Briefings in Bioinformatics** | Research | ~9.5 | 8-12周 | 生信顶刊 |
| **Nucleic Acids Research** | Database/Software | ~14.9 | 6-10周 | 如果能做数据库 |
| **Frontiers in Genetics** | Research | ~3.7 | 8-12周 | 开放获取 |
| **Computers in Biology and Medicine** | Research | ~6.7 | 8-12周 | 计算生物学 |

### 所需工作
1. 运行TCGA分析（1-2天）
2. 整理真实数据结果
3. 更新论文中的数字（样本量、p值、HR等）
4. 补充生物学讨论

### 论文结构
```
1. Introduction (~800词)
2. Materials and Methods (~1500词)
3. Results (~2000词)
   - 泛癌分析真实结果
   - 生存分析真实结果
   - 可能的生物学发现
4. Discussion (~800词)
5. Conclusion (~200词)
总计: ~5000词 + 表格/图
```

---

## 方案A详细指南（软件论文）

### Step 1: 准备图表（1-2天）

```r
# 运行图表生成脚本
source("inst/scripts/generate_paper_figures.R")

# 检查输出
list.files("figures")
```

### Step 2: 绘制Figure 1（半天）

使用BioRender或Figma，参考 `Figure1_design_guide.md`

### Step 3: 截取Shiny界面（1小时）

```r
# 启动Shiny应用
launch_death_app()

# 截图保存为Figure 6
```

### Step 4: 整理论文（1天）

1. 填写作者单位信息
2. 检查图表编号对应
3. 润色英文
4. 准备Cover Letter

### Step 5: 投稿

**推荐投稿顺序：**
1. Bioinformatics (Application Note) - 最权威
2. BMC Bioinformatics - 如果被拒
3. JOSS - 开源软件快速发表

---

## Cover Letter 模板（软件论文版）

```
Dear Editor,

We are pleased to submit our manuscript entitled "CellDeathAnalysis: 
An R Package for Comprehensive Analysis of Regulated Cell Death 
Pathways in Transcriptomic Data" for consideration as an Application 
Note in [Journal Name].

Regulated cell death (RCD) has emerged as a critical area in cancer 
research, particularly with recent discoveries of ferroptosis, 
cuproptosis, and disulfidptosis. However, researchers currently lack 
a unified computational tool for systematic analysis of these diverse 
pathways.

CellDeathAnalysis addresses this gap by providing:
• Curated gene sets for 14 cell death types (500+ genes)
• Multiple scoring methods (ssGSEA, GSVA, AUCell, Z-score)
• Integrated survival and enrichment analysis
• Single-cell RNA-seq support
• Machine learning prediction models
• Interactive Shiny application for non-programmers

The package is freely available on GitHub and will facilitate cell 
death research across the scientific community.

We believe this work is well-suited for [Journal Name] and will be 
of interest to researchers in cancer biology, bioinformatics, and 
computational biology.

Thank you for your consideration.

Sincerely,
Keran Sun
```

---

## 时间线对比

### 方案A：软件论文
```
Day 1-2:  运行generate_paper_figures.R，生成所有图
Day 3:    绘制Figure 1，截取Shiny截图
Day 4-5:  整理论文，润色英文
Day 6:    准备投稿材料，提交
---
总计: ~1周
```

### 方案B：完整研究论文
```
Day 1-3:  下载TCGA数据（可能需要多次尝试）
Day 4-5:  运行分析，处理报错
Day 6-7:  整理结果，更新论文数字
Day 8-10: 补充生物学讨论
Day 11-14: 润色，准备投稿
---
总计: ~2周
```

---

## 我的建议

### 如果您时间有限：选择方案A
- 软件论文同样有学术价值
- Bioinformatics Application Note 是生信领域认可度很高的发表形式
- 很多知名R包（如GSVA, clusterProfiler）都是以这种形式发表

### 如果您想要更大影响力：选择方案B
- 但需要确保TCGA分析能成功运行
- 可以先尝试5个癌症类型测试

### 折中方案
1. 先按方案A投稿软件论文（快速获得发表）
2. 后续有时间再用真实数据写一篇应用研究论文

---

## 文件清单

```
inst/manuscript/
├── manuscript_software_paper.md    # 方案A：软件论文（已完成）
├── manuscript_draft.md             # 方案B：研究论文模板
├── Figure1_design_guide.md         # Figure 1设计指南
└── submission_guide.md             # 本文件

inst/scripts/
├── generate_paper_figures.R        # 图表生成（使用示例数据）
├── TCGA_analysis.R                 # TCGA分析（需下载数据）
└── GEO_analysis.R                  # GEO分析（备选）
```

---

如有问题，请联系：s1214844197@163.com
