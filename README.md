# Auratus GWAS · 鲫鱼关联分析工作台

面向鲫鱼 (*Carassius auratus*) GWAS 下游分析的 Shiny 应用，按分析逻辑分为四步：

1. **群体表型比较**：ANOVA / Welch / Kruskal–Wallis，Tukey / Dunn 多重比较，异速校正
2. **LocusZoom 区域定位**：−log₁₀P、LD r²、基因结构，以及按距离排序的候选基因表
3. **基因型–表型关联**：合并相位杂合子，计算 MAF、HWE 精确检验，用加性模型估计 β 和偏 R²
4. **候选基因表达**：log₂ 尺度检验、跨基因 FDR 校正、Z-score 热图

在线版本：https://crafor18.shinyapps.io/Z_Y_ketang_data_analysis_app/

## 本地运行

```r
source("scripts/install_deps.R")
shiny::runApp(".", port = 3838)
```

需要 R ≥ 4.3。

## 目录结构

```text
app.R               页面导航与服务端入口
R/                  runApp() 自动加载
  00_setup.R        依赖、字体、示例文件路径
  utils_stats.R     基因型编码、HWE、加性模型、两两比较、表达检验
  utils_io.R        文件读取、列名识别、数据源与导出模块
  utils_theme.R     bslib 主题、配色、ggplot 出版主题
  utils_ui.R        页头、统计卡片、表格等 UI 组件
  mod_*.R           各页面模块（home / phenotype / locuszoom / genopheno / expression）
www/                样式 (styles.css) 与动效 (app.js)
docs/guide.md       应用内「指南」页
data/examples/      真实示例数据
tests/test-stats.R  统计函数测试：Rscript tests/test-stats.R
```

## 部署

凭据从环境变量读取，不写入代码：

```bash
export SHINYAPPS_TOKEN=...  SHINYAPPS_SECRET=...
Rscript rsconnect.R
```
