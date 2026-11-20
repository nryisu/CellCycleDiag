# CellCycleDiag 0.1.0

正式 R 包结构，六个公开函数。原稿、图和原分析脚本均保留。

## 安装

Windows 可直接安装旁边提供的 `CellCycleDiag_0.1.0.zip` 二进制包：

```r
install.packages("CellCycleDiag_0.1.0.zip", repos = NULL, type = "win.binary")
library(CellCycleDiag)
```

跨平台源码包是 `CellCycleDiag_0.1.0.tar.gz`：

```r
install.packages("CellCycleDiag_0.1.0.tar.gz", repos = NULL, type = "source")
```

先将 R 工作目录设为安装文件所在目录，或将文件名换为完整路径。
`CellCycleDiag_source.zip` 是供浏览、修改和上传 GitHub 的源码目录压缩包，
不是 Windows 安装包。

已有 S.Score、G2M.Score、theta 时，诊断和绘图不需要安装 Seurat 或 tricycle。
完整工作流需要另外安装：

```r
install.packages(c("Seurat", "BiocManager"))
BiocManager::install(c("tricycle", "MAST", "SingleCellExperiment",
                       "SummarizedExperiment"))
```

## 使用

```r
# seu 应已完成 QC、表达归一化；Seurat v5 的分层 data 应先合并。
# 此步骤默认仅对 human SYMBOL 使用 Seurat 的标准周期基因集。
seu <- add_cc_representation(seu, species = "human")

res <- diagnose_cc_imbalance(
  seu, group.by = "condition",
  ident.1 = "Control", ident.2 = "Treatment"
)
res
plot_cc_diagnostic(res)

sens <- run_cc_sensitivity(
  seu, group.by = "condition",
  ident.1 = "Control", ident.2 = "Treatment", donor.var = "donor"
)
sens$susceptibility
plot_cc_diagnostic(sens, type = "susceptibility")
saveRDS(sens, "cc_sensitivity.rds")
```

已有表示时直接从 `diagnose_cc_imbalance()` 开始即可，也可以传入 metadata data.frame。
mouse 或 ENSEMBL 输入须显式提供匹配 assay 行名的 `s.features`、`g2m.features`；
软件不会通过大小写转换猜测同源基因。替换已有表示需要 `overwrite = TRUE`。

## 与论文一致的定义

- D_score：S 和 G2M 两个经验 KS 距离的平均值。
- D_theta：共同原点、默认 12 箱的角度分布总变差距离；保留原脚本的右闭分箱约定。
- 三种 DE 均使用 MAST；Direct 包含指定的基线协变量；另两种模型分别增加
  S/G2M 分数或 sin(theta)/cos(theta)。提供 donor 时三种模型都控制 donor。
- 默认 DEG 为 Bonferroni 调整 P < 0.05、绝对 log2FC > 0.25、至少一组检测率 ≥ 0.1。
  梯度实验可设置 `logfc_cutoff = 0.10`。
- Jaccard distance、移除比例和新增基因数用于描述敏感性，不代表校正正确率。
- 两集合均为空时 Jaccard 为 NA；Direct 为空时 removed_fraction 为 NA。
- 缺失值、共线设计和模型失败会报错，不会被当作零 DEG。

这里的校正发生在 DE 模型中，不改写表达矩阵。Seurat 输出的 fold change
仍为边际 fold change；donor 固定效应不等于 pseudobulk 或独立生物学重复。
包不输出风险等级、推荐阈值或“是否应该校正”的结论。

## 验证和发布状态

本机已完成源码构建、安装、加载、帮助页、示例与核心测试检查；
100 个已保存的 mESC/hESC v2 抽样比较中，新旧指标实现数值一致。
本机 R 缺少 Seurat/MAST/tricycle，完整集成测试已编写但尚未运行。
详细范围见 VALIDATION.md。GitHub CI 配置已附带，尚未触发。

发布前必须在 DESCRIPTION 中将 `maintainer@example.org` 替换为真实维护者邮箱，
确认作者信息，并在有完整依赖的环境运行集成测试。当前没有创建或发布 GitHub 仓库。
