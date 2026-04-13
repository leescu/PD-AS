
saveRDS(sce_sub,file = "sce_sub.rds")
saveRDS(smc_sub,file = "smc_sub.rds")
anno <- read.csv("K:/WPY/cellranger/new_celltype.metadata.csv", stringsAsFactors = FALSE)
setwd("K:/WPY/cellranger")
data.A <- Read10X(data.dir = "A")
sce.A <- CreateSeuratObject(
  counts = data.A,
  project = "SampleA",
  min.cells = 3,
  min.features = 200
)

# B 样本
data.B <- Read10X(data.dir = "B")
sce.B <- CreateSeuratObject(
  counts = data.B,
  project = "SampleB",
  min.cells = 3,
  min.features = 200
)
colnames(sce.B) <- sub("-1$", "-2", colnames(sce.B))
sce <- merge(
  sce.A,
  y = sce.B,
  project = "A_B_scRNA"
)


sce
keep_cells <- intersect(colnames(sce), anno$Barcode)
length(keep_cells)
sce_sub <- subset(sce, cells = keep_cells)
meta0 <- sce_sub@meta.data
meta0$Barcode <- rownames(meta0)

anno2 <- as.data.frame(anno)
rownames(anno2) <- anno2$Barcode

meta_new <- merge(
  meta0,
  anno2,
  by = "Barcode",
  all.x = TRUE,
  sort = FALSE
)

rownames(meta_new) <- meta_new$Barcode
meta_new <- meta_new[colnames(sce_sub), , drop = FALSE]

sce_sub@meta.data <- meta_new
head(sce_sub@meta.data)
table(sce_sub$sampleid)
table(sce_sub$new_celltype)
sce_sub <- NormalizeData(sce_sub)
sce_sub <- FindVariableFeatures(sce_sub)
sce_sub <- ScaleData(sce_sub)
sce_sub <- RunPCA(sce_sub)
sce_sub <- FindNeighbors(sce_sub, dims = 1:20)
sce_sub <- FindClusters(sce_sub, resolution = 0.5)
sce_sub <- RunUMAP(sce_sub, dims = 1:20)
colors_0 <- c("#6A8EC9","#CC5B45") 
colors_1 <- c("#E64B35FF","#00A087FF","#3C5488FF","#7E6148FF","#D8D9DA") 
colors_2 <- c( "#CC5B45", "#EB7E60", "#FAA09C", "#F9B29C", "#FFCFD1", "#F5A216",
               "#EDC66A", "#FBDF9D", "#F9C89B", "#FBE3CD", "#B8FABF", "#98F4E0",
               "#B6E2DC", "#C6DCB9", "#3C5488", "#8FDBF3", "#8FB4DC", "#5CB0C3",
               "#6A8EC9", "#C6C3E1", "#EEC2E5", "#B46DA9", "#652884", "#7E6148",
               "#A6D854", "#FFD92F", "#E5C494", "#66C2A5", "#FC8D62", "#8DA0CB",
               "#E78AC3", "#A6D854", "#FFD92F", "#B3B3B3")
colors_3 <- c( "#CC5B45", "#FAA09C", "#FFCFD1", "#F5A216", "#B8FABF", "#3C5488",
               "#8FDBF3", "#5CB0C3", "#6A8EC9", "#B46DA9", "#652884", "#E78AC3",
               "#A6D854", "#FFD92F", "#B3B3B3")
p1<-DimPlot(
  sce_sub,
  reduction = "umap",
  group.by = "group",
  cols = colors_0, pt.size = 0.6
)
p1
ggsave("K:/WPY/Fig 1a.pdf", width = 15, height = 10, units = "cm")
p2<-DimPlot(
  sce_sub,
  reduction = "umap",group.by = "new_celltype",
  label = F,
  repel = TRUE,
  cols = colors_3,  pt.size = 0.25
)
p2
ggsave("K:/WPY/Fig 1b.pdf", width = 17, height = 10, units = "cm")
library(dplyr)
library(ggplot2)

# 统计各 group 内各细胞类型占比
df_plot <- sce_sub@meta.data %>%
  as.data.frame() %>%
  group_by(group, new_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(Frequency = n / sum(n) * 100)

# 设定横坐标顺序，可按总体丰度从高到低排序
celltype_order <- df_plot %>%
  group_by(new_celltype) %>%
  summarise(total_freq = sum(Frequency), .groups = "drop") %>%
  arrange(desc(total_freq)) %>%
  pull(new_celltype)

df_plot$new_celltype <- factor(df_plot$new_celltype, levels = celltype_order)

# 颜色
group_cols_fill <- c("Control" = "#CC5B45", "Periodontitis" = "#6A8EC9")
group_cols_line <- c("Control" = "#CC5B45", "Periodontitis" = "#6A8EC9")

# 作图
p <- ggplot(df_plot, aes(x = new_celltype, y = Frequency, fill = group, color = group)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.55,
    size = 1
  ) +
  geom_point(
    position = position_dodge(width = 0.75),
    size = 3
  ) +
  scale_fill_manual(values = group_cols_fill) +
  scale_color_manual(values = group_cols_line) +
  labs(x = "new_celltype", y = "Frequency") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )

p



# 1. 定义颜色映射（确保颜色数量与细胞类型匹配）
# 假设你的 celltype 顺序与 colors_3 一一对应
celltypes <- levels(df_plot$new_celltype)
names(colors_3) <- celltypes[1:length(colors_3)]

# 2. 准备绘图数据：计算每个扇区的起始和终止角度
df_nested <- df_plot %>%
  group_by(group) %>%
  arrange(desc(new_celltype)) %>% # 排序确保层叠顺序一致
  mutate(
    # 计算比例的累加值，用于确定扇区角度
    ymax = cumsum(Frequency),
    ymin = lag(ymax, default = 0),
    # 设置内外圈的半径位置
    # 内圈 (Control): 半径 2 到 3
    # 外圈 (Periodontitis): 半径 3.1 到 4.1
    rmin = ifelse(group == "Control", 2, 3.1),
    rmax = ifelse(group == "Control", 3, 4.1)
  ) %>%
  ungroup()


# 1. 明确定义细胞类型与颜色的对应关系
# 按照你上传的图例顺序手动列出
celltype_names <- c(
  "B_cell", "cDC", "cDC_pro", "Endothelial_cell", "Fibroblast", 
  "Mono_macrophage", "Neutrophil", "pDC", "Pericyte", "Plasma_cell", 
  "Platelets", "Schwann_cells", "Smooth_muscle_cells", "T_NK"
)

# 注意：如果你的数据里用的是下划线（如 B_cell），请把上面的空格改成下划线
# 确保 colors_3 的长度与 celltype_names 一致
names(colors_3) <- celltype_names

# 2. 强制转换数据中的因子顺序，使其与图例/颜色向量完全一致
df_nested$new_celltype <- factor(df_nested$new_celltype, levels = celltype_names)


# 2. 绘图：增加边框线层
p_final <- ggplot(df_nested) +
  # --- 核心数据层 ---
  geom_rect(aes(
    xmin = rmin, xmax = rmax, 
    ymin = ymin, ymax = ymax, 
    fill = new_celltype
  ), color = "white", size = 0.2) +
  
  # --- 增加两条外围装饰线 ---
  # 内圈 (Control) 的外边线
  annotate("rect", xmin = 3.0, xmax = 3.08, ymin = 0, ymax = 100, 
           fill = "#6A8EC9", alpha = 1) +
  
  # 外圈 (Periodontitis) 的外边线
  annotate("rect", xmin = 4.1, xmax = 4.18, ymin = 0, ymax = 100, 
           fill = "#CC5B45", alpha = 1) +
  
  # --- 坐标与样式 ---
  coord_polar(theta = "y") +
  scale_fill_manual(values = colors_3) + 
  xlim(0, 4.5) +
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right"
  ) +
  labs(subtitle = "Inner: Control (#CC5B45) | Outer: Periodontitis (#6A8EC9)")

# 显示图形
print(p_final)

ggsave("K:/WPY/Fig 1c.pdf", width = 17, height = 10, units = "cm")

## =========================
## 0) 参数区（按需改）
## =========================
kmer_min <- 1   # 保留 kmer > 1
nn_min   <- 3   # 保留 nn > 3（同一taxid出现的细胞数）

base_dir <- "L:/scRNA_microbe/SAHMI/functions/old"

samples <- list(
  A1 = list(
    report = file.path(base_dir, "A1", "A1.kraken.report.txt"),
    sckmer  = file.path(base_dir, "A1", "A1.sckmer.txt"),
    AB = "A"
  ),
  A2 = list(
    report = file.path(base_dir, "A2", "A2.kraken.report.txt"),
    sckmer  = file.path(base_dir, "A2", "A2.sckmer.txt"),
    AB = "A"
  ),
  B1 = list(
    report = file.path(base_dir, "B1", "B1.kraken.report.txt"),
    sckmer  = file.path(base_dir, "B1", "B1.sckmer.txt"),
    AB = "B"
  ),
  B2 = list(
    report = file.path(base_dir, "B2", "B2.kraken.report.txt"),
    sckmer  = file.path(base_dir, "B2", "B2.sckmer.txt"),
    AB = "B"
  )
)

## =========================
## 1) 读取 sckmer（兼容 header 被读成一个字段的情况）
## =========================
read_sckmer_robust <- function(path) {
  df <- read.table(path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  
  # 你的情况：names(df) 可能是 "barcode taxid kmer uniq" 一个字符串
  if (length(names(df)) == 1 && grepl("\\s+", names(df)[1])) {
    nm <- strsplit(trimws(names(df)[1]), "\\s+")[[1]]
    if (ncol(df) >= 4) {
      names(df)[1:4] <- nm[1:4]
    } else {
      stop("sckmer 文件未正确分列：", path)
    }
  }
  
  names(df) <- tolower(names(df))
  req <- c("barcode", "taxid", "kmer", "uniq")
  if (!all(req %in% names(df))) {
    stop("sckmer 列名不符合预期。当前列名：",
         paste(names(df), collapse = ", "),
         "\n文件：", path)
  }
  
  df$taxid <- as.integer(df$taxid)
  df$kmer  <- as.integer(df$kmer)
  df$uniq  <- as.integer(df$uniq)
  df
}

## =========================
## 2) 读取单个 lane 的 SAHMI，并 lane 内每个 barcode 取最强 taxon
## =========================
read_sahmi_lane <- function(report_file, sckmer_file, lane_id, AB,
                            kmer_min = 1, nn_min = 3) {
  
  report <- read.delim(report_file, header = FALSE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
  report$V8 <- trimws(report$V8)      # name
  report$V7 <- as.integer(report$V7)  # taxid
  
  kmer_data <- read_sckmer_robust(sckmer_file)
  
  x <- kmer_data %>%
    filter(kmer > kmer_min) %>%
    group_by(taxid) %>%
    mutate(nn = n()) %>%
    ungroup() %>%
    filter(nn > nn_min)
  
  x$name <- report$V8[match(x$taxid, report$V7)]
  x$barcode_1 <- paste0(x$barcode, "-1")
  x$lane <- lane_id
  x$AB <- AB
  
  # lane 内对同一 barcode 去重（取最强）
  meta_lane <- x %>%
    group_by(barcode_1) %>%
    arrange(desc(kmer), desc(uniq), desc(nn)) %>%
    slice(1) %>%
    ungroup()
  
  # 映射成 Seurat cell name 体系：A_ / B_
  meta_lane$cell_id <- paste0(AB, "_", meta_lane$barcode_1)
  
  as.data.frame(meta_lane)
}

## =========================
## 3) 读取四个 lane，合并，然后在 A / B 层面再次去重（关键）
## =========================
meta_list <- lapply(names(samples), function(lane_id) {
  message("Reading SAHMI lane: ", lane_id)
  read_sahmi_lane(
    report_file = samples[[lane_id]]$report,
    sckmer_file = samples[[lane_id]]$sckmer,
    lane_id     = lane_id,
    AB          = samples[[lane_id]]$AB,
    kmer_min    = kmer_min,
    nn_min      = nn_min
  )
})
meta_all <- do.call(rbind, meta_list)

# 在 A/B 层面，同一 cell_id（即 A_barcode_1）可能来自 A1 或 A2 -> 再去重一次
meta_all_clean <- meta_all %>%
  group_by(cell_id) %>%
  arrange(desc(kmer), desc(uniq), desc(nn)) %>%
  slice(1) %>%
  ungroup() %>%
  as.data.frame()

rownames(meta_all_clean) <- meta_all_clean$cell_id

# 检查是否仍重复（理论上不应重复）
stopifnot(!any(duplicated(rownames(meta_all_clean))))
meta_all_clean$Barcode <- ifelse(
  meta_all_clean$AB == "A",
     paste0(meta_all_clean$barcode, "-1"),
     paste0(meta_all_clean$barcode, "-2"))
rownames(meta_all_clean) <- meta_all_clean$Barcode
## =========================
common_cells <- intersect(colnames(sce_sub), rownames(meta_all_clean))
common_cells
message("Matched cells = ", length(common_cells))

sce_sub$microbe <- NA_character_
sce_sub@meta.data[common_cells, "microbe"]       <- meta_all_clean[common_cells, "name"]
sce_sub@meta.data[common_cells, "microbe_taxid"] <- meta_all_clean[common_cells, "taxid"]
sce_sub@meta.data[common_cells, "microbe_kmer"]  <- meta_all_clean[common_cells, "kmer"]
sce_sub@meta.data[common_cells, "microbe_uniq"]  <- meta_all_clean[common_cells, "uniq"]
sce_sub@meta.data[common_cells, "microbe_nn"]    <- meta_all_clean[common_cells, "nn"]
sce_sub@meta.data[common_cells, "microbe_lane"]  <- meta_all_clean[common_cells, "lane"]

sce_sub$microbe_pos <- ifelse(is.na(sce_sub$microbe), "No_microbe", "Microbe")

# 你的 Seurat cellname 已经是 A_ / B_，直接从 colnames 推断 AB
sce_sub$AB <- ifelse(grepl("^A_", colnames(sce_sub)), "A",
                     ifelse(grepl("^B_", colnames(sce_sub)), "B", NA))

## =========================
## 5) 结果检查 + 导出
## =========================
print(table(sce_sub$AB, useNA = "ifany"))
print(table(sce_sub$microbe_pos, useNA = "ifany"))

DimPlot(sce_sub, group.by = "microbe_pos")
DimPlot(sce_sub, group.by = "microbe_pos", split.by = "AB")

# 导出 microbe 计数（全体）
microbe_counts_all <- as.data.frame(table(sce_sub$microbe, useNA = "ifany"))
colnames(microbe_counts_all) <- c("microbe", "n")
write.csv(microbe_counts_all, file = "microbe_counts_all.csv", row.names = FALSE)

# 导出按 AB 分组计数
microbe_counts_byAB <- sce_sub@meta.data %>%
  mutate(microbe2 = ifelse(is.na(microbe), "NA", microbe)) %>%
  count(AB, microbe2, name = "n") %>%
  arrange(AB, desc(n))
write.csv(microbe_counts_byAB, file = "microbe_counts_byAB.csv", row.names = FALSE)

# 额外：看 A1/A2 合并时被“竞争掉”的数量（可选）
dup_before <- sum(duplicated(meta_all$cell_id))
message("Duplicates across lanes before A/B-level collapse = ", dup_before)
print(table(sce_sub$AB, sce_sub$microbe_pos,useNA = "ifany"))
DimPlot(
  sce_sub,
  reduction = "umap",
  group.by = "new_celltype",
  label = TRUE,
  repel = TRUE
) + NoLegend()


DimPlot(
  sce_sub,
  reduction = "umap",
  group.by = "group",
  cols = colors_0, pt.size = 1.1
)
DimPlot(
  sce_sub,
  reduction = "umap",group.by = "new_celltype",
  label = TRUE,
  repel = TRUE,
  cols = colors_3, pt.size = 1.1
)
DimPlot(sce_sub, group.by = "microbe_pos")
colors_a <- c( "#F5A216", "#652884")
DimPlot(sce_sub, group.by = "microbe_pos", cols = colors_a,pt.size = 0.3)
ggsave("K:/WPY/Fig 1d.pdf", width = 15, height = 10, units = "cm")


df_plot <- sce_sub@meta.data %>%
  as.data.frame() %>%
  group_by(microbe_pos, new_celltype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(microbe_pos) %>%
  mutate(Frequency = n / sum(n) * 100)

# 设定横坐标顺序，可按总体丰度从高到低排序
celltype_order <- df_plot %>%
  group_by(new_celltype) %>%
  summarise(total_freq = sum(Frequency), .groups = "drop") %>%
  arrange(desc(total_freq)) %>%
  pull(new_celltype)

df_plot$new_celltype <- factor(df_plot$new_celltype, levels = celltype_order)
# 颜色
group_cols_fill <- c("No_microbe" =  "#652884", "Microbe" = "#F5A216")
group_cols_line <- c("No_microbe" =  "#652884", "Microbe" ="#F5A216" )

# 作图
p <- ggplot(df_plot, aes(x = new_celltype, y = Frequency, fill =microbe_pos, color = microbe_pos)) +
  geom_bar(
    stat = "identity",
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.55,
    size = 1.5
  ) +
  geom_point(
    position = position_dodge(width = 0.75),
    size = 1.5
  ) +
  scale_fill_manual(values = group_cols_fill) +
  scale_color_manual(values = group_cols_line) +
  labs(x = "new_celltype", y = "Frequency") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 12)
  )

p
ggsave("K:/WPY/Fig 1e.pdf", width = 20, height = 10, units = "cm")


df_nested_final <- df_plot %>%
  group_by(microbe_pos) %>%
  # 按照细胞类型一致排序，确保内外圈扇区对齐（可选）
  arrange(microbe_pos, desc(new_celltype)) %>% 
  mutate(
    # 计算角度方向的累积位置 (0-100)
    ymax = cumsum(Frequency),
    ymin = lag(ymax, default = 0),
    # 分配半径区间：
    # 假设 Microbe 是内圈，其他（如 Non-Microbe）是外圈
    # 如果你想反过来，调换下方的 2/3 和 3.1/4.1
    rmin = ifelse(microbe_pos == "Microbe", 2, 3.1),
    rmax = ifelse(microbe_pos == "Microbe", 3, 4.1)
  ) %>%
  ungroup()

# 3. 绘图
p_nested <- ggplot(df_nested_final) +
  # --- 核心细胞类型分布层 ---
  geom_rect(aes(
    xmin = rmin, xmax = rmax, 
    ymin = ymin, ymax = ymax, 
    fill = new_celltype
  ), color = "white", size = 0.2) +
  
  # --- 增加两条组别外围线 ---
  # 内圈 Microbe 的外边线 (#CC5B45)
  annotate("rect", xmin = 3.0, xmax = 3.05, ymin = 0, ymax = 100, 
           fill = "#F5A216", alpha = 0.9) +
  
  # 外圈 (假设为另一组) 的外边线 (#6A8EC9)
  annotate("rect", xmin = 4.1, xmax = 4.15, ymin = 0, ymax = 100, 
           fill = "#652884", alpha = 0.9) +
  
  # --- 坐标系与颜色 ---
  coord_polar(theta = "y") +
  scale_fill_manual(values = colors_3) + 
  xlim(0, 4.5) + # 0控制中心孔洞大小
  theme_void() +
  theme(
    legend.title = element_blank(),
    legend.position = "right",
    plot.margin = unit(c(0,0,0,0), "cm")
  ) +
  labs(title = "Cell Type Distribution by Microbe Presence")
p_nested 
ggsave("K:/WPY/Fig 1f.pdf", width = 17, height = 10, units = "cm")


smc_sub <- subset(sce_sub, subset = new_celltype == "Smooth_muscle_cells")

# 2. 重新进行标准化和聚类（这一步至关重要，消除大群背景干扰）
smc_sub <- FindVariableFeatures(smc_sub, nfeatures = 2000)
smc_sub <- ScaleData(smc_sub)
smc_sub <- RunPCA(smc_sub)
smc_sub <- RunUMAP(smc_sub, dims = 1:15) # 维度可根据 ElbowPlot 调整
smc_sub <- FindNeighbors(smc_sub, dims = 1:15)

# 3. 尝试不同的分辨率，寻找潜在的亚群（如收缩型 vs. 转化型）
smc_sub <- FindClusters(smc_sub, resolution = 1.5) 
smc_sub <- FindClusters(smc_sub, resolution = 0.1) 


table(smc_sub$group)
Idents(smc_sub) <- "group"

# 2. 寻找 Cluster 0 的差异基因 (相对于其他所有 cluster)

group_diff_markers <- FindMarkers(smc_sub, 
                                  ident.1 = "Periodontitis", 
                                  ident.2 = "Control",
                                  logfc.threshold = 0.25, # 至少 1.2 倍差异
                                  min.pct = 0.1)          # 至少 10% 的细胞表达

# 3. 查看上调最显著的前 30 个基因 (即在牙周炎中升高的)
print(head(group_diff_markers %>% arrange(desc(avg_log2FC)), 30))

# 4. 查看下调最显著的前 30 个基因 (即在牙周炎中降低的)
print(head(group_diff_markers %>% arrange(avg_log2FC), 30))


library(ggplot2)
library(ggrepel)

# 标记上调、下调基因
group_diff_markers$gene <- rownames(group_diff_markers)
group_diff_markers$change <- ifelse(group_diff_markers$p_val_adj < 0.05 & abs(group_diff_markers$avg_log2FC) > 0.25, 
                                    ifelse(group_diff_markers$avg_log2FC > 0.25, "UP", "DOWN"), "NOT")

# 绘图
library(ggplot2)
library(ggrepel)
library(dplyr)
write.csv(plot_df, "K:/WPY/SMC_DEG_all.csv", row.names = FALSE)

# 处理 P 值极值
plot_df <- group_diff_markers
plot_df$logP <- -log10(plot_df$p_val_adj + 1e-305) # 防止 log(0)
plot_df$logP[plot_df$logP > 300] <- 300 # 统一截断在 300

# 筛选真正有意义的标注基因 (Log2FC > 1 且 logP 很大)
label_genes <- plot_df %>% 
  filter(abs(avg_log2FC) > 1 & logP > 50) %>% 
  arrange(desc(abs(avg_log2FC))) %>% 
  head(10) # 只选前 20 个

ggplot(plot_df, aes(x = avg_log2FC, y = logP, color = change)) +
  geom_point(alpha = 1, size = 1.5) +
  scale_color_manual(values = c("UP" = "#CC5B45", "DOWN" = "#6A8EC9", "NOT" = "grey80")) +
  # 关键：增加 max.overlaps 和标注设置
  geom_text_repel(data = label_genes, 
                  aes(label = gene),
                  size = 3.5, 
                  fontface = "bold",
                  box.padding = 0.5, 
                  point.padding = 0.3,
                  segment.color = "black",
                  max.overlaps = 50) + # 允许更多重叠尝试
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  theme_classic() +
  labs(x = "log2 (Fold Change)", 
       y = "-log10 (Adjusted P-value) [Capped at 300]")
ggsave("K:/WPY/Fig 1g.pdf", width = 3.8, height = 3, units = "in")


library(openxlsx)
library(clusterProfiler)
library(org.Mm.eg.db) # 小鼠数据库

# 1. 准备基因列表 (假设 genes_to_test 是你的差异基因 Symbol)
gene_convert <- bitr(genes_to_test, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)

# 2. 定义跑 GO 并格式化的函数
run_go_all <- function(ont) {
  res <- enrichGO(gene          = gene_convert$ENTREZID,
                  OrgDb         = org.Mm.eg.db,
                  ont           = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  readable      = TRUE)
  return(as.data.frame(res))
}

# 3. 运行 BP, CC, MF
go_bp <- run_go_all("BP")
go_cc <- run_go_all("CC")
go_mf <- run_go_all("MF")

# 4. 保存到 Excel 多个 Sheet
write.xlsx(list("GO_BP" = go_bp, "GO_CC" = go_cc, "GO_MF" = go_mf), 
           file = "K:/WPY/SMC_GO_Results.xlsx")

message("GO 结果已保存至 Excel!")
kegg_res <- enrichKEGG(gene         = gene_convert$ENTREZID,
                       organism     = 'mmu', # 如果是人则用 'hsa'
                       pvalueCutoff = 0.05)

# 2. 转换为可读的基因名 (Symbol)
kegg_res_readable <- setReadable(kegg_res, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

# 3. 绘制 KEGG 气泡图
p_kegg <- dotplot(kegg_res_readable, showCategory = 15) + 
  ggtitle("KEGG Pathway Enrichment") +
  theme_bw()

print(p_kegg)

# 4. 保存 KEGG 结果到之前的 Excel 中 (追加 Sheet)
wb <- loadWorkbook("K:/WPY/SMC_Cluster0_GO_Results.xlsx")
addWorksheet(wb, "KEGG_Pathway")
writeData(wb, "KEGG_Pathway", as.data.frame(kegg_res_readable))
saveWorkbook(wb, "K:/WPY/SMC_KEGG_Results.xlsx", overwrite = TRUE)


library(ggplot2)
library(dplyr)
library(stringr)
# 1. 手动筛选与 AS 最相关的 5 条 GO 通路 ID
go_as_ids <- c("GO:0048511", "GO:0007623", "GO:0045444", "GO:1901342", "GO:0030336")

# 2. 从 go_res 对象中提取数据并筛选
# 假设 go_res 是你运行 enrichGO 的结果对象
plot_go_data <- as.data.frame(go_res) %>%
  filter(ID %in% go_as_ids) %>%
  mutate(RichFactor = Count / as.numeric(sub(".*/", "", BgRatio)))  %>%
  # 关键步骤：将描述转换为每个单词首字母大写
  mutate(Description = str_to_title(Description))
# 3. 绘图
p_go_as <- ggplot(plot_go_data, aes(x = RichFactor, y = reorder(Description, RichFactor))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#A50F15", high = "#FEE5D9") + # 红色系
  labs(title = "Top 5 Atherosclerosis-related GO Terms",
       x = "Rich Factor", y = NULL, size = "Gene Count", color = "P-adjust") +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
    axis.text = element_text(size = 10, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_go_as)
ggsave("K:/WPY/Fig SA.pdf", width = 3.8, height = 3, units = "in")

# 1. 手动筛选与 AS 最相关的 5 条 KEGG 通路 ID (小鼠为 mmu)
# 1. 定义你想要的顺序（从上到下）
kegg_as_ids <- c("mmu04010", "mmu04710","mmu04310","mmu04611","mmu05418")

# 2. 筛选与处理数据
plot_kegg_data <- as.data.frame(kegg_res) %>%
  filter(ID %in% kegg_as_ids) %>%
  mutate(Description = case_when(
    ID == "mmu05418" ~ "Fluid shear stress and AS",
    TRUE ~ Description
  )) %>%
  mutate(Description = str_to_title(Description))%>%
  # 使用 parse_ratio 函数或手动计算
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    return(parts[1] / parts[2])
  }))

# --- 关键修改：锁定顺序 ---
# 我们根据 kegg_as_ids 的顺序来匹配 Description
# 使用 rev() 是因为 ggplot 坐标轴 0 点在下方，rev 确保第一个 ID 在最上面
target_order <- plot_kegg_data$Description[match(kegg_as_ids, plot_kegg_data$ID)]
plot_kegg_data$Description <- factor(plot_kegg_data$Description, levels = rev(target_order))

# 3. 绘图
p_kegg_as <- ggplot(plot_kegg_data, aes(x = GeneRatio, y = Description)) + # 去掉 reorder
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_size_continuous(range = c(5, 9)) +
  # 注意：不能同时使用 scale_color_manual 和 gradient，这里保留渐变
  scale_color_gradient(low = "#B22222", high = "#F8766D") +
  labs(title = "Top 5 Atherosclerosis-related KEGG Pathways",
       x = "Gene Ratio", y = NULL) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
    axis.text = element_text(size = 10, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_kegg_as)
ggsave("K:/WPY/Fig 1h.pdf", width = 7, height = 3.5, units = "in")




library(dplyr)
library(ggplot2)
library(scales)
microbe_data <- data.frame(Cell_ID = rownames(sce_sub@meta.data), 
                           microbe = sce_sub@meta.data[["microbe"]])

# 2. 保存为 CSV 文件
# 请根据你的实际路径修改 path
write.csv(microbe_data, 
          file = "K:/WPY/sce_sub_microbe_metadata.csv", 
          row.names = FALSE)


library(Seurat)
library(stringr)
library(dplyr)

# 1. 获取 meta.data 中 microbe 列的所有唯一值
all_microbes <- unique(as.character(sce_sub@meta.data$microbe))

# 2. 定义清洗逻辑
# 逻辑：
# - 属 (Genus): 取第一个单词，并去掉方括号
# - 种 (Species): 如果是标准双名法 (如 Klebsiella pneumoniae)，保留前两个词；
#                如果是 sp. 或只有属名，则种名记为 NA。
clean_df <- data.frame(raw_name = all_microbes) %>%
  mutate(
    # 处理属名：去掉方括号 [ ]，取第一个单词
    microbe_genus = str_remove_all(raw_name, "\\[|\\]") %>% word(1),
    
    # 处理种名：
    microbe_species = sapply(raw_name, function(x) {
      x_clean <- str_remove_all(x, "\\[|\\]")
      words <- str_split(x_clean, " ")[[1]]
      
      # 如果有两个及以上单词，且第二个单词不是 sp. 或类似占位符
      if (length(words) >= 2 && !str_detect(words[2], "(?i)sp|uncultured|sp\\.|R[0-9]+")) {
        return(paste(words[1], words[2]))
      } else {
        return(NA) # 无法确定具体种名
      }
    })
  )

# 3. 将清洗结果映射回 sce_sub 对象
# 使用 match 函数将 raw_name 对应的清洗结果填入 metadata
sce_sub@meta.data$microbe_genus <- clean_df$microbe_genus[match(sce_sub@meta.data$microbe, clean_df$raw_name)]
sce_sub@meta.data$microbe_species <- clean_df$microbe_species[match(sce_sub@meta.data$microbe, clean_df$raw_name)]

# 4. 验证清洗效果
print("前10行清洗结果预览：")
print(head(sce_sub@meta.data[, c("microbe", "microbe_genus", "microbe_species")], 10))

# 5. 保存结果
# 保存带清洗后信息的 metadata 为 CSV 方便复核
write.csv(sce_sub@meta.data, "K:/WPY/sce_sub_microbe_cleaned_v2.csv", row.names = TRUE)



# 1. 提取 Smooth_muscle_cells，去掉 NA
df_microbe <- sce_sub@meta.data %>%
  as.data.frame() %>%
  filter(new_celltype == "Smooth_muscle_cells") %>%
  filter(!is.na(microbe_genus), microbe_genus != "")

# 2. 取数量前15的菌
top15_microbes <- df_microbe %>%
  count(microbe_genus, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(microbe_genus)

# 3. 其余归为 Other
df_microbe <- df_microbe %>%
  mutate(microbe_plot = ifelse(microbe_genus %in% top15_microbes, microbe_genus, "Other"))

# 4. 计算每个 group 内各菌占比
df_plot <- df_microbe %>%
  group_by(group, microbe_plot) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(freq = n / sum(n))

# 5. 设置显示顺序：前15按总数降序，Other放最后
microbe_order <- df_microbe %>%
  mutate(microbe_plot = ifelse(microbe_genus %in% top15_microbes, microbe_genus, "Other")) %>%
  count(microbe_plot, sort = TRUE) %>%
  pull(microbe_plot)

microbe_order <- c(setdiff(microbe_order, "Other"), "Other")
df_plot$microbe_plot <- factor(df_plot$microbe_plot, levels = microbe_order)

# 6. 配色
n_microbe <- length(microbe_order)
microbe_cols <- setNames(c(colors_2[1:(n_microbe - 1)], "#B3B3B3"), microbe_order)

# 7. 作图
p <- ggplot(df_plot, aes(x = group, y = freq, fill = microbe_plot)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = microbe_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Group",
    y = "Proportion",
    fill = "Microbe"  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  
  )
p
ggsave("K:/WPY/Fig 1i.pdf", width = 6, height = 5, units = "in")

library(dplyr)
library(ggplot2)
library(scales)

# 1. 提取非NA的微生物信息
df_microbe <- sce_sub@meta.data %>%
  as.data.frame() %>%
  filter(!is.na(microbe_genus), microbe_genus != "")

# 2. 选总数前15的菌
top15_microbes <- df_microbe %>%
  count(microbe_genus, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(microbe_genus)

# 3. 其余归为 Other
df_microbe <- df_microbe %>%
  mutate(microbe_plot = ifelse(microbe_genus %in% top15_microbes, microbe_genus, "Other"))

# 4. 计算每个细胞类型内各种菌的比例
df_plot <- df_microbe %>%
  group_by(new_celltype, microbe_plot) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(new_celltype) %>%
  mutate(freq = n / sum(n))

# 5. 细胞类型顺序：按含菌细胞总数排序
celltype_order <- df_microbe %>%
  count(new_celltype, sort = TRUE) %>%
  pull(new_celltype)

df_plot$new_celltype <- factor(df_plot$new_celltype, levels = celltype_order)

# 6. 菌的顺序：前15按总数降序，Other放最后
microbe_order <- df_microbe %>%
  count(microbe_plot, sort = TRUE) %>%
  pull(microbe_plot)

microbe_order <- c(setdiff(microbe_order, "Other"), "Other")
df_plot$microbe_plot <- factor(df_plot$microbe_plot, levels = microbe_order)

# 7. 配色
n_microbe <- length(microbe_order)
microbe_cols <- setNames(c(colors_2[1:(n_microbe - 1)], "#B3B3B3"), microbe_order)

# 8. 作图
p <- ggplot(df_plot, aes(x = new_celltype, y = freq, fill = microbe_plot)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = microbe_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Cell type",
    y = "Proportion",
    fill = "Microbe"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

p
ggsave("K:/WPY/Fig 1j.pdf", width =11, height = 5, units = "in")



# 1. 提取 Smooth_muscle_cells，去掉 NA
df_microbe <- sce_sub@meta.data %>%
  as.data.frame() %>%
  filter(new_celltype == "Smooth_muscle_cells") %>%
  filter(!is.na(microbe_species), microbe_species != "")

# 2. 取数量前15的菌
top15_microbes <- df_microbe %>%
  count(microbe_species, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(microbe_species)

# 3. 其余归为 Other
df_microbe <- df_microbe %>%
  mutate(microbe_plot = ifelse(microbe_species %in% top15_microbes, microbe_species, "Other"))

# 4. 计算每个 group 内各菌占比
df_plot <- df_microbe %>%
  group_by(group, microbe_plot) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group) %>%
  mutate(freq = n / sum(n))

# 5. 设置显示顺序：前15按总数降序，Other放最后
microbe_order <- df_microbe %>%
  mutate(microbe_plot = ifelse(microbe_species %in% top15_microbes, microbe_species, "Other")) %>%
  count(microbe_plot, sort = TRUE) %>%
  pull(microbe_plot)

microbe_order <- c(setdiff(microbe_order, "Other"), "Other")
df_plot$microbe_plot <- factor(df_plot$microbe_plot, levels = microbe_order)

# 6. 配色
n_microbe <- length(microbe_order)
microbe_cols <- setNames(c(colors_2[1:(n_microbe - 1)], "#B3B3B3"), microbe_order)

# 7. 作图
p <- ggplot(df_plot, aes(x = group, y = freq, fill = microbe_plot)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = microbe_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Group",
    y = "Proportion",
    fill = "Microbe"  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
    
  )
p
ggsave("K:/WPY/Fig SB.pdf", width = 6, height = 5, units = "in")

library(dplyr)
library(ggplot2)
library(scales)

# 1. 提取非NA的微生物信息
df_microbe <- sce_sub@meta.data %>%
  as.data.frame() %>%
  filter(!is.na(microbe_species), microbe_species != "")

# 2. 选总数前15的菌
top15_microbes <- df_microbe %>%
  count(microbe_species, sort = TRUE) %>%
  slice_head(n = 20) %>%
  pull(microbe_species)

# 3. 其余归为 Other
df_microbe <- df_microbe %>%
  mutate(microbe_plot = ifelse(microbe_species %in% top15_microbes, microbe_species, "Other"))

# 4. 计算每个细胞类型内各种菌的比例
df_plot <- df_microbe %>%
  group_by(new_celltype, microbe_plot) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(new_celltype) %>%
  mutate(freq = n / sum(n))

# 5. 细胞类型顺序：按含菌细胞总数排序
celltype_order <- df_microbe %>%
  count(new_celltype, sort = TRUE) %>%
  pull(new_celltype)

df_plot$new_celltype <- factor(df_plot$new_celltype, levels = celltype_order)

# 6. 菌的顺序：前15按总数降序，Other放最后
microbe_order <- df_microbe %>%
  count(microbe_plot, sort = TRUE) %>%
  pull(microbe_plot)

microbe_order <- c(setdiff(microbe_order, "Other"), "Other")
df_plot$microbe_plot <- factor(df_plot$microbe_plot, levels = microbe_order)

# 7. 配色
n_microbe <- length(microbe_order)
microbe_cols <- setNames(c(colors_2[1:(n_microbe - 1)], "#B3B3B3"), microbe_order)

# 8. 作图
p <- ggplot(df_plot, aes(x = new_celltype, y = freq, fill = microbe_plot)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = microbe_cols) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    x = "Cell type",
    y = "Proportion",
    fill = "Microbe"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11)
  )

p
ggsave("K:/WPY/Fig SC.pdf", width =11, height = 5, units = "in")

library(Seurat)
library(ggplot2)
library(dplyr)

# --- 步骤 1: 标记阳性细胞 (只要包含其中一种菌就算阳性) ---
# 确保你已经运行了清洗 microbe 的代码，得到了 microbe_species 列

# 1. 标记红色复合体阳性细胞
sce_sub$Red_Complex_Status <- ifelse(
  sce_sub$microbe_genus %in% c("Porphyromonas", "Tannerella", "Treponema"),
  "Red+", "Negative"
)

# 2. 标记橙色复合体阳性细胞 (常见如 Fn, Pi, Pm)
sce_sub$Orange_Complex_Status <- ifelse(
  sce_sub$microbe_genus %in% c("Fusobacterium", "Prevotella", "Parvimonas","Campylobacter","Eubacterium"),
  "Orange+", "Negative"
)


# 提取所需的 metadata 和 UMAP 坐标
plot_data <- FetchData(sce_sub, vars = c("umap_1", "umap_2", "new_celltype", "Red_Complex_Status", "Orange_Complex_Status"))

# 信号：筛选出非阴性的细胞
bg_data <- plot_data
table(plot_data$Red_Complex_Status)
red_signal_data <- plot_data %>% filter(Red_Complex_Status == "Red+")
orange_signal_data <- plot_data %>% filter(Orange_Complex_Status == "Orange+")

# 例如：4 个亚群用 4 种颜色，Negative 设为浅灰
colors_3 <- c( "#CC5B45", "#FAA09C", "#FFCFD1", "#F5A216", "#B8FABF", "#3C5488",
               "#8FDBF3", "#5CB0C3", "#6A8EC9", "#B46DA9", "#652884", "#E78AC3",
               "#A6D854", "#FFD92F", "#B3B3B3")

p_complex_overlay <- ggplot() +
  
  # 图层 1 (底层): 绘制所有细胞的亚群图 (普通的 geom_point)
  geom_point(data = bg_data, 
             aes(x = umap_1, y = umap_2, color = new_celltype), 
             size = 0.8, 
             alpha = 0.6) + # 稍微透明一点，突出上层信号
  
  # 设置底层亚群颜色
  scale_color_manual(values = colors_3) +
  
  # 图层 2 (中层): 绘制橙色复合体阳性细胞 (带黑边、 shape=21)
  # 优先画橙色，防止被红色覆盖
  geom_point(data = orange_signal_data, 
             aes(x = umap_1, y = umap_2, fill = Orange_Complex_Status), 
             shape = 21,    # 关键：shape 21 是带边框的实心圆
             color = "black", # 边框颜色为黑色
             size =3.5,      # 信号点稍微做大一点，突出显示
             stroke = 0.2) +  # 边框的粗细
  
  # 图层 3 (顶层): 绘制红色复合体阳性细胞 (带黑边、 shape=21)
  geom_point(data = red_signal_data, 
             aes(x = umap_1, y = umap_2, fill = Red_Complex_Status), 
             shape = 21,    
             color = "black", 
             size = 3.5,      
             stroke = 0.2) +
  
  # 设置填充颜色 (上层红/橙底色)
  scale_fill_manual(values = c("Red+" = "#B22222", # 深红
                               "Orange+" = "#FF8C00")) + # 深橙
  
  # 主题设置
  theme_bw() + # 带框的主题
  theme(
    panel.grid = element_blank(), # 去掉网格线
    axis.text = element_text(color = "black"),
    legend.position = "right", # 图例放右边
    legend.title = element_blank(), # 去掉图例标题
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14)
  ) +
  labs(title = "Periodontal Complex Colonization across SMC Subtypes",
       x = "UMAP 1", y = "UMAP 2")
p_complex_overlay
ggsave("K:/WPY/Fig 1k.pdf", p_complex_overlay, width = 8, height = 6.5)


library(Seurat)
library(dplyr)

sce_smc <- subset(
  sce_sub,
  subset = group == "Periodontitis" & new_celltype == "Smooth_muscle_cells"
)

sce_smc$microbe_status <- ifelse(
  is.na(sce_smc$microbe) | sce_smc$microbe == "",
  "Microbe_Neg",
  "Microbe_Pos"
)

table(sce_smc$microbe_status)



Idents(sce_smc) <- sce_smc$microbe_status

deg_smc <- FindMarkers(
  sce_smc,
  ident.1 = "Microbe_Pos",
  ident.2 = "Microbe_Neg",
  logfc.threshold = 0,
  min.pct = 0.1,
  test.use = "wilcox"
)

deg_smc$gene <- rownames(deg_smc)
deg_smc <- deg_smc %>% arrange(p_val_adj, desc(avg_log2FC))

head(deg_smc)
write.csv(deg_smc, "K:/WPY/SMC_MicrobePos_vs_MicrobeNeg_DEG_all.csv", row.names = FALSE)


library(ggplot2)
deg_sig <- deg_smc %>%
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.25)

deg_up <- deg_smc %>%
  filter(p_val_adj < 0.05 & avg_log2FC > 0.25)

deg_down <- deg_smc %>%
  filter(p_val_adj < 0.05 & avg_log2FC < -0.25)

nrow(deg_sig)
nrow(deg_up)
nrow(deg_down)
deg_smc$group <- "NS"
deg_smc$group[deg_smc$p_val_adj < 0.05 & deg_smc$avg_log2FC > 0.25] <- "Up"
deg_smc$group[deg_smc$p_val_adj < 0.05 & deg_smc$avg_log2FC < -0.25] <- "Down"

ggplot(deg_smc, aes(x = avg_log2FC, y = -log10(p_val_adj), color = group)) +
  geom_point(size = 1) +
  scale_color_manual(values = c("Down" = "#6A8EC9", "NS" = "grey80", "Up" = "#CC5B45")) +
  theme_classic() +
  labs(x = "avg_log2FC", y = "-log10(adj.P)")


top10_up <- deg_up %>% slice_max(order_by = avg_log2FC, n = 10) %>% pull(gene)
top10_down <- deg_down %>% slice_min(order_by = avg_log2FC, n = 10) %>% pull(gene)
genes_show <- unique(c(top10_up, top10_down))
genes_show <- intersect(genes_show, rownames(sce_smc))

sce_smc <- ScaleData(sce_smc, features = genes_show)
table(sce_smc$microbe_status)

my_status_cols <- c("Microbe_Neg" = "#6A8EC9",
                    "Microbe_Pos" =  "#CC5B45")
p <- DoHeatmap(
  sce_smc,
  features = genes_show,
  group.by = "microbe_status",
  group.colors = my_status_cols
)
p
ggsave("K:/WPY/Fig SD.pdf", width = 3.8, height = 3, units = "in")




library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(ggplot2)
deg_smc$gene <- rownames(deg_smc)

deg_smc <- deg_smc %>%
  filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(p_val_adj, desc(avg_log2FC))

deg_up <- deg_smc %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)

deg_down <- deg_smc %>%
  filter(p_val_adj < 0.05, avg_log2FC < -0.25)

gene_up_symbol <- deg_up$gene
gene_down_symbol <- deg_down$gene

gene_up_entrez <- bitr(
  gene_up_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

gene_down_entrez <- bitr(
  gene_down_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)
ego_up <- enrichGO(
  gene = gene_up_symbol,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = gene_down_symbol,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
ekegg_up <- enrichKEGG(
  gene = gene_up_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

ekegg_down <- enrichKEGG(
  gene = gene_down_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(ego_up), "SMC_mouse_GO_BP_up.csv", row.names = FALSE)
write.csv(as.data.frame(ego_down), "SMC_mouse_GO_BP_down.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_up), "SMC_mouse_KEGG_up.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_down), "SMC_mouse_KEGG_down.csv", row.names = FALSE)
dotplot(ego_up, showCategory = 15) + ggtitle("Mouse SMC GO BP Up")
dotplot(ekegg_up, showCategory = 15) + ggtitle("Mouse SMC KEGG Up")




kegg_as_ids <- c("mmu05208", "mmu00190", "mmu04217", "mmu00480", "mmu00030")

# 2. 筛选数据
plot_kegg_data <- as.data.frame(ekegg_up) %>%
  filter(ID %in% kegg_as_ids) %>%
  # 将复杂的名称简写，方便排版
  mutate(Description = case_when(
    ID == "mmu05418" ~ "Fluid shear stress and AS",
    TRUE ~ Description
  ))%>%
  mutate(Description = str_to_title(Description))%>%
  # 使用 parse_ratio 函数或手动计算
  mutate(GeneRatio = sapply(GeneRatio, function(x) {
    parts <- as.numeric(unlist(strsplit(x, "/")))
    return(parts[1] / parts[2])
  }))
plot_kegg_data
# 3. 绘图
p_kegg_as <- ggplot(plot_kegg_data, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_size_continuous(range = c(4, 8)) +
  scale_color_manual(values = "#B22222") + # 如果只有上调，建议用统一深红或渐变
  scale_color_gradient(low = "#B22222", high = "#F8766D") +
  labs(title = "Top 5 Atherosclerosis-related KEGG Pathways",
       x = "Gene Ratio", y = NULL) +
  theme_bw() +
  theme(
    panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
    axis.text = element_text(size = 10, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

print(p_kegg_as)
ggsave("K:/WPY/Fig 1i.pdf", width = 8, height = 3.5, units = "in")


























library(ggplot2)
library(dplyr)
library(scales)
table(smc_sce_sub$group)
# 1. 定义复合体分类逻辑
# 红色复合体：Pg, Tf, Td
smc_sce_sub<-subset(sce_sub,new_celltype=="Smooth_muscle_cells"&group=="Periodontitis")
# 橙色复合体：Fn, Pi, Pm (根据你数据中的实际名称微调)
smc_sce_sub@meta.data <- smc_sce_sub@meta.data %>%
  mutate(Complex_Group = case_when(
    microbe_genus %in% c("Porphyromonas", "Tannerella", "Treponema") ~ "Red Complex",
    microbe_genus %in% c("Fusobacterium", "Prevotella", "Parvimonas","Campylobacter","Eubacterium") ~ "Orange Complex",
    is.na(microbe_genus) | microbe_genus == "None" ~ "Negative",
    TRUE ~ "Other Microbes" # 其他杂菌
  ))

library(Seurat)
library(dplyr)
library(ggplot2)

# 1. 设置分群标识


# 2. 差异分析 (推荐使用 MAST 算法，更适合处理单细胞中的稀疏数据和不平衡样本)
# 如果没有安装 MAST，可以使用默认的 wilcox
# 1. 合并所有 layers 
smc_sce_sub <- JoinLayers(smc_sce_sub)

# 2. 再次尝试运行差异分析
Idents(smc_sce_sub) <- "microbe_status"

deg_microbe <- FindMarkers(
  smc_sce_sub, 
  ident.1 = "Positive", 
  ident.2 = "Negative", 
  logfc.threshold = 0.1, 
  min.pct = 0.05, 
  only.pos = FALSE,
  assay = "RNA" # 确保指定 Assay
)

# 过滤出显著差异基因 (P_adj < 0.05)
deg_sig <- deg_microbe %>% filter(p_val_adj < 0.05)
write.csv(deg_sig, "K:/WPY/SMC_Microbe_Positive_vs_Negative_DEGs.csv")


group_diff_markers <- deg_microbe
group_diff_markers$gene <- rownames(group_diff_markers)

# 1. 标记上调、下调基因 (阈值设定为 P_adj < 0.05 且 |log2FC| > 0.25)
group_diff_markers$change <- case_when(
  group_diff_markers$p_val_adj < 0.05 & group_diff_markers$avg_log2FC > 0.25 ~ "UP",
  group_diff_markers$p_val_adj < 0.05 & group_diff_markers$avg_log2FC < -0.25 ~ "DOWN",
  TRUE ~ "NOT"
)

# 2. 处理 P 值极值与截断
plot_df <- group_diff_markers
plot_df$logP <- -log10(plot_df$p_val_adj + 1e-305) # 防止 log(0)
plot_df$logP[plot_df$logP > 300] <- 300 # 统一截断在 300

# 3. 筛选标注基因 (针对微生物正相关细胞，侧重标注上调明显的基因)
# 注意：如果基因太少，可以适当调低 abs(avg_log2FC) > 0.5
label_genes <- plot_df %>% 
  filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5) %>% 
  arrange(desc(abs(avg_log2FC))) %>% 
  head(10) # 标注前 15 个最显著的基因
# 4. 开始绘图
p_volcano_microbe <- ggplot(plot_df, aes(x = avg_log2FC, y = logP, color = change)) +
  # 背景散点层
  geom_point(alpha = 0.8, size = 1.2) +
  
  # 设置配色 (遵循你要求的红蓝配色)
  scale_color_manual(values = c("UP" = "#CC5B45", "DOWN" = "#6A8EC9", "NOT" = "grey85")) +
  
  # 关键：基因名称标注
  geom_text_repel(data = label_genes, 
                  aes(label = gene),
                  size = 3.5, 
                  fontface = "bold",
                  color = "black", # 标签颜色统一设为黑色更清晰
                  box.padding = 0.6, 
                  point.padding = 0.4,
                  segment.color = "grey30",
                  segment.size = 0.3,
                  max.overlaps = 100) + 
  
  # 辅助线
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  
  # 主题与标签
  theme_classic() +
  labs(
    title = "SMC: Microbe-Positive vs Microbe-Negative",
    x = "log2 (Fold Change)", 
    y = "-log10 (Adjusted P-value) [Capped at 300]"
  ) +
  theme(
    axis.text = element_text(size = 11, color = "black"),
    axis.title = element_text(size = 12, color = "black", face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "right"
  )

# 查看效果
print(p_volcano_microbe)
ggsave("K:/WPY/Fig SE.pdf", width = 3.8, height = 3, units = "in")

library(clusterProfiler)
library(org.Mm.eg.db)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggplot2)
library(dplyr)
library(stringr)
# 假设你的差异基因（Positive vs Negative）在 genes_to_test
gene_convert <- bitr(genes_to_test, fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Mm.eg.db)

# 运行 GO BP (通常 BP 最能反映生物学功能)
go_res <- enrichGO(gene          = gene_convert$ENTREZID,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   readable      = TRUE)
plot_go_data <- as.data.frame(go_res)

# 2. 按 p.adjust 排序（从小到大）
plot_go_data <- plot_go_data[order(plot_go_data$p.adjust), ]

# 3. 提取前 5 行（代替 slice）
if(nrow(plot_go_data) > 5) {
  plot_go_data <- plot_go_data[1:5, ]
}

# 4. 格式化描述（首字母大写）
plot_go_data$Description <- str_to_title(plot_go_data$Description)
write.csv(plot_go_data, "K:/WPY/SMC_Microbe_Positive_vs_Negative_go.csv")

# 5. 绘图
p_go_top5 <- ggplot(plot_go_data, aes(x = RichFactor, y = reorder(Description, RichFactor))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#A50F15", high = "#FEE5D9") + 
  labs(title = "Top 5 Enriched GO Terms", x = "Rich Factor", y = NULL, color = "P-adj") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
        axis.text = element_text(size = 10, color = "black", face = "bold"))

print(p_go_top5)
ggsave("K:/WPY/Fig SF.pdf", width = 3.8, height = 3, units = "in")
plot_kegg_data <- as.data.frame(kegg_res)
plot_kegg_data <- plot_kegg_data[order(plot_kegg_data$p.adjust), ]

# 2. 取前 5 行
write.csv(plot_kegg_data, "K:/WPY/SMC_Microbe_Positive_vs_Negative_kegg.csv")

if(nrow(plot_kegg_data) > 5) {
  plot_kegg_data <- plot_kegg_data[1:5, ]
}

# 3. 计算 GeneRatio 数值（因为 KEGG 默认是 "5/50" 这种字符串）
# 基础 R 的处理方式：
ratios <- strsplit(plot_kegg_data$GeneRatio, "/")
plot_kegg_data$RatioNum <- sapply(ratios, function(x) as.numeric(x[1]) / as.numeric(x[2]))

# 4. 格式化描述
plot_kegg_data$Description <- str_to_title(plot_kegg_data$Description)
# 修正缩写
plot_kegg_data$Description <- gsub("Tgf", "TGF", plot_kegg_data$Description)
plot_kegg_data$Description <- gsub("Mapk", "MAPK", plot_kegg_data$Description)

# 5. 绘图
p_kegg_top5 <- ggplot(plot_kegg_data, aes(x = RatioNum, y = reorder(Description, RatioNum))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "#B22222", high = "#F8766D") +
  labs(title = "Top 5 Enriched KEGG Pathways", x = "Gene Ratio", y = NULL, color = "P-adj") +
  theme_bw() +
  theme(panel.border = element_rect(colour = "black", fill=NA, linewidth=1.2),
        axis.text = element_text(size = 10, color = "black", face = "bold"))

print(p_kegg_top5)
ggsave("K:/WPY/Fig 1l.pdf", width = 8, height = 3.5, units = "in")













































library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(msigdbr)
library(ggplot2)
deg_smc$gene <- rownames(deg_smc)

deg_smc <- deg_smc %>%
  filter(!is.na(avg_log2FC), !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(p_val_adj, desc(avg_log2FC))

deg_up <- deg_smc %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.25)

deg_down <- deg_smc %>%
  filter(p_val_adj < 0.05, avg_log2FC < -0.25)

gene_up_symbol <- deg_up$gene
gene_down_symbol <- deg_down$gene

gene_up_entrez <- bitr(
  gene_up_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)

gene_down_entrez <- bitr(
  gene_down_symbol,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
)
ego_up <- enrichGO(
  gene = gene_up_symbol,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

ego_down <- enrichGO(
  gene = gene_down_symbol,
  OrgDb = org.Mm.eg.db,
  keyType = "SYMBOL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)
ekegg_up <- enrichKEGG(
  gene = gene_up_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)

ekegg_down <- enrichKEGG(
  gene = gene_down_entrez$ENTREZID,
  organism = "mmu",
  pvalueCutoff = 0.05
)
write.csv(as.data.frame(ego_up), "SMC_mouse_GO_BP_up.csv", row.names = FALSE)
write.csv(as.data.frame(ego_down), "SMC_mouse_GO_BP_down.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_up), "SMC_mouse_KEGG_up.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg_down), "SMC_mouse_KEGG_down.csv", row.names = FALSE)
dotplot(ego_up, showCategory = 15) + ggtitle("Mouse SMC GO BP Up")
dotplot(ekegg_up, showCategory = 15) + ggtitle("Mouse SMC KEGG Up")










deg_gsea <- deg_smc %>%
  filter(!is.na(avg_log2FC)) %>%
  group_by(gene) %>%
  slice_max(order_by = abs(avg_log2FC), n = 1) %>%
  ungroup()

geneList <- deg_gsea$avg_log2FC
names(geneList) <- deg_gsea$gene
geneList <- sort(geneList, decreasing = TRUE)
msig_h <- msigdbr(species = "Mus musculus", category = "H")
term2gene_h <- msig_h[, c("gs_name", "gene_symbol")]



gsea_h <- GSEA(
  geneList = geneList,
  TERM2GENE = term2gene_h,
  pvalueCutoff = 0.1,
  verbose = FALSE
)

gsea_h_res <- as.data.frame(gsea_h)
write.csv(gsea_h_res, "SMC_mouse_GSEA_Hallmark.csv", row.names = FALSE)



gsea_h_res %>%
  filter(ID %in% c(
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
    "HALLMARK_HYPOXIA",
    "HALLMARK_ANGIOGENESIS",
    "HALLMARK_APOPTOSIS",
    "HALLMARK_CHOLESTEROL_HOMEOSTASIS"
  )) %>%
  arrange(p.adjust)

