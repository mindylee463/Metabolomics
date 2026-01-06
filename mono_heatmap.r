###monoassociated
library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(openxlsx)

setwd('~/Desktop/Mei/monoassociate/')
ori_fm_data <- read_excel("Mei_LC-MS Metab Profile Raw data with IDs.xlsx", 
                          sheet = "Known_compound", col_names = FALSE)

rows_to_keep <- c(4:238)
metab <- ori_fm_data[rows_to_keep, ]
colnames(metab) <- as.character(ori_fm_data[3, ]) 

metab1 <- metab[colnames(metab)[c(1,12:19)]]

counts <- table(metab1$Metabolite.name)

# occurrencies < 2
keep_values <- names(counts[counts < 2])

# only keep these rows
metab_clean <- subset(metab1, `Metabolite.name` %in% keep_values)

###matrix
metab_clean <- metab_clean %>% 
  column_to_rownames(var = "Metabolite.name")

mat_metab <- t(metab_clean) 
dim(mat_metab)
mat_metab_num <- apply(mat_metab, 2, as.numeric) %>% as.data.frame()
rownames(mat_metab_num) <- rownames(mat_metab)

##heatmap
library(ComplexHeatmap)
library(circlize)
interest1<-c("FERULIC ACID",
             "M-COUMARIC ACID",
             "SYRINGIC ACID")
heatmap1<-mat_metab_num[,interest1]

group <- factor(c(rep("15% inWAX", 4), rep("15% cellulose", 4)))  # 修改为实际样本数量
names(group) <- rownames(heatmap1)

group_col <- c("15% cellulose" = "#92bfdb", "15% inWAX" = "#d73027")

# 4️⃣ 创建行注释对象
top_ha <- HeatmapAnnotation(
  Group = group,
  col = list(Group = group_col),
  show_annotation_name = FALSE,  # 不显示 title
  show_legend = FALSE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize=12, fontface="bold"),
    labels_gp = gpar(fontsize=10, fontface="bold")
  )
)
# 5️⃣ 行标准化（Z-score）
heatmap_z <- t(scale(heatmap1))
colnames(heatmap_z) <- c("1","2","3","4","1","2","3","4")
library(ComplexHeatmap)
library(circlize)
# 6️⃣ 绘制热图
pdf("heatmap1.pdf", width = 6.5, height = 3.5)
Heatmap(
  heatmap_z,
  name = "Z-score",
  rect_gp = gpar(col = "white", lwd = 0.1),
  col = colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,          # 代谢物聚类
  cluster_columns = FALSE,      # 样本顺序保持
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),       # 行名
  column_names_gp = gpar(fontsize = 11, fontface = "bold"), 
  column_split = group,         # ✅ 按样本分 block
  column_gap = unit(3, "mm"), border = TRUE, # ✅ 组间空隙
  top_annotation = top_ha       # 显示分组颜色条
)
dev.off() 

interest2 <- c("HYOCHOLIC ACID",
               "TAUROURSODEOXYCHOLIC ACID",
               "TAUROCHOLIC ACID",
               "CHOLIC ACID")

heatmap2<-mat_metab_num[,interest2]

heatmap_z2 <- t(scale(heatmap2))
colnames(heatmap_z2) <- c("1","2","3","4","1","2","3","4")
pdf("heatmap2.pdf", width = 6.5, height = 3.5)
Heatmap(
  heatmap_z2,
  name = "Z-score",
  rect_gp = gpar(col = "white", lwd = 0.1),
  col = colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,          # 代谢物聚类
  cluster_columns = FALSE,      # 样本顺序保持
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),       # 行名
  column_names_gp = gpar(fontsize = 11, fontface = "bold"), 
  column_split = group,         # ✅ 按样本分 block
  column_gap = unit(3, "mm"), border = TRUE, # ✅ 组间空隙
  top_annotation = top_ha       # 显示分组颜色条
)
dev.off() 


interest3 <- c("1-LINOLEOYL-SN-GLYCERO-3-PHOSPHORYLCHOLINE",
               "1-O-HEXADECYL-2-DEOXY-2-THIO-S-ACETYL-SN-GLYCERYL-3-PHOSPHORYLCHOLINE",
               "1-MYRISTOYL-SN-GLYCERO-3-PHOSPHOCHOLINE",
               "1-HEPTADECANOYL-SN-GLYCERO-3-PHOSPHOCHOLINE",
               "1-HEXADECYL-SN-GLYCERO-3-PHOSPHOCHOLINE",
               "1-PALMITOYL-SN-GLYCERO-3-PHOSPHOCHOLINE")

heatmap3<-mat_metab_num[,interest3]

heatmap_z3 <- t(scale(heatmap3))
colnames(heatmap_z3) <- c("1","2","3","4","1","2","3","4")
pdf("heatmap3.pdf", width = 8.7, height = 4.5)
Heatmap(
  heatmap_z3,
  name = "Z-score",
  rect_gp = gpar(col = "white", lwd = 0.1),
  col = colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,          # 代谢物聚类
  cluster_columns = FALSE,      # 样本顺序保持
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),       # 行名
  column_names_gp = gpar(fontsize = 11, fontface = "bold"), 
  column_split = group,         # ✅ 按样本分 block
  column_gap = unit(3, "mm"), border = TRUE, # ✅ 组间空隙
  top_annotation = top_ha       # 显示分组颜色条
)
dev.off() 


interest4 <- c("(2R)-3-HYDROXYISOVALEROYLCARNITINE",
               "(R)-BUTYRYLCARNITINE",
               "ACETYL-L-CARNITINE",
               "ADIPOYL-L-CARNITINE",
               "DECANOYLCARNITINE",
               "HEXANOYL-L-CARNITINE",
               "OCTANOYLCARNITINE")

heatmap4<-mat_metab_num[,interest4]

heatmap_z4 <- t(scale(heatmap4))
colnames(heatmap_z4) <- c("1","2","3","4","1","2","3","4")
pdf("heatmap4.pdf", width = 8.5, height = 4.5)
Heatmap(
  heatmap_z4,
  name = "Z-score",
  rect_gp = gpar(col = "white", lwd = 0.1),
  col = colorRamp2(c(-2,0,2), c("navy", "white", "firebrick3")),
  cluster_rows = FALSE,          # 代谢物聚类
  cluster_columns = FALSE,      # 样本顺序保持
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_names_rot = 0,
  row_names_gp = gpar(fontsize = 10, fontface = "bold"),       # 行名
  column_names_gp = gpar(fontsize = 11, fontface = "bold"), 
  column_split = group,         # ✅ 按样本分 block
  column_gap = unit(3, "mm"), border = TRUE, # ✅ 组间空隙
  top_annotation = top_ha       # 显示分组颜色条
)
dev.off() 
