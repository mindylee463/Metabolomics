library(readxl)
library(tibble)
library(dplyr)
library(tidyr)
library(patchwork)
library(ggplot2)
library(openxlsx)

setwd('~/Desktop/Mei/Jackson HFD mouse metabolomics data/10.10.25/Ttest_v2/')
ori_fm_data <- read_excel("../fecal metabolomics dry weight normalized.10.10.25.xlsx", 
                          sheet = "Known", col_names = FALSE)

rows_to_keep <- c(10:632)
metab <- ori_fm_data[rows_to_keep, ]
colnames(metab) <- as.character(ori_fm_data[9, ]) 

metab1 <- metab[colnames(metab)[c(1,16:30)]]

counts <- table(metab1$`Metabolite name`)

# occurrencies < 2
keep_values <- names(counts[counts < 2])

# only keep these rows
metab_clean <- subset(metab1, `Metabolite name` %in% keep_values)

###matrix
metab_clean <- metab_clean %>% 
  column_to_rownames(var = "Metabolite name")

mat_metab <- t(metab_clean) 
dim(mat_metab)
mat_metab_num <- apply(mat_metab, 2, as.numeric) %>% as.data.frame()
rownames(mat_metab_num) <- rownames(mat_metab)
sample_info <- data.frame(SampleID = rownames(mat_metab_num)) %>%
  mutate(Group = substr(SampleID, 1, 2)) %>% mutate(SampleType = 'Fecal_dwn') 


##test

test_with_log2fc <- function(mat, sample_info) {
  
  groups <- unique(sample_info$Group)
  comparisons <- combn(groups, 2, simplify = FALSE)  ## 生成所有两两组合
  
  res <- data.frame(Metabolite = character(),
                    Group1 = character(),
                    Group2 = character(),
                    log2FoldChange = numeric(),
                    p_value = numeric(),
                    stringsAsFactors = FALSE)
  
  for(cmp in comparisons){
    g1 <- cmp[1]
    g2 <- cmp[2]
    
    for(i in seq_len(ncol(mat))){
      # 取样本
      values1 <- mat[sample_info$SampleID[sample_info$Group == g1], i]
      values2 <- mat[sample_info$SampleID[sample_info$Group == g2], i]
      
      # 去掉 0
      values1_nonzero <- values1[values1 != 0]
      values2_nonzero <- values2[values2 != 0]
      
      # 如果去掉 0 后没有值或方差为 0，跳过
      if(length(values1_nonzero) == 0 || length(values2_nonzero) == 0 ||
         (var(values1_nonzero) == 0 && var(values2_nonzero) == 0)){
        next
      }
      
      # 先计算平均值，再算 log2FC（用原始数值）
      mean1 <- mean(values1_nonzero, na.rm = TRUE)
      mean2 <- mean(values2_nonzero, na.rm = TRUE)
      avgfc <- mean2 / mean1
      log2fc <- log2(mean2 / mean1)
      
      # 对剩余非零值做 log2 转换，再做 t-test
      values1_log2 <- log2(values1_nonzero)
      values2_log2 <- log2(values2_nonzero)
      
      test <- t.test(values1_log2, values2_log2, var.equal = TRUE, alternative = "two.sided")
      
      res <- rbind(res, data.frame(
        Metabolite = colnames(mat)[i],
        Group1 = g1,
        Group2 = g2,
        avgfc = avgfc,
        log2FoldChange = log2fc,
        p_value = test$p.value,
        log_pval = -log10(test$p.value),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  res$q_value <- p.adjust(res$p_value, method = "fdr")
  
  return(res)
}

result <- test_with_log2fc(mat_metab_num, sample_info)

result[result$Metabolite=="FERULATE",]

p_cutoff <- 0.05
fc_cutoff <- 0


result <- result %>%
  mutate(
    Significance = ifelse(p_value < p_cutoff & abs(log2FoldChange) > fc_cutoff, 
                          ifelse(log2FoldChange > 0, "Significant Up", "Significant Down"),
                          "Not significant")
  )

CFvsHF <- subset(result, Group1 == "HF" & Group2 == "CF")
WFvsHF <- subset(result, Group1 == "HF" & Group2 == "WF")
WFvsCF <- subset(result, Group1 == "CF" & Group2 == "WF")

wb <- createWorkbook()
addWorksheet(wb, "CFvsHF")
writeData(wb, "CFvsHF", CFvsHF)

addWorksheet(wb, "WFvsHF")
writeData(wb, "WFvsHF", WFvsHF)

addWorksheet(wb, "WFvsCF")
writeData(wb, "WFvsCF", WFvsCF)

saveWorkbook(wb, 'dry_weight_metabolite_Ttest_results20251105.xlsx', overwrite = TRUE)

create_volcano_plot <- function(data, 
                                p_cutoff = 0.05, 
                                fc_cutoff = 0.5,
                                top_n = 5,
                                title = "Volcano Plot",
                                x_limits = c(-10, 10),
                                x_breaks = seq(-10, 10, by = 2)) {
  
  if (!require(ggplot2)) install.packages("ggplot2"); library(ggplot2)
  if (!require(ggrepel)) install.packages("ggrepel"); library(ggrepel)
  if (!require(dplyr)) install.packages("dplyr"); library(dplyr)
  
  get_top_metabolites <- function(data, direction, n = 5) {
    result <- data %>%
      filter(Significance == paste("Significant", direction)) %>%
      arrange(p_value)
    
    if(nrow(result) > n) {
      return(head(result, n))
    } else {
      return(result)
    }
  }
  
  
  top_up <- get_top_metabolites(data, "Up", top_n)
  top_down <- get_top_metabolites(data, "Down", top_n)
  top_metabolites <- bind_rows(top_up, top_down)
  
  interesting_metabolites <- c("FERULATE")
  
  top_metabolites <- bind_rows(
    top_up,
    top_down,
    data %>% filter(Metabolite %in% interesting_metabolites)
  ) %>% distinct()
 
  volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(p_value))) +
    geom_point(aes(color = Significance), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = c(
      "Significant Up" = "red", 
      "Significant Down" = "blue",
      "Not significant" = "grey"
    )) +
    geom_hline(yintercept = -log10(p_cutoff), linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "black") +
    geom_point(
      data = data %>% filter(Metabolite %in% interesting_metabolites),
      color = "yellow", shape = 21, size = 3, fill = "yellow"
    )+
    geom_text_repel(
      data = top_metabolites,
      aes(label = Metabolite, color = 'black'),
      size = 1.5,
      box.padding = 0.6,
      point.padding = 0.3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    
    labs(
      title = title,
      subtitle = paste("Top", top_n, "up/down regulated metabolites labeled"),
      x = "log2 Fold Change", 
      y = "-log10(p-value)"
    ) +
    theme_bw() +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks
    )
  
  return(volcano_plot)
}

volcano1 <- create_volcano_plot(
  data = CFvsHF,
  p_cutoff = 0.05,           
  fc_cutoff = 1,             
  top_n = 20,                 
  title = "CF vs HF",
  x_limits = c(-10, 10),       
  x_breaks = seq(-10, 10, by = 2)  
)

volcano1
ggsave("volcano_plot_CF_vs_HF.pdf",  volcano1,  width = 10,  height = 8)

volcano2 <- create_volcano_plot(
  data = WFvsHF,
  p_cutoff = 0.05,           
  fc_cutoff = 1,             
  top_n = 20,                 
  title = "WF vs HF",
  x_limits = c(-10, 10),       
  x_breaks = seq(-10, 10, by = 2) 
)

volcano2
ggsave("volcano_plot_WF_vs_HF.pdf",  volcano2,  width = 10,  height = 8)

volcano3 <- create_volcano_plot(
  data = WFvsCF,
  p_cutoff = 0.05,           
  fc_cutoff = 1,            
  top_n = 20,                 
  title = "WF vs CF",
  x_limits = c(-10, 10),       
  x_breaks = seq(-10, 10, by = 2)  
)

volcano3
ggsave("volcano_plot_WF_vs_CF.pdf",  volcano3,  width = 10,  height = 8)


##Task4


library(UpSetR)

# Prepare a list of sets

met_sets <- list(
  CFvsHF_Up = CFvsHF[CFvsHF$Significance == "Significant Up", "Metabolite"],
  CFvsHF_Down = CFvsHF[CFvsHF$Significance == "Significant Down", "Metabolite"],
  WFvsHF_Up = WFvsHF[WFvsHF$Significance == "Significant Up", "Metabolite"],
  WFvsHF_Down = WFvsHF[WFvsHF$Significance == "Significant Down", "Metabolite"],
  WFvsCF_Up = WFvsCF[WFvsCF$Significance == "Significant Up", "Metabolite"],
  WFvsCF_Down = WFvsCF[WFvsCF$Significance == "Significant Down", "Metabolite"]
)


library(ggVennDiagram)
library(patchwork)

CFvsHF_metabolite = CFvsHF[CFvsHF$Significance != "Not significant", "Metabolite"]
WFvsHF_metabolite = WFvsHF[WFvsHF$Significance != "Not significant", "Metabolite"]
#WFvsCF_metabolite = WFvsCF[WFvsCF$Significance != "Not significant", "Metabolite"]

venn_list <- list(CFvsHF = CFvsHF_metabolite, 
                  WFvsHF = WFvsHF_metabolite
                  #WFvsCF = WFvsCF_metabolite
                  )

# Create ggplot Venn diagram
library(ggvenn)

p12_gg <- ggvenn(
  venn_list,
  fill_color = c("#7EA6D9", "#D4E6BC"),
  fill_alpha = 0.7,
  stroke_size = 0.8,
  stroke_color = "grey40",
  #   set_name = c("CF vs HF", "WF vs HF", "WF vs CF"),
  set_name_size = 5,
  set_name_color = "black",
  text_size = 4,
  show_percentage = FALSE
)
# Show plot
p12_gg
ggsave("Venn_plot_Ttest_20251105.pdf", plot = p12_gg, width = 8, height = 6, dpi = 300)
# Draw UpSet plot
pdf("UpsetR_plot_Ttest_20251105.pdf", width = 8, height = 6)
library(UpSetR)

upset(
  fromList(met_sets),
  sets = c("CFvsHF_Up","CFvsHF_Down","WFvsHF_Up","WFvsHF_Down"),
  nsets = 6, 
  order.by = "freq", 
  keep.order = TRUE,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = "#F5A889",
  main.bar.color = "#ACD6EC",
  text.scale = c(1.5, 1, 1.5, 1.5, 1.1, 1.5), # 控制字体大小
  number.angles = 0
)

dev.off()

##Table
library(dplyr)

sets <- names(met_sets)
all_elements <- unique(unlist(met_sets))

overlap_result <- data.frame(Element = all_elements, stringsAsFactors = FALSE)

for(s in sets){
  overlap_result[[s]] <- all_elements %in% met_sets[[s]]
}

print(overlap_result)

intersect_elements <- overlap_result %>%
  rowwise() %>%
  mutate(n_sets = sum(c_across(all_of(sets)))) %>%
  filter(n_sets > 1) %>%
  select(Element, n_sets, all_of(sets))

unique_elements <- overlap_result %>%
  rowwise() %>%
  mutate(n_sets = sum(c_across(all_of(sets)))) %>%
  filter(n_sets == 1) %>%
  select(Element, n_sets, all_of(sets))

all_met <- bind_rows(
  intersect_elements %>% mutate(Group = "Intersection"),
  unique_elements %>% mutate(Group = "Unique")
) 

matched_data <- result %>%
  filter(Metabolite %in% all_met$Element)

final_data <- matched_data %>%
  left_join(all_met, by = c("Metabolite" = "Element"))


head(final_data)

# Save as CSV
write.xlsx(final_data,
           file = "Table_WF_CF_HF_overlap.xlsx",
           sheetName = "Metabolite_Data",
           overwrite = TRUE)

###task5

function_metadata <- metab[,c(1:3,5)]

merged_table <- CFvsHF %>%
  left_join(function_metadata, 
            by = c("Metabolite" = "Metabolite name"))

write.table( merged_table, file = "CFvsHF_function_Ttest_20251105.txt", sep = "\t", row.names = FALSE,quote = FALSE)

merged_table <- WFvsHF %>%
  left_join(function_metadata, 
            by = c("Metabolite" = "Metabolite name"))

write.table( merged_table, file = "WFvsHF_function_Ttest_20251105.txt", sep = "\t", row.names = FALSE,quote = FALSE)

merged_table <- WFvsCF %>%
  left_join(function_metadata, 
            by = c("Metabolite" = "Metabolite name"))

write.table( merged_table, file = "WFvsCF_function_Ttest_20251105.txt", sep = "\t", row.names = FALSE,quote = FALSE)


###Task3
classification <- read.csv("../compund_classification.tsv",
                           header = TRUE,
                           stringsAsFactors = FALSE,
                           quote = "\"",
                           sep = "\t",
                           fileEncoding = "UTF-8")



get_class_counts <- function(compound_list, label) {
  classification %>%
    dplyr::filter(Compound %in% compound_list) %>%
    dplyr::count(Category, name = label) %>%
    dplyr::mutate(
      Percentage = !!rlang::sym(label) / sum(!!rlang::sym(label)) * 100,
      Label = paste0(
        Category, "\n", "(", sprintf("%.2f", Percentage), "%)"
      )
    )
}

cf_up   <- get_class_counts(met_sets$CFvsHF_Up,   "CFvsHF_Up")
cf_down <- get_class_counts(met_sets$CFvsHF_Down, "CFvsHF_Down")
wf_up   <- get_class_counts(met_sets$WFvsHF_Up,   "WFvsHF_Up")
wf_down <- get_class_counts(met_sets$WFvsHF_Down, "WFvsHF_Down")

# ==== 3. 合并为 Table 2 ====
table2 <- full_join(cf_up, cf_down, by = "Category") %>%
  full_join(., wf_up, by = "Category") %>%
  full_join(., wf_down, by = "Category") %>%
  replace_na(list(
    CFvsHF_Up = 0, CFvsHF_Down = 0,
    WFvsHF_Up = 0, WFvsHF_Down = 0
  )) %>%
  arrange(Category)

# 保存结果
write_csv(table2, "Table2_metabolite_class_summary_20251105.csv")

# ==== 4. 绘制饼图 ====

all_categories <- sort(unique(c(
  cf_up$Category,
  cf_down$Category,
  wf_up$Category,
  wf_down$Category
)))

# 生成统一颜色
base_colors <- brewer.pal(9, "Set3")

# 插值生成 9 个颜色
#colors <- colorRampPalette(base_colors)(length(all_categories))
colors <- base_colors
names(colors) <- all_categories

pdf("CF_vs_HF_&_WF_vs_HF_pie.pdf", width = 12, height = 10) 

par(mfrow = c(2, 2))
pie(cf_up$CFvsHF_Up,
    labels = cf_up$Label,
    col = colors[cf_up$Category],  # 用统一颜色
    main = "CF vs HF Up-regulated",
    border = "white")

# 绘图示例：WF Up
pie(wf_up$WFvsHF_Up,
    labels = wf_up$Label,
    col = colors[wf_up$Category],
    main = "WF vs HF Up-regulated",
    border = "white")

pie(cf_down$CFvsHF_Down,
    labels = cf_down$Label,
    col = colors[cf_down$Category],  # 用统一颜色
    main = "CF vs HF Down-regulated",
    border = "white")

# 绘图示例：WF Up
pie(wf_down$WFvsHF_Down,
    labels = wf_down$Label,
    col = colors[wf_down$Category],
    main = "WF vs HF Down-regulated",
    border = "white")

dev.off()


shared <- intersect(CFvsHF_metabolite, WFvsHF_metabolite)

# Only in each comparison
only_CF <- setdiff(CFvsHF_metabolite, WFvsHF_metabolite)
only_WF <- setdiff(WFvsHF_metabolite, CFvsHF_metabolite)



table3 <- data.frame(
  Metabolite = c(only_CF, only_WF, shared),
  Label = c(
    rep("Only CFvsHF", length(only_CF)),
    rep("Only WFvsHF", length(only_WF)),
    rep("Shared", length(shared))
  )
)
table3 <- merge(table3, classification[, c("Compound", "Category")],
                by.x = "Metabolite", by.y = "Compound", all.x = TRUE)

# Sort nicely
table3 <- table3[order(table3$Label, table3$Category, table3$Metabolite), ]

# Save as CSV
write.csv(table3, "Table3_CF_WF_overlap.csv", row.names = FALSE)


##heatmap
library(ComplexHeatmap)
library(circlize)
interest1<-c("FERULATE",
             "2-HYDROXYCINNAMIC ACID, PREDOMINANTLY TRANS",
             "ENTEROLACTONE",
             "ENTERODIOL",
             "7-O-METHYLERIODICTYOL",
             "3,4-DIHYDROXYACETOPHENONE",
             "JACEOSIDE",
             "NICOTINIC ACID")
heatmap1<-mat_metab_num[,interest1]

group <- factor(c(rep("HFD", 5), rep("HFD+cellulose", 5), rep("HFD+inWAX", 5)))  # 修改为实际样本数量
names(group) <- rownames(heatmap1)

group_col <- c("HFD" = "#4575b4", "HFD+cellulose" = "#92bfdb", "HFD+inWAX" = "#d73027")

# 4️⃣ 创建行注释对象
top_ha <- HeatmapAnnotation(
  Group = group,
  col = list(Group = group_col),
  show_annotation_name = FALSE,  # 不显示 title
  show_legend = FALSE,
  annotation_legend_param = list(
    title_gp = gpar(fontsize=12, fontface="bold"),
    labels_gp = gpar(fontsize=10)
  )
)
# 5️⃣ 行标准化（Z-score）
heatmap_z <- t(scale(heatmap1))
colnames(heatmap_z) <- c("1","2","3","4","5","1","2","3","4","5","1","2","3","4","5")
library(ComplexHeatmap)
library(circlize)
# 6️⃣ 绘制热图
pdf("heatmap1.pdf", width = 8.5, height = 3.5)
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


interest2 <- c("TYROSINE",
               "TRYPTOPHAN",
               "ASPARTATE",
               "GLYCINE",
               "HISTIDINE",
               "GLUTAMINE",
               "GLUTAMATE",
               "PHENYLALANINE",
               "N-ACETYL-L-GLUTAMINE",
               "N-ACETYLASPARTIC ACID")

heatmap2<-mat_metab_num[,interest2]

heatmap_z2 <- t(scale(heatmap2))
colnames(heatmap_z2) <- c("1","2","3","4","5","1","2","3","4","5","1","2","3","4","5")
pdf("heatmap2.pdf", width = 8, height = 4)
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

