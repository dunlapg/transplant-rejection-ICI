library(immunarch)
library(Seurat)
library(scRepertoire)
library(enrichR)
library(slingshot)
library(ggplot2)
library(gghighlight)
library(ggVennDiagram)
library(tibble)
library(tidyr)
library(RColorBrewer)
library(scales)
library(EnhancedVolcano)
library(cowplot)
library(ggpubr)
source("manuscript/AMP_manuscript_functions.R")

#load working directory
setwd("kidney_ICI")

#creatinine values over time
crvalues <- read.csv("cr_data_time.csv")
ggplot(crvalues, aes(Day, Cr)) + 
  geom_line() + geom_point() +
  scale_y_continuous(limits = c(0,4)) +
  theme_classic() +
  ylab("Cr (mg/dL)")
ggsave("analysis/figs/creatine_timecourse.png", width = 6, height = 3)

# Bulk TCR sequencing -----------------------------------------------------
#load TCR bulk data
bulktcr_kid <- repLoad(c("tcrb/data/GM-18-18046_KID_TCRB.tsv", "tcrb/data/GM-17-7547_pre_LN_TCRB.tsv", "tcrb/data/GM-19-7655_post_TUM_TCRB.tsv", "tcrb/data/GM-Pre_PBMC_180728_TCRB.tsv", "tcrb/data/GM-C1_PBMC_180816_TCRB.tsv", "tcrb/data/GM-C2_PBMC_180927_TCRB.tsv"))
bulktcr_kid$meta$sampleType <- c("Tissue", "Tissue", "Tissue", "Blood", "Blood", "Blood")
bulktcr_kid$meta$sampleShort <- c("Kidney", "LN Met", "Skin Met", "Pre Pembro", "Pembro 1", "Pembro 2")
bulktcr_kid$meta$sampleOrder <- c("Sample 1", "Sample 2", "Sample 3", "Sample 4", "Sample 5", "Sample 6")

bulktcr_cdr3_counts <- as.tibble(read.csv("analysis/data/bulk_uniqueCDR3_counts.csv"))
bulktcr_cdr3_freqs <- as.tibble(read.csv("analysis/data/bulk_uniqueCDR3_freqs.csv"))

bulktcr_cdr3_counts <- bulktcr_cdr3_counts[!grepl('\\*', bulktcr_cdr3_counts$CDR3aa), ]
bulktcr_cdr3_freqs <- bulktcr_cdr3_freqs[!grepl('\\*', bulktcr_cdr3_freqs$CDR3aa), ]

#adding in reactive TCRs from mcPAS
mcpas_cmv = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB", .pathology = "Cytomegalovirus (CMV)")
bulktcr_cmv<- dbAnnotate(bulktcr_kid$data, mcpas_cmv, "CDR3.aa", "CDR3.beta.aa")
mcpas_all = dbLoad("https://gitlab.com/immunomind/immunarch/raw/dev-0.5.0/private/McPAS-TCR.csv.gz", "mcpas", .species = "Human", .chain = "TRB")
bulktcr_mcpas_all <- dbAnnotate(bulktcr_kid$data, mcpas_all, "CDR3.aa", "CDR3.beta.aa")

#cross reference data to melanoma biopsies from pruessmann et al, 2020 (doi: 10.1038/s43018-019-0019-5)
bulktcr_pruessmann <- repLoad(c("tcrb/data/pruessmann_melanoma_biopsy_tcr"))
bulktcr_pruessmann <- pubRep(bulktcr_pruessmann$data, "aa", .verbose = T)
bulktcr_pruessmann$cloneSums <- rowSums(bulktcr_pruessmann[,3:201], na.rm = T)

#denote top kidney clones
kidney_cdr3_top25 <- bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa[1:25]

#add all as new columns
bulktcr_cdr3_counts <- bulktcr_cdr3_counts %>% 
  mutate(CMV_reactive = ifelse(CDR3aa %in% bulktcr_cmv$CDR3.aa, "Yes", "No")) %>%
  mutate(Kidney_Top25 = ifelse(CDR3aa %in% kidney_cdr3_top25, "Yes", "No")) %>%
  mutate(pruessmann_clone = if_else(CDR3aa %in% bulktcr_pruessmann$CDR3.aa, "Yes", "No"))
bulktcr_cdr3_freqs <- bulktcr_cdr3_freqs %>% 
  mutate(CMV_reactive = ifelse(CDR3aa %in% bulktcr_cmv$CDR3.aa, "Yes", "No")) %>%
  mutate(Kidney_Top25 = ifelse(CDR3aa %in% kidney_cdr3_top25, "Yes", "No")) %>%
  mutate(pruessmann_clone = if_else(CDR3aa %in% bulktcr_pruessmann$CDR3.aa, "Yes", "No"))


#bulktcrcolors <- c("#8DD3C7", "#FFED6F", "#FDB462", "#80B1D3", "#BEBADA", "#FB8072") #indiv colors for each sample

bulktcrcolors <- c("#8E44AD", "#6C3483", "#4A235A", "#C0392B", "#922B21", "#641E16") #colors hued by sample type
show_col(bulktcrcolors)

#number of clonotypes
exp_vol <- repExplore(bulktcr_kid$data, .method = "volume")

gather(bulktcr_cdr3_counts, sample, counts, kidney:pembro_2, factor_key=TRUE) %>% 
  filter(counts > 0) %>% 
  group_by(sample) %>% 
  count() %>% 
  ggplot(aes(x = sample, y = n, fill = sample)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_fill_manual(values = bulktcrcolors) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12)) +
  guides(fill=FALSE) +
  xlab(NULL) + ylab("Unique Clonotypes") + ylim(0,50000)
ggsave("analysis/figs/bulk_clonotypeCount_bySample.eps", width = 5, height = 4)

vis(exp_vol, .by = c("sampleType"), .meta = bulktcr_kid$meta, .test = FALSE) +
  scale_x_discrete(limits = c("Tissue", "Blood")) +
  scale_fill_manual(values = c("#D11F1F", "#8785BA")) +
  labs(subtitle = NULL) + 
  xlab(NULL) + 
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12)) +
  ggtitle(NULL)
ggsave("analysis/figs/bulk_clonotypeCount_bySampleType.eps", width = 2.5, height = 4)

#number of clones
exp_clones <- repExplore(bulktcr_kid$data, .method = "clones")

gather(bulktcr_cdr3_counts, sample, counts, kidney:pembro_2, factor_key=TRUE) %>% 
  filter(counts > 0) %>% 
  group_by(sample) %>% 
  tally(counts) %>% 
  ggplot(aes(x = sample, y = n, fill = sample)) + 
  geom_bar(stat = "identity", color = "black") +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_fill_manual(values = bulktcrcolors) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12)) +
  guides(fill=FALSE) +
  xlab(NULL) + ylab("Clones") + ylim(0,70000)
ggsave("analysis/figs/bulk_cellCount_bySample.eps", width = 5, height = 4)

vis(exp_clones, .by = c("sampleType"), .meta = bulktcr_kid$meta, .test = FALSE) +
  scale_x_discrete(limits = c("Tissue", "Blood")) +
  scale_fill_manual(values = c("#D11F1F", "#8785BA")) +
  labs(subtitle = NULL) + 
  xlab(NULL) + 
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12)) +
  ggtitle(NULL)
ggsave("analysis/figs/bulk_cellCount_bySampleType.eps", width = 2.5, height = 4)

#clonality per sample
imm_rare <- repClonality(bulktcr_kid$data, .method = "rare", .perc = T)
vis(imm_rare, .meta = bulktcr_kid$meta, .test = FALSE, .points = FALSE) +
  labs(title = NULL, subtitle = NULL) + 
  theme_classic() +
  scale_fill_brewer(palette = "Blues", name = "Clone Count", labels = c("1", "1 < X <= 3", "3 < X <= 10", "10 < X <= 30", "30 < X <= 100", "100 < X")) +
  scale_x_discrete(limits = bulktcr_kid$meta$Sample, labels = bulktcr_kid$meta$sampleShort) +
  ylab(NULL) + xlab(NULL) +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=11))  + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  scale_y_reverse(breaks=c(0, 0.25, 0.50, 0.75, 1), labels = c("100%", "75%", "50%", "25%", "0%"))
ggsave("analysis/figs/bulk_clonality_bySample_stacked.eps", width = 5, height = 3.5)

vis(imm_rare, .by= "sampleOrder", .meta = bulktcr_kid$meta, .test = FALSE, .points = FALSE) +
  labs(title = NULL, subtitle = NULL) + 
  theme_classic() +
  scale_fill_manual(values = bulktcrcolors, labels = bulktcr_kid$meta$sampleShort, name = NULL) +
  xlab("Clone Count") +
  theme(axis.text=element_text(size=12))
ggsave("analysis/figs/bulk_clonality_bySample_grouped.eps", width = 6, height = 3.5)

#clonal overlap
bulktcr_cdr3 <- list(kidney = bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa, ln = bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa, skin = bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa, pre_pembro = bulktcr_kid$data$`GM-Pre_PBMC_180728_TCRB`$CDR3.aa, pembro_1 = bulktcr_kid$data$`GM-C1_PBMC_180816_TCRB`$CDR3.aa, pembro_2 = bulktcr_kid$data$`GM-C2_PBMC_180927_TCRB`$CDR3.aa)

ggVennDiagram(bulktcr_cdr3[4:6], edge_size = 4, set_size = 4, label_alpha = 0, category.names = c("Pre Pembro", "Pembro C1", "Pembro C2")) + scale_colour_manual(values = bulktcrcolors[4:6]) + 
  scale_fill_gradient2(low = "white", high = "white") + 
  theme(legend.position = "none") 
ggsave("analysis/figs/bulk_bloodCollection_venn.eps", width = 4.5, height = 4.5)

bulktcr_kid_cloneOver10 <- bulktcr_kid
bulktcr_kid_cloneOver10$data$`GM-18-18046_KID_TCRB` <- bulktcr_kid_cloneOver10$data$`GM-18-18046_KID_TCRB` %>% filter(Clones > 10)
bulktcr_kid_cloneOver10$data$`GM-17-7547_pre_LN_TCRB` <- bulktcr_kid_cloneOver10$data$`GM-17-7547_pre_LN_TCRB` %>% filter(Clones > 10)
bulktcr_kid_cloneOver10$data$`GM-19-7655_post_TUM_TCRB` <- bulktcr_kid_cloneOver10$data$`GM-19-7655_post_TUM_TCRB` %>% filter(Clones > 10)
bulktcr_kid_cloneOver10$data$`GM-Pre_PBMC_180728_TCRB` <- bulktcr_kid_cloneOver10$data$`GM-Pre_PBMC_180728_TCRB` %>% filter(Clones > 10)
bulktcr_kid_cloneOver10$data$`GM-C1_PBMC_180816_TCRB` <- bulktcr_kid_cloneOver10$data$`GM-C1_PBMC_180816_TCRB` %>% filter(Clones > 10)
bulktcr_kid_cloneOver10$data$`GM-C2_PBMC_180927_TCRB` <- bulktcr_kid_cloneOver10$data$`GM-C2_PBMC_180927_TCRB` %>% filter(Clones > 10)

bulktcr_cloneOver10_cdr3 <- list(kidney = bulktcr_kid_cloneOver10$data$`GM-18-18046_KID_TCRB`$CDR3.aa, ln = bulktcr_kid_cloneOver10$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa, skin = bulktcr_kid_cloneOver10$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa, pre_pembro = bulktcr_kid_cloneOver10$data$`GM-Pre_PBMC_180728_TCRB`$CDR3.aa, pembro_1 = bulktcr_kid_cloneOver10$data$`GM-C1_PBMC_180816_TCRB`$CDR3.aa, pembro_2 = bulktcr_kid_cloneOver10$data$`GM-C2_PBMC_180927_TCRB`$CDR3.aa)

ggVennDiagram(bulktcr_cloneOver10_cdr3[4:6], edge_size = 4, set_size = 4, label_alpha = 0) + scale_colour_manual(values = bulktcrcolors[4:6]) + 
  scale_fill_gradient2(low = "white", high = "white") + 
  theme(legend.position = "none")

imm_ovPub <- repOverlap(bulktcr_kid$data, .method = "public", .verbose = F)
vis(imm_ovPub) + scale_y_discrete(labels=rev(bulktcr_kid$meta$sampleShort)) + 
  scale_x_discrete(labels=bulktcr_kid$meta$sampleShort) + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
  guides(fill = guide_colorbar(title="# Clonotypes")) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  ggtitle(NULL) + xlab(NULL) + ylab(NULL)
ggsave("analysis/figs/bulk_overlap_shared.eps", width = 6, height = 6)

imm_ovMor <- repOverlap(bulktcr_kid$data, .method = "morisita", .verbose = F)
vis(imm_ovMor) + scale_y_discrete(labels=rev(bulktcr_kid$meta$sampleShort)) + 
  scale_x_discrete(labels=bulktcr_kid$meta$sampleShort) + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
  guides(fill = guide_colorbar(title="Morisita \nIndex")) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text=element_text(size=12),
        legend.text =element_text(size=10)) +
  ggtitle(NULL) + xlab(NULL) + ylab(NULL)
ggsave("analysis/figs/bulk_overlap_morisita.eps", width = 6, height = 6)

#clonal diversity
div_invsimp <- repDiversity(bulktcr_kid$data, "inv.simp")

ggplot(div_invsimp, aes(x=Sample, y=Value, fill=Sample)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = bulktcrcolors) +
  labs(subtitle = NULL) + 
  ylab("Diversity Value (Inverse Simpson)") + 
  xlab(NULL) + 
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12)) +
  scale_x_discrete(labels=bulktcr_kid$meta$sampleShort) 
ggsave("analysis/figs/bulk_diversity_invsimp.eps", width = 5, height = 3.5)

divsimp_tissue <- div_invsimp[1:3,]
divsimp_blood <- div_invsimp[4:6,]

divsimp_tissue <- ggplot(divsimp_tissue, aes(x=Sample, y=Value, fill=Sample)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = bulktcrcolors[1:3]) +
  labs(subtitle = NULL) + 
  ylab("Diversity Value (Inverse Simpson)") + 
  xlab(NULL) + 
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=bulktcr_kid$meta$sampleShort[1:3]) 

divsimp_blood <- ggplot(divsimp_blood, aes(x=Sample, y=Value, fill=Sample)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = bulktcrcolors[4:6]) +
  labs(subtitle = NULL) + 
  ylab("Diversity Value (Inverse Simpson)") + 
  xlab(NULL) + 
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=12)) +
  scale_x_discrete(labels=bulktcr_kid$meta$sampleShort[4:6]) 

diversityplots <- align_plots(divsimp_tissue, divsimp_blood, align="hv", axis="tblr")
p1x <- ggdraw(diversityplots[[1]])
p2x <- ggdraw(diversityplots[[2]])

save_plot("analysis/figs/bulk_diversity_invsimp_tissue.eps", p1x, base_width = 2.4, base_height = 3.5)
save_plot("analysis/figs/bulk_diversity_invsimp_blood.eps", p2x, base_width = 2.4, base_height = 3.5)


div_invsimp_pbmc <- div_invsimp[4:6,]
div_invsimp_pbmc$Date <- c("2018-07-28", "2018-08-16", "2018-09-27")
div_invsimp_pbmc$Date <- as.Date(div_invsimp_pbmc$Date)
ggplot(div_invsimp_pbmc, aes(x=Date, y=Value, group = 1)) +
  geom_line(size = 1) + theme_classic() + 
  ylim(0,3000) + 
  ylab("Diversity Value (Inverse Simpson)") + 
  geom_point(aes(color=Sample), size = 5) + 
  scale_color_manual(labels = c("Pre Pembro", "After Pembro C1", "After Pembro C2"), values = bulktcrcolors[4:6]) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") 
ggsave("analysis/figs/bulk_diversity_invsimp_pbmcTimline.eps", width = 7, height = 3.5)

div_truediv <- repDiversity(bulktcr_kid$data, "div")

ggplot(div_truediv, aes(x=Sample, y=Value, fill=Sample)) +
  geom_bar(stat = "identity", color = "black") +
  scale_fill_manual(values = bulktcrcolors) +
  labs(subtitle = NULL) + 
  ylab("Diversity Index") + 
  xlab(NULL) + 
  guides(fill=FALSE) +
  theme_classic() +
  theme(axis.text.x  = element_text(angle=45, hjust = 1)) +
  scale_x_discrete(labels=bulktcr_kid$meta$sampleShort) 
ggsave("analysis/figs/bulk_diversity_true.eps", width = 5, height = 3.5)

div_truediv_pbmc <- div_truediv[4:6,]
div_truediv_pbmc$Date <- c("2018-07-28", "2018-08-16", "2018-09-27")
div_truediv_pbmc$Date <- as.Date(div_truediv_pbmc$Date)
ggplot(div_truediv_pbmc, aes(x=Date, y=Value, group = 1)) +
  geom_line(size = 1) + theme_classic() + 
  ylim(0,500) + 
  ylab("Diversity Index") + 
  geom_point(aes(color=Sample), size = 5) + 
  scale_color_manual(labels = c("Pre Pembro", "After Pembro C1", "After Pembro C2"), values = bulktcrcolors[4:6]) +
  scale_x_date(date_labels = "%m-%d", date_breaks = "1 week") 
ggsave("analysis/figs/bulk_diversity_true_pbmcTimline.eps", width = 7, height = 3.5)

#expansion and contraction plots
bulktcr_cdr3_counts %>% 
  mutate(expansion = pembro_1/pre_pembro) %>%
  mutate(type = ifelse(expansion > 2, "Expanded", ifelse(expansion < 0.5, "Contracted", "Stable"))) %>%
  mutate(type = ifelse(is.na(type), "Absent", type)) %>%
  select(CDR3aa, pembro_1, pre_pembro, expansion, type) %>%
  ggplot(aes(x=pre_pembro, y=pembro_1, color = type)) + 
  geom_point(size = 1, position = position_jitter(width = 0.05, height = 0.05)) + 
  geom_abline() + 
  scale_color_manual(values = c("#7fbf7b","#e0e0e0","#a1cbe2", "#b3b3b3"), breaks  = c("Expanded", "Stable", "Contracted", "Absent")) +
  theme_classic() +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000), limits = c(0, 2000)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000), limits = c(0, 2000)) +
  xlab("Pre Pembro Clone Count") + ylab("Pembro 1 Clone Count") + labs(color=NULL) + theme(legend.position = c(0.82, 0.20), 
                                                                                           legend.title=element_blank(), 
                                                                                           legend.box.margin = margin(0, 0, 0, 0), 
                                                                                           legend.text = element_text(margin = margin(r = 0, unit = "pt"), size = 8), 
                                                                                           legend.background = element_blank(),
                                                                                           legend.box.background = element_blank(),
                                                                                           legend.key.height = unit(3, "mm")) +
  guides(color = guide_legend(override.aes = list(size=2)))
ggsave("analysis/figs/bulk_expansion_scatter_pembro1_vs_pre.eps", width = 3.1, height = 3)

bulktcr_cdr3_counts %>% 
  mutate(expansion = pembro_2/pre_pembro) %>%
  mutate(type = ifelse(expansion > 2, "Expanded", ifelse(expansion < 0.5, "Contracted", "Stable"))) %>%
  mutate(type = ifelse(is.na(type), "Absent", type)) %>%
  select(CDR3aa, pembro_2, pre_pembro, expansion, type) %>%
  ggplot(aes(x=pre_pembro, y=pembro_2, color = type)) + 
  geom_point(size = 1, position = position_jitter(width = 0.05, height = 0.05)) + 
  geom_abline() + 
  scale_color_manual(values = c("#7fbf7b","#e0e0e0","#a1cbe2", "#b3b3b3"), breaks  = c("Expanded", "Stable", "Contracted", "Absent")) +
  theme_classic() +
  scale_x_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000), limits = c(0, 2000)) + 
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), breaks = c(0, 1, 10, 100, 1000), limits = c(0, 2000)) +
  xlab("Pre Pembro Clone Count") + ylab("Pembro 2 Clone Count") + labs(color=NULL) + theme(legend.position = c(0.82, 0.20), 
                                                                                           legend.title=element_blank(), 
                                                                                           legend.box.margin = margin(0, 0, 0, 0), 
                                                                                           legend.text = element_text(margin = margin(r = 0, unit = "pt"), size = 8), 
                                                                                           legend.background = element_blank(),
                                                                                           legend.box.background = element_blank(),
                                                                                           legend.key.height = unit(3, "mm"))+
  guides(color = guide_legend(override.aes = list(size=2)))
ggsave("analysis/figs/bulk_expansion_scatter_pembro2_vs_pre.eps", width = 3.1, height = 3)

#identifying expanded TCRs
jitter <- position_jitter(width = 0.1, height = 0.1)

#clones expanded with pembro in blood
pembro_bloodexpanded_all_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2)

pembro_bloodexpanded_existing_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2 & is.finite(expansion))

pembro_bloodexpanded_emerged_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(!is.finite(expansion)) %>% na.omit()

#clones expanded with pembro in blood, present in kidney
pembro_bloodexpanded_all_kidneypresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2) %>% filter(kidney > 0)

pembro_bloodexpanded_existing_kidneypresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2 & is.finite(expansion)) %>% filter(kidney > 0)

pembro_bloodexpanded_emerged_kidneypresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(!is.finite(expansion)) %>% na.omit() %>% filter(kidney > 0)

#clones expanded with pembro in blood, present in kidney, absent in tumor
pembro_bloodexpanded_all_kidneypresent_tumorabsent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2) %>% filter(kidney > 0)%>% filter(ln == 0 & skin == 0)

pembro_bloodexpanded_existing_kidneypresent_tumorabsent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2 & is.finite(expansion)) %>% filter(kidney > 0) %>% filter(ln == 0 & skin == 0)

pembro_bloodexpanded_emerged_kidneypresent_tumorabsent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(!is.finite(expansion)) %>% na.omit() %>% filter(kidney > 0) %>% filter(ln == 0 & skin == 0)

#clones expanded with pembro in blood, present in tumor
pembro_bloodexpanded_all_lnpresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2) %>% filter(ln > 0)

pembro_bloodexpanded_existing_lnpresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion > 2 & is.finite(expansion)) %>% filter(ln > 0)

pembro_bloodexpanded_emerged_lnpresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(!is.finite(expansion)) %>% na.omit() %>% filter(ln > 0)

#clones expanded across tumors (only pre-existing) after pembro 
pembro_tumorexpanded_and_preexisting_all_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = skin/ln) %>% filter(expansion > 2 & is.finite(expansion)) 

#clones expanded across tumors after pembro (only pre-existing) and present in kidney
pembro_tumorexpanded_and_preexisting_kidneypresent_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = skin/ln) %>% filter(expansion > 2 & is.finite(expansion)) %>% filter(kidney > 0)

#clones expanded across tumors (only pre-existing) and expanded in blood after pembro 
pembro_tumorexpanded_and_preexisting_bloodexpanded_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion_tumor = skin/ln) %>% filter(expansion_tumor > 2 & is.finite(expansion_tumor)) %>%
  mutate(expansion_blood = pembro_1/pre_pembro) %>% filter(expansion_blood > 2)

#clones kidney only (tumor absent), blood any
pembro_kidneypresent_tumorabsentboth_clones <- bulktcr_cdr3_freqs %>% filter(kidney > 0) %>% filter(ln == 0 & skin == 0)

#clones tumor only (kidney absent), blood any
pembro_tumorpresentboth_kidneyabsent_clones <- bulktcr_cdr3_freqs %>% filter(ln > 0 & skin > 0) %>% filter(kidney == 0)

#clones ln tumor only (kidney absent), blood any 
pembro_lntumorpresent_kidneyabsent_clones <- bulktcr_cdr3_freqs %>% filter(ln > 0) %>% filter(kidney == 0)

#clones skin tumor only (kidney absent), blood any 
pembro_skintumorpresent_kidneyabsent_clones <- bulktcr_cdr3_freqs %>% filter(skin > 0) %>% filter(kidney == 0)

#clones present in kidney and both blood, blood any 
pembro_tumorpresentboth_kidneypresent_clones <- bulktcr_cdr3_freqs %>% filter(ln > 0 & skin > 0) %>% filter(kidney > 0)

#clones stable with pembro in blood
pembro_bloodstable_all_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion < 2 & expansion > 0.5)

#clones contract with pembro in blood
pembro_bloodcontract_all_clones <- bulktcr_cdr3_freqs %>% 
  mutate(expansion = pembro_1/pre_pembro) %>% filter(expansion <= 0.5)



# Single-cell TCR and transcriptome analysis of MLR -----------------------
#load 10x data
MLR_sample1_data <- Read10X_h5("10x_MLR/data/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
MLR_sample1 <- CreateSeuratObject(MLR_sample1_data, project = "10xMLR")

MLR_sample1_TCR_data <- read.csv("10x_MLR/data/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
TCR_list <- list(MLR_sample1_TCR_data)
MLR_sample1_TCR <- combineTCR(TCR_list, samples = "Patient", ID = "1", cells = "T-AB", filterMulti = T)
MLR_sample1_TCR$Patient_1$barcode <- substring(MLR_sample1_TCR$Patient_1$barcode, 11)

MLR_sample1 <- combineExpression(MLR_sample1_TCR$Patient_1, MLR_sample1, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=3, Medium=10, Large=30, Hyperexpanded=100))

MLR_sample1$cloneType <- as.character(MLR_sample1$cloneType)
MLR_sample1$cloneType <- replace_na(MLR_sample1$cloneType, "No TCR")
Idents(MLR_sample1) <- MLR_sample1$cloneType

slot(MLR_sample1, "meta.data")$cloneType <- factor(slot(MLR_sample1, "meta.data")$cloneType, 
                                                   levels = c("Hyperexpanded (30 < X <= 100)", 
                                                              "Large (10 < X <= 30)", 
                                                              "Medium (3 < X <= 10)", 
                                                              "Small (1 < X <= 3)", 
                                                              "Single (0 < X <= 1)",
                                                              "No TCR"))

current.clone.names <- c("Hyperexpanded (30 < X <= 100)", "Large (10 < X <= 30)", "Medium (3 < X <= 10)", "Small (1 < X <= 3)", "Single (0 < X <= 1)", "No TCR")
new.clone.names <- c("Hyperexpanded (30 < X)", "Large (10 < X <= 30)", "Medium (3 < X <= 10)", "Small (1 < X <= 3)", "Single (0 < X <= 1)", "No TCR information")
MLR_sample1@active.ident <- plyr::mapvalues(x = MLR_sample1@active.ident, from = current.clone.names, to = new.clone.names)
MLR_sample1$cloneType <- MLR_sample1@active.ident

slot(MLR_sample1, "meta.data")$cloneType <- factor(slot(MLR_sample1, "meta.data")$cloneType, 
                                                   levels = c("Hyperexpanded (30 < X)", 
                                                              "Large (10 < X <= 30)", 
                                                              "Medium (3 < X <= 10)", 
                                                              "Small (1 < X <= 3)", 
                                                              "Single (0 < X <= 1)", 
                                                              "No TCR information"))

MLR_sample1$clonal <-ifelse(MLR_sample1$Frequency > 1, "Yes", "No")

#QC
MLR_sample1[["percent.mito"]] <- PercentageFeatureSet(MLR_sample1, pattern = "^MT-")
MLR_sample1[["percent.ribo"]] <- PercentageFeatureSet(MLR_sample1, pattern = "^RP[SL]")

Idents(MLR_sample1) <- MLR_sample1$orig.ident
VlnPlot(MLR_sample1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) +
  theme(legend.position = "none") &
  scale_x_discrete(labels = "MLR") &
  xlab(NULL)
ggsave("analysis/figs/sc_QCmetrics_preSubset.eps", width = 6, height = 3)

MLR_preQC_cellnum <- MLR_sample1@meta.data %>% group_by(orig.ident) %>% count()
MLR_preQC_cellnum <- MLR_preQC_cellnum$n #1587

MLR_preQC_VDJcellnum <- MLR_sample1@meta.data %>% filter(cloneType != "No TCR information") %>% group_by(orig.ident) %>% count()
MLR_preQC_VDJcellnum <- MLR_preQC_VDJcellnum$n #1327

MLR_sample1_subset <- subset(MLR_sample1, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mito < 25)

MLR_postQC_cellnum <- MLR_sample1_subset@meta.data %>% group_by(orig.ident) %>% count()
MLR_postQC_cellnum <- MLR_postQC_cellnum$n #1505

MLR_postQC_VDJcellnum <- MLR_sample1_subset@meta.data %>% filter(cloneType != "No TCR information") %>% group_by(orig.ident) %>% count()
MLR_postQC_VDJcellnum <- MLR_postQC_VDJcellnum$n #1311

MLR_postQC_clonalNonclonalnums <- MLR_sample1_subset@meta.data %>% group_by(clonal) %>% count()

MLR_cellnums <-  data.frame(QC = c("preQC", "postQC"), Nums = c(MLR_preQC_cellnum, MLR_postQC_cellnum))

ggplot(MLR_cellnums, aes(x = QC, y = Nums, label = Nums)) + 
  geom_col(fill = "#C0C0C0") + 
  geom_text(position = position_dodge(width = 1), vjust = 2, size = 4) +
  scale_x_discrete(limits = c("preQC", "postQC"), labels = c("Pre QC", "Post QC")) + 
  xlab(NULL) + ylab("Number of Cells") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("analysis/figs/sc_MLR_QC_barplot.eps", width = 2, height = 2.5)


#dimensionality reduction
MLR_sample1_subset <- NormalizeData(MLR_sample1_subset, normalization.method = "LogNormalize", scale.factor = 10000)
MLR_sample1_subset <- FindVariableFeatures(MLR_sample1_subset, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(MLR_sample1_subset), 10)

all.genes <- rownames(MLR_sample1_subset)
MLR_sample1_subset <- ScaleData(MLR_sample1_subset, features = all.genes, vars.to.regress = "percent.mito")
MLR_sample1_subset <- RunPCA(MLR_sample1_subset, features = VariableFeatures(object = MLR_sample1_subset))

MLR_sample1_subset <- FindNeighbors(MLR_sample1_subset, dims = 1:30)
MLR_sample1_subset <- FindClusters(MLR_sample1_subset, resolution = 0.7)

MLR_sample1_subset <- RunUMAP(MLR_sample1_subset, dims = 1:40)

MLR_sample1_subset$seurat_clusters <- paste0("C", MLR_sample1_subset$seurat_clusters)
MLR_sample1_subset$orig_cluster_order <- MLR_sample1_subset$seurat_clusters
slot(MLR_sample1_subset, "meta.data")$orig_cluster_order <- factor(slot(MLR_sample1_subset, "meta.data")$orig_cluster_order, 
                                                                   levels = c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"))

Idents(MLR_sample1_subset) <- MLR_sample1_subset$seurat_clusters
current.cluster.names <- c("C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8")
new.cluster.names <- c("C4", "C6", "C3", "C1", "C5", "C7", "C2", "C9", "C8")
MLR_sample1_subset@active.ident <- plyr::mapvalues(x = MLR_sample1_subset@active.ident, from = current.cluster.names, to = new.cluster.names)
MLR_sample1_subset$seurat_clusters <- MLR_sample1_subset@active.ident

slot(MLR_sample1_subset, "meta.data")$seurat_clusters <- factor(slot(MLR_sample1_subset, "meta.data")$seurat_clusters, 
                                                                levels = c( "C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"))
Idents(MLR_sample1_subset) <- MLR_sample1_subset$seurat_clusters

#singlecellcolors <- c("#e76f51", "#f4a261", "#152418", "#2a4830", "#3f6d48", "#549160", "#69b578", "#78bc86", "#96cba1")

singlecellcolors_set2 <- c(colorRampPalette(brewer.pal(9, "Set2"))(9))
singlecellcolors_set3 <- c(colorRampPalette(brewer.pal(9, "Set3"))(9))
show_col(c(singlecellcolors_set2))
show_col(c(singlecellcolors_set3))

singlecellcolors <- c(singlecellcolors_set3[1], singlecellcolors_set2[4],singlecellcolors_set3[c(9,4:8)], singlecellcolors_set2[7])

DimPlot(MLR_sample1_subset, group.by = "seurat_clusters") + 
  scale_color_manual(values = singlecellcolors) +
  ggtitle(NULL) + xlab("UMAP 1") + ylab("UMAP 2")
ggsave("analysis/figs/sc_UMAP_MLR_byCluster.eps", width = 5, height = 4)

MLRclusterUMAP <- DimPlot(MLR_sample1_subset, group.by = "seurat_clusters") + 
  scale_color_manual(values = singlecellcolors) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

MLRclusterUMAPlabeled <- DimPlot(MLR_sample1_subset, group.by = "seurat_clusters", label = T, repel = T, label.box = T) + 
  scale_color_manual(values = singlecellcolors) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

MLR_10xcluster_breakdown <- as.data.frame(table(MLR_sample1_subset$seurat_clusters))
ggplot(MLR_10xcluster_breakdown, aes(x = 2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  guides(fill = guide_legend(reverse = TRUE)) +
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = singlecellcolors) +
  theme_void() +
  xlim(0.5, 2.5)  
ggsave("analysis/figs/sc_MLR_piechart_byCluster.eps", width = 4, height = 3)

MLR10x_markers_tcells <- FindAllMarkers(object=MLR_sample1_subset, only.pos = FALSE, logfc.threshold = 0.25, min.pct = 0.20)
write.csv(MLR10x_markers_tcells, "analysis/data/sc_MLR_RNA_cluster_markers.csv")

MLR10x_markers_tcells_pos <- FindAllMarkers(object=MLR_sample1_subset, only.pos = T, logfc.threshold = 0.25, min.pct = 0.20)

DotPlot(MLR_sample1_subset, features = c("CD3E", "CD8A", "CD4", "ITGAE", "MKI67", "PCNA", "IL2RA", "GZMB", "GNLY", "PRF1", "NKG7", "KLRB1",  "CTLA4", "PDCD1", "IL7R", "CD69", "SELL", "LEF1", "FOXP3", "IFI44L", "IFIT1", "MALAT1"), group.by = "seurat_clusters") + 
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))  + RotatedAxis() +
  theme(legend.title=element_text(size=11),legend.text=element_text(size=10)) +
  scale_y_discrete(limits = rev(levels(MLR_sample1_subset))) +
  xlab(NULL) +
  ylab("Cluster") 
ggsave("analysis/figs/sc_dotplot_Clusters.eps", width = 8, height = 3.5)

StackedVlnPlot(MLR_sample1_subset, features = c("CD3E", "CD8A", "CD4", "ITGAE", "MKI67", "PCNA", 
                                                "IL2RA", "IFNG", "GZMB", "PRF1", "NKG7",  "GNLY", "KLRB1",
                                                "CTLA4", "PDCD1", "IL7R", "CD69", "SELL", "LEF1", 
                                                "FOXP3", "IFI44L", "IFIT1", "MALAT1"), cols = singlecellcolors)
ggsave("analysis/figs/sc_stackedVln_Clusters.eps", width = 4, height = 6)

top10mlr_deg_pos <- MLR10x_markers_tcells_pos %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n=10)
DoHeatmap(MLR_sample1_subset, features = top10mlr_deg_pos$gene, group.colors = singlecellcolors, size = 4) + 
  scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.y = element_text(size = 6))
ggsave("analysis/figs/sc_heatmap_Clusters.eps", width = 5, height = 7)

MLR_CD4 <- FeaturePlot(MLR_sample1_subset, features = c("CD4")) +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds"))[2:9], guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), breaks = c(0,1,2,3), limits = c(0,3)) +
  ggtitle("CD4") + 
  theme(legend.position = "bottom", legend.justification = "center",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

MLR_CD8 <- FeaturePlot(MLR_sample1_subset, features = c("CD8A")) +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds"))[2:9], guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")) +
  ggtitle("CD8A") + 
  theme(legend.position = "bottom", legend.justification = "center",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

MLR_CD4 | MLR_CD8
ggsave("analysis/figs/sc_UMAP_CD4_CD8.eps", width = 8, height = 4)

MLR_CD4
ggsave("analysis/figs/sc_UMAP_CD4.eps", width = 4, height = 4)

MLR_CD8
ggsave("analysis/figs/sc_UMAP_CD8.eps", width = 4, height = 4)



#proliferation GO signature
prolif_genelist <- read.csv("10x_MLR/proliferation_geneset.txt", header = F)
prolif_genelist <- prolif_genelist$V1
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(prolif_genelist), name = "Proliferation_sig", ctrl = 1000)

VlnPlot(MLR_sample1_subset, features = "Proliferation_sig1", group.by = "seurat_clusters") + 
  scale_fill_manual(values = singlecellcolors) +
  ggtitle(NULL) + xlab(NULL) +
  scale_y_continuous(breaks=c(0,0.25), labels=c("Min","Max"), limits = c(0,0.25)) +
  theme(legend.position = "none") + ylab("Proliferation Signature") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) 
ggsave("analysis/figs/sc_Violin_proliferation_MLR_byCluster.eps", width = 5, height = 3)

MLR_sample1_subset@meta.data %>%
  ggplot(aes(x = seurat_clusters, y = Proliferation_sig1, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("Proliferation Signature") + 
  scale_fill_manual(values = singlecellcolors) +
  scale_y_continuous(breaks=c(0,0.25), labels=c("Min","Max"), limits = c(0,0.25)) 
ggsave("analysis/figs/sc_boxplot_proliferation_MLR_byCluster.eps", width = 5, height = 3)

FeaturePlot(MLR_sample1_subset, features = "Proliferation_sig1") +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds")), 
                        breaks=c(0,0.25), labels=c("Min","Max"), na.value = "white",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits = c(0,0.25)) + ggtitle("Proliferation Signature") +
  theme(legend.position = "bottom", legend.justification = "center",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())
ggsave("analysis/figs/sc_UMAP_proliferation_MLR.eps", width = 4, height = 4)

#activation GO signature
active_genelist <- read.csv("10x_MLR/activation_geneset.csv", header = F)
active_genelist <- active_genelist$V1
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(active_genelist), name = "Activation_sig", ctrl = 1000)

VlnPlot(MLR_sample1_subset, features = "Activation_sig1", group.by = "seurat_clusters") + 
  scale_fill_manual(values = singlecellcolors) +
  ggtitle(NULL) + xlab(NULL) +
  scale_y_continuous(breaks=c(0,0.30), labels=c("Min","Max"), limits = c(0,0.30)) 
ggsave("analysis/figs/sc_Violin_activation_MLR_byCluster.eps", width = 5, height = 3)

FeaturePlot(MLR_sample1_subset, features = "Activation_sig1") +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds")), 
                        breaks=c(0,0.30), labels=c("Min","Max"), na.value = "white",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits = c(0,0.30)) + ggtitle(NULL) + 
  theme(legend.position = "bottom", legend.justification = "center",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())
ggsave("analysis/figs/sc_UMAP_activation_MLR.eps", width = 4, height = 4)

#mito and ribo signatures
FeaturePlot(MLR_sample1_subset, features = "percent.mito") +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds")), 
                        breaks=c(0,25), labels=c("Min","Max"), na.value = "white",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits = c(0,25)) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank()) +
  ggtitle("% Mito Genome Reads") 
ggsave("analysis/figs/sc_UMAP_mito_MLR.eps", width = 4, height = 3.5)

VlnPlot(MLR_sample1_subset, features = "percent.mito", group.by = "seurat_clusters") + 
  scale_fill_manual(values = singlecellcolors) +
  ggtitle(NULL) + xlab(NULL) 
ggsave("analysis/figs/sc_Violin_mito_MLR_byCluster.eps", width = 5, height = 3)

FeaturePlot(MLR_sample1_subset, features = "percent.ribo") +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds")), 
                        breaks=c(0,50), labels=c("Min","Max"), na.value = "white",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits = c(0,50)) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank()) +
  ggtitle("% Ribo Protein Gene Reads") 
ggsave("analysis/figs/sc_UMAP_ribo_MLR.eps", width = 4, height = 3.5)

#Additional signatures (from Pai, Chow, et al- see "10x_MLR/pai_chow_genesets_full.xlsx" for more information)
gene_sigs <- read.csv("10x_MLR/pai_chow_genesets.csv", header = T)
virus_specific_genelist <- gene_sigs$Virus.Specific[-which(gene_sigs$Virus.Specific == "")]
influenza_TIL_genelist <- gene_sigs$Influenza.TIL[-which(gene_sigs$Influenza.TIL == "")]
tumor_specific_genelist <- gene_sigs$Tumor.Specific
tumor_reactive_genelist <- gene_sigs$Tumor.reactivity[-which(gene_sigs$Tumor.reactivity == "")]
prolif_peripheralcd8_genelist <- gene_sigs$Proliferating.peripheral.CD8.[-which(gene_sigs$Proliferating.peripheral.CD8. == "")]
progenitor_genelist <- gene_sigs$Progenitor.score[-which(gene_sigs$Progenitor.score == "")]
cd8_exhaustion_genelist <- gene_sigs$CD8..exhaustion[-which(gene_sigs$CD8..exhaustion == "")]
proliferation_genelist <- gene_sigs$Proliferation[-which(gene_sigs$Proliferation == "")]

MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(virus_specific_genelist), name = "Virus_specific_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(influenza_TIL_genelist), name = "Influenza_TIL_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(tumor_specific_genelist), name = "Tumor_specific_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(tumor_reactive_genelist), name = "Tumor_reactive_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(prolif_peripheralcd8_genelist), name = "Proliferating_peripheralCD8_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(progenitor_genelist), name = "Progenitor_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(cd8_exhaustion_genelist), name = "CD8_exhaustion_sig", ctrl = 1000)
MLR_sample1_subset <- AddModuleScore(MLR_sample1_subset, list(proliferation_genelist), name = "Proliferation_Li_sig", ctrl = 1000)

VlnPlot(MLR_sample1_subset, features = c("Virus_specific_sig1", "Influenza_TIL_sig1", "Tumor_specific_sig1", "Tumor_reactive_sig1", "Proliferating_peripheralCD8_sig1", "Progenitor_sig1", "CD8_exhaustion_sig1", "Proliferation_sig1", "Activation_sig1", "Proliferation_Li_sig1"), group.by = "seurat_clusters") & 
  scale_fill_manual(values = singlecellcolors) & xlab(NULL)
ggsave("analysis/figs/sc_Violin_signatures_byCluster.eps", width = 12, height = 9)

MLR_sample1_subset@meta.data %>%
  filter(seurat_clusters %in% c("C1", "C2")) %>%
  ggplot(aes(x = seurat_clusters, y = Virus_specific_sig1, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("Viral Specific\n (Oliviera et al, 2021)") + 
  scale_fill_manual(values = singlecellcolors) +
  scale_y_continuous(breaks=c(-0.3, 0.3), labels=c("Min","Max"), limits = c(-0.3, 0.3)) 
ggsave("analysis/figs/sc_boxplot_C1C2_viral.eps", width = 2, height = 2.5)

MLR_sample1_subset@meta.data %>%
  filter(seurat_clusters %in% c("C1", "C2")) %>%
  ggplot(aes(x = seurat_clusters, y = Influenza_TIL_sig1, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("Influenza Specific\n (Caushi et al, 2021)") + 
  scale_fill_manual(values = singlecellcolors) +
  scale_y_continuous(breaks=c(-0.25, 0.55), labels=c("Min","Max"), limits = c(-0.25, 0.55)) 
ggsave("analysis/figs/sc_boxplot_C1C2_influenza.eps", width = 2, height = 2.5)

VlnPlot(MLR_sample1_subset, features = "Proliferation_Li_sig1", group.by = "seurat_clusters") + 
  scale_fill_manual(values = singlecellcolors) +
  ggtitle(NULL) + xlab(NULL) +
  scale_y_continuous(breaks=c(-0.25, 1.35), labels=c("Min","Max"), limits = c(-0.25, 1.35)) +
  theme(legend.position = "none") + ylab("Proliferation Signature") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1)) 
ggsave("analysis/figs/sc_Violin_proliferation_LiAmit_MLR_byCluster.eps", width = 5, height = 3)

MLR_sample1_subset@meta.data %>%
  ggplot(aes(x = seurat_clusters, y = Proliferation_Li_sig1, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("Proliferation Signature") + 
  scale_fill_manual(values = singlecellcolors) +
  scale_y_continuous(breaks=c(-0.25, 1.35), labels=c("Min","Max"), limits = c(-0.25, 1.35)) 
ggsave("analysis/figs/sc_boxplot_proliferation_LiAmit_MLR_byCluster.eps", width = 5, height = 3)

ProlifLiUMAP <- FeaturePlot(MLR_sample1_subset, features = "Proliferation_Li_sig1") +
  scale_color_gradientn(colours = (brewer.pal(n = 9, name = "Reds")), 
                        breaks=c(-0.25, 1.35), labels=c("Min","Max"), na.value = "white",
                        guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"), limits = c(-0.25, 1.35)) + ggtitle("Proliferation Signature") +
  theme(legend.position = "bottom", legend.justification = "center",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

ProlifLiUMAP
ggsave("analysis/figs/sc_UMAP_proliferation_LiAmit_MLR.eps", width = 4, height = 4)

#cell cycle analysis
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
MLR_sample1_subset <- CellCycleScoring(MLR_sample1_subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
slot(MLR_sample1_subset, "meta.data")$Phase <- factor(slot(MLR_sample1_subset, "meta.data")$Phase, 
                                                      levels = c("G1", "S", "G2M"))

MLRcellcycleUMAP <- DimPlot(MLR_sample1_subset, group.by = "Phase") +
  scale_color_manual(values = c("#f9dbbd", "#d81159","#0496ff"), labels = c("G1", "S", "G2/M")) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

MLRcellcycleUMAP 
ggsave("analysis/figs/sc_UMAP_cellcycle_phase_MLR.eps", width = 4.5, height = 3.5)

sscore <- FeaturePlot(MLR_sample1_subset, features = c("S.Score")) + ggtitle("S Score") + xlab("UMAP 1") + ylab("UMAP 2")
g2mscore <- FeaturePlot(MLR_sample1_subset, features = c("G2M.Score")) + ggtitle("G2/M Score") + xlab("UMAP 1") + ylab("UMAP 2")

sscore | g2mscore
ggsave("analysis/figs/sc_UMAP_cellcycle_scores_MLR.eps", width = 7, height = 3.5)

VlnPlot(MLR_sample1_subset, features = "S.Score", group.by = "seurat_clusters") + 
  scale_fill_manual(values = singlecellcolors) +
  ggtitle("S Phase Score") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "none") 
ggsave("analysis/figs/sc_Violin_Sphase_MLR_byCluster.eps", width = 5.5, height = 3)

MLR_sample1_subset@meta.data %>%
  ggplot(aes(x = seurat_clusters, y = S.Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("S Phase Score") + 
  scale_fill_manual(values = singlecellcolors) +
  scale_y_continuous(breaks=c(-0.37, 0.69), labels=c("Min","Max"), limits = c(-0.37, 0.69)) 
ggsave("analysis/figs/sc_boxplot_Sphase_MLR_byCluster.eps", width = 5, height = 3)


VlnPlot(MLR_sample1_subset, features = "G2M.Score", group.by = "seurat_clusters") + 
  scale_fill_manual(values = singlecellcolors) +
  ggtitle("G2/M Phase Score") + xlab(NULL) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        legend.position = "none") 
ggsave("analysis/figs/sc_Violin_G2Mphase_MLR_byCluster.eps", width = 5.5, height = 3)

MLR_sample1_subset@meta.data %>%
  ggplot(aes(x = seurat_clusters, y = G2M.Score, fill = seurat_clusters)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("G2/M Phase Score") + 
  scale_fill_manual(values = singlecellcolors) +
  scale_y_continuous(breaks=c(-0.34, 1.04), labels=c("Min","Max"), limits = c(-0.34, 1.04)) 
ggsave("analysis/figs/sc_boxplot_G2Mphase_MLR_byCluster.eps", width = 5, height = 3)

#pathway analysis
MLR10x_markers_tcells_DEGsforpathway <- FindAllMarkers(object=MLR_sample1_subset, only.pos = TRUE, logfc.threshold = 0.50, min.pct = 0.30)

dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015", "KEGG_2021_Human")

enriched_C1 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C1", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C1 <- lapply(enriched_C1, function(x) 
  cbind(x, Cluster = "C1"))
enriched_C2 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C2", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C2 <- lapply(enriched_C2, function(x) 
  cbind(x, Cluster = "C2"))
enriched_C3 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C3", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C3 <- lapply(enriched_C3, function(x) 
  cbind(x, Cluster = "C3"))
enriched_C4 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C4", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C4 <- lapply(enriched_C4, function(x) 
  cbind(x, Cluster = "C4"))
enriched_C5 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C5", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C5 <- lapply(enriched_C5, function(x) 
  cbind(x, Cluster = "C5"))
enriched_C6 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C6", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C6 <- lapply(enriched_C6, function(x) 
  cbind(x, Cluster = "C6"))
enriched_C7 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C7", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C7 <- lapply(enriched_C7, function(x) 
  cbind(x, Cluster = "C7"))
enriched_C8 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C8", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C8 <- lapply(enriched_C8, function(x) 
  cbind(x, Cluster = "C8"))
enriched_C9 <- enrichr(c(MLR10x_markers_tcells_DEGsforpathway$gene[grep("C9", MLR10x_markers_tcells_DEGsforpathway$cluster)]), dbs)
enriched_C9 <- lapply(enriched_C9, function(x) 
  cbind(x, Cluster = "C9"))

GObiologicalprocess_top5 <- rbind(enriched_C1$GO_Biological_Process_2015[1:5,], 
                                  enriched_C2$GO_Biological_Process_2015[1:5,], enriched_C3$GO_Biological_Process_2015[1:5,],
                                  enriched_C4$GO_Biological_Process_2015[1:5,], enriched_C5$GO_Biological_Process_2015[1:5,],
                                  enriched_C6$GO_Biological_Process_2015[1:5,], enriched_C7$GO_Biological_Process_2015[1:5,],
                                  enriched_C8$GO_Biological_Process_2015[1:5,], enriched_C9$GO_Biological_Process_2015[1:5,])
GObiologicalprocess_top5$TermCluster <- paste0(GObiologicalprocess_top5$Term, GObiologicalprocess_top5$Cluster)
GObiologicalprocess_top5$TermCluster <- factor(GObiologicalprocess_top5$TermCluster,levels = c(GObiologicalprocess_top5$TermCluster))

ggplot(GObiologicalprocess_top5, aes(x=TermCluster, y=-log(Adjusted.P.value))) +
  geom_segment(aes(x=TermCluster, xend=TermCluster, y=0, yend=-log(Adjusted.P.value), color = Cluster)) +
  geom_point( aes(color=Cluster), size=4) +
  scale_color_manual(values = singlecellcolors) + 
  scale_x_discrete(labels = rev(GObiologicalprocess_top5$Term), limits = rev(levels(GObiologicalprocess_top5$TermCluster))) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab(NULL) + ylab("-log(Adjusted P Value)") +
  theme_classic() +
  coord_flip()
ggsave("analysis/figs/sc_UMAP_pathway_biologicalProcess_byCluster_LFC.50_pct30.eps", width = 12, height = 7)


GOkegg_top5 <- rbind(enriched_C1$KEGG_2021_Human[1:5,], enriched_C2$KEGG_2021_Human[1:5,], 
                     enriched_C3$KEGG_2021_Human[1:5,], enriched_C4$KEGG_2021_Human[1:5,],
                     enriched_C5$KEGG_2021_Human[1:5,], enriched_C6$KEGG_2021_Human[1:5,],
                     enriched_C7$KEGG_2021_Human[1:5,], enriched_C8$KEGG_2021_Human[1:5,],
                     enriched_C9$KEGG_2021_Human[1:5,]) %>%
  filter(!is.na(Term))
GOkegg_top5$TermCluster <- paste0(GOkegg_top5$Term, GOkegg_top5$Cluster)
GOkegg_top5$TermCluster <- factor(GOkegg_top5$TermCluster,levels = c(GOkegg_top5$TermCluster))

ggplot(GOkegg_top5, aes(x=TermCluster, y=-log(Adjusted.P.value))) +
  geom_segment(aes(x=TermCluster, xend=TermCluster, y=0, yend=-log(Adjusted.P.value), color = Cluster)) +
  geom_point( aes(color=Cluster), size=4) +
  scale_color_manual(values = singlecellcolors) + 
  scale_x_discrete(labels = rev(GOkegg_top5$Term), limits = rev(levels(GOkegg_top5$TermCluster))) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab(NULL) + ylab("-log(Adjusted P Value)") +
  theme_classic() +
  coord_flip()
ggsave("analysis/figs/sc_UMAP_pathway_KEGG_byCluster_LFC.50_pct30.eps", width = 12, height = 7)

#slingshot trajectory inference
#on selected clusters
MLR_sample1_subset_slingshotremove <- subset(MLR_sample1_subset, idents = c("C7","C2","C9","C8"), invert = TRUE)

kidney_MLR_TI_subset_start4 <- slingshot(Embeddings(MLR_sample1_subset_slingshotremove, "umap"), clusterLabels = MLR_sample1_subset_slingshotremove$seurat_clusters, 
                                         start.clus = "C4", stretch = 0, omega = TRUE)
kidney_MLR_TI_subset_start6 <- slingshot(Embeddings(MLR_sample1_subset_slingshotremove, "umap"), clusterLabels = MLR_sample1_subset_slingshotremove$seurat_clusters, 
                                         start.clus = "C6", stretch = 0, omega = TRUE)
kidney_MLR_TI_subset_unbiasedstart <- slingshot(Embeddings(MLR_sample1_subset_slingshotremove, "umap"), clusterLabels = MLR_sample1_subset_slingshotremove$seurat_clusters, 
                                                stretch = 0, omega = TRUE)

setEPS()
postscript("analysis/figs/sc_slingshot_subset_start4.eps", width = 5, height = 4)
par(mar=c(5,4.5,1,1))
plot(Embeddings(MLR_sample1_subset, "umap"), col = singlecellcolors[Idents(MLR_sample1_subset)], pch = 16, cex = 0.65, axes = F)
box(bty="l")
axis(1)
axis(2, las = 2)
lines(kidney_MLR_TI_subset_start4, lwd = 3, col = 'black')
dev.off()

setEPS()
postscript("analysis/figs/sc_slingshot_subset_start6.eps", width = 5, height = 4)
par(mar=c(5,4.5,1,1))
plot(Embeddings(MLR_sample1_subset, "umap"), col = singlecellcolors[Idents(MLR_sample1_subset)], pch = 16, cex = 0.65, axes = F)
box(bty="l")
axis(1)
axis(2, las = 2)
lines(kidney_MLR_TI_subset_start6, lwd = 3, col = 'black')
dev.off()

setEPS()
postscript("analysis/figs/sc_slingshot_subset_unbiasedstart.eps", width = 5, height = 4)
par(mar=c(5,4.5,1,1))
plot(Embeddings(MLR_sample1_subset, "umap"), col = singlecellcolors[Idents(MLR_sample1_subset)], pch = 16, cex = 0.65, axes = F)
box(bty="l")
axis(1)
axis(2, las = 2)
lines(kidney_MLR_TI_subset_unbiasedstart, lwd = 3, col = 'black')
dev.off()

#DEG C1 vs C2
C1vsC2markers_wilcox <- FindMarkers(MLR_sample1_subset, ident.1 = "C1", ident.2 = "C2", min.pct = 0.5, test.use = "wilcox")
EnhancedVolcano(C1vsC2markers_wilcox, lab = as.character(rownames(C1vsC2markers_wilcox)), x = 'avg_log2FC', y = 'p_val_adj', title = "C1 vs C2", subtitle = NULL, legendPosition = 'right', legendLabels = NULL, legendIconSize = NULL, legendLabSize = NULL)
ggsave("analysis/figs/sc_Volcano_C1vsC2_wilcox.pdf", width = 24, height = 14)

#clonality
cloneType_byCluster_table <- table(MLR_sample1_subset$cloneType, MLR_sample1_subset$seurat_clusters)

MLRclonalityUMAP <- DimPlot(MLR_sample1_subset, group.by = "cloneType") + 
  scale_color_brewer(palette = "Blues", direction = -1) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

MLRclonalityUMAP
ggsave("analysis/figs/sc_UMAP_MLR_byClonality.eps", width = 6.5, height = 4)

MLR_sample1_subset_cluster_clones <- expression2List(MLR_sample1_subset, group = "seurat_clusters")

MLR_clust_morisita <- clonalOverlap(MLR_sample1_subset_cluster_clones, cloneCall = "aa", method = "morisita", exportTable = T)
MLR_clust_morisita$names <- NULL
MLR_clust_morisita <- as.matrix(MLR_clust_morisita)
MLR_clust_morisita <- Matrix::forceSymmetric(MLR_clust_morisita,uplo="U")
MLR_clust_morisita <- as.matrix(MLR_clust_morisita)
class(MLR_clust_morisita) <- c("immunr_ov_matrix", "matrix", "array")
vis(MLR_clust_morisita) + 
  scale_y_discrete(labels=rev(levels(MLR_sample1_subset$seurat_clusters))) + 
  scale_x_discrete(labels=levels(MLR_sample1_subset$seurat_clusters)) + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
  guides(fill = guide_colorbar(title="Morisita \nIndex")) +
  ggtitle(NULL) + xlab(NULL) + ylab(NULL) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
        axis.text=element_text(size=12)) 
ggsave("analysis/figs/sc_cluster_overlap_morisita.eps", width = 8, height = 8)

MLR_sample1_subset@meta.data %>%
  group_by(seurat_clusters, cloneType) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=Percent, fill=cloneType)) +
  geom_col(color = "black", size = 0.5) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        legend.text = element_text(size=12))
ggsave("analysis/figs/sc_MLR_ClusterClonality_byCluster.eps", width = 7.5, height = 4)

c1clones <- MLR_sample1_subset@meta.data %>% filter(seurat_clusters == "C1") %>% select(CTaa)
c2clones <- MLR_sample1_subset@meta.data %>% filter(seurat_clusters == "C2") %>% select(CTaa)

c1clones <-c1clones$CTaa
c2clones <-c2clones$CTaa

c1clones <- na.omit(c1clones)
c2clones <- na.omit(c2clones)


#plotting same dim UMAPs
MLR_UMAPs <- align_plots(MLRclusterUMAP, MLRclusterUMAPlabeled, MLR_CD4, MLR_CD8, ProlifLiUMAP, MLRcellcycleUMAP, MLRclonalityUMAP, align="hv", axis="tblr")
p1x <- ggdraw(MLR_UMAPs[[1]])
p2x <- ggdraw(MLR_UMAPs[[2]])
p3x <- ggdraw(MLR_UMAPs[[3]])
p4x <- ggdraw(MLR_UMAPs[[4]])
p5x <- ggdraw(MLR_UMAPs[[5]])
p6x <- ggdraw(MLR_UMAPs[[6]])
p7x <- ggdraw(MLR_UMAPs[[7]])

save_plot("analysis/figs/sc_UMAP_MLR_byCluster_aligned.eps", p1x, base_width = 5.5, base_height = 4)
save_plot("analysis/figs/sc_UMAP_MLR_byCluster_labeled_aligned.eps", p2x, base_width = 5.5, base_height = 4)
save_plot("analysis/figs/sc_UMAP_CD4_aligned.eps", p3x, base_width = 5.5, base_height = 4)
save_plot("analysis/figs/sc_UMAP_CD8_aligned.eps", p4x, base_width = 5.5, base_height = 4)
save_plot("analysis/figs/sc_UMAP_proliferation_LiAmit_MLR_aligned.eps", p5x, base_width = 5.5, base_height = 4)
save_plot("analysis/figs/sc_UMAP_cellcycle_phase_MLR_aligned.eps", p6x, base_width = 5.5, base_height = 4)
save_plot("analysis/figs/sc_UMAP_MLR_byClonality_aligned.eps", p7x, base_width = 5.5, base_height = 4)


#examining location of top MLR clones
MLR_clonotypes <- as_tibble(MLR_sample1_subset@meta.data)
MLR_clonotypes <- MLR_clonotypes %>% 
  mutate(CTaa_split = CTaa) %>%
  dplyr::count(CTaa, CTgene, CTaa_split, sort = T) %>%
  drop_na() %>%
  mutate(percent = n/sum(n)*100) %>%
  separate(CTaa_split, into = c("alpha", "beta"), sep="_(?=[^_]+$)")
MLR_clonotypes <- MLR_clonotypes[-c(5),] #removing a double alpha
write.csv(MLR_clonotypes, "analysis/data/sc_MLR_topClonotypes.csv")

Idents(MLR_sample1_subset) <- MLR_sample1_subset$CTaa

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[1])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 1", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 1") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_01.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[2])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 2", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 2") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_02.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[3])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 3", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 3") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_03.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[4])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 4", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 4") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_04.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[5])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 5", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 5") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_05.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[6])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 6", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 6") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_06.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[7])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 7", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 7") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_07.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[8])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 8", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 8") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_08.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[9])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 9", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 9") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_09.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[10])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 10", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 10") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_10.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[11])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 11", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 11") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_11.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[12])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 12", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 12") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_12.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[13])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 13", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 13") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_13.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[14])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 14", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 14") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_14.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[15])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 15", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 15") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_15.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[16])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 16", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 16") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_16.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[17])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 17", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 17") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_17.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[18])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 18", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 18") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_18.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[19])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 19", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 19") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_19.eps", width = 3, height = 3)

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_clonotypes$CTaa[20])) +
  scale_color_manual(name = NULL, limits=c("Group_1", "Unselected"), labels = c("Clonotype 20", "All Others"), values = c("#073f7d", "#e6e5e3")) +
  ggtitle("Clonotype 20") +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_UMAP_MLR_10x_top_clone_20.eps", width = 3, height = 3)

MLR_sample1_subset@meta.data %>%
  separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% 
  filter(beta != "NA" & alpha != "NA")


MLR_subset_metadata <- MLR_sample1_subset@meta.data

MLR_subset_metadata_alphabeta <- MLR_subset_metadata %>% 
  separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% 
  filter(beta != "NA") #note, consider also alpha != "NA"

Idents(MLR_sample1_subset) <- MLR_sample1_subset$seurat_clusters

saveRDS(MLR_sample1_subset, file = "analysis/data/MLR_10x_seuratObj.rds")


# Combining bulk and 10X MLR data ---------------------------------------------

#plotting bulk frequencies of MLR clones (colored using cluster, or pre-dominate cluster)
MLR_tcr1 <- c(MLR_clonotypes$beta[1])
tc1 <- trackClonotypes(bulktcr_kid$data, MLR_tcr1, .col = "aa")
p1 <- vis(tc1, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.0004, 0.0008)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        axis.text=element_text(size=14)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 1") 

p1
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_01.png", width = 4, height = 3)

MLR_tcr2 <- c(MLR_clonotypes$beta[2])
tc2 <- trackClonotypes(bulktcr_kid$data, MLR_tcr2, .col = "aa")
p2 <- vis(tc2, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.0001, 0.0002)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 2") 

p2
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_02.png", width = 4, height = 3)

MLR_tcr3 <- c(MLR_clonotypes$beta[3])
tc3 <- trackClonotypes(bulktcr_kid$data, MLR_tcr3, .col = "aa")
p3 <- vis(tc3, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.00001, 0.00002)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 3") 

p3
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_03.png", width = 4, height = 3)

MLR_tcr4 <- c(MLR_clonotypes$beta[4])
tc4 <- trackClonotypes(bulktcr_kid$data, MLR_tcr4, .col = "aa")
p4 <- vis(tc4, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.00004, 0.00008)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 4") 

p4
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_04.png", width = 4, height = 3)

MLR_tcr5 <- c(MLR_clonotypes$beta[5])
tc5 <- trackClonotypes(bulktcr_kid$data, MLR_tcr5, .col = "aa")
p5 <- vis(tc5, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.003, 0.006)) +
  scale_fill_manual(values = singlecellcolors[2]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 5") 

p5
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_05.png", width = 4, height = 3)

MLR_tcr6 <- c(MLR_clonotypes$beta[6])
tc6 <- trackClonotypes(bulktcr_kid$data, MLR_tcr6, .col = "aa")
p6 <- vis(tc6, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.01, 0.02)) +
  scale_fill_manual(values = singlecellcolors[2]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 6") 

p6
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_06.png", width = 4, height = 3)

MLR_tcr7 <- c(MLR_clonotypes$beta[7])
tc7 <- trackClonotypes(bulktcr_kid$data, MLR_tcr7, .col = "aa")
p7 <- vis(tc7, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.01, 0.02)) +
  scale_fill_manual(values = singlecellcolors[2]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 7") 

p7
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_07.png", width = 4, height = 3)

MLR_tcr8 <- c(MLR_clonotypes$beta[8])
tc8 <- trackClonotypes(bulktcr_kid$data, MLR_tcr8, .col = "aa")
p8 <- vis(tc8, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), limits = c(0,0.0002), breaks = c(0,0.0001,0.0002)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 8") 

p8
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_08.png", width = 4, height = 3)

MLR_tcr9 <- c(MLR_clonotypes$beta[9])
tc9 <- trackClonotypes(bulktcr_kid$data, MLR_tcr9, .col = "aa")
p9 <- vis(tc9, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0, 0.00005, 0.0001)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 9") 

p9
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_09.png", width = 4, height = 3)

MLR_tcr10 <- c(MLR_clonotypes$beta[10])
tc10 <- trackClonotypes(bulktcr_kid$data, MLR_tcr10, .col = "aa")
p10 <- vis(tc10, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0,0.0004, 0.0008)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 10") 

p10
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_10.png", width = 4, height = 3)

MLR_tcr11 <- c(MLR_clonotypes$beta[11])
tc11 <- trackClonotypes(bulktcr_kid$data, MLR_tcr11, .col = "aa")
p11 <- vis(tc11, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0,0.00004, 0.00008)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 11") 

p11
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_11.png", width = 4, height = 3)

MLR_tcr12 <- c(MLR_clonotypes$beta[12])
tc12 <- trackClonotypes(bulktcr_kid$data, MLR_tcr12, .col = "aa")
p12 <- vis(tc12, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0,0.0004, 0.0008)) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 12") 

p12
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_12.png", width = 4, height = 3)

MLR_tcr13 <- c(MLR_clonotypes$beta[13])
tc13 <- trackClonotypes(bulktcr_kid$data, MLR_tcr13, .col = "aa")
p13 <- vis(tc13, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%")) +
  scale_fill_manual(values = singlecellcolors[5]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ylim(0,0.01) +
  ggtitle("Clonotype 13")

p13
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_13.png", width = 4, height = 3)

MLR_tcr14 <- c(MLR_clonotypes$beta[14])
tc14 <- trackClonotypes(bulktcr_kid$data, MLR_tcr14, .col = "aa")
p14 <- vis(tc14, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%")) +
  scale_fill_manual(values = singlecellcolors[5]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ylim(0,0.01) +
  ggtitle("Clonotype 14") 

p14
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_14.png", width = 4, height = 3)

MLR_tcr15 <- c(MLR_clonotypes$beta[15])
tc15 <- trackClonotypes(bulktcr_kid$data, MLR_tcr15, .col = "aa")
p15 <- vis(tc15, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0,0.003,0.006)) +
  scale_fill_manual(values = singlecellcolors[2]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 15") 

p15
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_15.png", width = 4, height = 3)

MLR_tcr16 <- c(MLR_clonotypes$beta[16])
tc16 <- trackClonotypes(bulktcr_kid$data, MLR_tcr16, .col = "aa")
p16 <- vis(tc16, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%")) +
  scale_fill_manual(values = singlecellcolors[1]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 16") 

p16
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_16.png", width = 4, height = 3)

MLR_tcr17 <- c(MLR_clonotypes$beta[17])
tc17 <- trackClonotypes(bulktcr_kid$data, MLR_tcr17, .col = "aa")
p17 <- vis(tc17, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%")) +
  scale_fill_manual(values = singlecellcolors[2]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 17") 

p17
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_17.png", width = 4, height = 3)

MLR_tcr18 <- c(MLR_clonotypes$beta[18])
tc18 <- trackClonotypes(bulktcr_kid$data, MLR_tcr18, .col = "aa")
p18 <- vis(tc18, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_fill_manual(values = singlecellcolors[5]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ylim(0,0.01) +
  ggtitle("Clonotype 18") 

p18
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_18.png", width = 4, height = 3)

MLR_tcr19 <- c(MLR_clonotypes$beta[19])
tc19 <- trackClonotypes(bulktcr_kid$data, MLR_tcr19, .col = "aa")
p19 <- vis(tc19, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%"), breaks = c(0,0.00001)) +
  scale_fill_manual(values = singlecellcolors[5]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 19") 

p19
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_19.png", width = 4, height = 3)

MLR_tcr20 <- c(MLR_clonotypes$beta[20])
tc20 <- trackClonotypes(bulktcr_kid$data, MLR_tcr20, .col = "aa")
p20 <- vis(tc20, .plot = "smooth") +
  theme_classic() +
  scale_x_discrete(labels = bulktcr_kid$meta$sampleShort) +
  scale_y_continuous(labels = function(x) paste0(100*x, "%")) +
  scale_fill_manual(values = singlecellcolors[6]) +
  theme(axis.text.x = element_text(angle=45, hjust = 1),
        axis.text=element_text(size=14)) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle("Clonotype 20") 

p20
ggsave("analysis/figs/bulk_alluvial_MLR_10x_top_clone_20.png", width = 4, height = 3)

bulkpercentplots <- align_plots(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, align="hv", axis="tblr")
p1x <- ggdraw(bulkpercentplots[[1]])
p2x <- ggdraw(bulkpercentplots[[2]])
p3x <- ggdraw(bulkpercentplots[[3]])
p4x <- ggdraw(bulkpercentplots[[4]])
p5x <- ggdraw(bulkpercentplots[[5]])
p6x <- ggdraw(bulkpercentplots[[6]])
p7x <- ggdraw(bulkpercentplots[[7]])
p8x <- ggdraw(bulkpercentplots[[8]])
p9x <- ggdraw(bulkpercentplots[[9]])
p10x <- ggdraw(bulkpercentplots[[10]])
p11x <- ggdraw(bulkpercentplots[[11]])
p12x <- ggdraw(bulkpercentplots[[12]])
p13x <- ggdraw(bulkpercentplots[[13]])
p14x <- ggdraw(bulkpercentplots[[14]])
p15x <- ggdraw(bulkpercentplots[[15]])
p16x <- ggdraw(bulkpercentplots[[16]])
p17x <- ggdraw(bulkpercentplots[[17]])
p18x <- ggdraw(bulkpercentplots[[18]])
p19x <- ggdraw(bulkpercentplots[[19]])
p20x <- ggdraw(bulkpercentplots[[20]])

save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_01_aligned.png", p1x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_02_aligned.png", p2x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_03_aligned.png", p3x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_04_aligned.png", p4x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_05_aligned.png", p5x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_06_aligned.png", p6x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_07_aligned.png", p7x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_08_aligned.png", p8x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_09_aligned.png", p9x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_10_aligned.png", p10x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_11_aligned.png", p11x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_12_aligned.png", p12x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_13_aligned.png", p13x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_14_aligned.png", p14x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_15_aligned.png", p15x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_16_aligned.png", p16x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_17_aligned.png", p17x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_18_aligned.png", p18x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_19_aligned.png", p19x, base_width = 4, base_height = 3)
save_plot("analysis/figs/bulk_alluvial_MLR_10x_top_clone_20_aligned.png", p20x, base_width = 4, base_height = 3)

#adding putative CD4/CD8 annotation to bulk TCRs that match MLR clones
MLR_subset_counts <- MLR_sample1_subset@assays$RNA@data
MLR_subset_counts <- as.data.frame(MLR_subset_counts)
MLR_subset_counts <- t(MLR_subset_counts)
MLR_subset_counts <- as_tibble(MLR_subset_counts, rownames = NA)

MLR_CD4_CD8_counts <- cbind(MLR_subset_counts$CD4, MLR_subset_counts$CD8A)
rownames(MLR_CD4_CD8_counts) <- rownames(MLR_subset_counts)
colnames(MLR_CD4_CD8_counts) <- c("CD4", "CD8")

MLR_meta_CD4_CD8_counts <- merge(MLR_subset_metadata,  MLR_CD4_CD8_counts, by = 0)

MLR_clone_CD4CD8 <- MLR_meta_CD4_CD8_counts %>%
  group_by(CTaa) %>%
  summarise(meanCD4 = mean(CD4), meanCD8 = mean(CD8), n = n()) %>% 
  mutate(type = ifelse(meanCD4 > meanCD8, "CD4", ifelse(meanCD4 < meanCD8, "CD8", "Other"))) %>%
  separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>%
  filter(beta != "NA" & alpha != "NA") %>%
  select(beta, type, n) %>% arrange(desc(n)) %>%
  rename(CDR3aa = beta)

bulktcr_cdr3_counts <- unique(dplyr::left_join(bulktcr_cdr3_counts, MLR_clone_CD4CD8[, c("CDR3aa", "type")], by="CDR3aa"))
bulktcr_cdr3_freqs <- unique(dplyr::left_join(bulktcr_cdr3_freqs, MLR_clone_CD4CD8[, c("CDR3aa", "type")], by="CDR3aa"))

MLR_subset_metadata %>% 
  separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% 
  filter(beta != "NA" & alpha != "NA") %>% 
  add_count(seurat_clusters, name = "cluster_totaln") %>% 
  add_count(seurat_clusters, beta, name = "clonotypesize") %>% 
  select(beta, seurat_clusters, cluster_totaln, clonotypesize) %>%
  mutate(found_kid = ifelse((beta %in% bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa), "yes", "no")) %>% 
  mutate(found_ln = ifelse((beta %in% bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa), "yes", "no")) %>% 
  mutate(found_skin = ifelse((beta %in% bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa), "yes", "no")) %>% 
  mutate(found_pbl0 = ifelse((beta %in% bulktcr_kid$data$`GM-Pre_PBMC_180728_TCRB`$CDR3.aa), "yes", "no")) %>% 
  mutate(found_pbl1 = ifelse((beta %in% bulktcr_kid$data$`GM-C1_PBMC_180816_TCRB`$CDR3.aa), "yes", "no")) %>% 
  mutate(found_pbl2 = ifelse((beta %in% bulktcr_kid$data$`GM-C2_PBMC_180927_TCRB`$CDR3.aa), "yes", "no")) %>% 
  count(beta, seurat_clusters, cluster_totaln, found_kid, found_ln, found_skin, found_pbl0, found_pbl1, found_pbl2, name = "clonotypesize") %>%
  gather(bulk_sample, found, found_kid:found_pbl2, factor_key=TRUE) %>% 
  group_by(seurat_clusters, bulk_sample, found, cluster_totaln) %>%
  summarise(test = sum(clonotypesize)) %>%
  mutate(freqofclone = test / cluster_totaln) %>%
  filter(found == "yes") %>%
  complete(seurat_clusters, bulk_sample, fill = list(freqofclone = 0)) %>%
  distinct(.keep_all=TRUE) %>%
  ggplot(aes(fill=bulk_sample, y = freqofclone, x = seurat_clusters)) +
  geom_bar(position = "dodge", stat = "identity", colour="black") +
  scale_fill_manual(values = bulktcrcolors, labels = bulktcr_kid$meta$sampleShort, name = "Bulk Sample") +
  theme_classic() +
  ylab("Frequency of cells in cluster found in sample") + xlab(NULL) + 
  scale_y_continuous(expand = c(0, 0)) 
ggsave("analysis/figs/sc_clusters_bulk_representation.eps", width = 8, height = 4)


Idents(MLR_sample1_subset) <- MLR_sample1_subset$barcode

MLR_kidney_clones <- MLR_subset_metadata_alphabeta %>%
  filter(beta %in% bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa) 

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_kidney_clones$barcode)) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[1])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Kidney\n Matched Clones")
ggsave("analysis/figs/sc_UMAP_MLR_kidneyclones.eps", width = 4, height = 4)

MLR_ln_clones <- MLR_subset_metadata_alphabeta %>%
  filter(beta %in% bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa) 

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_ln_clones$barcode)) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[2])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("LN Met\n Matched Clones")
ggsave("analysis/figs/sc_UMAP_MLR_lnclones.eps", width = 4, height = 4)

MLR_skin_clones <- MLR_subset_metadata_alphabeta %>%
  filter(beta %in% bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa) 

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_skin_clones$barcode)) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[3])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Skin Met\n Matched Clones")
ggsave("analysis/figs/sc_UMAP_MLR_skinclones.eps", width = 4, height = 4)

MLR_prepembro_clones <- MLR_subset_metadata_alphabeta %>%
  filter(beta %in% bulktcr_kid$data$`GM-Pre_PBMC_180728_TCRB`$CDR3.aa) 

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_prepembro_clones$barcode)) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[4])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Pre Pembro\n Matched Clones")
ggsave("analysis/figs/sc_UMAP_MLR_prepembroclones.eps", width = 4, height = 4)

MLR_pembro1_clones <- MLR_subset_metadata_alphabeta %>%
  filter(beta %in% bulktcr_kid$data$`GM-C1_PBMC_180816_TCRB`$CDR3.aa) 

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_pembro1_clones$barcode)) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[5])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Pembro 1\n Matched Clones")
ggsave("analysis/figs/sc_UMAP_MLR_pembro1clones.eps", width = 4, height = 4)

MLR_pembro2_clones <- MLR_subset_metadata_alphabeta %>%
  filter(beta %in% bulktcr_kid$data$`GM-C2_PBMC_180927_TCRB`$CDR3.aa) 

DimPlot(MLR_sample1_subset, cells.highlight = WhichCells(MLR_sample1_subset, idents = MLR_pembro2_clones$barcode)) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[6])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Pembro 2\n Matched Clones")
ggsave("analysis/figs/sc_UMAP_MLR_pembro2clones.eps", width = 4, height = 4)

Idents(MLR_sample1_subset) <- MLR_sample1_subset$CTaa



#getting clone counts
MLRtotalclones_percluster <- MLR_subset_metadata_alphabeta %>% 
  group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "totalclones")

MLRkidneyclones_percluster <- MLR_subset_metadata_alphabeta %>% 
  filter(beta %in% bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa) %>% group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "kidney")

MLRlnclones_percluster <- MLR_subset_metadata_alphabeta %>% 
  filter(beta %in% bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa) %>% group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "ln")

MLRskinclones_percluster <- MLR_subset_metadata_alphabeta %>% 
  filter(beta %in% bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa) %>% group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "skin")

MLRprepembroclones_percluster <- MLR_subset_metadata_alphabeta %>% 
  filter(beta %in% bulktcr_kid$data$`GM-Pre_PBMC_180728_TCRB`$CDR3.aa) %>% group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "prepembro")

MLRpembro1clones_percluster <- MLR_subset_metadata_alphabeta %>% 
  filter(beta %in% bulktcr_kid$data$`GM-C1_PBMC_180816_TCRB`$CDR3.aa) %>% group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "pembro1")

MLRpembro2clones_percluster <- MLR_subset_metadata_alphabeta %>% 
  filter(beta %in% bulktcr_kid$data$`GM-C2_PBMC_180927_TCRB`$CDR3.aa) %>% group_by(alpha, beta, seurat_clusters) %>% distinct(alpha, beta, seurat_clusters) %>% group_by(seurat_clusters) %>% count(name = "pembro2")

MLRclonesbulkmatch_percluster <- left_join(MLRtotalclones_percluster, MLRkidneyclones_percluster, by='seurat_clusters') %>%
  left_join(., MLRlnclones_percluster, by='seurat_clusters') %>%
  left_join(., MLRskinclones_percluster, by='seurat_clusters') %>%
  left_join(., MLRprepembroclones_percluster, by='seurat_clusters') %>%
  left_join(., MLRpembro1clones_percluster, by='seurat_clusters') %>%
  left_join(., MLRpembro2clones_percluster, by='seurat_clusters') %>% 
  replace(is.na(.), 0)

p1 <- MLRclonesbulkmatch_percluster %>% 
  mutate(cluster1 = ifelse(seurat_clusters == "C1", "C1", "Other")) %>% 
  pivot_longer(kidney:pembro2, names_to = "tissue", values_to = "clones_matched") %>% 
  group_by(cluster1, tissue) %>% 
  add_tally(clones_matched, name = "clones_matched_tissue_sum_C1vsother") %>% 
  add_tally(totalclones, name = "total_clones_C1vsother") %>% 
  mutate (percent_matched = clones_matched_tissue_sum_C1vsother / total_clones_C1vsother) %>% 
  select(cluster1, tissue, percent_matched) %>% 
  distinct() %>%
  filter(tissue %in% c("kidney", "ln", "skin")) %>%
  ggplot(aes(fill=cluster1, y=percent_matched, x=tissue)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(limits = c(0,1), expand = c(0,0), labels = function(x) paste0(100*x, "%")) +
  scale_x_discrete(limits = c("kidney", "ln", "skin"), labels = c("Kidney", "LN Met", "Skin Met")) + 
  theme_classic() + 
  scale_fill_manual(values = c(singlecellcolors[1], "grey"), name = NULL) + 
  xlab(NULL) + ylab("Clones Matched to Bulk TCRseq") + theme(axis.text.x = element_text(angle=45, hjust = 1),
                                                             axis.text=element_text(size=10),
                                                             axis.title = element_text(size=10),
                                                             legend.text = element_text(size=9),
                                                             legend.key.size = unit(0.5, "cm")) 

p2 <- MLRclonesbulkmatch_percluster %>% 
  mutate(cluster1 = ifelse(seurat_clusters == "C1", "C1", "Other")) %>% 
  pivot_longer(kidney:pembro2, names_to = "tissue", values_to = "clones_matched") %>% 
  group_by(cluster1, tissue) %>% 
  add_tally(clones_matched, name = "clones_matched_tissue_sum_C1vsother") %>% 
  add_tally(totalclones, name = "total_clones_C1vsother") %>% 
  mutate (percent_matched = clones_matched_tissue_sum_C1vsother / total_clones_C1vsother) %>% 
  select(cluster1, tissue, percent_matched) %>% 
  distinct() %>%
  filter(tissue %in% c("prepembro", "pembro1", "pembro2")) %>%
  ggplot(aes(fill=cluster1, y=percent_matched, x=tissue)) + 
  geom_bar(position="dodge", stat="identity") + scale_y_continuous(limits = c(0,1), expand = c(0,0), labels = function(x) paste0(100*x, "%")) +
  scale_x_discrete(limits = c("prepembro", "pembro1", "pembro2"), labels = c("Pre Pembro", "Pembro 1", "Pembro 2")) + 
  theme_classic() + 
  scale_fill_manual(values = c(singlecellcolors[1], "grey"), name = NULL) + 
  xlab(NULL) + ylab("Clones Matched to Bulk TCRseq") + theme(axis.text.x = element_text(angle=45, hjust = 1),
                                                             axis.text=element_text(size=10),
                                                             axis.title = element_text(size=10),
                                                             legend.text = element_text(size=9),
                                                             legend.key.size = unit(0.5, "cm")) 

mlrclonesbloodline <- MLRclonesbulkmatch_percluster %>% 
  mutate(cluster1 = ifelse(seurat_clusters == "C1", "C1", "Other")) %>% 
  pivot_longer(kidney:pembro2, names_to = "tissue", values_to = "clones_matched") %>% 
  group_by(cluster1, tissue) %>% 
  add_tally(clones_matched, name = "clones_matched_tissue_sum_C1vsother") %>% 
  add_tally(totalclones, name = "total_clones_C1vsother") %>% 
  mutate (percent_matched = clones_matched_tissue_sum_C1vsother / total_clones_C1vsother) %>% 
  select(cluster1, tissue, percent_matched) %>% 
  distinct() %>%
  filter(tissue %in% c("prepembro", "pembro1", "pembro2")) 

mlrclonesbloodline$cluster1 <- factor(mlrclonesbloodline$cluster1, levels = c("Other", "C1"))

p3 <- ggplot(mlrclonesbloodline, aes(color=cluster1, y=percent_matched, x=tissue, group = cluster1)) + 
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0,1), expand = c(0,0), labels = function(x) paste0(100*x, "%")) +
  scale_x_discrete(limits = c("prepembro", "pembro1", "pembro2"), labels = c("Pre Pembro", "Pembro 1", "Pembro 2")) + 
  theme_classic() + 
  scale_color_manual(values = c(singlecellcolors[1], "grey"), name = NULL, limits = c("C1", "Other")) + 
  xlab(NULL) + ylab("Clones Matched to Bulk TCRseq") + theme(axis.text.x = element_text(angle=45, hjust = 1),
                                                             axis.text=element_text(size=10),
                                                             axis.title = element_text(size=10),
                                                             legend.text = element_text(size=9),
                                                             legend.key.size = unit(0.5, "cm")) 


both3 <- align_plots(p1, p2, p3, align="hv", axis="tblr")
p1x <- ggdraw(both3[[1]])
p2x <- ggdraw(both3[[2]])
p3x <- ggdraw(both3[[3]])

save_plot("analysis/figs/MLRvsbulk_C3clonematch_barplot_tissue.eps", p1x, base_width = 4.5, base_height = 2.8)
save_plot("analysis/figs/MLRvsbulk_C3clonematch_barplot_pbmc.eps", p2x, base_width = 4.5, base_height = 2.8)
save_plot("analysis/figs/MLRvsbulk_C3clonematch_line_pbmc.eps", p3x, base_width = 4.5, base_height = 2.8)

MLR_subset_metadata %>% 
  separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% 
  filter(beta != "NA") %>% filter(beta %in% bulktcr_kid$data$`GM-Pre_PBMC_180728_TCRB`$CDR3.aa) %>% filter(seurat_clusters == "C1")

MLRpembro2expansion <- MLR_subset_metadata_alphabeta %>%
  group_by(beta, seurat_clusters) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  ungroup() %>%
  mutate(freq = n / sum(n)) %>%
  left_join(., bulktcr_kid$data$`GM-C2_PBMC_180927_TCRB`, by = c("beta" = "CDR3.aa")) %>%
  select(beta, seurat_clusters, n, freq, Proportion) %>%
  mutate(MLRvspembro = freq / Proportion)
MLRpembro2expansion$MLRvspembro[is.infinite(MLRpembro2expansion$MLRvspembro)] <- NA  
MLRpembro2expansion %>% na.omit(MLRvspembro) %>% 
  ggplot(aes(x = seurat_clusters, y = MLRvspembro, size = n, color = seurat_clusters)) +
  geom_point(position = position_jitter(width = 0.15, height = 0.15)) + scale_color_manual(values = singlecellcolors) +   
  scale_x_discrete(drop = FALSE) +
  labs(size = "MLR Clone \nCount", color = "Cluster", y = "MLR freq / Blood freq", x = NULL) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 1, 10, 100, 1000), limits = c(NA, 1000)) + theme_classic() + guides(color = "none")  
ggsave("analysis/figs/MLRvspembro2expansion_byCluster.eps", width = 5, height = 3)

MLRpembro1expansion <- MLR_subset_metadata_alphabeta %>%
  group_by(beta, seurat_clusters) %>% 
  tally() %>% 
  arrange(desc(n)) %>% 
  ungroup() %>%
  mutate(freq = n / sum(n)) %>%
  left_join(., bulktcr_kid$data$`GM-C1_PBMC_180816_TCRB`, by = c("beta" = "CDR3.aa")) %>%
  select(beta, seurat_clusters, n, freq, Proportion) %>%
  mutate(MLRvspembro = freq / Proportion)
MLRpembro1expansion$MLRvspembro[is.infinite(MLRpembro1expansion$MLRvspembro)] <- NA  
MLRpembro1expansion %>% na.omit(MLRvspembro) %>% 
  ggplot(aes(x = seurat_clusters, y = MLRvspembro, size = n, color = seurat_clusters)) +
  geom_point(position = position_jitter(width = 0.15, height = 0.15)) +  scale_color_manual(values = singlecellcolors) +   
  scale_x_discrete(drop = FALSE) +
  labs(size = "MLR Clone \nCount", color = "Cluster", y = "MLR freq / Pembro 1 freq", x = NULL) +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), 
                     breaks = c(0, 1, 10, 100, 1000), limits = c(NA, 5000)) + theme_classic() + guides(color = "none")  
ggsave("analysis/figs/MLRvspembro1expansion_byCluster.eps", width = 5, height = 3)

# Single-cell TCR and transcriptome analysis of non-naive CD8 post-pembro -----------------------
#load 10x data
#pembro1
pembro1_nonnaiveCD8_data <- Read10X_h5("nonNaiveCD8_scRNA/data/pembro1/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
pembro1_nonnaiveCD8 <- CreateSeuratObject(pembro1_nonnaiveCD8_data, project = "Pembro1")
pembro1_nonnaiveCD8$timepoint <- "Pembro1"

pembro1_nonnaiveCD8_TCR_data <- read.csv("nonNaiveCD8_scRNA/data/pembro1/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
TCR_list <- list(pembro1_nonnaiveCD8_TCR_data)
pembro1_nonnaiveCD8_TCR <- combineTCR(TCR_list, samples = "Pembro", ID = "1", cells = "T-AB", filterMulti = T)
pembro1_nonnaiveCD8_TCR$Pembro_1$barcode <- substring(pembro1_nonnaiveCD8_TCR$Pembro_1$barcode, 10)

pembro1_nonnaiveCD8 <- combineExpression(pembro1_nonnaiveCD8_TCR$Pembro_1, pembro1_nonnaiveCD8, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=3, Medium=30, Large=100, Hyperexpanded=2000))

pembro1_nonnaiveCD8$cloneType <- as.character(pembro1_nonnaiveCD8$cloneType)
pembro1_nonnaiveCD8$cloneType <- replace_na(pembro1_nonnaiveCD8$cloneType, "No TCR")
Idents(pembro1_nonnaiveCD8) <- pembro1_nonnaiveCD8$cloneType

slot(pembro1_nonnaiveCD8, "meta.data")$cloneType <- factor(slot(pembro1_nonnaiveCD8, "meta.data")$cloneType, 
                                                           levels = c("Hyperexpanded (100 < X <= 2000)", 
                                                                      "Large (30 < X <= 100)", 
                                                                      "Medium (3 < X <= 30)", 
                                                                      "Small (1 < X <= 3)", 
                                                                      "Single (0 < X <= 1)",
                                                                      "No TCR"))

current.clone.names <- c("Hyperexpanded (100 < X <= 2000)", "Large (30 < X <= 100)", "Medium (3 < X <= 30)", "Small (1 < X <= 3)", "Single (0 < X <= 1)", "No TCR")
new.clone.names <- c("Hyperexpanded (100 < X)", "Large (30 < X <= 100)", "Medium (3 < X <= 30)", "Small (1 < X <= 3)", "Single (0 < X <= 1)", "No TCR information")
pembro1_nonnaiveCD8@active.ident <- plyr::mapvalues(x = pembro1_nonnaiveCD8@active.ident, from = current.clone.names, to = new.clone.names)
pembro1_nonnaiveCD8$cloneType <- pembro1_nonnaiveCD8@active.ident

slot(pembro1_nonnaiveCD8, "meta.data")$cloneType <- factor(slot(pembro1_nonnaiveCD8, "meta.data")$cloneType, 
                                                           levels = c("Hyperexpanded (100 < X)", 
                                                                      "Large (30 < X <= 100)", 
                                                                      "Medium (3 < X <= 30)", 
                                                                      "Small (1 < X <= 3)", 
                                                                      "Single (0 < X <= 1)", 
                                                                      "No TCR information"))

pembro1_nonnaiveCD8$clonal <- ifelse(pembro1_nonnaiveCD8$Frequency > 1, "Yes", "No")

#pembro2
pembro2_nonnaiveCD8_data <- Read10X_h5("nonNaiveCD8_scRNA/data/pembro2/RNA/filtered_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE)
pembro2_nonnaiveCD8 <- CreateSeuratObject(pembro2_nonnaiveCD8_data, project = "Pembro2")
pembro2_nonnaiveCD8$timepoint <- "Pembro2"

pembro2_nonnaiveCD8_TCR_data <- read.csv("nonNaiveCD8_scRNA/data/pembro2/TCR/filtered_contig_annotations.csv", stringsAsFactors = F)
TCR_list <- list(pembro2_nonnaiveCD8_TCR_data)
pembro2_nonnaiveCD8_TCR <- combineTCR(TCR_list, samples = "Pembro", ID = "2", cells = "T-AB", filterMulti = T)
pembro2_nonnaiveCD8_TCR$Pembro_2$barcode <- substring(pembro2_nonnaiveCD8_TCR$Pembro_2$barcode, 10)

pembro2_nonnaiveCD8 <- combineExpression(pembro2_nonnaiveCD8_TCR$Pembro_2, pembro2_nonnaiveCD8, cloneCall="aa", proportion = FALSE, cloneTypes=c(Single=1, Small=3, Medium=30, Large=100, Hyperexpanded=2000))

pembro2_nonnaiveCD8$cloneType <- as.character(pembro2_nonnaiveCD8$cloneType)
pembro2_nonnaiveCD8$cloneType <- replace_na(pembro2_nonnaiveCD8$cloneType, "No TCR")
Idents(pembro2_nonnaiveCD8) <- pembro2_nonnaiveCD8$cloneType

slot(pembro2_nonnaiveCD8, "meta.data")$cloneType <- factor(slot(pembro2_nonnaiveCD8, "meta.data")$cloneType, 
                                                           levels = c("Hyperexpanded (100 < X <= 2000)", 
                                                                      "Large (30 < X <= 100)", 
                                                                      "Medium (3 < X <= 30)", 
                                                                      "Small (1 < X <= 3)", 
                                                                      "Single (0 < X <= 1)",
                                                                      "No TCR"))

current.clone.names <- c("Hyperexpanded (100 < X <= 2000)", "Large (30 < X <= 100)", "Medium (3 < X <= 30)", "Small (1 < X <= 3)", "Single (0 < X <= 1)", "No TCR")
new.clone.names <- c("Hyperexpanded (100 < X)", "Large (30 < X <= 100)", "Medium (3 < X <= 30)", "Small (1 < X <= 3)", "Single (0 < X <= 1)", "No TCR information")
pembro2_nonnaiveCD8@active.ident <- plyr::mapvalues(x = pembro2_nonnaiveCD8@active.ident, from = current.clone.names, to = new.clone.names)
pembro2_nonnaiveCD8$cloneType <- pembro2_nonnaiveCD8@active.ident

slot(pembro2_nonnaiveCD8, "meta.data")$cloneType <- factor(slot(pembro2_nonnaiveCD8, "meta.data")$cloneType, 
                                                           levels = c("Hyperexpanded (100 < X)", 
                                                                      "Large (30 < X <= 100)", 
                                                                      "Medium (3 < X <= 30)", 
                                                                      "Small (1 < X <= 3)", 
                                                                      "Single (0 < X <= 1)", 
                                                                      "No TCR information"))

pembro2_nonnaiveCD8$clonal <- ifelse(pembro2_nonnaiveCD8$Frequency > 1, "Yes", "No")

#merge
nonnaiveCD8_Seurat <- merge(pembro1_nonnaiveCD8, y = c(pembro2_nonnaiveCD8), project = "PostPembro_nonNaiveCD8_RNA_TCR")

slot(nonnaiveCD8_Seurat, "meta.data")$cloneType <- factor(slot(nonnaiveCD8_Seurat, "meta.data")$cloneType, 
                                                          levels = c("Hyperexpanded (100 < X)", 
                                                                     "Large (30 < X <= 100)", 
                                                                     "Medium (3 < X <= 30)", 
                                                                     "Small (1 < X <= 3)", 
                                                                     "Single (0 < X <= 1)", 
                                                                     "No TCR information"))


#QC and initial subset
nonnaiveCD8_Seurat[["percent.mito"]] <- PercentageFeatureSet(nonnaiveCD8_Seurat, pattern = "^MT-")
nonnaiveCD8_Seurat[["percent.ribo"]] <- PercentageFeatureSet(nonnaiveCD8_Seurat, pattern = "^RP[SL]")

Idents(nonnaiveCD8_Seurat) <- nonnaiveCD8_Seurat$timepoint

timepoint_cols <- c("#D9514EFF", "#2DA8D8FF")
VlnPlot(nonnaiveCD8_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3) +
  theme(legend.position = "none") &
  xlab(NULL)
a <- VlnPlot(nonnaiveCD8_Seurat, features = c("nFeature_RNA"), pt.size = 0, cols = timepoint_cols) +
  theme(legend.position = "none") + ggtitle("Number of Genes") +
  xlab(NULL)
b <- VlnPlot(nonnaiveCD8_Seurat, features = c("percent.mito"),  pt.size = 0, cols = timepoint_cols) +
  theme(legend.position = "none") + ggtitle("% Mitochondrial Reads") +
  xlab(NULL)

a | b

preqc <- table(nonnaiveCD8_Seurat$timepoint)

nonNaiveCD8_preQC_cellnum <- nonnaiveCD8_Seurat@meta.data %>% group_by(timepoint) %>% count() #pembro1 7704, pembro2 8479

nonNaiveCD8_preQC_VDJcellnum <- nonnaiveCD8_Seurat@meta.data %>% filter(cloneType != "No TCR information") %>% group_by(timepoint) %>% count()  #pembro1 6103, pembro2 7159


nonnaiveCD8_Seurat_subset <- subset(nonnaiveCD8_Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 15)


nonNaiveCD8_postQC_cellnum <- nonnaiveCD8_Seurat_subset@meta.data %>% group_by(timepoint) %>% count()  #pembro1 6403, pembro2 7783

nonNaiveCD8_postQC_VDJcellnum <- nonnaiveCD8_Seurat_subset@meta.data %>% filter(cloneType != "No TCR information") %>% group_by(timepoint) %>% count()  #pembro1 5958, pembro2 7076

nonNaiveCD8_postQC_clonalNonclonalnums <- nonnaiveCD8_Seurat_subset@meta.data %>% group_by(clonal, timepoint) %>% count()

postqc <- table(nonnaiveCD8_Seurat_subset$timepoint)
filtered <- preqc - postqc

nonnaive_cellnums <- data.frame(Timepoint = c("Pembro1", "Pembro2", "Pembro1", "Pembro2"), QC = c("preQC", "preQC", "postQC", "postQC"), Nums = c(preqc, postqc))
nonnaive_cellnums$timeqc <- paste0(nonnaive_cellnums$Timepoint, "_", nonnaive_cellnums$QC)
nonnaive_cellnums$QC <- factor(nonnaive_cellnums$QC, levels = c("preQC", "postQC"))
nonnaive_cellnums$timeqc <- factor(nonnaive_cellnums$timeqc, levels = c("Pembro1_preQC", "Pembro1_postQC", "Pembro2_preQC", "Pembro2_postQC"))

ggplot(nonnaive_cellnums, aes(x = Timepoint, y = Nums, label = Nums, fill = timeqc)) + geom_bar(position="dodge", stat="identity") +
  geom_text(position = position_dodge(width = 1), vjust = 2, size = 4) +
  scale_fill_manual(values = c("#831e1c", timepoint_cols[1], "#1f75b2", timepoint_cols[2])) + 
  xlab(NULL) + ylab("Number of Cells") +
  theme_classic() +
  theme(legend.position = "none")
ggsave("analysis/figs/sc_nonnaiveCD8_QC_barplot.eps", width = 3, height = 2.5)

#cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
nonnaiveCD8_Seurat_subset <- CellCycleScoring(nonnaiveCD8_Seurat_subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
slot(nonnaiveCD8_Seurat_subset, "meta.data")$Phase <- factor(slot(nonnaiveCD8_Seurat_subset, "meta.data")$Phase, 
                                                             levels = c("G1", "S", "G2M"))

#spliting CDR3s to separate alpha and beta
nonnaiveCD8_Seurat_subset$barcode_timepoint <- rownames(nonnaiveCD8_Seurat_subset@meta.data)

nonnaiveCD8_Seurat_subset$beta <- sub(".*_", "", nonnaiveCD8_Seurat_subset$CTaa) 
nonnaiveCD8_Seurat_subset$beta <- ifelse(is.na(nonnaiveCD8_Seurat_subset$beta), "NA", nonnaiveCD8_Seurat_subset$beta)

nonnaiveCD8_Seurat_subset$alpha <- sub("_.*", "", nonnaiveCD8_Seurat_subset$CTaa) 
nonnaiveCD8_Seurat_subset$alpha <- ifelse(is.na(nonnaiveCD8_Seurat_subset$alpha), "NA", nonnaiveCD8_Seurat_subset$alpha)

nonnaive_topclone_list <- nonnaiveCD8_Seurat_subset@meta.data  %>% group_by(beta) %>% count() %>% filter(beta != "NA") %>% arrange(desc(n)) %>% na.omit()

clone1cells <- nonnaiveCD8_Seurat_subset@meta.data %>% filter(beta == nonnaive_topclone_list$beta[1]) %>% select(barcode_timepoint)
clone2cells <- nonnaiveCD8_Seurat_subset@meta.data %>% separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% filter(beta == nonnaive_topclone_list$beta[2]) %>% select(barcode_timepoint)
clone3cells <- nonnaiveCD8_Seurat_subset@meta.data %>% separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% filter(beta == nonnaive_topclone_list$beta[3]) %>% select(barcode_timepoint)
clone4cells <- nonnaiveCD8_Seurat_subset@meta.data %>% separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% filter(beta == nonnaive_topclone_list$beta[4]) %>% select(barcode_timepoint)
clone5cells <- nonnaiveCD8_Seurat_subset@meta.data %>% separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% filter(beta == nonnaive_topclone_list$beta[5]) %>% select(barcode_timepoint)

#dimensionality reduction with percent mito regress and harmonized by timepoint
nonnaiveCD8_Seurat_subset <- NormalizeData(nonnaiveCD8_Seurat_subset, normalization.method = "LogNormalize", scale.factor = 10000)
nonnaiveCD8_Seurat_subset <- FindVariableFeatures(nonnaiveCD8_Seurat_subset, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(nonnaiveCD8_Seurat_subset)
nonnaiveCD8_Seurat_subset <- ScaleData(nonnaiveCD8_Seurat_subset, vars.to.regress = c("percent.mito"))
nonnaiveCD8_Seurat_subset <- RunPCA(nonnaiveCD8_Seurat_subset, features = VariableFeatures(object = nonnaiveCD8_Seurat_subset))

dims_use = 1:30
nonnaiveCD8_Seurat_subset <- RunHarmony(nonnaiveCD8_Seurat_subset, group.by.vars=c("timepoint"))
nonnaiveCD8_Seurat_subset <- RunUMAP(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindNeighbors(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindClusters(object=nonnaiveCD8_Seurat_subset, resolution=0.15, verbose=FALSE)

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = TRUE, logfc.threshold = 0.50, min.pct = 0.30)

a <- DimPlot(nonnaiveCD8_Seurat_subset, label = T, repel = T) +
  scale_color_brewer(palette = "Paired") + xlab("UMAP 1") + ylab("UMAP 2")

b <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "timepoint", shuffle = T, cols = timepoint_cols) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(NULL) 

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$beta
c <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[1]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[2]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[3]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[4])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Group_2", "Group_3", "Group_4", "Unselected"), 
                     labels = c("Clone #1", "Clone #2", "Clone #3", "Clone #4", "Other"),
                     values = c(brewer.pal(4, "Set2"), "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$seurat_clusters

d <- FeaturePlot(nonnaiveCD8_Seurat_subset, features = "GZMB") + 
  scale_color_distiller(palette = "Reds", direction = 1, name = "GZMB") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(NULL)

wrap_plots(a, b, c, d)

#dimensionality reduction with removing VDJ as variable, percent mito regress, and harmonized by timepoint
"%!in%" <- Negate("%in%")
unwanted_genes <- "^TRBV*|^TRBD*|^TRBJ*|^TRDV*|^TRDD*|^TRDJ*|^TRAV*|^TRAJ*|^TRGV*|^TRGJ*"
unwanted_genes <- grep(pattern = unwanted_genes, x = nonnaiveCD8_Seurat_subset[["RNA"]]@var.features, value = T)
nonnaiveCD8_Seurat_subset@assays$RNA@var.features <- nonnaiveCD8_Seurat_subset@assays$RNA@var.features[nonnaiveCD8_Seurat_subset@assays$RNA@var.features %!in% unwanted_genes]

nonnaiveCD8_Seurat_subset <- ScaleData(nonnaiveCD8_Seurat_subset, vars.to.regress = c("percent.mito"))
nonnaiveCD8_Seurat_subset <- RunPCA(nonnaiveCD8_Seurat_subset, features = VariableFeatures(object = nonnaiveCD8_Seurat_subset))

dims_use = 1:30
nonnaiveCD8_Seurat_subset <- RunHarmony(nonnaiveCD8_Seurat_subset, group.by.vars=c("timepoint"))
nonnaiveCD8_Seurat_subset <- RunUMAP(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindNeighbors(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindClusters(object=nonnaiveCD8_Seurat_subset, resolution=0.2, verbose=FALSE)

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = TRUE, logfc.threshold = 0.50, min.pct = 0.30)

a <- DimPlot(nonnaiveCD8_Seurat_subset, label = T, repel = T) +
  scale_color_brewer(palette = "Paired") + xlab("UMAP 1") + ylab("UMAP 2")

b <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "timepoint", shuffle = T, cols = timepoint_cols) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(NULL) 

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$beta
c <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[1]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[2]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[3]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[4])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Group_2", "Group_3", "Group_4", "Unselected"), 
                     labels = c("Clone #1", "Clone #2", "Clone #3", "Clone #4", "Other"),
                     values = c(brewer.pal(4, "Set2"), "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$seurat_clusters

d <- FeaturePlot(nonnaiveCD8_Seurat_subset, features = "GZMB") + 
  scale_color_distiller(palette = "Reds", direction = 1, name = "GZMB") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(NULL)

wrap_plots(a, b, c, d)

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$beta
e <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[1])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #1", "Other"),
                     values = c(brewer.pal(4, "Set2")[1], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

f <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[2])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #2", "Other"),
                     values = c(brewer.pal(4, "Set2")[2], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

g <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[3])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #3", "Other"),
                     values = c(brewer.pal(4, "Set2")[3], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

h <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[4])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #4", "Other"),
                     values = c(brewer.pal(4, "Set2")[4], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

wrap_plots(e, f, g, h)
Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$seurat_clusters

#dimensionality reduction with removing VDJ as variable, percent mito regress, and harmonized by timepoint and CDR3
nonnaiveCD8_Seurat_subset <- ScaleData(nonnaiveCD8_Seurat_subset, vars.to.regress = c("percent.mito"))
nonnaiveCD8_Seurat_subset <- RunPCA(nonnaiveCD8_Seurat_subset, features = VariableFeatures(object = nonnaiveCD8_Seurat_subset))

dims_use = 1:20
nonnaiveCD8_Seurat_subset <- RunHarmony(nonnaiveCD8_Seurat_subset, group.by.vars=c("timepoint", "CTaa"))
nonnaiveCD8_Seurat_subset <- RunUMAP(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindNeighbors(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindClusters(object=nonnaiveCD8_Seurat_subset, resolution=0.15, verbose=FALSE)

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = TRUE, logfc.threshold = 0.50, min.pct = 0.30)

a <- DimPlot(nonnaiveCD8_Seurat_subset, label = T, repel = T) +
  scale_color_brewer(palette = "Paired") + xlab("UMAP 1") + ylab("UMAP 2")

b <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "timepoint", shuffle = T, cols = timepoint_cols) + xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(NULL) 

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$beta
c <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[1]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[2]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[3]),
                                                               WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[4])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Group_2", "Group_3", "Group_4", "Unselected"), 
                     labels = c("Clone #1", "Clone #2", "Clone #3", "Clone #4", "Other"),
                     values = c(brewer.pal(4, "Set2"), "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$seurat_clusters

d <- FeaturePlot(nonnaiveCD8_Seurat_subset, features = "GZMB") + 
  scale_color_distiller(palette = "Reds", direction = 1, name = "GZMB") +
  xlab("UMAP 1") + ylab("UMAP 2") + ggtitle(NULL)

wrap_plots(a, b, c, d)

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$beta
e <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[1])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #1", "Other"),
                     values = c(brewer.pal(4, "Set2")[1], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

f <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[2])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #2", "Other"),
                     values = c(brewer.pal(4, "Set2")[2], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

g <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[3])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #3", "Other"),
                     values = c(brewer.pal(4, "Set2")[3], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

h <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_topclone_list$beta[4])), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Clone #4", "Other"),
                     values = c(brewer.pal(4, "Set2")[4], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2")

wrap_plots(e, f, g, h)
Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$seurat_clusters


#dimensionality reduction with removing VDJ as variable, percent mito regress, and harmonized by timepoint and CDR3 (reduced theta)
nonnaiveCD8_Seurat_subset <- ScaleData(nonnaiveCD8_Seurat_subset, vars.to.regress = c("percent.mito"))
nonnaiveCD8_Seurat_subset <- RunPCA(nonnaiveCD8_Seurat_subset, features = VariableFeatures(object = nonnaiveCD8_Seurat_subset))

dims_use = 1:20
nonnaiveCD8_Seurat_subset <- RunHarmony(nonnaiveCD8_Seurat_subset, group.by.vars=c("timepoint", "CTaa"), theta = c(2,0.1))
nonnaiveCD8_Seurat_subset <- RunUMAP(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindNeighbors(object=nonnaiveCD8_Seurat_subset, reduction="harmony", dims=dims_use, verbose=FALSE)
nonnaiveCD8_Seurat_subset <- FindClusters(object=nonnaiveCD8_Seurat_subset, resolution=0.2, verbose=FALSE)

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = FALSE, logfc.threshold = 0.25, min.pct = 0.20)

StackedVlnPlot(nonnaiveCD8_Seurat_subset, features = c("GZMB", "FGFBP2", "GNLY", "PRF1", "GZMK", "CD69", "HLA-DRA", "CCR7", "LEF1", "TRAV1-2", "SLC4A10", "CXCR3", "MALAT1", "MKI67"), cols = brewer.pal(10, "Paired"))

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$seurat_clusters
current.cluster.names <- c(0:6)
new.cluster.names <- c("CTL", "EM", "CM", "MAIT", "ZNF683+", "Mito-hi", "Cycling")
nonnaiveCD8_Seurat_subset@active.ident <- plyr::mapvalues(x = nonnaiveCD8_Seurat_subset@active.ident, from = current.cluster.names, to = new.cluster.names)
nonnaiveCD8_Seurat_subset$clusters_names <- nonnaiveCD8_Seurat_subset@active.ident

slot(nonnaiveCD8_Seurat_subset, "meta.data")$clusters_names <- factor(slot(nonnaiveCD8_Seurat_subset, "meta.data")$clusters_names, 
                                                                      levels = c("CTL", "EM", "ZNF683+", "CM", "MAIT", "Cycling", "Mito-hi"))
Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$clusters_names

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = FALSE, logfc.threshold = 0.25, min.pct = 0.20)
write.csv(nonNaiveCD8_cluster_markers, "analysis/data/sc_nonNaive_RNA_cluster_markers.csv")

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = TRUE, logfc.threshold = 0.50, min.pct = 0.30)

StackedVlnPlot(nonnaiveCD8_Seurat_subset, features = c("GZMB", "FGFBP2", "GNLY", "PRF1", "GZMK", "CD69", "HLA-DRA", "ZNF683", "CXCR3", "LEF1", "CCR7", "TRAV1-2", "SLC4A10", "MKI67", "MALAT1"), cols = brewer.pal(10, "Paired"))
ggsave("analysis/figs/sc_nonnaiveCD8_stackedVln_Clusters.eps", width = 4, height = 5.5)

top10nonnaiveCD8_deg_pos <- nonNaiveCD8_cluster_markers %>% group_by(cluster) %>% slice_head(n=10)
DoHeatmap(nonnaiveCD8_Seurat_subset, features = c(top10nonnaiveCD8_deg_pos$gene[1:20], "ZNF683", top10nonnaiveCD8_deg_pos$gene[21:70]), group.colors = brewer.pal(10, "Paired"), size = 4) + 
  scale_fill_distiller(palette = "RdBu") + 
  theme(axis.text.y = element_text(size = 6))
ggsave("analysis/figs/sc_nonnaiveCD8_heatmap_Clusters.eps", width = 7, height = 7)

cd8cols <- c(colorRampPalette(brewer.pal(10, "Paired"))(10))

nonnaiveclusterUMAP <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "clusters_names") + 
  scale_color_manual(values = cd8cols) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

nonnaiveclusterUMAPPlabeled <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "clusters_names", label = T, repel = T, label.box = T) + 
  scale_color_manual(values = cd8cols) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = "CASSQDPGPGNTIYF"),
                                                          WhichCells(nonnaiveCD8_Seurat_subset, idents = "CASSLGRQTPTGELFF"),
                                                          WhichCells(nonnaiveCD8_Seurat_subset, idents = "CASSLTFASGGLRRSTDTQYF"),
                                                          WhichCells(nonnaiveCD8_Seurat_subset, idents = "CASSLYGQGGSNQPQHF")), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Group_2", "Group_3", "Group_4", "Unselected"), 
                     labels = c("Clone #1", "Clone #2", "Clone #3", "Clone #4", "Other"),
                     values = c(brewer.pal(4, "Set2"), "gray90"))

nonnaiveCD8_Seurat_subset@meta.data %>%
  ggplot(aes(x = clusters_names, y = percent.mito, fill = clusters_names)) +
  geom_boxplot(outlier.size = 0.1) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("Mito Reads (%)") + 
  scale_fill_manual(values = cd8cols) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

nonNaiveCD8_cluster_markers <- FindAllMarkers(object=nonnaiveCD8_Seurat_subset, only.pos = TRUE, logfc.threshold = 0.4, min.pct = 0.30)
dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015", "KEGG_2021_Human", "Azimuth_Cell_Types_2021")

enriched_CM <- enrichr(c(nonNaiveCD8_cluster_markers$gene[grep("CM", nonNaiveCD8_cluster_markers$cluster)]), dbs)
enriched_CM <- lapply(enriched_CM, function(x) 
  cbind(x, Cluster = "CM"))
enriched_EM <- enrichr(c(nonNaiveCD8_cluster_markers$gene[grep("EM", nonNaiveCD8_cluster_markers$cluster)]), dbs)
enriched_EM <- lapply(enriched_EM, function(x) 
  cbind(x, Cluster = "EM"))
enriched_CXCR3 <- enrichr(c(nonNaiveCD8_cluster_markers$gene[grep("CXCR3+", nonNaiveCD8_cluster_markers$cluster)]), dbs)
enriched_CXCR3 <- lapply(enriched_CXCR3, function(x) 
  cbind(x, Cluster = "CXCR3+"))
enriched_CTL <- enrichr(c(nonNaiveCD8_cluster_markers$gene[grep("CTL", nonNaiveCD8_cluster_markers$cluster)]), dbs)
enriched_CTL <- lapply(enriched_CTL, function(x) 
  cbind(x, Cluster = "CTL"))
enriched_Cycling <- enrichr(c(nonNaiveCD8_cluster_markers$gene[grep("Cycling", nonNaiveCD8_cluster_markers$cluster)]), dbs)
enriched_Cycling <- lapply(enriched_Cycling, function(x) 
  cbind(x, Cluster = "Cycling"))
enriched_Mitohi <- enrichr(c(nonNaiveCD8_cluster_markers$gene[grep("Mito-hi", nonNaiveCD8_cluster_markers$cluster)]), dbs)
enriched_Mitohi <- lapply(enriched_Mitohi, function(x) 
  cbind(x, Cluster = "Mito-hi"))

GObiologicalprocess_top5 <- rbind(enriched_CTL$GO_Biological_Process_2015[1:5,], 
                                  enriched_EM$GO_Biological_Process_2015[1:5,], enriched_CXCR3$GO_Biological_Process_2015[1:5,],
                                  enriched_CM$GO_Biological_Process_2015[1:5,], enriched_Cycling$GO_Biological_Process_2015[1:5,],
                                  enriched_Mitohi$GO_Biological_Process_2015[1:5,])
GObiologicalprocess_top5$TermCluster <- paste0(GObiologicalprocess_top5$Term, GObiologicalprocess_top5$Cluster)
GObiologicalprocess_top5$TermCluster <- factor(GObiologicalprocess_top5$TermCluster,levels = c(GObiologicalprocess_top5$TermCluster))

ggplot(GObiologicalprocess_top5, aes(x=TermCluster, y=-log(Adjusted.P.value))) +
  geom_segment(aes(x=TermCluster, xend=TermCluster, y=0, yend=-log(Adjusted.P.value), color = Cluster)) +
  geom_point( aes(color=Cluster), size=4) +
  scale_color_manual(values = cd8cols) + 
  scale_x_discrete(labels = rev(GObiologicalprocess_top5$Term), limits = rev(levels(GObiologicalprocess_top5$TermCluster))) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  xlab(NULL) + ylab("-log(Adjusted P Value)") +
  theme_classic() +
  coord_flip()


#Symphony
#Watson et al 2021  
watson_seurat_cd8_blood <- readRDS("nonNaiveCD8_scRNA/data/other_datasets/watson_seurat_cd8_blood.rds")
Idents(watson_seurat_cd8_blood) <- watson_seurat_cd8_blood$assignment_final
watson_seurat_cd8_blood <- subset(watson_seurat_cd8_blood, idents =  c("Unknown", "gd", "Naive"), invert = T)
watson_seurat_cd8_blood$assignment_final <- Idents(watson_seurat_cd8_blood)

watson_seurat_cd8_blood$assignment_final <- factor(watson_seurat_cd8_blood$assignment_final, levels = c("Effector", "EM", "CM", "MAIT", "Mitotic"))
Idents(watson_seurat_cd8_blood) <- watson_seurat_cd8_blood$assignment_final

watson_seurat_cd8_blood <- NormalizeData(watson_seurat_cd8_blood, normalization.method = "LogNormalize", scale.factor = 10000)

#Symphony on repertoire T cell data
symphT_reference <- symphony::buildReference(
  watson_seurat_cd8_blood@assays$RNA@data,
  watson_seurat_cd8_blood@meta.data,
  vars = c('donor'),         # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = FALSE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'new_cluster_name', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20                     # number of PCs
)

queryT <- mapQuery(nonnaiveCD8_Seurat_subset@assays$RNA@data,             # query gene expression (genes x cells)
                   nonnaiveCD8_Seurat_subset@meta.data,        # query metadata (cells x attributes)
                   symphT_reference,             # Symphony reference object
                   vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                   do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                   do_umap = FALSE)        # project query cells into reference UMAP

queryT <- knnPredict(queryT, symphT_reference, train_labels = symphT_reference$meta_data$assignment_final,
                     k = 10, save_as = 'cell_type_10_nn', confidence = TRUE)

nonnaiveCD8_Seurat_subset$symphony_narrow_type_watson <- queryT$meta_data$cell_type_10_nn

# Cell type mappings heatmap
query_immune = queryT$meta_data
query_immune$symphony_narrow_type = droplevels(as.factor(query_immune$cell_type_10_nn))
res_immune = symphony:::evaluate(query_immune$clusters_names, query_immune$cell_type_10_nn)
Conf_immune_watson = res_immune$Conf / rowSums(res_immune$Conf)

dev.off()
#setEPS(width = 6.5, height = 7.5)
#postscript("CD8_watsonetal_symphony.eps")
pheatmap::pheatmap(Conf_immune_watson, cluster_rows = F, cluster_cols = F, angle_col = c("45"), cellwidth = 20, cellheight = 20, color=colorRampPalette(c("white", "#962f23", "#5e0e04"))(1000), border_color = "white")
dev.off()


#pauken et al 2021
pauken_seurat_cd8_blood <- readRDS("nonNaiveCD8_scRNA/data/other_datasets/pauken_seurat_cd8_blood.rds")

DefaultAssay(pauken_seurat_cd8_blood) <- "RNA"

Idents(pauken_seurat_cd8_blood) <- pauken_seurat_cd8_blood$Seurat_clusters

current.clust <- c(0:8)
new.clust <- c("Naive/CM", "Effector/EM/Exhausted Precursor", "Effector/EM", "CM", "Effector/EM", "Naive/CM", "Naive/CM", "CM", "Exhausted Precursor")
pauken_seurat_cd8_blood@active.ident <- plyr::mapvalues(x = pauken_seurat_cd8_blood@active.ident, from = current.clust, to = new.clust)
pauken_seurat_cd8_blood$cluster_mappings <- pauken_seurat_cd8_blood@active.ident
Idents(pauken_seurat_cd8_blood) <- pauken_seurat_cd8_blood$cluster_mappings

pauken_seurat_cd8_blood$cluster_mappings <- factor(pauken_seurat_cd8_blood$cluster_mappings, levels = c("Exhausted Precursor", "Effector/EM/Exhausted Precursor", "Effector/EM", "CM", "Naive/CM" ))
Idents(pauken_seurat_cd8_blood) <- pauken_seurat_cd8_blood$cluster_mappings

pauken_seurat_cd8_blood <- NormalizeData(pauken_seurat_cd8_blood, normalization.method = "LogNormalize", scale.factor = 10000)

#Symphony on repertoire T cell data
symphT_reference <- symphony::buildReference(
  pauken_seurat_cd8_blood@assays$RNA@data,
  pauken_seurat_cd8_blood@meta.data,
  vars = c('Sample'),         # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = FALSE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'new_cluster_name', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20                     # number of PCs
)

queryT <- mapQuery(nonnaiveCD8_Seurat_subset@assays$RNA@data,             # query gene expression (genes x cells)
                   nonnaiveCD8_Seurat_subset@meta.data,        # query metadata (cells x attributes)
                   symphT_reference,             # Symphony reference object
                   vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                   do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                   do_umap = FALSE)        # project query cells into reference UMAP

queryT <- knnPredict(queryT, symphT_reference, train_labels = symphT_reference$meta_data$cluster_mappings,
                     k = 10, save_as = 'cell_type_10_nn', confidence = TRUE)

nonnaiveCD8_Seurat_subset$symphony_narrow_type_pauken <- queryT$meta_data$cell_type_10_nn

# Cell type mappings heatmap
query_immune = queryT$meta_data
query_immune$symphony_narrow_type = droplevels(as.factor(query_immune$cell_type_10_nn))
res_immune = symphony:::evaluate(query_immune$clusters_names, query_immune$cell_type_10_nn)
Conf_immune_pauken = res_immune$Conf / rowSums(res_immune$Conf)

pheatmap::pheatmap(Conf_immune_pauken, cluster_rows = F, cluster_cols = F, angle_col = c("45"), cellwidth = 20, cellheight = 20, color=colorRampPalette(c("white", "#962f23", "#5e0e04"))(1000), border_color = "white")


#luoma et al 2022
luoma_seurat_cd8_blood <- readRDS("nonNaiveCD8_scRNA/data/other_datasets/luoma_blood//luoma_seurat_cd8_blood.rds")

Idents(luoma_seurat_cd8_blood) <- luoma_seurat_cd8_blood$CellType_ID

luoma_seurat_cd8_blood <- subset(luoma_seurat_cd8_blood, idents = 7, invert = T)

current.clust <- c(1:6,8:9)
new.clust <- c("GZMK+", "GZMBhi", "CCR7+", "FGFBP2hi", "IL7R+", "CD38+", "LTB+", "KLRB1+")
luoma_seurat_cd8_blood@active.ident <- plyr::mapvalues(x = luoma_seurat_cd8_blood@active.ident, from = current.clust, to = new.clust)
luoma_seurat_cd8_blood$cluster_mappings <- luoma_seurat_cd8_blood@active.ident
Idents(luoma_seurat_cd8_blood) <- luoma_seurat_cd8_blood$cluster_mappings

luoma_seurat_cd8_blood$cluster_mappings <- factor(luoma_seurat_cd8_blood$cluster_mappings, levels = c("FGFBP2hi", "GZMBhi", "GZMK+", "IL7R+", "CCR7+", "LTB+", "KLRB1+", "CD38+"))
Idents(luoma_seurat_cd8_blood) <- luoma_seurat_cd8_blood$cluster_mappings

luoma_seurat_cd8_blood <- NormalizeData(luoma_seurat_cd8_blood, normalization.method = "LogNormalize", scale.factor = 10000)
luoma_seurat_cd8_blood$Patient_Stage <- paste0(luoma_seurat_cd8_blood$Patient_ID, "-", luoma_seurat_cd8_blood$Stage)

#Symphony on repertoire T cell data
symphT_reference <- symphony::buildReference(
  luoma_seurat_cd8_blood@assays$RNA@data,
  luoma_seurat_cd8_blood@meta.data,
  vars = c('Patient_Stage'),         # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = FALSE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'new_cluster_name', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20                     # number of PCs
)

queryT <- mapQuery(nonnaiveCD8_Seurat_subset@assays$RNA@data,             # query gene expression (genes x cells)
                   nonnaiveCD8_Seurat_subset@meta.data,        # query metadata (cells x attributes)
                   symphT_reference,             # Symphony reference object
                   vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                   do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                   do_umap = FALSE)        # project query cells into reference UMAP

queryT <- knnPredict(queryT, symphT_reference, train_labels = symphT_reference$meta_data$cluster_mappings,
                     k = 10, save_as = 'cell_type_10_nn', confidence = TRUE)

nonnaiveCD8_Seurat_subset$symphony_narrow_type_luoma <- queryT$meta_data$cell_type_10_nn

# Cell type mappings heatmap
query_immune = queryT$meta_data
query_immune$symphony_narrow_type = droplevels(as.factor(query_immune$cell_type_10_nn))
res_immune = symphony:::evaluate(query_immune$clusters_names, query_immune$cell_type_10_nn)
Conf_immune_luoma = res_immune$Conf / rowSums(res_immune$Conf)

pheatmap::pheatmap(Conf_immune_luoma, cluster_rows = F, cluster_cols = F, angle_col = c("45"), cellwidth = 20, cellheight = 20, color=colorRampPalette(c("white", "#962f23", "#5e0e04"))(1000), border_color = "white")


#li et al 2022
li_seurat_cd8_blood <- readRDS("nonNaiveCD8_scRNA/data/other_datasets/li_datasets/li_seurat_CD8_blood.rds")

li_seurat_cd8_blood$clusters_names <- Idents(li_seurat_cd8_blood)
DefaultAssay(li_seurat_cd8_blood) <- "RNA"

Idents(li_seurat_cd8_blood) <- li_seurat_cd8_blood$clusters_names

li_seurat_cd8_blood <- subset(li_seurat_cd8_blood, idents = "naive CD8" , invert = T)

current.clust <- c("effector CD8", "GZMK+ CD8", "memory CD8",  "KIR+ effector CD8", "IFN-stim CD8", "MAIT", "proliferating CD8")
new.clust <- c("Effector", "GZMK+", "Memory", "KIR+", "IFN-Stim", "MAIT", "Proliferating")
li_seurat_cd8_blood@active.ident <- plyr::mapvalues(x = li_seurat_cd8_blood@active.ident, from = current.clust, to = new.clust)
li_seurat_cd8_blood$clusters_names <- li_seurat_cd8_blood@active.ident
Idents(li_seurat_cd8_blood) <- li_seurat_cd8_blood$clusters_names

li_seurat_cd8_blood$clusters_names <- factor(li_seurat_cd8_blood$clusters_names, levels = c("Effector", "GZMK+", "Memory", "KIR+", "IFN-Stim", "MAIT", "Proliferating"))
Idents(li_seurat_cd8_blood) <- li_seurat_cd8_blood$clusters_names

li_seurat_cd8_blood <- NormalizeData(li_seurat_cd8_blood, normalization.method = "LogNormalize", scale.factor = 10000)

#Symphony on repertoire T cell data
symphT_reference <- symphony::buildReference(
  li_seurat_cd8_blood@assays$RNA@data,
  li_seurat_cd8_blood@meta.data,
  vars = c('patient'),         # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = FALSE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'new_cluster_name', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20                     # number of PCs
)

queryT <- mapQuery(nonnaiveCD8_Seurat_subset@assays$RNA@data,             # query gene expression (genes x cells)
                   nonnaiveCD8_Seurat_subset@meta.data,        # query metadata (cells x attributes)
                   symphT_reference,             # Symphony reference object
                   vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                   do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                   do_umap = FALSE)        # project query cells into reference UMAP

queryT <- knnPredict(queryT, symphT_reference, train_labels = symphT_reference$meta_data$clusters_names,
                     k = 10, save_as = 'cell_type_10_nn', confidence = TRUE)

nonnaiveCD8_Seurat_subset$symphony_narrow_type_li <- queryT$meta_data$cell_type_10_nn

# Cell type mappings heatmap
query_immune = queryT$meta_data
query_immune$symphony_narrow_type = droplevels(as.factor(query_immune$cell_type_10_nn))
res_immune = symphony:::evaluate(query_immune$clusters_names, query_immune$cell_type_10_nn)
Conf_immune_li = res_immune$Conf / rowSums(res_immune$Conf)

pheatmap::pheatmap(Conf_immune_li, cluster_rows = F, cluster_cols = F, angle_col = c("45"), cellwidth = 20, cellheight = 20, color=colorRampPalette(c("white", "#962f23", "#5e0e04"))(1000), border_color = "white")

dev.off()
setEPS(width = 9, height = 4)
postscript("analysis/figs/sc_nonnaiveCD8_symphony_ref.eps")
pheatmap::pheatmap(cbind(Conf_immune_watson, Conf_immune_luoma, Conf_immune_li), cluster_rows = F, cluster_cols = F, angle_col = c("45"), cellwidth = 20, cellheight = 20, 
                   color=colorRampPalette(c("white", "#962f23", "#5e0e04"))(1000),border_color = "white", gaps_col = c(5,13))
dev.off()

nonnaiveCD8_Seurat_subset@meta.data %>%
  group_by(clusters_names, timepoint) %>%
  dplyr::count() %>%
  group_by(timepoint) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill = timepoint)) +
  geom_bar(position="dodge", stat="identity") +
  scale_fill_manual(values = timepoint_cols, name = NULL) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  ylab("Percent cluster of\nnon-naive CD8 T cells") +
  xlab(NULL) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size = 9, angle = 45, hjust = 1),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10)) 
ggsave("analysis/figs/sc_nonnaiveCD8_Clusters_byTimepoint_split.eps", width = 4, height = 2.5)

nonnaiveCD8_Seurat_subset@meta.data %>%
  group_by(clusters_names, timepoint) %>%
  dplyr::count() %>%
  group_by(timepoint) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=timepoint,y=Percent, fill=clusters_names)) +
  geom_col()+
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_x_discrete(limits = c("Pembro2", "Pembro1"), labels = c("Pembro 2", "Pembro 1")) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm'),
        legend.position = "none") + coord_flip()
ggsave("analysis/figs/sc_nonnaiveCD8_Clusters_byTimepoint.eps", width = 3.5, height = 1)

nonnaivetimepointUMAP <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "timepoint", shuffle = T, cols = timepoint_cols) + 
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

DimPlot(nonnaiveCD8_Seurat_subset, group.by = "timepoint", split.by = "timepoint", shuffle = T, cols = timepoint_cols, ncol = 1) + 
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
ggsave("analysis/figs/sc_nonnaiveCD8_UMAP_splitbyTimepoint.eps", width = 2.5, height = 5)

nonnaiveCD8_Seurat_subset$timepoint_cluster <- paste0(nonnaiveCD8_Seurat_subset$timepoint, "_", nonnaiveCD8_Seurat_subset$clusters_names)
Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$timepoint_cluster

nonnaivepembro1UMAP <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_CTL"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_EM"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_ZNF683+"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_CM"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_MAIT"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_Cycling"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro1_Mito-hi")), 
                               cols.highlight = rev(cd8cols), sizes.highlight = 0.05, cols = "gray80") +
  scale_color_manual(breaks = c("Group_1", "Group_2", "Group_3", "Group_4", "Group_5", "Group_6", "Group_7", "Unselected"), 
                     labels = c("CTL", "EM", "ZNF683+", "CM", "MAIT", "Cycling", "Mito-hi", "Pembro2"), 
                     values = c(cd8cols[1:7], "gray90")) +
  ggtitle("Pembro 1") + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(legend.title = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12))


nonnaivepembro2UMAP <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = list(WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_CTL"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_EM"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_ZNF683+"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_CM"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_MAIT"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_Cycling"),
                                                                                 WhichCells(nonnaiveCD8_Seurat_subset, idents = "Pembro2_Mito-hi")), 
                               cols.highlight = rev(cd8cols), sizes.highlight = 0.05, cols = "gray80") +
  scale_color_manual(breaks = c("Group_1", "Group_2", "Group_3", "Group_4", "Group_5", "Group_6", "Group_7", "Unselected"), 
                     labels = c("CTL", "EM", "ZNF683+", "CM", "MAIT", "Cycling", "Mito-hi", "Pembro2"), 
                     values = c(cd8cols[1:7], "gray90")) +
  ggtitle("Pembro 2") + xlab("UMAP 1") + ylab("UMAP 2") + 
  theme(legend.title = element_text(size = 10),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12))


#clonality/clonal comparisons
nonnaive_cloneType_byCluster_table <- table(nonnaiveCD8_Seurat_subset$cloneType, nonnaiveCD8_Seurat_subset$clusters_names)

nonnaiveclonalityUMAP <- DimPlot(nonnaiveCD8_Seurat_subset, group.by = "cloneType") + 
  scale_color_brewer(palette = "Blues", direction = -1) +
  ggtitle(NULL) +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank())

nonnaiveCD8_Seurat_subset_cluster_clones <- expression2List(nonnaiveCD8_Seurat_subset, group = "clusters_names")

nonnaive_clust_morisita <- clonalOverlap(nonnaiveCD8_Seurat_subset_cluster_clones, cloneCall = "aa", method = "morisita", exportTable = T)
nonnaive_clust_morisita$names <- NULL
nonnaive_clust_morisita <- as.matrix(nonnaive_clust_morisita)
nonnaive_clust_morisita <- Matrix::forceSymmetric(nonnaive_clust_morisita,uplo="U")
nonnaive_clust_morisita <- as.matrix(nonnaive_clust_morisita)
class(nonnaive_clust_morisita) <- c("immunr_ov_matrix", "matrix", "array")
vis(nonnaive_clust_morisita) + 
  scale_y_discrete(limits = rev(levels(nonnaiveCD8_Seurat_subset$clusters_names)), labels=rev(levels(nonnaiveCD8_Seurat_subset$clusters_names))) + 
  scale_x_discrete(limits = levels(nonnaiveCD8_Seurat_subset$clusters_names), labels=levels(nonnaiveCD8_Seurat_subset$clusters_names)) + 
  theme(axis.text.x=element_text(angle=90,hjust=0.95,vjust=0.2)) +
  guides(fill = guide_colorbar(title="Morisita \nIndex")) +
  ggtitle(NULL) + xlab(NULL) + ylab(NULL) +
  scale_fill_distiller(palette = "Blues", direction = 1) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1),
        axis.text=element_text(size=12)) 
ggsave("analysis/figs/sc_nonnaiveCD8_cluster_overlap_morisita.eps", width = 7, height = 7)

nonnaiveCD8_Seurat_subset@meta.data %>%
  group_by(clusters_names, cloneType) %>%
  count() %>%
  group_by(clusters_names) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=cloneType)) +
  geom_col(color = "black", size = 0.5) +
  scale_fill_brewer(palette = "Blues", direction = -1) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) +
  xlab(NULL) +
  ylab(NULL) +
  labs(fill=NULL) +
  theme_classic() +
  theme(axis.text=element_text(size=12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size=11))



#comparing clone proportions between bulkTCR and non-naive scTCR
bulkpembro1clonefreqs <- bulktcr_cdr3_freqs %>% select(CDR3aa, pembro_1)
colnames(bulkpembro1clonefreqs) <- c("beta", "pembro1_bulk")

nonnaiveCD8pembro1clonefreqs <- nonnaiveCD8_Seurat_subset@meta.data %>% filter(timepoint == "Pembro1") %>% filter(beta != "NA") %>%
  group_by(beta) %>% count() %>% ungroup()  %>% mutate(freq = n/sum(n)) %>% select(beta, freq)
colnames(nonnaiveCD8pembro1clonefreqs) <- c("beta", "pembro1_sc")

left_join(bulkpembro1clonefreqs, nonnaiveCD8pembro1clonefreqs, by = "beta") %>% filter(pembro1_bulk > 0 & pembro1_sc > 0)  %>% ggplot(aes(x = pembro1_bulk, y = pembro1_sc)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = bulktcrcolors[4], linetype = "dashed") + 
  stat_cor(label.y.npc = "top", label.x.npc = 0.2, size = 4) +
  xlab("Pembro 1 Bulk TCRseq") + ylab("Pembro 1 Non-Naive CD8 TCRseq") +
  theme_classic() 

bulkpembro2clonefreqs <- bulktcr_cdr3_freqs %>% select(CDR3aa, pembro_2)
colnames(bulkpembro2clonefreqs) <- c("beta", "pembro2_bulk")

nonnaiveCD8pembro2clonefreqs <- nonnaiveCD8_Seurat_subset@meta.data %>% filter(timepoint == "Pembro2") %>% filter(beta != "NA") %>%
  group_by(beta) %>% count() %>% ungroup()  %>% mutate(freq = n/sum(n)) %>% select(beta, freq)
colnames(nonnaiveCD8pembro2clonefreqs) <- c("beta", "pembro2_sc")

left_join(bulkpembro2clonefreqs, nonnaiveCD8pembro2clonefreqs, by = "beta") %>% filter(pembro2_bulk > 0 & pembro2_sc > 0)  %>% ggplot(aes(x = pembro2_bulk, y = pembro2_sc)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = bulktcrcolors[4], linetype = "dashed") + 
  stat_cor(label.y.npc = "top", label.x.npc = 0.2, size = 4) +
  xlab("Pembro 2 Bulk TCRseq") + ylab("Pembro 2 Non-Naive CD8 TCRseq") +
  theme_classic() 


#clone matching
Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$barcode

#kidney bulk TCR
nonnaive_kidney_clones <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa) 
DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_kidney_clones$barcode), sizes.highlight = 0.1) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[1])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Kidney\n Matched Clones")

nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(kidney = ifelse(beta %in% bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa, "Kidney", "No")) %>%
  na.omit(beta) %>%
  distinct(beta, clusters_names, .keep_all=T) %>%
  group_by(clusters_names, kidney) %>%
  dplyr::count() %>% 
  group_by(clusters_names) %>% mutate(Percent=100*n/sum(n)) %>%
  filter(kidney == "Kidney") %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) + xlab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10)) 

nonnaive_kidney_match_percent <- nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(tissue = ifelse(beta %in% bulktcr_kid$data$`GM-18-18046_KID_TCRB`$CDR3.aa, "Kidney", "No")) %>%
  na.omit(beta) %>%
  distinct(beta, clusters_names, .keep_all=T) %>%
  group_by(clusters_names, tissue) %>%
  dplyr::count() %>% 
  group_by(clusters_names) %>% mutate(Percent=100*n/sum(n)) %>%
  filter(tissue == "Kidney") %>%
  ungroup() %>% select(clusters_names, tissue, Percent) 

#ln bulk TCR
nonnaive_ln_clones <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa) 
DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_ln_clones$barcode), sizes.highlight = 0.1) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[2])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("LN\n Matched Clones")

nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(LN = ifelse(beta %in% bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa, "LN", "No")) %>%
  na.omit(beta) %>%
  distinct(beta, clusters_names, .keep_all=T) %>%
  group_by(clusters_names, LN) %>%
  dplyr::count() %>% 
  group_by(clusters_names) %>% mutate(Percent=100*n/sum(n)) %>%
  filter(LN == "LN") %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) + xlab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10)) 

nonnaive_ln_match_percent <- nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(tissue = ifelse(beta %in% bulktcr_kid$data$`GM-17-7547_pre_LN_TCRB`$CDR3.aa, "LN", "No")) %>%
  na.omit(beta) %>%
  distinct(beta, clusters_names, .keep_all=T) %>%
  group_by(clusters_names, tissue) %>%
  dplyr::count() %>% 
  group_by(clusters_names) %>% mutate(Percent=100*n/sum(n)) %>%
  filter(tissue == "LN") %>%
  ungroup()  %>% select(clusters_names, tissue, Percent)

#skin bulk TCR
nonnaive_skin_clones <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa) 
DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = nonnaive_skin_clones$barcode), sizes.highlight = 0.1) +
  scale_color_manual(values = c("#ECECEC", bulktcrcolors[3])) +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  ggtitle("Skin\n Matched Clones")

nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(skin = ifelse(beta %in% bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa, "Skin", "No")) %>%
  na.omit(beta) %>%
  distinct(beta, clusters_names, .keep_all=T) %>%
  group_by(clusters_names, skin) %>%
  dplyr::count() %>% 
  group_by(clusters_names) %>% mutate(Percent=100*n/sum(n)) %>%
  filter(skin == "Skin") %>%
  ungroup() %>%
  ggplot(aes(x=clusters_names,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  ggtitle(NULL) + xlab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10)) 

nonnaive_skin_match_percent <- nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(tissue = ifelse(beta %in% bulktcr_kid$data$`GM-19-7655_post_TUM_TCRB`$CDR3.aa, "Skin", "No")) %>%
  na.omit(beta) %>%
  distinct(beta, clusters_names, .keep_all=T) %>%
  group_by(clusters_names, tissue) %>%
  dplyr::count() %>% 
  group_by(clusters_names) %>% mutate(Percent=100*n/sum(n)) %>%
  filter(tissue == "Skin") %>%
  ungroup() %>% select(clusters_names, tissue, Percent)

rbind(nonnaive_kidney_match_percent, nonnaive_ln_match_percent, nonnaive_skin_match_percent) %>%
  ggplot(aes(x = clusters_names, y = Percent, fill = tissue)) +
  geom_col(position = "dodge2") +
  scale_fill_manual(values = bulktcrcolors[1:3], name = "Tissue") +
  theme_classic() + xlab(NULL) + ylab("Percent Matching Clones") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#MLR matching
c1clones <- MLR_subset_metadata %>% filter(seurat_clusters == "C1") %>% separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% na.omit(beta) %>% distinct(beta)
allo_clones_nonnaive <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% c1clones$beta) 

c1matchumap <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = allo_clones_nonnaive$barcode), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("C1 Matching Clones", "Other"),
                     values = c(singlecellcolors[1], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("MLR C1 (Alloreactive)\n Matched Clones") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 11))

allo_clones_nonnaive %>%
  group_by(clusters_names, timepoint) %>%
  dplyr::count() %>%
  group_by(timepoint) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=timepoint,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(limits = c("Pembro2", "Pembro1")) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm'),
        legend.position = "none") + coord_flip()
ggsave("analysis/figs/sc_nonnaiveCD8_C1match_byTimepoint.eps", width = 3.5, height = 1)

c2clones <- MLR_subset_metadata %>% filter(seurat_clusters == "C2") %>% separate(CTaa, into = c("alpha", "beta"), sep="_(?=[^_]+$)") %>% na.omit(beta) %>% distinct(beta)
c2_clones_nonnaive <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% c2clones$beta) %>% filter(beta != "NA")

c2matchumap <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = na.omit(c2_clones_nonnaive$barcode)), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("C2 Clones", "Other"),
                     values = c(singlecellcolors[2], "gray90")) + xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("MLR C2\n Matched Clones") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 11))

c2_clones_nonnaive %>%
  group_by(clusters_names, timepoint) %>%
  dplyr::count() %>%
  group_by(timepoint) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=timepoint,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(limits = c("Pembro2", "Pembro1")) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm'),
        legend.position = "none") + coord_flip()
ggsave("analysis/figs/sc_nonnaiveCD8_C2match_byTimepoint.eps", width = 3.5, height = 1)

nonnaive_UMAP_MLRmatch <- align_plots(c1matchumap, c2matchumap, align="hv", axis="tblr")
p1x <- ggdraw(nonnaive_UMAP_MLRmatch[[1]])
p2x <- ggdraw(nonnaive_UMAP_MLRmatch[[2]])

save_plot("analysis/figs/sc_nonnaive_UMAP_MLRc1_match.eps", p1x, base_width = 3.5, base_height = 3.7)
save_plot("analysis/figs/sc_nonnaive_UMAP_MLRc2_match.eps", p2x, base_width = 3.5, base_height = 3.7)

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$clusters_names

nonnaiveCD8_Seurat_subset$MLRmatch <- ifelse(nonnaiveCD8_Seurat_subset$barcode %in% allo_clones_nonnaive$barcode, "C1", 
                                             ifelse(nonnaiveCD8_Seurat_subset$barcode %in% c2_clones_nonnaive$barcode, "C2", "No Match"))
nonnaiveCD8_Seurat_subset$MLRmatchcluster <- paste0(nonnaiveCD8_Seurat_subset$clusters_names, "_", nonnaiveCD8_Seurat_subset$MLRmatch)

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$clusters_names
nonnaiveCD8_Seurat_subset_cxcr3 <- subset(nonnaiveCD8_Seurat_subset, idents = c("ZNF683+"))

Idents(nonnaiveCD8_Seurat_subset_cxcr3) <- nonnaiveCD8_Seurat_subset_cxcr3$MLRmatch

current.clone.names <- c("C1", "C2", "No Match")
new.clone.names <- c("C1 Match", "C2 Match", "No Match")
nonnaiveCD8_Seurat_subset_cxcr3@active.ident <- plyr::mapvalues(x = nonnaiveCD8_Seurat_subset_cxcr3@active.ident, from = current.clone.names, to = new.clone.names)
nonnaiveCD8_Seurat_subset_cxcr3$MLRmatch <- nonnaiveCD8_Seurat_subset_cxcr3@active.ident

slot(nonnaiveCD8_Seurat_subset_cxcr3, "meta.data")$MLRmatch <- factor(slot(nonnaiveCD8_Seurat_subset_cxcr3, "meta.data")$MLRmatch, 
                                                                      levels = c("C1 Match", "C2 Match", "No Match"))
Idents(nonnaiveCD8_Seurat_subset_cxcr3) <- nonnaiveCD8_Seurat_subset_cxcr3$MLRmatch

StackedVlnPlot(nonnaiveCD8_Seurat_subset_cxcr3, features = c("ZNF683", "CXCR3", "HLA-DRA", "HLA-DRB5", "TCF7", "CCL5", "GZMA", "GZMH", "GNLY"), cols = c(singlecellcolors[1:2], "gray70"))
ggsave("analysis/figs/sc_nonnaiveCD8_cxcr3_MLRmatch_violin.eps", width = 3, height = 4.5)



#comparing pre pembro to pembro 1 in bulk TCR
Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$barcode

emerged_clones_nonnaive <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% pembro_bloodexpanded_emerged_clones$CDR3aa) 

emergedclonesmatchumap <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = emerged_clones_nonnaive$barcode), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Emerged", "Other"),
                     values = c("#073f7d", "gray90")) + xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Emerged (Expanded)\n Clones") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15))

existing_clones_nonnaive <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% pembro_bloodexpanded_existing_clones$CDR3aa) 

existingclonesmatchumap <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = existing_clones_nonnaive$barcode), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Existing", "Other"),
                     values = c("#073f7d", "gray90")) + xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Existing (Expanded)\n Clones") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15))

stable_clones_nonnaive <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% pembro_bloodstable_all_clones$CDR3aa) 

stableclonesmatchumap <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = stable_clones_nonnaive$barcode), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Stable", "Other"),
                     values = c("#073f7d", "gray90")) + xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Stable\n Clones") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15))


contract_clones_nonnaive <- nonnaiveCD8_Seurat_subset@meta.data %>%
  filter(beta %in% pembro_bloodcontract_all_clones$CDR3aa) 

contractclonesmatchumap <- DimPlot(nonnaiveCD8_Seurat_subset, cells.highlight = WhichCells(nonnaiveCD8_Seurat_subset, idents = contract_clones_nonnaive$barcode), sizes.highlight = 0.5) +
  scale_color_manual(limits = c("Group_1", "Unselected"), 
                     labels = c("Contracted", "Other"),
                     values = c("#073f7d", "gray90")) + xlab("UMAP 1") + ylab("UMAP 2") +
  ggtitle("Contracted\n Clones") +
  theme(panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 15))

nonnaive_expandcontract_UMAPs <- align_plots(emergedclonesmatchumap, existingclonesmatchumap, stableclonesmatchumap, contractclonesmatchumap, align="hv", axis="tblr")
p1x <- ggdraw(nonnaive_expandcontract_UMAPs[[1]])
p2x <- ggdraw(nonnaive_expandcontract_UMAPs[[2]])
p3x <- ggdraw(nonnaive_expandcontract_UMAPs[[3]])
p4x <- ggdraw(nonnaive_expandcontract_UMAPs[[4]])

save_plot("analysis/figs/sc_nonnaive_UMAP_expandedEmerged_aligned.eps", p1x, base_width = 3.5, base_height = 3.7)
save_plot("analysis/figs/sc_nonnaive_UMAP_expandedExisting_aligned.eps", p2x, base_width = 3.5, base_height = 3.7)
save_plot("analysis/figs/sc_nonnaive_UMAP_stable_aligned.eps", p3x, base_width = 3.5, base_height = 3.7)
save_plot("analysis/figs/sc_nonnaive_UMAP_contracted_aligned.eps", p4x, base_width = 3.5, base_height = 3.7)


nonnaiveCD8_Seurat_subset@meta.data %>%
  mutate(class = ifelse(barcode %in% emerged_clones_nonnaive$barcode, "Emerged", 
                        ifelse(barcode %in% existing_clones_nonnaive$barcode, "Existing",
                               ifelse(barcode %in% stable_clones_nonnaive$barcode, "Stable",
                                      ifelse(barcode %in% contract_clones_nonnaive$barcode, "Contracted", "Other"))))) %>%
  filter(class != "Other") %>%
  group_by(clusters_names, class) %>%
  dplyr::count() %>%
  group_by(class) %>%
  mutate(Percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=class,y=Percent, fill=clusters_names)) +
  geom_col() +
  scale_fill_manual(values = brewer.pal(10, "Paired"), name = "Cell Type", drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  scale_x_discrete(limits = rev(c("Emerged", "Existing", "Stable", "Contracted"))) +
  ggtitle(NULL) +
  labs(fill = "Sample") +
  xlab(NULL) + ylab(NULL) +
  theme_classic() + 
  theme(axis.text.x = element_text(size = 9),
        legend.text = element_text(size=9),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.35, 'cm'),
        legend.position = "none") + coord_flip()
ggsave("analysis/figs/sc_nonnaiveCD8_expandContract.eps", width = 3.5, height = 1.5)


VlnPlot(nonnaiveCD8_Seurat_subset, features = c("CXCR3", "ZNF683", "CD38", "PDCD1"),
        group.by = "clusters_names", split.by = "timepoint", 
        cols = timepoint_cols, split.plot = T, ncol = 1) & 
  xlab(NULL)

cd8signatures <- read.csv("cd8_signatures.csv")
nonnaiveCD8_Seurat_subset <- AddModuleScore(nonnaiveCD8_Seurat_subset, list(cd8signatures$Activation), name = "Activation_sig", ctrl = 1000)

nonnaiveCD8_Seurat_subset@meta.data %>%
  ggplot(aes(x = clusters_names, y = Activation_sig1, fill = clusters_names)) +
  geom_boxplot(outlier.size = 0.5) + theme_classic() + theme(legend.position = "none", axis.text=element_text(size=10)) + 
  xlab((NULL)) + ylab("Activation Signature") + 
  scale_fill_manual(values = cd8cols) +
  scale_y_continuous(breaks=c(0,0.4), labels=c("Min","Max"), limits = c(0,0.4)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 9))
ggsave("analysis/figs/sc_nonnaiveCD8_activationSig.eps", width = 3.5, height = 2)

#plotting same dim UMAPs
nonnaive_UMAPs <- align_plots(nonnaiveclusterUMAP, nonnaiveclusterUMAPPlabeled, nonnaivetimepointUMAP, nonnaivepembro1UMAP, nonnaivepembro2UMAP, nonnaiveclonalityUMAP, align="hv", axis="tblr")
p1x <- ggdraw(nonnaive_UMAPs[[1]])
p2x <- ggdraw(nonnaive_UMAPs[[2]])
p3x <- ggdraw(nonnaive_UMAPs[[3]])
p4x <- ggdraw(nonnaive_UMAPs[[4]])
p5x <- ggdraw(nonnaive_UMAPs[[5]])
p6x <- ggdraw(nonnaive_UMAPs[[6]])

save_plot("analysis/figs/sc_nonnaive_UMAP_byCluster_aligned.eps", p1x, base_width = 5.5, base_height = 3.5)
save_plot("analysis/figs/sc_nonnaive_UMAP_byCluster_labeled_aligned.eps", p2x, base_width = 5.5, base_height = 3.5)
save_plot("analysis/figs/sc_nonnaive_UMAP_byTimepoint_aligned.eps", p3x, base_width = 5.5, base_height = 3.5)
save_plot("analysis/figs/sc_nonnaive_UMAP_pembro1_aligned.eps", p4x, base_width = 5.5, base_height = 3.5)
save_plot("analysis/figs/sc_nonnaive_UMAP_pembro2_aligned.eps", p5x, base_width = 5.5, base_height = 3.5)
save_plot("analysis/figs/sc_nonnaive_UMAP_clonality_MLR_aligned.eps", p6x, base_width = 5.5, base_height = 3.5)

Idents(nonnaiveCD8_Seurat_subset) <- nonnaiveCD8_Seurat_subset$clusters_names
CD8prolif <- subset(nonnaiveCD8_Seurat_subset, subset = clusters_names == "Cycling")
CD8_exceptprolif <- subset(nonnaiveCD8_Seurat_subset, subset = clusters_names == "Cycling", invert = T)

symphT_reference <- symphony::buildReference(
  CD8_exceptprolif@assays$RNA@data,
  CD8_exceptprolif@meta.data,
  vars = c("timepoint"),         # variables to integrate over
  K = 100,                   # number of Harmony clusters
  verbose = TRUE,            # verbose output
  do_umap = FALSE,            # can set to FALSE if want to run umap separately later
  do_normalize = FALSE,      # set to TRUE if input counts are not normalized yet
  vargenes_method = 'vst',   # method for variable gene selection ('vst' or 'mvp')
  #vargenes_groups = 'clusters_names', # metadata column specifying groups for variable gene selection 
  topn = 2000,               # number of variable genes to choose per group
  d = 20                     # number of PCs
)

queryT <- mapQuery(CD8prolif$RNA@data,             # query gene expression (genes x cells)
                   CD8prolif@meta.data,        # query metadata (cells x attributes)
                   symphT_reference,             # Symphony reference object
                   vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                   do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                   do_umap = FALSE)        # project query cells into reference UMAP

queryT <- knnPredict(queryT, symphT_reference, train_labels = symphT_reference$meta_data$clusters_names,
                     k = 10, save_as = 'cell_type_10_nn', confidence = TRUE)

CD8_cluster_breakdown_exceptprolif <- as.data.frame(table(CD8_exceptprolif$clusters_names))
ggplot(CD8_cluster_breakdown_exceptprolif, aes(x = 2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  guides(fill = guide_legend(reverse = TRUE)) +
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = cd8cols) +
  theme_void() +
  xlim(0.5, 2.5) +
  theme(legend.key.size = unit(0.35, 'cm'),
        legend.text = element_text(size = 9))

CD8_reClusterProlif_cluster_breakdown <- as.data.frame(table(queryT$meta_data$cell_type_10_nn))
ggplot(CD8_reClusterProlif_cluster_breakdown, aes(x = 2, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", color = "white") +
  coord_polar(theta = "y", start = 0) +
  guides(fill = guide_legend(reverse = TRUE)) +
  guides(fill=guide_legend(title=NULL)) +
  scale_fill_manual(values = cd8cols) +
  theme_void() +
  xlim(0.5, 2.5) +
  theme(legend.key.size = unit(0.35, 'cm'),
        legend.text = element_text(size = 9))  
ggsave("analysis/figs/sc_nonnaiveCD8_prolif_reCluster.eps", width = 3, height = 2)


#trajectory on selected clusters
nonnaive_subset_slingshotremove <- subset(nonnaiveCD8_Seurat_subset, idents = c("Cycling", "Mito-hi", "MAIT"), invert = TRUE)

nonnaive_TI_subset_unbiasedstart <- slingshot(Embeddings(nonnaive_subset_slingshotremove, "umap"), clusterLabels = nonnaive_subset_slingshotremove$clusters_names, 
                                              stretch = 0, omega = T)

setEPS()
postscript("analysis/figs/nonnaive_sc_slingshot_subset_unbiased.eps", width = 5, height = 4)
par(mar=c(5,4.5,1,1))
plot(Embeddings(nonnaiveCD8_Seurat_subset, "umap"), col = cd8cols[Idents(nonnaiveCD8_Seurat_subset)], pch = 16, cex = 0.5, axes = F)
box(bty="l")
axis(1)
axis(2, las = 2)
lines(nonnaive_TI_subset_unbiasedstart, lwd = 2, col = 'black', type = 'lineages')
dev.off()

nc <- 3
pt <- slingPseudotime(nonnaive_TI_subset_unbiasedstart)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(nonnaive_TI_subset_unbiasedstart), col = colors, pch = 16, cex = 0.5, main = i)
  lines(nonnaive_TI_subset_unbiasedstart, lwd = 2, col = 'black', type = 'lineages')
}

# Final system output -----------------------

writeLines(capture.output(sessionInfo()), "analysis/data/sessionInfo.txt")

save.image("analysis/data/kidneyICI_projectAnalysis.RData")
