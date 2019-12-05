#Libraries of required packages-------------------------------------------------------
library(DESeq2)
library(ggplot2)
library(vsn)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
library(calibrate)
library(msmsTests)
library(xlsx)
library(magrittr)

#Build DESeq Matrix---------------------------------------------------
dds_final <- DESeqDataSetFromMatrix(countData = Count_adult2, 
                                          colData = data_adult,
                                          design = ~ Gender + Genotype + Gender:Genotype) 
dds_final

dds_final$Gender = relevel(dds_final$Gender, "Male")
dds_final$Gender
dds_final$Genotype = relevel(dds_final$Genotype, "WT")
dds_final$Genotype

#DESeq---------------------------------------------
design(dds_final) <- ~ Gender + Genotype + Gender:Genotype
dds_final <- DESeq(dds_final)
resultsNames(dds_final)

###The effect of gender on WT mice. (WT: Male vs. Female)----------------
res = results(dds_final, contrast=c("Gender","Female","Male"), alpha = 0.05)
ix = which.min(res$padj) # most significant
res <- res[order(res$padj),] # sort
summary(res)

#plot of the most significant gene in above comparison. ie. Xist
#barplot(assay(dds_final)[ix,],las=2, main=rownames(dds_final)[ ix  ]  )

###The effect of gender in HET mice. (HET: Male vs. Female)
res_HET <- results(dds_final, list( c("Gender_Female_vs_Male","GenderFemale.GenotypeHET") ), alpha = 0.05)
ix_HET = which.min(res_HET$padj) # most significant
res_HET <- res_HET[order(res_HET$padj),] # sort
summary(res_HET)

###the difference between HET and WT w/o gender effect. (Male: WT vs. HET)
res_male = results(dds_final, contrast=c("Genotype","HET","WT"), alpha = 0.05)
ix_male = which.min(res_male$padj) # most significant
res_male <- res_male[order(res_male$padj),] # sort
summary(res_male)

barplot(assay(dds_final)[ix_male,],las=2, main=rownames(dds_final)[ix_male]  )

###the difference between HET and WT w/ gender effect. (Female: WT vs. HET)
res_female = results(dds_final, list( c("Genotype_HET_vs_WT","GenderFemale.GenotypeHET") ), alpha = 0.05)
ix_female = which.min(res_female$padj) # most significant
res_female <- res_female[order(res_female$padj),] # sort
summary(res_female)

###The different response in genotypes to gender effect. (interaction term)
# is the effect of gender different across WT and HET?
res_inter = results(dds_final, name="GenderFemale.GenotypeHET", alpha = 0.05)
ix_inter = which.min(res_inter$padj) # most significant
res_inter <- res_inter[order(res_inter$padj),] # sort
summary(res_inter)

barplot(assay(dds_final)[ix_inter,],las=2, main=rownames(dds_final)[ix_inter]  )


#output files-------------------------
write.xlsx(as.data.frame(res),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/res_WTmaleVfemale.xlsx")
write.xlsx(as.data.frame(res_HET),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/res_HETmaleVfemale.xlsx")
write.xlsx(as.data.frame(res_male),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/res_MaleWTvHET.xlsx")
write.xlsx(as.data.frame(res_female),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/res_FemaleWTvHET.xlsx")
write.xlsx(as.data.frame(res_inter),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/res_interaction.xlsx")

#get raw data of each gene list
result_res_new = na.omit(res)
result_res_new <- as.data.frame(result_res_new, header = T)
result_res_new_05 <- subset(result_res_new, padj < 0.05)
result_res_names <- rownames(result_res_new_05)
result_res_raw <- Count_adult2[result_res_names, ]
write.xlsx(as.data.frame(result_res_raw),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/raw_WTmaleVfemale.xlsx")

result_resHET_new = na.omit(res_HET)
result_resHET_new <- as.data.frame(result_resHET_new, header = T)
result_resHET_new_05 <- subset(result_resHET_new, padj < 0.05)
result_resHET_names <- rownames(result_resHET_new_05)
result_resHET_raw <- Count_adult2[result_resHET_names, ]
write.xlsx(as.data.frame(result_resHET_raw),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/raw_HETmaleVfemale.xlsx")

result_resmale_new = na.omit(res_male)
result_resmale_new <- as.data.frame(result_resmale_new, header = T)
result_resmale_new_05 <- subset(result_resmale_new, padj < 0.05)
result_resmale_names <- rownames(result_resmale_new_05)
result_resmale_raw <- Count_adult2[result_resmale_names, ]
write.xlsx(as.data.frame(result_resmale_raw),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/raw_MaleWTvHET.xlsx")

raw_FemaleWTvHET <- Count_adult2[res_female$padj < 0.05, ]
raw_FemaleWTvHET <- na.omit(raw_FemaleWTvHET)
write.xlsx(as.data.frame(raw_FemaleWTvHET),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/raw_FemaleWTvHET.xlsx")

raw_interaction <- Count_adult2[res_inter$padj < 0.05, ]
raw_interaction <- na.omit(raw_interaction)
write.xlsx(as.data.frame(raw_interaction),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/raw_interaction.xlsx")

list_interest <- c("Stx1b", "Vamp2", "Bet1", "Itgam", "Ppp1cb", "Ppp2cb", "Akt1", "Gnai1")
raw_interest <- Count_adult2[list_interest, ]
write.xlsx(as.data.frame(raw_interest),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/interest_raw.xlsx")

#get the gene counts for qPCR from raw data.----------------------------------
qPCR_list <- c("Slc3a1","Fam120a","Cadm2","Ankrd10","Hipk3","Rn45s","1700058G18Rik","Lnp","Map4k5","Bmi1","Tspan7","Spock3", "Ptprn2")
qPCR_refList <- c("Gapdh","Actb","Akt1","Xist")
qPCR_gene_counts <- Count_adult2[qPCR_list, ]
qPCR_ref_counts <- Count_adult2[qPCR_refList, ]
write.xlsx(as.data.frame(qPCR_gene_counts),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/inter_gene_raw_2.xlsx")
write.xlsx(as.data.frame(qPCR_ref_counts),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/20190530/ref_gene_raw.xlsx")

akt_isoforms <- c("Akt1", "Akt2", "Akt3")
Akt_raw <- Count_adult2[akt_isoforms, ]
#Akt1_raw <- Count_adult2["Akt1", ]
write.xlsx(as.data.frame(Akt_raw),
           file="C:/Users/Tzu-Ping Lee/Desktop/RNA_analysis/DESeq2 analysis/Akt_raw.xlsx")

## Data visualization------------------------------------- 
################################
# PCA and heatmap, clusteing analysis with pcaExplorer
library(pcaExplorer)
pcaExplorer(dds = dds_final)

##Volcano Plots##-----------------------------------
library(ggrepel)

#axis titles
x_title <- expression(log["2"]("Fold Change"))
y_title <- expression(-log["10"](paste(italic("P"), "-value")))

#WT only
volcano.input.WT <- data.frame(gene = row.names(res),
                               pvalue = res$padj, 
                               lfc = res$log2FoldChange)
volcano.input.WT <- na.omit(volcano.input.WT)

#coord_cartesian(xlim=c(-3,3), ylim=c(0, 8)) [this is for adding axis limits]
plot.WT <- ggplot(volcano.input.WT, aes(x = lfc, y = -log10(pvalue))) +
  geom_point() + 
  labs(title = "WT♂ vs. WT♀", y = y_title, x = x_title) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(xintercept = c(0), linetype = "dotted") + 
  coord_cartesian(xlim=c(-13,10), ylim=c(0, 50)) + 
  theme_classic(base_size = 25)
plot.WT 

#subset p<0.05 genes for labelling, may use this to require genes with desired criteria of pvalue and fold change
volcano.WT.sig <- subset(volcano.input.WT, pvalue < 0.05) #%>% subset(abs(lfc) < 12)

plot.WT + 
  geom_text_repel(aes(label=gene), size=6, data=volcano.WT.sig) + 
  geom_point(data = volcano.WT.sig, color = "red") +
  theme(plot.title = element_text(hjust = 0.5))

write.csv(volcano.WT.sig, file="C:/Users/Tzu-Ping Lee/Desktop/論文/Graphs/volcano_WTp0.05.csv")

#HET only
volcano.input.HET <- data.frame(gene = row.names(res_HET),
                               pvalue = res_HET$padj, 
                               lfc = res_HET$log2FoldChange)
volcano.input.HET <- na.omit(volcano.input.HET)

HETonly_title <- expression(paste(italic("Akt1")^"+/-", " ♂", " vs. ", italic("Akt1")^"+/-", " ♀"))

volcano.HET.sig <- subset(volcano.input.HET, pvalue < 0.05) %>% subset(abs(lfc) > log2(2))
rownames(volcano.HET.sig) <- volcano.HET.sig$gene
volcano.HET.sig$gene <- NULL
HET_gene <- c("Hspa1l", "Foxp3")
HET_color_gene <- volcano.HET.sig[HET_gene,]

plot.HET <- ggplot(volcano.input.HET, aes(x = lfc, y = -log10(pvalue))) +
  geom_point() + 
  labs(title = HETonly_title, y = y_title, x = x_title) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(1, -1), linetype = "dotted") +
  coord_cartesian(xlim=c(-10,10), ylim=c(0, 50)) + 
  theme_classic(base_size = 25) + 
  geom_point(data = volcano.HET.sig, color = "red") + 
  geom_point(data = HET_color_gene, color = "blue") + 
  geom_label_repel(aes(label=row.names(HET_color_gene)), size=8, data=HET_color_gene, 
                   box.padding = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(plot.title = element_text(hjust = 0.5))
plot.HET

write.csv(volcano.HET.sig, file="volcano_HETlogFC1.csv")

#Male only
volcano.input.male <- data.frame(gene = row.names(res_male),
                                pvalue = res_male$padj, 
                                lfc = res_male$log2FoldChange)
volcano.input.male <- na.omit(volcano.input.male)

Maleonly_title <- expression(paste("WT", " ♂", " vs. ", italic("Akt1")^"+/-", " ♂"))

volcano.male.sig <- subset(volcano.input.male, pvalue < 0.05) %>% subset(abs(lfc) > log2(1.3))
rownames(volcano.male.sig) <- volcano.male.sig$gene
volcano.male.sig$gene <- NULL
male_gene <- c("Flt1")
male_color_gene <- volcano.male.sig[male_gene,]

plot.male <- ggplot(volcano.input.male, aes(x = lfc, y = -log10(pvalue))) +
  geom_point() + 
  labs(title = Maleonly_title, y = y_title, x = x_title) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(log2(1.3), -log2(1.3)), linetype = "dotted") +
  coord_cartesian(xlim=c(-4,4), ylim=c(0, 5)) + 
  theme_classic(base_size = 25) + 
  #geom_text_repel(aes(label=gene), size=3, data=volcano.male.sig) + 
  geom_label_repel(aes(label=row.names(male_color_gene)), size=8, data=male_color_gene, 
                   box.padding = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  geom_point(data = volcano.male.sig, color = "red") +
  theme(plot.title = element_text(hjust = 0.5))
plot.male

write.csv(volcano.male.sig, file="volcano_maleFC1.3.csv")

#Female Only
volcano.input.female <- data.frame(gene = row.names(res_female),
                                 pvalue = res_female$padj, 
                                 lfc = res_female$log2FoldChange)
volcano.input.female <- na.omit(volcano.input.female)

Femaleonly_title <- expression(paste("WT", " ♀", " vs. ", italic("Akt1")^"+/-", " ♀"))

volcano.female.sig <- subset(volcano.input.female, pvalue < 0.05) %>% subset(abs(lfc) > log2(1.3))
rownames(volcano.female.sig) <- volcano.female.sig$gene
volcano.female.sig$gene <- NULL
female_gene <- c("Tspan7")
female_color_gene <- volcano.female.sig[female_gene,]

plot.female <- ggplot(volcano.input.female, aes(x = lfc, y = -log10(pvalue))) +
  geom_point() + 
  labs(title = Femaleonly_title, y = y_title, x = x_title) + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(log2(1.3), -log2(1.3)), linetype = "dotted") +
  coord_cartesian(xlim=c(-4,4), ylim=c(0, 15)) + 
  theme_classic(base_size = 25) + 
  geom_point(data = volcano.female.sig, color = "red") + 
  geom_point(data = female_color_gene, color = "blue") + 
  geom_label_repel(aes(label=row.names(female_color_gene)), size=8, data=female_color_gene, 
                   box.padding = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50') +
  theme(plot.title = element_text(hjust = 0.5))
plot.female

write.csv(volcano.female.sig, file="volcano_femaleFC1.3.csv")

#combine the graphs into single image.
library(ggpubr)
library(cowplot)

plots_comb <- ggarrange(plot.HET, plot.male, plot.female, ncol = 3, nrow = 1)
ggsave("vol_plots_comb", plot = plots_comb, device = "pdf" ,path = "/20190716_VolPlots",
       scale = 1, width = 36, height = 12.23, units = c("cm"),
       dpi = 900, limitsize = TRUE) #ggsave error not resolved.... 20190716

#interactions res_inter
volcano.input.inter <- data.frame(gene = row.names(res_inter),
                                   pvalue = res_inter$padj, 
                                   lfc = res_inter$log2FoldChange)
volcano.input.inter <- na.omit(volcano.input.inter)

plot.inter <- ggplot(volcano.input.inter, aes(x = lfc, y = -log10(pvalue))) +
  geom_point() + 
  labs(title = "Interaction", y = "-log10(P-value)", x = "Log2(Fold Change)") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") + geom_vline(xintercept = c(1, -1), linetype = "dotted") +
  theme_classic()
plot.inter

volcano.inter.sig <- subset(volcano.input.inter, pvalue < 0.05) #%>% subset(abs(lfc) > 1)

plot.inter + 
  geom_text_repel(aes(label=gene), size=2, data=volcano.inter.sig) + 
  geom_point(data = volcano.inter.sig, color = "red")
