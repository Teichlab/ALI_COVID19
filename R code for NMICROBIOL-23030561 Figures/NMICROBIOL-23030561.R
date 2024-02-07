# Code for NMICROBIOL-23030561 Main Figures written by Claire M Smith c.m.smith@ucl.ac.uk

#load main object
obj <- readRDS("seurat_ALI_raw_counts.rds")



#=====================================================================================
#
#  Code chunk Figure 1d
#
#=====================================================================================

library(ggrepel)
library(ggpubr)

df$age_group <- factor(df$age_group, levels = c('Paed', 'Adult', 'Elderly'))
df$entry_factor <- factor(df$entry_factor, levels = c('TMPRSS2', 'ACE2', 'shACE2'))

my_colors2 <- c("#85c5c9", "#b09a6f", "#70b06f")            # Create vector of colors
names(my_colors2) <- levels(factor(c(levels(df$age_group), levels(df$age_group)))) # Extract all levels of both data


p <- df2 %>%
  filter(treatment == "mock") %>%
  ggplot(aes(x=age_group, y=ratio.GAPDH , fill = age_group)) +
  geom_point( shape = 21, size = 3) +
  geom_boxplot(alpha=0.5)+
  labs(title = 'Entry Factor Protein')+
  xlab("Age group")+
  theme(legend.position = "right") + theme(legend.direction = "vertical") +
  scale_fill_manual(name = "Group", values = my_colors2)  +
  facet_wrap(~entry_factor, scales = "free_y")+ theme_bw(base_size = 16)+ theme(legend.position = "none")
p

my_comparisons <- list( c("Paed", "Adult"), c("Adult", "Elderly"), c("Paed", "Elderly"))
p +stat_compare_means(comparisons = my_comparisons, method = 't.test', hide.ns = FALSE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+ # Add pairwise comparisons p-value
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))

#=====================================================================================
#
#  Code chunk Figure 1f
#
#=====================================================================================

obj <- SetIdent(obj, value = obj@meta.data$time)
levels(obj)

SARS_expression = GetAssayData(object = obj, assay = "RNA", slot = "data")["VIRAL-SARS-CoV2",]
pos_ids = names(which(SARS_expression>0))
pos_cells = subset(obj,cells=pos_ids)
pos_cells@meta.data$age_group <- factor(pos_cells@meta.data$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))

vln_df = data.frame(Expression_Level =pos_cells[["RNA"]]@data[ "VIRAL-SARS-CoV2", ], donor = pos_cells$donor_id, age_group = pos_cells$age_group, time = pos_cells$time)

df_mean <- vln_df %>%
  group_by(age_group, time, donor) %>%
  summarize(ave2 = n(), ave = mean(Expression_Level)) %>%
  ungroup()

df_mean$age_group <- factor(df_mean$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
vln_df$age_group <- factor(vln_df$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))

my_colors <- c("#85c5c9", "#b09a6f", "#70b06f")
names(my_colors) <- levels(factor(c(levels(vln_df$age_group), levels(vln_df$age_group)))) # Extract all levels of both data

vln <- ggplot(df_mean, aes(x = age_group, y = ave, fill = age_group)) +
  geom_jitter(data = vln_df,
              mapping = aes(x = age_group, y = Expression_Level),
              color="grey", alpha = 0.1) +
  geom_boxplot()+
  scale_fill_manual(name = "Group", values = my_colors)  +
  geom_jitter(alpha = 1, color = "red") +
  ylab("Expression level/donor")+
  xlab("Age group")+
  theme_bw(base_size = 20) + theme(legend.position = "none") + theme(legend.direction = "vertical")+
  facet_wrap(~time, ncol = 2)
vln

my_comparisons <- list( c("Paediatric", "Adult"), c("Adult", "Elderly"), c("Paediatric", "Elderly"))
vln+   stat_compare_means(comparisons = my_comparisons)+ # Add pairwise comparisons p-value
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.15)))

#=====================================================================================
#
#  Code chunk Figure 1g
#
#=====================================================================================

##load source data as df2
df2$age_group <- factor(df2$age_group, levels = c('Elderly', 'Adult', 'Paediatric'))
df2$timepoint <- factor(df2$timepoint, levels = c('4h', '24h', '72h'))
df2$fraction_1000cells <- as.numeric(df2$fraction_1000cells)
df2$read_percell <- as.numeric(df2$read_percell)

#order
df2$annotation_v5 <- factor(df2$annotation_v5 , levels = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1', 'Basaloid-like 2', 'Cycling basal', 'Hillock', 'Ionocyte', 'Goblet 1', 'Goblet 2',
                                                           'Goblet 2 BPIFA1+', 'Goblet 2 PLAU+', 'Goblet 2 inflammatory', 'Secretory', 'Secretory 2', 'Secretory 3', 'Secretory 4', 'Ciliated 1', 'Ciliated 2', 'Deutorosomal',  'Squamous', 'Transit epi', 'Transit epi 2'))
###as dot plot
fp <-  df2 %>%
  filter(timepoint != "4h")%>%
  filter(Infected_count >0)%>%
  filter(treatment == "SARS")%>%
  ggplot(aes(x = annotation_v5, y=age_group, fill = read_percell, size = fraction_1000cells)) +
  geom_point(shape = 21)+
  scale_fill_gradientn(colours = c("white", "dark red", "yellow"), limits = c(0,2), oob = scales::squish, name = 'v.reads/cell')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  labs(title = 'Viral reads/cell type') +
  ylab('')+
  facet_wrap(~timepoint, ncol = 1, strip.position = "top")
fp



#=====================================================================================
#
#  Code chunk Figure 1k
#
#=====================================================================================

##load source data as df4

df4$Age_group <- factor(df4$Age_group, levels = c('Elderly', 'Adult', 'Paed'))
df4$Site <- factor(df4$Site, levels = c('Intracellular', 'Extracellular'))

df_mean <- df4 %>%
  group_by(Age_group, Protein, Site) %>%
  summarize(ave = mean(value), FC= mean(Fold.Change)) %>%
  ungroup()

my_colors2 <- c("#85c5c9", "#b09a6f", "#70b06f")            # Create vector of colors
names(my_colors2) <- levels(factor(c(levels(df4$Age_group), levels(df4$Age_group)))) # Extract all levels of both data

dotplot <- df_mean %>%
  filter(ave >0)%>%
  filter(Site == "Intracellular" | Site == "Extracellular") %>%
  ggplot(aes(x=Protein, y = Age_group, fill = Age_group, size = ave)) +
  geom_point(shape = 21) +
  cowplot::theme_cowplot() +
  labs(title="", y = "")+
  scale_fill_manual(name = "", values = my_colors2)  +
  theme(axis.ticks = element_blank()) +
  labs(title = 'SARS-CoV-2 Protein Expression')+
  scale_size_area(max_size = 10, name = 'Abundance')+
  facet_wrap(~Site, ncol=3)

dotplot


#=====================================================================================
#
#  Code chunk Figure 3g
#
#=====================================================================================
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggpubr)

##load source data as df2
df2$age_group <- factor(df2$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
df2$time <- factor(df2$time , levels = c('72h', '24h', '4h'))
df2$treatment <- factor(df2$treatment , levels = c('mock', 'SARS'))
df2$log.value <- log10(df2$value)

df_mean <- df2 %>%
  group_by(age_group, variable, treatment, time) %>%
  summarize(average = mean(value), n=n(), sd=sd(value))

bar2 <-  df2 %>%
  ggplot(aes(x = treatment, y = value,  fill= time, group = treatment)) +
  geom_point(color = "black", shape =21, size = 2)+
  geom_col(data = df_mean, group = "treatment",
           mapping = aes (x=treatment, y = average),alpha = 0.5,  color = "black",
           position=position_stack())+
  xlab('Condition') +
  ylab('pg/ml') +
  labs(title = "IFN Protein")+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  theme_bw(base_size = 24)+
  theme(legend.position = "right")+
  facet_grid(variable~age_group , scales = "free_y")+
  stat_compare_means(comparisons=my_comparisons, method = "t.test",label.y = Inf, vjust = 3, hide.ns = FALSE,   paired = TRUE, symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns")))+ # Add pairwise comparisons p-value
  scale_y_continuous(expand = expansion(mult = c(0.02, .5)))
bar2

bar2+ scale_fill_grey(start = 0, end = 1)




#=====================================================================================
#
#  Code chunk Figure 4b
#
#=====================================================================================

#df2 <- source data for 4b
df_mean <- df2 %>%
  group_by(Age_Group) %>%
  summarize(sgRPTL = (sgRPTL), n= n(), mean = median(sgRPTL)) %>%
  ungroup()

pp <- df2 %>%
  ggplot(mapping = aes(x = Age_Group, y = sgRPTL, fill = Age_Group)) +
  geom_boxplot(alpha = 0.6) +
  geom_point(alpha = 0.2, size =2, stat = "identity", position=position_jitter(width = .2))+
  geom_point(data = df_mean,
             mapping = aes(x = Age_Group, y = mean),
             color="red") +
  geom_line(data = df_mean,
            mapping = aes(x = Age_Group, y = mean, group =1))+
  scale_fill_manual(name = "Group", values = my_colors)  +
  theme_bw(base_size = 18) + theme(legend.position = "none") + theme(legend.direction = "vertical")+
  scale_y_continuous(trans='log2')+
  labs(title = 'Total sgRNA')+
  ylab('sgRPTL') +
  xlab('Age Group')

pp

#=====================================================================================
#
#  Code chunk Figure 4c
#
#=====================================================================================

#df2 <- source data
my_colors <- c("#85c5c9", "#b09a6f", "#70b06f")            # Create vector of colors
names(my_colors) <- levels(factor(c(levels(df2$Age_Group), levels(df2$Age_Group)))) # Extract all levels of both data

df2$orf <- factor(df2$orf, levels = c("S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF8", "N"))
df2$Age_Group <- as.factor(df2$Age_Group)
df2$Age_Group <- factor(df2$Age_Group, levels = c("Paed", "Adult", "Elderly"))

df_mean2 <- df2 %>%
  group_by(sample, Donor, Age_Group, orf) %>%
  summarize(Count_sgRPTL = sum(sgRPTL)) %>%
  ungroup()

gg <- df2 %>%
  filter(Age_Group != "Adult") %>%
  ggbarplot(x = "orf" ,y= "sgRPTL" ,fill= "Age_Group", add = "mean_se",position = position_dodge(0.8))+
  geom_point(aes(x = orf ,y= sgRPTL ,fill= Age_Group), position = position_dodge(0.9), shape = 21)+
  scale_fill_manual(name = "Group", values = my_colors)  +
  ylab('sgRPTL') +
  xlab('Gene') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  coord_flip()+
  theme_bw(base_size = 18) +
  theme(legend.position = "none")
gg


#=====================================================================================
#
#  Code chunk Figure 5a
#
#====================================================================================

#df <- read source data
df$age_group <- factor(df$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
df$timepoint <- factor(df$timepoint, levels = c('4h', '24h', '72h'))
df$log2 <- as.numeric(df$log2)

##blue, orange, gree
my_colors <- c("#85c5c9", "#b09a6f", "#70b06f")            # Create vector of colors
names(my_colors) <- levels(factor(c(levels(df$age_group), levels(df$age_group)))) # Extract all levels of both data

#order
df$annotation_v5 <- factor(df$annotation_v5 , levels = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1', 'Basaloid-like 2', 'Cycling basal', 'Hillock'))

df_sum <- df %>%
  group_by(age_group, treatment, annotation_v5) %>%
  summarize(sum = sum(per_1000), FC = mean(log2)) %>%
  ungroup()


fp <-  df_sum %>%
  ggplot(aes(x = annotation_v5, y = sum, fill = FC, color = treatment)) +
  geom_col(position = position_dodge(0.9))+
  scale_color_manual(values = c("black", "dark red"))+
  scale_fill_gradient2(low = 'white', high = 'yellow', mid = 'red',  midpoint = 0, limit = c(0, 2), oob = scales::squish, name = 'FC from mock')+
  ylab('Frequency/1000 cells') +
  labs(title = 'KRT5hi cells all TPs')+
  theme_bw(base_size = 14) + theme(legend.position = "right") + theme(legend.direction = "vertical") +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~age_group, ncol = 1, strip.position = "right")

fp

#=====================================================================================
#
#  Code chunk Figure 5c
#
#====================================================================================

#df <- read source data 5c
library(ggrepel)

p <- df %>%
  ggplot(aes(x=log2FoldChange, y=-log10p_val, col=diffexpressed, label=Protein)) +
  geom_point(color = ifelse(df$Protein == "", "black", "red")) +
  theme_minimal(base_size = 18) + theme(legend.position = "right") + theme(legend.direction = "vertical") +
  geom_label_repel(color = ifelse(df$Protein == "", "grey", "red"), max.overlaps = 300)+
  scale_color_manual(values=c("blue", "black", "red")) +
  labs(title = 'Elderly Secretome DE:SARS v Mock')+
  xlim(-5,5)+
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=(4), col="red")

p
p + geom_text_repel(aes(label=delabel), max.overlaps = 10) +
  theme(strip.background = element_rect(colour = "black", fill = "white"), strip.text.x = element_text(colour = "black", face = "bold"))
