# Code for NMICROBIOL-23030561 Extended Figures written by Claire M Smith c.m.smith@ucl.ac.uk


#=====================================================================================
#
#  Code chunk Extended Figure 1c
#
#=====================================================================================

#df <- load source data

my_levels <- c('Paediatric', 'Adult', 'Elderly')
df$age_group <- factor(df$age_group , levels = my_levels)

dp2 <-  df %>%
  ggplot(aes(x = treatment, y = perc, fill=phase)) +
  geom_col(position = position_stack())+
  xlab('Treatment') +
  ylab('% Total') +
  labs(title = 'Cell Cycle Phase')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~age_group, ncol = 3, strip.position = "top")

dp2

#=====================================================================================
#
#  Code chunk Extended Figure 1g
#
#=====================================================================================

#df <- load source data from Ext Fig 1g
df$age_group <- factor(df$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
my_colors <- c("#85c5c9", "#b09a6f", "#70b06f")            # Create vector of colors
names(my_colors) <- levels(factor(c(levels(df$age_group), levels(df$age_group)))) # Extract all levels of both data

#order
df$annotation_v5 <- factor(df$annotation_v5 , levels = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1', 'Basaloid-like 2', 'Cycling basal', 'Hillock', 'Secretory', 'Secretory 2', 'Secretory 3', 'Secretory 4', 'Transit epi', 'Transit epi 2', 'Goblet 1', 'Goblet 2',
                                                         'Goblet 2 BPIFA1+', 'Goblet 2 PLAU+', 'Goblet 2 inflammatory', 'Ciliated 1', 'Ciliated 2', 'Deutorosomal', 'Ionocyte', 'Squamous'))
plot1 <-  df %>%
  ggplot(aes(x = annotation_v5, y = per_1000, fill = age_group)) +
  geom_col()+
  #  geom_col(position = position_dodge(0.8))+
  scale_fill_manual(name = "Group", values = my_colors)  +
  # ylab('Frequency/1000 cells') +
  labs(title = 'in vitro nasal ALI cultures')+
  ylab('Frequency/1000 cells') +
  xlab('') +
  theme_bw(base_size = 14) + theme(legend.position = "none") +
  facet_wrap(~age_group, ncol = 1, strip.position = "right")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
plot1



#=====================================================================================
#
#  Code chunk Extended Figure 3d
#
#=====================================================================================

#df2 <- load source data from Ext Fig 3d
df2$age_group <- factor(df2$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
df2$timepoint <- factor(df2$timepoint, levels = c('4h', '24h', '72h'))
df2$Infection <- factor(df2$Infection, levels = c('SARS-CoV2', 'no viral reads'))
df2$fraction <- as.numeric(df2$fraction)
df2$fraction_1000cells <- as.numeric(df2$fraction_1000cells)

df_sum <- df2 %>%
  filter(Infection == "SARS-CoV2") %>%
  group_by(age_group, timepoint, cell_domain) %>%
  summarize(sum = sum(fraction_1000cells)/10, reads = mean(read_percell)) %>%
  ungroup()

df_sum$age_group <- factor(df_sum$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
df_sum$cell_domain <- factor(df_sum$cell_domain , levels = c('KRT5hi', 'SCGB1A1hi', 'Ciliated/Other'))

my_colors2 <- c("#fde725","#440d54", "#1e908c")            # Create vector of colors
names(my_colors2) <- levels(factor(c(levels(df_sum$cell_domain), levels(df_sum$cell_domain)))) # Extract all levels of both data

pie <-  df_sum %>%
  filter(reads >0) %>%
  ggplot(aes(x = "", y = sum, fill = cell_domain, alpha = reads)) +
  geom_col(size=.2)+
  guides(fill = guide_legend(title = "Cell domain"), alpha= guide_legend(title = "Viral reads/cell")) +
  ylab('') +
  xlab('Percent cells with viral reads') +
  labs(title = 'Viral reads/cell domain')+
  scale_alpha_continuous(range = c(0.3, 1))+
  scale_fill_manual(name = "Group", values = my_colors2)  +
  facet_wrap(~age_group+timepoint, ncol = 2)

pie +  coord_polar(theta = "y")



#=====================================================================================
#
#  Code chunk Extended Figure 3e
#
#=====================================================================================

#df2 <- load source data from Ext Fig 3d
df2$age_group <- factor(df2$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
df2$timepoint <- factor(df2$timepoint, levels = c('4h', '24h', '72h'))
df2$Infection <- factor(df2$Infection, levels = c('SARS-CoV2', 'no viral reads'))
df2$fraction_1000cells <- as.numeric(df2$fraction_1000cells)

df2$annotation_v5 <- factor(df2$annotation_v5 , levels = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1', 'Basaloid-like 2', 'Cycling basal', 'Hillock', 'Secretory', 'Secretory 2', 'Secretory 3', 'Secretory 4', 'Transit epi', 'Transit epi 2', 'Goblet 1', 'Goblet 2',
                                                         'Goblet 2 BPIFA1+', 'Goblet 2 PLAU+', 'Goblet 2 inflammatory', 'Ciliated 1', 'Ciliated 2', 'Deutorosomal', 'Ionocyte', 'Squamous'))

df2$age_group <- factor(df2$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))

fp <-  df2 %>%
  filter(fraction_1000cells > 1)%>%
  ggplot(aes(y = fraction_1000cells, x=annotation_v5, fill = Infection)) +
  geom_col(position = position_stack())+
  scale_fill_manual(values = c("dark red", "light grey"))  +
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  ylab('Percent Infected cells') +
  facet_wrap(~age_group+timepoint, ncol = 3, strip.position = "top")

fp




#=====================================================================================
#
#  Code chunk Extended Figure 3g
#
#=====================================================================================

#df3.long <- load source data from Ext Fig 3gh
df3.long$age_group <- factor(df3.long$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))

##blue, orange, gree
my_colors <- c("#450c54","#1e908c","orange")

names(my_colors) <- levels(factor(df3.long$annotation_v5 , levels = c('Basal', 'Gob/Sec','Ciliated/Other')))
df3.long$time <- factor(df3.long$time , levels = c('4h', '24h', '72h'))
df3.long$treatment <- factor(df3.long$treatment , levels = c('mock', 'COVID_19'))
df3.long$age_group <- factor(df3.long$age_group , levels = c('Paediatric', 'Adult', 'Elderly'))

library(ggpubr)
library(ggpmisc)

correlation1 <-  df3.long %>%
  filter(treatment == "COVID_19") %>%
  filter(VIRAL.SARS.CoV2 > 0 & VIRAL.SARS.CoV2 < 5)%>%
  filter(value > 0)%>%
  ggplot(aes(x = VIRAL.SARS.CoV2, y = value)) +
  geom_point()+
  xlim(0,4)+
  geom_smooth(method = "lm", formula= y~x, se=FALSE, size=1)+
  stat_cor(method = "pearson", label.x = 0.5, label.y = 4) +
  facet_wrap(~variable, ncol =  2)+
  theme_bw(base_size = 16) + theme(legend.position = "none") + theme(legend.direction = "vertical")
correlation1

#=====================================================================================
#
#  Code chunk Extended Figure 3h
#
#=====================================================================================

#see Ext Fig 3g for data processing

correlation2 <-  df3.long %>%
  filter(treatment == "COVID_19") %>%
  filter(VIRAL.SARS.CoV2 > 0 & VIRAL.SARS.CoV2 < 5)%>%
  filter(value > 0)%>%
  ggplot(aes(x = VIRAL.SARS.CoV2, y = value, color = cell_group,)) +
  geom_point()+
  scale_color_manual(name = "Cell Type", values = my_colors)  +
  xlab('SARS-CoV2 Expression') +
  ylab('Expression Level') +
  geom_smooth(method = "lm", formula= y~x, se=FALSE, size=1)+
  labs(title = "Entry Factors")+
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"))) +
  facet_wrap(~variable, scales = "free_y", ncol =  2)+
  theme_bw(base_size = 16) + theme(legend.position = "none") + theme(legend.direction = "vertical")
correlation2

proportionPlot <- df3.long %>%
  filter(treatment == "COVID_19") %>%
  filter(value > 0)%>%
  ggplot(aes(x = age_group, fill = cell_group,)) +
  geom_bar(position = position_stack())+
  scale_fill_manual(name = "Cell Type", values = my_colors)  +
  scale_x_discrete(guide = guide_axis(angle = 45))+
  ylab('Cell Count') +
  xlab('Age') +
  facet_wrap(~variable, scales = "free_y", nrow = 2)+
  theme_bw(base_size = 16) + theme(legend.position = "none")
proportionPlot

correlation2+ proportionPlot +plot_layout(byrow=1, widths = c(0.9, 0.2), guides = "auto")


#=====================================================================================
#
#  Code chunk Extended Figure 5e
#
#=====================================================================================

#df <- source data Ext Fig 5e
df$age_group <- factor(df$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))

#order
df$annotation_v5 <- factor(df$annotation_v5 , levels = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1', 'Basaloid-like 2', 'Cycling basal', 'Hillock', 'Secretory', 'Secretory 2', 'Secretory 3', 'Secretory 4', 'Transit epi', 'Transit epi 2', 'Goblet 1', 'Goblet 2',
                                                         'Goblet 2 BPIFA1+', 'Goblet 2 PLAU+', 'Goblet 2 inflammatory', 'Ciliated 1', 'Ciliated 2', 'Deutorosomal', 'Ionocyte', 'Squamous'))

plot2 <-  df %>%
  ggplot(aes(x = annotation_v5, y = per_1000, fill = treatment)) +
  geom_col(position = position_dodge())+
  ylab('Frequency/1000 cells') +
  theme(legend.direction = "vertical") +  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~age_group, ncol = 1, strip.position = "right")

plot2

#=====================================================================================
#
#  Code chunk Extended Figure 5f
#
#=====================================================================================

#df <- source data Ext Fig 5f
df$age_group <- factor(df$age_group, levels = c('Paediatric', 'Adult', 'Elderly'))
df$timepoint2 <- factor(df$timepoint2, levels = c('Mock', '4h', '24h', '72h'))
df$fold.change <- as.numeric(df$fold.change)
df$log2FC <- as.numeric(df$log2FC)

#order
df$annotation_v5 <- factor(df$annotation_v5 , levels = c('Basal 1', 'Basal 2', 'Basal|EMT1', 'Basal|EMT2', 'Basaloid-like 1', 'Basaloid-like 2', 'Cycling basal', 'Hillock', 'Secretory', 'Secretory 2', 'Secretory 3', 'Secretory 4', 'Transit epi', 'Transit epi 2', 'Goblet 1', 'Goblet 2',
                                                         'Goblet 2 BPIFA1+', 'Goblet 2 PLAU+', 'Goblet 2 inflammatory', 'Ciliated 1', 'Ciliated 2', 'Deutorosomal', 'Ionocyte', 'Squamous'))
fct_rev(df$annotation_v5)

dp <-  df %>%
  ggplot(aes(x = timepoint2, y = annotation_v5, fill = fold.change, size = per_1000)) +
  geom_point(shape=21)+
  xlab('Time') +
  scale_fill_gradient2(low = 'dark blue', high = 'dark red', mid = 'white',  midpoint = 1, limit = c(0, 2), oob = scales::squish, name = 'FC from mock')+
  theme(axis.text.x = element_text(angle = 45, hjust=1))+
  facet_wrap(~age_group, ncol = 3, strip.position = "top")
dp

#=====================================================================================
#
#  Code chunk Extended Figure 6b
#
#=====================================================================================
#install.packages("ggcorrplot")
library(ggcorrplot)
# create a correlation matrix for df_subset
df_subset <- df1 %>%
  filter(age_group == "Paediatric")%>%
  filter(treatment == "COVID_19")%>%
  subset(select = c(IFNA1, IFNB1, IFNG, IFNL1, IFNL2, IL6, ISG15, MX1, IRF7, IFITM3, RSAD2, DDX58, IFIH1, VIRAL.SARS.CoV2))

cor_mat <- round(cor(df_subset), 1)
# Compute a matrix of correlation p-values
p_mat <- cor_pmat(df_subset, p.method = "pearson")
p_mat[is.na(p_mat)] <- 1

# plot the correlation matrix with p-values using ggcorrplot
corr_plot <- ggcorrplot(cor(df_subset),
                        type = "lower",
                        title = "Paediatric COVID19 Correlation",
                        method = "square",
                        outline.color = "black",
                        insig = "pch",
                        lab = TRUE,
                        digits = 1)
corr_plot

#=====================================================================================
#
#  Code chunk Extended Figure 6d
#
#=====================================================================================

df$diffexpressed <- "NO"
df$diffexpressed[df$log2FoldChange > 0.6 & df$log10p_val  < 0.05] <- "UP"
df$diffexpressed[df$log2FoldChange < -0.6 & df$log10p_val  < 0.05] <- "DOWN"
df$delabel <- NA
df$delabel[df$diffexpressed != "NO"] <- df$Protein[df$diffexpressed != "NO"]

mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")

p <- df %>%
  ggplot(aes(x=log2FoldChange, y=-log10p_val , col=diffexpressed, fill = age_group, label=delabel)) +
  geom_point( shape = 21, size = 3) +
  theme_minimal(base_size = 18) + theme(legend.position = "right") + theme(legend.direction = "vertical") +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=c("blue", "black", "red")) +
  scale_fill_manual(name = "Group", values = my_colors2)  +
  labs(title = 'Unique Secreted Proteins - Paediatric',
       subtitle = 'DE:SARS v Mock 72h')+
  xlim(-12,12)+
  ylim(4,15)+
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=(4), col="red")

p


#=====================================================================================
#
#  Code chunk Extended Figure 7g
#
#=====================================================================================

df2$Region_Start <- as.numeric(df2$Region_Start)
df2$Age_Group <- factor(df2$Age_Group, levels = c('Paed', 'Adult', 'Elderly'))

library(ggridges)
library(ggplot2)

###Histogram
rp <- df2 %>%
  filter(Reference.allele == "No") %>%
  ggplot(aes(x=Region_Start, fill = Age_Group)) +
  geom_histogram(binwidth = 15000, colour = "black")+
  labs(
    title = '',
    y = "Count Mutations",
    x = 'Region in Genome') +
  facet_wrap(~Age_Group, ncol = 1)
rp

###Polyline
rp2 <- df2 %>%
  filter(Reference.allele == "No") %>%
  ggplot(aes(x=Region_Start, colour = Age_Group)) +
  geom_freqpoly(binwidth = 100, alpha = 1, size = 0.2)+
  labs(
    title = 'Mutations position in genome',
    y = "Count",
    x = 'Region in Genome')+
  facet_wrap(~Age_Group, ncol = 1)

rp2

rp +rp2


#=====================================================================================
#
#  Code chunk Extended Figure 8a
#
#=====================================================================================

dotplot <- df4 %>%
  ggplot(aes(x=Treatment, y = value, color = Treatment)) +
  geom_point() +
  geom_boxplot()+
  ylab('Abundance ITGB6') +
  theme(axis.ticks = element_blank()) +
  facet_wrap(~Age_group)

dotplot
