install.packages("devtools")
devtools::install_github("AndersenLab/easyXpress",force=TRUE)
library(tidyverse)
library(dplyr)
library(broom)
library(ggplot2)
library(lme4)
library(multcompView)
library(cowplot)
library(ggpubr)
library(glue)
library(rio)

setwd("/projects/b1059/projects/JB/2024_bt")

#Read in prepared data
bt1<- read.csv("data/20240206_Assay1processed.csv")
bt2<- read.csv("data/20240206_Assay2processed.csv")


##################
######Assay1######
##################

#Divide into single and double deletions
sbt1 <- bt1 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")) %>%
  dplyr::mutate(strain=factor(strain, levels=c("N2","VC364","ECA3257","ECA3260","ECA3595","ECA3600","ECA882","ECA3275","ECA3584")))

dbt1 <- bt1 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628","ECA3582","ECA3583")))


####Single####

#MainFig1-Treated
sstatst <- sbt1 %>%
  dplyr::filter(concentration_um==30)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

stF1<- sstatst %>%
  dplyr::filter(group2 %in% c("N2","ECA882","ECA3257","VC364","ECA3595","ECA3275"))

Fig1 <- sbt1 %>%
  dplyr::filter(concentration_um == "30") %>% 
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","VC364","ECA3595","ECA3275")) %>%
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","VC364" = "grey", "ECA3257" = "grey" , "ECA3595" = "grey",  "ECA3275" = "grey")) +
  scale_x_discrete(breaks = c(c("N2","ECA882","VC364","ECA3257","ECA3595","ECA3275")),
                   labels = c("WT","ben-1", "tbb-1","mec-7","tbb-4","tbb-6")) +
  stat_pvalue_manual(stF1, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
Fig1
ggsave("plots/Fig3.png", plot = Fig1, device = "png", width = 8, height = 4, units = "in",dpi=300)

#SFig1-Full Trt
SFig1 <- sbt1 %>%
  dplyr::filter(concentration_um == "30") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","VC364" = "grey", "ECA3257" = "grey","ECA3260"="grey" , "ECA3595" = "grey","ECA3600"="grey",  "ECA3275" = "grey","ECA3584" = "grey")) +
  scale_x_discrete(breaks = c(c("N2","ECA882","VC364","ECA3257","ECA3260","ECA3595","ECA3600","ECA3275","ECA3584")),
                   labels = c("WT","ben-1(ean64)", "tbb-1(gk207)","mec-7(ean253)","mec-7(ean254)","tbb-4(ean263)","tbb-4(ean268)","tbb-6(ean255)","tbb-6(ean261)")) +
  stat_pvalue_manual(sstatst, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
SFig1
ggsave("plots/Sfig3.png", plot = SFig1, device = "png", width = 8, height = 4, units = "in",dpi=300)


#Controls
sstatsc <- sbt1 %>%
  dplyr::filter(concentration_um==0)%>%
  aov(median_wormlength_um ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig2 <- sbt1 %>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","VC364" = "grey", "ECA3257" = "grey","ECA3260"="grey" , "ECA3595" = "grey","ECA3600"="grey",  "ECA3275" = "grey","ECA3584" = "grey")) +
  scale_x_discrete(breaks = c(c("N2","ECA882","VC364","ECA3257","ECA3260","ECA3595","ECA3600","ECA3275","ECA3584")),
                   labels = c("WT","ben-1(ean64)" ,"tbb-1(gk207)","mec-7(ean253)","mec-7(ean254)","tbb-4(ean263)","tbb-4(ean268)","tbb-6(ean255)","tbb-6(ean261)")) +
  stat_pvalue_manual(sstatsc, label = "p.adj.signif", y.position = c(900), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig2
ggsave("plots/sfig4.png", plot = sfig2, device = "png", width = 8, height = 5, units = "in",dpi=300)

####Double####

#Control
dstatsc <- dbt1 %>%
  dplyr::filter(concentration_um==0)%>%
  aov(median_wormlength_um ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig7 <- dbt1 %>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","ECA3693"= "grey","ECA3694"= "grey","ECA3581"="grey","ECA3585"= "grey","ECA3628"= "grey",  "ECA3582"= "grey","ECA3583"= "grey")) +
  scale_x_discrete(breaks = c(c("N2" ,"ECA882" ,"ECA3693","ECA3694","ECA3581","ECA3585","ECA3628",  "ECA3582","ECA3583")),
                   labels = c("WT","ben-1(ean64)", "ben-1(ean64) tbb-1(ean285)","ben-1(ean64) tbb-1(ean286)", "ben-1(ean64); mec-7(ean258)",  "ben-1(ean64); tbb-4(ean262)","ben-1(ean64); tbb-4(ean282)", "ben-1(ean64); tbb-6(ean259)","ben-1(ean64); tbb-6(ean260)")) +
  stat_pvalue_manual(dstatsc, label = "p.adj.signif", y.position = c(900), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig7
ggsave("plots/SFig9.png", plot = sfig7, device = "png", width = 8, height = 5, units = "in",dpi=300)

#Treated
#levels for stats to work as needed
dbt1_2 <- dbt1 %>%
 dplyr::mutate(strain=factor(strain,levels=c("ECA882","N2","ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628","ECA3582","ECA3583")))


dstatst <- dbt1_2 %>%
  dplyr::filter(concentration_um==30)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="ECA882")
dstatst <- dstatst %>%
  dplyr::mutate(group2=factor(group2,levels=c("N2","ECA882","ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628","ECA3582","ECA3583")))


sfig8 <- dbt1 %>%
  dplyr::filter(concentration_um == "30") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","ECA3693"= "grey","ECA3694"= "grey","ECA3581"="grey","ECA3585"= "grey","ECA3628"= "grey",  "ECA3582"= "grey","ECA3583"= "grey")) +
  scale_x_discrete(breaks = c(c("N2" ,"ECA882" ,"ECA3693","ECA3694","ECA3581","ECA3585","ECA3628",  "ECA3582","ECA3583")),
                   labels = c("WT","ben-1(ean64)", "ben-1(ean64) tbb-1(ean285)","ben-1(ean64) tbb-1(ean286)", "ben-1(ean64); mec-7(ean258)",  "ben-1(ean64); tbb-4(ean262)","ben-1(ean64); tbb-4(ean282)", "ben-1(ean64); tbb-6(ean259)","ben-1(ean64); tbb-6(ean260)")) +
  stat_pvalue_manual(dstatst, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig8
ggsave("plots/SFig10.png", plot = sfig8, device = "png", width = 8, height = 5, units = "in",dpi=300)






##################
######Assay2######
##################

#Divide into single and double deletions
sbt2 <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")) %>%
  dplyr::mutate(strain=factor(strain, levels=c("N2","VC364","ECA3257","ECA3260","ECA3595","ECA3600","ECA882","ECA3275","ECA3584")))

dbt2 <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628","ECA3582","ECA3583")))


####Single####

#Sfig3-Treated
sstatst2 <- sbt2 %>%
  dplyr::filter(concentration_um==30)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig3 <- sbt2 %>%
  dplyr::filter(concentration_um == "30") %>% 
 ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","VC364" = "grey", "ECA3257" = "grey","ECA3260"="grey" , "ECA3595" = "grey","ECA3600"="grey",  "ECA3275" = "grey","ECA3584" = "grey")) +
  scale_x_discrete(breaks = c(c("N2","ECA882","VC364","ECA3257","ECA3260","ECA3595","ECA3600","ECA3275","ECA3584")),
                   labels = c("WT","ben-1(ean64)" ,"tbb-1(gk207)","mec-7(ean253)","mec-7(ean254)","tbb-4(ean263)","tbb-4(ean268)","tbb-6(ean255)","tbb-6(ean261)")) +
  stat_pvalue_manual(sstatst2, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig3
ggsave("plots/SFig5.png", plot = sfig3, device = "png", width = 8, height = 4, units = "in",dpi=300)

#Controls-sfig4
sstatsc2 <- sbt2 %>%
  dplyr::filter(concentration_um==0)%>%
  aov(median_wormlength_um ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig4 <- sbt2 %>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","VC364" = "grey", "ECA3257" = "grey","ECA3260"="grey" , "ECA3595" = "grey","ECA3600"="grey",  "ECA3275" = "grey","ECA3584" = "grey")) +
  scale_x_discrete(breaks = c(c("N2","ECA882","VC364","ECA3257","ECA3260","ECA3595","ECA3600","ECA3275","ECA3584")),
                   labels = c("WT","ben-1(ean64)" ,"tbb-1(gk207)","mec-7(ean253)","mec-7(ean254)","tbb-4(ean263)","tbb-4(ean268)","tbb-6(ean255)","tbb-6(ean261)")) +
  stat_pvalue_manual(sstatsc2, label = "p.adj.signif", y.position = c(900), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig4
ggsave("plots/SFig6.png", plot = sfig4, device = "png", width = 8, height = 5, units = "in",dpi=300)

####Double####

#Control
dstatsc2 <- dbt2 %>%
  dplyr::filter(concentration_um==0)%>%
  aov(median_wormlength_um ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig5 <- dbt2 %>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","ECA3693"= "grey","ECA3694"= "grey","ECA3726"="grey","ECA3581"="grey","ECA3585"= "grey","ECA3628"= "grey",  "ECA3582"= "grey","ECA3583"= "grey")) +
  scale_x_discrete(breaks = c(c("N2" ,"ECA882" ,"ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628",  "ECA3582","ECA3583")),
                   labels = c("WT","ben-1(ean64)", "ben-1(ean64) tbb-1(ean285)","ben-1(ean64) tbb-1(ean286)", "ben-1(ean64); mec-7(ean257)","ben-1(ean64); mec-7(ean258)",  "ben-1(ean64); tbb-4(ean262)","ben-1(ean64); tbb-4(ean282)", "ben-1(ean64); tbb-6(ean259)","ben-1(ean64); tbb-6(ean260)")) +
  stat_pvalue_manual(dstatsc2, label = "p.adj.signif", y.position = c(900), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig5
ggsave("plots/SFig8.png", plot = sfig5, device = "png", width = 8, height = 5, units = "in",dpi=300)

#Treated-full
#change levels for appropriate stats
dbt2_2 <- dbt2 %>%
  dplyr::mutate(strain=factor(strain,levels=c("ECA882","N2","ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628","ECA3582","ECA3583")))


dstatst2 <- dbt2_2 %>%
  dplyr::filter(concentration_um==30)%>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="ECA882")
dstatst2 <- dstatst2 %>%
  dplyr::mutate(group2=factor(group2,levels=c("N2","ECA882","ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628","ECA3582","ECA3583")))

sfig6 <- dbt2 %>%
  dplyr::filter(concentration_um == "30") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","ECA3693"= "grey","ECA3694"= "grey","ECA3726"="grey","ECA3581"="grey","ECA3585"= "grey","ECA3628"= "grey",  "ECA3582"= "grey","ECA3583"= "grey")) +
  scale_x_discrete(breaks = c(c("N2" ,"ECA882" ,"ECA3693","ECA3694","ECA3726","ECA3581","ECA3585","ECA3628",  "ECA3582","ECA3583")),
                   labels = c("WT","ben-1(ean64)", "ben-1(ean64) tbb-1(ean285)","ben-1(ean64) tbb-1(ean286)", "ben-1(ean64); mec-7(ean257)","ben-1(ean64); mec-7(ean258)",  "ben-1(ean64); tbb-4(ean262)","ben-1(ean64); tbb-4(ean282)", "ben-1(ean64); tbb-6(ean259)","ben-1(ean64); tbb-6(ean260)")) +
  stat_pvalue_manual(dstatst2, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
sfig6
ggsave("plots/SFig7.png", plot = sfig6, device = "png", width = 8, height = 5, units = "in",dpi=300)

#Treated-Fig2
ds2 <- dstatst2 %>%
  dplyr::filter(group2 %in% c("N2","ECA882","ECA3726","ECA3693","ECA3585","ECA3582"))

Fig2 <- dbt2 %>%
  dplyr::filter(concentration_um == "30") %>% 
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3693","ECA3585","ECA3582")) %>%
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized animal length") +
  scale_fill_manual(values = c("N2" = "#FFA500","ECA882" = "red","ECA3693"= "grey","ECA3726"="grey","ECA3585"= "grey", "ECA3582"= "grey")) +
  scale_x_discrete(breaks = c(c("N2" ,"ECA882" ,"ECA3693","ECA3726","ECA3585",  "ECA3582")),
                   labels = c("WT","ben-1", "ben-1 tbb-1", "ben-1; mec-7",  "ben-1; tbb-4", "ben-1; tbb-6")) +
  stat_pvalue_manual(ds2, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, margin = unit(c(0, 0, 0, 0), units = "in"), face = "italic"),
        plot.background = element_rect(fill="white"),
        strip.background = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
  )
Fig2
ggsave("plots/Fig4.png", plot = Fig2, device = "png", width = 8, height = 5, units = "in",dpi=300)
