# Install required R packages
install.packages(c("tidyverse", "dplyr", "broom", "ggplot2", "lme4", "multcompView", 
                   "cowplot", "ggpubr", "glue", "rio", "rstatix"))

# Load the installed packages
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
library(rstatix)


setwd("/projects/b1059/projects/JB/2024_test")

#Read in prepared data
bt2<- read.csv("data/20240113_Assay2processed.csv")
bt1<- read.csv("data/20240113_Assay1processed.csv")

####Assay 2####
#single deletions#
single <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")) %>%
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")))
  

cont <- single %>%
  dplyr::filter(concentration_um==0)

statscont <- cont %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

#control#
sfig1 <- single %>%
 dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3257" = "pink","ECA3260"="pink" ,"VC364" = "lightgreen", "ECA3595" = "#AED6F1","ECA3600"="#AED6F1", "ECA3275" = "#D2B4DE","ECA3584" = "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584"),
                   labels = c("WT", "ben-1", "mec-7","mec-7", "tbb-1", "tbb-4","tbb-4", "tbb-6","tbb-6")) +
  stat_pvalue_manual(statscont2, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig1
ggsave("plots/sfig1.png", plot = sfig1, device = "png", width = 8, height = 4, units = "in",dpi=300)


#Treated#
#main#
singtrt <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","VC364","ECA3595","ECA3275")) %>%
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","ECA3257","VC364","ECA3595","ECA3275")))



trt <- singtrt %>%
  dplyr::filter(concentration_um==50)

statstrt <- trt %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

#trt#
Fig1 <- singtrt %>%
  dplyr::filter(concentration_um == "50") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3257" = "pink","VC364" = "lightgreen", "ECA3595" = "#AED6F1", "ECA3275" = "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3257","VC364","ECA3595","ECA3275"),
                   labels = c("WT", "ben-1", "mec-7", "tbb-1", "tbb-4","tbb-6")) +
  stat_pvalue_manual(statstrt, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

Fig1
ggsave("plots/Fig1.png", plot = Fig1, device = "png", width = 8, height = 4, units = "in",dpi=300)




#supp trt#
singtrt2 <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")) %>%
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")))

trt2 <- singtrt2 %>%
  dplyr::filter(concentration_um==50)

statstrt2 <- trt2 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig2 <- singtrt2 %>%
  dplyr::filter(concentration_um == "50") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3257" = "pink","ECA3260"="pink" ,"VC364" = "lightgreen", "ECA3595" = "#AED6F1","ECA3600"="#AED6F1", "ECA3275" = "#D2B4DE","ECA3584" = "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584"),
                   labels = c("WT", "ben-1", "mec-7","mec-7", "tbb-1", "tbb-4","tbb-4", "tbb-6","tbb-6")) +
  stat_pvalue_manual(statstrt2, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig2
ggsave("plots/sfig2.png", plot = sfig2, device = "png", width = 8, height = 4, units = "in",dpi=300)



######Double#####
double <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")))

cont1 <- double %>%
  dplyr::filter(concentration_um==0)

statscont3 <- cont1 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

#control#
sfig3<- double %>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3726"="pink","ECA3581"="pink","ECA3693"= "lightgreen","ECA3694"= "lightgreen","ECA3585"= "#AED6F1","ECA3628"= "#AED6F1","ECA3582"= "#D2B4DE","ECA3583"= "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583"),
                   labels = c("WT", "ben-1", "ben1\nmec-7","ben1\nmec-7", "ben1\ntbb-1","ben1\ntbb-1", "ben1\ntbb-4","ben1\ntbb-4", "ben1\ntbb-6","ben1\ntbb-6")) +
  stat_pvalue_manual(statscont3, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig3
ggsave("plots/sfig3.png", plot = sfig3, device = "png", width = 8, height = 4, units = "in",dpi=300)


#Treated#
#main#
dbtrt <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3693","ECA3585","ECA3582")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3726","ECA3693","ECA3585","ECA3582")))


trt3 <- dbtrt %>%
  dplyr::filter(concentration_um==50)

statstrt3 <- trt3 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

#trt#
Fig2 <- dbtrt %>%
  dplyr::filter(concentration_um == "50") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3726"="pink","ECA3693"= "lightgreen","ECA3585"= "#AED6F1","ECA3582"= "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3726","ECA3693","ECA3585","ECA3582"),
                   labels = c("WT", "ben-1", "ben-1\nmec-7", "ben-1\ntbb-1", "ben-1\ntbb-4", "ben-1\ntbb-6")) +
  stat_pvalue_manual(statstrt3, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

Fig2
ggsave("plots/Fig2.png", plot = Fig2, device = "png", width = 8, height = 4, units = "in",dpi=300)




#supp trt#
double2 <- bt2 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")))

trt4 <- double2 %>%
  dplyr::filter(concentration_um==50)

statstrt4 <- trt4 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig4 <- double2 %>%
  dplyr::filter(concentration_um == "50") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3726"="pink","ECA3581"="pink","ECA3693"= "lightgreen","ECA3694"= "lightgreen","ECA3585"= "#AED6F1","ECA3628"= "#AED6F1","ECA3582"= "#D2B4DE","ECA3583"= "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583"),
                   labels = c("WT", "ben-1", "ben1\nmec-7","ben1\nmec-7", "ben1\ntbb-1","ben1\ntbb-1", "ben1\ntbb-4","ben1\ntbb-4", "ben1\ntbb-6","ben1\ntbb-6")) +
  stat_pvalue_manual(statstrt4, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig4
ggsave("plots/sfig4.png", plot = sfig4, device = "png", width = 8, height = 4, units = "in",dpi=300)



##########Assay 1 All supplement#############

sa1 <- bt1 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")) %>%
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")))
test <- bt1 %>%
  dplyr::filter(strain=="ECA3726")

contsa1 <- sa1 %>%
  dplyr::filter(concentration_um==0)

statsconta1 <- contsa1 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

#control#
sfig5 <- sa1%>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3257" = "pink","ECA3260"="pink" ,"VC364" = "lightgreen", "ECA3595" = "#AED6F1","ECA3600"="#AED6F1", "ECA3275" = "#D2B4DE","ECA3584" = "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584"),
                   labels = c("WT", "ben-1", "mec-7","mec-7", "tbb-1", "tbb-4","tbb-4", "tbb-6","tbb-6")) +
  stat_pvalue_manual(statsconta1, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig5
ggsave("plots/sfig5.png", plot = sfig5, device = "png", width = 8, height = 4, units = "in",dpi=300)


#supp trt#
singa1trt2 <- bt1 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")) %>%
  dplyr::mutate(strain=factor(strain, levels = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584")))

trtsa1 <- singa1trt2 %>%
  dplyr::filter(concentration_um==50)

statstrta1 <- trtsa1 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig6 <- singa1trt2 %>%
  dplyr::filter(concentration_um == "50") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3257" = "pink","ECA3260"="pink" ,"VC364" = "lightgreen", "ECA3595" = "#AED6F1","ECA3600"="#AED6F1", "ECA3275" = "#D2B4DE","ECA3584" = "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3257","ECA3260","VC364","ECA3595","ECA3600","ECA3275","ECA3584"),
                   labels = c("WT", "ben-1", "mec-7","mec-7", "tbb-1", "tbb-4","tbb-4", "tbb-6","tbb-6")) +
  stat_pvalue_manual(statstrta1, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig6
ggsave("plots/sfig6.png", plot = sfig6, device = "png", width = 8, height = 4, units = "in",dpi=300)


##double Assay1##

doublea1 <- bt1 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")))

contas12 <- doublea1 %>%
  dplyr::filter(concentration_um==0)

statscontas3 <- contas12  %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

#control#
sfig7<- doublea1 %>%
  dplyr::filter(concentration_um == "0") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3726"="pink","ECA3581"="pink","ECA3693"= "lightgreen","ECA3694"= "lightgreen","ECA3585"= "#AED6F1","ECA3628"= "#AED6F1","ECA3582"= "#D2B4DE","ECA3583"= "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583"),
                   labels = c("WT", "ben-1", "ben1\nmec-7","ben1\nmec-7", "ben1\ntbb-1","ben1\ntbb-1", "ben1\ntbb-4","ben1\ntbb-4", "ben1\ntbb-6","ben1\ntbb-6")) +
  stat_pvalue_manual(statscontas3, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig7
ggsave("plots/sfig7.png", plot = sfig7, device = "png", width = 8, height = 4, units = "in",dpi=300)

#supp trt#
double2as <- bt1 %>%
  dplyr::filter(strain %in% c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")) %>%
  dplyr::mutate(strain=factor(strain,levels=c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583")))

trtas4 <- double2as %>%
  dplyr::filter(concentration_um==50)

statstrtas4 <- trtas4 %>%
  aov(median_wormlength_um_reg_delta ~ strain, data = .)%>%
  rstatix::tukey_hsd() %>%
  dplyr::filter(group1=="N2")

sfig8 <- double2as %>%
  dplyr::filter(concentration_um == "50") %>% 
  ggplot +
  aes(x = strain, y = median_wormlength_um_reg_delta) +
  geom_jitter(width = 0.1, size = 0.3) +
  geom_boxplot(aes(fill = strain), alpha = 0.8, outlier.shape = NA) +
  xlab(" ") +
  ylab("Normalized Animal Length") +
  scale_fill_manual(values = c("N2" = "#FFA500", "ECA882" = "red", "ECA3726"="pink","ECA3581"="pink","ECA3693"= "lightgreen","ECA3694"= "lightgreen","ECA3585"= "#AED6F1","ECA3628"= "#AED6F1","ECA3582"= "#D2B4DE","ECA3583"= "#D2B4DE")) +
  scale_x_discrete(breaks = c("N2","ECA882","ECA3726","ECA3581","ECA3693","ECA3694","ECA3585","ECA3628","ECA3582","ECA3583"),
                   labels = c("WT", "ben-1", "ben1\nmec-7","ben1\nmec-7", "ben1\ntbb-1","ben1\ntbb-1", "ben1\ntbb-4","ben1\ntbb-4", "ben1\ntbb-6","ben1\ntbb-6")) +
  stat_pvalue_manual(statstrtas4, label = "p.adj.signif", y.position = c(90), xmax = "group2", remove.bracket = TRUE) +
  cowplot::theme_cowplot(12) +
  theme(plot.background = element_rect(fill = "white"),
        legend.position = "None",
        plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        axis.text.x = element_text(size = 12, face = "italic"))

sfig8
ggsave("plots/sfig8.png", plot = sfig8, device = "png", width = 8, height = 4, units = "in",dpi=300)

