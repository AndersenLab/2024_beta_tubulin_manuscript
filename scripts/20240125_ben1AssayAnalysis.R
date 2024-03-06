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

today <- format(Sys.Date(), "%Y%m%d")
setwd("/projects/b1059/projects/JB/2024_bt")


##########################################
##############Assay 1#####################
##########################################
today <- format(Sys.Date(), "%Y%m%d")
###HTA Data Processing and Cleaning####
# figure directory 
figure_dir <- "/projects/b1059/projects/JB/2024_bt/plots"
dirs <- "/projects/b1059/projects/JB/2024_bt/data/Assay1" #data and design match
data <- "20230803_TCben1_Analysis-20230810.RData"

bt1 <- easyXpress::readXpress(filedir = dirs,
                               design = T,
                               rdafile = "20230803_TCben1_Analysis-20230810.RData") 

ms<- easyXpress::modelSelection(bt1$raw_data)

# Use the edgeOF fucntion
ef <- easyXpress::edgeOF(data = ms)

# Use the clusterOF function
cf <- easyXpress::clusterOF(data = ef)
#Check flags
c1 <- easyXpress::checkOF(data = cf, drug, concentration_um)
c1$p
ggsave(filename = "plots/Assay1flags.png",c1$p, device=png)
#CheckObjecs
c2 <- easyXpress::checkObjs(data = cf, OF = 'filter', drug, concentration_um)
c2
ggsave(filename = "plots/Assay1obje.png",c2, device=png)
##Make photo arrays to check models
cm <- cf%>%
  mutate(
    i.dir = "/projects/b1059/projects/Skyler/cellprofiler-nf/projects/20230803_TCben1/Analysis-20230810/processed_images/",
    w.lab = paste(drug, strain, sep = "_"))

# Now we can run the checkModels function
cm_out <- easyXpress::checkModels(data = cm ,
                                  # the grouping vars (...), make a plot for each.
                                  drug,
                                  proc.img.dir = "i.dir",
                                  well.label = "w.lab",
                                  # save in the repo you cloned
                                  out.dir = glue::glue("{figure_dir}/{today}_checkModelsassay1"))


# add the user variable that will be converted to an object flag
u = cm %>%
  dplyr::mutate(user = dplyr::case_when(drug == "Albendazole" &
                                          model == "MDHD" ~ "junk",
                                        drug == "Albendazole" &
                                          worm_length_um < 165 ~ "junk",
                                        drug == "DMSO" &
                                          model == "MDHD" ~ "junk",
                                        drug == "DMSO" &
                                          worm_length_um < 165 ~ "junk",
                                        TRUE ~ NA_character_))
# Run the userOF function and specify user variable as the flag
uf <- easyXpress::userOF(data = u, user)

# Check the object data again to see if the bimodal distributions are resolved.
easyXpress::checkObjs(data = uf, OF = "filter", drug, concentration_um)
#Apply classifier
cl <- easyXpress::classifierOF(data = uf)
#Remove Outliers
o <- easyXpress::outlierOF(data = cl)
#Use the checkObjs() function again to check the effect of filtering all the object flags. 
z<- easyXpress::checkObjs(data = o, OF = 'filter', drug, concentration_um)
z
ggsave(filename = "plots/Assay1flagsrm.png",z, device=png)
#Use the checkOF function again to see how much data is being flagged.
co2 <- easyXpress::checkOF(data = o, drug, concentration_um)
co2$p
ggsave(filename="plots/assay1presumflag.png",co2$p,device=png)
#Filter all flaggs and summarize
proc.objs <- easyXpress::filterOF(o, rmVars = T)
raw.wells <- easyXpress::summarizeWells(data = o)
#Check titering
tf <- easyXpress::titerWF(data = raw.wells,
                          Metadata_Experiment, bleach, strain, drug,                          doseR = F)
tf$p
ggsave(filename = "plots/Assay1titer.png",tf$p, device=png)
#Check wells for too many/few objects
n <- easyXpress::nWF(data = tf$d, drug, concentration_um, max = 30, min = 5)
n$p
ggsave(filename = "plots/Assay1wellobjects.png",n$p, device=png)
#Check well outliers
ow <- easyXpress::outlierWF(data = n$d,
                            Metadata_Experiment, bleach, drug,
                            concentration_um, strain) %>%
  dplyr::mutate(assay_bleach = paste(Metadata_Experiment, bleach, sep = "_"))
#Check how many wells are flagged
cw1 <- easyXpress::checkWF(data = ow, drug, concentration_um)

cw1$p
ggsave(filename = "plots/Assay1wellflags.png",cw1$p, device=png)
#Filter bad wells
fw <- easyXpress::filterWF(data = ow, rmVars = T)
#Check how many wells are retained

cb <- easyXpress::checkBalance(data = fw, drug, concentration_um,
                               design = bt1$design %>%
                                 dplyr::mutate(assay_bleach =
                                                 paste(Metadata_Experiment, bleach)),
                               x = assay_bleach)
poo<- cb$p +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
  ggplot2::geom_hline(yintercept = 0.75, linetype = 2, color = "red")
ggsave(filename = "plots/Assay1wellretained.png",poo, device=png)
#set all to ABZ 0 and 30
zz <- fw %>%
  mutate(drug = ifelse(drug == "DMSO", "Albendazole", drug))

zz<-zz %>% mutate(genotype = case_when(
  strain == "N2" ~ "WT",
  strain == "ECA882" ~ "ben-1(ean64)",
  strain == "VC364" ~ "tbb-1(gk207)",
  strain == "ECA3257" ~ "mec-7(ean253)",
  strain == "ECA3260" ~ "mec-7(ean254)",
  strain == "ECA3595" ~ "tbb-4(ean263)",
  strain == "ECA3600" ~ "tbb-4(ean268)",
  strain == "ECA3275" ~ "tbb-6(ean255)",
  strain == "ECA3584" ~ "tbb-6(ean261)",
  strain == "ECA3582" ~ "ben-1(ean64); tbb-6(ean259)",
  strain == "ECA3583" ~ "ben-1(ean64); tbb-6(ean260)",
  strain == "ECA3726" ~ "ben-1(ean64); mec-7(ean257)",
  strain == "ECA3694" ~ "ben-1(ean64) tbb-1(ean286)",
  strain == "ECA3693" ~ "ben-1(ean64) tbb-1(ean285)",
  strain == "ECA3581" ~ "ben-1(ean64); mec-7(ean258)",
  strain == "ECA3585" ~ "ben-1(ean64); tbb-4(ean262)",
  strain == "ECA3628" ~ "ben-1(ean64); tbb-4(ean282)",
  TRUE ~ NA_character_
))

#Check effects
ce1 <- chekEff(data = zz,genotype,
                            x = concentration_um,
                            y = median_wormlength_um,
                            fill = assay_bleach,
                            scales = "free_x")
ce1 <- ce1 + theme(strip.text = element_text(face = "italic"))


ggsave(filename = "plots/SuppFig2.png",ce1, device=png)
##Regression,Delta, and final check
# Regress the effect of independent bleaches for each drug using regEff()
reg <- easyXpress::regEff(data = zz,
                          drug,
                          d.var = median_wormlength_um,
                          c.var = assay_bleach)

# Look at the regression coefficients in the diagnostic plot p2
reg$p2
ggsave(filename = "plots/Assay1regco.png",reg$p2, device=png)
del<- easyXpress::delta(data = reg$d,
                        assay_bleach, drug, strain, # group with ...
                        WF = "filter",
                        doseR = TRUE,
                        vars = "median_wormlength_um_reg")
del
# check the finalized data
a<- easyXpress::checkEff(data = del, drug, strain, x = concentration_um,
                         y = median_wormlength_um_reg_delta,
                         fill = assay_bleach,
                         scales = "free_x")
ggsave(filename = "plots/Assay1finaleff.png",a, device=png)
#Final Output
write.csv(del, "data/20240206_Assay1processed.csv")










##############################################
###################Assay2#####################
##############################################
today <- format(Sys.Date(), "%Y%m%d")
###HTA Data Processing and Cleaning####
# figure directory 
figure_dir <- "/projects/b1059/projects/JB/2024_bt/plots"
dirs <- "/projects/b1059/projects/JB/2024_bt/data/Assay2" #data and design match
data <- "20230929_TCben1_Analysis-20231106.RData"

bt2 <- easyXpress::readXpress(filedir = dirs,
                              design = T,
                              rdafile = "20230929_TCben1_Analysis-20231106.RData") 


ms<- easyXpress::modelSelection(bt2$raw_data)

# Use the edgeOF fucntion
ef <- easyXpress::edgeOF(data = ms)

# Use the clusterOF function
cf <- easyXpress::clusterOF(data = ef)
#Check flags
c1 <- easyXpress::checkOF(data = cf, drug, concentration_um)
c1$p
ggsave(filename = "plots/Assay2checkflags.png",c1$p, device=png)
#CheckObjecs
c2 <- easyXpress::checkObjs(data = cf, OF = 'filter', drug, concentration_um)
c2
ggsave(filename = "plots/Assay2checkonjs.png",c2, device=png)
##Make photo arrays to check models
cm <- cf%>%
  mutate(
    i.dir = "/projects/b1059/projects/Skyler/cellprofiler-nf/projects/20230929_TCben1/Analysis-20231106/processed_images/",
    w.lab = paste(drug, strain, sep = "_"))

# Now we can run the checkModels function
cm_out <- easyXpress::checkModels(data = cm ,
                                  # the grouping vars (...), make a plot for each.
                                  drug,
                                  proc.img.dir = "i.dir",
                                  well.label = "w.lab",
                                  # save in the repo you cloned
                                  out.dir = glue::glue("{figure_dir}/{today}_Assay2checkModels"))


# add the user variable that will be converted to an object flag
u = cm %>%
  dplyr::mutate(user = dplyr::case_when(drug == "Albendazole" &
                                          model == "MDHD" ~ "junk",
                                        drug == "Albendazole" &
                                          worm_length_um < 165 ~ "junk",
                                        drug == "DMSO" &
                                          model == "MDHD" ~ "junk",
                                        drug == "DMSO" &
                                          worm_length_um < 165 ~ "junk",
                                        TRUE ~ NA_character_))
# Run the userOF function and specify user variable as the flag
uf <- easyXpress::userOF(data = u, user)


# Check the object data again to see if the bimodal distributions are resolved.
easyXpress::checkObjs(data = uf, OF = "filter", drug, concentration_um)
#Apply classifier
cl <- easyXpress::classifierOF(data = uf)
#Remove Outliers
o <- easyXpress::outlierOF(data = cl)
#Use the checkObjs() function again to check the effect of filtering all the object flags. 
z<- easyXpress::checkObjs(data = o, OF = 'filter', drug, concentration_um)
z
ggsave(filename = "plots/Assay2filtobjc.png",z, device=png)
#Use the checkOF function again to see how much data is being flagged.
co2 <- easyXpress::checkOF(data = o, drug, concentration_um)
co2$p
ggsave(filename="plots/assay2presumflag.png",co2$p,device=png)
#Filter all flaggs and summarize
proc.objs <- easyXpress::filterOF(o, rmVars = T)
raw.wells <- easyXpress::summarizeWells(data = o)
#Check titering
tf <- easyXpress::titerWF(data = raw.wells,
                          Metadata_Experiment, bleach, strain, drug,                          doseR = F)
tf$p
ggsave(filename = "plots/Assay2titer.png",tf$p, device=png)
#Check wells for too many/few objects
n <- easyXpress::nWF(data = tf$d, drug, concentration_um, max = 30, min = 5)
n$p
ggsave(filename = "plots/Assay2wellobjects.png",n$p, device=png)
#Check well outliers
ow <- easyXpress::outlierWF(data = n$d,
                            Metadata_Experiment, bleach, drug,
                            concentration_um, strain) %>%
  dplyr::mutate(assay_bleach = paste(Metadata_Experiment, bleach, sep = "_"))
#Check how many wells are flagged
cw1 <- easyXpress::checkWF(data = ow, drug, concentration_um)

cw1$p
ggsave(filename = "plots/Assay2wellsflagged.png",cw1$p, device=png)
#Filter bad wells
fw <- easyXpress::filterWF(data = ow, rmVars = T)


cb <- easyXpress::checkBalance(data = fw, drug, concentration_um,
                               design = bt2$design %>%
                                 dplyr::mutate(assay_bleach =
                                                 paste(Metadata_Experiment, bleach)),
                               x = assay_bleach)
poo2 <- cb$p +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) +
  ggplot2::geom_hline(yintercept = 0.75, linetype = 2, color = "red")
ggsave(filename = "plots/Assay2wellsretained.png",poo2, device=png)
#set all to ABZ 0 and 30
zz <- fw %>%
  mutate(drug = ifelse(drug == "DMSO", "Albendazole", drug))
zz<-zz %>% mutate(genotype = case_when(
  strain == "N2" ~ "WT",
  strain == "ECA882" ~ "ben-1(ean64)",
  strain == "VC364" ~ "tbb-1(gk207)",
  strain == "ECA3257" ~ "mec-7(ean253)",
  strain == "ECA3260" ~ "mec-7(ean254)",
  strain == "ECA3595" ~ "tbb-4(ean263)",
  strain == "ECA3600" ~ "tbb-4(ean268)",
  strain == "ECA3275" ~ "tbb-6(ean255)",
  strain == "ECA3584" ~ "tbb-6(ean261)",
  strain == "ECA3582" ~ "ben-1(ean64); tbb-6(ean259)",
  strain == "ECA3583" ~ "ben-1(ean64); tbb-6(ean260)",
  strain == "ECA3726" ~ "ben-1(ean64); mec-7(ean257)",
  strain == "ECA3694" ~ "ben-1(ean64) tbb-1(ean286)",
  strain == "ECA3693" ~ "ben-1(ean64) tbb-1(ean285)",
  strain == "ECA3581" ~ "ben-1(ean64); mec-7(ean258)",
  strain == "ECA3585" ~ "ben-1(ean64); tbb-4(ean262)",
  strain == "ECA3628" ~ "ben-1(ean64); tbb-4(ean282)",
  TRUE ~ NA_character_
))
#Check effects
ce1 <- chekEff(data = zz,genotype,
               x = concentration_um,
               y = median_wormlength_um,
               fill = assay_bleach,
               scales = "free_x")

ce1
ggsave(filename = "plots/Suppfig1.png",ce1, device=png)

##Regression,Delta, and final check
# Regress the effect of independent bleaches for each drug using regEff()
reg <- easyXpress::regEff(data = zz,
                          drug,
                          d.var = median_wormlength_um,
                          c.var = assay_bleach)

# Look at the regression coefficients in the diagnostic plot p2
reg$p2
ggsave(filename = "plots/Assay2bleachf.png",reg$p2, device=png)
del<- easyXpress::delta(data = reg$d,
                        assay_bleach, drug, strain, # group with ...
                        WF = "filter",
                        doseR = TRUE,
                        vars = "median_wormlength_um_reg")
del
# check the finalized data
a<- easyXpress::checkEff(data = del, drug, strain, x = concentration_um,
                         y = median_wormlength_um_reg_delta,
                         fill = assay_bleach,
                         scales = "free_x")
ggsave(filename = "plots/Assay2finaleffs.png",a ,device=png)
#Final Output
write.csv(del, "data/20240206_Assay2processed.csv")





# ####COMBINED PROCESSED####
# 
# 
# bt1 <- read.csv("data/20240113_Assay1processed.csv")
# bt2 <- read.csv("data/20240113_Assay2processed.csv")
# 
# fw_assay1_clean <- bt1 %>%
#   filter(!is.na(drug))  # Remove rows where drug is NA
# 
# fw_assay2_clean <- bt2 %>%
#   filter(!is.na(drug))  # Remove rows where drug is NA
# 
# HTA_combined_ABZ <- bind_rows(fw_assay1_clean, fw_assay2_clean)
# 
# # use the checkEff function
# ce1_assay_ABZ <- easyXpress::checkEff(data = HTA_combined_ABZ,
#                                       drug, strain,
#                                       x = concentration_um,
#                                       y = median_wormlength_um,
#                                       fill = assay_bleach,
#                                       scales = "free_x")
# 
# # look at the plot
# ce1_assay_ABZ
# 
# # Assuming your data frame is named del_assay_HTA
# HTA_combined_ABZ <- HTA_combined_ABZ %>%
#   mutate(assay_bleach = ifelse(assay_bleach == "TCben1_1", "1_1", assay_bleach))
# 
# HTA_combined_ABZ <- HTA_combined_ABZ %>%
#   mutate(assay_bleach = ifelse(assay_bleach == "TCben1_2", "1_2", assay_bleach))
# # Assuming your data frame is named del_assay_HTA
# HTA_combined_ABZ <- HTA_combined_ABZ %>%
#   mutate(assay_bleach = ifelse(assay_bleach == "TCben1_3", "1_3", assay_bleach))
# 
# HTA_combined_ABZ <- HTA_combined_ABZ %>%
#   mutate(assay_bleach = ifelse(assay_bleach == "ben1_1", "2_1", assay_bleach))
# # Assuming your data frame is named del_assay_HTA
# HTA_combined_ABZ <- HTA_combined_ABZ %>%
#   mutate(assay_bleach = ifelse(assay_bleach == "ben1_2", "2_2", assay_bleach))
# 
# HTA_combined_ABZ <- HTA_combined_ABZ %>%
#   mutate(assay_bleach = ifelse(assay_bleach == "ben1_3", "2_3", assay_bleach))
# 
# 
# zz <- HTA_combined_ABZ %>%
#   mutate(drug = ifelse(drug == "DMSO", "Albendazole", drug))
# 
# # regEff
# # Regress the effect of independent bleaches for each drug using regEff()
# reg_assay <- easyXpress::regEff(data = HTA_combined_ABZ, #drop
#                                 drug,
#                                 d.var = median_wormlength_um,
#                                 c.var = assay_bleach)
# 
# # Look at the regression coefficients in the diagnostic plot p2
# reg_assay$p2
# 
# 
# # delta 
# # use the delta() function
# del_assay <- easyXpress::delta(data = reg_assay$d,
#                                assay_bleach, drug, strain, # group with ...
#                                WF = "filter",
#                                doseR = TRUE,
#                                vars = "median_wormlength_um_reg")
# 
# View(del_assay)
# 
# # check the finalized data
# easyXpress::checkEff(data = del_assay, drug, strain, x = concentration_um,
#                      y = median_wormlength_um_reg_delta,
#                      fill = assay_bleach,
#                      scales = "free_x")
# 
# 
# write.csv(del_assay,"data/20240131_combinedben1.csv")