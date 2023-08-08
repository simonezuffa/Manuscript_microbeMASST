# Set working directory
setwd("~/PycharmProjects/Manuscript_microbeMASST")

# Packages
library(tidyverse)
library(ggpubr)
library(UpSetR)
library(Spectra)
library(msdata)
library(MsBackendMgf)
library(taxize)
library(VennDiagram)
library(ComplexUpset)
library(patchwork)

# Data obtained from classical molecular networking (https://gnps.ucsd.edu/ProteoSAFe/status.jsp?task=893fd89b52dc4c07a292485404f97723)
data_3dmap <- read.delim("MSV000079949/feature_table_MSV000079949.tsv") %>% 
  as.data.frame() %>% column_to_rownames("FeatureID") %>% t() %>% as.data.frame() %>% rownames_to_column("SampleID")
data_3dmap$SampleID <- gsub(".mzXML", "", data_3dmap$SampleID)

metadata_3dmap <- data_3dmap %>% dplyr::select(SampleID) %>% 
  separate(SampleID, c("info1", "info2", "info3", "info4", "info5"), remove = FALSE)

mn_annotation <- read.delim("MSV000079949/mn_annotation_MSV000079949.tsv")
mn_annotation$X.Scan. <- as.character(mn_annotation$X.Scan.)

mol_network <- read.delim("MSV000079949/network.txt") # singletons (ComponentIndex = -1) are present in the feature table

# Extract MSMS only present in SPF animals from the different organs/biofluids of interest

############################
# GASTRO-INTESTINAL REGION #
############################

cecum <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Cecum")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

cecum_check <- cecum %>% dplyr::mutate(SumAll = rowSums(cecum)) %>%
  dplyr::mutate(GF = rowSums(cecum[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(cecum[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 4911 features found in cecum

colon <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Colon")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

colon_check <- colon %>% dplyr::mutate(SumAll = rowSums(colon)) %>%
  dplyr::mutate(GF = rowSums(colon[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(colon[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 5717 features found in colon

duodenum <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Duo")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean))  %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

duodenum_check <- duodenum %>% dplyr::mutate(SumAll = rowSums(duodenum)) %>%
  dplyr::mutate(GF = rowSums(duodenum[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(duodenum[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 5881 features found in duodenum

ileum <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Ile")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean))  %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

ileum_check <- ileum %>% dplyr::mutate(SumAll = rowSums(ileum)) %>%
  dplyr::mutate(GF = rowSums(ileum[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(ileum[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 5475 features found in ileum

jejunum <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Jeju")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean))  %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

jejunum_check <- jejunum %>% dplyr::mutate(SumAll = rowSums(jejunum)) %>%
  dplyr::mutate(GF = rowSums(jejunum[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(jejunum[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 5981 features found in jejunum

stool <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Stool")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean))  %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

stool_check <- stool %>% dplyr::mutate(SumAll = rowSums(stool)) %>%
  dplyr::mutate(GF = rowSums(stool[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(stool[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 1189 features found in stool

stomach <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Stomach")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

stomach_check <- stomach %>% dplyr::mutate(SumAll = rowSums(stomach)) %>%
  dplyr::mutate(GF = rowSums(stomach[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(stomach[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 3298 features found in stomach

esophagus <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Eso")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

esophagus_check <- esophagus %>% dplyr::mutate(SumAll = rowSums(esophagus)) %>%
  dplyr::mutate(GF = rowSums(esophagus[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(esophagus[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 1617 features found in esophagus

mouth <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Mouth")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

mouth_check <- mouth %>% dplyr::mutate(SumAll = rowSums(mouth)) %>%
  dplyr::mutate(GF = rowSums(mouth[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(mouth[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 875 features found in mouth

# Extract metabolites exclusively present in SPF animals

cecum_spf <- cecum_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Cecum = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Cecum)

colon_spf <- colon_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Colon = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Colon)

duodenum_spf <- duodenum_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Duodenum = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Duodenum)

ileum_spf <- ileum_check %>% dplyr::filter(GF == 0 & SPF > 0) %>% 
  dplyr::mutate(Ileum = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Ileum)

jejunum_spf <- jejunum_check %>% dplyr::filter(GF == 0 & SPF > 0) %>% 
  dplyr::mutate(Jejunum = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Jejunum)

stool_spf <- stool_check %>% dplyr::filter(GF == 0 & SPF > 0) %>% 
  dplyr::mutate(Stool = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Stool)

stomach_spf <- stomach_check %>% dplyr::filter(GF == 0 & SPF > 0) %>% 
  dplyr::mutate(Stomach = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Stomach)

esophagus_spf <- esophagus_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Esophagus = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Esophagus)

mouth_spf <- mouth_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Mouth = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Mouth)

# Check if the feature are retained between regions
features_gi <- mouth_spf %>%
  dplyr::full_join(esophagus_spf) %>%
  dplyr::full_join(stomach_spf) %>%
  dplyr::full_join(duodenum_spf) %>%
  dplyr::full_join(jejunum_spf) %>%
  dplyr::full_join(ileum_spf) %>%
  dplyr::full_join(cecum_spf) %>%
  dplyr::full_join(colon_spf) %>%
  dplyr::full_join(stool_spf) %>%
  dplyr::mutate(Region = str_count(paste(Mouth, Esophagus, Stomach, Duodenum, Jejunum, Ileum, Cecum, Colon, Stool, sep = " "), "Yes")) %>%
  arrange(desc(Region)) %>%
  left_join(mn_annotation, by = c("Feature" = "X.Scan."))


################################################################################
# ORGANS of interest (Blood, Brain, Heart, Kid, Liver, Lung, Pancreas, Spleen) #
################################################################################

blood <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Blood")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

blood_check <- blood %>% dplyr::mutate(SumAll = rowSums(blood)) %>%
  dplyr::mutate(GF = rowSums(blood[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(blood[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 1098 features found in blood

brain <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Brain")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

brain_check <- brain %>% dplyr::mutate(SumAll = rowSums(brain)) %>%
  dplyr::mutate(GF = rowSums(brain[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(brain[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 2392 features found in brain

heart <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Heart")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

heart_check <- heart %>% dplyr::mutate(SumAll = rowSums(heart)) %>%
  dplyr::mutate(GF = rowSums(heart[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(heart[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 1953 features found in heart

liver <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Liver")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

liver_check <- liver %>% dplyr::mutate(SumAll = rowSums(liver)) %>%
  dplyr::mutate(GF = rowSums(liver[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(liver[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 2621 features found in liver

kidney <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Kid")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

kidney_check <- kidney %>% dplyr::mutate(SumAll = rowSums(kidney)) %>%
  dplyr::mutate(GF = rowSums(kidney[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(kidney[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 2008 features found in kidney

lung <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Lung")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

lung_check <- lung %>% dplyr::mutate(SumAll = rowSums(lung)) %>%
  dplyr::mutate(GF = rowSums(lung[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(lung[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 2039 features found in lung

pancreas <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Pancreas")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

pancreas_check <- pancreas %>% dplyr::mutate(SumAll = rowSums(pancreas)) %>%
  dplyr::mutate(GF = rowSums(pancreas[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(pancreas[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 1821 features found in pancreas

spleen <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Spleen")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

spleen_check <- spleen %>% dplyr::mutate(SumAll = rowSums(spleen)) %>%
  dplyr::mutate(GF = rowSums(spleen[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(spleen[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 1715 features found in spleen

# Extract metabolites exclusively present in SPF animals

blood_spf <- blood_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Blood = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Blood)

brain_spf <- brain_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Brain = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Brain)

heart_spf <- duodenum_check %>% dplyr::filter(GF == 0 & SPF > 0) %>% 
  dplyr::mutate(Heart = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Heart)

pancreas_spf <- pancreas_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Pancreas = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Pancreas)

liver_spf <- liver_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Liver = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Liver)

kidney_spf <- kidney_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Kidney = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Kidney)

lung_spf <- lung_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Lung = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Lung)

spleen_spf <- spleen_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Spleen = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Spleen)

# Check if the feature are retained between organs
features_organs <- blood_spf %>%
  dplyr::full_join(brain_spf) %>%
  dplyr::full_join(heart_spf) %>%
  dplyr::full_join(pancreas_spf) %>%
  dplyr::full_join(liver_spf) %>%
  dplyr::full_join(kidney_spf) %>%
  dplyr::full_join(lung_spf) %>%
  dplyr::full_join(spleen_spf) %>%
  dplyr::mutate(Region = str_count(paste(Blood, Brain, Heart, Pancreas, Liver, Kidney, Lung, Spleen, sep = " "), "Yes")) %>%
  arrange(desc(Region)) %>%
  left_join(mn_annotation, by = c("Feature" = "X.Scan."))


#############################################################
# FEMALE REPRODUCTIVE TRACT (Vagina, Uterus, Ovary, Cervix) #
#############################################################

vagina <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Vagina")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

vagina_check <- vagina %>% dplyr::mutate(SumAll = rowSums(vagina)) %>%
  dplyr::mutate(GF = rowSums(vagina[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(vagina[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 979 features found in vagina

cervix <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Cervix")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

cervix_check <- cervix %>% dplyr::mutate(SumAll = rowSums(cervix)) %>%
  dplyr::mutate(GF = rowSums(cervix[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(cervix[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 674 features found in cervix

uterus <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Uterus")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

uterus_check <- uterus %>% dplyr::mutate(SumAll = rowSums(uterus)) %>%
  dplyr::mutate(GF = rowSums(uterus[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(uterus[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 2054 features found in uterus

ovary <- data_3dmap %>% dplyr::filter(str_detect(SampleID, pattern = "Ovary")) %>% 
  left_join((metadata_3dmap %>% dplyr::select(SampleID, info1))) %>% 
  column_to_rownames("SampleID") %>% group_by(info1) %>%
  summarise(across(everything(), mean)) %>% 
  column_to_rownames("info1") %>%
  t() %>% as.data.frame()

ovary_check <- ovary %>% dplyr::mutate(SumAll = rowSums(ovary)) %>%
  dplyr::mutate(GF = rowSums(ovary[, 1:4])) %>%
  dplyr::mutate(SPF = rowSums(ovary[, 5:8])) %>%
  dplyr::filter(SumAll > 0) # 934 features found in ovary

# Extract metabolites exclusively present in SPF animals

vagina_spf <- vagina_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Vagina = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Vagina)

cerix_spf <- cervix_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Cervix = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Cervix)

uterus_spf <- uterus_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Uterus = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Uterus)

ovary_spf <- ovary_check %>% dplyr::filter(GF == 0 & SPF > 0) %>%
  dplyr::mutate(Ovary = "Yes") %>% rownames_to_column("Feature") %>%
  dplyr::select(Feature, Ovary)

# Check if the feature are retained between female reproductive tract
features_female <- vagina_spf %>%
  dplyr::full_join(cerix_spf) %>%
  dplyr::full_join(uterus_spf) %>%
  dplyr::full_join(ovary_spf) %>%
  dplyr::mutate(Region = str_count(paste(Vagina, Cervix, Uterus, Ovary, sep = " "), "Yes")) %>%
  arrange(desc(Region)) %>%
  left_join(mn_annotation, by = c("Feature" = "X.Scan."))


########################################
# Agglomerate all features of interest #
########################################
features_interest <- features_gi %>% 
  full_join(features_organs, by = "Feature") %>% 
  full_join(features_female, by = "Feature") %>%
  distinct(Feature, .keep_all = TRUE)

features_interest_select <- features_interest %>% dplyr::select(Feature, Mouth, Esophagus, Stomach, Duodenum, Jejunum, Ileum, Cecum, Colon, Stool,
                                                                Blood, Brain, Heart, Pancreas, Liver, Kidney, Lung, Spleen, Vagina, Cervix, Uterus, Ovary) %>%
  dplyr::mutate(Region = str_count(paste(Mouth, Esophagus, Stomach, Duodenum, Jejunum, Ileum, Cecum, Colon, Stool,
                                         Blood, Brain, Heart, Pancreas, Liver, Kidney, Lung, Spleen, Vagina, Cervix, Uterus, Ovary, sep = " "), "Yes")) %>%
  arrange(desc(Region)) %>%
  left_join(mn_annotation, by = c("Feature" = "X.Scan."))

#write_csv(x = features_interest_select, file = "MSMS_SPF_Location.csv")



# List of MSMS of interest found only in SPF mice
feature_microbial_possible <- data.frame(Feature = features_interest_select$Feature) %>% arrange(Feature)

# Check network and remove feature that have an edge with a feature only present in GF animals and that have the same parent ion mass
mol_network$CLUSTERID1 <- as.character(mol_network$CLUSTERID1)
mol_network$CLUSTERID2 <- as.character(mol_network$CLUSTERID2)

mol_network_expand <- mol_network %>% 
  dplyr::mutate(CLUSTERID1_colonized = case_when(CLUSTERID1 %in% feature_microbial_possible$Feature ~ "Yes",
                                                 TRUE ~ "No")) %>% 
  dplyr::mutate(CLUSTERID2_colonized = case_when(CLUSTERID2 %in% feature_microbial_possible$Feature ~ "Yes",
                                                 TRUE ~ "No")) %>%
  dplyr::mutate(Matching = case_when(CLUSTERID1_colonized == "Yes" & CLUSTERID2_colonized == "Yes" ~ "No",
                                     CLUSTERID1_colonized == "No" & CLUSTERID2_colonized == "No" ~ "No",
                                     TRUE ~ "Yes"))

matching_02 <- mol_network_expand %>% dplyr::filter(Matching == "Yes") %>% dplyr::filter(abs(DeltaMZ) < 0.02)
features_to_remove <- data.frame(Features = c(matching_02$CLUSTERID1, matching_02$CLUSTERID2)) %>% distinct(Features)

reduntant_spf <- mol_network_expand %>% dplyr::filter(Matching == "No") %>% dplyr::filter(abs(DeltaMZ) < 0.02) %>%
  dplyr::filter(CLUSTERID1_colonized == "Yes") %>% dplyr::filter(CLUSTERID1 != CLUSTERID2) # possible redundant MSMS within SPF

######################################################
# Filter MGF file obtained from molecular networking #
######################################################

# Following code takes time (run only if necessary)
dda <- Spectra("MSV000079949/Quinn_Mouse_Organs_GFvSPF_MSV000079949_MN/METABOLOMICS-SNETS-V2-893fd89b-download_clustered_spectra-main.mgf", source = MsBackendMgf()) # the mgf file is too big to be uploaded into GitHub. You can downlaod it from the link to the MN job
features_interest_select$Feature <- as.numeric(features_interest_select$Feature)
dda_filtered <- dda[features_interest_select$Feature] # put vector with indices to filter specific scans
#export(dda_filtered, MsBackendMgf(), file = "microbemasst_10380.mgf", exportTitle = FALSE)

# The filtered mgf file (10380 MS/MS spectra) is run in batch with microbeMASST (takes approximately 2 hours)
# Filter for the matching edges with GF spectra will be applied later


######################################################
# READ ALL RESULTS (.tsv) OBTAINED WITH microbeMASST #
######################################################

# Output of microbeMASST is too big for GitHub (21GB). Here the code used to generate 
# the microbemasst_output_MSV000079949.csv table that can be read to move on with analysis

# function to read in a tsv and add the file name as a column
#customized_read_tsv <- function(file){
#  read_tsv(file, col_types = cols(.default = "c")) %>%
#    mutate(fileName = file)}

#microbemasst_output <- list.files(path = "microbeMASST_MSV000079949", full.names = TRUE, pattern = "counts_microbe.tsv") %>% # list all the files
#  lapply(customized_read_tsv) %>% # read them all in with our custom function
#  reduce(bind_rows) %>% # stack them all on top of each other
#  dplyr::select(fileName, Taxaname_file,	Taxaname_alternative,	Taxa_NCBI) # select the correct columns

#write_csv(microbemasst_output, "microbemasst_output_MSV000079949.csv")
microbemasst_output <- read_csv("microbemasst_output_MSV000079949.csv")

microbemasst_output$fileName <- gsub("microbeMASST_MSV000079949/mgf__", "", microbemasst_output$fileName)
microbemasst_output$fileName <- gsub("_counts_microbe.tsv", "", microbemasst_output$fileName)

# Check how many unique MS/MS had a match in microbeMASST
microbemasst_output_distinct <- microbemasst_output %>% distinct(fileName, .keep_all = TRUE) # 3785 spectra were found with microbeMASST of the 10380 that were searched

# Remove MSMS matching to GF MSMS
microbemasst_output_distinct_filter <- microbemasst_output_distinct %>% dplyr::filter(!(fileName %in% features_to_remove$Features)) # 3262 MSMS spectra after filtering
microbemasst_output_distinct_filter_annotated <- microbemasst_output_distinct_filter %>% left_join(mn_annotation, by = c("fileName" = "X.Scan."))

#write_csv(microbemasst_output_distinct_filter_annotated, "microbemasst_3262.csv")

dda_3262 <- dda[as.numeric(microbemasst_output_distinct_filter_annotated$fileName)] # put vector with indices to filter specific scans
#export(dda_3262, MsBackendMgf(), file = "microbemasst_3262.mgf", exportTitle = FALSE)

# Check MSMS found also in human cell lines
fetures_human <- microbemasst_output %>% dplyr::filter(Taxaname_file == "Homo sapiens") %>% distinct(fileName, .keep_all = TRUE) # 968 spectra were found in human cell lines
fetures_human_annotated <- fetures_human %>% left_join(mn_annotation, by = c("fileName" = "X.Scan."))

#write_csv(fetures_human_annotated, "microbemasst_human_968.csv")

# Remove spectra that were found in human cell lines
unique_microbial <- microbemasst_output %>% dplyr::filter(!(fileName %in% fetures_human$fileName)) %>%
  dplyr::filter(!(fileName %in% features_to_remove$Features)) %>% dplyr::filter(!(Taxa_NCBI == "blank" | Taxa_NCBI == "qc")) %>%
  distinct(fileName, .keep_all = TRUE) %>% left_join(mn_annotation, by = c("fileName" = "X.Scan.")) %>%
  dplyr::filter(is.na(Compound_Name) | !(str_detect(Compound_Name, pattern = regex("ampi|amitr|MS_Contaminant|Sodium", ignore_case = TRUE)))) 

# 2425 unique MS/MS spectra that only had match to microbial cultures

#write_csv(unique_microbial, "microbemasst_2425.csv")

dda_2425 <- dda[as.numeric(unique_microbial$fileName)] # put vector with indices to filter specific scans
#export(dda_2425, MsBackendMgf(), file = "microbemasst_2425.mgf", exportTitle = FALSE)

taxon_apperance <- microbemasst_output %>% dplyr::filter(fileName %in% unique_microbial$fileName) %>% 
  group_by(Taxa_NCBI) %>% summarise(n = n()) %>% arrange(desc(n)) %>% dplyr::filter(!(Taxa_NCBI %in% c("blank", "qc")))

# Classify nodes tree
Sys.setenv(ENTREZ_KEY = "enter_here_your_api_key")

microbe_lineage <- classification(taxon_apperance$Taxa_NCBI, db = 'ncbi', batch_size = 4) # it can return HTTP error - if so re run it
microbe_lineage_long <- do.call(rbind, microbe_lineage)

microbe_lineage_long_rooted <- microbe_lineage_long %>% 
  dplyr::mutate(rank_fix = case_when(name == "cellular organisms" ~ "root",
                                     str_detect(name, "unclassified ") ~ "unclassified",
                                     TRUE ~ rank))

microbe_lineage_wide <- microbe_lineage_long_rooted %>% rownames_to_column("Sample") %>% 
  dplyr::filter(rank_fix %in% c("root", "no rank", "superkingdom", "kingdom", "phylum",
                                "class", "order", "family", "genus", "subgenus",
                                "species", "subspecies", "strain", "varietas")) %>%
  dplyr::filter(!(str_detect(name, "incerta")))

microbe_lineage_wide$Sample <- gsub("\\..*", "", microbe_lineage_wide$Sample)
microbe_lineage_wide_final <- microbe_lineage_wide %>% pivot_wider(id_cols = "Sample", names_from = "rank_fix", values_from = "name")

# Merge with 
taxon_apperance_info <- taxon_apperance %>% left_join(microbe_lineage_wide_final, by = c("Taxa_NCBI" = "Sample"))



##############################################################
# Create Venn diagram of features between bacteria and fungi #
##############################################################
unique_microbial_info <- microbemasst_output %>% 
  dplyr::filter(!(fileName %in% fetures_human$fileName)) %>%
  dplyr::filter(!(fileName %in% features_to_remove$Features)) %>%
  left_join(mn_annotation, by = c("fileName" = "X.Scan.")) %>%
  dplyr::filter(is.na(Compound_Name) | !(str_detect(Compound_Name, pattern = regex("ampi|amitr|MS_Contaminant|Sodium", ignore_case = TRUE)))) %>%
  left_join(microbe_lineage_wide_final, by = c("Taxa_NCBI" = "Sample")) %>% 
  group_by(fileName) %>% distinct(superkingdom, .keep_all = TRUE) %>%
  dplyr::filter(!(Taxa_NCBI == "blank" | Taxa_NCBI == "qc")) %>%
  dplyr::select(1:6, 52:65) # 2425 unique MS/MS found in 3082 superkindgom cultures cause some are produced by both bacteria and fungi

list_ven <- list(
  Bacteria = (unique_microbial_info %>% dplyr::filter(superkingdom == "Bacteria"))$fileName,
  Fungi = (unique_microbial_info %>% dplyr::filter(superkingdom == "Eukaryota"))$fileName)

# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

display_venn(list_ven, fill = c("#999999", "#E69F00"), cex = 1)



################################################################################################
# Generate an MGF file of these features of interest and run it with SIRIUS and also make a MN #
################################################################################################
foi_microbial <- as.numeric(distinct(unique_microbial_info, fileName, .keep_all = TRUE)$fileName)
dda_interest_sirius <- dda[foi_microbial]
#export(dda_interest_sirius, MsBackendMgf(), file = "Fetures_interest_MSV000079949_SIRIUS.mgf", exportTitle = FALSE)



########################
# Import SIRIUS output #
########################
sirius <- read_delim("SIRIUS_MSV000079949/canopus_compound_summary.tsv")
sirius$id <- gsub(".+_", "", sirius$id)
sirius_fix <- sirius %>% dplyr::select(id, `NPC#class`, `NPC#pathway`, `ClassyFire#most specific class`)
colnames(sirius_fix) <- c("fileName", "NPC_Class", "NPC_Pathway", "ClassyFire_Class")

# Features found in bacterial, fungi, only bacteria, and only fungi
bacterial <- unique_microbial_info %>% dplyr::filter(superkingdom == "Bacteria")
fungal <- unique_microbial_info %>% dplyr::filter(superkingdom == "Eukaryota")

bacterial_only <- unique_microbial_info %>% dplyr::filter(fileName %in% bacterial$fileName) %>% 
  dplyr::filter(!(fileName %in% fungal$fileName)) %>% left_join(sirius, by = c("fileName" = "id"))
fungal_only <- unique_microbial_info %>% dplyr::filter(fileName %in% fungal$fileName) %>% 
  dplyr::filter(!(fileName %in% bacterial$fileName)) %>% left_join(sirius, by = c("fileName" = "id"))

# Add SIRIUS info to the full table
unique_microbial_sirius <- unique_microbial %>% dplyr::select(fileName, Compound_Name, Adduct, Precursor_MZ, 
                                                              ExactMass, LibraryQualityString) %>% left_join(sirius_fix)

#write.csv(unique_microbial_sirius, "microbemasst_2425_sirius.csv")

#######################
# GENERATE UPSET PLOT #
#######################

# Check all bacterial features
bacterial_features <- microbemasst_output %>% 
  dplyr::filter(fileName %in% bacterial$fileName) %>%
  left_join(microbe_lineage_wide_final, by = c("Taxa_NCBI" = "Sample")) %>%
  dplyr::filter(superkingdom == "Bacteria") %>%
  group_by(phylum) %>%
  distinct(fileName, .keep_all = TRUE) %>%
  left_join(mn_annotation, by = c("fileName" = "X.Scan.")) %>%
  dplyr::filter(!(is.na(phylum)))

# I have remove MSMS spectra that do not have a phylum assosciated cause samples were generically categorized as bacteria or fungi 
# This leaves 1920 MSMS out of the 2330 that were classified as bacteria

summary(as.factor(bacterial_features$phylum))

list_ven_bacteria_small <- list(
  Actinomycetota = (bacterial_features %>% dplyr::filter(phylum == "Actinomycetota"))$fileName,
  Bacillota = (bacterial_features %>% dplyr::filter(phylum == "Bacillota"))$fileName,
  Bacteroidota = (bacterial_features %>% dplyr::filter(phylum == "Bacteroidota"))$fileName,
  Cyanobacteriota = (bacterial_features %>% dplyr::filter(phylum == "Cyanobacteriota"))$fileName,
  Fusobacteriota = (bacterial_features %>% dplyr::filter(phylum == "Fusobacteriota"))$fileName,
  Pseudomonadota = (bacterial_features %>% dplyr::filter(phylum == "Pseudomonadota"))$fileName,
  Verrucomicrobiota = (bacterial_features %>% dplyr::filter(phylum == "Verrucomicrobiota"))$fileName)

upset_plot_bacterial <- UpSetR::upset(fromList(list_ven_bacteria_small), nsets = 7, nintersects = NA, point.size = 1.5, line.size = 1, text.scale = 0.5)

# Use package complexupset to plot sirius annotations
complex_bacteria <- bacterial_features %>% dplyr::filter(phylum %in% c("Actinomycetota", "Bacillota", "Bacteroidota", "Cyanobacteriota", 
                                                                       "Fusobacteriota", "Pseudomonadota", "Verrucomicrobiota")) %>%
  dplyr::select(fileName, phylum) %>% dplyr::mutate(Presence = 1) %>%
  pivot_wider(id_cols = fileName, names_from = phylum, values_from = Presence) %>%
  mutate_all(~replace_na(., 0)) %>%
  left_join(sirius_fix) %>% 
  dplyr::mutate(adjusted = case_when(NPC_Class == "Penicillins" ~ ClassyFire_Class, 
                                     NPC_Class == "Cephalosporins" ~ ClassyFire_Class,
                                     NPC_Class == "Erythromycins" ~ ClassyFire_Class,
                                     TRUE ~ NPC_Class))

top_annotaions <- complex_bacteria %>%
  group_by(adjusted) %>% 
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(30) %>% 
  left_join(sirius_fix, by = c("adjusted" = "NPC_Class")) %>% distinct(adjusted, .keep_all = TRUE)

# create new column using case_when function
complex_bacteria <- complex_bacteria %>%
  mutate(final = case_when(
    adjusted %in% top_annotaions$adjusted ~ adjusted,
    TRUE ~ "wOther" ))

phyla <- c("Actinomycetota", "Bacillota", "Bacteroidota", "Cyanobacteriota", "Fusobacteriota", "Pseudomonadota", "Verrucomicrobiota")

upset_bac <- upset(complex_bacteria, phyla,
      base_annotations = list('Intersection size'= intersection_size(counts = FALSE, mapping = aes(fill = NPC_Pathway))),
      width_ratio = 0.1, height_ratio = 0.3,
      sort_intersections = "ascending",
      sort_sets = "descending",
      sort_intersections_by=c('degree', 'cardinality')) & scale_fill_manual(values = c("#005f73", "#0a9396", "#94d2bd", "#e9d8a6",
                                                                                       "#ee9b00", "#ca6702", "#bb3e03", "#ae2012"))

upset_bac + plot_layout(guides='collect') &
  theme(legend.position='none')

# The two upset plots are combined using Illustrator to obtain final figure in the manuscript

# Check all fungal features
fungal_features <- microbemasst_output %>% 
  dplyr::filter(fileName %in% fungal$fileName) %>%
  left_join(microbe_lineage_wide_final, by = c("Taxa_NCBI" = "Sample")) %>%
  dplyr::filter(superkingdom == "Eukaryota") %>%
  group_by(phylum) %>%
  distinct(fileName, .keep_all = TRUE) %>%
  left_join(mn_annotation, by = c("fileName" = "X.Scan.")) %>%
  dplyr::filter(!(is.na(phylum)))

# I have remove MSMS spectra that do not have a phylum associated cause samples were generically categorized as bacteria or fungi 
# This leaves 628 MSMS out of the 752 that were classified as bacteria

summary(as.factor(fungal_features$phylum))

list_fungal_phylum <- list(
  Ascomycota = (fungal_features %>% dplyr::filter(phylum == "Ascomycota"))$fileName,
  Bacillariophyta = (fungal_features %>% dplyr::filter(phylum == "Bacillariophyta"))$fileName,
  Basidiomycota = (fungal_features %>% dplyr::filter(phylum == "Basidiomycota"))$fileName,
  Euglenozoa = (fungal_features %>% dplyr::filter(phylum == "Euglenozoa"))$fileName,
  Mucoromycota = (fungal_features %>% dplyr::filter(phylum == "Mucoromycota"))$fileName,
  Zoopagomycota = (fungal_features %>% dplyr::filter(phylum == "Zoopagomycota"))$fileName)

upset_plot_fungal <- UpSetR::upset(fromList(list_fungal_phylum), nsets = 6, nintersects = NA, point.size = 1.5, line.size = 1, text.scale = 0.5)

# Use package complexupset to plot sirius annotations
complex_fungal <- fungal_features %>% dplyr::filter(phylum %in% c("Ascomycota", "Bacillariophyta", "Basidiomycota", "Euglenozoa", 
                                                                       "Mucoromycota", "Zoopagomycota")) %>%
  dplyr::select(fileName, phylum) %>% dplyr::mutate(Presence = 1) %>%
  pivot_wider(id_cols = fileName, names_from = phylum, values_from = Presence) %>%
  mutate_all(~replace_na(., 0)) %>%
  left_join(sirius_fix) %>% 
  dplyr::mutate(adjusted = case_when(NPC_Class == "Penicillins" ~ ClassyFire_Class, 
                                     NPC_Class == "Cephalosporins" ~ ClassyFire_Class,
                                     NPC_Class == "Erythromycins" ~ ClassyFire_Class,
                                     TRUE ~ NPC_Class))

top_annotaions_fungal <- complex_fungal %>%
  group_by(adjusted) %>% 
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(30)

# create new column using case_when function
complex_fungal <- complex_fungal %>%
  mutate(final = case_when(
    adjusted %in% top_annotaions_fungal$adjusted ~ adjusted,
    TRUE ~ "wOther"
  ))

phyla_fungi <- c("Ascomycota", "Bacillariophyta", "Basidiomycota", "Euglenozoa", "Mucoromycota", "Zoopagomycota")

upset_fun <- upset(complex_fungal, phyla_fungi,
      base_annotations = list('Intersection size'= intersection_size( counts = FALSE, mapping = aes(fill = NPC_Pathway))),
      width_ratio = 0.1, height_ratio = 0.3,
      sort_intersections = "ascending",
      sort_sets = "descending",
      sort_intersections_by=c('degree', 'cardinality')) & scale_fill_manual(values = c("#005f73", "#0a9396", "#94d2bd", "#e9d8a6",
                                                                                       "#ee9b00", "#ca6702", "#bb3e03", "#ae2012"))

upset_fun + plot_layout(guides='collect') &
  theme(legend.position='none')


####################################################################################
# Extract MS/MS spectra that are also observed in the Karin dataset (MSV000080918) #
####################################################################################

# I am going to use the same output of microbemasst but in this case I am looking for al matches
# that I am going to filter for a specific dataset (MSV000080918)
#list_microbemasst_output <- list.files(path = "microbeMASST_MSV000079949", full.names = TRUE, pattern = "matches.tsv") %>% as.data.frame()
#microbemasst_matches <- data.frame(Feature = list_microbemasst_output$.) %>% 
#  dplyr::filter(str_detect(Feature, pattern = paste(unique_microbial_info$fileName, collapse = "|")))

#microbemasst_output_bacterial_foi <- microbemasst_matches$Feature %>%
#  lapply(customized_read_tsv) %>%
#  reduce(bind_rows) %>%
#  dplyr::select(fileName, USI)

#microbial_karin <- microbemasst_output_bacterial_foi %>% dplyr::filter(str_detect(USI, pattern = "MSV000080918")) 

#write_csv(microbial_karin, "microbemasst_output_karin.csv")
microbial_karin <- read_csv("microbemasst_output_karin.csv")

karin_split <- microbial_karin %>% separate(USI, into = c("one", "two", "three", "four", "five"), sep = ":") %>%
  separate(fileName, into = c("six", "seven", "eight", "nine"), sep = "_") %>%
  dplyr::select(nine, three, five)

colnames(karin_split) <- c("Feature", "Sample", "Scan")

karin_split_filter <- karin_split %>% dplyr::filter(Feature %in% unique_microbial_info$fileName)

# Get metadata file
karin_metadata <- read_csv("metadata_MSV000080918.csv")
karin_metadata$SampleID <- gsub("_mzML_1", "", karin_metadata$SampleID)

# Total number of features observed in karin dataset --> 621 spectra
features_karin <- karin_split_filter %>% distinct(Feature)

features_karin_whole <- unique_microbial %>% dplyr::filter(fileName %in% features_karin$Feature) %>% distinct(fileName, .keep_all = TRUE) %>%
  dplyr::select(fileName, SpectrumID, Compound_Name, Adduct, Precursor_MZ, ExactMass, LibraryQualityString) %>% left_join(sirius_fix)

#write_csv(features_karin_whole, "microbemasst_621.csv")

# How many were still observed in animals receiving the antibiotics cocktail? --> 109 spectra
features_antibiotics <- karin_split_filter %>% left_join(karin_metadata, by = c("Sample" = "SampleID")) %>%
  dplyr::filter(str_detect(ATTRIBUTE_antibiotics, pattern = "Yes")) %>% distinct(Feature, .keep_all = TRUE)

feature_karin_microbial <- features_karin %>% dplyr::filter(!(Feature %in% features_antibiotics$Feature)) # 512 MSMS spectra

# Generate list with annotation and sirius putative classification
unique_microbial_karin_sirius <- feature_karin_microbial %>% left_join(mn_annotation, by = c("Feature" = "X.Scan.")) %>%
  dplyr::select(Feature, Compound_Name) %>% left_join(sirius_fix, by = c("Feature" = "fileName"))

#write.csv(unique_microbial_karin_sirius, "microbemasst_512.csv")

# Filter the initial list for the features obtained from the karin dataset
unique_microbial_info_karin <- unique_microbial_karin_sirius %>% left_join(unique_microbial_info, by = c("Feature" = "fileName"))

list_karin <- list(
  Bacteria = (unique_microbial_info_karin %>% dplyr::filter(superkingdom == "Bacteria"))$Feature,
  Fungi = (unique_microbial_info_karin %>% dplyr::filter(superkingdom == "Eukaryota"))$Feature)

# Venn diagram
display_venn(list_karin, fill = c("#999999", "#E69F00"), cex = 1)

# Generate the bacterial UpSet plot
bacterial_features_karin <- microbemasst_output %>% 
  dplyr::filter(fileName %in% unique_microbial_info_karin$Feature) %>%
  left_join(microbe_lineage_wide_final, by = c("Taxa_NCBI" = "Sample")) %>%
  dplyr::filter(superkingdom == "Bacteria") %>%
  group_by(phylum) %>%
  distinct(fileName, .keep_all = TRUE) %>%
  left_join(mn_annotation, by = c("fileName" = "X.Scan.")) %>%
  dplyr::filter(!(is.na(phylum)))

summary(as.factor(bacterial_features_karin$phylum))

list_ven_bacteria_karin <- list(
  Actinomycetota = (bacterial_features_karin %>% dplyr::filter(phylum == "Actinomycetota"))$fileName,
  Bacillota = (bacterial_features_karin %>% dplyr::filter(phylum == "Bacillota"))$fileName,
  Bacteroidota = (bacterial_features_karin %>% dplyr::filter(phylum == "Bacteroidota"))$fileName,
  Cyanobacteriota = (bacterial_features_karin %>% dplyr::filter(phylum == "Cyanobacteriota"))$fileName,
  Fusobacteriota = (bacterial_features_karin %>% dplyr::filter(phylum == "Fusobacteriota"))$fileName,
  Pseudomonadota = (bacterial_features_karin %>% dplyr::filter(phylum == "Pseudomonadota"))$fileName,
  Verrucomicrobiota = (bacterial_features_karin %>% dplyr::filter(phylum == "Verrucomicrobiota"))$fileName)

upset_plot_karin <- UpSetR::upset(fromList(list_ven_bacteria_karin), nsets = 7, nintersects = NA, point.size = 1.5, line.size = 1, text.scale = 0.5)

# Use package complexupset to plot sirius annotations
complex_bacteria_karin <- bacterial_features_karin %>% dplyr::filter(phylum %in% c("Actinomycetota", "Bacillota", "Bacteroidota", "Cyanobacteriota", 
                                                                       "Fusobacteriota", "Pseudomonadota", "Verrucomicrobiota")) %>%
  dplyr::select(fileName, phylum) %>% dplyr::mutate(Presence = 1) %>%
  pivot_wider(id_cols = fileName, names_from = phylum, values_from = Presence) %>%
  mutate_all(~replace_na(., 0)) %>%
  left_join(sirius_fix) %>% 
  dplyr::mutate(adjusted = case_when(NPC_Class == "Penicillins" ~ ClassyFire_Class, 
                                     NPC_Class == "Cephalosporins" ~ ClassyFire_Class,
                                     NPC_Class == "Erythromycins" ~ ClassyFire_Class,
                                     TRUE ~ NPC_Class))

top_annotaions_karin <- complex_bacteria_karin %>%
  group_by(adjusted) %>% 
  summarize(count = n()) %>%
  arrange(desc(count)) %>%
  top_n(30)

# create new column using case_when function
complex_bacteria_karin <- complex_bacteria_karin %>%
  mutate(final = case_when(
    adjusted %in% top_annotaions_karin$adjusted ~ adjusted,
    TRUE ~ "wOther"
  ))

upset_bac_karin <- upset(complex_bacteria_karin, phyla,
      base_annotations = list('Intersection size'= intersection_size( counts = FALSE, mapping = aes(fill = NPC_Pathway))),
      width_ratio = 0.1, height_ratio = 0.3,
      sort_intersections = "ascending",
      sort_sets = "descending",
      sort_intersections_by=c('degree', 'cardinality')) & scale_fill_manual(values = c("#005f73", "#0a9396", "#94d2bd", "#e9d8a6",
                                                                                       "#ee9b00", "#ca6702", "#bb3e03", "#ae2012"))

upset_bac_karin + plot_layout(guides='collect') &
  theme(legend.position='none')

# What about the spectra that where not found in the karin dataset? maybe because the spectra were not initially found in stools?
location_karin <- features_interest_select %>% dplyr::filter(Feature %in% unique_microbial_info_karin$Feature) # not very interesting
location_not_karin <- features_interest_select %>% dplyr::filter(Feature %in% unique_microbial_info$fileName) %>%
  dplyr::filter(!(Feature %in% unique_microbial_info_karin$Feature))

# Create mgf of 512
foi_microbial_karin <- as.numeric(distinct(location_karin, Feature, .keep_all = TRUE)$Feature)
dda_interest_karin <- dda[foi_microbial_karin]
#export(dda_interest_karin, MsBackendMgf(), file = "microbemasst_512.mgf", exportTitle = FALSE)


#############################################################################
# Check location of the MS/MS spectra in the different body sites for human #
#############################################################################
redu_metadata <- read.delim("redu_metadata.tsv")
redu_metadata_interest <- redu_metadata %>% dplyr::filter(str_detect(NCBITaxonomy, pattern = "9606")) # 58 datasets human samples
summary(factor(redu_metadata_interest$UBERONBodyPartName))

redu_metadata_body <- redu_metadata %>% dplyr::filter(str_detect(NCBITaxonomy, pattern = "9606")) %>%
  dplyr::filter(UBERONBodyPartName != "not applicable") %>% dplyr::select(filename, ATTRIBUTE_DatasetAccession, UBERONBodyPartName, NCBITaxonomy, DOIDCommonName)
redu_metadata_body$filename <- str_extract(redu_metadata_body$filename, "[^/]*$")
redu_metadata_body$filename <- gsub(".mzML", "", redu_metadata_body$filename)
redu_metadata_body$filename <- gsub(".mzXML", "", redu_metadata_body$filename)
redu_metadata_body$UBERONBodyPartName <- gsub("blood serum", "blood plasma", redu_metadata_body$UBERONBodyPartName)

microbial_location <- microbemasst_output_bacterial_foi %>% 
  dplyr::filter(str_detect(USI, pattern = paste(redu_metadata_body$ATTRIBUTE_DatasetAccession, collapse = "|")))

#write_csv(microbial_location, "microbemasst_output_location.csv")
microbial_location <- read_csv("microbemasst_output_location.csv")

microbial_location_split <- microbial_location %>% separate(USI, into = c("one", "two", "three", "four", "five"), sep = ":") %>%
  separate(fileName, into = c("six", "seven", "eight", "nine"), sep = "_") %>%
  dplyr::select(nine, two, three, five)

colnames(microbial_location_split) <- c("Feature", "MassIVE", "Sample", "Scan")

microbial_location_split_filter <- microbial_location_split %>% dplyr::filter(Feature %in% location_karin$Feature) # 455 MSMS in 45 datasets

# Baloon Plot
microbial_location_split_human_baloon <- microbial_location_split_filter %>% 
  left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>%
  dplyr::filter(!is.na(DOIDCommonName)) %>%
  dplyr::filter(!(DOIDCommonName == "not applicable" | UBERONBodyPartName == "not applicable" | UBERONBodyPartName == "not specified" | DOIDCommonName == "ML import: not available")) %>%
  #dplyr::filter(!(DOIDCommonName == "not specified" | DOIDCommonName == "not collected" | DOIDCommonName == "disease NOS")) %>%
  dplyr::mutate(DOIDCommonName = case_when(DOIDCommonName %in% c("not specified", "not collected", "disease NOS") ~ "Not Specified",
                                           TRUE ~ DOIDCommonName)) %>%
  dplyr::mutate(Disease = case_when(
    DOIDCommonName %in% c("ulcerative colitis", "inflammatory bowel disease", "Crohn's disease") ~ "IBD",
    DOIDCommonName == "diabetes mellitus" ~ "T2D",
    DOIDCommonName %in% c("acquired immunodeficiency syndrome", "human immunodeficiency virus infectious disease") ~ "HIV",
    DOIDCommonName %in% c("sleep deprivation", "circadian rhythm disorders") ~ "SD",
    DOIDCommonName == "psoriasis" ~ "Psoriasis",
    DOIDCommonName == "mental depression" ~ "Depression",
    DOIDCommonName == "primary bacterial infectious disease" ~ "PBI",
    DOIDCommonName == "obesity" ~ "Obesity",
    DOIDCommonName == "no disease reported" ~ "Healthy",
    DOIDCommonName == "Kawasaki disease" ~ "KD",
    DOIDCommonName == "ischemic stroke" ~ "IS",
    DOIDCommonName == "Influenza Virus" ~ "Influenza",
    DOIDCommonName == "dental caries" ~ "Dental Caries",
    DOIDCommonName == "COVID-19" ~ "COVID-19",
    DOIDCommonName == "Chagas disease" ~ "CD",
    DOIDCommonName == "Alzheimer's disease" ~ "AD",
    TRUE ~ DOIDCommonName)) %>% 
  dplyr::mutate(Location = case_when(
    str_detect(UBERONBodyPartName, "skin") ~ "Skin",
    str_detect(UBERONBodyPartName, "blood") ~ "Blood",
    UBERONBodyPartName == "anal region" ~ "Anal Region",
    UBERONBodyPartName == "cerebrospinal fluid" ~ "CSF",
    UBERONBodyPartName == "brain" ~ "Brain",
    UBERONBodyPartName == "milk" ~ "Milk",
    UBERONBodyPartName == "oral cavity" ~ "Oral Cavity",
    UBERONBodyPartName == "feces" ~ "Stool",
    UBERONBodyPartName == "heart" ~ "Heart",
    UBERONBodyPartName == "vagina" ~ "Vagina",
    UBERONBodyPartName == "lung" ~ "Lung",
    UBERONBodyPartName == "nasal cavity" ~ "Skin",
    UBERONBodyPartName == "saliva" ~ "Saliva",
    UBERONBodyPartName == "stomach" ~ "Stomach",
    UBERONBodyPartName == "urine" ~ "Urine",
    TRUE ~ UBERONBodyPartName)) %>%
  group_by(Disease, Location) %>%
  distinct(Feature, .keep_all = TRUE) %>%
  summarise(UniqueFeatures = n())

microbial_location_split_human_baloon$Location <- factor(microbial_location_split_human_baloon$Location, 
                                                         levels = c("Brain", "CSF", "Oral Cavity", "Saliva", "Skin", "Blood", "Urine", 
                                                                    "Anal Region", "Vagina", "Milk", "Stool")) %>% fct_rev()

microbial_location_split_human_baloon$Disease <- factor(microbial_location_split_human_baloon$Disease, 
                                                        levels = c("Healthy", "COVID-19", "HIV", "PBI", "Psoriasis", 
                                                                   "SD", "Depression", "AD", "IS", "Dental Caries",
                                                                   "KD", "IBD", "T2D", "Obesity", "Not Specified")) 

Baloon_plot <- microbial_location_split_human_baloon %>% 
  ggballoonplot(x = "Location", y = "Disease", size = "UniqueFeatures", fill = "UniqueFeatures",
                title = "Unique Detection in Human Samples") + 
  scale_fill_viridis_c(direction = -1) +
  theme(axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7),
        plot.title = element_text(size = 8)) + coord_flip()

#ggsave(plot = Baloon_plot, filename = "Baloon_plot.svg", device = "svg", dpi = "retina", width = 4, height = 2.5)

# List 
microbial_location_split_human <- microbial_location_split_filter %>% 
  left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% 
  group_by(UBERONBodyPartName) %>% 
  #distinct(Feature, .keep_all = TRUE) %>%
  summarise(n = n()) %>% 
  dplyr::arrange(desc(n)) #figure 2e



# Distribution

microbial_location_split_brain <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "brain") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  dplyr::filter(!(phylum %in% c("Cyanobacteriota", "Deinococcota", "Myxococcota"))) %>%
  group_by(phylum) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>% 
  dplyr::filter(!(is.na(phylum))) %>%
  ggdonutchart(x = "count", fill = "phylum", title = "Brain") + 
  scale_fill_viridis_d() + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("Phyla"))

microbial_location_split_blood <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "blood plasma") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  dplyr::filter(!(phylum %in% c("Cyanobacteriota", "Deinococcota", "Myxococcota"))) %>%
  group_by(phylum) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>% 
  dplyr::filter(!(is.na(phylum))) %>%
  ggdonutchart(x = "count", fill = "phylum", title = "Blood") + 
  scale_fill_viridis_d() + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("Phyla"))

microbial_location_split_feces <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "feces") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  dplyr::filter(!(phylum %in% c("Cyanobacteriota", "Deinococcota", "Myxococcota"))) %>%
  group_by(phylum) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>% 
  dplyr::filter(!(is.na(phylum))) %>%
  ggdonutchart(x = "count", fill = "phylum", title = "Stool") + 
  scale_fill_viridis_d() + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("Phyla"))

microbial_location_split_urine <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "urine") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  dplyr::filter(!(phylum %in% c("Cyanobacteriota", "Deinococcota", "Myxococcota"))) %>%
  group_by(phylum) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>% 
  dplyr::filter(!(is.na(phylum))) %>%
  ggdonutchart(x = "count", fill = "phylum", title = "Urine") + 
  scale_fill_viridis_d() + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("Phyla"))

collapsed_microbial <- ggarrange(microbial_location_split_brain, microbial_location_split_blood, microbial_location_split_feces, microbial_location_split_urine,
                                 nrow = 1, common.legend = TRUE, legend = "right")

#ggsave(plot = collapsed_microbial, filename = "microbial_donut.svg", device = "svg", dpi = "retina", width = 8, height = 2)



microbial_location_split_brain_class <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "brain") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  distinct(Feature, .keep_all = TRUE) %>%
  group_by(NPC_Pathway) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>%
  replace(is.na(.), "Unknownw") %>%
  rbind(data.frame(NPC_Pathway = "Shikimates and Phenylpropanoids", count = 0)) %>%
  rbind(data.frame(NPC_Pathway = "Terpenoids", count = 0)) %>%
  ggdonutchart(x = "count", fill = "NPC_Pathway", title = "Brain") + 
  scale_fill_viridis_d(option="magma") + theme(legend.position = "right") + 
  theme(axis.text = element_blank()) + guides(fill = guide_legend("NPC_Pathway"))

microbial_location_split_blood_class <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "blood plasma") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  distinct(Feature, .keep_all = TRUE) %>%
  group_by(NPC_Pathway) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>%
  replace(is.na(.), "Unknownw") %>%
  ggdonutchart(x = "count", fill = "NPC_Pathway", title = "BLood") + 
  scale_fill_viridis_d(option="magma") + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("NPC_Pathway"))

microbial_location_split_feces_class <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "feces") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  distinct(Feature, .keep_all = TRUE) %>%
  group_by(NPC_Pathway) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>%
  replace(is.na(.), "Unknownw") %>%
  ggdonutchart(x = "count", fill = "NPC_Pathway", title = "Feces") + 
  scale_fill_viridis_d(option="magma") + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("NPC_Pathway"))

microbial_location_split_urine_class <- microbial_location_split_filter %>% left_join(redu_metadata_body, by = c("Sample" = "filename")) %>%
  dplyr::filter(NCBITaxonomy == "9606|Homo sapiens") %>% dplyr::filter(UBERONBodyPartName == "urine") %>% distinct(Feature) %>%
  left_join(unique_microbial_sirius, by = c("Feature" = "fileName")) %>%
  left_join(bacterial_features, by = c("Feature" = "fileName")) %>%
  distinct(Feature, .keep_all = TRUE) %>%
  group_by(NPC_Pathway) %>%
  summarise(count = n()) %>% arrange(desc(count)) %>%
  replace(is.na(.), "Unknownw") %>%
  rbind(data.frame(NPC_Pathway = "Shikimates and Phenylpropanoids", count = 0)) %>%
  ggdonutchart(x = "count", fill = "NPC_Pathway", title = "Urine") + 
  scale_fill_viridis_d(option="magma") + theme(legend.position="right") + theme(axis.text = element_blank()) + guides(fill = guide_legend("NPC_Pathway"))

collapsed_microbial_class <- ggarrange(microbial_location_split_brain_class, microbial_location_split_blood_class, 
                                       microbial_location_split_feces_class, microbial_location_split_urine_class,
                                       nrow = 1, common.legend = TRUE, legend = "right")

#ggsave(plot = collapsed_microbial_class, filename = "microbial_class_donut.svg", device = "svg", dpi = "retina", width = 8, height = 2)
