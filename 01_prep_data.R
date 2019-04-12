#################################################
## Allen Roberts
## March 2019
## Senegal Y Chromosome Analysis
#################################################

rm(list = ls())

set.seed(0)

library(tidyverse)
library(readxl)
library(lubridate)

## Combined data file
dat <- read_excel("../../yc/data/c.20180626_Senegal and Seattle Combined Juin 2018.xlsx")

## Number of FSW screened
nrow(dat) ## 350

## Yc results correction - discordant samples
dat$ychromom3[dat$h00IdNumber == 3056] <- "P"

## Site
dat$site <- substr(as.character(dat$h00IdNumber), start = 1, stop = 1)
dat$site[dat$site == "2"] <- "Pikine"
dat$site[dat$site == "3"] <- "Mbao"
dat$site[dat$site == "4"] <- "Rifisque"
dat$site[dat$site == "5"] <- "Diamniado"

## Age
dat$age <- dat$F00_scr_q01AgeDuSujet
dat$age_cat <- cut(dat$age, c(17, 29, 39, 49, Inf),  labels = c("<30", "30-39", "40-49", "50+"))

## Registered FSW
dat$registered <- factor(dat$F00_scr_q10aTSEnregistree, levels = c(2, 1), labels = c("No", "Yes"))

## Ethnic Group
dat$ethnic_group <- factor(dat$F00_scr_q03Ethnie, levels = c(1:8, 99), labels = c("Wolof", "Poulaar", "Serere", "Mandingue/Bambara", "Diola", "Soninke", "Maure", "Manjak", "Other")) ## Note that Geoff's tables has alternative names and spellings for these groups
dat$ethnic_group[dat$F00_scr_q03EthnieAutre == "MANDJAQUE"] <- "Manjak"
dat$ethnic_group_cat <- factor(c("Wolof", "Poulaar", "Serere", "Mandingue", "Other", "Other", "Other", "Other", "Other")[dat$ethnic_group], levels = c("Wolof", "Poulaar", "Serere", "Mandingue", "Other"))

## Education
dat$edu <- ifelse(is.na(dat$F00_scr_q09Scolarite), 0, dat$F00_scr_q09Scolarite)
dat$edu <- factor(dat$edu, levels = c(0, 1, 2), labels = c("None", "Primary", "Secondary"))

sex_question_stems <- c("NbreClieDernSemaine", "UtiliseUnPresevatif", "FreqUtilPreservatif", "UtilisePvatifDepuis", "PartenaireSexPrinc", "NbrePartMoisDernier", "FreqUtilervatifPart", "UtilisePartenDepuis", "PreservaPartenaires")
sex_dat <- dat[, c("h00IdNumber", "site", grep(paste(sex_question_stems, collapse = "|"), names(dat), value = TRUE))]

## Tidy up
sex_dat <- sex_dat %>% 
  rename(id = h00IdNumber) %>%
  gather(key, value, -id, -site) %>%
  mutate(visit = gsub(".*_(.*)\\_.*","\\1", key, perl=T)) %>% 
  mutate(variable = gsub(".*(\\d{2})", "", key)) %>%
  mutate(visit = replace(visit, visit == "j0", "scr")) %>% ## Note that this is to facilitate merging with Yc
  select(-key) %>%
  spread(key = variable, value = value)

## Clean up specific responses
sex_dat$NbreClieDernSemaine[sex_dat$NbreClieDernSemaine %in% c("NE CONNAIT PAS", "REFUS")] <- NA
sex_dat$NbreClieDernSemaine[sex_dat$NbreClieDernSemaine %in% c("O", "OO")] <- "0"
sex_dat$NbreClieDernSemaine <- as.integer(sex_dat$NbreClieDernSemaine)

## Ychromosome 
grep("ychrom", names(dat), value = TRUE)
yc_dat <- dat[, c("h00IdNumber", "site", grep("ychrom", names(dat), value = TRUE))]

## Tidy up
yc_dat <- yc_dat %>% 
  rename(id = h00IdNumber) %>%
  gather(key, value, -id, -site) %>%
  filter(!is.na(value)) %>%
  mutate(type = ifelse(grepl("date", key), "date", "value")) %>%
  mutate(visit = ifelse(grepl("date", key), substr(key, start = nchar("ychromo") + 1, stop = nchar(key) - nchar("date")), substr(key, start = nchar("ychromo") + 1, stop = nchar(key)))) %>%
  select(-key) %>%
  spread(key = type, value = value)

## Fix mislabeled ID 
yc_dat$id[yc_dat$id == 4059] <- 4052

## Remove uncertain value
yc_dat <- yc_dat[yc_dat$value != "U", ]

## Merge sexual behavior and y-chromosome results
sex_dat <- left_join(sex_dat, yc_dat, by = c("id", "visit", "site"))

## New variables
sex_dat$yc_result <- ifelse(sex_dat$value == "P", 1, 0)

## Visit as a factor variable
yc_dat$visit_factor <- factor(yc_dat$visit, levels = c("scr", "m1", "m3", "m6", "m9", "m12"))
sex_dat$visit_factor <- factor(sex_dat$visit, levels = c("scr", "m1", "m3", "m6", "m9", "m12"))

## Create condom binaries
## Sex in prior month
sex_dat$always_condom_client <- ifelse(sex_dat$FreqUtilPreservatif == "toujours", 1, 0)
sex_dat$not_always_condom_client <- 1 - sex_dat$always_condom_client
sex_dat$always_condom_main <- ifelse(sex_dat$FreqUtilervatifPart == "toujours", 1, 0)
sex_dat$not_always_condom_main <- 1 - sex_dat$always_condom_main
sex_dat$always_condom_both <- ifelse(sex_dat$always_condom_client == 1 & sex_dat$always_condom_main == 1, 1, 0)
sex_dat$not_always_condom_both <- 1 - sex_dat$always_condom_both

## Most recent sex
sex_dat$used_condom_recent_client <- ifelse(sex_dat$UtiliseUnPresevatif == "o", 1, 0)
sex_dat$not_used_condom_recent_client <- 1 - sex_dat$used_condom_recent_client
sex_dat$used_condom_recent_main <- ifelse(sex_dat$PreservaPartenaires == "o", 1, 0)
sex_dat$not_used_condom_recent_main <- 1 - sex_dat$used_condom_recent_main
sex_dat$used_condom_recent_both <- ifelse(sex_dat$used_condom_recent_client == 1 & sex_dat$used_condom_recent_main == 1, 1, 0)
sex_dat$not_used_condom_recent_both <- 1 - sex_dat$used_condom_recent_both

## Reported sexual partners
## Reports sex with a main partner in the last month
sex_dat$reports_main_partner <- ifelse(sex_dat$NbrePartMoisDernier > 0, 1, 0)
sex_dat$reports_main_partner[sex_dat$PartenaireSexPrinc == "n"] <- 0

## Reports sex with a client in the last week
sex_dat$reports_client <- ifelse(sex_dat$NbreClieDernSemaine > 0, 1, 0)

## Reports either
sex_dat$reports_either <- ifelse(sex_dat$reports_client == 1 | sex_dat$reports_main_partner == 1, 1, 0)

## STI detection
sti_question_stems <- c("RaisonPPrelevement", "EstDispost", "RaisonPat", "Resultatst", "TypeTest", "Specificst", "atestSyphilis", "aResultatRPR", "bResultatTPHA")
sti_dat <- dat[, c("h00IdNumber", "site", grep(paste(sti_question_stems, collapse = "|"), names(dat), value = TRUE))]

## Tidy up
sti_dat <- sti_dat %>% 
  rename(id = h00IdNumber) %>%
  gather(key, value, -id, -site) %>%
  mutate(visit = gsub(".*_(.*)\\_.*","\\1", key, perl=T)) %>% 
  mutate(variable = gsub(".*(\\d{2})", "", key)) %>%
  mutate(visit = replace(visit, visit == "j0", "scr")) %>% ## Note that this is to facilitate merging with Yc
  select(-key) %>%
  spread(key = variable, value = value)

sti_dat$visit_factor <- factor(sti_dat$visit, levels = c("scr", "m1", "m3", "m6", "m9", "m12"))

## Drop individuals who weren't enrolled - based on j0 visit date
enrolled_ids <- unique(dat$h00IdNumber[!is.na(dat$F00_j0_h00VisitDate)])
dat <- dat[dat$h00IdNumber %in% enrolled_ids, ]
sex_dat <- sex_dat[sex_dat$id %in% enrolled_ids, ]
sti_dat <- sti_dat[sti_dat$id %in% enrolled_ids, ]
yc_dat <- yc_dat[yc_dat$id %in% enrolled_ids, ]

