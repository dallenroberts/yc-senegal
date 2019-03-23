#################################################
## Allen Roberts
## January 2019
## Senegal Y Chromosome Analysis
#################################################

rm(list = ls())

set.seed(0)

library(tidyverse)
library(readxl)
library(lubridate)
library(binom)
library(gridExtra)
library(cowplot)
library(geepack)
library(sandwich)
library(lmtest)

## Combined data file
dat <- read_excel("../data/c.20180626_Senegal and Seattle Combined Juin 2018.xlsx")

## ID: seems to be h00IdNumber: _ - _ _ _
length(unique(dat$h00IdNumber))
nrow(dat)
grep("h00Id", names(dat), value = TRUE) ## Why are there different h00Id numbers? 

## Site
dat$site <- substr(as.character(dat$h00IdNumber), start = 1, stop = 1)
dat$site[dat$site == "2"] <- "Pikine"
dat$site[dat$site == "3"] <- "Mbao"
dat$site[dat$site == "4"] <- "Rifisque"
dat$site[dat$site == "5"] <- "Diamniado"

## Create "presence" indicator based on having a valid visit date - seems like there is a presence indicator in the  data dictionary but I can't find it in this dataset
## Note that read_excel has a problem parsing some of the dates and gives the excel date verison (number of days since 1899-12-30) as a character for some of these. Could mess with this if I wanted to later:
## https://github.com/tidyverse/readxl/issues/118

grep("VisitDate", names(dat), value = TRUE) ## There's a "F21_m12_h00VisitDate" that I'm not sure what it corresponds to

for(visit in c("scr", "j0", "j7", "m1", "m3", "m6", "m9", "m12")) {
  
  ## Get the right variable
  if(visit != "m12") {
    date_var <- names(dat)[grepl("VisitDate", names(dat)) & grepl(visit, names(dat))]
  } else {
    date_var <- names(dat)[grepl("VisitDate", names(dat)) & grepl(visit, names(dat)) & grepl("F01", names(dat))]
  }

  dat[, paste0("presence_", visit)] <- ifelse(is.na(dat[, date_var]), 0, 1)
}

## Missingness Table
presence_dat <- dat[, c("h00IdNumber", "site", grep("presence", names(dat), value = TRUE))]

presence_dat <- presence_dat %>%
  rename(id = h00IdNumber) %>%
  gather(key, value, -id, -site) %>% ## Craetes a warning about attributes that I think doesn't matter
  mutate(visit = gsub("presence_", "", key), variable = "presence") %>%
  select(-key) %>%
  spread(key = variable, value = value)
           
presence_dat %>%
  group_by(visit) %>%
  summarise(present = sum(presence), absent = sum(1 - presence), pct = mean(presence))

## Sexual behavior questions
  ## NbreClieDernSemaine: Number of clients in the last week
  ## UtiliseUnPresevatif: Did you use a condom the last time that you had sex (vaginal or anal) with a client?
  ## FreqUtilPreservatif: Frequecy of condom use with clients in the last month
  ## UtilisePvatifDepuis: For those who report always using a condom with clients, for how long have they been consistently using one?
  ## PartenaireSexPrinc: Has a main partner in the last six months 
  ## NbrePartMoisDernier: Number of main partners last month
  ## PreservaPartenaires: Did you use a condom the last time you had sex with a main partner?
  ## FreqUtilervatifPart: Frequecy of condom use with main partners in the last month
  ## UtilisePartenDepuis: For those who report always using a condom with main partners, for how long have they been consistently using one?


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

## Clean up
sex_dat$NbreClieDernSemaine[sex_dat$NbreClieDernSemaine %in% c("NE CONNAIT PAS", "REFUS")] <- NA
sex_dat$NbreClieDernSemaine[sex_dat$NbreClieDernSemaine %in% c("O", "OO")] <- 0

## View(sex_dat) ## Some clients seem to be missing all of the sexual behavior data - is there a trigger somewhere?

## Create condom binaries
sex_dat$always_condom_client <- ifelse(sex_dat$FreqUtilPreservatif == "toujours", 1, 0)
sex_dat$always_condom_main <- ifelse(sex_dat$FreqUtilervatifPart == "toujours", 1, 0)

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

## Need to construct date-type variable - note that the collection dates don't seem to be systematically dmy or mdy. Should cross-reference with visit date.

#################################################
## Geoff's tables
#################################################
## 165 samples from 131 FSW tested for Yc; 164 gave valid results
length(unique(yc_dat$id)) ## I get 132 different women, not 131
geoff_unique_pids <- read_excel("../data/geoff_unique_patient_ids.xlsx")
unique(yc_dat$id)[!unique(yc_dat$id) %in% geoff_unique_pids$id] ## 4008 - this is the undefined result, which Geoff must have removed

nrow(yc_dat) ## I get 164 results, not 165
table(yc_dat$value) ## 1 U - I assume this is the undefined one
yc_dat[yc_dat$value == "U",] ## This is the ID highlighted in yellow in the spreadsheet

geoff_ids <- read_excel("../data/geoff_ids.xlsx")
geoff_ids <- geoff_ids %>%
  mutate(visit = ifelse(visit == "screen", "scr", paste0("m", visit)))
geoff_ids[duplicated(geoff_ids[c("id", "visit")]),] ## 2 swab run for this patient - one was positive and one was negative. Geoff thinks that the positive swab is correct. Overall dataset then has the wrong result.

## Make those corrections here
yc_dat$value[yc_dat$id == 3056] <- "P"
dat$ychromom3[dat$h00IdNumber == 3056] <- "P"

## Remove uncertain value
yc_dat <- yc_dat[yc_dat$value != "U", ]

############################################################################################
## Table 1: Demographics and baseline characteristics of FSW that had Yc testing
############################################################################################
## Baseline characterstics dataset - note that this does not include woman for whom swab result was undefined
base_dat <- dat[dat$h00IdNumber %in% unique(yc_dat$id), ]

## Number of FSW tested for yc
nrow(base_dat) ## 131

## Age (median, IQR)
length(which(is.na(base_dat$F00_scr_q01AgeDuSujet)))
quantile(base_dat$F00_scr_q01AgeDuSujet, probs = c(0.25, 0.5, 0.75))
IQR(base_dat$F00_scr_q01AgeDuSujet)

## Born in Senegal (%), and birth locations of individuals born outside of Senegal
mean(base_dat$F00_scr_q04Nationalite == 1)
base_dat$F00_scr_q04NationaliteAutre[base_dat$F00_scr_q04Nationalite == 2]

## Registered vs. unregistered FSW (%)
mean(base_dat$F00_scr_q10aTSEnregistree == 1)

## Site (N, %)
table(base_dat$site)
table(base_dat$site)/nrow(base_dat)

## Ethnic Group, and number missing
table(base_dat$F00_scr_q03Ethnie, base_dat$F00_scr_q03EthnieAutre, exclude = NULL)

base_dat$ethnic_group <- factor(base_dat$F00_scr_q03Ethnie, levels = c(1:8, 99), labels = c("Wolof", "Poulaar", "Serere", "Mandingue/Bambara", "Diola", "Soninke", "Maure", "Manjak", "Other")) ## Note that Geoff's tables has alternative names and spellings for these groups
base_dat$ethnic_group[base_dat$F00_scr_q03EthnieAutre == "MANDJAQUE"] <- "Manjak"
table(base_dat$ethnic_group)
table(base_dat$ethnic_group)/length(which(!is.na(base_dat$ethnic_group)))

## Self-reported number of sex partners in the prior week (Median, range), and number missing - I think this is really number of clients, not number of partners, right?
quantile(as.numeric(base_dat$F14_j0_q03NbreClieDernSemaine), probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
min(as.numeric(base_dat$F14_j0_q03NbreClieDernSemaine), na.rm = TRUE)
max(as.numeric(base_dat$F14_j0_q03NbreClieDernSemaine), na.rm = TRUE)
mean(as.numeric(base_dat$F14_j0_q03NbreClieDernSemaine), na.rm = TRUE)
length(which(is.na(base_dat$F14_j0_q03NbreClieDernSemaine)))

base_dat$num_clients_cat <- cut(as.integer(base_dat$F14_j0_q03NbreClieDernSemaine), breaks = c(-Inf, 0, 2, Inf ), labels = c("0", "1-2", "3+"))
base_dat %>%
  group_by(num_clients_cat) %>%
  summarise(n = n()) %>%
  filter(!is.na(num_clients_cat)) %>%
  mutate(pct = n/sum(n))

## Self-reported had a "Main Partner" - this is in the last six months
table(base_dat$F14_j0_q11PartenaireSexPrinc, exclude = NULL)
sum(base_dat$F14_j0_q11PartenaireSexPrinc == "o", na.rm = TRUE)
mean(base_dat$F14_j0_q11PartenaireSexPrinc == "o", na.rm = TRUE)

# sum(base_dat$F14_j0_q13NbrePartMoisDernier > 0, na.rm = TRUE)
# mean(base_dat$F14_j0_q13NbrePartMoisDernier > 0, na.rm = TRUE)
# 
# table(base_dat$F14_j0_q11PartenaireSexPrinc, base_dat$F14_j0_q12NbrePartenaireSex, exclude = NULL, deparse.level = 2)
# 
# table(base_dat$F14_j0_q12NbrePartenaireSex, exclude = NULL) ## Number of main partners - this isn't always 1. Missings could be people who report no main partner
# table(base_dat$F14_j0_q12NbrePartenaireSex[base_dat$F14_j0_q11PartenaireSexPrinc == "o"], exclude = NULL)
# median(base_dat$F14_j0_q12NbrePartenaireSex[base_dat$F14_j0_q11PartenaireSexPrinc == "o"], na.rm = TRUE)
# min(base_dat$F14_j0_q12NbrePartenaireSex[base_dat$F14_j0_q11PartenaireSexPrinc == "o"], na.rm = TRUE)
# max(base_dat$F14_j0_q12NbrePartenaireSex[base_dat$F14_j0_q11PartenaireSexPrinc == "o"], na.rm = TRUE)

## Self-reported condom use with clients in the last month for vaginal or anal sex, and number missing
table(base_dat$F14_j0_q05FreqUtilPreservatif, exclude = NULL)
table(base_dat$F14_j0_q05FreqUtilPreservatif, exclude = NULL)/length(which(!is.na(base_dat$F14_j0_q05FreqUtilPreservatif)))
table(base_dat$F14_j0_q05FreqUtilPreservatif)/length(which(!is.na(base_dat$F14_j0_q05FreqUtilPreservatif)))

## Self-reported condom use with main partner in the last month for vaginal or anal sex, and number missing
table(base_dat$F14_j0_q15FreqUtilervatifPart[base_dat$F14_j0_q11PartenaireSexPrinc == "o"], exclude = NULL)
table(base_dat$F14_j0_q15FreqUtilervatifPart[base_dat$F14_j0_q11PartenaireSexPrinc == "o"], exclude = NULL)/length(which(!is.na(base_dat$F14_j0_q15FreqUtilervatifPart[base_dat$F14_j0_q11PartenaireSexPrinc == "o"])))
table(base_dat$F14_j0_q15FreqUtilervatifPart[base_dat$F14_j0_q11PartenaireSexPrinc == "o"])/length(which(!is.na(base_dat$F14_j0_q15FreqUtilervatifPart[base_dat$F14_j0_q11PartenaireSexPrinc == "o"])))


############################################################################################
## Table 2: Detection of Yc in FSW vaginal swabs
############################################################################################
## FSW with detectable Yc
yc_dat %>%
  group_by(id) %>%
  summarise(any_pos = any(value == "P")) %>%
  ungroup() %>%
  summarise(num_fsw = n_distinct(id), num_fsw_pos = sum(any_pos), pct_fsw_pos = mean(any_pos)) %>%
  mutate(se = sqrt(pct_fsw_pos * (1-pct_fsw_pos)/num_fsw)) %>%
  mutate(lower = pct_fsw_pos - 1.96*se, upper = pct_fsw_pos + 1.96*se) 

## Wilson interval
binom.confint(x = 32, n = 131, methods = "wilson")

## Vaginal swabs with detectable Yc
yc_dat %>%
  summarise(num_swabs = n(), num_swabs_pos = sum(value == "P"), pct_swabs_pos = mean(value == "P")) %>%
  mutate(se = sqrt(pct_swabs_pos * (1-pct_swabs_pos)/num_swabs)) %>%
  mutate(lower = pct_swabs_pos - 1.96*se, upper = pct_swabs_pos + 1.96*se)

## Wilson interval
binom.confint(x = 35, n = 163, methods = "wilson")

## Swab positivity by visit month
yc_dat %>%
  mutate(visit_group = ifelse(visit %in% c("m1", "m3"), "m1-3", visit)) %>%
  mutate(visit_group = factor(visit_group, levels = c("scr", "m1-3", "m6", "m9", "m12"))) %>%
  group_by(visit_group) %>%
  summarise(num_swabs = n(), num_swabs_pos = sum(value == "P"), pct_swabs_pos = mean(value == "P")) %>%
  mutate(se = sqrt(pct_swabs_pos * (1-pct_swabs_pos)/num_swabs)) %>%
  mutate(lower = pct_swabs_pos - 1.96*se, upper = pct_swabs_pos + 1.96*se)

## No significant difference between Baseline and PrEP use for detection of Y-chromosomes (p>0.05 by Fisher Exact Test)
cont_table_pre_post <- as.matrix(yc_dat %>%
  mutate(visit_group = ifelse(visit == "scr", "pre-PrEP", "post-PrEP")) %>%
  mutate(visit_group = factor(visit_group, levels = c("pre-PrEP", "post-PrEP"))) %>%
  group_by(visit_group) %>%
  summarise(num_swabs_pos = sum(value == "P"), num_swabs_neg = sum(value == "N")) %>%
  select(num_swabs_pos, num_swabs_neg))
row.names(cont_table_pre_post) <-c("Pre-PrEP", "Post-PrEP")
cont_table_pre_post
fisher.test(cont_table_pre_post)

## Sites:  Diamniadio (n=40),  Mbao (n=43), Pikine (n=41), Rufisque (n=40) 
yc_dat %>%
  group_by(site) %>%
  summarise(num = n())

############################################################################################
## Table 3: Longitudinal detection of Yc in FSW vaginal swabs
############################################################################################
## FSW with longitudinal sampling
yc_dat %>%
  add_count(id) %>%
  filter(n > 1) %>%
  group_by(id) %>%
  summarise(always_detect = all(value == "P"), never_detect = all(value == "N")) %>%
  summarise(n_fsw = n(), n_always_detect = sum(always_detect), pct_always_detect = mean(always_detect), n_never_detect = sum(never_detect), pct_never_detect = mean(never_detect))
  
## Need to verify that n_fsw is different by one because of the discordant result in Geoff's data - hard to do because not sure how he calculated it

## Different initial result than subsequent result - use visit numbers to establish temporality
yc_dat %>%
  add_count(id) %>%
  filter(n > 1) %>%
  mutate(visit_num = ifelse(visit == "scr", 0, as.numeric(gsub("m", "", visit)))) %>%
  group_by(id) %>%
  summarise(initial_pos = value[visit_num == min(visit_num)] == "P", fu_pos = any(value[visit_num != min(visit_num)] == "P"), all_same = (all(value == "P") | all(value == "N"))) %>%
  mutate(n_fsw = n()) %>%
  filter(all_same == FALSE) %>%
  ungroup() %>%
  summarise(n_init_neg_later_pos = sum(initial_pos == FALSE), pct_init_neg_later_pos = sum(initial_pos == FALSE)/unique(n_fsw), n_init_pos_later_neg = sum(initial_pos == TRUE), pct_init_pos_later_neg = sum(initial_pos == TRUE)/unique(n_fsw))


############################################################################################
## Table 4: Self-reported condom use in the last month prior to vaginal swab in those positive for Yc
############################################################################################
pos_swabs <- yc_dat[yc_dat$value == "P",]
nrow(pos_swabs) 
length(unique(pos_swabs$id))
## 35 swabs from 32 women

neg_swabs <- yc_dat[yc_dat$value == "N",]
nrow(neg_swabs)
length(unique(neg_swabs$id))
## 128 swabs from 108 women

# Merge together with sex dat
yc_dat <- left_join(yc_dat, sex_dat, by = c("id", "visit", "site"))

## Condom usage with clients by detection status
yc_dat %>%
  group_by(value, FreqUtilPreservatif) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n)) %>%
  ungroup() %>%
  summarise(sum(n))

## Missings removed
yc_dat %>%
  group_by(value, FreqUtilPreservatif) %>%
  summarise(n = n()) %>%
  filter(!is.na(FreqUtilPreservatif)) %>%
  mutate(pct = n/sum(n))
  
## Condom usage with main partner by detection status
yc_dat %>%
  filter(NbrePartMoisDernier > 0 | is.na(NbrePartMoisDernier)) %>%
  group_by(value, FreqUtilervatifPart) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>% ## Drop missings
  filter(NbrePartMoisDernier > 0 | is.na(NbrePartMoisDernier)) %>%
  group_by(value, FreqUtilervatifPart) %>%
  summarise(n = n()) %>%
  filter(!is.na(FreqUtilervatifPart)) %>%
  mutate(pct = n/sum(n))

yc_dat %>% 
  filter(NbrePartMoisDernier == 0) %>% 
  group_by(value, FreqUtilervatifPart) %>%
  summarise(n = n())

## Yc positivity by condom status - use collapsed categories
yc_dat %>%
  group_by(always_condom_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbreClieDernSemaine > 0) %>%
  group_by(always_condom_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))


yc_dat %>%
  group_by(always_condom_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(PartenaireSexPrinc == "o") %>%
  group_by(always_condom_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat$always_condom_both <- ifelse(yc_dat$always_condom_client == 1 & yc_dat$always_condom_main == 1, 1, 0)
yc_dat %>%
  group_by(always_condom_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbrePartMoisDernier > 0 & NbreClieDernSemaine > 0) %>%
  group_by(always_condom_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))


## Yc positivity by condom use at last sex act
## With clients
yc_dat %>%
  group_by(UtiliseUnPresevatif, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbreClieDernSemaine > 0) %>%
  group_by(UtiliseUnPresevatif, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## With main partners
yc_dat %>%
  group_by(PreservaPartenaires, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(PreservaPartenaires, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## Both
yc_dat$used_condom_last_both <- ifelse(yc_dat$PreservaPartenaires == "o" & yc_dat$UtiliseUnPresevatif == "o", 1, 0)
yc_dat %>%
  group_by(used_condom_last_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(used_condom_last_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  filter(NbrePartMoisDernier > 0) %>%
  filter(NbreClieDernSemaine > 0) %>%
  group_by(used_condom_last_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

############################################################################################
## Sexual behavior explorations
############################################################################################
## Detection among FSW reporting constant condom use vs. never condom use
yc_dat <- yc_dat %>%
  mutate(always_condom = ifelse(FreqUtilervatifPart == "toujours" & FreqUtilPreservatif == "toujours", 1, 0)) %>%
  mutate(never_condom = ifelse(FreqUtilervatifPart == "jamais" & FreqUtilPreservatif == "jamais", 1, 0))

## Always
yc_dat %>% 
  group_by(always_condom) %>%
  summarise(n = n(), num_pos = sum(value == "P"), num_neg = sum(value == "N")) %>%
  mutate(pct_pos = num_pos/n)

## Never
yc_dat %>% 
  group_by(never_condom) %>%
  summarise(n = n(), num_pos = sum(value == "P"), num_neg = sum(value == "N")) %>%
  mutate(pct_pos = num_pos/n)

## Partnership reporting
yc_dat <- yc_dat %>%
  mutate(reports_main_partner = ifelse(NbrePartMoisDernier > 0, 1, 0)) %>%
  mutate(reports_main_partner = replace(reports_main_partner, PartenaireSexPrinc == "n", 0)) %>%
  mutate(reports_client = ifelse(NbreClieDernSemaine > 0, 1, 0)) %>%
  mutate(reports_either = ifelse(reports_client == 1 | reports_main_partner == 1, 1, 0))

## How many of those sampled reported having a main partner in the last month? How many had detectable Yc?
yc_dat %>%
  group_by(reports_main_partner) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  group_by(reports_main_partner, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## How many of those sampled reported having a client in the last 7 days?
yc_dat %>%
  group_by(reports_client) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  group_by(reports_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## How many of those sampled reported either having a main partner in the last month or having a client in the last 7 days?
yc_dat %>%
  group_by(reports_either) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

yc_dat %>%
  group_by(reports_either, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))


## How has condom use changed over (study) time?
## Create visit groupings
sex_dat$visit_factor <- factor(sex_dat$visit, levels = c("scr", "m1", "m3", "m6", "m9", "m12"))
sex_dat$visit_group <- factor(c("scr", "m1-3", "m1-3", "m6", "m9", "m12")[sex_dat$visit_factor], levels = c("scr", "m1-3", "m6", "m9", "m12"))

## Proportion NOT always using a condom 
condom_dat <- sex_dat %>%
  select(id, visit_factor, always_condom_client, always_condom_main) %>%
  gather(key = "measure", value = "value", -id, -visit_factor) %>%
  mutate(measure = fct_recode(measure, "Main Partner" = "always_condom_main", "Clients" = "always_condom_client")) %>%
  group_by(measure, visit_factor) %>%
  summarise(pct = mean(1-value, na.rm = TRUE), num_inconsistent = sum(1-value, na.rm = TRUE), n = sum(!is.na(value))) %>%
  mutate(visit_num = c(0, 1, 3, 6, 9, 12))

condom_dat[, c("lower", "upper")] <- binom.confint(condom_dat$num_inconsistent, condom_dat$n, methods = "wilson")[, c("lower", "upper")]


condom_use_time_line_plot <- ggplot(condom_dat, aes(x = visit_num, y = pct * 100, group = measure)) +
  geom_line(aes(color = measure)) +
  # geom_text(aes(label = paste0(floor(pct*100), "%")), vjust = -1, size = 2.5) +
  # geom_text(aes(label = paste0("n=", n)), vjust = 1.5, size = 2.5) +
  labs(x = "Months after enrollment", y = "% reporting \n inconsistent condom use") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  theme_classic() +
  ## ggtitle("Proportion reporting inconsistent condom use among all FSW") +
  theme(plot.title = element_text(size = 8),
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "bottom",
        axis.text = element_text(size = 6), 
        axis.title=element_text(size=7), 
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))
  
# jpeg(file = "../condom_use_time.jpg", height = 4, width = 5, unit = 'in', res = 500)
# print(condom_use_time_plot)
# dev.off()

condom_use_time_point_plot <- ggplot(condom_dat, aes(x = visit_num, y = pct * 100, group = measure)) +
  geom_point(aes(color = measure)) +
  geom_errorbar(aes(ymin = 100*lower, ymax = 100*upper, colour = measure), width = 0.25) +
  # geom_text(aes(label = paste0(floor(pct*100), "%")), vjust = -1, size = 2.5) +
  # geom_text(aes(label = paste0("n=", n)), vjust = 1.5, size = 2.5) +
  labs(x = "Months after enrollment", y = "% reporting \n inconsistent condom use") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  theme_classic() +
  ## ggtitle("Proportion reporting inconsistent condom use among all FSW") +
  theme(plot.title = element_text(size = 8),
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "bottom",
        axis.text = element_text(size = 6), 
              axis.title=element_text(size=7),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

## Among women included in Yc sample
sex_dat_yc <- sex_dat[sex_dat$id %in% unique(yc_dat$id), ]

condom_dat_yc <- sex_dat_yc %>%
  select(id, visit_factor, always_condom_client, always_condom_main) %>%
  gather(key = "measure", value = "value", -id, -visit_factor) %>%
  mutate(measure = fct_recode(measure, "Clients" = "always_condom_client", "Main Partner" = "always_condom_main")) %>%
  group_by(measure, visit_factor) %>%
  summarise(pct = mean(value, na.rm = TRUE), n = sum(value, na.rm = TRUE))

condom_use_time_plot_yc <- ggplot(condom_dat_yc, aes(x = visit_factor, y = pct * 100, group = measure)) +
  geom_line(aes(color = measure)) +
  geom_text(aes(label = paste0(floor(pct*100), "%")), vjust = -1, size = 2.5) +
  geom_text(aes(label = paste0("n=", n)), vjust = 1.5, size = 2.5) +
  labs(x = "Months after enrollment", y = "% reporting always using condoms") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme_classic() +
  ggtitle("Proportion reporting always condom use among FSW included in Yc sample") +
  theme(plot.title = element_text(size = 8))


jpeg(file = "../condom_use_time_yc.jpg", height = 4, width = 5, unit = 'in', res = 500)
print(condom_use_time_plot_yc)
dev.off()

## By site
condom_dat_site <- sex_dat %>%
  select(id, visit_factor, site, always_condom_client, always_condom_main) %>%
  gather(key = "measure", value = "value", -id, -visit_factor, -site) %>%
  mutate(measure = fct_recode(measure, "Clients" = "always_condom_client", "Main Partner" = "always_condom_main")) %>%
  group_by(measure, site, visit_factor) %>%
  summarise(pct = mean(value, na.rm = TRUE), n = sum(value, na.rm = TRUE))

condom_use_time_site_plot <- ggplot(condom_dat_site, aes(x = visit_factor, y = pct * 100, group = site)) +
  geom_line(aes(color = site)) +
  ## geom_text(aes(label = paste0("n=", n)), vjust = -1, size = 2.5) +
  labs(x = "Months after enrollment", y = "% reporting always using condoms") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  facet_wrap(~measure) +
  theme_bw() +
  ggtitle("Proportion reporting always condom use among all FSW") +
  theme(plot.title = element_text(size = 12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(colour='black', fill=NA),
        panel.border = element_rect(fill = NA, color = "black"))

jpeg(file = "../condom_use_time_site.jpg", height = 5, width = 7, unit = 'in', res = 500)
print(condom_use_time_site_plot)
dev.off()


## Number of clients over time
clients_dat <- sex_dat %>%
  select(id, visit_factor, NbreClieDernSemaine) %>%
  mutate(number_clients = as.numeric(NbreClieDernSemaine)) %>%
  group_by(visit_factor) %>%
  summarise(mean = mean(number_clients, na.rm = TRUE), median = median(number_clients, na.rm = TRUE), sd = sd(number_clients, na.rm = TRUE), n = n()) %>%
  # gather(key = "measure", value = "value", -visit_factor) %>%
  # filter(measure != "n") %>%
  mutate(visit_num =c(0, 1, 3, 6, 9, 12)) %>%
  mutate(lower = mean - (1.96*sd/sqrt(n)), upper = mean + (1.96*sd/sqrt(n)))



clients_time_line_plot <- ggplot(data = clients_dat, aes(x = visit_num, y = mean)) +
  geom_line() +
  theme_classic() +
  labs(x = "Months after enrollment", y = "Reported number of \n clients in past week") +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, by = 1)) +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3))  +
  theme(axis.text = element_text(size = 6), 
        axis.title=element_text(size=7))

clients_time_point_plot <- ggplot(data = clients_dat, aes(x = visit_num, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  theme_classic() +
  labs(x = "Months after enrollment", y = "Reported number of \n clients in past week") +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  theme(axis.text = element_text(size = 6), 
        axis.title=element_text(size=7))
              
# jpeg(file = "../clients_time.jpg", height = 4, width = 5, unit = 'in', res = 500)
# print(clients_time_plot)
# dev.off()

clients_site_dat <- sex_dat %>%
  select(id, visit_factor, site, NbreClieDernSemaine) %>%
  mutate(number_clients = as.numeric(NbreClieDernSemaine)) %>%
  group_by(site, visit_factor) %>%
  summarise(mean = mean(number_clients, na.rm = TRUE), median = median(number_clients, na.rm = TRUE), n = n()) %>%
  gather(key = "measure", value = "value", -visit_factor, -site) %>%
  filter(measure != "n")

clients_time_site_plot <- ggplot(data = clients_site_dat, aes(x = visit_factor, y = value, group = site)) +
  geom_line(aes(colour = site)) +
  theme_classic() +
  facet_wrap(~measure) +
  labs(x = "Months after enrollment", y = "Reported number of clients in past week") +
  scale_y_continuous(limits = c(0, 8), breaks = seq(0, 8, by = 2))

jpeg(file = "../clients_time_site.jpg", height = 5, width = 7, unit = 'in', res = 500)
print(clients_time_site_plot)
dev.off()


## Proportion with Yc detected over time
yc_time_dat <- yc_dat %>%
  mutate(visit = factor(visit, levels = c("scr", "m1", "m3", "m6", "m9", "m12"), labels = c("0", "1", "3", "6", "9", "12"))) %>%
  group_by(visit) %>%
  summarise(num_swabs = n(), num_swabs_pos = sum(value == "P"), pct_swabs_pos = mean(value == "P")) %>%
  mutate(visit_num = as.numeric(as.character(visit)))

yc_time_dat[, c("lower", "upper")] <- binom.confint(yc_time_dat$num_swabs_pos, yc_time_dat$num_swabs, methods = "wilson")[, c("lower", "upper")]

yc_time_line_plot <- ggplot(data = yc_time_dat, aes(x = visit_num, y = 100*pct_swabs_pos)) +
  geom_line() +
  theme_classic() +
  labs(x = "Months after enrollment", y = "% of swabs with \n Yc detected") +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme(axis.text = element_text(size = 6), 
        axis.title=element_text(size=7))

yc_time_point_plot <- ggplot(data = yc_time_dat, aes(x = visit_num, y = 100*pct_swabs_pos)) +
  geom_point() +
  geom_errorbar(aes(ymin = 100*lower, ymax = 100*upper), width = 0.2) +
  theme_classic() +
  labs(x = "Months after enrollment", y = "% of swabs with \n Yc detected") +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme(axis.text = element_text(size = 6), 
        axis.title=element_text(size=7))


## Combined plots:
# grid.arrange(condom_use_time_line_plot, clients_time_line_plot, yc_time_line_plot, 
#              nrow = 3,
#              heights = c(2.75, 2, 2))
# 
# grid.arrange(condom_use_time_point_plot, clients_time_point_plot, yc_time_point_plot, 
#              nrow = 3,
#              heights = c(2.75, 2, 2))

jpeg(file = "../risk_comp_line.jpg", height = 4.5, width = 4, unit = 'in', res = 500)
print(plot_grid(condom_use_time_line_plot, clients_time_line_plot, yc_time_line_plot,nrow = 3, align = "v", rel_heights = c(1.1, 1, 1), labels = c("(A)", "(B)", "(C)"), label_size = 6, label_y = 0.2))
dev.off()


jpeg(file = "../risk_comp_point.jpg", height = 4.5, width = 4, unit = 'in', res = 500)
print(plot_grid(condom_use_time_point_plot, clients_time_point_plot, yc_time_point_plot, 
          nrow = 3, align = "v", rel_heights = c(1.1, 1, 1), labels = c("(A)", "(B)", "(C)"), label_size = 6, label_y = 0.2))
dev.off()


# plots <- list(condom_use_time_line_plot, clients_time_line_plot, yc_time_line_plot)
# grobs <- list()
# widths <- list()
# 
# for (i in 1:length(plots)){
#   grobs[[i]] <- ggplotGrob(plots[[i]])
#   widths[[i]] <- grobs[[i]]$widths[2:5]
# }
# maxwidth <- do.call(grid::unit.pmax, widths)
# for (i in 1:length(grobs)){
#   grobs[[i]]$widths[2:5] <- as.list(maxwidth)
# }
# do.call("grid.arrange", c(grobs, ncol = 1))

## Get legends to align
# plot_grid(
#   condom_use_time_line_plot, 
#   clients_time_line_plot +
#     geom_line(aes(color = "Test")) +
#     scale_color_manual(values = NA) +
#     theme(legend.text = element_blank()
#           , legend.title = element_blank()),
#   yc_time_line_plot + 
#     geom_line(aes(color = "Test")) +
#     scale_color_manual(values = NA) +
#     theme(legend.text = element_blank()
#           , legend.title = element_blank())
#   , align = "hv"
#   , ncol = 1
# )


## Questions about number of partners vary. Last 7 days for clients, last month for main partners
# Quelle est votre fréquence d’utilisation des préservatifs durant ce dernier mois lors de vos rapports
# sexuels (par voie vaginale ou anale) avec vos clients?

# Combien de partenaires sexuels principaux avez-vous eu le mois dernier.

## Questions about condom use are all in the last month
# Quelle est votre fréquence d’utilisation des préservatifs durant ce dernier mois lors de vos rapports
# sexuels (par voie vaginale ou anale) avec vos clients?

# Quelle était votre fréquence d’utilisation des préservatifs le mois dernier lorsque vous aviez des
# rapports sexuels (par voie vaginale ou anale) votre (vos) partenaire(s) principal (aux)?
  
## There is also a question about how long they have constantly used condoms for (among those who responded that they always use condoms)
# Depuis combien de temps utilisez-vous constamment les préservatifs (toujours) lorsque vous faites
# des rapports sexuels (par voie vaginale ou anale) avec vos clients?

## Condom use with clients

## Condom use with main partner

## No "main partner": N = 4 (14.4%)

############################################################################################
## Missingness
############################################################################################
presence_dat_merge <- presence_dat[presence_dat$visit != "scr", ] ## Need the merge to match up with Yc, which has scr for what I think really should be j0. Not sure about this and need to check
presence_dat_merge$visit[presence_dat_merge$visit == "j0"] <- "scr"

yc_dat <- left_join(yc_dat, presence_dat_merge, by = c("id", "site", "visit"))

yc_dat[yc_dat$presence == 0,] ## Note that all of the variables are missing when presence is missing. But the fact that they have a swab means that they showed up!

## Missingness by visit
yc_dat %>%
  group_by(visit) %>%
  summarise(present = sum(presence), absent = sum(1 - presence), total = n(), pct = mean(1-presence))

## Note that 1 person missing m12 (or just didn't get a valid visit date...). Also, 9 people missing "scr" (or just didn't get a valid visit date), but that could really be j0

## Missingness by site
yc_dat %>%
  group_by(site) %>%
  summarise(present = sum(presence), absent = sum(1 - presence), total = n(), pct = mean(1-presence))

## Missingness by calendar period

## Number missing all sexual behavior questions
length(which(apply(yc_dat[, which(names(yc_dat) == "FreqUtilervatifPart"):which(names(yc_dat) == "UtilisePvatifDepuis")], 1, function(x) sum(is.na(x))) == length(which(names(yc_dat) == "FreqUtilervatifPart"):which(names(yc_dat) == "UtilisePvatifDepuis"))))
## 13 individuals - so there are 3 who have valid visit dates but are missing all sexual behavior variables

############################################################################################
## Repeated Swab Agreement
############################################################################################
## Number of swabs per woman
yc_dat %>%
  add_count(id) %>%
  group_by(n) %>%
  select(id, n) %>%
  distinct(id, n) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  mutate(pct = count/(sum(count)))

long_dat <- yc_dat %>%
  add_count(id) %>%
  filter(n > 1) %>% 
  mutate(visit_num = ifelse(visit == "scr", 0, as.numeric(gsub("m", "", visit)))) %>%
  arrange(id, visit_num) %>%
  group_by(id) %>%
  mutate(swab_num = rank(visit_num))
  

## 2x2 table for all women who had only 2 swabs
two_swabs <- long_dat %>% 
                filter(n == 2) %>%
                select(id, swab_num, value) %>%
                mutate(swab_num = paste("swab", swab_num, sep = "_")) %>%
                spread(key = swab_num, value = value)

table(two_swabs$swab_1, two_swabs$swab_2, deparse.level = 2)

## Pre-study vs. post-study for these individuals
long_dat %>% 
  filter(n == 2) %>%
  select(id, visit, value) %>%
  group_by(visit) %>%
  summarise(n = n(), num_pos = sum(value == "P"), num_neg = sum(value == "N"))

## One person had 3 swabs - all were negative
long_dat %>% 
  filter(n == 3) %>%
  select(id, swab_num, value)


############################################################################################
## Predictors of Yc detection
############################################################################################
pred_dat <- base_dat

## Baseline Characteristics
pred_dat$id <- pred_dat$h00IdNumber

## Age
pred_dat$age <- pred_dat$F00_scr_q01AgeDuSujet
pred_dat$age_quart <- cut(pred_dat$age, c(17, 32, 40, 45, 57), labels = c("18-32", "33-40", "41-45", "45-57"))
pred_dat$age_cat <- cut(pred_dat$age, c(17, 29, 39, 49, Inf),  labels = c("<30", "30-39", "40-49", "50+"))

## Registered FSW
pred_dat$registered <- factor(pred_dat$F00_scr_q10aTSEnregistree, levels = c(2, 1), labels = c("No", "Yes"))

## Ethnic Group
pred_dat$ethnic_group_cat <- factor(c("Wolof", "Poulaar", "Serere", "Mandingue", "Other", "Other", "Other", "Other", "Other")[pred_dat$ethnic_group], levels = c("Wolof", "Poulaar", "Serere", "Mandingue", "Other"))

## Marital status
pred_dat$married <- factor(pred_dat$F00_scr_q07StatutMatrimonial, levels = seq(1,5), labels = c("Single", "Married", "Separated", "Divorced", "Widowed"))
pred_dat$married_cat <- factor(c("Never Married", "Currently/Previously Married", "Never Married", "Currently/Previously Married", "Currently/Previously Married")[pred_dat$married], levels = c("Never Married", "Currently/Previously Married"))

## Education
table(pred_dat$F00_scr_q08EtudeFrancaise, exclude = NULL)
table(pred_dat$F00_scr_q09Scolarite, exclude = NULL)
table(pred_dat$F00_scr_q09ScolariteAutre, exclude = NULL)

pred_dat$edu <- ifelse(is.na(pred_dat$F00_scr_q09Scolarite), 0, pred_dat$F00_scr_q09Scolarite)
pred_dat$edu <- factor(pred_dat$edu, levels = c(0, 1, 2), labels = c("None", "Primary", "Secondary"))

## Merge to sexual behavior data
var_names <- c("id", "age", "age_cat", "age_quart", "registered", "ethnic_group", "ethnic_group_cat", "married", "married_cat", "edu")

pred_dat <- left_join(yc_dat, pred_dat[, var_names], by = "id")

## Site
pred_dat$site <- factor(pred_dat$site)

## Positive Yc
pred_dat$pos_yc <- ifelse(pred_dat$value == "P", 1, 0)

## Sample random result for women with multiple swabs
pred_dat_rand <- pred_dat %>% 
  group_by(id) %>%
  sample_n(1)

## Age
pred_dat_rand %>%
  group_by(age_quart) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)

pred_dat_rand %>%
  group_by(age_cat) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)

chisq.test(pred_dat_rand$age_quart, pred_dat_rand$pos_yc)
chisq.test(pred_dat_rand$age_cat, pred_dat_rand$pos_yc)
fisher.test(pred_dat_rand$age_cat, pred_dat_rand$pos_yc)

summary(glm(pos_yc ~ age, data = pred_dat_rand, family = "binomial"))
summary(glm(pos_yc ~ age_quart, data = pred_dat_rand, family = "binomial"))
summary(glm(pos_yc ~ age_cat, data = pred_dat_rand, family = "binomial"))

## Registered FSW
pred_dat_rand %>%
  group_by(registered) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)

chisq.test(pred_dat_rand$registered, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ registered, data = pred_dat_rand, family = "binomial"))


## Ethnic Group
pred_dat_rand %>%
  group_by(ethnic_group_cat) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)

fisher.test(pred_dat_rand$ethnic_group_cat[!is.na(pred_dat_rand$ethnic_group_cat)], pred_dat_rand$pos_yc[!is.na(pred_dat_rand$ethnic_group_cat)])
summary(glm(pos_yc ~ ethnic_group_cat, data = pred_dat_rand, family = "binomial"))

## Marital status
pred_dat_rand %>%
  group_by(married_cat) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
chisq.test(pred_dat_rand$married_cat, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ married_cat, data = pred_dat_rand, family = "binomial"))

pred_dat_rand %>%
  group_by(married) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)

## Education
pred_dat_rand %>%
  group_by(edu) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$edu, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ edu, data = pred_dat_rand, family = "binomial"))

## Site
pred_dat_rand %>%
  group_by(site) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
chisq.test(pred_dat_rand$site, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ site, data = pred_dat_rand, family = "binomial"))

## Condom usage with clients in last month
pred_dat_rand %>%
  group_by(FreqUtilPreservatif) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$FreqUtilPreservatif, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ FreqUtilPreservatif, data = pred_dat_rand, family = "binomial"))

## Collapsed condom categories
pred_dat_rand %>%
  group_by(always_condom_client) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$always_condom_client, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ always_condom_client, data = pred_dat_rand, family = "binomial"))


## Condom usage with main partner in last month
pred_dat_rand %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(FreqUtilervatifPart) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$FreqUtilervatifPart[pred_dat_rand$NbrePartMoisDernier != 0], pred_dat_rand$pos_yc[pred_dat_rand$NbrePartMoisDernier != 0])
summary(glm(pos_yc ~ FreqUtilervatifPart, data = pred_dat_rand[pred_dat_rand$NbrePartMoisDernier != 0, ], family = "binomial"))

## Collapsed condom category
pred_dat_rand %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_main) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$always_condom_main[pred_dat_rand$NbrePartMoisDernier != 0], pred_dat_rand$pos_yc[pred_dat_rand$NbrePartMoisDernier != 0])
summary(glm(pos_yc ~ always_condom_main, data = pred_dat_rand[pred_dat_rand$NbrePartMoisDernier != 0, ], family = "binomial"))

## Number of clients in last 7 days
pred_dat_rand$NbreClieDernSemaine <- as.numeric(pred_dat_rand$NbreClieDernSemaine)
table(pred_dat_rand$NbreClieDernSemaine)
summary(glm(pos_yc ~ NbreClieDernSemaine, data = pred_dat_rand, family = "binomial"))
exp(coefficients(glm(pos_yc ~ NbreClieDernSemaine, data = pred_dat_rand, family = "binomial")))
exp(confint(glm(pos_yc ~ NbreClieDernSemaine, data = pred_dat_rand, family = "binomial")))

pred_dat_rand$num_clients_cat <- cut(pred_dat_rand$NbreClieDernSemaine, breaks = c(-Inf, 0, 1, 2, 3, Inf ), labels = c("0", "1", "2", "3", "> 3"))
pred_dat_rand %>%
  group_by(num_clients_cat) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$num_clients_cat, pred_dat_rand$pos_yc)

## Reports a main partner in the last month
pred_dat_rand %>%
  group_by(reports_main_partner) %>%
  summarise(num_pos = sum(pos_yc), n = n()) %>%
  mutate(pct = num_pos/n)
fisher.test(pred_dat_rand$reports_main_partner, pred_dat_rand$pos_yc)
summary(glm(pos_yc ~ reports_main_partner, data = pred_dat_rand, family = "binomial"))

############################################################################################
## Table S1: Demographics and baseline characterstics of all FSW in PrEP Demo Project
############################################################################################
## Number of FSW 
nrow(dat) ## 131

## Age (median, IQR)
length(which(is.na(dat$F00_scr_q01AgeDuSujet)))
quantile(dat$F00_scr_q01AgeDuSujet, probs = c(0.25, 0.5, 0.75))
IQR(dat$F00_scr_q01AgeDuSujet)

## Born in Senegal (%), and birth locations of individuals born outside of Senegal
mean(dat$F00_scr_q04Nationalite == 1)
dat$F00_scr_q04NationaliteAutre[dat$F00_scr_q04Nationalite == 2]

## Registered vs. unregistered FSW (%)
mean(dat$F00_scr_q10aTSEnregistree == 1, na.rm = TRUE)
length(which(is.na(dat$F00_scr_q10aTSEnregistree)))

## Site (N, %)
table(dat$site)
table(dat$site)/nrow(dat)

## Ethnic Group, and number missing
table(dat$F00_scr_q03Ethnie, dat$F00_scr_q03EthnieAutre, exclude = NULL)

dat$ethnic_group <- factor(dat$F00_scr_q03Ethnie, levels = c(1:8, 99), labels = c("Wolof", "Poulaar", "Serere", "Mandingue/Bambara", "Diola", "Soninke", "Maure", "Manjak", "Other")) ## Note that Geoff's tables has alternative names and spellings for these groups
dat$ethnic_group[dat$F00_scr_q03EthnieAutre == "MANDJAQUE"] <- "Manjak"
table(dat$ethnic_group, exclude = NULL)
table(dat$ethnic_group)/length(which(!is.na(dat$ethnic_group)))

## Self-reported number of sex partners in the prior week (Median, range), and number missing - I think this is really number of clients, not number of partners, right?
dat$F14_j0_q03NbreClieDernSemaine[dat$F14_j0_q03NbreClieDernSemaine %in% c("NE CONNAIT PAS", "REFUS")] <- NA
length(which(is.na(dat$F14_j0_q03NbreClieDernSemaine)))

dat$num_clients_cat <- cut(as.integer(dat$F14_j0_q03NbreClieDernSemaine), breaks = c(-Inf, 0, 2, Inf ), labels = c("0", "1-2", "3+"))
dat %>%
  group_by(num_clients_cat) %>%
  summarise(n = n()) %>%
  filter(!is.na(num_clients_cat)) %>%
  mutate(pct = n/sum(n))

## Self-reported had a "Main Partner" - this is in the last six months
table(dat$F14_j0_q11PartenaireSexPrinc, exclude = NULL)
sum(dat$F14_j0_q11PartenaireSexPrinc == "o", na.rm = TRUE)
mean(dat$F14_j0_q11PartenaireSexPrinc == "o", na.rm = TRUE)

## Self-reported condom use with clients in the last month for vaginal or anal sex, and number missing
table(dat$F14_j0_q05FreqUtilPreservatif, exclude = NULL)
table(dat$F14_j0_q05FreqUtilPreservatif)/length(which(!is.na(dat$F14_j0_q05FreqUtilPreservatif)))

## Self-reported condom use with main partner in the last month for vaginal or anal sex, and number missing
table(dat$F14_j0_q15FreqUtilervatifPart[dat$F14_j0_q11PartenaireSexPrinc == "o"], exclude = NULL)
table(dat$F14_j0_q15FreqUtilervatifPart[dat$F14_j0_q11PartenaireSexPrinc == "o"])/length(which(!is.na(dat$F14_j0_q15FreqUtilervatifPart[dat$F14_j0_q11PartenaireSexPrinc == "o"])))


############################################################################################
## Frequency of STI detection
############################################################################################
sti_question_stems <- c("RaisonPPrelevement", "EstDispost", "RaisonPat", "Resultatst", "TypeTest", "Specificst", "atestSyphilis", "aResultatRPR", "bResultatTPHA")
sti_dat <- dat[, c("h00IdNumber", "site", grep(paste(sti_question_stems, collapse = "|"), names(dat), value = TRUE))]

## Visit date is messing something up
##grep("VisitDate", names(dat), value = TRUE)

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

## Gonorrhea
table(sti_dat$EstDispostGonorrhee, exclude = NULL)
table(sti_dat$ResultatstGonorrhee, exclude = NULL)
mean(sti_dat$ResultatstGonorrhee == "positif", na.rm = TRUE)


## At baseline
table(sti_dat$EstDispostGonorrhee[sti_dat$visit == "scr"], exclude = NULL)
table(sti_dat$ResultatstGonorrhee[sti_dat$visit == "scr"], exclude = NULL)
mean(sti_dat$ResultatstGonorrhee[sti_dat$visit == "scr"] == "positif", na.rm = TRUE)
## Only 3 gonorrhea results, all negative...

## Post baseline
table(sti_dat$EstDispostGonorrhee[sti_dat$visit != "scr"], exclude = NULL)
table(sti_dat$ResultatstGonorrhee[sti_dat$visit != "scr"], exclude = NULL)
mean(sti_dat$ResultatstGonorrhee[sti_dat$visit != "scr"] == "positif", na.rm = TRUE)
## 8/118 (6.8%)

## First test
sti_dat %>%
  arrange(id, visit_factor) %>%
  group_by(id) %>%
  summarise(first_result = ResultatstGonorrhee[which(!is.na(ResultatstGonorrhee))[1]]) %>%
  ungroup() %>%
  summarise(n = sum(!is.na(first_result)), num_pos = sum(first_result == "positif", na.rm = TRUE), num_neg = sum(first_result == "negatif", na.rm = TRUE), pct_pos = mean(first_result == "positif", na.rm = TRUE))

## Follow-up gonorrhea tests
sti_dat %>%
  arrange(id, visit_factor) %>%
  filter(!is.na(ResultatstGonorrhee)) %>%
  group_by(id) %>%
  mutate(test_num = 1:n()) %>%
  mutate(first_result = ResultatstGonorrhee[which(test_num == 1)]) %>%
  mutate(prev_result = lag(ResultatstGonorrhee)) %>%
  filter(prev_result == "negatif") %>%
  ungroup() %>%
  summarise(n = sum(!is.na(ResultatstGonorrhee)), num_pos = sum(ResultatstGonorrhee == "positif", na.rm = TRUE), num_neg = sum(ResultatstGonorrhee == "negatif", na.rm = TRUE), pct_pos = mean(ResultatstGonorrhee == "positif", na.rm = TRUE))
  

## Incident gonorrhea post baseline
ids_gon <- unique(na.omit((sti_dat$id[sti_dat$visit == "scr" & sti_dat$ResultatstGonorrhee == "negatif"])))
table(sti_dat$EstDispostGonorrhee[sti_dat$visit != "scr" & sti_dat$id %in% ids_gon], exclude = NULL)
table(sti_dat$ResultatstGonorrhee[sti_dat$visit != "scr" & sti_dat$id %in% ids_gon], exclude = NULL)
mean(sti_dat$ResultatstGonorrhee[sti_dat$visit != "scr" & sti_dat$id %in% ids_gon] == "positif", na.rm = TRUE)
## Only 2 f/u results for people who tested negative at baseline, both of which were negative

## Chlamydia
table(sti_dat$EstDispostChlamydia, exclude = NULL)
table(sti_dat$ResultatstChlamydia, exclude = NULL)
mean(sti_dat$ResultatstChlamydia == "positif", na.rm = TRUE)

## At baseline
table(sti_dat$EstDispostChlamydia[sti_dat$visit == "scr"], exclude = NULL)
table(sti_dat$ResultatstChlamydia[sti_dat$visit == "scr"], exclude = NULL)
mean(sti_dat$ResultatstChlamydia[sti_dat$visit == "scr"] == "positif", na.rm = TRUE)
## Only 3 chlamydia results, all negative...

## Post baseline
table(sti_dat$EstDispostChlamydia[sti_dat$visit != "scr"], exclude = NULL)
table(sti_dat$ResultatstChlamydia[sti_dat$visit != "scr"], exclude = NULL)
mean(sti_dat$ResultatstChlamydia[sti_dat$visit != "scr"] == "positif", na.rm = TRUE)
## 8/120 (6.7%)

## First test
sti_dat %>%
  arrange(id, visit_factor) %>%
  group_by(id) %>%
  summarise(first_result = ResultatstChlamydia[which(!is.na(ResultatstChlamydia))[1]]) %>%
  ungroup() %>%
  summarise(n = sum(!is.na(first_result)), num_pos = sum(first_result == "positif", na.rm = TRUE), num_neg = sum(first_result == "negatif", na.rm = TRUE), pct_pos = mean(first_result == "positif", na.rm = TRUE))

## Follow-up chlamydia tests
sti_dat %>%
  arrange(id, visit_factor) %>%
  filter(!is.na(ResultatstChlamydia)) %>%
  group_by(id) %>%
  mutate(test_num = 1:n()) %>%
  mutate(first_result = ResultatstChlamydia[which(test_num == 1)]) %>%
  mutate(prev_result = lag(ResultatstChlamydia)) %>%
  filter(prev_result == "negatif") %>%
  ungroup() %>%
  summarise(n = sum(!is.na(ResultatstChlamydia)), num_pos = sum(ResultatstChlamydia == "positif", na.rm = TRUE), num_neg = sum(ResultatstChlamydia == "negatif", na.rm = TRUE), pct_pos = mean(ResultatstChlamydia == "positif", na.rm = TRUE))

## Incident chlamydia post baseline
ids_chl <- unique(na.omit((sti_dat$id[sti_dat$visit == "scr" & sti_dat$ResultatstChlamydia == "negatif"])))
table(sti_dat$EstDispostChlamydia[sti_dat$visit != "scr" & sti_dat$id %in% ids_chl], exclude = NULL)
table(sti_dat$ResultatstChlamydia[sti_dat$visit != "scr" & sti_dat$id %in% ids_chl], exclude = NULL)
mean(sti_dat$ResultatstChlamydia[sti_dat$visit != "scr" & sti_dat$id %in% ids_chl] == "positif", na.rm = TRUE)
## Only 3 f/u results for people who tested negative at baseline, both of which were negative

## Syphilis
table(sti_dat$atestSyphilisRPR, exclude = NULL)
table(sti_dat$aResultatRPR, exclude = NULL)
mean(sti_dat$aResultatRPR == "positif", na.rm = TRUE)

table(sti_dat$atestSyphilisTPHA, exclude = NULL)
table(sti_dat$bResultatTPHA, exclude = NULL)
mean(sti_dat$bResultatTPHA == "positif", na.rm = TRUE)

table(sti_dat$aResultatRPR, sti_dat$bResultatTPHA, exclude = NULL, deparse.level = 2)
21/521
21/24

## RPR +, then TPPA +. Exclude id that had positive RPR at baseline
ids_pos_base_rpr <- unique(na.omit(sti_dat$id[sti_dat$aResultatRPR == "positif" & sti_dat$visit == "scr"])) ## Just one to exclude

table(sti_dat$aResultatRPR[(!sti_dat$id %in% ids_pos_base_rpr)])
## 22/523 (4.2%) of RPRs post-baseline were positive 
table(sti_dat$bResultatTPHA[(!sti_dat$id %in% ids_pos_base_rpr) & sti_dat$aResultatRPR == "positif"])
## Of 22 positive post-baseline RPRs, 21 (96%) of concurrent TPHAs were also positive

## Baseline syphilis
## TPHA
table(sti_dat$atestSyphilisTPHA[sti_dat$visit == "scr"], exclude = NULL)
table(sti_dat$bResultatTPHA[sti_dat$visit == "scr"], exclude = NULL)
mean(sti_dat$bResultatTPHA[sti_dat$visit == "scr"] == "positif", na.rm = TRUE)
## 7/54 had positve TPHA at baseline

## Post-baseline TPHA
table(sti_dat$atestSyphilisTPHA[sti_dat$visit != "scr"], exclude = NULL)
table(sti_dat$bResultatTPHA[sti_dat$visit != "scr"], exclude = NULL)
mean(sti_dat$bResultatTPHA[sti_dat$visit != "scr"] == "positif", na.rm = TRUE)
## 56/467 had positve TPHA post baseline

## RPR
table(sti_dat$atestSyphilisRPR[sti_dat$visit == "scr"], exclude = NULL)
table(sti_dat$aResultatRPR[sti_dat$visit == "scr"], exclude = NULL)
mean(sti_dat$aResultatRPR[sti_dat$visit == "scr"] == "positif", na.rm = TRUE)
## 1/54 were positive at baseline

table(sti_dat$aResultatRPR[sti_dat$visit == "scr"], sti_dat$bResultatTPHA[sti_dat$visit == "scr"], exclude = NULL, deparse.level = 2)

## Post-baseline RPR
table(sti_dat$atestSyphilisRPR[sti_dat$visit != "scr"], exclude = NULL)
table(sti_dat$aResultatRPR[sti_dat$visit != "scr"], exclude = NULL)
mean(sti_dat$aResultatRPR[sti_dat$visit != "scr"] == "positif", na.rm = TRUE)
## 23/472 were positive after baseline

table(sti_dat$aResultatRPR[sti_dat$visit != "scr"], sti_dat$bResultatTPHA[sti_dat$visit != "scr"], exclude = NULL, deparse.level = 2)

## Incident RPR
## First test
sti_dat %>%
  arrange(id, visit_factor) %>%
  group_by(id) %>%
  summarise(first_result = aResultatRPR[which(!is.na(aResultatRPR))[1]]) %>%
  ungroup() %>%
  summarise(n = sum(!is.na(first_result)), num_pos = sum(first_result == "positif", na.rm = TRUE), num_neg = sum(first_result == "negatif", na.rm = TRUE), pct_pos = mean(first_result == "positif", na.rm = TRUE))


## Including TPHA
sti_dat %>%
  arrange(id, visit_factor) %>%
  group_by(id) %>%
  summarise(first_result = aResultatRPR[which(!is.na(aResultatRPR))[1]], tpha = bResultatTPHA[which(!is.na(aResultatRPR))[1]]) %>%
  ungroup() %>%
  filter(first_result == "positif") %>%
  summarise(n = sum(!is.na(tpha)), num_pos = sum(tpha == "positif", na.rm = TRUE), num_neg = sum(tpha == "negatif", na.rm = TRUE), pct_pos = mean(tpha == "positif", na.rm = TRUE))


## Follow-up RPR tests
sti_dat %>%
  arrange(id, visit_factor) %>%
  filter(!is.na(aResultatRPR)) %>%
  group_by(id) %>%
  mutate(test_num = 1:n()) %>%
  mutate(first_result = aResultatRPR[which(test_num == 1)]) %>%
  mutate(prev_result = lag(aResultatRPR)) %>%
  filter(prev_result == "negatif") %>%
  ungroup() %>%
  summarise(n = sum(!is.na(aResultatRPR)), num_pos = sum(aResultatRPR == "positif", na.rm = TRUE), num_neg = sum(aResultatRPR == "negatif", na.rm = TRUE), pct_pos = mean(aResultatRPR == "positif", na.rm = TRUE))

## Including TPHA
sti_dat %>%
  arrange(id, visit_factor) %>%
  filter(!is.na(aResultatRPR)) %>%
  group_by(id) %>%
  mutate(test_num = 1:n()) %>%
  mutate(first_result = aResultatRPR[which(test_num == 1)]) %>%
  mutate(prev_result = lag(aResultatRPR)) %>%
  filter(prev_result == "negatif") %>%
  ungroup() %>%
  filter(aResultatRPR == "positif") %>%
  summarise(n = sum(!is.na(bResultatTPHA)), num_pos = sum(bResultatTPHA == "positif", na.rm = TRUE), num_neg = sum(bResultatTPHA == "negatif", na.rm = TRUE), pct_pos = mean(bResultatTPHA == "positif", na.rm = TRUE))
  

## Incident RPR after baseline
ids_rpr <- unique(na.omit(sti_dat$id[sti_dat$aResultatRPR == "negatif" & sti_dat$visit == "scr"]))
length(ids_rpr)
table(sti_dat$atestSyphilisRPR[sti_dat$visit != "scr" & sti_dat$id %in% ids_rpr], exclude = NULL)
table(sti_dat$aResultatRPR[sti_dat$visit != "scr" & sti_dat$id %in% ids_rpr], exclude = NULL)
mean(sti_dat$aResultatRPR[sti_dat$visit != "scr" & sti_dat$id %in% ids_rpr] == "positif", na.rm = TRUE)
## Of 53 women with negative RPR at baseline, 3 out of 81 subsequent tests were positive for RPR
length(unique(na.omit(sti_dat$id[sti_dat$aResultatRPR == "positif" & sti_dat$visit != "scr" & sti_dat$id %in% ids_rpr])))
length(unique(na.omit(sti_dat$id[sti_dat$aResultatRPR %in% c("positif", "negatif") & sti_dat$visit != "scr" & sti_dat$id %in% ids_rpr])))
## Of 53 women with negative RPR at baseline, 44 were tested at least once subsequently with RPR. Of these, two (4.6%) had positive RPRs.

## Incident TPHA
## First test
sti_dat %>%
  arrange(id, visit_factor) %>%
  group_by(id) %>%
  summarise(first_result = bResultatTPHA[which(!is.na(bResultatTPHA))[1]]) %>%
  ungroup() %>%
  summarise(n = sum(!is.na(first_result)), num_pos = sum(first_result == "positif", na.rm = TRUE), num_neg = sum(first_result == "negatif", na.rm = TRUE), pct_pos = mean(first_result == "positif", na.rm = TRUE))

## Follow-up TPHA tests
sti_dat %>%
  arrange(id, visit_factor) %>%
  filter(!is.na(bResultatTPHA)) %>%
  group_by(id) %>%
  mutate(test_num = 1:n()) %>%
  mutate(first_result = bResultatTPHA[which(test_num == 1)]) %>%
  mutate(prev_result = lag(bResultatTPHA)) %>%
  filter(prev_result == "negatif") %>%
  ungroup() %>%
  summarise(n = sum(!is.na(bResultatTPHA)), num_pos = sum(bResultatTPHA == "positif", na.rm = TRUE), num_neg = sum(bResultatTPHA == "negatif", na.rm = TRUE), pct_pos = mean(bResultatTPHA == "positif", na.rm = TRUE))


## Incident TPHA - post baseline
ids_tpha <- unique(na.omit(sti_dat$id[sti_dat$bResultatTPHA == "negatif" & sti_dat$visit == "scr"]))
length(ids_tpha)
table(sti_dat$atestSyphilisTPHA[sti_dat$visit != "scr" & sti_dat$id %in% ids_tpha], exclude = NULL)
table(sti_dat$bResultatTPHA[sti_dat$visit != "scr" & sti_dat$id %in% ids_tpha], exclude = NULL)
mean(sti_dat$bResultatTPHA[sti_dat$visit != "scr" & sti_dat$id %in% ids_tpha] == "positif", na.rm = TRUE)
## Of 47 women with negative TPHA at baseline, 1 out of 71 (1.4%) subsequent tests were positive for TPHA
length(unique(na.omit(sti_dat$id[sti_dat$bResultatTPHA == "positif" & sti_dat$visit != "scr" & sti_dat$id %in% ids_tpha])))
length(unique(na.omit(sti_dat$id[sti_dat$bResultatTPHA %in% c("positif", "negatif") & sti_dat$visit != "scr" & sti_dat$id %in% ids_tpha])))
## Of 47 women with negative TPHA at baseline, 39 were tested at least once subsequently with TPHA. Of these, one (2.6%) had positive TPHA.


## Number with available results
sti_dat %>%
  group_by(visit_factor) %>%
  summarise(n_gon = sum(!is.na(ResultatstGonorrhee)),
            n_chl = sum(!is.na(ResultatstChlamydia)),
            n_rpr = sum(!is.na(aResultatRPR)),
            n_tpha = sum(!is.na(bResultatTPHA)),
            n_both = sum(!is.na(aResultatRPR) & !is.na(bResultatTPHA))) %>%
  mutate(pct_gon = n_gon/350,
         pct_chl = n_chl/350,
         pct_rpr = n_rpr/350,
         pct_tpha = n_tpha/350,
         pct_both = n_both/350)
                         






## Merge with yc_dat
yc_dat <- left_join(yc_dat, sti_dat, by = c("id", "visit", "site"))

## Gonorrhea
table(yc_dat$EstDispostGonorrhee, exclude = NULL)
table(yc_dat$ResultatstGonorrhee, exclude = NULL)
mean(yc_dat$ResultatstGonorrhee == "positif", na.rm = TRUE)

table(yc_dat$ResultatstGonorrhee, yc_dat$value, deparse.level = 2)

## Chlamydia
table(yc_dat$EstDispostChlamydia, exclude = NULL)
table(yc_dat$ResultatstChlamydia, exclude = NULL)
mean(yc_dat$ResultatstChlamydia == "positif", na.rm = TRUE)

table(yc_dat$ResultatstChlamydia, yc_dat$value, deparse.level = 2)

## Syphilis
table(yc_dat$atestSyphilisRPR, exclude = NULL)
table(yc_dat$aResultatRPR, exclude = NULL)
mean(yc_dat$aResultatRPR == "positif", na.rm = TRUE)

table(yc_dat$aResultatRPR, yc_dat$value, deparse.level = 2)

table(yc_dat$atestSyphilisTPHA, exclude = NULL)
table(yc_dat$bResultatTPHA, exclude = NULL)
mean(yc_dat$bResultatTPHA == "positif", na.rm = TRUE)

table(yc_dat$bResultatTPHA, yc_dat$value, deparse.level = 2)

table(yc_dat$aResultatRPR, yc_dat$bResultatTPHA, exclude = NULL, deparse.level = 2)

## Sample
# DatePrelevement
# RaisonPPrelevement

## Gonorrhea
## grep("Gonorrhee", names(dat), value = TRUE)
# EstDispostGonorrhee
# RaisonPatGonorrhee
# ResultatstGonorrhee
# TypeTestGonorrhee
# SpecificstGonorrhee

## Chlamydia
# EstDispostChlamydia
# RaisonPatChlamydia
# ResultatstChlamydia
# SpecificstChlamydia
# TypeTestChlamydia

## Syphilis
# EstDispoestSyphilis
# atestSyphilisTPHA
# aResultatRPR
# atestSyphilisRPR
# bResultatTPHA

############################################################################################
## Extensions
############################################################################################
## Frequency of misreporting, either direction

## Predictors of misreporting - either direction? Overreporting condom use?

