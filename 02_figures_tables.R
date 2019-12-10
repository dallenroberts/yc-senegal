library(tidyverse)
library(gridExtra)
library(cowplot)
library(binom)

############################################################################################
## Results in text
############################################################################################
## Number of Yc samples
nrow(yc_dat) ## 154

## Number of women tested for Yc
length(unique(yc_dat$id)) ## 121

## Percentage of women with positive Yc result
yc_woman_dat <- yc_dat %>%
  group_by(id) %>%
  summarise(has_positive = any(value == "P")) %>%
  ungroup() %>%
  summarise(num_positive = sum(has_positive), pct = sum(has_positive)/n(), n = n())
yc_woman_dat
binom.confint(yc_woman_dat$num_positive, yc_woman_dat$n, methods = "wilson")

## Percentage of swabs with positive Yc result
yc_swab_dat <- yc_dat %>%
  summarise(num_positive = sum(value == "P"), pct = sum(value == "P")/n(), n = n())
yc_swab_dat
binom.confint(yc_swab_dat$num_positive, yc_swab_dat$n, methods = "wilson")

## Number of Yc samples that were positive
length(which(yc_dat$value == "P")) ## 34

## Percentage of Yc samples that were positive
length(which(yc_dat$value == "P"))/nrow(yc_dat) ## 22

## STI testing results
## Percent who were either positive at baseline or became positive over the course of the study
sti_dat %>%
  group_by(id) %>%
  summarise(pos_gon = any(ResultatstGonorrhee == "positif", na.rm = TRUE), tested_gon = any(!is.na(ResultatstGonorrhee)),
            pos_chl = any(ResultatstChlamydia == "positif", na.rm = TRUE), tested_chl = any(!is.na(ResultatstChlamydia)),
            pos_tpha = any(bResultatTPHA == "positif", na.rm = TRUE), tested_tpha = any(!is.na(bResultatTPHA))) %>%
  ungroup() %>%
  select(-id) %>%
  colSums

## Gonorrhea: 8/106 = 7.54%
## Chlamyida: 8/106 = 7.54%
## TPHA: 34/221 = 15.4%

## Testing at baseline
sti_dat %>%
  filter(visit == "scr") %>%
  summarise(pos_gon = sum(ResultatstGonorrhee == "positif", na.rm = TRUE), tested_gon = sum(!is.na(ResultatstGonorrhee)),
            pos_chl = sum(ResultatstChlamydia == "positif", na.rm = TRUE), tested_chl = sum(!is.na(ResultatstChlamydia)),
            pos_tpha = sum(bResultatTPHA == "positif", na.rm = TRUE), tested_tpha = sum(!is.na(bResultatTPHA)))

############################################################################################
## Table 1: Demographics and baseline characterstics of FSW enrolled in PrEP Demo Project
############################################################################################
## Number of FSW enrolled
nrow(dat) ## 267 women enrolled

## Age (median, IQR)
length(which(is.na(dat$age)))
quantile(dat$age, probs = c(0.25, 0.5, 0.75))
IQR(dat$age)

## Born in Senegal (%), and birth locations of individuals born outside of Senegal
length(which(dat$F00_scr_q04Nationalite == 1))
mean(dat$F00_scr_q04Nationalite == 1)
dat$F00_scr_q04NationaliteAutre[dat$F00_scr_q04Nationalite == 2]

## Registered vs. unregistered FSW (%)
length(which(dat$registered == "Yes"))
mean(dat$registered == "Yes", na.rm = TRUE)
length(which(is.na(dat$registered)))

## Site (N, %)
table(dat$site)
table(dat$site)/nrow(dat)

## Ethnic Group, and number missing
table(dat$F00_scr_q03Ethnie, dat$F00_scr_q03EthnieAutre, exclude = NULL)
table(dat$ethnic_group, exclude = NULL)
table(dat$ethnic_group)/length(which(!is.na(dat$ethnic_group)))

## Self-reported number of clients in the prior week (Median, range), and number missing
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
table(dat$F14_j0_q15FreqUtilervatifPart[which(dat$F14_j0_q11PartenaireSexPrinc == "o")], exclude = NULL)
table(dat$F14_j0_q15FreqUtilervatifPart[which(dat$F14_j0_q11PartenaireSexPrinc == "o")])/length(which(!is.na(dat$F14_j0_q15FreqUtilervatifPart[which(dat$F14_j0_q11PartenaireSexPrinc == "o")])))

## Confidence in condom use the next time having sex
## With main partner
dat$ConfianceservParten<- factor(dat$F14_j0_q18ConfianceservParten, levels = c("pasdutout", "moins", "confiante", "tres"), labels = c("Not at all", "Less confident", "Confident", "Very confident"))
table(dat$ConfianceservParten, exclude = NULL) 
table(dat$ConfianceservParten[which(dat$F14_j0_q11PartenaireSexPrinc == "o")], exclude = NULL) 
table(dat$ConfianceservParten[which(dat$F14_j0_q11PartenaireSexPrinc == "o")])/length(which(!is.na(dat$ConfianceservParten[which(dat$F14_j0_q11PartenaireSexPrinc == "o")])))

## With client
dat$ConfiancePreserv <- factor(dat$F14_j0_q08ConfiancePreserv, levels = c("pasdutout", "moins", "confiante", "tres"), labels = c("Not at all", "Less confident", "Confident", "Very confident"))
table(dat$ConfiancePreserv, exclude = NULL)
table(dat$ConfiancePreserv) /length(which(!is.na(dat$ConfiancePreserv)))

############################################################################################
## Figure 1: Trends in self-reported sexual behavior over time since enrollment 
############################################################################################

## Proportion NOT always using a condom 
condom_dat <- sex_dat %>%
  select(id, visit_factor, always_condom_client, always_condom_main) %>%
  gather(key = "measure", value = "value", -id, -visit_factor) %>%
  mutate(measure = fct_recode(measure, "Main Partner" = "always_condom_main", "Clients" = "always_condom_client")) %>%
  group_by(measure, visit_factor) %>%
  summarise(pct = mean(1-value, na.rm = TRUE), num_inconsistent = sum(1-value, na.rm = TRUE), n = sum(!is.na(value))) %>%
  mutate(visit_num = c(0, 1, 3, 6, 9, 12))

condom_dat[, c("lower", "upper")] <- binom.confint(condom_dat$num_inconsistent, condom_dat$n, methods = "wilson")[, c("lower", "upper")]

## Proportion of participants reporting inconsistent condom use over the last month at baseline and at 12 months
## With clients
condom_dat[condom_dat$measure == "Clients" & condom_dat$visit_factor %in% c("scr", "m12"), ]

## With main partner
condom_dat[condom_dat$measure == "Main Partner" & condom_dat$visit_factor %in% c("scr", "m12"), ]

condom_use_time_point_plot <- ggplot(condom_dat, aes(x = visit_num, y = pct * 100, group = measure)) +
  geom_point(aes(color = measure, shape = measure)) +
  geom_errorbar(aes(ymin = 100*lower, ymax = 100*upper, colour = measure), width = 0.25) +
  labs(x = "Months after enrollment", y = "% reporting \n inconsistent condom use") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  theme_classic() +
  scale_color_manual(values = c("gray60", "black")) +
  scale_shape_manual(values = c(16, 4)) +
  ## ggtitle("Proportion reporting inconsistent condom use among all FSW") +
  theme(plot.title = element_text(size = 8),
        legend.title=element_blank(),
        legend.text = element_text(size = 6),
        legend.position = "bottom",
        axis.text = element_text(size = 6), 
        axis.title=element_text(size=7),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-10,-10,-10,-10))

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

## Number of clients reported at baseline and at month 12
clients_dat[clients_dat$visit_factor %in% c("scr", "m12"), ]

clients_time_point_plot <- ggplot(data = clients_dat, aes(x = visit_num, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.25) +
  theme_classic() +
  labs(x = "Months after enrollment", y = "Reported number of \n clients in past week") +
  scale_y_continuous(limits = c(0, 6), breaks = seq(0, 6, by = 1)) +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  theme(axis.text = element_text(size = 6), 
        axis.title=element_text(size=7))

## Combined plot
fig1 <- plot_grid(condom_use_time_point_plot, clients_time_point_plot,
          nrow = 2, align = "v", rel_heights = c(1.1, 1), labels = c("(A)", "(B)"), label_size = 6, label_y = 0.2)
jpeg(file = "../../yc/fig1.jpg", height = 4.5, width = 4, unit = 'in', res = 1000)
print(fig1)
dev.off()

############################################################################################
## Table 2: Yc detection vs. self-reported behavior
############################################################################################
## Reported sexual partners
## Client in the last week
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(reports_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## Main partner in the last month
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(reports_main_partner, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## Both
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(reports_either, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## Reported condom use in the prior month
## With client
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(always_condom_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## With main partner - only among participants reporting at least one main partner in the last month
sex_dat %>%
  filter(!is.na(value)) %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## With both - only among participants reporting at least one main partner in the last month
sex_dat %>%
  filter(!is.na(value)) %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## Reported condom use at most recent sex
## With clients
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(used_condom_recent_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## With main partners - only among participants reporting at least one main partner in the last month
sex_dat %>%
  filter(!is.na(value)) %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(used_condom_recent_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

## With both - only among participants reporting at least one main partner in the last month
sex_dat %>%
  filter(!is.na(value)) %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(used_condom_recent_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

############################################################################################
## Figure 2: Yc detection over time
############################################################################################
yc_time_dat <- yc_dat %>%
  mutate(visit = factor(visit, levels = c("scr", "m1", "m3", "m6", "m9", "m12"), labels = c("0", "1", "3", "6", "9", "12"))) %>%
  group_by(visit) %>%
  summarise(num_swabs = n(), num_swabs_pos = sum(value == "P"), pct_swabs_pos = mean(value == "P")) %>%
  mutate(visit_num = as.numeric(as.character(visit)))

yc_time_dat[, c("lower", "upper")] <- binom.confint(yc_time_dat$num_swabs_pos, yc_time_dat$num_swabs, methods = "wilson")[, c("lower", "upper")]

yc_time_dat

fig2 <- ggplot(data = yc_time_dat, aes(x = visit_num, y = 100*pct_swabs_pos)) +
  geom_point() +
  geom_errorbar(aes(ymin = 100*lower, ymax = 100*upper), width = 0.2) +
  theme_classic() +
  labs(x = "Months after enrollment", y = "% of swabs with \n Yc-DNA detected") +
  scale_x_continuous(limits = c(-0.25, 12.25), breaks = seq(0, 12, by = 3)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  theme(axis.text = element_text(size = 8), 
       axis.title=element_text(size=8))

jpeg(file = "../../yc/fig2.jpg", height = 2.5, width = 4, unit = 'in', res = 1000)
print(fig2)
dev.off()

############################################################################################
## Table S1: Demographics and baseline characteristics of FSW that had Yc testing
############################################################################################
## Baseline characterstics dataset - note that this does not include woman for whom swab result was undefined
base_dat_yc <- dat[dat$h00IdNumber %in% unique(yc_dat$id), ]

## Number of FSW tested for yc
nrow(base_dat_yc) ## 121

## Age (median, IQR)
length(which(is.na(base_dat_yc$age)))
quantile(base_dat_yc$age, probs = c(0.25, 0.5, 0.75))
IQR(base_dat_yc$age)

## Born in Senegal (%), and birth locations of individuals born outside of Senegal
length(which(base_dat_yc$F00_scr_q04Nationalite == 1))
mean(base_dat_yc$F00_scr_q04Nationalite == 1)
base_dat_yc$F00_scr_q04NationaliteAutre[base_dat_yc$F00_scr_q04Nationalite == 2]

## Registered vs. unregistered FSW (%)
length(which(base_dat_yc$registered == "Yes"))
mean(base_dat_yc$registered == "Yes")

## Site (N, %)
table(base_dat_yc$site)
table(base_dat_yc$site)/nrow(base_dat_yc)

## Ethnic Group, and number missing
table(base_dat_yc$F00_scr_q03Ethnie, base_dat_yc$F00_scr_q03EthnieAutre, exclude = NULL)
table(base_dat_yc$ethnic_group, exclude = NULL)
table(base_dat_yc$ethnic_group)/length(which(!is.na(base_dat_yc$ethnic_group)))

## Self-reported number of clients in the prior week (Median, range), and number missing 
base_dat_yc$num_clients_cat <- cut(as.integer(base_dat_yc$F14_j0_q03NbreClieDernSemaine), breaks = c(-Inf, 0, 2, Inf ), labels = c("0", "1-2", "3+"))
base_dat_yc %>%
  group_by(num_clients_cat) %>%
  summarise(n = n()) %>%
  filter(!is.na(num_clients_cat)) %>%
  mutate(pct = n/sum(n))
length(which(is.na(base_dat_yc$num_clients_cat)))

## Self-reported had a "Main Partner" - in the last six months
table(base_dat_yc$F14_j0_q11PartenaireSexPrinc, exclude = NULL)
sum(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o", na.rm = TRUE)
mean(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o", na.rm = TRUE)

## Self-reported condom use with clients in the last month for vaginal or anal sex, and number missing
table(base_dat_yc$F14_j0_q05FreqUtilPreservatif, exclude = NULL)
table(base_dat_yc$F14_j0_q05FreqUtilPreservatif, exclude = NULL)/length(which(!is.na(base_dat_yc$F14_j0_q05FreqUtilPreservatif)))
table(base_dat_yc$F14_j0_q05FreqUtilPreservatif)/length(which(!is.na(base_dat_yc$F14_j0_q05FreqUtilPreservatif)))

## Self-reported condom use with main partner in the last month for vaginal or anal sex, and number missing
table(base_dat_yc$F14_j0_q15FreqUtilervatifPart[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")], exclude = NULL)
table(base_dat_yc$F14_j0_q15FreqUtilervatifPart[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")], exclude = NULL)
table(base_dat_yc$F14_j0_q15FreqUtilervatifPart[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")])/length(which(!is.na(base_dat_yc$F14_j0_q15FreqUtilervatifPart[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")])))

## Confidence in condom use the next time having sex
## With main partner
# base_dat_yc$ConfianceservParten<- factor(base_dat_yc$F14_j0_q18ConfianceservParten, levels = c("pasdutout", "moins", "confiante", "tres"), labels = c("Not at all", "Less confident", "Confident", "Very confident"))
table(base_dat_yc$ConfianceservParten[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")], exclude = NULL) 
table(base_dat_yc$ConfianceservParten[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")])/length(which(!is.na(base_dat_yc$ConfianceservParten[which(base_dat_yc$F14_j0_q11PartenaireSexPrinc == "o")])))

## With client
# base_dat_yc$ConfiancePreserv <- factor(base_dat_yc$F14_j0_q08ConfiancePreserv, levels = c("pasdutout", "moins", "confiante", "tres"), labels = c("Not at all", "Less confident", "Confident", "Very confident"))
table(base_dat_yc$ConfiancePreserv, exclude = NULL)
table(base_dat_yc$ConfiancePreserv) /length(which(!is.na(base_dat_yc$ConfiancePreserv)))
