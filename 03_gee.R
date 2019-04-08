library(tidyverse)
library(geepack)
library(sandwich)
library(lmtest)


############################################################################################
## Trends in sexual behavior over time
############################################################################################
gee_dat <- sex_dat
gee_dat$visit_num <- c(0, 1, 3, 6, 9, 12)[match(gee_dat$visit, c("scr", "m1", "m3", "m6", "m9", "m12"))]
gee_dat$yc_result <- ifelse(gee_dat$value == "P", 1, 0)

## Number of clients
mod_num_clients <- geeglm(NbreClieDernSemaine ~ visit_num, data = gee_dat, family = "gaussian", id = id, corstr = 'exchangeable')
summary(mod_num_clients)
cc_num_clients <- coef(summary(mod_num_clients))
citab_num_clients <- with(as.data.frame(cc_num_clients),
                 cbind(lwr=Estimate-1.96*Std.err,
                       upr=Estimate+1.96*Std.err))
rownames(citab_num_clients) <- rownames(cc_num_clients)
## beta  = -0.135, 95% confidence interval: [-0.135, 0.0341]
## p-value for trend = 0.24

## Condom use with clients
mod_condom_clients <- geeglm(always_condom_client ~ visit_num, data = gee_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable') 
summary(mod_condom_clients)
cc_condom_clients <- coef(summary(mod_condom_clients))
citab_condom_clients <- with(as.data.frame(cc_condom_clients),
                 cbind(est = exp(Estimate),
                       lwr=exp(Estimate-1.96*Std.err),
                       upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_clients) <- rownames(cc_condom_clients)
## OR: 1.06, 95% CI: [0.997, 1.12]
## p-value: 0.06

## Condom use with main partner
mod_condom_main <- geeglm(always_condom_main ~ visit_num, data = gee_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable') 
summary(mod_condom_main)
cc_condom_main <- coef(summary(mod_condom_main))
citab_condom_main <- with(as.data.frame(cc_condom_main),
                             cbind(est = exp(Estimate),
                                   lwr=exp(Estimate-1.96*Std.err),
                                   upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_main) <- rownames(cc_condom_main)
## OR: 1.01, 95% CI: [0.98, 1.04]
## p-value: 0.49


## Yc detection
mod_yc <- geeglm(yc_result ~ visit_num, data = gee_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable') 
summary(mod_yc)
cc_yc <- coef(summary(mod_yc))
citab_yc <- with(as.data.frame(cc_yc),
              cbind(est = exp(Estimate),
                    lwr=exp(Estimate-1.96*Std.err),
                    upr=exp(Estimate+1.96*Std.err)))
rownames(citab_yc) <- rownames(cc_yc)
## OR = 0.966, 95% confidence interval: [0.89, 1.05]
## p-value for trend = 0.4

############################################################################################
## Predictors of Yc detection
############################################################################################
## Baseline predictors
pred_dat <- base_dat_yc
pred_dat$id <- pred_dat$h00IdNumber
pred_dat <- left_join(yc_dat, pred_dat[, c("id", "age_cat", "registered", "ethnic_group_cat", "edu")], by = "id")
pred_dat$yc_result <- ifelse(pred_dat$value == "P", 1, 0)

## Age group
pred_dat %>%
  group_by(age_cat, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_age <- geeglm(yc_result ~ age_cat, data = pred_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable') 
summary(mod_age)
cc_age <- coef(summary(mod_age))
citab_age <- with(as.data.frame(cc_age),
                 cbind(est = exp(Estimate),
                       lwr=exp(Estimate-1.96*Std.err),
                       upr=exp(Estimate+1.96*Std.err)))
rownames(citab_age) <- rownames(cc_age)

## Registered
pred_dat %>%
  group_by(registered, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_reg <- geeglm(yc_result ~ registered, data = pred_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable')
summary(mod_reg)
cc_reg <- coef(summary(mod_reg))
citab_reg <- with(as.data.frame(cc_reg),
                  cbind(est = exp(Estimate),
                        lwr=exp(Estimate-1.96*Std.err),
                        upr=exp(Estimate+1.96*Std.err)))
rownames(citab_reg) <- rownames(cc_reg)

## Ethnic Group
pred_dat %>%
  group_by(ethnic_group_cat, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_eg <- geeglm(yc_result ~ ethnic_group_cat, data = na.omit(pred_dat[, c("id", "yc_result", "ethnic_group_cat")]), family = binomial(link = "logit"), id = id, corstr = 'exchangeable')
summary(mod_eg)
cc_eg <- coef(summary(mod_eg))
citab_eg <- with(as.data.frame(cc_eg),
                  cbind(est = exp(Estimate),
                        lwr=exp(Estimate-1.96*Std.err),
                        upr=exp(Estimate+1.96*Std.err)))
rownames(citab_eg) <- rownames(cc_eg)

## Education
pred_dat %>%
  group_by(edu, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_edu <- geeglm(yc_result ~ edu, data = pred_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable') 
summary(mod_edu)
cc_edu <- coef(summary(mod_edu))
citab_edu <- with(as.data.frame(cc_edu),
                  cbind(est = exp(Estimate),
                        lwr=exp(Estimate-1.96*Std.err),
                        upr=exp(Estimate+1.96*Std.err)))
rownames(citab_edu) <- rownames(cc_edu)

## Site
pred_dat %>%
  group_by(site, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_site <- geeglm(yc_result ~ site, data = pred_dat, family = binomial(link = "logit"), id = id, corstr = 'exchangeable') 
summary(mod_site)
cc_site <- coef(summary(mod_site))
citab_site <- with(as.data.frame(cc_site),
                  cbind(est = exp(Estimate),
                        lwr=exp(Estimate-1.96*Std.err),
                        upr=exp(Estimate+1.96*Std.err)))
rownames(citab_site) <- rownames(cc_site)

################################################################################################
## Sexual behavior predictors

## Sex with client in last week
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(reports_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_client <- geeglm(yc_result ~ reports_client, data = na.omit(sex_dat[, c("id", "yc_result", "reports_client")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_client)
cc_client <- coef(summary(mod_client))
citab_client <- with(as.data.frame(cc_client),
                   cbind(est = exp(Estimate),
                         lwr=exp(Estimate-1.96*Std.err),
                         upr=exp(Estimate+1.96*Std.err)))
rownames(citab_client) <- rownames(cc_client)
## OR = 0.999, 95% CI: [0.416, 2.401], p = 0.998

## Sex with main partner in last month
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(reports_main_partner, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_main <- geeglm(yc_result ~ reports_main_partner, data = na.omit(sex_dat[, c("id", "yc_result", "reports_main_partner")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_main)
cc_main <- coef(summary(mod_main))
citab_main <- with(as.data.frame(cc_main),
                     cbind(est = exp(Estimate),
                           lwr=exp(Estimate-1.96*Std.err),
                           upr=exp(Estimate+1.96*Std.err)))
rownames(citab_main) <- rownames(cc_main)
## OR = 2.104, 95% CI: [0.836, 5.299], p = 0.11

## Both
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(reports_either, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

mod_main_client <- geeglm(yc_result ~ reports_either, data = na.omit(sex_dat[, c("id", "yc_result", "reports_either")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_main_client)
cc_main_client <- coef(summary(mod_main_client))
citab_main_client <- with(as.data.frame(cc_main_client),
                   cbind(est = exp(Estimate),
                         lwr=exp(Estimate-1.96*Std.err),
                         upr=exp(Estimate+1.96*Std.err)))
rownames(citab_main_client) <- rownames(cc_main_client)
## OR = 0.807, 95% CI: [0.236, 2.76], p = 0.733

## Number of clients in last week
mod_num_clients <- geeglm(yc_result ~ NbreClieDernSemaine, data = na.omit(sex_dat[, c("id", "yc_result", "NbreClieDernSemaine")]), family = gaussian, id = id, corstr = "exchangeable")

summary(mod_num_clients)
cc_num_clients <- coef(summary(mod_num_clients))
citab_num_clients <- with(as.data.frame(cc_num_clients),
                   cbind(est = exp(Estimate),
                         lwr=exp(Estimate-1.96*Std.err),
                         upr=exp(Estimate+1.96*Std.err)))
rownames(citab_num_clients) <- rownames(cc_num_clients)
## OR = 1.01, 95% CI: [0.994, 1.02], p = 0.31

## Condom use with clients
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(always_condom_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

# 
# mod_condom_client <- geeglm(yc_result ~ always_condom_client, data = na.omit(sex_dat[, c("id", "yc_result", "always_condom_client")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")
mod_condom_client <- geeglm(yc_result ~ not_always_condom_client, data = na.omit(sex_dat[, c("id", "yc_result", "not_always_condom_client")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_condom_client)
cc_condom_client <- coef(summary(mod_condom_client))
citab_condom_client <- with(as.data.frame(cc_condom_client),
                   cbind(est = exp(Estimate),
                         lwr=exp(Estimate-1.96*Std.err),
                         upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_client) <- rownames(cc_condom_client)
## OR = 1.337, 95% CI: [0.1696, 10.54], p = 0.78
## OR = 0.748, 95% CI: [0.0948, 5.896], p = 0.78

## Condom use with main partner
sex_dat %>%
  filter(!is.na(value)) %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

# mod_condom_main <- geeglm(yc_result ~ always_condom_main, data = na.omit(sex_dat[sex_dat$NbrePartMoisDernier > 0, c("id", "yc_result", "always_condom_main")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")
mod_condom_main <- geeglm(yc_result ~ not_always_condom_main, data = na.omit(sex_dat[sex_dat$NbrePartMoisDernier > 0, c("id", "yc_result", "not_always_condom_main")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_condom_main)
cc_condom_main <- coef(summary(mod_condom_main))
citab_condom_main <- with(as.data.frame(cc_condom_main),
                            cbind(est = exp(Estimate),
                                  lwr=exp(Estimate-1.96*Std.err),
                                  upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_main) <- rownames(cc_condom_main)
## OR = 1.28, 95% CI: [0.411, 3.96], p = 0.67
## OR = 0.78, 95% CI: [0.25, 2.43], p = 0.67

## Both
sex_dat %>%
  filter(!is.na(value)) %>%
  filter(NbrePartMoisDernier > 0) %>%
  group_by(always_condom_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

# mod_condom_both <- geeglm(yc_result ~ always_condom_both, data = na.omit(sex_dat[sex_dat$NbrePartMoisDernier > 0, c("id", "yc_result", "always_condom_both")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")
mod_condom_both <- geeglm(yc_result ~ not_always_condom_both, data = na.omit(sex_dat[sex_dat$NbrePartMoisDernier > 0, c("id", "yc_result", "not_always_condom_both")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_condom_both)
cc_condom_both <- coef(summary(mod_condom_both))
citab_condom_both <- with(as.data.frame(cc_condom_both),
                          cbind(est = exp(Estimate),
                                lwr=exp(Estimate-1.96*Std.err),
                                upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_both) <- rownames(cc_condom_both)
## OR = 0.9600, 95% CI: [0.350, 2.637], p = 0.71
## OR = 0.80, 95% CI: [0.25, 2.55], p = 0.71

## Condom use at most recent sex with client
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(used_condom_recent_client, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

# mod_condom_client_recent <- geeglm(yc_result ~ used_condom_recent_client, data = na.omit(sex_dat[, c("id", "yc_result", "used_condom_recent_client")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")
mod_condom_client_recent <- geeglm(yc_result ~ not_used_condom_recent_client, data = na.omit(sex_dat[, c("id", "yc_result", "not_used_condom_recent_client")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_condom_client_recent)
cc_condom_client_recent <- coef(summary(mod_condom_client_recent))
citab_condom_client_recent <- with(as.data.frame(cc_condom_client_recent),
                            cbind(est = exp(Estimate),
                                  lwr=exp(Estimate-1.96*Std.err),
                                  upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_client_recent) <- rownames(cc_condom_client_recent)
## OR = 0.571, 95% CI: [0.0499, 6.53], p = 0.65
## OR = 1.75, 95% CI: [0.15, 20.03], p = 0.65

## Condom use at most recent sex with main partner
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(used_condom_recent_main, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

# mod_condom_main_recent <- geeglm(yc_result ~ used_condom_recent_main, data = na.omit(sex_dat[, c("id", "yc_result", "used_condom_recent_main")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")
mod_condom_main_recent <- geeglm(yc_result ~ not_used_condom_recent_main, data = na.omit(sex_dat[, c("id", "yc_result", "not_used_condom_recent_main")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_condom_main_recent)
cc_condom_main_recent <- coef(summary(mod_condom_main_recent))
citab_condom_main_recent <- with(as.data.frame(cc_condom_main_recent),
                                   cbind(est = exp(Estimate),
                                         lwr=exp(Estimate-1.96*Std.err),
                                         upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_main_recent) <- rownames(cc_condom_main_recent)
## OR = 0.73, 95% CI: [0.28, 1.90], p = 0.51
## OR = 1.38, 95% CI: [0.53, 3.59], p = 0.51

## Condom use at most recent sex with both
sex_dat %>%
  filter(!is.na(value)) %>%
  group_by(used_condom_recent_both, value) %>%
  summarise(n = n()) %>%
  mutate(pct = n/sum(n))

# mod_condom_both_recent <- geeglm(yc_result ~ used_condom_recent_both, data = na.omit(sex_dat[, c("id", "yc_result", "used_condom_recent_both")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")
mod_condom_both_recent <- geeglm(yc_result ~ not_used_condom_recent_both, data = na.omit(sex_dat[, c("id", "yc_result", "not_used_condom_recent_both")]), family = binomial(link = "logit"), id = id, corstr = "exchangeable")

summary(mod_condom_both_recent)
cc_condom_both_recent <- coef(summary(mod_condom_both_recent))
citab_condom_both_recent <- with(as.data.frame(cc_condom_both_recent),
                                 cbind(est = exp(Estimate),
                                       lwr=exp(Estimate-1.96*Std.err),
                                       upr=exp(Estimate+1.96*Std.err)))
rownames(citab_condom_both_recent) <- rownames(cc_condom_both_recent)
## OR = 1.27, 95% CI: [0.50, 3.21], p = 0.62

