
## Test for trends 
gee_dat <- yc_dat ## Sex dat
gee_dat$visit_num <- c(0, 1, 3, 6, 9, 12)[match(gee_dat$visit, c("scr", "m1", "m3", "m6", "m9", "m12"))]

gee_dat$yc_result <- ifelse(gee_dat$value == "P", 1, 0)

mod_glm <- glm(yc_result ~ visit_num, data = gee_dat, family = binomial(link = "log"))
summary(mod_glm)

mod_gee <- geeglm(yc_result ~ visit_num, data = gee_dat, family = binomial(link = "log"), id = id, corstr = 'exchangeable') 
summary(mod_gee)

# mod_ind_fe <- glm(yc_result ~ visit_num + factor(id), data = gee_dat, family = binomial(link = "log"))
# summary(mod_ind_fe)

# mod_robust <- glm(yc_result ~ visit_num, data = gee_dat, family = "binomial")
# coeftest(mod_robust, vcov = vcovHC(mod_robust, type="HC1"))

test <- na.omit(gee_dat[, c("yc_result", "reports_main_partner", "id")])
mod_gee <- geeglm(yc_result ~ reports_main_partner, data = test, family = binomial(link = "log"), id = id, corstr = 'exchangeable')
summary(mod_gee)
exp(coefficients(mod_gee))
