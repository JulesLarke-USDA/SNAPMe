# Author: Jules Larke
# Date: 041223
# Purpose: Create linear models to examine correlations in snapme and bitesnap nutrient vectors across covariates
# covariates: age, BMI and education and Caucasian/Non-Caucasian.

library(ggplot2)
library(dplyr)
library(boot)

dat <- read.csv('../input/snapme_bitesnap_nutrients.csv')

cal <- boot(dat,function(data,indices)
  summary(lm(KCAL ~ calories.cal. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
cal$t0

fit = summary(lm(KCAL ~ calories.cal. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
cal_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

cal_df = data.frame(cal$t0, t(quantile(cal$t,c(0.025,0.975))), pval, 'Energy (Kcal)')
colnames(cal_df)[1] = 'R2'
colnames(cal_df)[2] = 'lower'
colnames(cal_df)[3] = 'upper'
colnames(cal_df)[4] = 'pval'
colnames(cal_df)[5] = 'name'

# protein
pro <- boot(dat,function(data,indices)
  summary(lm(PROT ~ protein.g. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
pro$t0

fit = summary(lm(PROT ~ protein.g. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
pro_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

pro_df = data.frame(pro$t0, t(quantile(pro$t,c(0.025,0.975))), pval, 'Protein (g)')
colnames(pro_df)[1] = 'R2'
colnames(pro_df)[2] = 'lower'
colnames(pro_df)[3] = 'upper'
colnames(pro_df)[4] = 'pval'
colnames(pro_df)[5] = 'name'

?boot
# carbs
carb <- boot(dat,function(data,indices)
  summary(lm(log(CARB) ~ log(totalCarb.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
carb$t0

fit = summary(lm(log(CARB) ~ log(totalCarb.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
carb_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

carb_df = data.frame(carb$t0, t(quantile(carb$t,c(0.025,0.975))), pval, 'Carbohydrate (g)')
colnames(carb_df)[1] = 'R2'
colnames(carb_df)[2] = 'lower'
colnames(carb_df)[3] = 'upper'
colnames(carb_df)[4] = 'pval'
colnames(carb_df)[5] = 'name'

# total fat
tfat <- boot(dat,function(data,indices)
  summary(lm(TFAT ~ totalFat.g. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
tfat$t0

fit = summary(lm(TFAT ~ totalFat.g. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
tfat_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

tfat_df = data.frame(tfat$t0, t(quantile(tfat$t,c(0.025,0.975))), pval, 'Total Fat (g)')
colnames(tfat_df)[1] = 'R2'
colnames(tfat_df)[2] = 'lower'
colnames(tfat_df)[3] = 'upper'
colnames(tfat_df)[4] = 'pval'
colnames(tfat_df)[5] = 'name'

# sat fat
satfat <- boot(dat,function(data,indices)
  summary(lm(sqrt(SFAT) ~ sqrt(saturatedFat.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
satfat$t0

fit = summary(lm(sqrt(SFAT) ~ sqrt(saturatedFat.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)

qqnorm(residuals(fit),
       ylab="Sample Quantiles for residuals")
qqline(residuals(fit),
       col="red")

hist(fit$residuals)

pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
sfat_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

satfat_df = data.frame(satfat$t0, t(quantile(satfat$t,c(0.025,0.975))), pval, 'Saturated Fat (g)')
colnames(satfat_df)[1] = 'R2'
colnames(satfat_df)[2] = 'lower'
colnames(satfat_df)[3] = 'upper'
colnames(satfat_df)[4] = 'pval'
colnames(satfat_df)[5] = 'name'

# MUFA
mufa <- boot(dat,function(data,indices)
  summary(lm(sqrt(MFAT) ~ sqrt(monounsaturatedFat.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
mufa$t0

fit = summary(lm(sqrt(MFAT) ~ sqrt(monounsaturatedFat.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)

qqnorm(residuals(fit),
       ylab="Sample Quantiles for residuals")
qqline(residuals(fit),
       col="red")

hist(fit$residuals)

pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
mufa_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

mufa_df = data.frame(mufa$t0, t(quantile(mufa$t,c(0.025,0.975))), pval, 'Monounsaturated Fat (g)')
colnames(mufa_df)[1] = 'R2'
colnames(mufa_df)[2] = 'lower'
colnames(mufa_df)[3] = 'upper'
colnames(mufa_df)[4] = 'pval'
colnames(mufa_df)[5] = 'name'

# PUFA
pufa <- boot(dat,function(data,indices)
  summary(lm(sqrt(PFAT) ~ sqrt(polyunsaturatedFat.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
pufa$t0

fit = summary(lm(sqrt(PFAT) ~ sqrt(polyunsaturatedFat.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)

qqnorm(residuals(fit),
       ylab="Sample Quantiles for residuals")
qqline(residuals(fit),
       col="red")

hist(fit$residuals)

pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
pufa_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

pufa_df = data.frame(pufa$t0, t(quantile(pufa$t,c(0.025,0.975))), pval, 'Polyunsaturated Fat (g)')
colnames(pufa_df)[1] = 'R2'
colnames(pufa_df)[2] = 'lower'
colnames(pufa_df)[3] = 'upper'
colnames(pufa_df)[4] = 'pval'
colnames(pufa_df)[5] = 'name'

# water
water <- boot(dat,function(data,indices)
  summary(lm(log(dat$MOIS + sqrt(dat$MOIS^2 + 1)) ~ log(dat$water.ml. + sqrt(dat$water.ml.^2 + 1)) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
water$t0

dat$trans_MOIS = log(dat$MOIS + sqrt(dat$MOIS^2 + 1))
dat$trans_water = log(dat$water.ml. + sqrt(dat$water.ml.^2 + 1))

fit = summary(lm(trans_MOIS ~ trans_water + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)

qqnorm(residuals(fit),
       ylab="Sample Quantiles for residuals")
qqline(residuals(fit),
       col="red")

hist(fit$residuals)

pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
water_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

water_df = data.frame(water$t0, t(quantile(water$t,c(0.025,0.975))), pval, 'Water (mL)')
colnames(water_df)[1] = 'R2'
colnames(water_df)[2] = 'lower'
colnames(water_df)[3] = 'upper'
colnames(water_df)[4] = 'pval'
colnames(water_df)[5] = 'name'
# generalized log transformation

# Sugar
sugar <- boot(dat,function(data,indices)
  summary(lm(sqrt(SUGR) ~ sqrt(sugars.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
sugar$t0

fit = summary(lm(sqrt(SUGR) ~ sqrt(sugars.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
sugar_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

sugar_df = data.frame(sugar$t0, t(quantile(sugar$t,c(0.025,0.975))), pval, 'Sugar (g)')
colnames(sugar_df)[1] = 'R2'
colnames(sugar_df)[2] = 'lower'
colnames(sugar_df)[3] = 'upper'
colnames(sugar_df)[4] = 'pval'
colnames(sugar_df)[5] = 'name'


# fiber
fiber <- boot(dat,function(data,indices)
  summary(lm(log(FIBE) ~ log(dietaryFiber.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
fiber$t0

fit = summary(lm(log(FIBE) ~ log(dietaryFiber.g.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
fiber_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

fiber_df = data.frame(fiber$t0, t(quantile(fiber$t,c(0.025,0.975))), pval, 'Fiber, total dietary (g)')
colnames(fiber_df)[1] = 'R2'
colnames(fiber_df)[2] = 'lower'
colnames(fiber_df)[3] = 'upper'
colnames(fiber_df)[4] = 'pval'
colnames(fiber_df)[5] = 'name'

# Iron
iron <- boot(dat,function(data,indices)
  summary(lm(sqrt(IRON) ~ sqrt(iron.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
iron$t0

fit = summary(lm(sqrt(IRON) ~ sqrt(iron.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
iron_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

iron_df = data.frame(iron$t0, t(quantile(iron$t,c(0.025,0.975))), pval, 'Iron (mg)')
colnames(iron_df)[1] = 'R2'
colnames(iron_df)[2] = 'lower'
colnames(iron_df)[3] = 'upper'
colnames(iron_df)[4] = 'pval'
colnames(iron_df)[5] = 'name'

# magnesium
magn <- boot(dat,function(data,indices)
  summary(lm(log(MAGN) ~ log(magnesium.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
magn$t0

fit = summary(lm(log(MAGN) ~ log(magnesium.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
magn_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

magn_df = data.frame(magn$t0, t(quantile(magn$t,c(0.025,0.975))), pval, 'Magnesium (mg)')
colnames(magn_df)[1] = 'R2'
colnames(magn_df)[2] = 'lower'
colnames(magn_df)[3] = 'upper'
colnames(magn_df)[4] = 'pval'
colnames(magn_df)[5] = 'name'

# sodium
sodi <- boot(dat,function(data,indices)
  summary(lm(sqrt(SODI) ~ sqrt(sodium.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
sodi$t0

fit = summary(lm(sqrt(SODI) ~ sqrt(sodium.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
sodi_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

sodi_df = data.frame(sodi$t0, t(quantile(sodi$t,c(0.025,0.975))), pval, 'Sodium (mg)')
colnames(sodi_df)[1] = 'R2'
colnames(sodi_df)[2] = 'lower'
colnames(sodi_df)[3] = 'upper'
colnames(sodi_df)[4] = 'pval'
colnames(sodi_df)[5] = 'name'
# sqrt transform

# zinc
zinc <- boot(dat,function(data,indices)
  summary(lm(sqrt(ZINC) ~ sqrt(zinc.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
zinc$t0

fit = summary(lm(sqrt(ZINC) ~ sqrt(zinc.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
zinc_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

zinc_df = data.frame(zinc$t0, t(quantile(zinc$t,c(0.025,0.975))), pval, 'Zinc (mg)')
colnames(zinc_df)[1] = 'R2'
colnames(zinc_df)[2] = 'lower'
colnames(zinc_df)[3] = 'upper'
colnames(zinc_df)[4] = 'pval'
colnames(zinc_df)[5] = 'name'
# sqrt transform

# copper
copp <- boot(dat,function(data,indices)
  summary(lm(log(COPP) ~ log(copper.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
copp$t0

fit = summary(lm(log(COPP) ~ log(copper.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
copp_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

copp_df = data.frame(copp$t0, t(quantile(copp$t,c(0.025,0.975))), pval, 'Copper (mg)')
colnames(copp_df)[1] = 'R2'
colnames(copp_df)[2] = 'lower'
colnames(copp_df)[3] = 'upper'
colnames(copp_df)[4] = 'pval'
colnames(copp_df)[5] = 'name'
# log transform

# Selenium
sele <- boot(dat,function(data,indices)
  summary(lm(sqrt(SELE) ~ sqrt(selenium.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
sele$t0

fit = summary(lm(sqrt(SELE) ~ sqrt(selenium.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
sele_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

sele_df = data.frame(sele$t0, t(quantile(sele$t,c(0.025,0.975))), pval, 'Selenium (mcg)')
colnames(sele_df)[1] = 'R2'
colnames(sele_df)[2] = 'lower'
colnames(sele_df)[3] = 'upper'
colnames(sele_df)[4] = 'pval'
colnames(sele_df)[5] = 'name'
# sqrt transform


# thiamin
thia <- boot(dat,function(data,indices)
  summary(lm(log(VB1) ~ log(thiamin.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
thia$t0

fit = summary(lm(log(VB1) ~ log(thiamin.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
thia_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

thia_df = data.frame(thia$t0, t(quantile(thia$t,c(0.025,0.975))), pval, 'Thiamin (mg)')
colnames(thia_df)[1] = 'R2'
colnames(thia_df)[2] = 'lower'
colnames(thia_df)[3] = 'upper'
colnames(thia_df)[4] = 'pval'
colnames(thia_df)[5] = 'name'
# log transform

# riboflavin
ribo <- boot(dat,function(data,indices)
  summary(lm(sqrt(VB2) ~ sqrt(riboflavin.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
ribo$t0

fit = summary(lm(sqrt(VB2) ~ sqrt(riboflavin.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
ribo_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

ribo_df = data.frame(ribo$t0, t(quantile(ribo$t,c(0.025,0.975))), pval, 'Riboflavin (mg)')
colnames(ribo_df)[1] = 'R2'
colnames(ribo_df)[2] = 'lower'
colnames(ribo_df)[3] = 'upper'
colnames(ribo_df)[4] = 'pval'
colnames(ribo_df)[5] = 'name'
# sqrt transform


# niacin
niac <- boot(dat,function(data,indices)
  summary(lm(NIAC ~ niacin.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
niac$t0

fit = summary(lm(NIAC ~ niacin.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
niac_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

niac_df = data.frame(niac$t0, t(quantile(niac$t,c(0.025,0.975))), pval, 'Niacin (mg)')
colnames(niac_df)[1] = 'R2'
colnames(niac_df)[2] = 'lower'
colnames(niac_df)[3] = 'upper'
colnames(niac_df)[4] = 'pval'
colnames(niac_df)[5] = 'name'

# vitamin B6
b6 <- boot(dat,function(data,indices)
  summary(lm(VB6 ~ vitaminB6.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
b6$t0

fit = summary(lm(VB6 ~ vitaminB6.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
b6_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

b6_df = data.frame(b6$t0, t(quantile(b6$t,c(0.025,0.975))), pval, 'Vitamin B6 (mg)')
colnames(b6_df)[1] = 'R2'
colnames(b6_df)[2] = 'lower'
colnames(b6_df)[3] = 'upper'
colnames(b6_df)[4] = 'pval'
colnames(b6_df)[5] = 'name'

# folate
folate <- boot(dat,function(data,indices)
  summary(lm(FF ~ folate.mcg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
folate$t0

fit = summary(lm(FF ~ folate.mcg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
folate_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

folate_df = data.frame(folate$t0, t(quantile(folate$t,c(0.025,0.975))), pval, 'Folate, food (mcg)')
colnames(folate_df)[1] = 'R2'
colnames(folate_df)[2] = 'lower'
colnames(folate_df)[3] = 'upper'
colnames(folate_df)[4] = 'pval'
colnames(folate_df)[5] = 'name'

# B12
b12 <- boot(dat,function(data,indices)
  summary(lm(sqrt(VB12) ~ sqrt(vitaminB12.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
b12$t0

fit = summary(lm(sqrt(VB12) ~ sqrt(vitaminB12.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
b12_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

b12_df = data.frame(b12$t0, t(quantile(b12$t,c(0.025,0.975))), pval, 'Vitamin B12 (mcg)')
colnames(b12_df)[1] = 'R2'
colnames(b12_df)[2] = 'lower'
colnames(b12_df)[3] = 'upper'
colnames(b12_df)[4] = 'pval'
colnames(b12_df)[5] = 'name'

# vitamin A
vita <- boot(dat,function(data,indices)
  summary(lm(log(VARA) ~ log(vitaminA.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
vita$t0

fit = summary(lm(log(VARA) ~ log(vitaminA.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
vita_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

vita_df = data.frame(vita$t0, t(quantile(vita$t,c(0.025,0.975))), pval, 'Vitamin A, (mcg [RAE])')
colnames(vita_df)[1] = 'R2'
colnames(vita_df)[2] = 'lower'
colnames(vita_df)[3] = 'upper'
colnames(vita_df)[4] = 'pval'
colnames(vita_df)[5] = 'name'

# beta cryptoxanthin
cryp <- boot(dat,function(data,indices)
  summary(lm(log(dat$CRYP + sqrt(dat$CRYP^2 + 1)) ~ log(dat$betaCryptoxanthin.mcg. + sqrt(dat$betaCryptoxanthin.mcg.^2 + 1)) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
cryp$t0


dat$trans_CRYP = log(dat$CRYP + sqrt(dat$CRYP^2 + 1))
dat$trans_betaCryptoxanthin = log(dat$betaCryptoxanthin.mcg. + sqrt(dat$betaCryptoxanthin.mcg.^2 + 1))

fit = summary(lm(trans_CRYP ~ trans_betaCryptoxanthin + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
cryp_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

cryp_df = data.frame(cryp$t0, t(quantile(cryp$t,c(0.025,0.975))), pval, 'Cryptoxanthin, beta (mcg)')
colnames(cryp_df)[1] = 'R2'
colnames(cryp_df)[2] = 'lower'
colnames(cryp_df)[3] = 'upper'
colnames(cryp_df)[4] = 'pval'
colnames(cryp_df)[5] = 'name'


# lycopene
dat$trans_LYCO = dat$LYCO^(1/3) # cube root xform
dat$trans_lycopene = dat$lycopene.mcg.^(1/3)

lyco <- boot(dat,function(data,indices)
  summary(lm(trans_LYCO ~ trans_lycopene + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
lyco$t0

fit = summary(lm(trans_LYCO ~ trans_lycopene + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
lyco_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

lyco_df = data.frame(lyco$t0, t(quantile(lyco$t,c(0.025,0.975))), pval, 'Lycopene (mcg)')
colnames(lyco_df)[1] = 'R2'
colnames(lyco_df)[2] = 'lower'
colnames(lyco_df)[3] = 'upper'
colnames(lyco_df)[4] = 'pval'
colnames(lyco_df)[5] = 'name'

# vitamin E
vite <- boot(dat,function(data,indices)
  summary(lm(log(ATOC) ~ log(vitaminE.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
vite$t0

fit = summary(lm(log(ATOC) ~ log(vitaminE.mg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
vite_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

vite_df = data.frame(vite$t0, t(quantile(vite$t,c(0.025,0.975))), pval, 'Vitamin E, alpha-tocopherol (mg)')
colnames(vite_df)[1] = 'R2'
colnames(vite_df)[2] = 'lower'
colnames(vite_df)[3] = 'upper'
colnames(vite_df)[4] = 'pval'
colnames(vite_df)[5] = 'name'

# vitamin K
vk <- boot(dat,function(data,indices)
  summary(lm(log(VK) ~ log(vitaminK.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
vk$t0

fit = summary(lm(log(VK) ~ log(vitaminK.mcg.) + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
vk_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

vk_df = data.frame(vk$t0, t(quantile(vk$t,c(0.025,0.975))), pval, 'Vitamin K, phylloquinone (mcg)')
colnames(vk_df)[1] = 'R2'
colnames(vk_df)[2] = 'lower'
colnames(vk_df)[3] = 'upper'
colnames(vk_df)[4] = 'pval'
colnames(vk_df)[5] = 'name'


# Cholesterol
chole <- boot(dat,function(data,indices)
  summary(lm(CHOLE ~ cholesterol.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
chole$t0

fit = summary(lm(CHOLE ~ cholesterol.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
chole_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

chole_df = data.frame(chole$t0, t(quantile(chole$t,c(0.025,0.975))), pval, 'Cholesterol (mg)')
colnames(chole_df)[1] = 'R2'
colnames(chole_df)[2] = 'lower'
colnames(chole_df)[3] = 'upper'
colnames(chole_df)[4] = 'pval'
colnames(chole_df)[5] = 'name'

# Choline
choline <- boot(dat,function(data,indices)
  summary(lm(CHOLN ~ choline.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level,data[indices,]))$adj.r.squared,R=1000)
choline$t0

fit = summary(lm(CHOLN ~ choline.mg. + scale(Age, scale = F) + scale(bmi, scale = F) + ethnic_group + edu_level, data = dat))
shapiro.test(fit$residuals)
hist(fit$residuals)
pval = broom::glance(fit)$p.value
pval = p.adjust(pval, method = 'BH', n = 43)
choline_label = paste0("italic(p) ==", formatC(pval, format = "e", digits = 2))

choline_df = data.frame(choline$t0, t(quantile(choline$t,c(0.025,0.975))), pval, 'Choline (mg)')
colnames(choline_df)[1] = 'R2'
colnames(choline_df)[2] = 'lower'
colnames(choline_df)[3] = 'upper'
colnames(choline_df)[4] = 'pval'
colnames(choline_df)[5] = 'name'

nutrient_df <- rbind(b12_df, b6_df, cal_df, carb_df, chole_df, choline_df, copp_df, cryp_df, fiber_df, folate_df, iron_df, lyco_df, magn_df, mufa_df, niac_df, pro_df, pufa_df, ribo_df, satfat_df, sele_df, sodi_df, sugar_df, tfat_df, thia_df, vita_df, vite_df, vk_df, water_df, zinc_df)
nutrient_df <- nutrient_df %>% arrange(desc(R2))

nutrient_plot <- ggplot(nutrient_df, aes(x = R2, y = forcats::fct_reorder(name, R2))) +
  geom_point(size = 5) +
  geom_segment(aes(x= lower, xend = upper, yend = name), size = 1, linetype = 1, arrow = arrow(ends = "both", angle = 90, length = unit(0.05, "in"))) +
  ylab("Nutrient") +
  xlab("Adjusted R-squared") +
  theme_classic() +
  annotate(
    "text",
    x = 1.1,
    y = 1,
    label = folate_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 2,
    label = lyco_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 3,
    label = carb_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 4,
    label = vita_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 5,
    label = sodi_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 6,
    label = b6_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 7,
    label = copp_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 8,
    label = cal_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 9,
    label = vite_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 10,
    label = b12_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 11,
    label = thia_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 12,
    label = sele_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 13,
    label = sfat_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 14,
    label = magn_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 15,
    label = fiber_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 16,
    label = pufa_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 17,
    label = tfat_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 18,
    label = iron_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 19,
    label = water_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 20,
    label = mufa_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 21,
    label = zinc_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 22,
    label = pro_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 23,
    label = vk_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 24,
    label = ribo_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 25,
    label = niac_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 26,
    label = cryp_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 27,
    label = sugar_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 28,
    label = choline_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 29,
    label = chole_label,
    parse = T
  ) +
  coord_cartesian(xlim = c(0, 1), # This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  #geom_vline(xintercept= 0.5, linetype = "dashed") +
  theme(plot.margin = unit(c(1,5,1,1), "lines"),
    axis.text.x  = element_text(size = 10, color = 'black'),
    axis.text.y = element_text(size = 10, color = 'black'),
    axis.title.y = element_text(size = 10),
    axis.title.x = element_text(size = 10),
    legend.position = "none"
  )
nutrient_plot

## Nutrients that could not be fit to linear model as residuals were non-normally distributed after transformation
# Use partial spearman rank correlation

library(RVAideMemoire)
# Alcohol
dat_spread <- dat %>% mutate(value = 1) %>% tidyr::spread(edu_level, value,  fill = 0)
dat_spread <- dat_spread %>% mutate(value = 1) %>% tidyr::spread(ethnic_group, value,  fill = 0) 

alc = RVAideMemoire::pcor.test(dat_spread$ALC, dat_spread$alcohol.g., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
alc$p.value = p.adjust(alc$p.value, method = 'BH', n = 43)
alc_label = paste0("italic(p) ==", formatC(alc$p.value, format = "e", digits = 2))

alc_df = data.frame(alc$estimate, alc$conf.int[1], alc$conf.int[2], alc$p.value, 'Alcohol (g)')
colnames(alc_df)[1] = 'rho'
colnames(alc_df)[2] = 'lower'
colnames(alc_df)[3] = 'upper'
colnames(alc_df)[4] = 'pval'
colnames(alc_df)[5] = 'name'

# calcium
calc = RVAideMemoire::pcor.test(dat_spread$CALC, dat_spread$calcium.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
calc$p.value = p.adjust(calc$p.value, method = 'BH', n = 43)
calc_label = paste0("italic(p) ==", formatC(calc$p.value, format = "e", digits = 2))

calc_df = data.frame(calc$estimate, calc$conf.int[1], calc$conf.int[2], calc$p.value, 'Calcium (mg)')
colnames(calc_df)[1] = 'rho'
colnames(calc_df)[2] = 'lower'
colnames(calc_df)[3] = 'upper'
colnames(calc_df)[4] = 'pval'
colnames(calc_df)[5] = 'name'
# theobromine
theo = RVAideMemoire::pcor.test(dat_spread$THEO, dat_spread$theobromine.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
theo$p.value = p.adjust(theo$p.value, method = 'BH', n = 43)
theo_label = paste0("italic(p) ==", formatC(theo$p.value, format = "e", digits = 2))

theo_df = data.frame(theo$estimate, theo$conf.int[1], theo$conf.int[2], theo$p.value, 'Theobromine (mg)')
colnames(theo_df)[1] = 'rho'
colnames(theo_df)[2] = 'lower'
colnames(theo_df)[3] = 'upper'
colnames(theo_df)[4] = 'pval'
colnames(theo_df)[5] = 'name'

# caffeine
caff = RVAideMemoire::pcor.test(dat_spread$CAFF, dat_spread$caffeine.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
caff$p.value = p.adjust(caff$p.value, method = 'BH', n = 43)
caff_label = paste0("italic(p) ==", formatC(caff$p.value, format = "e", digits = 2))

caff_df = data.frame(caff$estimate, caff$conf.int[1], caff$conf.int[2], caff$p.value, 'Caffeine (mg)')
colnames(caff_df)[1] = 'rho'
colnames(caff_df)[2] = 'lower'
colnames(caff_df)[3] = 'upper'
colnames(caff_df)[4] = 'pval'
colnames(caff_df)[5] = 'name'

# phosphorus
phos = RVAideMemoire::pcor.test(dat_spread$PHOS, dat_spread$phosphorus.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
phos$p.value = p.adjust(phos$p.value, method = 'BH', n = 43)
phos_label = paste0("italic(p) ==", formatC(phos$p.value, format = "e", digits = 2))

phos_df = data.frame(phos$estimate, phos$conf.int[1], phos$conf.int[2], phos$p.value, 'Phosphorus (mg)')
colnames(phos_df)[1] = 'rho'
colnames(phos_df)[2] = 'lower'
colnames(phos_df)[3] = 'upper'
colnames(phos_df)[4] = 'pval'
colnames(phos_df)[5] = 'name'

# potassium
pota = RVAideMemoire::pcor.test(dat_spread$POTA, dat_spread$potassium.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
pota$p.value = p.adjust(pota$p.value, method = 'BH', n = 43)
pota_label = paste0("italic(p) ==", formatC(pota$p.value, format = "e", digits = 2))

pota_df = data.frame(pota$estimate, pota$conf.int[1], pota$conf.int[2], pota$p.value, 'Potassium (mg)')
colnames(pota_df)[1] = 'rho'
colnames(pota_df)[2] = 'lower'
colnames(pota_df)[3] = 'upper'
colnames(pota_df)[4] = 'pval'
colnames(pota_df)[5] = 'name'

# vit c
vitc = RVAideMemoire::pcor.test(dat_spread$VC, dat_spread$vitaminC.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
vitc$p.value = p.adjust(vitc$p.value, method = 'BH', n = 43)
vitc_label = paste0("italic(p) ==", formatC(vitc$p.value, format = "e", digits = 2))

vitc_df = data.frame(vitc$estimate, vitc$conf.int[1], vitc$conf.int[2], vitc$p.value, 'Vitamin C (mg)')
colnames(vitc_df)[1] = 'rho'
colnames(vitc_df)[2] = 'lower'
colnames(vitc_df)[3] = 'upper'
colnames(vitc_df)[4] = 'pval'
colnames(vitc_df)[5] = 'name'

# folic acid
folic = RVAideMemoire::pcor.test(dat_spread$FA, dat_spread$folicAcid.mcg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
folic$p.value = p.adjust(folic$p.value, method = 'BH', n = 43)
folic_label = paste0("italic(p) ==", formatC(folic$p.value, format = "e", digits = 2))

folic_df = data.frame(folic$estimate, folic$conf.int[1], folic$conf.int[2], folic$p.value, 'Folic acid (mcg)')
colnames(folic_df)[1] = 'rho'
colnames(folic_df)[2] = 'lower'
colnames(folic_df)[3] = 'upper'
colnames(folic_df)[4] = 'pval'
colnames(folic_df)[5] = 'name'

# Retinol
ret = RVAideMemoire::pcor.test(dat_spread$RET, dat_spread$retinol.mcg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
ret$p.value = p.adjust(ret$p.value, method = 'BH', n = 43)
ret_label = paste0("italic(p) ==", formatC(ret$p.value, format = "e", digits = 2))

ret_df = data.frame(ret$estimate, ret$conf.int[1], ret$conf.int[2], ret$p.value, 'Retinol (mcg)')
colnames(ret_df)[1] = 'rho'
colnames(ret_df)[2] = 'lower'
colnames(ret_df)[3] = 'upper'
colnames(ret_df)[4] = 'pval'
colnames(ret_df)[5] = 'name'

# beta Carotene
bcar = RVAideMemoire::pcor.test(dat_spread$BCAR, dat_spread$betaCarotene.mcg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
bcar$p.value = p.adjust(bcar$p.value, method = 'BH', n = 43)
bcar_label = paste0("italic(p) ==", formatC(bcar$p.value, format = "e", digits = 2))

bcar_df = data.frame(bcar$estimate, bcar$conf.int[1], bcar$conf.int[2], bcar$p.value, 'Beta-carotene (mcg)')
colnames(bcar_df)[1] = 'rho'
colnames(bcar_df)[2] = 'lower'
colnames(bcar_df)[3] = 'upper'
colnames(bcar_df)[4] = 'pval'
colnames(bcar_df)[5] = 'name'

# alpha carotene
acar = RVAideMemoire::pcor.test(dat_spread$ACAR, dat_spread$alphaCarotene.mcg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
acar$p.value = p.adjust(acar$p.value, method = 'BH', n = 43)
acar_label = paste0("italic(p) ==", formatC(acar$p.value, format = "e", digits = 2))

acar_df = data.frame(acar$estimate, acar$conf.int[1], acar$conf.int[2], acar$p.value, 'Alpha-carotene (mcg)')
colnames(acar_df)[1] = 'rho'
colnames(acar_df)[2] = 'lower'
colnames(acar_df)[3] = 'upper'
colnames(acar_df)[4] = 'pval'
colnames(acar_df)[5] = 'name'

# Lutein + zeaxanthin
lz = RVAideMemoire::pcor.test(dat_spread$LZ, dat_spread$luteinAndZeaxanthin.mcg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
lz$p.value = p.adjust(lz$p.value, method = 'BH', n = 43)
lz_label = paste0("italic(p) ==", formatC(lz$p.value, format = "e", digits = 2))

lz_df = data.frame(lz$estimate, lz$conf.int[1], lz$conf.int[2], lz$p.value, 'Lutein + zeaxanthin (mcg)')
colnames(lz_df)[1] = 'rho'
colnames(lz_df)[2] = 'lower'
colnames(lz_df)[3] = 'upper'
colnames(lz_df)[4] = 'pval'
colnames(lz_df)[5] = 'name'

# n-3 fatty acids Note: bitesnap estimate is actually in g (not mg which is in the variable name), so on the same scale as the ASA24 estimate
omega3 = RVAideMemoire::pcor.test(dat_spread$Omega.3.fatty.acids..g., dat_spread$omega3FattyAcids.mg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
omega3$p.value = p.adjust(omega3$p.value, method = 'BH', n = 43)
omega3_label = paste0("italic(p) ==", formatC(omega3$p.value, format = "e", digits = 2))

omega3_df = data.frame(omega3$estimate, omega3$conf.int[1], omega3$conf.int[2], omega3$p.value, 'Omega-3s (EPA + DPA + DHA) (g)')
colnames(omega3_df)[1] = 'rho'
colnames(omega3_df)[2] = 'lower'
colnames(omega3_df)[3] = 'upper'
colnames(omega3_df)[4] = 'pval'
colnames(omega3_df)[5] = 'name'

#Vitamin D
vitd = RVAideMemoire::pcor.test(dat_spread$VITD, dat_spread$vitaminD.mcg., dat_spread[,c("Age", "bmi", "Caucasian", "Highschool graduate", "MS, MD, DDS, JD, PhD, EdD", "Some college or AS")], conf.level = 0.95, nrep = 1000, method = "spearman")
vitd$p.value = p.adjust(vitd$p.value, method = 'BH', n = 43)
vitd_label = paste0("italic(p) ==", formatC(vitd$p.value, format = "e", digits = 2))

vitd_df = data.frame(vitd$estimate, vitd$conf.int[1], vitd$conf.int[2], vitd$p.value, 'Vitamin D (D2 + D3) (mcg)')
colnames(vitd_df)[1] = 'rho'
colnames(vitd_df)[2] = 'lower'
colnames(vitd_df)[3] = 'upper'
colnames(vitd_df)[4] = 'pval'
colnames(vitd_df)[5] = 'name'

nutrient_df2 <- rbind(alc_df, calc_df, theo_df, caff_df, phos_df, pota_df, vitc_df, folic_df, ret_df, bcar_df, acar_df, lz_df, omega3_df, vitd_df)
nutrient_df2 <- nutrient_df2 %>% arrange(desc(rho))

nutrient_plot2 <- ggplot(nutrient_df2, aes(x = rho, y = forcats::fct_reorder(name, rho))) +
  geom_point(size = 5) +
  geom_segment(aes(x= lower, xend = upper, yend = name), size = 1, linetype = 1, arrow = arrow(ends = "both", angle = 90, length = unit(0.05, "in"))) +
  ylab("Nutrient") +
  xlab("Spearman's rho") +
  theme_classic() +
  annotate(
    "text",
    x = 1.1,
    y = 1,
    label = calc_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 2,
    label = phos_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 3,
    label = folic_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 4,
    label = ret_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 5,
    label = pota_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 6,
    label = theo_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 7,
    label = acar_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 8,
    label = bcar_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 9,
    label = lz_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 10,
    label = omega3_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 11,
    label = vitd_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 12,
    label = vitc_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 13,
    label = alc_label,
    parse = T
  ) +
  annotate(
    "text",
    x = 1.1,
    y = 14,
    label = caff_label,
    parse = T
  ) +
  coord_cartesian(xlim = c(0, 1), # This focuses the x-axis on the range of interest
                  clip = 'off') +   # This keeps the labels from disappearing
  #geom_vline(xintercept= 0.5, linetype = "dashed") +
  theme(plot.margin = unit(c(1,5,1,1), "lines"),
        axis.text.x  = element_text(size = 10, color = 'black'),
        axis.text.y = element_text(size = 10, color = 'black'),
        axis.title.y = element_text(size = 10),
        axis.title.x = element_text(size = 10),
        legend.position = "none"
  )
nutrient_plot2

ggsave("figure4A.png",
       plot = nutrient_plot,
       width = 8.5,
       height = 5.5,
       dpi = 1000)

ggsave("figure4B.png",
       plot = nutrient_plot2,
       width = 8.5,
       height = 3,
       dpi = 1000)
