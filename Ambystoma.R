## Load in data ##
# Read in csv file of data---------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/rangerjohn99/Ambystoma/main/Ambystoma_final.csv")
Ambystoma_final <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(Ambystoma_final)

# Make variables factors---------------------------------
Ambystoma_final$larval <- as.factor(Ambystoma_final$larval)
Ambystoma_final$atlas_fused <- as.factor(Ambystoma_final$atlas_fused) 
Ambystoma_final$offset <- as.factor(Ambystoma_final$offset)

Ambystoma_final$Sex[Ambystoma_final$Sex=='']=NA
Ambystoma_final$Sex = droplevels(Ambystoma_final$Sex)
str(Ambystoma_final)

Ambystoma_final$species <-gsub("_UTEP.*","", Ambystoma_final$specimen) #makes species column
Ambystoma_final$species <- as.factor(Ambystoma_final$species)


# ANALYSES #---------------------------------
library(car)
library(performance)
library(olsrr)
library(afex)

# ANOVA ---------------------------------
# Check model assumptions
hist(Ambystoma_final$trunk, main="Histogram of Trunk Count", xlab="Number of Trunk verts")  

model <- lm(trunk~log(SVL_P)+larval+Sex+offset,data=Ambystoma_final)
ols_test_normality(model)
plot(model)
is_norm <- check_normality(model)
plot(is_norm)
plot(is_norm, type = "qq")
plot(is_norm, type = "qq", detrend = TRUE)

model <- aov(trunk~log(SVL_P)+larval+Sex+offset,data=Ambystoma_final)
check_model(model)

myaov <- aov(trunk~log(SVL_P)+larval+Sex+offset+species,data=Ambystoma_final,contrasts=list(larval=contr.sum,Sex=contr.sum,offset=contr.sum))  

Anova(myaov, type=3)

summary.lm(myaov)$adj.r.squared

# Ordinal Logistic Regression---------------------------------
library(MASS)
#Ordering the dependent variable
Ambystoma_final$trunk = factor(Ambystoma_final$trunk, levels = c("12", "13", "14", "15", "16", "17"), ordered = TRUE) 
Ambystoma_final$caudosacral = factor(Ambystoma_final$caudosacral, levels = c("1", "2", "3"), ordered = TRUE) 

#Exploratory data analysis
#Summarizing the data
summary(Ambystoma_final)
library(ggplot2)
ggplot(Ambystoma_final, aes(x = trunk, y = SVL_P, fill = Sex)) +facet_wrap(. ~ species, ncol = 3)+   geom_dotplot(binaxis = "y",stackdir = "center") +   theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# Making the models
model_fit <- polr(trunk~log(SVL_P)+larval+Sex+offset+species,data=Ambystoma_final, Hess = TRUE)
summary(model_fit)

model_fit_nosex <- polr(trunk~log(SVL_P)+larval+offset+Sex,data=Ambystoma_final, Hess = TRUE)
summary(model_fit_nosex)

# Comparing models
library(AICcmodavg)

models <- list(model_fit, model_fit_nosex)
model.names <- c('model_fit', 'model_fit_nosex')
aictab(cand.set = models, modnames = model.names)

# Significance of coefficients and intercepts
summary_table <- coef(summary(model_fit))
pval <- pnorm(abs(summary_table[, "t value"]),lower.tail = FALSE)* 2
summary_table <- cbind(summary_table, "p value" = round(pval,3))
summary_table

# centrum model

Ambystoma_final$centrum_trunk <- as.factor(Ambystoma_final$centrum_trunk)
Ambystoma_final$centrum_sacral <- as.factor(Ambystoma_final$centrum_sacral) 
Ambystoma_final$atlas_fused <- as.factor(Ambystoma_final$atlas_fused)

# 8th centrum model

centrum_model <- lm(centrum_trunk~log(SVL_P)+larval+Sex+species,data=Ambystoma_final)
ols_test_normality(centrum_model)
plot(centrum_model)
iscentrum8_norm <- check_normality(centrum_model)
plot(iscentrum8_norm)
plot(iscentrum8_norm, type = "qq")
plot(iscentrum8_norm, type = "qq", detrend = TRUE)

centrum_model <- aov(centrum_trunk~log(SVL_P)+larval+Sex+species,data=Ambystoma_final)
check_model(centrum_model)

centrumaov <- aov(centrum_trunk~log(SVL_P)+larval+Sex++species,data=Ambystoma_final,contrasts=list(larval=contr.sum,Sex=contr.sum))  

Anova(centrumaov, type=3)

summary.lm(centrumaov)$adj.r.squared

# caudosacral centrum model

caudosacral_model <- lm(centrum_sacral~log(SVL_P)+larval+Sex+species,data=Ambystoma_final)
ols_test_normality(caudosacral_model)
plot(caudosacral_model)
iscaudosacral_norm <- check_normality(caudosacral_model)
plot(iscaudosacral_norm)
plot(iscaudosacral_norm, type = "qq")
plot(iscaudosacral_norm, type = "qq", detrend = TRUE)

caudosacral_model <- aov(centrum_trunk~log(SVL_P)+larval+Sex+species,data=Ambystoma_final)
check_model(caudosacral_model)

caudosacralaov <- aov(centrum_sacral~log(SVL_P)+larval+Sex+offset+species,data=Ambystoma_final,contrasts=list(larval=contr.sum,Sex=contr.sum))  

Anova(caudosacralaov, type=3)

summary.lm(caudosacralaov)$adj.r.squared

# atlas model

atlas_model <- lm(atlas_fused~log(SVL_P)+larval+Sex+species,data=Ambystoma_final)
ols_test_normality(atlas_model)
plot(atlas_model)
isatlas_norm <- check_normality(atlas_model)
plot(isatlas_norm)
plot(isatlas_norm, type = "qq")
plot(isatlas_norm, type = "qq", detrend = TRUE)

atlas_model <- aov(atlas_fused~log(SVL_P)+larval+Sex+species,data=Ambystoma_final)
check_model(atlas_model)

atlasaov <- aov(atlas_fused~log(SVL_P)+larval+Sex+species,data=Ambystoma_final,contrasts=list(larval=contr.sum,Sex=contr.sum))  

Anova(atlasaov, type=3)

summary.lm(atlasaov)$adj.r.squared

# SVL_P MF ANOVA

svl_model <- aov(log(SVL_P)~Sex+species,data=Ambystoma_final)
check_model(svl_model)

svlaov <- aov(log(SVL_P)~Sex+species,data=Ambystoma_final,contrasts=list(Sex=contr.sum)) 

Anova(svlaov, type=3)

summary.lm(svlaov)$adj.r.squared
