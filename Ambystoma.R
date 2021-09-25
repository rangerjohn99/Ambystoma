
## Load in data ##
# Read in csv file of data---------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/rangerjohn99/Ambystoma/main/Ambystoma_final_v2.csv")
Ambystoma_final <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(Ambystoma_final)

# Make variables factors---------------------------------
Ambystoma_final$larval = factor(Ambystoma_final$larval, levels = c("Adult", "Larval"))
levels(Ambystoma_final$larval)
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
Ambystoma_final$centrum_trunk = factor(Ambystoma_final$centrum_trunk, levels = c("Closed", "Open"))
levels(Ambystoma_final$centrum_trunk)

Ambystoma_final$centrum_sacral = factor(Ambystoma_final$centrum_sacral, levels = c("Closed", "Open"))
levels(Ambystoma_final$centrum_sacral)

# 8th centrum model

centrum_model <- glm(centrum_trunk~log(SVL_P)+larval+Sex+species,data=Ambystoma_final, family = 'binomial')
summary(centrum_model)
exp(centrum_model$coefficients)


# sacral centrum model

sacral_model <- glm(centrum_sacral~log(SVL_P)+larval+Sex+species,data=Ambystoma_final, family = 'binomial')
summary(sacral_model)
exp(sacral_model$coefficients)

# atlas model

atlas_model <- glm(atlas_fused~log(SVL_P)+larval+Sex+species,data=Ambystoma_final,family = 'binomial')
summary(atlas_model)
exp(atlas_model$coefficients)

# SVL_P MF ANOVA

svl_model <- aov(log(SVL_P)~Sex+species,data=Ambystoma_final)
check_model(svl_model)

svlaov <- aov(log(SVL_P)~Sex+species,data=Ambystoma_final,contrasts=list(Sex=contr.sum)) 

Anova(svlaov, type=3)

summary.lm(svlaov)$adj.r.squared

# test individual species for sexual dimorphism (n=10)

# maculatum

maculatum <- subset(Ambystoma_final, species == "Ambystoma_maculatum", select = c("SVL_P", "Sex"))

svl_model_maculatum <- aov(log(SVL_P)~Sex,data=maculatum)
check_model(svl_model_maculatum)

svlaov_maculatum <- aov(log(SVL_P)~Sex,data=maculatum,contrasts=list(Sex=contr.sum)) 

Anova(svlaov_maculatum, type=3)

summary.lm(svlaov_maculatum)$adj.r.square

# opacum

opacum <- subset(Ambystoma_final, species == "Ambystoma_opacum", select = c("SVL_P", "Sex"))

svl_model_opacum <- aov(log(SVL_P)~Sex,data=opacum)
check_model(svl_model_opacum)

svlaov_opacum <- aov(log(SVL_P)~Sex,data=opacum,contrasts=list(Sex=contr.sum)) 

Anova(svlaov_opacum, type=3)

summary.lm(svlaov_opacum)$adj.r.square

# tigrinum

tigrinum <- subset(Ambystoma_final, species == "Ambystoma_tigrinum", select = c("SVL_P", "Sex"))

svl_model_tigrinum <- aov(log(SVL_P)~Sex,data=tigrinum)
check_model(svl_model_tigrinum)

svlaov_tigrinum <- aov(log(SVL_P)~Sex,data=tigrinum,contrasts=list(Sex=contr.sum)) 

Anova(svlaov_tigrinum, type=3)

summary.lm(svlaov_tigrinum)$adj.r.square
