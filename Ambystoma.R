## Load in data ##
# Read in csv file of data---------------------------------
library(curl)
f2 <- curl("https://raw.githubusercontent.com/rangerjohn99/Ambystoma/main/Ambystoma_final.csv")
Ambystoma_final <- read.csv(f2, header = TRUE, sep = ",", stringsAsFactors = TRUE) # this is a matrix of each specimen  
head(Ambystoma_final)


Ambystoma_final$larval <- as.factor(Ambystoma_final$larval)
Ambystoma_final$atlas_fused <- as.factor(Ambystoma_final$atlas_fused) 
Ambystoma_final$offset <- as.factor(Ambystoma_final$offset)

Ambystoma_final$Sex[Ambystoma_final$Sex=='']=NA
Ambystoma_final$Sex = droplevels(Ambystoma_final$Sex)
str(Ambystoma_final)

Ambystoma_final$species <-gsub("_UTEP.*","", Ambystoma_final$specimen)

library(car)
library(lme4)

# need to deal with getting a species variable
# not sure how to deal with NA and missing values (why does sex not have NA values?)
model <- lmer(trunk~log(SVL_P)+larval+Sex+offset,data=Ambystoma_final)
plot(check_distribution(model))

myaov <- aov(trunk~log(SVL_P)+larval+Sex+offset,data=Ambystoma_final,contrasts=list(larval=contr.sum,Sex=contr.sum,offset=contr.sum))  

Anova(myaov, type=3)

summary.lm(myaov)$adj.r.squared
