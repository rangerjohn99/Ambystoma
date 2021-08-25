Ambystoma_final$larval <- as.factor(Ambystoma_final$larval)
Ambystoma_final$atlas_fused <- as.factor(Ambystoma_final$atlas_fused) 
Ambystoma_final$offset <- as.factor(Ambystoma_final$offset)

gsub("_UTEP.*","", Ambystoma_final$specimen)

library(car)

myaov <- aov(trunk~log(SVL_P)+larval+Sex+offset,data=Ambystoma_final,contrasts=list(larval=contr.sum,Sex=contr.sum,offset=contr.sum))  

Anova(myaov, type=3)
