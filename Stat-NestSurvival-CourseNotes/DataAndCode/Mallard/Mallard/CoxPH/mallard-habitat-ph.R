# Analysis of mallard data illustrating proportional hazard model and compare to logistic exposure model
# How to fit a nest level categorical variable (habitat)

# 2019-05-01 CJS Initial code

# This is the mallard data that ships with RMark

library(AICcmodavg)
library(emmeans)
library(ggplot2)
library(plyr)
library(readxl)
library(survival)
library(survminer)

source(file.path("..","..","CoxPH-additional-functions.R"))

# The dataframe must contain the following fields with the following names
#
#    FirstFound: day the nest was first found
#    LastPresent: last day that a chick was present in the nest
#    LastChecked: last day the nest was checked
#    Fate: fate of the nest; 0=hatch an
#    Freq: number of nests with this data
#   
# Also contains the following fields
#   Robel    - 	Measurement of Robel pole of visibility of nest
#   PpnGrass - proportion of grassland cover on the 10.4 km2 study site t
#   AgeFound - Age of nest when found
#   AgeDay1  - Age of nest on day 1 of study (can be negative)
#   Habitat	- N=Native; P=Planted; W=Wetland; R=roadside right of way
#
# In this example, the multiple visits to a nest have been collapsed
# to a single record for each nest.
# In more complex examples, you may have multple records per nest
# as shown in the mallard example.
#

malldata <- readxl::read_excel(file.path("..","mallard.xlsx"), 
                               sheet="mallard")
head(malldata)

malldata2 <- expand.nest.data(malldata)

malldata2$Habitat <- factor(malldata2$Habitat)  # Convert categorical variable to be factors

# Fit a particular model
# This is a model with S varying by habitat type
mod.hab <-  glm(Survive~Habitat,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.hab)

# fit a null models
mod.null <-  glm(Survive~1,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.null)

AICcmodavg::aictab( list(mod.hab, mod.null))



###------------------------------------------------------------
### [proportional hazard model

malldata3 <- expand.nest.data.ph(malldata)
head(malldata3[,c("FirstFound","LastPresent","LastChecked","Fate","Start","End","Fail","Surv")], n=100)


malldata3$Habitat <- factor(malldata3$Habitat)  # Convert categorical variable to be factors

mod.hab.ph <-  coxph(Surv~Habitat,data=malldata3)
summary(mod.hab.ph)

car::Anova(mod.hab.ph, type=3)

mod.hab.ph.emmo <- emmeans::emmeans(mod.hab.ph, ~Habitat)
summary(pairs(mod.hab.ph.emmo), infer=TRUE)
summary(pairs(mod.hab.ph.emmo, type="response"), infer=TRUE)


# baseline cumulative hazard
cumhaz <- basehaz(mod.hab.ph)
# estimate change in cumulative hazard and plot
cumhaz$deltaHaz <- c(NA,diff(cumhaz$hazard))
basehaz.plot <- ggplot(data=cumhaz, aes(x=time, y=deltaHaz))+
  ggtitle("Estimated baseline hazard function over time")+
  geom_point()+
  geom_smooth(se=FALSE)
basehaz.plot
ggsave(basehaz.plot,
    file=file.path("..","..","..","..","MyStuff","Images","mallard-ph-basehazard.png"), h=4, w=6, units="in", dpi=300)


# test the assumption of proportionality
# look for approximate parallelism of the curveys
cox.zph(mod.hab.ph)
ggcoxzph(cox.zph(mod.hab.ph))

png(file=file.path("..","..","..","..","MyStuff","Images","mallard-ph-hab-testprop.png"), h=4, w=6, units="in", res=300)
ggcoxzph(cox.zph(mod.hab.ph))
dev.off()


# estimated survival curves for the two habitats
pred.data <- expand.grid( Habitat=unique(malldata3$Habitat))
pred.survival <- survfit(mod.hab.ph, newdata=pred.data, se.fit=TRUE)
plot(pred.survival)


mod.null.ph <- coxph(Surv~1, data=malldata3)


AICcmodavg::aictab( list(mod.hab.ph, mod.null.ph))
