# Analysis of mallard data using logistic exposure model vs proportional hazard model
# How to fit a nest level continous covariate (Robel height

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

malldata2$Habitat <- factor(malldata2$Habitat)  # RMark wants categorical variable to be factors



# Fit a particular model
# This is a model with S varying by robel height
mod.rob <- glm(Survive~Robel,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.rob)

# fit a null models
mod.null <-  glm(Survive~1,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.null)

AICcmodavg::aictab( list(mod.rob, mod.null))


#-------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------

malldata3 <- expand.nest.data.ph(malldata)


mod.rob.ph <-  coxph(Surv~Robel,data=malldata3)
summary(mod.rob.ph)


# test the assumption of proportionality
# look for approximate parallelism of the curveys
cox.zph(mod.rob.ph)
ggcoxzph(cox.zph(mod.rob.ph))

png(file=file.path("..","..","..","..","MyStuff","Images","mallard-ph-rob-testprop.png"), h=4, w=6, units="in", res=300)
ggcoxzph(cox.zph(mod.rob.ph))
dev.off()


# baseline cumulative hazard
cumhaz <- basehaz(mod.rob.ph)
# estimate change in cumulative hazard and plot
cumhaz$deltaHaz <- c(NA,diff(cumhaz$hazard))
ggplot(data=cumhaz, aes(x=time, y=deltaHaz))+
  ggtitle("Estimated baseline hazard function over time")+
  geom_point()+
  geom_smooth(se=FALSE)


# estimated hazard curves at three levels of the Robel
pred.data <- expand.grid( Robel=seq(min(malldata3$Robel),max(malldata3$Robel), length.out=3))
pred.survival <- survfit(mod.rob.ph, newdata=pred.data, se.fit=TRUE)
plot(pred.survival)



mod.null.ph <- coxph(Surv~1, data=malldata3)


AICcmodavg::aictab( list(mod.rob.ph, mod.null.ph))


