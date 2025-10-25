# Analysis of mallard data using logistic exposure model
# How to fit a nest level continous covariate (Robel height

# 2019-05-01 CJS Initial code

# This is the mallard data that ships with RMark

library(AICcmodavg)
library(emmeans)
library(ggplot2)
library(plyr)
library(readxl)

source(file.path("..","..","logistic-exposure-model.R"))

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
rm(malldata)

malldata2$Habitat <- factor(malldata2$Habitat)  # RMark wants categorical variable to be factors



# Fit a particular model
# This is a model with S varying by robel height
mod.rob <- glm(Survive~Robel,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.rob)

# You can make predictions just as before
pred.data <- data.frame(Robel=seq(min(malldata2$Robel), max(malldata2$Robel), length.out=50))
logit.dsr.pred.rob <- predict(mod.rob, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.rob <- cbind(pred.data, logit.dsr=logit.dsr.pred.rob$fit, logit.dsr.se=logit.dsr.pred.rob$se.fit)
dsr.rob$dsr <- expit(dsr.rob$logit.dsr)
dsr.rob$dsr.se <- dsr.rob$logit.dsr.se* dsr.rob$dsr * (1-dsr.rob$dsr)
dsr.rob


# we plot the results

dsr.rob.plot <- ggplot(data=dsr.rob, aes(x=Robel, y=dsr))+
  ggtitle("DSR by Robel height for mallards", subtitle="Logistic exposure model")+
  geom_line()+
  geom_ribbon(aes(ymin=expit(logit.dsr-1.95*logit.dsr.se), 
                  ymax=expit(logit.dsr+1.95*logit.dsr.se)),alpha=0.1)+
  ylab("DSR (95% ci)")
dsr.rob.plot 
ggsave(dsr.rob.plot,
       file=file.path("..","..","..","..","MyStuff","Images","mallard-dsr-robel-le.png"), h=4, w=6, units="in", dpi=300)





#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# fit a null models
mod.null <-  glm(Survive~1,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.null)

AICcmodavg::aictab( list(mod.rob, mod.null))


