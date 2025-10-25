# Analysis of mallard data illustrating logistic exposure models
# How to fit a nest level categorical variable (habitat)

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
# This is a model with S varying by habitat type
mod.hab <-  glm(Survive~Habitat,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.hab)


# You can make predictions just as before
pred.data <- data.frame(Habitat=unique(malldata2$Habitat))
logit.dsr.pred.hab <- predict(mod.hab, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.hab <- cbind(pred.data, logit.dsr=logit.dsr.pred.hab$fit, logit.dsr.se=logit.dsr.pred.hab$se.fit)
dsr.hab$dsr <- expit(dsr.hab$logit.dsr)
dsr.hab$dsr.se <- dsr.hab$logit.dsr.se* dsr.hab$dsr * (1-dsr.hab$dsr)
dsr.hab



# extract the logit(DSR) for each habitat using emmeans
mod.hab.emmo <- emmeans::emmeans(mod.hab, ~Habitat)
dsr.logit <- CLD(mod.hab.emmo)
dsr.logit

# Compute the multiple comparison on the ordinary scale
# We need to update the transformation so we can get answers on the original scale as well
mod.hab.rg <- update(ref_grid(mod.hab, at=list(exposure=1)), tran = logexp())

mod.hab.emmo2 <- emmeans::emmeans(mod.hab.rg, ~Habitat)
CLD(mod.hab.emmo2)
dsr <- CLD(mod.hab.emmo2, type="response") # on DSR scale
dsr

# pairwise effects on logit scale
summary(pairs(mod.hab.emmo2),infer=TRUE)

# pairwise effects on the DSR scale 
mod.hab.emmo3 <- regrid(emmeans(mod.hab, ~Habitat), transform="log")
summary(pairs(mod.hab.emmo3, type="response"), infer=TRUE)

dsr.hab.plot <- ggplot(data=dsr, aes(x=Habitat, y=prob))+
  ggtitle("DSR in different habitats for mallards", subtitle="Logistic exposure model")+
  geom_point()+
  geom_errorbar(aes(ymin=asymp.LCL, ymax=asymp.UCL), width=.1)+
  ylab("DSR (95% ci)")
dsr.hab.plot 
ggsave(dsr.hab.plot,
       file=file.path("..","..","..","..","MyStuff","Images","mallard-dsr-habitat-le.png"), h=4, w=6, units="in", dpi=300)






# you can get the nest survival probability in the usual way by using the estimated DSRs.
# we want to do this for each habitat 
ns <- plyr::ddply(dsr, "Habitat", function(x, ndays){
  #browser()
  ns    <- x$prob^ndays
  ns.se <- x$SE* ndays * x$prob^(ndays-1) 
  data.frame(ns=ns, ns.se=ns.se)
}, ndays=20)
ns


#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# fit a null models
mod.null <-  glm(Survive~1,
         family=binomial(link=logexp(malldata2$Exposure)),
         data=malldata2)
summary(mod.null)

AICcmodavg::aictab( list(mod.hab, mod.null))
