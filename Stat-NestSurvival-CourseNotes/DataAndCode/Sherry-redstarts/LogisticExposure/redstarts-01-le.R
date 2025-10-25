# Analysis of Redstart data using logisitic exposure models

# 2019-05-01 CJS Initial code

# Sherry TW, Wilson S, Hunter S, Holmes RT (2015) 
# Impacts of nest predators and weather on reproductive success and 
# population limitation in a long-distance migratory songbird. 
# Journal of Avian Biology 46(6): 559-569. https://doi.org/10.1111/jav.00536

# Data from
# Sherry TW, Wilson S, Hunter S, Holmes RT (2015) 
# Data from: Impacts of nest predators and weather on reproductive success and 
# population limitation in a long-distance migratory songbird. 
# Dryad Digital Repository. https://doi.org/10.5061/dryad.73870


library(AICcmodavg)
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
# In this example, the multiple visits to a nest have been collapsed
# to a single record for each nest.
# In more complex examples, you may have multple records per nest
# as shown in the mallard example.
#

reddata <- readxl::read_excel(file.path("..","Sherry.xlsx"), 
                               sheet="NestData")
head(reddata)

# We expand the data to generate the effective sample size
reddata2 <- expand.nest.data(reddata)
head(reddata2)
rm(reddata)



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Baffled nests only

dim(reddata2)
reddata2.baf <- reddata2[ reddata2$BaffleStatus=="Y",]
dim(reddata2.baf)


# Fit a particular model
# This is a model with S constant over time (closest to Mayfield method)
mod.baf <-  glm(Survive~1,
         family=binomial(link=logexp(reddata2.baf$Exposure)),
         data=reddata2.baf)
summary(mod.baf)


# Convert the logit(DSR) to DSR
DSR <- expit(coef(mod.baf))
DSR.se <- arm::se.coef(mod.baf)*DSR*(1-DSR)
cat("DSR ", DSR, "(SE ", DSR.se, ")\n")

# Find confidence intervals by taking expit of confit of coefficients
expit(confint(mod.baf))

# make a data frame for model averaging
dsr.dot <- data.frame(Day=1:39, logit.dsr=coef(mod.baf), logit.dsr.sr=arm::se.coef(mod.baf),
                                dsr=DSR, dsr.se=DSR.se)

# Compute the nest survival for 20 days
days <- 20
NS <- DSR^days
NS.se <- DSR.se * days * DSR^(days-1)
cat("NS ", days," days ", NS, "(SE ", NS.se, ")\n")



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# control nests only

dim(reddata2)
reddata2.cntl <- reddata2[ reddata2$BaffleStatus=="N" & !reddata2$Year %in% c(1983, 1984, 1990),]
dim(reddata2.cntl)



# 3. Fit a particular model
# This is a model with S constant over time (closest to Mayfield method)
mod.cntl <-  glm(Survive~1,
         family=binomial(link=logexp(reddata2.cntl$Exposure)),
         data=reddata2.cntl)
summary(mod.cntl)


# Convert the logit(DSR) to DSR
DSR <- expit(coef(mod.cntl))
DSR.se <- arm::se.coef(mod.cntl)*DSR*(1-DSR)
cat("DSR ", DSR, "(SE ", DSR.se, ")\n")

# Find confidence intervals by taking expit of confit of coefficients
expit(confint(mod.cntl))

# Compute the nest survival for 20 days
days <- 20
NS <- DSR^days
NS.se <- DSR.se * days * DSR^(days-1)
cat("NS ", days," days ", NS, "(SE ", NS.se, ")\n")


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# linear effect of date (relative to 27 May) 

dim(reddata2)

# 3. Fit a particular model
# This is a model with linear effect of date (relative to 27 May)
mod.T <-  glm(Survive~Day,
            family=binomial(link=logexp(reddata2$Exposure)),
            data=reddata2)
summary(mod.T)


# predict the survival probability as a function of time since 21 May
pred.data <- data.frame(Day=1:60)
logit.dsr.pred.linear <- predict(mod.T, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.linear <- cbind(pred.data, logit.dsr=logit.dsr.pred.linear$fit, logit.dsr.se=logit.dsr.pred.linear$se.fit)
head(dsr.linear)
dsr.linear$dsr <- expit(dsr.linear$logit.dsr)
dsr.linear$dsr.se <- dsr.linear$logit.dsr.se* dsr.linear$dsr * (1-dsr.linear$dsr)
head(dsr.linear)


plotdata <- dsr.linear
# add approximate ci on dsr
plotdata$dsr.lcl <- expit(plotdata$logit.dsr - 1.96*plotdata$logit.dsr.se)
plotdata$dsr.ucl <- expit(plotdata$logit.dsr + 1.96*plotdata$logit.dsr.se)

ggplot(data=plotdata, aes(x=Day, y=dsr, group=1))+
  ggtitle("Linear effect on DSR")+
  geom_line()+
  geom_ribbon(aes(ymin=dsr.lcl, ymax=dsr.ucl), alpha=0.2)+
  ylim(0.5,1)

# To compute the nest survival, you also need the covariance between all of the estimates
# This is NOT produced by the predict function and so we need to do this by hand (groan)!

# for what values do you want predictions?
days=1:20
pred.matrix <- data.frame(intecept=1, data=days)
head(pred.matrix)

# compute the predictions and variance covariance matrix
pred.logit <- as.vector( as.matrix(pred.matrix) %*% coef(mod.T))
pred.logit.vcov <- as.matrix(pred.matrix) %*% vcov(mod.T) %*% t(as.matrix(pred.matrix))

# use msm package to get estimates and vcov after expit()
# create a list of functions to back transform
library(msm)
funs <- paste("~1/(1+exp(-x",days,"))",sep="")
funs <- lapply(funs, as.formula)
head(funs)

# find the inverse predictions and vcov
pred.dsr     <- expit(pred.logit)
pred.dsr
pred.dsr.vcov <- msm::deltamethod(funs, mean=pred.logit, cov=pred.logit.vcov, ses=FALSE)

# now find the product of the dsr over the 20 time periods
funs <- as.formula( paste("~", paste(paste("x",days,sep=""),collapse=" * "), collapse=""))
funs
ns <- prod(pred.dsr)
ns.se <- msm::deltamethod(funs, mean=pred.dsr, cov=pred.dsr.vcov, ses=TRUE)
cat("NS ", min(days)," to ",max(days),"days:", NS, "(SE ", NS.se, ")\n")

