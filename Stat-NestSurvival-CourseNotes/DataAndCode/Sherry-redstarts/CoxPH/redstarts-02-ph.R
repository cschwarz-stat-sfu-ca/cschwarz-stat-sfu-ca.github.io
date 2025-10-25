# Analysis of Shelly redstart data using the proportional hazard models
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
library(survival)

source(file.path("..","..","CoxPH-additional-functions.R"))

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
reddata2 <- expand.nest.data.ph(reddata)
head(reddata2)
rm(reddata)


# as noted in the paper, we remove nests from years with no baffeling
dim(reddata2)
reddata2 <- reddata2[ !(reddata2$BaffleStatus=="N" & reddata2$Year %in% c(1983, 1984, 1990)),]
dim(reddata2)

# create factor variable for year
reddata2$YearF <- factor(reddata2$Year)

# create factor for baffle
reddata2$BaffleStatus <- factor(reddata2$BaffleStatus)

#--------------------------------------------
#  Set up the set of model to fit
model.list.csv <- textConnection(
" S
 ~DBH
 ~DBH+YearF
 ~DBH+YearF+BaffleStatus
")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

model.fits <- plyr::dlply(model.list, c("model.number","S"), function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- coxph(formula=as.formula(paste("Surv", eval(x$S))),data=input.data)
  fit
  
},input.data=reddata2)

# examine results for a particular model
model.number <-3
names(model.fits[[model.number]])
model.fits[[model.number]]$formula

summary(model.fits[[model.number]])


# Model comparision and averaging
# collect models and make AICc table
AICcmodavg::aictab(model.fits)


# look at impact of baffle status (averaged over all years)
fit.emmo <- emmeans::emmeans(model.fits[[3]], ~BaffleStatus)
pairs(fit.emmo)
summary(pairs(fit.emmo, type="response"), infer=TRUE)
