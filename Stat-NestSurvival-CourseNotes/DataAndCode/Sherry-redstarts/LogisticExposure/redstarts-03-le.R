# Analysis of Redstart data illustrating logistic exposure models
# Reproduce Table 1b.

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
reddata <- as.data.frame(reddata)


# get the yearly date and merge with the nest data
annual.data <- readxl::read_excel(file.path("..","Sherry.xlsx"),
                                  sheet="AnnualCovariates", skip=9)
head(annual.data)

dim(reddata)
reddata <- merge(reddata, annual.data, all.x=TRUE)
dim(reddata)

# any missing data?
sum(!complete.cases(reddata))

# as noted in the paper, we remove nests with baffeling
dim(reddata)
reddata <- reddata[ !(reddata$BaffleStatus=="Y"),]
dim(reddata)

# expand the data
reddata2 <- expand.nest.data(reddata)
rm(reddata)

# create factor variable for year
reddata2$YearF <- factor(reddata2$Year)

#-------------------------------------
# Set up the set of model to fit
model.list.csv <- textConnection(
  " S
  ~DBH
  ~DBH+Predators+ MayTemp                   +JuneRain  +NestAge
  ~DBH+Predators+ MayTemp          +MayRain +JuneRain  +NestAge
  ~DBH+Predators+ MayTemp+JuneTemp          +JuneRain  +NestAge
  ~DBH+Predators+ MayTemp                   +JuneRain  +NestAge + Density
  ~DBH+Predators+ MayTemp                   +JuneRain  
  ~DBH+Predators+ MayTemp          +MayRain            +NestAge 
  ~DBH+Predators+ MayTemp                              +NestAge 
  ~DBH+Predators+ MayTemp          +MayRain            +NestAge + Density
  ")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list


model.fits <- plyr::dlply(model.list, c("S","model.number"), function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- glm(formula=as.formula(paste("Survive", eval(x$S))),
            family=binomial(link=logexp(input.data$Exposure)),
            data=input.data)
  fit
  
},input.data=reddata2)


# examine individual model results
model.number <-1
names(model.fits[[model.number]])
model.fits[[model.number]]$formula

summary(model.fits[[model.number]])


# Model comparision and averaging
# collect models and make AICc table
AICcmodavg::aictab(model.fits)


# Estimates of beta from the top model
summary(model.fits[[2]])

