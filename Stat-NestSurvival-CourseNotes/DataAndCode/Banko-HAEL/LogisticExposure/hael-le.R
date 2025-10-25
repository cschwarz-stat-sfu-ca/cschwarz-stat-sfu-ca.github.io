# Analysis of Hawaii Elepaio nesting data using the logistic exposure model
# 2019-05-01 CJS Initial code

# Paul C Banko, Kelly A Jaenecke, Robert W Peck, Kevin W Brinck, 
# Increased nesting success of Hawaii Elepaio in response to the removal 
#  of invasive black rats, 
# The Condor: Ornithological Applications, Volume 121, Issue 1, 1 February 2019, duz003,
# https://doi.org/10.1093/condor/duz003

# Data from
# Banko, P.C., Jaenecke, K.A., Peck, R.W., and
# Brinck, K.W., 2018, Hawaii Volcanoes National Park Elepaio
# nest monitoring and black rat mark recapture data 2015–2017:
# US Geological Survey data release, https://doi.org/10.5066/P93TOM58

# In Hawai‘i and other oceanic islands with few native land mammals, 
# black rats (Rattus rattus) are among the most damaging invasive 
# vertebrate species to native forest bird populations and habitats, 
# due to their arboreal behavior and generalist foraging habitats and habitat use. 
# We evaluated the nesting response of Hawai‘i ‘Elepaio (Chasiempis sandwichensis; Monarchidae), 
# a generalist insectivore, to the removal of black rats using rodenticide 
# in a before-after-control-impact study in high and low, mesic montane habitat 
# recovering from long-term damage from introduced ungulates and weeds. 
# We monitored nesting and rat activity during 2015–2016 before applying diphacinone 
# bait in 2017 to remove rats from two 700 x 700-m treatment plots that were paired 
# with two non-treatment plots. We continued monitoring through July 2017. 
# This data release consists of three tabular datasets including nest monitoring data, 
# rat capture data, and spatial data for all field plots.


library(AICcmodavg)
library(ggplot2)
library(plyr)
library(readxl)

source(file.path("..","..", "logistic-exposure-model.R"))

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

haeldata <- read.csv(file.path("..","HAEL_NestMonitoring_HAVO_2015-2017.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)
head(haeldata)


# We expand the data to generate the effective sample size
haeldata2 <- expand.nest.data(haeldata)
head(haeldata2)
rm(haeldata)


# create factor variable for year
xtabs(~year, data=haeldata2, exclude=NULL, na.action=na.pass)
haeldata2$yearF <- factor(haeldata2$year)


#--------------------------------------------
#  Set up the set of model to fit
model.list.csv <- textConnection(
" S
 ~NestAge + tmt
 ~          tmt
 ~          tmt + elev
")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

model.fits <- plyr::dlply(model.list, c("model.number","S"), function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- glm(formula=as.formula(paste("Survive", eval(x$S))),
            family=binomial(link=logexp(input.data$Exposure)),
            data=input.data)
  fit
  
},input.data=haeldata2)

# examine results for a particular model
model.number <-3
names(model.fits[[model.number]])
model.fits[[model.number]]$formula

summary(model.fits[[model.number]])


# Model comparision and averaging
# collect models and make AICc table
AICcmodavg::aictab(model.fits)

# results differ slight from paper because of approximation used for 
# NestAge in final intervals where the mid-point is used, but the individual
# ages are used in MARK during likelihood constructions. All of the results
# for the other models match perfectly.

