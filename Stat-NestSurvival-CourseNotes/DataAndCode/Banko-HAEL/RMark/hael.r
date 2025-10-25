# Analysis of Hawaii Elepaio nesting data using the nest survival model of RMark
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

library(ggplot2)
library(readxl)
library(RMark)


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


# create factor variable for year and treatment and elevation
xtabs(~year, data=haeldata, exclude=NULL, na.action=na.pass)
haeldata$yearF <- factor(haeldata$year)
haeldata$tmt   <- factor(haeldata$tmt )
haeldata$elev  <- factor(haeldata$elev)

# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
hael.proc <- process.data(haeldata, model="Nest", group=c("tmt","elev"), nocc=max(haeldata$LastChecked))
hael.proc


# 2. Examine and/or modify the ddl. (Not done here)
hael.ddl <- make.design.data(hael.proc)
str(hael.ddl)
hael.ddl 
  

# 3. Set up the set of model to fit
model.list.csv <- textConnection(
" S
  ~NestAge + tmt
  ~          tmt
  ~          tmt + elev
  ")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="Nest",
                     model.parameters=list(
                       S   =list(formula=as.formula(eval(x$S)))
                     )
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=hael.proc, input.ddl=hael.ddl)


# examine individula model results
model.number <-3
names(model.fits[[model.number]])
model.fits[[model.number]]$model.name

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "S", se=TRUE)


# Model comparision and averaging
# collect models and make AICc table

model.set <- RMark::collect.models( type="Nest")
model.set

names(model.set)
model.set$model.table


# model averaged values
get.real(model.set[[1]], "S", se=TRUE)
get.real(model.set[[2]], "S", se=TRUE)
get.real(model.set[[3]], "S", se=TRUE)

S.ma <- RMark::model.average(model.set, param="S")
head(S.ma)

cleanup(ask=FALSE)
