# Analysis of Redstart data illustrating basic RMark features
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

# create factor variable for year
reddata$YearF <- factor(reddata$Year)

# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
red.proc <- process.data(reddata, model="Nest", group=c("YearF"), nocc=max(reddata$LastChecked))
red.proc


# 2. Examine and/or modify the ddl. (Not done here)
red.ddl <- make.design.data(red.proc)
str(red.ddl)
red.ddl 
  

# 3. Set up the set of model to fit
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
  
},input.data=red.proc, input.ddl=red.ddl)


# examine individula model results
model.number <-1
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


# Estimates of beta from the top model
model.fits[[2]]$results$beta

# there is no automatic function in \RMark to model average beta
# parameters but view
# http://www.phidot.org/forum/viewtopic.php?t=996&postdays=0&postorder=asc&start=0
# I would put the beta parameters into the derived parameters area and the use the 
# additional RMark functions to model average derived parameters

cleanup(ask=FALSE)
