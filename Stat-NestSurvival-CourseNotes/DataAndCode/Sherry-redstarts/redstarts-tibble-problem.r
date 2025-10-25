# Analysis of Redstart data illustrating basic RMark features

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
reddata <- as.data.frame(reddata) # Rmark gets upset with tibbles
head(reddata)

range(reddata$FirstFound)
range(reddata$LastPresent)
range(reddata$LastChecked)


# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# Baffled nests only

dim(reddata)
reddata.baf <- reddata[ reddata$BaffleStatus=="Y",]
dim(reddata.baf)



# 1. Process the data.
# The nocc variable is the data at which hatching occurs
red.baf.proc <- process.data(reddata.baf, model="Nest", nocc=max(reddata.baf$LastChecked))
red.baf.proc


# 2. Examine and/or modify the ddl. (Not done here)
red.baf.ddl <- make.design.data(red.baf.proc)
str(red.baf.ddl)
red.baf.ddl 
  

# 3. Fit a particular model
# This is a model with S constant over time (closest to Mayfield method)
mod.baf <-  RMark::mark(red.baf.proc, ddl=red.baf.ddl,
                          model="Nest",
                          model.parameters=list(
                            S   =list(formula=~1)
                          )
  )
summary(mod.baf)


