# Analysis of killdeer data illustrating basic RMark features

# 2019-05-01 CJS Initial code

# This is the killdeer data that ships with RMark
# which has been saved as a *.csv file

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

killdata <- readxl::read_excel(file.path("..","Killdeer.xlsx"), 
                               sheet="killdeer")
head(killdata)


# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
kill.proc <- process.data(killdata, model="Nest", nocc=40)
kill.proc


# 2. Examine and/or modify the ddl. (Not done here)
kill.ddl <- make.design.data(kill.proc)
str(kill.ddl)
kill.ddl 
  

# 3. Fit a particular model
# This is a model with S constant over time (closest to Mayfield method)
mod.res <-  RMark::mark(kill.proc, ddl=kill.ddl,
                          model="Nest",
                          model.parameters=list(
                            S   =list(formula=~1)
                          )
  )
summary(mod.res)


# Look the objects returned in more details
names(mod.res)
names(mod.res$results)

# look at estimates on beta and original scale
mod.res$results$beta  # on the logit scale

mod.res$results$real# on the regular 0-1 scale for each site

# derived variabldes is the nest survival probability over the 40 (nocc) days 
names(mod.res$results$derived)

mod.res$results$derived$"S Overall Survival"



# alternatively
get.real(mod.res, "S", se=TRUE)


# Notice that the nest survival is the product of the individual S
prod(get.real(mod.res, "S", se=TRUE)$estimate)


