# Analysis of mallard data illustrating basic RMark features
# How to fit a nest level continous covariate (Robel height

# 2019-05-01 CJS Initial code

# This is the mallard data that ships with RMark

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
malldata <- as.data.frame(malldata)


# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable the number of sampling occasions in the data
mall.proc <- process.data(malldata, model="Nest", 
                          nocc=max(malldata$LastChecked))
mall.proc


# 2. Examine and/or modify the ddl. (Not done here)
mall.ddl <- make.design.data(mall.proc)
str(mall.ddl)
mall.ddl 
  

# 3. Fit a particular model
# This is a model with S varying by robel height
mod.rob <-  RMark::mark(mall.proc, ddl=mall.ddl,
                          model="Nest",
                          model.parameters=list(
                            S   =list(formula=~Robel)
                          )
  )
summary(mod.rob)

# estimates of DSR are for average value of Robel

# Look the objects returned in more details
names(mod.rob)
names(mod.rob$results)

# look at estimates on beta and original scale
mod.rob$results$beta  # on the logit scale

mod.rob$results$real# on the regular 0-1 scale for each habitat

# derived variabldes is the nest survival probability over the (nocc) days 
names(mod.rob$results$derived)

mod.rob$results$derived$"S Overall Survival"


# we need to use covariate predictions to get estimated DSR at different Robel heights
head(get.real(mod.rob, "S", se=TRUE))
# because the DSR depends on the Robel value and NOT the day, we can predict at day =1
# corresonding to index.all.diff=1

range(malldata$Robel)
pred.data <- data.frame(Robel=seq(min(malldata$Robel), max(malldata$Robel), length.out=50),
                        index=1)
head(pred.data)

# we plot the results
# because the DSR is the same for all days, we use the time=1 values
plotdata <- covariate.predictions(mod.rob, data=pred.data)$estimates
head(plotdata)

dsr.rob.plot <- ggplot(data=plotdata, aes(x=Robel, y=estimate))+
  ggtitle("DSR by Robel height for mallards")+
  geom_line()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl),alpha=0.1)+
  ylab("DSR (95% ci)")
dsr.rob.plot 
ggsave(dsr.rob.plot,
       file=file.path("..","..","..","..","MyStuff","Images","mallard-dsr-robel.png"), h=4, w=6, units="in", dpi=300)





#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# fit a null models
# 3. Fit a particular model
# This is a model with S varying by habitat type
mod.null <-  RMark::mark(mall.proc, ddl=mall.ddl,
                        model="Nest",
                        model.parameters=list(
                          S   =list(formula=~1)
                        )
)
summary(mod.null)

collect.models(type="Nest")

cleanup(ask=FALSE)
