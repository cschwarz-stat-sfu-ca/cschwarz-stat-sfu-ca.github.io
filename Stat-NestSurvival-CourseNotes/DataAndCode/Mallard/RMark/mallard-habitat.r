# Analysis of mallard data illustrating basic RMark features
# How to fit a nest level categorical variable (habitat)

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

malldata$Habitat <- factor(malldata$Habitat)  # RMark wants categorical variable to be factors

# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable the number of sampling occasions in the data
mall.proc <- process.data(malldata, model="Nest", 
                          group="Habitat", 
                          nocc=max(malldata$LastChecked))
mall.proc


# 2. Examine and/or modify the ddl. (Not done here)
mall.ddl <- make.design.data(mall.proc)
str(mall.ddl)
mall.ddl 
  

# 3. Fit a particular model
# This is a model with S varying by habitat type
mod.hab <-  RMark::mark(mall.proc, ddl=mall.ddl,
                          model="Nest",
                          model.parameters=list(
                            S   =list(formula=~Habitat)
                          )
  )
summary(mod.hab)


# Look the objects returned in more details
names(mod.hab)
names(mod.hab$results)

# look at estimates on beta and original scale
mod.hab$results$beta  # on the logit scale

mod.hab$results$real# on the regular 0-1 scale for each habitat

# derived variabldes is the nest survival probability over the (nocc) days 
names(mod.hab$results$derived)

mod.hab$results$derived$"S Overall Survival"



# alternatively
get.real(mod.hab, "S", se=TRUE)

# we plot the results
# because the DSR is the same for all days, we use the time=1 values
plotdata <- get.real(mod.hab, "S", se=TRUE)
head(plotdata)

plotdata <- plotdata[plotdata$time==1,]
plotdata

dsr.hab.plot <- ggplot(data=plotdata, aes(x=Habitat, y=estimate))+
  ggtitle("DSR in different habitats for mallards")+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  ylab("DSR (95% ci)")
dsr.hab.plot 
ggsave(dsr.hab.plot,
       file=file.path("..","..","..","..","MyStuff","Images","mallard-dsr-habitat.png"), h=4, w=6, units="in", dpi=300)


# you can get the nest survival probability in the usual way by using the estimated DSRs.
# we want to do this for each habitat 
ns <- plyr::ldply(unique(malldata$Habitat), function(habitat, ndays){
  daily.dsr <- get.real(mod.hab, "S", se=TRUE, vcv=TRUE, expand=TRUE)
  # this has all habitats together; select only habitat of interest
  select <- daily.dsr$estimates$Habitat==habitat
  daily.dsr$estimates <- daily.dsr$estimates[select,]
  daily.dsr$vcv.real  <- daily.dsr$vcv.real[select,select]
  #browser()
  ns    <- prod(daily.dsr$estimates$estimate[1:ndays])
  ns.se <- deltamethod.special("prod", 
                      daily.dsr$estimates$estimate[1:ndays],
                      daily.dsr$vcv.real[1:ndays, 1:ndays])
  data.frame(Habitat=habitat, ns=ns, ns.se=ns.se)
}, ndays=20)
ns


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