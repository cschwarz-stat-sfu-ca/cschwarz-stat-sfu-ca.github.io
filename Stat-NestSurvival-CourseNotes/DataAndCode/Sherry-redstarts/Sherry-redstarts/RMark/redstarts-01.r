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


# Look the objects returned in more details
names(mod.baf)
names(mod.baf$results)

# look at estimates on beta and original scale
mod.baf$results$beta  # on the logit scale

mod.baf$results$real# on the regular 0-1 scale for each site

# derived variabldes is the nest survival probability over the 40 (nocc) days 
names(mod.baf$results$derived)

mod.baf$results$derived$"S Overall Survival"


# alternatively
get.real(mod.baf, "S", se=TRUE)


# Notice that the nest survival is the product of the individual S
prod(get.real(mod.baf, "S", se=TRUE)$estimate)


# problem is that the derived parameter is computed over 61 days (nocc)
# we only want the probability computed over 20 days.
# Use the deltamethod.special function

# for all 61 days
prod(get.real(mod.baf, "S", se=TRUE)$estimate)
deltamethod.special("prod", 
                    get.real(mod.baf, "S", se=TRUE)$estimate,
                    get.real(mod.baf, "S", se=TRUE, vcv=TRUE, expand=TRUE)$vcv.real)

# for 20 days
prod(get.real(mod.baf, "S", se=TRUE)$estimate[1:20])
deltamethod.special("prod", 
                    get.real(mod.baf, "S", se=TRUE)$estimate[1:20],
                    get.real(mod.baf, "S", se=TRUE, vcv=TRUE, expand=TRUE)$vcv.real[1:20, 1:20])



#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# control nests only

dim(reddata)
reddata.cntl <- reddata[ reddata$BaffleStatus=="N" & !reddata$Year %in% c(1983, 1984, 1990),]
dim(reddata.cntl)



# 1. Process the data.
# The nocc variable is the data at which hatching occurs
red.cntl.proc <- process.data(reddata.cntl, model="Nest", nocc=max(reddata.cntl$LastChecked))
red.cntl.proc


# 2. Examine and/or modify the ddl. (Not done here)
red.cntl.ddl <- make.design.data(red.cntl.proc)
str(red.cntl.ddl)
red.cntl.ddl 


# 3. Fit a particular model
# This is a model with S constant over time (closest to Mayfield method)
mod.cntl <-  RMark::mark(red.cntl.proc, ddl=red.cntl.ddl,
                        model="Nest",
                        model.parameters=list(
                          S   =list(formula=~1)
                        )
)
summary(mod.cntl)


# Look the objects returned in more details
names(mod.cntl)
names(mod.cntl$results)

# look at estimates on beta and original scale
mod.cntl$results$beta  # on the logit scale

mod.cntl$results$real# on the regular 0-1 scale for each site

# derived variabldes is the nest survival probability over the 40 (nocc) days 
names(mod.cntl$results$derived)

mod.cntl$results$derived$"S Overall Survival"


# alternatively
get.real(mod.cntl, "S", se=TRUE)


# Notice that the nest survival is the product of the individual S
prod(get.real(mod.cntl, "S", se=TRUE)$estimate)


# problem is that the derived parameter is computed over 61 days (nocc)
# we only want the probability computed over 20 days.
# Use the deltamethod.special function

# for all 61 days
prod(get.real(mod.cntl, "S", se=TRUE)$estimate)
deltamethod.special("prod", 
                    get.real(mod.cntl, "S", se=TRUE)$estimate,
                    get.real(mod.cntl, "S", se=TRUE, vcv=TRUE, expand=TRUE)$vcv.real)

# for 20 days
prod(get.real(mod.cntl, "S", se=TRUE)$estimate[1:20])
deltamethod.special("prod", 
                    get.real(mod.cntl, "S", se=TRUE)$estimate[1:20],
                    get.real(mod.cntl, "S", se=TRUE, vcv=TRUE, expand=TRUE)$vcv.real[1:20, 1:20])


#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
# linear effect of date (relative to 27 May) 

dim(reddata)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
red.all.proc <- process.data(reddata, model="Nest", nocc=max(reddata$LastChecked))
red.all.proc


# 2. Examine and/or modify the ddl. (Not done here)
red.all.ddl <- make.design.data(red.all.proc)
str(red.all.ddl)
red.all.ddl 


# 3. Fit a particular model
# This is a model with linear effect of date (relative to 27 May)
mod.T <-  RMark::mark(red.all.proc, ddl=red.all.ddl,
                         model="Nest",
                         model.parameters=list(
                           S   =list(formula=~Time)
                         )
)
summary(mod.T)


# Look the objects returned in more details
names(mod.T)
names(mod.T$results)

# look at estimates on beta and original scale
mod.T$results$beta  # on the logit scale

mod.T$results$real# on the regular 0-1 scale for each site

# derived variabldes is the nest survival probability over the 40 (nocc) days 
names(mod.T$results$derived)

mod.T$results$derived$"S Overall Survival"


# alternatively
get.real(mod.T, "S", se=TRUE)

plotdata <- get.real(mod.T, "S", se=TRUE)
ggplot(data=plotdata, aes(x=time, y=estimate, group=1))+
  ggtitle("Linear effect on DSR")+
  geom_line()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)+
  ylim(0.5,1)


# Notice that the nest survival is the product of the individual S
prod(get.real(mod.T, "S", se=TRUE)$estimate)


# problem is that the derived parameter is computed over 61 days (nocc)
# we only want the probability computed over 20 days but starting when?.
# Use the deltamethod.special function

# for 20 days at start of experiment/
prod(get.real(mod.T, "S", se=TRUE)$estimate[1:20])
deltamethod.special("prod", 
                    get.real(mod.T, "S", se=TRUE)$estimate[1:20],
                    get.real(mod.T, "S", se=TRUE, vcv=TRUE, expand=TRUE)$vcv.real[1:20, 1:20])


# for 20 days at middle of experiment/
prod(get.real(mod.T, "S", se=TRUE)$estimate[21:40])
deltamethod.special("prod", 
                    get.real(mod.T, "S", se=TRUE)$estimate[21:40],
                    get.real(mod.T, "S", se=TRUE, vcv=TRUE, expand=TRUE)$vcv.real[21:40, 21:40])



cleanup(ask=FALSE)
