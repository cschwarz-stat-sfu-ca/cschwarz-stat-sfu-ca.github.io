# Analysis of killdeer data illustrating basic RMark features
# Incorportating an age variable

# 2019-05-01 CJS Initial code

# This is the killdeer data that ships with RMark
# I added arbitrary nest ages.

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
#    AgeDay1: age of nest at day 1 of study (can be negative)
#
# In this example, the multiple visits to a nest have been collapsed
# to a single record for each nest.
# In more complex examples, you may have multple records per nest
# as shown in the mallard example.
#

killdata <- readxl::read_excel(file.path("..","Killdeer.xlsx"), 
                               sheet="killdeer-age")

killdata <- as.data.frame(killdata)  # sometimes doesn't play nice with tibbles
head(killdata)

# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
kill.proc <- process.data(killdata, model="Nest", nocc=max(killdata$LastChecked))
kill.proc


# 2. Examine and/or modify the ddl. (Not done here)
# The Age variable here is time since the start of the study which is not useful.
kill.ddl <- make.design.data(kill.proc)
str(kill.ddl)
kill.ddl 
  

# 3. Fit a particular model
# This is a model with S a function of nest age
mod.res <-  RMark::mark(kill.proc, ddl=kill.ddl,
                          model="Nest",
                          model.parameters=list(
                            S   =list(formula=~NestAge)
                          )
  )
summary(mod.res)


# Look the objects returned in more details
names(mod.res)
names(mod.res$results)

# look at estimates on beta and original scale
mod.res$results$beta  # on the logit scale

mod.res$results$real# on the regular 0-1 scale for each day AVERAGE OVER NEST AGE


# The real estimates are not useful because the value for each day is averaged over the
# nest ages at that time.
# For example, the beta values are shown above and the
# logit(DSR) for nest 1 day old is
logit_DSR_1 = sum( mod.res$results$beta$estimate *c(1,1))
logit_DSR_1
# and estimate of survival of nest 1 day old is
1/(1+exp(-logit_DSR_1))

# compared to 
head(mod.res$results$real)

# The average nest age at day 1 is
average_nest_age_1 = mean(killdata$AgeDay1)
average_nest_age_1
logit_DSR_1_avg = sum( mod.res$results$beta$estimate *c(1,average_nest_age_1))
logit_DSR_1_avg
# and estimate of survival of nest on day 1 at average age of nests is
1/(1+exp(-logit_DSR_1_avg))  # now matches the real estimates



# You want to predict survival as a function of nest age.
# This is a bit tricker

# First get the all.diff.index values for each day of the study.
get.real(mod.res, param="S", se=TRUE)
# we see that all.diff.index==1 is for nest survival on day 1 of the study

# we will predict then the DSR for day 1 at various ages
pred.ages <- data.frame(NestAge1=1:20, index=1)
covariate.predictions(mod.res, data=pred.ages )$estimates[1:10,]

# because the model does NOT include a term for day effects, these predictions will be the same
# for all days of the study
pred.ages <- data.frame(NestAge4=1:20, index=4)
covariate.predictions(mod.res, data=pred.ages )$estimates[1:10,]

plotdata <-covariate.predictions(mod.res, data=pred.ages )$estimates
dsr.age <- ggplot(data=plotdata, aes(x=NestAge4, y=estimate, group=1))+
  ggtitle("Survival as a function of nest age")+
  geom_line()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.1)
dsr.age
ggsave(dsr.age,
       file=file.path("..","..","..","..","MyStuff","Images","killdear-age.png"),h=4, w=6, units="in", dpi=300)
 

#_----------------------------------------------------------------------
#--------------------------------------------------------------------------
# Compare to the null model

mod.null <-  RMark::mark(kill.proc, ddl=kill.ddl,
                        model="Nest",
                        model.parameters=list(
                          S   =list(formula=~1)
                        )
)

collect.models(type="Nest")

cleanup(ask=FALSE)