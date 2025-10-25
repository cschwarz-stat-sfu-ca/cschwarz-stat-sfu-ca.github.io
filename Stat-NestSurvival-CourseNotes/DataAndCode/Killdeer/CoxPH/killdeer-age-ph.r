# Analysis of killdeer data using the proportional hazerds model
# 2019-05-01 CJS Initial code

# This is the killdeer data that ships with RMark

library(AICcmodavg)
library(ggplot2)
library(plyr)
library(readxl)
library(survival)
library(survminer)

source(file.path("..","..","CoxPH-additional-functions.R"))  # need to expand the nest data
# The dataframe must contain the following fields with the following names
#
#    FirstFound: day the nest was first found
#    LastPresent: last day that a chick was present in the nest
#    LastChecked: last day the nest was checked
#    Fate: fate of the nest; 0=hatch an
#    Freq: number of nests with this data
#    AgeDay1 - age at day 1 so that age of best at each day can be imputed
#
# In this example, the multiple visits to a nest have been collapsed
# to a single record for each nest.
# In more complex examples, you may have multple records per nest
# as shown in the mallard example.
#

killdata <- readxl::read_excel(file.path("..","Killdeer.xlsx"), 
                               sheet="killdeer-age")
head(killdata)

# There are no covariates except nest age
killdata2 <- expand.nest.data.ph(killdata)
head(killdata2)
killdata2

#-------------------------------
# Caution -- Age and Time are very highly correlated!
ggplot(data=killdata2, aes(x=Start, y=NestAge))+
  geom_point()



fit.Sdot.ph <- coxph(Surv ~1,   data=killdata2)
summary(fit.Sdot.ph)

#-----------------
fit.SNestAge.ph <- coxph(Surv~NestAge, data=killdata2)
summary(fit.SNestAge.ph)

cox.zph(fit.SNestAge.ph)
ggcoxzph(cox.zph(fit.SNestAge.ph))

png(file=file.path("..","..","..","..","MyStuff","Images","killdeer-ph-age-testprop.png"), h=4, w=6, units="in", res=300)
ggcoxzph(cox.zph(fit.SNestAge.ph))
dev.off()


# baseline cumulative hazard
cumhaz <- basehaz(fit.SNestAge.ph)
# estimate change in cumulative hazard and plot
cumhaz$deltaHaz <- c(NA,diff(cumhaz$hazard))
ggplot(data=cumhaz, aes(x=time, y=deltaHaz))+
  ggtitle("Estimated baseline hazard function over time")+
  geom_point()+
  geom_smooth(se=FALSE)


# estimated hazard curves at three ages 
pred.data <- expand.grid( NestAge=seq(min(killdata2$NestAge),max(killdata2$NestAge), length.out=3))
pred.survival <- survfit(fit.SNestAge.ph, newdata=pred.data, se.fit=TRUE)
plot(pred.survival)




AICcmodavg::aictab(list(fit.Sdot.ph, fit.SNestAge.ph))


# Now to fit the logistic exposure model
killdata3 <- expand.nest.data(killdata)
fit.Sdot <- glm(Survive~1,
         family=binomial(link=logexp(killdata3$Exposure)),
         data=killdata3)
fit.SNestAge <- glm(Survive~NestAge,
         family=binomial(link=logexp(killdata3$Exposure)),
         data=killdata3)
AICcmodavg::aictab(list(fit.Sdot, fit.SNestAge))



