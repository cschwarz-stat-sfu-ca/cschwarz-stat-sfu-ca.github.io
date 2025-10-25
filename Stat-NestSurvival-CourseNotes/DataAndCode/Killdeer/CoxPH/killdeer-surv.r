# Analysis of killdeer data using survival models
# Cannot use Cox proportional hazard models because no groups here
# 2019-05-01 CJS Initial code

# This is the killdeer data that ships with RMark

library(AICcmodavg)
library(ggplot2)
library(plyr)
library(readxl)
library(survival)

source(file.path("..","CoxPH-additional-functions.R"))

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

killdata <- readxl::read_excel("Killdeer.xlsx", 
                               sheet="killdeer")
head(killdata)

# We expand the data to generate the effective sample size
killdata2 <- expand.nest.data(killdata)
head(killdata2)


# Now to fit the logistic exposure model
fit.Sdot <- glm(Survive~1,
         family=binomial(link=logexp(killdata2$Exposure)),
         data=killdata2)
summary(fit.Sdot)
-2*logLik(fit.Sdot)

# Convert the logit(DSR) to DSR
DSR <- expit(coef(fit.Sdot))
DSR.se <- arm::se.coef(fit.Sdot)*DSR*(1-DSR)
cat("DSR ", DSR, "(SE ", DSR.se, ")\n")

# Find confidence intervals by taking expit of confit of coefficients
expit(confint(fit.Sdot))

# make a data frame for model averaging
dsr.dot <- data.frame(Day=1:39, logit.dsr=coef(fit.Sdot), logit.dsr.sr=arm::se.coef(fit.Sdot),
                                dsr=DSR, dsr.se=DSR.se)

# Compute the nest survival
days <- 39
NS <- DSR^days
NS.se <- DSR.se * days * DSR^(days-1)
cat("NS ", days," days ", NS, "(SE ", NS.se, ")\n")

# fit a survival model with a constant hazard which is equivalent to exponential survival model
# We need to set up the intervals
killdata2$start <- killdata2$Day
select <-  killdata2$LastChecked > killdata2$LastPresent & killdata2$Survive==0
killdata2$start[ select] <- killdata2$LastPresent[ select]
killdata2$stop   <- killdata2$Day+1
killdata2$stop[ select] <- killdata2$LastChecked[select]

head(killdata2, n=20)
tail(killdata2)

killdata2$Surv <- Surv(time=killdata2$start, time2=killdata2$stop, event=1-killdata2$Survive, type="counting")
head(killdata2, n=20)

fit <- survreg(Surv ~ 1, dist='exp', data=killdata2)





#--------------------------------------------------------------
# Linear in date
# This is a close approximation as the final interval uses the midpoint of the interval 
# rather than the individual daily survival probabilities computed with 
# the actual date

head(killdata2)

fit.linear <- glm(Survive~Day,
         family=binomial(link=logexp(killdata2$Exposure)),
         data=killdata2)
summary(fit.linear)
-2*logLik(fit.linear)

# We now want to predict the DSR for days 1..39
# In this case predict() will work on the logit scale because the exposure=1 default is ok, 
# but we need to expit the results
pred.data <- data.frame(Day=1:39)
logit.dsr.pred.linear <- predict(fit.linear, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.linear <- cbind(pred.data, logit.dsr=logit.dsr.pred.linear$fit, logit.dsr.se=logit.dsr.pred.linear$se.fit)
head(dsr.linear)
dsr.linear$dsr <- expit(dsr.linear$logit.dsr)
dsr.linear$dsr.se <- dsr.linear$logit.dsr.se* dsr.linear$dsr * (1-dsr.linear$dsr)
head(dsr.linear)

# plot this in the usual way
killdeer.linear.le <- ggplot(data=dsr.linear, aes(x=Day, y=dsr, group=1))+
  ggtitle("Estimated linear fit on logit scale", subtitle="Logistic exposure model")+
  geom_line()+
  geom_ribbon(aes(ymin=expit(logit.dsr-1.96*logit.dsr.se),
                  ymax=expit(logit.dsr+1.96*logit.dsr.se)), alpha=.1)+
  ylim(0.5, 1)
killdeer.linear.le
ggsave(killdeer.linear.le,
       file=file.path("..","..","..","MyStuff","Images","killdeer-linear-S-le.png"),h=4, w=6, units="in", dpi=300)

#--------------------------------------------------------------
# Quadratic in date
# This is a close approximation as the final interval uses the midpoint of the interval 
# rather than the individual daily survival probabilities computed with 
# the actual date
killdata2$Day2 <- (killdata2$Day-20)^2

fit.quad <- glm(Survive~Day+Day2,
         family=binomial(link=logexp(killdata2$Exposure)),
         data=killdata2)
summary(fit.quad)
-2*logLik(fit.quad)

# We now want to predict the DSR for days 1..39
# In this case predict() will work on the logit scale because the exposure=1 default is ok, 
# but we need to expit the results
pred.data <- data.frame(Day=1:39)
pred.data$Day2 <- (pred.data$Day-20)^2  # we need to match coding above

logit.dsr.pred.quad <- predict(fit.quad, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.quad  <- cbind(pred.data, logit.dsr=logit.dsr.pred.quad$fit, logit.dsr.se=logit.dsr.pred.quad$se.fit)
head(dsr.quad )
dsr.quad $dsr <- expit(dsr.quad $logit.dsr)
dsr.quad $dsr.se <- dsr.quad $logit.dsr.se* dsr.quad$dsr * (1-dsr.quad$dsr)
head(dsr.quad)

# plot this in the usual way
killdeer.quad.le <- ggplot(data=dsr.quad, aes(x=Day, y=dsr, group=1))+
  ggtitle("Estimated quadratic fit on logit scale", subtitle="Logistic exposure model")+
  geom_line()+
  geom_ribbon(aes(ymin=expit(logit.dsr-1.96*logit.dsr.se),
                  ymax=expit(logit.dsr+1.96*logit.dsr.se)), alpha=.1)+
  ylim(0.5, 1)
killdeer.quad.le
ggsave(killdeer.quad.le,
       file=file.path("..","..","..","MyStuff","Images","killdeer-quad-S-le.png"),h=4, w=6, units="in", dpi=300)



#----------------------------------------------------------------------------
# DSR varies in the first and second half of the study
# need to define the part of the study
# some ambiguity if an interval straddles the breakpoint

killdata2$studyhalf <- car::recode(killdata2$Day,
                                  " lo:20='1st'; 21:hi='2nd'")
head(killdata2)

fit.2part <- glm(Survive~studyhalf,
         family=binomial(link=logexp(killdata2$Exposure)),
         data=killdata2)
summary(fit.2part)
-2*logLik(fit.2part)

# We now want to predict the DSR for days 1..39
# In this case predict() will work on the logit scale because the exposure=1 default is ok, 
# but we need to expit the results
pred.data <- data.frame(Day=1:39)
pred.data$studyhalf <- car::recode(pred.data$Day,
                                  " lo:20='1st'; 21:hi='2nd'")
logit.dsr.pred.2part <- predict(fit.2part, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.2part <- cbind(pred.data, logit.dsr=logit.dsr.pred.2part$fit, logit.dsr.se=logit.dsr.pred.2part$se.fit)
head(dsr.2part)
dsr.2part$dsr <- expit(dsr.2part$logit.dsr)
dsr.2part$dsr.se <- dsr.2part$logit.dsr.se* dsr.2part$dsr * (1-dsr.2part$dsr)
head(dsr.2part)
tail(dsr.2part)

# plot this in the usual way
killdeer.2part.le <- ggplot(data=dsr.2part, aes(x=Day, y=dsr, group=1))+
  ggtitle("Estimated 2part fit on logit scale", subtitle="Logistic exposure model")+
  geom_line()+
  geom_ribbon(aes(ymin=expit(logit.dsr-1.96*logit.dsr.se),
                  ymax=expit(logit.dsr+1.96*logit.dsr.se)), alpha=.1)+
  ylim(0.5, 1)
killdeer.2part.le
ggsave(killdeer.2part.le,
       file=file.path("..","..","..","MyStuff","Images","killdeer-2part-S-le.png"),h=4, w=6, units="in", dpi=300)


