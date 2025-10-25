# Analysis of killdeer data illustrating logistic exposure models
# Incorportating an age variable

# 2019-05-01 CJS Initial code

# This is the killdeer data that ships with RMark
# I added arbitrary nest ages.

library(AICcmodavg)
library(ggplot2)
library(plyr)
library(readxl)

source(file.path("..","..","logistic-exposure-model.R"))


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

killdata2 <- expand.nest.data(killdata)
head(killdata2)
print(killdata2[,c(1,2,3,4,7,8,11)], row.names=FALSE)


# Fit a particular model
# This is a model with S a function of nest age
mod.age <-  glm(Survive~NestAge,
         family=binomial(link=logexp(killdata2$Exposure)),
         data=killdata2)
summary(mod.age)

# predict surival as a function of nest age
pred.data <- data.frame(NestAge=1:10)
logit.dsr.pred.age <- predict(mod.age, newdata=pred.data, se.fit=TRUE)  

# put these together in a data frame
dsr.age <- cbind(pred.data, logit.dsr=logit.dsr.pred.age$fit, logit.dsr.se=logit.dsr.pred.age$se.fit)
head(dsr.age)
dsr.age$dsr <- expit(dsr.age$logit.dsr)
dsr.age$dsr.se <- dsr.age$logit.dsr.se* dsr.age$dsr * (1-dsr.age$dsr)
head(dsr.age)

# plot this in the usual way
killdeer.age.le <- ggplot(data=dsr.age, aes(x=NestAge, y=dsr, group=1))+
  ggtitle("DSR as a function of nest age", subtitle="Logistic exposure model")+
  geom_line()+
  geom_ribbon(aes(ymin=expit(logit.dsr-1.96*logit.dsr.se),
                  ymax=expit(logit.dsr+1.96*logit.dsr.se)), alpha=.1)+
  ylim(0.5, 1)
killdeer.age.le
ggsave(killdeer.age.le,
       file=file.path("..","..","..","..","MyStuff","Images","killdeer-age-le.png"),h=4, w=6, units="in", dpi=300)


#_----------------------------------------------------------------------
#--------------------------------------------------------------------------
# Compare to the null model

mod.null <-  glm(Survive~1,
         family=binomial(link=logexp(killdata2$Exposure)),
         data=killdata2)

AICcmodavg::aictab(list(mod.age, mod.null))

