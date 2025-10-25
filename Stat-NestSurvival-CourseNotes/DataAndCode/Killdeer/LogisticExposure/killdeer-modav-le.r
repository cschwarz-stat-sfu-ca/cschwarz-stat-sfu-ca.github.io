# Analysis of killdeer data illustrating model averaging
# and logistic exposure models

# 2019-05-01 CJS Initial version

# This is the killdeer data that ships with RMark

library(AICcmodavg)
library(ggplot2)
library(plyr)
library(readxl)

source(file.path("..","..","logistic-exposure-model.R"))

# The dataframe must contain the following fields with the following names
#
#   FirstFound: day the nest was first found
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

# We expand the data to generate the effective sample size
killdata2 <- expand.nest.data(killdata)
head(killdata2)

# Add the Day2 term and the first/second half variables
killdata2$Day2 <- (killdata2$Day-20)^2
killdata2$studyhalf <- car::recode(killdata2$Day,
                                  " lo:20='1st'; 21:hi='2nd'")

#  Set up the set of model to fit
model.list.csv <- textConnection(
" S
 ~1
 ~Day
 ~Day+Day2
 ~studyhalf
")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

model.fits <- plyr::dlply(model.list, c("S","model.number"), function(x,input.data){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- glm(formula=as.formula(paste("Survive", eval(x$S))),
            family=binomial(link=logexp(input.data$Exposure)),
            data=input.data)
  fit
  
},input.data=killdata2)


# examine individual model results
model.number <-1
names(model.fits[[model.number]])
model.fits[[model.number]]$formula

summary(model.fits[[model.number]])


# Model comparision and averaging
# collect models and make AICc table
AICcmodavg::aictab(model.fits)


# compute the DSR and model average the daily values
# we need to set up the prediction matrix with variables that match the models

pred.data      <- data.frame(Day=1:39)
pred.data$Day2 <- (pred.data$Day-20)^2  # we need to match coding above
pred.data$studyhalf <- car::recode(pred.data$Day,
                                  " lo:20='1st'; 21:hi='2nd'")

dsr.indiv <- plyr::ldply(model.fits, function(x,pred.data){
   # get the predictions on the logit scale and then back transform
    logit.dsr.pred <- predict(x, newdata=pred.data, se.fit=TRUE)  

    # put these together in a data frame
    dsr <- cbind(pred.data, logit.dsr=logit.dsr.pred$fit, logit.dsr.se=logit.dsr.pred$se.fit)
    dsr$dsr   <- expit(dsr$logit.dsr)
    dsr$dsr.se <- dsr$logit.dsr.se* dsr$dsr * (1-dsr$dsr)
    dsr
},pred.data=pred.data)

plotdata.indiv <- dsr.indiv
head(plotdata.indiv)

# do the model averaging for each day's DSR
# Extract the logL and K number of parameters for each model in the model set.
# The model number is used to index the values
model.info <- plyr::ldply(model.fits, function(x){
    #browser()
    logL <- logLik(x)
    K    <- length(coef(x))
    nobs <- nrow(x$data)
    data.frame(logL=logL, K=K, nobs=nobs)
})
model.info

dsr.ma <- plyr::ddply(dsr.indiv, c("Day"), function(x, model.info){
   # merge the model information with the estimates
   x <-merge(x, model.info)
   # get the model averaged values
   #browser()
   ma <- AICcmodavg::modavgCustom(x$logL, x$K, modnames=x$S, nobs=x$nobs,
                                  estimate=x$dsr, se=x$dsr.se)
   data.frame(dsr=ma$Mod.avg.est, dsr.se=ma$Uncond.SE, dsr.lcl=ma$Lower.CL, dsr.ucl=ma$Upper.CL)
},model.info=model.info)
head(dsr.ma)

plotdata.ma  <- dsr.ma
plotdata.ma$S <- "Model avg"

# model average the predicted nest suvival probability
final.fit <- ggplot(data=plotdata.indiv, aes(x=Day, y=dsr, color=S))+
  ggtitle("Comparing the model fits", subtitle="Fit using logistic exposure models")+
  geom_line(aes(group=S))+
  geom_line(data=plotdata.ma, color='black', size=2, group=1)+
  theme(legend.justification=c(0,0), legend.position=c(0,0))
final.fit
ggsave(final.fit,
       file=file.path("..","..","..","..","MyStuff","Images","killdeer-modavg-le.png"),h=4, w=6, units="in", dpi=300)


