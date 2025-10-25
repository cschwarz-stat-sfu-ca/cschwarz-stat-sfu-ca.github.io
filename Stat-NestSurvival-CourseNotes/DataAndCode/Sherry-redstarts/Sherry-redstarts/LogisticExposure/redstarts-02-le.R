# Analysis of Shelly redstart data using the logistic exposure model
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
#
# In this example, the multiple visits to a nest have been collapsed
# to a single record for each nest.
# In more complex examples, you may have multple records per nest
# as shown in the mallard example.
#

reddata <- readxl::read_excel(file.path("..","Sherry.xlsx"), 
                               sheet="NestData")
head(reddata)

# We expand the data to generate the effective sample size
reddata2 <- expand.nest.data(reddata)
head(reddata2)
rm(reddata)


# as noted in the paper, we remove nests from years with no baffeling
dim(reddata2)
reddata2 <- reddata2[ !(reddata2$BaffleStatus=="N" & reddata2$Year %in% c(1983, 1984, 1990)),]
dim(reddata2)

# create factor variable for year
reddata2$YearF <- factor(reddata2$Year)

# create factor for baffle
reddata2$BaffleStatus <- factor(reddata2$BaffleStatus)

#--------------------------------------------
#  Set up the set of model to fit
model.list.csv <- textConnection(
" S
 ~DBH
 ~DBH+YearF
 ~DBH+YearF+BaffleStatus
")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

model.fits <- plyr::dlply(model.list, c("S","model.number"), function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- glm(formula=as.formula(paste("Survive", eval(x$S))),
            family=binomial(link=logexp(input.data$Exposure)),
            data=input.data)
  fit
  
},input.data=reddata2)

# examine results for a particular model
model.number <-1
names(model.fits[[model.number]])
model.fits[[model.number]]$formula

summary(model.fits[[model.number]])


# Model comparision and averaging
# collect models and make AICc table
AICcmodavg::aictab(model.fits)


#----------------------------------------------------------------------


# compute the DSR and model average the daily values
# we need to set up the prediction matrix with variables that match the models
mean.dbh= mean(reddata2$DBH)
mean.dbh

pred.data      <- expand.grid(DBH=mean.dbh, 
                             YearF=as.factor(unique(as.numeric(as.character(reddata2$YearF)))),
                             BaffleStatus=unique(as.character(reddata2$BaffleStatus)))
pred.data                            

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

dsr.ma <- plyr::ddply(dsr.indiv, c("YearF","BaffleStatus"), function(x, model.info){
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


# plot of bafflestatus by year (model averages())
final.est <- ggplot(data=plotdata.ma, aes(x=as.numeric(as.character(YearF)), y=dsr, color=BaffleStatus, shape=BaffleStatus))+
  ggtitle("Effects of baffle status by year", subtitle="Logistic exposure model")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  geom_errorbar(aes(ymin=dsr.lcl, ymax=dsr.ucl),width=.1, position=position_dodge(w=0.2))+
  ylab("DSR (95% ci)")+xlab("Year")+
  scale_x_continuous(breaks=1980:2020)+
  theme(legend.position=c(0,0), legend.justification=c(0,0))
final.est
ggsave(final.est,
       file=file.path("..","..","..","..","MyStuff","Images","sherry-baffle-year-le.png"), h=4, w=6, units="in", dpi=300)


#-------------------------------------------------------------
# DSR by DBH for 1995 N
pred.data      <- expand.grid(DBH=seq(min(reddata2$DBH),max(reddata2$DBH), length.out=50), 
                             YearF=as.factor(1985),
                             BaffleStatus="N")
pred.data                            

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

dsr.ma <- plyr::ddply(dsr.indiv, c("DBH"), function(x, model.info){
   # merge the model information with the estimates
   x <-merge(x, model.info)
   # get the model averaged values
   #browser()
   ma <- AICcmodavg::modavgCustom(x$logL, x$K, modnames=x$S, nobs=x$nobs,
                                  estimate=x$dsr, se=x$dsr.se)
   data.frame(dsr=ma$Mod.avg.est, dsr.se=ma$Uncond.SE, dsr.lcl=ma$Lower.CL, dsr.ucl=ma$Upper.CL)
},model.info=model.info)
head(dsr.ma)

dbh.effect <- ggplot(data=dsr.ma, aes(x=DBH, y=dsr))+
  ggtitle("Effect of DBH for 1985 N", subtitle="Logistic exposure model")+
  geom_line()+
  geom_ribbon(aes(ymin=dsr.lcl, ymax=dsr.ucl), alpha=0.2)

dbh.effect
ggsave(dbh.effect, 
       file=file.path("..","..","..","..","MyStuff","Images","sherry-dbh-le.png"), h=4, w=6, units="in", dpi=300)

