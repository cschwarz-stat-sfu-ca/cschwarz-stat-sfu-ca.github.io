# Nest survival model using JAGS 
# Only fixed effects at this time
#
# 2019-06-30 CHJS First Edition
#

library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)

options(width=200)

source(file.path("..","..","jags-nest-survival-fixed-effects.r"))


# The input dataframe must contain the following fields with the following names
#
#    NestID: id code of the nest (alpha numeric)
#    FirstFound: day the nest was first found
#    LastPresent: last day that a chick was present in the nest
#    LastChecked: last day the nest was checked
#    Fate: fate of the nest; 0 = success; 1=fail
#    AgheDay1 = age of the nest on day 1 (if you are fitting age of nest models) 
#
# You could also have a nest level covariates, survey level covariates, and
# next x survey time covariates as well

nestdata <- readxl::read_excel(file.path("..","Killdeer.xlsx"), 
                               sheet="killdeer-age")
nestdata <- plyr::rename(nestdata, c("id"="NestId"))
head(nestdata)

# Unfortunately, JAGS cannot deal with alpha numeric code and 
# so we need to convert the alphanumberic NestID to numeric codes
# by declaring NestId as a factor and extracting the level values
nestdata$NestId.num <- as.numeric(factor(nestdata$NestId))

# We must create a file with every combination of next x day nest was "active"
# being every day from FirstCound to LastChecked-1

nesttime <- plyr::adply(nestdata, 1, function(x){
     nesttime <- expand.grid(NestId.num=x$NestId.num, 
                             Day=x$FirstFound:(x$LastChecked-1),
                             stringsAsFactors=FALSE)
     nesttime
})


# Extract the nest level covariates (including AgeNest1)
# The next level covariates should be indexed using NestId
# If AgeNest1 variable is present then the age of the nest is computed
#
nest.covariates <- NULL

if( !is.null(nest.covariates)){
   nesttime <- merge(nesttime, nest.covariates, by="NestId")
}


# Extract any survey time covariates such as time, time^2, early/late
# weather covariates ect.
# All of these covariates will affect all nests simultaneouls
nesttime $Day2 <- (nesttime$Day-20)^2  # day^2 for quadratic trends
nesttime $Period <- car::recode(nesttime$Day,
                    paste("lo:", (max(nesttime$Day)+min(nesttime$Day))/2, "='Early';",
                          "else='Late'"))
xtabs(~Period+Day, data=nesttime, exclude=NULL, na.action=na.pass)


# if there is a AgeDay1 variable, we compute the nest age for each time for each nest
if( !is.null(nesttime$AgeDay1)){
   nesttime$NestAge <- nesttime$AgeDay1 + nesttime$Day -1
}
head(nesttime)


# Add any next x day survey covariates to the nesttime data
#

# there is nothing here for this example


# Set up the design matrix for the fixed effects
fe.design <- model.matrix(  ~ 1, data=nesttime)

head(fe.design)


# Finally, load and run the JAGS model
fitted.model <- jags.nest.survival.fixed.effects(
         nestdata=nestdata,    # nest data
         nesttime=nesttime,    # daily nest values with nest, time, nest x time covariates
         fe.design=fe.design,  # fixed effects design matrix
         init.seed=12321312)   # initial seed)  

# the nesttime dataframe has the estimated DSR for every combination of NestId.num and Day

# in this case, we fit a S~1 model, so all rows of nesttime have the same estimated DSR
fitted.model$nesttime[1,]



# the results list has lots of other stuff
results <- fitted.model$results
names(results)
names(results$BUGSoutput)

# we can also look at the beta estimates
# in this case this is the logit DSR which is the same for all nest x days
results$BUGSoutput$summary[ grepl("beta", row.names(results$BUGSoutput$summary)),,drop=FALSE]



#######################################

# get the full summary table
results$BUGSoutput$summary
results$BUGSoutput$summary[grepl("beta",rownames(results$BUGSoutput$summary)),
                           c("mean", "sd", "2.5%","97.5%","Rhat", "n.eff")]

# get just the means
results$BUGSoutput$mean
results$BUGSoutput$mean$parm

# the results$BUGSoutput$sims.array is a 3-d object [iterations, chains, variables]
# that can be used for posterior plots, diagnostics etc
parm <- "beta"
dim(results$BUGSoutput$sims.array)
#results$BUGSoutput$sims.array[1:5,,]
results$BUGSoutput$sims.array[1:5,1,parm,drop=FALSE]


# the results$BUGSoutput$sims.matrix is a 2-d object [iterations, variables] with chains stacked
# on top of each other
dim(results$BUGSoutput$sims.matrix)
#results$BUGSoutput$sims.matrix[1:5,]
results$BUGSoutput$sims.matrix[1:5,parm,drop=FALSE]


# make a posterior density plot of the parameter
plotdata <- data.frame(parm=results$BUGSoutput$sims.matrix[,parm], stringsAsFactors=FALSE)
head(plotdata)
postplot.parm <- ggplot2::ggplot( data=plotdata, aes(x=parm, y=..density..))+
  geom_histogram(alpha=0.3)+
  geom_density()+
  ggtitle(paste("Posterior density plot for ", parm))
postplot.parm
#ggsave(plot=postplot.parm, file='...posterior.png', h=4, w=6, units="in", dpi=300)



# make a trace plot (notice we use the sims.array here)
plotdata <- data.frame(results$BUGSoutput$sims.array[,,parm,drop=FALSE], stringsAsFactors=FALSE)
plotdata$iteration <- 1:nrow(plotdata)
head(plotdata)

# convert from wide to long format
plotdata2 <- reshape2::melt(data=plotdata, 
                            id.vars="iteration",
                            measure.vars=paste("X",1:results$BUGSoutput$n.chains,".",parm,sep=""),
                            variable.name="chain",
                            value.name=parm)
head(plotdata2)
traceplot.parm <- ggplot2::ggplot(data=plotdata2, aes_string(x="iteration", y=parm, color="chain"))+
  ggtitle("Trace plot")+
  geom_line(alpha=.2)
traceplot.parm
#ggsave(plot=traceplot.parm, file='....-trace-parm.png', h=4, w=6, units="in", dpi=300)


# autocorrelation plot
# First compute the autocorrelation plot
acf.parm <-acf( results$BUGSoutput$sims.matrix[,parm,drop=FALSE], plot=FALSE)
acf.parm
acfplot.parm <- ggplot(data=with(acf.parm, data.frame(lag, acf)), aes(x = lag, y = acf)) +
  ggtitle(paste("Autocorrelation plot for ",parm))+
  geom_hline(aes(yintercept = 0)) +
  geom_segment(aes(xend = lag, yend = 0))
acfplot.parm
#ggsave(plot=acfplot.parm, file="...-acf-parm.png",h=4, w=6, units="in", dpi=300)


# You can save the results information to an r-object for retreival later using a load# command

#save(list=c("results"), file="results.Rdata")
#load(file="results.Rdata")

