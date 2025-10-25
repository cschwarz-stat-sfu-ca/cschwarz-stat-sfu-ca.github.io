# Nest survival model using JAGS 
# Fit a nest level continuous variable (robel height)

#
# 2019-06-28 CHJS First Edition
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

nestdata <- readxl::read_excel(file.path("..","mallard.xlsx"), 
                               sheet="mallard")
# we need to add a NestId variables
nestdata$NestId <- paste("Nest", 1:nrow(nestdata),sep=".")
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

# none added here

# if there is a AgeDay1 variable, we compute the nest age for each time for each nest
if( !is.null(nesttime$AgeDay1)){
   nesttime$NestAge <- nesttime$AgeDay1 + nesttime$Day -1
}
head(nesttime)


# Add any next x day survey covariates to the nesttime data
#

# there is nothing here for this example


# Set up the design matrix for the fixed effects
fe.design <- model.matrix(  ~ Robel, data=nesttime)

head(fe.design)


# Finally, the actual call to JAGS
fitted.model <- jags.nest.survival.fixed.effects(
         nestdata=nestdata,    # nest data
         nesttime=nesttime,    # daily nest values with nest, time, nest x time covariates
         fe.design=fe.design,  # fixed effects design matrix
         init.seed=12321312)   # initial seed)  

# the results list has lots of other stuff
results <- fitted.model$results

# the nesttime dataframe has the estimated DSR for every combination of NestId.num and Day

# in this case, we fit a S~Robel model, so all rows of nesttime have the same estimated DSR
plotdata <- plyr::ddply(fitted.model$nesttime, "Robel", function(x){ x[1,]})
head(plotdata)
ggplot(data=plotdata, aes(x=Robel, y=mean))+
  ggtitle("Estimated DSR as a function of Robel height")+
  geom_line(group=1)+
  geom_ribbon(aes(ymin=X2.5., ymax=X97.5.), alpha=0.1)+
  ylab("DSR (95% ci)")


