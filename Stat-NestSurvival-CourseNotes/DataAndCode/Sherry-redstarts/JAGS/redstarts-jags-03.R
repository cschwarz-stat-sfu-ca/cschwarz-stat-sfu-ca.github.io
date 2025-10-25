# Nest survival model using JAGS 
# Fit a nest level categorical variable (habitat)

# Analysis of Shelly redstart data using a Bayesian model

# Sherry TW, Wilson S, Hunter S, Holmes RT (2015) 
# Impacts of nest predators and weather on reproductive success and 
# population limitation in a long-distance migratory songbird. 
# Journal of Avian Biology 46(6): 559-569. https://doi.org/10.1111/jav.00536

# Data from
# Sherry TW, Wilson S, Hunter S, Holmes RT (2015) 
# Data from: Impacts of nest predators and weather on reproductive success and 
# population limitation in a long-distance migratory songbird. 
# Dryad Digital Repository. https://doi.org/10.5061/dryad.73870

#
# 2019-06-28 CHJS First Edition
#

library("R2jags")  # used for call to JAGS
library(coda)
library(ggplot2)
library(reshape2)

options(width=200)

source(file.path("..","..","jags-nest-survival-fixed-effects.r"))


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

nestdata <- readxl::read_excel(file.path("..","Sherry.xlsx"), 
                              sheet="NestData")
head(nestdata)


# as noted in the paper, we remove nests from years with no baffeling
dim(nestdata)
nestdata <- nestdata[ !(nestdata$BaffleStatus=="N" & nestdata$Year %in% c(1983, 1984, 1990)),]
dim(nestdata)

# create factor variable for year because it numeric
nestdata$YearF <- factor(nestdata$Year)


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
fe.design <- model.matrix(  ~ DBH+BaffleStatus+YearF, data=nesttime)

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
# evaluated at the covariates for that nest and time

# This is problematic as DBH will vary among the nests.
# We want the DSR at the mean DBH for each baffle status x year combination
# Consequently, we need to create a design matrix for every baffle status x year and
# the mean DBH and then multiply by the estimated betas 

pred.df <- unique(fitted.model$nesttime[, c("BaffleStatus","YearF")])
pred.df$DBH <- mean(nestdata$DBH)

# Notice that the model must MATCH exactly in the order of terms etc as used to generate
# the fe.design matrix seen earlier
pred.design <- model.matrix( ~ DBH+BaffleStatus+YearF, data=pred.df)

head(fe.design)
head(pred.design)

# extract the posterior sample of the beta as a matrix
# the results$BUGSoutput$sims.matrix is a 2-d object [iterations, variables] with chains stacked
# on top of each other
dim(  results$BUGSoutput$sims.matrix)
select <- grepl("^beta", colnames(results$BUGSoutput$sims.matrix))
colnames(results$BUGSoutput$sims.matrix)[select]
beta.matrix <- results$BUGSoutput$sims.matrix[,select]

# create posterior sample of DSR (first on the logit scale) and then back transformed
DSR <- beta.matrix %*% t(pred.design)
DSR <- 1/(1+exp(-DSR))
dim(DSR)

pred.df$mean <- apply(DSR,2,mean)
pred.df$sd   <- apply(DSR,2,sd)
pred.df$lcl  <- apply(DSR,2, quantile, prob=.025)
pred.df$ucl  <- apply(DSR,2,  quantile, prob=.975)

# now we are ready to plot
pred.df$Year <- as.numeric(as.character(pred.df$YearF))
# in this case, we fit a S~Habitat model, so all rows of nesttime have the same estimated DSR

ggplot(data=pred.df, aes(x=Year, y=mean, color=BaffleStatus))+
  ggtitle("Estimated DSR by Baffle Status and Year",
          subtitle="Evaluated at the mean DBH")+
  geom_line( position=position_dodge(w=0.2))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.1, position=position_dodge(w=0.2))+
  ylab("DSR (95% ci)")+
  theme(legend.position=c(0,0),legend.justification=c(0,0))


