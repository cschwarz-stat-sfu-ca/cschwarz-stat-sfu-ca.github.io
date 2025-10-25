# Downloaded from
#   https://www.researchgate.net/post/Does_anybody_have_code_for_running_a_logistic_exposure_nest_survival_model_in_R
#



# Example of R code to fit daily nest survival (DSR) models
# using logistic exposure
# load MASS library if not already
library(MASS)
#####################################################
# Example data: Vesper sparrow in western US.
# nest = Nest ID
# spnum = (VESP) 9
# year = year :)
# int = nest check interval. NOTE! the nest discovery date is not included. 
# t = number of days since last check (related to note above. Can't have zeros)
# fate = 1: success, 0: failed
# nage = nest age (in days)
# sday = season day (from April 15)
# treeht = mean tree height around nest
# sageht = mean sagebrush height around nest
# grassht = mean grass height around nest
exdat <- structure(list(nest = c(1L, 1L, 1L, 1L, 1L, 12L, 12L, 12L, 12L, 16L, 16L, 18L, 18L, 18L, 18L, 18L, 19L, 19L, 19L, 13L, 13L, 13L, 14L, 14L, 15L, 15L, 15L, 17L, 17L, 17L, 17L, 17L), spnum = structure(c(9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L, 9L), .Label = c("BRSP", "BTSP", "DUFL", "GRFL", "GTTO", "SAGS", "sath", "SATH", "VESP" ), class = "factor"), year = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("2006", "2007", "2008", "2009", "2010", "2011", "2013"), class = "factor"), int = c(2L, 3L, 4L, 5L, 6L, 2L, 3L, 4L, 5L, 2L, 3L, 2L, 3L, 4L, 5L, 6L, 2L, 3L, 4L, 2L, 3L, 4L, 2L, 3L, 2L, 3L, 4L, 2L, 3L, 4L, 5L, 6L), t = c(2L, 3L, 3L, 4L, 3L, 3L, 2L, 4L, 4L, 2L, 5L, 2L, 4L, 2L, 2L, 3L, 3L, 5L, 4L, 2L, 5L, 3L, 2L, 5L, 2L, 2L, 3L, 4L, 4L, 3L, 2L, 5L), fate = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L, 1L, 0L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 0L), nage = c(14L, 17L, 20L, 24L, 27L, 15L, 17L, 21L, 25L, 19L, 24L, 16L, 20L, 22L, 24L, 27L, 21L, 26L, 30L, 5L, 10L, 13L, 18L, 23L, 22L, 24L, 27L, 10L, 14L, 17L, 19L, 24L), sday = c(153L, 156L, 159L, 163L, 166L, 171L, 173L, 177L, 181L, 187L, 192L, 200L, 204L, 206L, 208L, 211L, 158L, 163L, 167L, 145L, 150L, 153L, 158L, 163L, 185L, 187L, 190L, 190L, 194L, 197L, 199L, 204L ), treeht = c(0, 0, 0, 0, 0, 44, 44, 44, 44, 72.8, 72.8, 64.3, 64.3, 64.3, 64.3, 64.3, 0, 0, 0, 51.5, 51.5, 51.5, 37, 37, 56.5, 56.5, 56.5, 21, 21, 21, 21, 21), sageht = c(64, 64, 64, 64, 64, 41, 41, 41, 41, 22.3, 22.3, 44.7, 44.7, 44.7, 44.7, 44.7, 25, 25, 25, 36, 36, 36, 80.3, 80.3, 22.2, 22.2, 22.2, 23.3, 23.3, 23.3, 23.3, 23.3), grassht = c(24, 24, 24, 24, 24, 17.9, 17.9, 17.9, 17.9, 26.7, 26.7, 15.9, 15.9, 15.9, 15.9, 15.9, 13.1, 13.1, 13.1, 17.6, 17.6, 17.6, 11.7, 11.7, 19.3, 19.3, 19.3, 19.7, 19.7, 19.7, 19.7, 19.7)), .Names = c("nest", "spnum", "year", "int", "t", "fate", "nage", "sday", "treeht", "sageht", "grassht"), class = "data.frame", row.names = c(NA, 32L))

# Convert to the summary format as required by MARK
# The dataframe must contain the following fields with the following names
#
#    FirstFound: day the nest was first found
#    LastPresent: last day that a chick was present in the nest
#    LastChecked: last day the nest was checked
#    Fate: fate of the nest; 0=hatch an
#    Freq: number of nests with this data


exdat

write.csv(exdat, "sparrow-original.csv", row.names=FALSE)

exdat$Fate <- 1- exdat$Fate  # conver to failure

nest.data <- plyr::ddply(exdat, "nest", function(x){
  # create a nest record for each visit
  x <- x[ order(x$sday),]  # sort by season day
  FirstFound = x$sday[1] - x$t[1]
  LastChecked= max(x$sday)
  LastPresent= x$sday[ nrow(x)] - as.numeric(x$t[nrow(x)]*(x$fate[nrow(x)]==0)) # back out a failure
  AgeDay1 = min(x$nage)-min(x$sday)+1
  #browser()
  data.frame(FirstFound=FirstFound,
             LastPresent=LastPresent,
             LastChecked=LastChecked,
             Fate=1-min(x$fate),
             Freq=1,
             AgeDay1=AgeDay1,
             treeht = x$treeht[1],
             sageht = x$sageht[1],
             grassht= x$grassht[1])
})
nest.data
write.csv(nest.data, "sparrow-summary.csv", row.names=FALSE)


#####################################################
# Next, we need to source a custom link function. Courtesy of Ben Bolker.
# Note that this version of the logexp link function will allow you to use
# lme4 (glmer) to fit mixed effect models as well. The example below just # uses glm to fit models of all fixed effects.
logexp <- function(exposure = 1) {
linkfun <- function(mu) qlogis(mu^(1/exposure))
## FIXME: is there some trick we can play here to allow
## evaluation in the context of the 'data' argument?
linkinv <- function(eta) plogis(eta)^exposure
logit_mu_eta <- function(eta) {
ifelse(abs(eta)>30,.Machine$double.eps,
exp(eta)/(1+exp(eta))^2)
## OR .Call(stats:::C_logit_mu_eta, eta, PACKAGE = "stats")
}
mu.eta <- function(eta) {
exposure * plogis(eta)^(exposure-1) *
logit_mu_eta(eta)
}
valideta <- function(eta) TRUE
link <- paste("logexp(", deparse(substitute(exposure)), ")",
sep="")
structure(list(linkfun = linkfun, linkinv = linkinv,
mu.eta = mu.eta, valideta = valideta,
name = link),
class = "link-glm")
}
# Now fit models... ####################################################### Fit a model of DSR as a function of tree height
mfit1 <- glm(fate ~ treeht, family=binomial(link=logexp(exdat$t)), data=exdat)
summary(mfit1)

# Effect not significant with sample data
####################################################### Fit a model of DSR as a quadratic function of nest age
mfit2 <- glm(fate ~ nage + I(nage^2),
family=binomial(link=logexp(exdat$t)), data=exdat)
summary(mfit2)
# These effects also not significant with sample data
####################################################### Have fun!