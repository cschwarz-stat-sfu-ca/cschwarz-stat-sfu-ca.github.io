# logistic exposure model of Shaffer 2004 Auk
# Taken from
# https://rpubs.com/bbolker/logregexp
# https://www.researchgate.net/post/Does_anybody_have_code_for_running_a_logistic_exposure_nest_survival_model_in_R
# http://www.perrywilliams.us/wp-content/uploads/2018/03/Crimmins2016factors.pdf

# define the modification for the logit() link to account for exposure
logexp <- function(exposure = 1)
{
  linkfun <- function(mu) qlogis(mu^(1/exposure))
  ## FIXME: is there some trick we can play here to allow
  ##   evaluation in the context of the 'data' argument?
  linkinv <- function(eta)  plogis(eta)^exposure
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


logit <- function(x){ log(x/(1-x))}
expit <- function(x){1/(1+exp(-x))}


# expand the data to generate the effective sample size
# See the GIM, Chapter 17, page 17-8
# We assume that the variable names on the nest record match those
# required by RMark, but these can be changed
# Returns an expanded dataset with variables Exposure, Survive, Day, NestAge

expand.nest.data <- function( nestdata, 
                              FirstFound="FirstFound", 
                              LastChecked="LastChecked",
                              LastPresent="LastPresent",
                              Fate="Fate",
                              AgeDay1="AgeDay1"){
  require(plyr)
  nestdata2 <- plyr::adply(nestdata, 1, function(x){
    # expand for days when nest is known to be alive
    # create data frame because tibbles can cause problems
    x <- as.data.frame(x)
    #browser()
    orig.x <- x
    x <- x[rep(1,x[,LastPresent]-x[,FirstFound]),]  
    #browser()
    if(nrow(x)>0){
       x$Day <- x[,FirstFound][1]:(x[,LastPresent][1]-1)
       x$Exposure<- 1
       x$Survive <- 1
       #browser()
       if(AgeDay1 %in% names(orig.x)) x$NestAge <- orig.x[,AgeDay1] + x$Day -1 
    }
    
    # now for the final record where the nest fails in interval
    # We define the date of the last interval as the midpoint of the interval
    # We define the nestage in the last interval as the age in the midpoint of the interval as well.
    #browser()
    if(orig.x[,LastChecked] > orig.x[,LastPresent]){
       x2 <- orig.x
       x2$Exposure <- orig.x[,LastChecked] - orig.x[,LastPresent]
       x2$Survive  <- 1-orig.x[,Fate]
       x2$Day      <- x2[,LastPresent] + x2$Exposure/2
       if(AgeDay1 %in% names(orig.x)) x2$NestAge  <- orig.x[,AgeDay1] + x2$Day -1
       x <- rbind(x, x2)
       #browser()
    }
    # remove any records with 0 exposure
    x <- x[ !x$Exposure==0,]
    # return the expanded data set
    x
  })
  nestdata2}

  
