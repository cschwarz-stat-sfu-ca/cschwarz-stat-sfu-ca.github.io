# RMark additional functions

# 2018-02-15 Software contributed by Carl Schwarz (cschwarz.stat.sfu.ca@gmail.com)



##################################################################################
# Model averaging derived parameters
# http://www.phidot.org/forum/viewtopic.php?f=21&t=2865&p=9588&hilit=model+averaging+derived+estimates#p9588


RMark.model.average.derived <- function(model.set, parameter, estimate="estimate", se="se", alpha=0.25){
  # model average derived parameters
  # like RMark, alpha=.025 gives a 95% confidence interval
  # The aic table has all of the model information, but also the aic table at the end
  
  require(plyr)
  require(boot)
  
  # Extract the estimates and standard errors
  est    <- plyr::laply(model.set[-length(model.set)], function(x){x$results$derived[[parameter]][,estimate]})
  est.se <- plyr::laply(model.set[-length(model.set)], function(x){x$results$derived[[parameter]][,se      ]})
  AICc   <- plyr::laply(model.set[-length(model.set)], function(x){x$results$AICc}) 
 
  # Get the model averaged values
  ma.est <- data.frame(model.average(list(estimate=est, se=est.se, AICc=AICc)))
  wald.params <- c("lambda Rate of Change","log odds lambda")
  if (!(parameter %in% wald.params)){
    ma.est$est.logit <- logit(ma.est$estimate)
    ma.est$est.logit.se <- ma.est$se/(ma.est$estimate*(1-ma.est$estimate))
    ma.est$ci.logit.lower <- ma.est$est.logit - qnorm(.975)*ma.est$est.logit.se
    ma.est$ci.logit.upper <- ma.est$est.logit + qnorm(.975)*ma.est$est.logit.se
    ma.est$lcl <- inv.logit(ma.est$ci.logit.lower)
    ma.est$ucl <- inv.logit(ma.est$ci.logit.upper)
    keep.params <- c("estimate","se","lcl","ucl")
    ma.est <- ma.est[,which(names(ma.est) %in% keep.params)]
  } 
  else {
    ma.est$lcl <- ma.est$estimate - qnorm(.975)*ma.est$se
    ma.est$ucl <- ma.est$estimate + qnorm(.975)*ma.est$se 
  }
  ma.est
  
}
