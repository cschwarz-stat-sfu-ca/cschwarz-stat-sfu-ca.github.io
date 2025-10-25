# Jags model for Nest survival using fixed effects and random effects
# The fixed effects are passed to this code using the fe.design matrix.
# A single random effect can be passed using the   re.z1.design matrix.
# If you have more than one random effect, create a re.z2.design matrix and 
# modify the code as needed
#
# 2019-06-29 CJS First edition
#

# This function takes the input data sets, creates some temporary data sets, and then calls JAGS
# Returns the BUGS output and the updated nesttime data frame with the estimated DSR for each 
# combination of nest x day

jags.nest.survival.random.effects <- function(
         nestdata,   # nest data
         nesttime,   # daily nest values with nest, time, nest x time covariates
         fe.design,  # fixed effects design matrix
         re.z1.design,# random effects design matrix.
         init.seed=round(1000000*runif(1)),  # initial seed 
         model.file="model.txt",  # text for the jags file
         debug=FALSE   # if you want to debug the function in case of problems
         ){
  if(debug)browser()

  set.seed(init.seed)

# You must create 2 data frames which are then used to run the JAGS models

#
# nestdata -  Survival intervals collected on the nest data
#
#  This data frame must contain the following variables
#    NestId.num - nest id coded using a number. Usually these are the factor levels of NestId
#                 which was alpha-numeric
#    FirstFound - day number when the nest was first found
#    LastPresent - day number when the nest was last know to be "alive"
#    LastChecked - day nunber of last check on the nest
#       if LastPresent > FirstFound the nest survived for 1 or more days
#       if LastChecked > LastPresent the nest failed somewhere between LastChecked and LastPresent
#       if LastChecked = LastPresent nest survived upto final nest check
# 
#  You can have multiple lines for each nest to account for covariates that change over time
#  It is assumed that freq=1 for all records, i.e. individual nest data rather than grouped data

if(is.null(nestdata)){
   stop("Missing nestdata frame")
}
if(is.null(nestdata$NestId.num)){
   stop("Missing numeric NestId.num variable in nestdata")
}
if(is.null(nestdata$FirstFound)){
   stop("Missing FirstFound variable in nestdata")
}
if(is.null(nestdata$LastPresent)){
   stop("Missing LastPresent variable in nestdata")
}
if(is.null(nestdata$LastChecked)){
   stop("MIssing LastChecked variable in nestdata")
}

# We convert the nest data into two new data frames
#    nest.survive
#    nest.fail
# that contain information on intervals where the nest is known to be alive
# and intervals in which the nest failed

# Remove any records where the intervals are not at least 1 day long

nestdata.survive <- nestdata[,c("NestId.num","FirstFound","LastPresent")]
nestdata.survive <- nestdata.survive[nestdata.survive$LastPresent > nestdata.survive$FirstFound,]
head(nestdata.survive)

nestdata.fail <- nestdata[, c("NestId.num","LastPresent","LastChecked")]
nestdata.fail <- nestdata.fail[ nestdata.fail$LastChecked > nestdata.fail$LastPresent,]
head(nestdata.fail)

# The following variables are then passed to JAGS automatically
#    Nnest.survive - number of nest info records with survival intervals
#    NS.NestId.num - id number of the nest
#    NS.FirstFound - when the nest was first found
#    NS.LastPresent- day when the nest was still alive

# Failure intervals
#   
#    Nnest.fail    - number of nest info records with failure intervals
#    NF.NestId.num: id code of the nest (converted from the NestId alpha numeric code) as a factor
#    NF.LastPresent: last day that a chick was present in the nest
#    NF.LastChecked: last day the nest was checked


# You also need to create the nest x time data matrix which is used to create
# the fixed and random effect design matrices. This matrix should contain
# the NestId.num and Day variables plus any nest level, survey level, or nest x survey level
# covariates.
# 

if(is.null(nesttime)){
   stop("MIssing nesttime data frame")
}
if(is.null(nesttime$NestId.num)){
   stop("Missing NestId.num varible from nesttime data frame")
}
if(is.null(nesttime$Day)){
   stop("missing Day variable from nesttime data frame")
}

if(length(setdiff(nesttime$NestId.num, nestdata$NestId.num))>0){
   stop("NestID.num in nesttime data frame missing fvrom NestID.num in nestdata data frame")
}
if(length(setdiff(nestdata$NestId.num, nesttime$NestId.num ))>0){
  stop("NestID.num in nesdata data frame missing fvrom NestID.num in nesttime data frame")
}

# The following variables are passed to JAGS
#    Nnesttime  - number of nesttime rows
#    NT.NestId.num - nest id for this row of the design matrix
#    NT.Day        - day number for this nest




# You will have created a fixed effects design matrix
if(is.null(fe.design)){
  stop("Missing fixed effects design matrix - fe.design")
}

# You will have created a random effects design matrix
# using somethink like re.z1.design <- model.matrix(  ~ -1 + YearF, data=nesttime)
if(is.null(re.z1.design)){
  stop("Missing random effects design matrix - re.z1.design)")
}




# The model file.
# The cat() command is used to save the model to the working directory.
# Notice that you CANNOT have any " (double quotes) in the bugs code
# between the start and end of the cat("...",) command.

cat(file="model.txt", "
    ############################################################

  var NS.NestId.num [Nnest.survive],
      NS.FirstFound [Nnest.survive],
      NS.LastPresent[Nnest.survive],

      NF.NestId.num [Nnest.fail],
      NF.LastPresent[Nnest.fail],
      NF.LastChecked[Nnest.fail],

      NT.NestId.num[Nnesttime],
      NT.Day       [Nnesttime],
  
      fe.design  [Nnesttime, Nbeta],
      beta       [Nbeta], 
      
      re.z1.design[Nnesttime, Nz1],
      z1          [Nz1]

data {
   for(i in 1:Nnest.fail){
      one.f[i] <- 0  # for the zero's trick
   }
   for(i in 1:Nnest.survive){
      one.s[i] <- 1
   }
}

model {


   # compute the estimated prob(survival) for each nest x day value
   # from the fixed effects design matrix and beta values and random effects design matrix and z1 values
   for(i in 1:Nnesttime){
       logit(S[NT.NestId.num[i], NT.Day[i]]) <- inprod( fe.design   [i, 1:Nbeta],  beta[1:Nbeta]) +
                                                inprod( re.z1.design[i, 1:Nz1],    z1  [1:Nz1])
   }
 
   # distribution of random effects
   for(i in 1:Nz1){
      z1[i] ~ dnorm(0, z1tau)
   }
   
   # prior on the variance of the random effect
   # see http://www.stat.columbia.edu/~gelman/research/published/taumain.pdf
   z1tau <- 1/(z1sd^2)
   z1sd ~ dunif(0,1000)
   
   
   # compute the contribution from each survival portion of the nest record
   for(i in 1:Nnest.survive){
      p1[i] <-    prod( S[NS.NestId.num[i], NS.FirstFound[i]:(NS.LastPresent[i]-1)])
      # use the one's trick
      one.s[i] ~ dbern(p1[i])
   }

   # compute the contribution from each failure portion of the nest record
   for(i in 1:Nnest.fail){
      p2[i] <-    prod( S[NF.NestId.num[i], NF.LastPresent[i]:(NF.LastChecked[i]-1)])
      # use the one's trick
      one.f[i] ~ dbern(p2[i])
   }

   # prior distributions for the betas.
   for(i in 1:Nbeta){
      beta[i] ~ dnorm(0, .001)  # 
   }




    }
    ") # End of the model


# datalist
data.list <-list(
  Nnest.survive   = nrow(nestdata.survive),
  NS.NestId.num   = nestdata.survive$NestId.num,
  NS.FirstFound   = nestdata.survive$FirstFound,
  NS.LastPresent  = nestdata.survive$LastPresent,
  
  Nnest.fail      = nrow(nestdata.fail),
  NF.NestId.num   = nestdata.fail$NestId.num,
  NF.LastPresent  = nestdata.fail$LastPresent,
  NF.LastChecked  = nestdata.fail$LastChecked,
  
  Nnesttime    = nrow(nesttime),
  NT.NestId.num = nesttime$NestId.num,
  NT.Day        = nesttime$Day,
  
  fe.design     = fe.design,
  Nbeta         = ncol(fe.design),
   
  re.z1.design  = re.z1.design,
  Nz1           = ncol(re.z1.design)
)

# check the list
data.list


# Next create the initial values.
# If you are using more than one chain, you need to create a function
# that returns initial values for each chain.

init.list <- list(
  list(beta=rep(0, ncol(fe.design)), z1=rep(.1, ncol(re.z1.design))),
  list(beta=rep(0, ncol(fe.design)), z1=rep(.1, ncol(re.z1.design))),
  list(beta=rep(0, ncol(fe.design)), z1=rep(.1, ncol(re.z1.design)))
)



# Next create the list of parameters to monitor.
# The deviance is automatically monitored.
# 
monitor.list <- c("beta","S", "z1", "z1sd")

results <- jags( 
  data      =data.list,   # list of data variables
  inits     =init.list,   # list/function for initial values
  parameters=monitor.list,# list of parameters to monitor
  model.file="model.txt",  # file with bugs model
  n.chains=3,
  n.iter  =5000,          # total iterations INCLUDING burn in
  n.burnin=2000,          # number of burning iterations
  n.thin=2,               # how much to thin
  DIC=TRUE,               # is DIC to be computed?
  working.dir=getwd()    # store results in current working directory
)


# Extract the S[nest, day] and merge with the nesttime dataframe

names(results)
names(results$BUGSoutput)

row.names(results$BUGSoutput$summary)
select <- grepl("S[", row.names(results$BUGSoutput$summary), fixed=TRUE)
summary.select <- as.data.frame( results$BUGSoutput$summary[select,])
summary.select$NestId.num <- as.numeric( substr(row.names(summary.select),
                                               1+regexpr("[", row.names(summary.select),fixed=TRUE),
                                              -1+regexpr(",", row.names(summary.select),fixed=TRUE)))
summary.select$Day         <- as.numeric( substr(row.names(summary.select),
                                                1+regexpr(",", row.names(summary.select),fixed=TRUE),
                                               -1+regexpr("]", row.names(summary.select),fixed=TRUE)))
names(summary.select)<- make.names(names(summary.select))
dim(nesttime)
nesttime<- merge(nesttime, summary.select, by=c("NestId.num","Day"),all=TRUE)
dim(nesttime)

# return the bugs object and the updated nesttime matrix
list(results=results, nesttime=nesttime, init.seed=init.seed)
}