# Analysis of Redstart data illustrating basic RMark features
# Reproduce Table 1a.

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


library(ggplot2)
library(readxl)
library(RMark)


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
reddata <- as.data.frame(reddata)

# as noted in the paper, we remove nests from years with no baffeling
dim(reddata)
reddata <- reddata[ !(reddata$BaffleStatus=="N" & reddata$Year %in% c(1983, 1984, 1990)),]
dim(reddata)

# create factor variable for year
reddata$YearF <- factor(reddata$Year)

# create factor for baffle
reddata$BaffleStatus <- factor(reddata$BaffleStatus)

# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
red.proc <- process.data(reddata, model="Nest", group=c("YearF","BaffleStatus"), nocc=max(reddata$LastChecked))
red.proc


# 2. Examine and/or modify the ddl. (Not done here)
red.ddl <- make.design.data(red.proc)
str(red.ddl)
red.ddl 
  

# 3. Set up the set of model to fit
model.list.csv <- textConnection(
  " S
  ~1+DBH
  ~1+DBH+YearF
  ~1+DBH+YearF+BaffleStatus
  ")

model.list <- read.csv(model.list.csv, header=TRUE, as.is=TRUE, strip.white=TRUE)
model.list$model.number <- 1:nrow(model.list)
model.list

# fit the models
myobs <- ls()
myobs <- myobs[ grepl("m...",myobs,fixed=TRUE)]
cat("Removing ", myobs, "\n")
rm(list=myobs)

model.fits <- plyr::dlply(model.list, "model.number", function(x,input.data, input.ddl){
  cat("\n\n***** Starting ", unlist(x), "\n")
  
  fit <- RMark::mark(input.data, ddl=input.ddl,
                     model="Nest",
                     model.parameters=list(
                       S   =list(formula=as.formula(eval(x$S)))
                     )
                     #,brief=TRUE,output=FALSE, delete=TRUE
                     #,invisible=TRUE,output=TRUE  # set for debugging
  )
  mnumber <- paste("m...",formatC(x$model.number, width = 3, format = "d", flag = "0"),sep="")
  assign( mnumber, fit, envir=.GlobalEnv)
  #browser()
  fit
  
},input.data=red.proc, input.ddl=red.ddl)


# examine individula model results
model.number <-1
names(model.fits[[model.number]])
model.fits[[model.number]]$model.name

summary(model.fits[[model.number]])
model.fits[[model.number]]$results$real
model.fits[[model.number]]$results$beta
model.fits[[model.number]]$results$derived

get.real(model.fits[[model.number]], "S", se=TRUE)


# Model comparision and averaging
# collect models and make AICc table

model.set <- RMark::collect.models( type="Nest")
model.set

names(model.set)
model.set$model.table


# model averaged values
get.real(model.set[[1]], "S", se=TRUE)
get.real(model.set[[2]], "S", se=TRUE)
get.real(model.set[[3]], "S", se=TRUE)

S.ma <- RMark::model.average(model.set, param="S")
head(S.ma)

# plot the DSR without baffles for each year
# because there are no time trends we can use the estimates at time=1
# This are at the average value of dbh
plotdata <- RMark::model.average(model.set, param="S", vcv=TRUE)$estimates
plotdata <- plotdata[ plotdata$time==1,]
plotdata

final.est <- ggplot(data=plotdata, aes(x=as.numeric(as.character(YearF)), y=estimate, color=BaffleStatus, shape=BaffleStatus))+
  ggtitle("Effects of baffle status by year")+
  geom_point(position=position_dodge(w=0.2))+
  geom_line(position=position_dodge(w=0.2))+
  geom_errorbar(aes(ymin=lcl, ymax=ucl),width=.1, position=position_dodge(w=0.2))+
  ylab("DSR (95% ci)")+xlab("Year")+
  scale_x_continuous(breaks=1980:2020)+
  theme(legend.position=c(0,0), legend.justification=c(0,0))
final.est
ggsave(final.est,
       file=file.path("..","..","..","..","MyStuff","Images","sherry-baffle-year.png"), h=4, w=6, units="in", dpi=300)


# impact of dbh
head(get.real(model.set[[1]], param="S", se=TRUE))

pred.data <- data.frame(DBH=seq(min(reddata$DBH), max(reddata$DBH), length.out=50),
                        BaffleStatus="N",
                        YearF=1985,
                        index=1)
plotdata <- covariate.predictions(model.set, data=pred.data)$estimates
head(plotdata)
tail(plotdata)

dbh.effect <- ggplot(data=plotdata, aes(x=DBH, y=estimate))+
  ggtitle("Effect of DBH for 1985 N")+
  geom_line()+
  geom_ribbon(aes(ymin=lcl, ymax=ucl), alpha=0.2)

dbh.effect
ggsave(dbh.effect, 
       file=file.path("..","..","..","..","MyStuff","Images","sherry-dbh.png"), h=4, w=6, units="in", dpi=300)

cleanup(ask=FALSE)
