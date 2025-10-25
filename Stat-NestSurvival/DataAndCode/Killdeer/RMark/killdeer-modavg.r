# Analysis of killdeer data illustrating model averaging

# 2019-05-01 CJS Initial version

# This is the killdeer data that ships with RMark
# which has been saved as a *.csv file

library(ggplot2)
library(RMark)

source(file.path("..","..","RMark.additional.functions.r"))

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


# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# 1. Process the data.
# The nocc variable is the data at which hatching occurs
kill.proc <- process.data(killdata, model="Nest", nocc=40)
kill.proc


# 2. Examine and/or modify the ddl. Here you could standardize covariates
kill.ddl <- make.design.data(kill.proc)
kill.ddl$S$Time2 <- (kill.ddl$S$Time-20)^2
kill.ddl 
 

# 3. Set up the set of model to fit
model.list.csv <- textConnection(
" S
 ~1
 ~Time
 ~Time+Time2
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
  
},input.data=kill.proc, input.ddl=kill.ddl)


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
S.ma




# plot the estimates from the model plus the model averaging
length(model.set)
plotdata.indiv <- plyr::ldply(model.set[-length(model.set)], function(x){
   est <- get.real(x, "S", se=TRUE)
   est$Model <- x$model.name
   est
})

plotdata.ma  <- S.ma
plotdata.ma$Model <- "Model avg"


# model average the predicted nest suvival probability
final.fit <- ggplot(data=plotdata.indiv, aes(x=time, y=estimate, color=Model))+
  ggtitle("Comparing the model fits")+
  geom_line(aes(group=Model))+
  geom_line(data=plotdata.ma, color='black', size=2, group=1)+
  theme(legend.justification=c(0,0), legend.position=c(0,0))
final.fit
ggsave(final.fit,
       file=file.path("..","..","..","..","MyStuff","Images","killdeer-modavg.png"),h=4, w=6, units="in", dpi=300)



# model average the derived parameter of the oveall nest survival
model.fits[[1]]$results$derived$"S Overall Survival"
model.fits[[2]]$results$derived$"S Overall Survival"
model.fits[[3]]$results$derived$"S Overall Survival"

RMark.model.average.derived(model.set, param="S Overall Survival")




# cleanup
cleanup(ask=FALSE)


