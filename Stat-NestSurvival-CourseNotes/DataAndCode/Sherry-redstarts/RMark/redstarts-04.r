# Analysis of Redstart data illustrating basic RMark features
# Reproduce Figure 2 where the Nest survival is computed for each year.

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


# get the yearly date and merge with the nest data
annual.data <- readxl::read_excel(file.path("..","Sherry.xlsx"),
                                  sheet="AnnualCovariates", skip=9)
head(annual.data)

dim(reddata)
reddata <- merge(reddata, annual.data, all.x=TRUE)
dim(reddata)

# any missing data?
sum(!complete.cases(reddata))

# as noted in the paper, we remove nests with baffeling
dim(reddata)
reddata <- reddata[ !(reddata$BaffleStatus=="Y"),]
dim(reddata)

# create factor variable for year
reddata$YearF <- factor(reddata$Year)

# what are the parameters of the model
# There is only one parameter, the daily survival probality (S)
setup.parameters("Nest", check=TRUE)


# Now we look through each year of the study to compute the DSR using a S(.) model
# and compute the nest survival probability
fits <- plyr::dlply(reddata, c("Year","MayTemp"), function(x){
   cat("Fitting ", x$Year[1], x$MayTemp[1], "\n")
  
   # 1. Process the data.
   # The nocc variable is the data at which hatching occurs
   red.proc <- process.data(x, model="Nest", group=c("YearF"), nocc=max(x$LastChecked))
   red.proc


   # 2. Examine and/or modify the ddl. (Not done here)
   red.ddl <- make.design.data(red.proc)
   str(red.ddl)
   red.ddl 
 
   # 3. Fit the dot model 
   fit <- RMark::mark(red.proc, ddl=red.ddl,
                     model="Nest",
                     model.parameters=list(
                       S   =list(formula=~1))
                     )
   print(summary(fit))
   
   # 4. Estimate the nest success for 20 days
   dsr <-  get.real(fit, "S", se=TRUE, vcv=TRUE, expand=TRUE)

   ns    <- prod(dsr$estimates$estimate[1:20])
   ns.se <- deltamethod.special("prod", 
                                dsr$estimates$estimate[1:20],
                                dsr$vcv.real[1:20, 1:20])
   list(red.proc=red.proc, red.ddl=red.ddl, fit=fit, dsr=dsr, ns=ns, ns.se=ns.se)
  
})


# Extract the nest success and create the plot
ns<- plyr::ldply(fits, function(x){
  data.frame(ns=x$ns, ns.se=x$ns.se)
})
ns$lcl <- ns$ns - 1.95*ns$ns.se
ns$ucl <- ns$ns + 1.96*ns$ns.se
ns

ns.vs.temp <- ggplot(data=ns, aes(x=MayTemp, y=ns))+
  ggtitle("Nest success vs. May temp", subtitle="S(.) model fit to each year separately")+
  geom_point()+
  geom_errorbar(aes(ymin=lcl, ymax=ucl), width=0.01)+
  geom_smooth(method="lm")+
  ylim(0,1)
ns.vs.temp
ggsave(ns.vs.temp,
       file=file.path("..","..","..","..","MyStuff","Images","sherry-ns-vs-temp.png"),  h=4, w=6, units="in", dpi=300)

cleanup(ask=FALSE)
