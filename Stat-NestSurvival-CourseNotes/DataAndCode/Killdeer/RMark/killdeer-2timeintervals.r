# Analysis of killdeer data illustrating basic RMark features
# Fit a model with DSR differing in first and second half of the study

# 2019-05-01 CJS Code produced.

# This is the killdeer data that ships with RMark
# which has been saved as a *.csv file

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


# 2. Examine and/or modify the ddl. 
# Here we create a variable for first/second half of the study
kill.ddl <- make.design.data(kill.proc)

kill.ddl$S$studyhalf <- car::recode(kill.ddl$S$Time,
                                  " lo:20='1st'; 21:hi='2nd'")
str(kill.ddl)
kill.ddl 
  

# 3. Fit a particular model
# This is a model with S linear over time.
# Notice the use of Time vs time.
mod.res <-  RMark::mark(kill.proc, ddl=kill.ddl,
                          model="Nest",
                          model.parameters=list(
                            S   =list(formula=~studyhalf)
                          )
  )
summary(mod.res)


# Look the objects returned in more details
names(mod.res)
names(mod.res$results)

# look at estimates on beta and original scale
mod.res$results$beta  # on the logit scale

mod.res$results$real# on the regular 0-1 scale for each site

# derived variabldes is the nest survival probability over the 40 (nocc) days 
names(mod.res$results$derived)

mod.res$results$derived$"S Overall Survival"



# alternatively
get.real(mod.res, "S", se=TRUE)


# Notice that the nest survival is the product of the individual S
prod(get.real(mod.res, "S", se=TRUE)$estimate)

mod.res$design.matrix

# plot the trend over time
plotdata <- get.real(mod.res, "S", se=TRUE)

plot.res <- ggplot(data=plotdata, aes(x=time, y=estimate, group=1))+
  ggtitle("Model with two groups of S")+
  geom_line()+
  geom_ribbon(aes(x=time, ymin=lcl, ymax=ucl), alpha=0.2)+
  #geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
  ylim(0,1)
plot.res
ggsave(plot.res,
       file=file.path("..","..","..","..","MyStuff","Images","killdeer-2groups-S.png"), h=4, w=6, units="in", dpi=300)


