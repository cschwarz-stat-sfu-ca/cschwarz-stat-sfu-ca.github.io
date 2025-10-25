# Analysis of sparrpw data illustrating logistic exposure models
# Incorportating an age variable

# This was taken from 
# https://www.researchgate.net/post/Does_anybody_have_code_for_running_a_logistic_exposure_nest_survival_model_in_R
# which appears to be the data set used by Ben Bolker

# 2019-05-01 CJS Initial code

library(AICcmodavg)
library(ggplot2)
library(plyr)
library(readxl)

source(file.path("..","..","logistic-exposure-model.R"))


# The dataframe must contain the following fields with the following names
#
#    FirstFound: day the nest was first found
#    LastPresent: last day that a chick was present in the nest
#    LastChecked: last day the nest was checked
#    Fate: fate of the nest; 0=hatch an
#    Freq: number of nests with this data
#    AgeDay1: age of nest at day 1 of study (can be negative)
#
# In this example, the multiple visits to a nest have been collapsed
# to a single record for each nest.
# In more complex examples, you may have multple records per nest
# as shown in the mallard example.
#

spardata <- read.csv(file.path("..","sparrow-summary.csv"), header=TRUE, as.is=TRUE, strip.white=TRUE)

spardata <- as.data.frame(spardata)  # sometimes doesn't play nice with tibbles
head(spardata)

spardata2 <- expand.nest.data(spardata)
head(spardata2)


# Fit a particular model
# This is a model with S a function of treeheight
mod.treeht <-  glm(Survive~treeht,
         family=binomial(link=logexp(spardata2$Exposure)),
         data=spardata2)
summary(mod.treeht)


# This is a model with S a function of nest age
# This will differ from the original results because
# we are using the individual day points with individual nest ages
# rather than the nest age at the start of the interval
mod.age.quad <-  glm(Survive~NestAge+I(NestAge^2),
         family=binomial(link=logexp(spardata2$Exposure)),
         data=spardata2)
summary(mod.age.quad)

