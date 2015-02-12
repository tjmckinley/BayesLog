#load data
library(MASS)
bwt <- with(birthwt,
{
    race <- factor(race, labels = c("white", "black", "other"))
    ptd <- factor(ptl > 0)
    ftv <- factor(ftv)
    levels(ftv)[-(1:2)] <- "2+"
    data.frame(low = factor(low), age, lwt, race, smoke = factor(smoke > 0),
        ptd, ht = factor(ht > 0), ui = factor(ui > 0), ftv)
})

#source in necessary functions
source("analysisfunctions.R")

#convert data frame into correct format for use
#in Bayesian model
bwt <- createLinear(bwt, "low")

#run mcmc
bwt.mcmc <- run.mcmc(bwt, "low")

