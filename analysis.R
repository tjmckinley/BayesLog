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

#run mcmc
bwt.mcmc <- run.mcmc(bwt, "low")
#bwt.mcmc1 <- run.mcmc(bwt, "low", varselect = T)

