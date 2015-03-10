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
#bwt.mcmc <- run.mcmc(bwt, "low")
bwt.mcmc1 <- run.mcmc(bwt, "low", varselect = F, random = "locrand", nchains = 1)

summary.varselect(bwt.mcmc1)

#check against BMA package
library(BMA)
bwt.bma <- bic.glm(low ~ ., data = bwt, strict = FALSE, OR = 1000, glm.family = "binomial", factor.type = TRUE)

#check against stepAIC
bwt.aic <- glm(low ~ ., data = bwt, family = "binomial")
bwt.aic <- stepAIC(bwt.aic, scope = list(upper = ~ ., lower = ~ 1))
