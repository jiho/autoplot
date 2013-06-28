#
#     Simple regression example for mcmc objects
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

## Generate data

set.seed(12345)
n <- 50                         # nb of points
x <- runif(n, min=0, max=100)   # x coordinate
# random deviation around a linear relationship
a <- 2
b <- 10
y <- a * x + b + rnorm(n, mean=0, sd=10)

dd <- data.frame(x, y)
dl <- list(n=n, x=x, y=y)

library("ggplot2")
ggplot(dd) + geom_point(aes(x=x, y=y))


## Regression model with a frequentist approach
mF <- lm(y ~ x)
summary(mF)
plot(mF, 1:2)
shapiro.test(mF$residuals)
ggplot(dd) + geom_point(aes(x=x, y=y)) + geom_path(aes(x=x, y=predict(mF)), colour="red")


## Regression model in JAGS
library("rjags")            # interface R with JAGS
library("coda")             # extract data from the BUGS run

# compile and adapt the model in JAGS
mJ <- jags.model("demo-mcmc.bug", data=dl, n.chains=3)
# additional burn-in
update(mJ, 9000)
# run the model while recording some data
cJ <- coda.samples(mJ, c("a", "b", "sigmaR"), n.iter=10000, thin=15)

summary(cJ)

# default plot : densities
source("autoplot-mcmc.R")
autoplot(cJ)

# specify plots
autoplot(cJ, plot="trace")
autoplot(cJ, plot="autocor")
autoplot(cJ, plot="dens")     # abbreviations are OK
autoplot(cJ, plot="histo")    # but they must be inambiguous
autoplot(cJ, plot="histogram")
autoplot(cJ, plot="history")


## Regression model in WinBUGS
library("R2WinBUGS")

mW <- bugs(data=d, inits=NULL, parameters.to.save=c("a", "b", "sigmaR"),
           model.file="demo-mcmc.bug", n.chains=3, n.iter=20000, n.burnin=10000, n.thin=15)
# NB: the arguments "bugs.directory", "program", "working.directory", and optionally everything related to WINE should be set according to your own setup
cW <- as.mcmc.list(mW)

summary(cW)

# default plot : densities
source("autoplot-mcmc.R")
autoplot(cW)

# specify plots
autoplot(cW, plot="trace")
autoplot(cW, plot="autocor")
autoplot(cW, plot="dens")     # abbreviations are OK
autoplot(cW, plot="histo")    # but they must be inambiguous
autoplot(cW, plot="histogram")
autoplot(cW, plot="history")
