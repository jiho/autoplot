#
#     Data extraction (fortify)
#     and plotting in ggplot2 (autplot)
#     for mcmc objects
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

# Extract the data from an mcmc object into a data.frame
fortify.mcmc <- function(x) {
    # extract attributes, inluding iteraton number
    mcpar <- attr(x, "mcpar")
    # convert into data.frame
    x <- as.data.frame(x)
    # add iteration number
    x$iter <- seq(mcpar[1], mcpar[2], mcpar[3])
    return(x)
}

# Extract the data from each mcmc element of an mcmc.list and store it in a single data.frame
fortify.mcmc.list <- function(x) {
    # extract each element of the list
    for (i in 1:length(x)) {
        x[[i]] <- fortify(x[[i]])
        # identify the chain
        x[[i]]$chain <- i
    }
    # concatenate them all in a single data.frame
    x <- do.call("rbind", x)
    x$chain <- factor(x$chain)
    return(x)
}

# Generic (to be included in ggplot2 eventually)
autoplot <- function(X, ...) { UseMethod("autoplot")}

# Pre-made plots for the results of an mcmc analysis (mcmc.list class)
autoplot.mcmc.list <- function(x, plot=c("density", "histogram", "trace", "history", "autocor", "quantiles")) {

    # check type of plot
    plot <- match.arg(plot)

    # convert data to a data.frame
    xd <- fortify(x)

    # convert it into the tall format for plotting
    require("reshape")
    xM <- melt(xd, id.vars=c("iter", "chain"))

    # detect when there are several chains, to colour them differently
    severalChains <- (length(unique(xd$chain)) > 1)

    # Density distribution of recorded variables (estimated by kernel)
    if (plot == "density") {
        p <- ggplot(xM) + geom_area(aes(x=value), stat="density") + facet_wrap(~variable, scale="free")
    }

    # Histogram of recorded variables
    if (plot == "histogram") {
        p <- ggplot(xM) + geom_histogram(aes(x=value)) + facet_wrap(~variable, scale="free")
    }

    # Latest values in the mcmc iterations
    if (plot == "trace") {
        # get the last 200 records
        xMt <- xM[xM$iter > xd[max(1,nrow(xd)-200),"iter"],]
        if (severalChains) {
            p <- ggplot(xMt) + geom_path(aes(x=iter, y=value, colour=chain)) + facet_wrap(~variable, scale="free")
        } else {
            p <- ggplot(xMt) + geom_path(aes(x=iter, y=value)) + facet_wrap(~variable, scale="free")
        }
    }

    # All values in the mcmc iterations
    if (plot == "history") {
        if (severalChains) {
            p <- ggplot(xM) + geom_path(aes(x=iter, y=value, colour=chain)) + facet_wrap(~variable, scale="free")
        } else {
            p <- ggplot(xM) + geom_path(aes(x=iter, y=value)) + facet_wrap(~variable, scale="free")
        }
    }

    # Autocorrelation function
    if (plot == "autocor") {
        require("plyr")
        # compute ACF per chain and variable
        xACF <- ddply(xM, ~chain+variable, function(y) {
            a <- acf(y$value, plot=FALSE)
            # compute confidence interval
            # (same for all since we use as many points)
            # from plot.acf source :
            ci <- 0.95      # percentage of the confidence interval
            cl <- qnorm((1 + ci)/2)/sqrt(a$n.used)
            data.frame(lag=a$lag, ACF=a$acf, cl)
        })

        if (severalChains) {
            p <- ggplot(xACF) + geom_bar(aes(x=lag, y=ACF, fill=chain), stat="identity", position="dodge") + facet_wrap(~variable) + geom_hline(aes(yintercept=c(cl, -cl)), linetype="dotted")
        } else {
            p <- ggplot(xACF) + geom_bar(aes(x=lag, y=ACF), stat="identity", position="dodge") + geom_hline(aes(yintercept=c(cl, -cl)), linetype="dotted") + facet_wrap(~variable)
        }
    }

    # Evolution of quantiles in time for recorded variables
    if (plot == "quantiles") {
        warning("Not implemented yet")
        p <- NULL
    }

    return(p)
}
