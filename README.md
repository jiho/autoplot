# autoplot

[ggplot2](http://ggplot2.org/) is a versatile plotting package which allows to produce almost any kind of plot from data stored as a data.frame, by combining unit elements. However, it requires the user to design the plot entirely, from scratch. Many R functions for statistical analyses (linear models, factorial analyses, etc.) output objects of a given class and allow to easily plot classic diagnostics using plot(), by defining a specialized method for this generic function. This package aims at reproducing this functionality in ggplot2 while benefiting from its improved versatility. It provides two sets of methods : (i) the fortify() methods extract data from the original object and format it as a data.frame; (ii) the autoplot() methods use these data.frames and leverage ggplot2 to produce diagnostic plots.

## Installation

The package is under development and is not on CRAN yet. To install it, the simplest way it therefore

    # install.packages("devtools")
    library("devtools")
    install_github("autoplot", "jiho")

## Development

All suggestions and improvements are welcome. Suggestions are preferred under the form of Github "issues" (see the tab above). Improvements should be submitted under the form of "pull-requests" (again, see the tab above).

The best way to create a pull-request (according to me) is:

1. fork the repository on GitHub's website. It creates your very own copy of `autoplot`

2. clone *my version* the project

        git clone https://github.com/jiho/autoplot.git

3. create a branch

        git branch myfeature
        git checkout myfeature

4. make the changes you want and commit them to the branch

        <code, code>
        git commit
    
5. push the branch to *your* github repository

        git remote add mine https://github.com/<yourname>/autoplot.git
        git push myfeature mine

6. create the pull request from your branch to my master branch, on GitHub's website.

This way, you can easily update the master branch to get my latest changes and base your next contribution on it (which is necessary for things not to break):

    # update the master branch to match mine
    git checkout master
    git pull
    
    # create a new branch
    git branch myotherfeature
    git checkout myotherfeature
    
    # create a new change
    <code, code>
    git commit
    
    # push it online
    git push myotherfeature mine

Technically, the documentation is written with [roxygen](http://roxygen.org/), unit tests are performed with [testthat](https://github.com/hadley/testthat).
