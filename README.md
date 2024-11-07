# Sharp Variance Bounds

This package provides a consistent estimate of the sharp variance bound for general regression-adjusted estimators of the average treatment effect in two-arm randomized controlled experiments. This variance estimate can be used to construct the asymptotically narrowest conservative Wald-type confidence interval with the nominal coverage for general regression adjusted estimators for the treatment effect. For more details, see [(Mikhaeil and Green, 2024)](https://arxiv.org/abs/2411.00191).

## Installation 
The package can be installed with:

`
#install.packages("remotes")
remotes::install_github("JonasMikhaeil/SharpVarianceBounds")
`
## Usage

`
library(sharpvar)
population <-create_dataset(1000,0.5,1)
Y1 <- population$Y1
Y0 <- population$Y0
x <- population$X

dat <- randomize(Y0,Y1,x,0.5)
mod <- lm(Y_obs ~ treatment*I(X-mean(X)) , data=dat)

sharp_var <-sharpvar(mod$residuals,dat$treatment,upper=TRUE)
`

See vignette for examples.  
