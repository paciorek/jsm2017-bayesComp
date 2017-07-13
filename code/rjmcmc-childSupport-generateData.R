## unfortunately we realized after writing about the child support enforcement example from Gelman and Hill "Data Analysis Using Regression and Multilevel/Hierarchical Models" that that particular dataset was not publicly available.

## This file simulates data of the same form to illustrate the methodology. In this case we simulate the data such that there is only a small relationship between enforcement and outcome. Also note that for simplicity I leave out the 'mom race' variable.

set.seed(1)
m <- 20
cities <- data.frame(id = letters[1:m], enforce = rnorm(m), benefit = rnorm(m))

n <- 1367
persons <- data.frame(ID = 1:n, dad_age = runif(n, 17, 40),
                      city_ID = sample(1:m, n, replace = TRUE))

cities_true_randomEffect <- rnorm(m, 0, .2)
benefit_coef <- .5
enforcement_coef <- .2
city_effect <- enforcement_coef * cities$enforce + benefit_coef *
    cities$benefit + cities_true_randomEffect

age_coef <- 0.01
expit <- function(x) exp(x)/(1+exp(x))
true_prob <- expit(city_effect[persons$city_ID] + age_coef*persons$dad_age)

persons$support <- rbinom(n, 1, true_prob)
