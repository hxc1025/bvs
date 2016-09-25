#bvs
R package (Bayesian Variable Selection using Spike-and-Slab Prior for Linear and Generalized Linear Regressions)

##Maintainer
Xichen Huang, the University of Illinois at Urbana-Champaign

##Description
Variational Bayesian algorithms for fitting the Bayesian linear and logistic regression with the spike-and-slab prior (Normal and Point Mass).

##Install the package
```{r}
# Need R package "devtools"
# install.packages("devtools")
# Install "bvs"
library("devtools")
install_github("hxc1025/bvs")
library('bvs')
```

##Examples
```{r}
library('bvs')
```

###Bayesian Linear Regression
####Generate simulation data of Fan and Li (2001) example 1.
```{r}
# The default of sim.data.1 is (n = 40, sigma2 = 3, beta = c(3, 1.5, 0, 0, 2, 0, 0, 0), rho = 0.5)
set.seed(100)
train.data = sim.data.1()
x = train.data$x
y = train.data$y
```

####Fit Bayesian linear regression.
```{r}
# Algorithm 1 (A slight modification of Carbonetto and Stephens 2012)
bvs.fit.1 = BVS.linear.1(x, y, v1=1)
# Returns a list containing mu, phi, sigma_j_2, etc. 
round(with(bvs.fit.1, cbind(mu, sigma_j_sq, phi)), 2)

# Algorithm 2 (Batch-wise updating)
bvs.fit.2 = BVS.linear.2(x, y, v1=1)
round(with(bvs.fit.2, cbind(mu, sigma_j_sq, phi)), 2)
```

####Cross-validation For v1
```{r}
# Given a sequence of v1, perform cross-validation to select v1.
# Can be used for algorithm 1
bvs.cv.1 = BVS.linear.cv(x, y, v1.all=10^seq(-2,2,0.1), FUN=BVS.linear.1, my.seed=200)
bvs.cv.1$v1.min
# Or algorithm 2
bvs.cv.2 = BVS.linear.cv(x, y, v1.all=10^seq(-2,2,0.1), FUN=BVS.linear.2, my.seed=200)
bvs.cv.2$v1.min
```

####Prediction
```{r}
# Generate test data
set.seed(300)
test.data = sim.data.1()
x.test = test.data$x
y.test = test.data$y
# Use the above bvs.fit.2 
y.hat = predict.BVS.linear(bvs.fit.2, x.test, type=1)
# Type can be 1, 2, 3. For more details, use ?predict.BVS.linear
```

###Bayesian logistic regression
####Generate data
```{r}
n = 100; p = 20; 
beta.0 = c(3,2,1)
beta.0 = c(beta.0, rep(0, p-length(beta.0)))
set.seed(100)
x = matrix(rnorm(n*p), n, p)
y = x%*%beta.0 + rnorm(n, sd=1)
y = ifelse(y>0, 1, 0)
# Model fit
bvs.fit = BVS.logistic(x, y, v1=1)
round(head(cbind(beta.0, with(bvs.fit, cbind(mu, sigma2, phi))), 8), 2)
```

####Cross-validation for v1
```{r}
bvs.cv = BVS.logistic.cv(x, y, v1.list=10^seq(-2,2,0.1), my.seed=400)
bvs.cv$v1.best
# Also return the best model. To access, use bvs.cv$BVS
```

####Prediction on a new data set
```{r}
# Test data
set.seed(500)
x.test = matrix(rnorm(n*p), n, p)
y.test = x.test%*%beta.0 + rnorm(n, sd=1)
y.test = ifelse(y.test>0, 1, 0)
# Prediction on bvs.fit
y.hat.1 = predict.BVS.logistic(bvs.fit, x.test, type=1)$y
table(y.hat.1, y.test)
# Prediction given bvs.cv
y.hat.2 = predict.BVS.logistic(bvs.cv$BVS, x.test, type=1)$y
table(y.hat.2, y.test)
# Increase accuracy by 1%. 
```
