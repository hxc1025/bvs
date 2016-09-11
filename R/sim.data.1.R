sim.data.1 <-
function(n = 40, sigma2 = 3, beta = c(3, 1.5, 0, 0, 2, 0, 0, 0), rho = 0.5, Sigma = NULL){
    p = length(beta)
    if(is.null(Sigma)){
        Sigma = outer(1:p, 1:p, function(u,v){rho^(abs(u-v))})
    }
    
    x = rMVN(n, mu = rep(0,p), Sigma = Sigma)
    y = rnorm(n, x%*%beta, sigma2)
    
    list(x = x, y = y)
}
