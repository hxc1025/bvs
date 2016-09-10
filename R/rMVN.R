rMVN <-
function(n = 1, mu, Sigma, L = NULL, U = NULL){

    p = length(mu)

    if(is.null(L) && is.null(U)){
        L = t(chol(Sigma))
    }else if(!is.null(L)){
        L = L
    }else if(!is.null(U)){
        L = t(U)
    }

    x = mu + L %*% matrix(rnorm(n*p), p, n)
    
    t(x)
    
}
