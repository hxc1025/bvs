BVS.linear.Gibbs <-
function(x, y, v1 = 1, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 10^5, thres = 10^-3, normalize = FALSE, burn.in = 1000){

    x = as.matrix(x)
    y = as.vector(y)

    if(normalize){
        x.norm = Norm(x)
        n = x.norm$n
        p = x.norm$p
        x = x.norm$x
        x.mean = x.norm$x.mean
        x.sd = x.sd
        remove.var = x.norm$remove.var
    }else{
        n = nrow(x)
        p = ncol(x)
        x.mean = NULL
        x.sd = NULL
        remove.var = NULL
    }

    y.mean = mean(y)
    y = y - y.mean
    
    # initial value
    sigma2 = 1/rgamma(1, nu/2, (nu*lambda)/2)
    theta = rbeta(1, a0, b0)
    gamma = rbinom(p, 1, theta)
    beta = rnorm(p, 0, sqrt(v1*sigma2))*gamma
    res = Inf
    Ent = ent(gamma)
    
    # auxiliary
    kappa = drop(y - x%*%beta)
    xx = colSums(x^2)
    mu = rep(0, p)
    sigma_j_sq = rep(0, p)
    phi = rep(0, p)

    post.save = list(
        # alpha = rep(NA, iter.max), 
        beta = matrix(NA, p, iter.max), 
        gamma = matrix(NA, p, iter.max),
        sigma2 = rep(NA, iter.max), 
        theta = rep(NA, iter.max)
    )

    t = 0
    while(t<iter.max && res>thres){
        t = t + 1
        
        # beta and gamma
        for(j in 1:p){
            kappa = kappa + x[,j]*beta[j]
            mu[j] = (x[,j] %*% kappa) / (xx[j] + 1/v1)
            sigma_j_sq[j] = sigma2/(xx[j] + 1/v1)
            BF = log(theta/(1-theta) / sqrt(1+v1*xx[j])) + mu[j]^2/sigma_j_sq[j]/2
            phi[j] = 1/(1+exp(-BF))
            gamma[j] = rbinom(1, 1, phi[j])
            if(gamma[j]==1){
                beta[j] = rnorm(1, mu[j], sqrt(sigma_j_sq[j]))
            }else{
                beta[j] = 0
            }
            kappa = kappa - x[,j]*beta[j]
        }
        
        # sigma2
        sigma2 = 1/rgamma(1, (n+nu+sum(gamma))/2, (sum(kappa^2)+sum(beta[which(gamma!=0)]^2)/v1+nu*lambda)/2) 
                
        # theta
        theta = rbeta(1, sum(gamma)+a0, p-sum(gamma)+b0)    
                
        # Save the results
        post.save$beta[,t] = beta
        post.save$gamma[,t] = gamma
        post.save$sigma2[t] = sigma2 
        post.save$theta[t] = theta
        
        # Entropy
        if((t %%  1000 == 0) && (t>burn.in)){
            Ent.old = Ent
            Ent = ent(rowMeans(post.save$gamma[, (t-burn.in):t]))
            res = abs(Ent-Ent.old)
        }
    }    
    
    post.save$beta = post.save$beta[,1:t]
    post.save$gamma = post.save$gamma[,1:t]
    post.save$sigma2 = post.save$sigma2[1:t] 
    post.save$theta = post.save$theta[1:t]
    
    list(post.save = post.save, x.mean = x.mean, x.sd = x.sd, remove.var = remove.var, y.mean = y.mean)
    
}
