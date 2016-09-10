BVS.linear.1 <-
function(x, y, v1 = 1, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 100, thres = 10^-3, normalize = TRUE, trace = FALSE, param.ini = list(mu=NULL, phi=NULL, sigma_j_sq=NULL, sigma2=NULL, theta=NULL), keep.data = FALSE){

    x = as.matrix(x)
    y = as.vector(y)
    
    if(normalize){
        x.norm = Norm(x)
        x = x.norm$x
        x.mean = x.norm$x.mean
        x.sd = x.norm$x.sd
        remove.var = x.norm$remove.var
    }else{
        x.mean = NULL
        x.sd = apply(x,2,sd)
        remove.var = which(x.sd==0)
        if(length(remove.var)==0){
            remove.var = NULL
        }else{
            x = x[,-remove.var]
        }        
        x.sd = NULL
    }

    n = nrow(x)
    p = ncol(x)
    
    y.mean = mean(y)
    y = y - y.mean

    # initial value
    if(any(is.null(param.ini$mu))){
        mu = rep(0, p)
    }else{
        mu = param.ini$mu
    }
    if(is.null(param.ini$theta)){
        theta = a0/(a0+b0)
    }else{
        theta = param.ini$theta
    }
    if(is.null(param.ini$sigma2)){
        sigma2 = var(y)
    }else{
        sigma2 = param.ini$sigma2
    }
    if(any(is.null(param.ini$phi))){
        phi = rep(theta, p)
    }else{
        phi = param.ini$phi
    }
    if(any(is.null(param.ini$sigma_j_sq))){
        sigma_j_sq = rep(v1*sigma2, p)
    }else{
        sigma_j_sq = param.ini$sigma_j_sq
    }
    
    # auxiliary parameters and variables
    kappa = y - x%*%(mu*phi)
    xx = colSums(x^2)
    
    t = 0
    res = Inf
    Ent = rep(Inf, p)

    if(trace){
        post.save = list(
            mu = matrix(NA, p, iter.max), 
            sigma_j_sq = matrix(NA, p, iter.max), 
            phi = matrix(NA, p, iter.max),
            sigma2 = rep(NA, iter.max), 
            theta = rep(NA, iter.max),
            Ent = matrix(NA, p, iter.max)
        )
    }else{
        post.save = NULL
    }
    
    while(t<iter.max && res>thres){
        t = t + 1

        # mu, sigma_j_sq, and phi
        for(j in 1:p){
            kappa = kappa + x[,j]*(mu[j]*phi[j])
            sigma_j_sq[j] = sigma2/(xx[j] + 1/v1)
            mu[j] = x[,j] %*% kappa / (xx[j] + 1/v1)
            log.odds = log(theta/(1-theta)) + log(sqrt(sigma_j_sq[j]/sigma2/v1)) + mu[j]^2/sigma_j_sq[j]/2
            phi[j] = 1/(1+exp(-log.odds))
            kappa = kappa - x[,j]*(mu[j]*phi[j])
        }

        # theta
        theta = (sum(phi)+a0-1)/(p+a0+b0-2)

        # sigma2
        E.kappa.2 = sum(kappa^2) + sum(xx*(phi*(1-phi)*mu^2 + phi*sigma_j_sq))
        sigma2 = (E.kappa.2 + sum(phi*(mu^2+sigma_j_sq))/v1 + nu*lambda)/(n + sum(phi) + nu + 2)

        # Entropy
        Ent.old = Ent
        Ent = log(phi^phi) + log((1-phi)^(1-phi))
        
        # Stopping Criterion
        res = max(abs(Ent - Ent.old))
        
        # Save the results
        if(trace){
            post.save$mu[,t] = mu
            post.save$sigma_j_sq[,t] = sigma_j_sq
            post.save$phi[,t] = phi
            post.save$sigma2[t] = sigma2 
            post.save$theta[t] = theta
            post.save$Ent[,t] = Ent
        }
    }
    
    if(trace){
        post.save$mu = post.save$mu[,1:t]
        post.save$sigma_j_sq = post.save$sigma_j_sq[,1:t]
        post.save$phi = post.save$phi[,1:t]
        post.save$sigma2 = post.save$sigma2[1:t] 
        post.save$theta= post.save$theta[1:t]
        post.save$Ent = post.save$Ent[,1:t]
    }
    
    if(!keep.data){
        x = NULL
        y = NULL
    }
    
    list(mu = mu, sigma_j_sq = sigma_j_sq, phi = phi, sigma2 = sigma2, theta = theta, post.save = post.save, x.mean = x.mean, x.sd = x.sd, remove.var = remove.var, y.mean = y.mean, x = x, y = y)
}
