BVS.logistic <-
function(x, y, m=NULL, v1=1, a0=1, b0=1, iter.max=10^3, thres=10^-3, normalize=TRUE, trace=FALSE, param.ini=list(alpha=NULL, mu=NULL, phi=NULL, sigma2=NULL, theta=NULL, w.E = NULL), keep.data = FALSE, phi.c = 4){

    x = as.matrix(x)
    y = as.vector(y)
    n = length(y)
    
    if(is.null(m)){
        m = rep(1, n)
    }
        
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
    
    p = ncol(x)
    
    # initial value
    if(is.null(param.ini$alpha)){
        alpha = 0
    }else{
        alpha = param.ini$alpha
    }
    if(any(is.null(param.ini$w.E))){
        w.E = rep(0.25, n)
    }else{
        w.E = param.ini$w.E
    }
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
    if(any(is.null(param.ini$phi))){
        phi = rep(theta, p)
    }else{
        phi = param.ini$phi
    }
    if(any(is.null(param.ini$sigma2))){
        sigma2 = rep(v1, p)
    }else{
        sigma2 = param.ini$sigma2
    }
    
    # initialize 
    b.m = phi*mu
    b.D = mu^2*phi*(1-phi)+sigma2*phi
    Ent = rep(Inf, p)
    phi.flag = rep(TRUE, p)    
    phi.logit = rep(NA, p)
    
    if(trace){
        post.save = list(
            alpha = rep(NA, iter.max), 
            mu = matrix(NA, p, iter.max), 
            sigma2 = matrix(NA, p, iter.max), 
            w.E = matrix(NA, n, iter.max), 
            phi = matrix(NA, p, iter.max), 
            theta = rep(NA, n, iter.max)
        )
    }else{
        post.save = NULL
    }
    
    res = Inf

    t = 0
    
    while(t<iter.max & res > thres){
        t = t + 1
        
        # w
        w.c = sqrt((alpha + as.vector(x%*%b.m))^2 + colSums(t(x^2)*b.D))
        w.E = m/(2*w.c)*tanh(w.c/2)
        
        # Working residual
        kappa = y - m/2 - (alpha + x%*%b.m)*w.E
        
        # beta and gamma, need to iterative over j=1:p
        for(j in 1:p){
            kappa = kappa + x[,j]*phi[j]*mu[j]*w.E
            sigma2[j] = 1/(sum(w.E*x[,j]^2)+1/v1)
            mu[j] = sum(kappa*x[,j]) * sigma2[j]
            if(phi.flag[j]){
                phi.logit[j] = log(theta/(1-theta)) + log(sigma2[j]/v1)/2 + mu[j]^2/sigma2[j]/2
                phi[j] = 1/(1+exp(-phi.logit[j]))
                phi.flag[j] = (abs(phi.logit[j])<phi.c)
            }
            kappa = kappa - x[,j]*phi[j]*mu[j]*w.E
        }

        b.m = phi*mu
        b.D = mu^2*phi*(1-phi)+phi*sigma2        
        
        # alpha
        kappa = kappa + alpha*w.E
        alpha = sum(kappa)/sum(w.E)
        kappa = kappa - alpha*w.E
        
        # theta
        theta = (sum(phi)+a0-1)/(p+a0+b0-2)

        # Stopping Criterion
        Ent.old = Ent
        Ent = -log(phi^phi)-log((1-phi)^(1-phi))
        res = max(abs(Ent - Ent.old))
        
        # save the result
        if(trace){
            post.save$alpha[t] = alpha
            post.save$mu[,t] = mu
            post.save$sigma2[,t] = sigma2
            post.save$phi[,t] = phi
            post.save$w.E[,t] = w.E        
            post.save$theta[t] = theta    
        }
    }

    if(trace){
        post.save$alpha = post.save$alpha[1:t]
        post.save$mu = post.save$mu[,1:t]
        post.save$sigma2 = post.save$sigma2[,1:t]
        post.save$phi = post.save$phi[,1:t]
        post.save$w.E = post.save$w.E[,1:t]
        post.save$theta= post.save$theta[1:t]
    }
    
    if(!keep.data){
        x = NULL
        y = NULL
    }
    
    list(mu = drop(mu), sigma2 = drop(sigma2), phi = drop(phi), w.E = drop(w.E), alpha = drop(alpha), theta = drop(theta), post.save = post.save, x.mean = x.mean, x.sd = x.sd, remove.var = remove.var, x = x, y = y)
    
}
