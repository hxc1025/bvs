BVS.linear.2 <-
function(x, y, v1=1, a0=2, b0=2, lambda=1, nu=1, iter.max=100, thres=10^-3, normalize=TRUE, phi.cut=0.01, trace=FALSE, param.ini=list(mu=NULL, phi=NULL, sigma_j_sq=NULL, sigma2=NULL, theta=NULL), use.lambda_n1=FALSE, keep.data=FALSE){

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
        phi = rep(1, p)
    }else{
        phi = param.ini$phi
    }
    if(any(is.null(param.ini$sigma_j_sq))){
        sigma_j_sq = rep(v1*sigma2, p)
    }else{
        sigma_j_sq = param.ini$sigma_j_sq
    }

    # auxiliary parameters and variables
    xx = t(x) %*% x
    xy = drop(t(x) %*% y)
    update.id = 1:p
    B = xx; diag(B) = 0; 
    
    q = p
    
    if(n<=p || use.lambda_n1==TRUE){
        x.svd = svd(x)        
        if(!require("MASS")){
            library("MASS")
        }
        S.sq = diag(length(x.svd$d))
        diag(S.sq) = (x.svd$d)^2
        S.sq.inv = ginv(S.sq)
        s = sum(diag(S.sq.inv) * x.svd$d^2 != 0)
        F = x.svd$v[,1:s]
        S.sq.inv = diag(s)
        diag(S.sq.inv) = (x.svd$d[1:s])^(-2)
        if(use.lambda_n1==TRUE){
            lambda_n1 = x.svd$d[s]^2
        }
    }
    
    t = 0
    res = Inf
    Ent = Inf
    
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
                
        # mu        
        if(n>q){
            if(t == 1){
                A.inv = solve(t(xx*phi) + diag(diag(xx)*(1-phi)+1/v1))
                mu = A.inv %*% t(x) %*% y
                mu = drop(mu)
            }else{
                C.inv = diag(length(update.id))
                diag(C.inv) = 1/(phi[update.id] - phi.old[update.id])
                U = B[,update.id]
                V = matrix(0, q, p); V[,update.id] = diag(q)
                A.inv = A.inv - A.inv%*%U%*%solve(C.inv+V%*%A.inv%*%U, V%*%A.inv)
                mu = A.inv %*% t(x) %*% y
                mu = drop(mu)
            }
        }else{
            Delta.inv = diag(1/(diag(xx)*(1-phi)+1/v1))
            V = F*phi
            A.inv =  Delta.inv - Delta.inv%*%F%*%ginv(S.sq.inv+t(V)%*%Delta.inv%*%F)%*%t(V)%*%Delta.inv
            mu = A.inv %*% t(x) %*% y
            mu = drop(mu)
        }
        
        # sigma_j_sq
        if(use.lambda_n1==TRUE){
            sigma_j_sq = rep(sigma2/(lambda_n1 + 1/v1), p)
        }else{
            sigma_j_sq = sigma2/(diag(xx) + 1/v1)
        }
           
        # phi
        phi.old = phi
        if(t == 1){
            phi.logit = log(theta/(1-theta)) + log(sigma_j_sq/v1/sigma2)/2+mu^2*(diag(xx)+1/v1)/sigma2/2
            phi = 1/(1+exp(-phi.logit))
            phi.temp = rbind(phi, 1-phi)
            update.id = which(apply(phi.temp, 2, min) > phi.cut)
        }else{
            phi.temp = rbind(phi, 1-phi)
            update.id = which(apply(phi.temp, 2, min) > phi.cut)
            phi.logit[update.id] = log(theta/(1-theta)) + log(sigma_j_sq[update.id]/v1/sigma2)/2 + mu[update.id]^2*(diag(xx)[update.id]+1/v1)/sigma2/2
            phi[update.id] = 1/(1+exp(-phi.logit[update.id]))
        }
        
        # theta
        theta = (sum(phi)+a0-1)/(p+a0+b0-2)

        # sigma2
        kappa = y - x%*%(mu*phi)
        E.kappa.2 = sum(kappa^2) + sum(diag(xx)*(phi*(1-phi)*mu^2 + phi*sigma_j_sq))
        sigma2 = (E.kappa.2 + sum(phi*(mu^2+sigma_j_sq))/v1 + nu*lambda)/(n + sum(phi) + nu + 2)
        
        # Entropy
        Ent.old = Ent
        Ent = log(phi^phi) + log((1-phi)^(1-phi))
        
        # Stopping Criterion
        res = max(abs(Ent - Ent.old))
        
        q = length(update.id)
        
        # Save the results
        if(trace){
            post.save$mu[,t] = mu
            post.save$sigma_j_sq[,t] = sigma_j_sq
            post.save$phi[,t] = phi
            post.save$sigma2[t] = sigma2 
            post.save$theta[t] = theta
            post.save$Ent[,t] = Ent
        }
        
        if(q==0){
            break
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
    
    list(mu = unname(drop(mu)), sigma_j_sq = unname(drop(sigma_j_sq)), phi = unname(drop(phi)), sigma2 = sigma2, theta = theta, post.save = post.save, x.mean = x.mean, x.sd = x.sd, remove.var = remove.var,y.mean = y.mean, x = x, y = y)
}
