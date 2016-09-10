BVS.logistic.cv.2 <-
function(x, y, m=NULL, v1.list=NULL, fold=5, a0=1, b0=1, iter.max=20, thres=10^-3, normalize=FALSE, param.ini = list(alpha=NULL, mu=NULL, phi=NULL, sigma2=NULL, theta=NULL, w.E = NULL), LogLoss = FALSE, use.mean = TRUE, type=NULL){
    
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
    if(any(is.null(param.ini$alpha))){
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
        sigma2 = 1
    }else{
        sigma2 = param.ini$sigma2
    }
    
    # fit all the BVS models and extract the unique models according to the selected variables. 
    n.v1 = length(v1.list)
    candidate.model.all = matrix(NA, p, n.v1)
    BVS.all = vector("list", n.v1)
    for(i in 1:n.v1){
        BVS.all[[i]] = BVS.logistic(x = x, y = y, m = m, normalize = FALSE, v1 = v1.list[i], iter.max = iter.max, thres = thres, param.ini = param.ini, a0 = a0, b0 = b0)
        candidate.model.all[,i] = (BVS.all[[i]]$phi>0.5)
    }
    candidate.model = unique(candidate.model.all, MARGIN=2)
    if(all(colSums(candidate.model)==0)){
        stop("Use smaller v1 to include variables!")
    }
    candidate.model = candidate.model[,which(colSums(candidate.model)>0), drop = FALSE]
    
    n.candidate = ncol(candidate.model)
    v1.id = rep(NA, n.candidate)
    for(j in 1:n.candidate){
        v1.id[j] = tail(which(colSums(candidate.model.all - candidate.model[,j]) == 0), 1)
    }
    
    # cross-validation
    if(1 %in% type){
        cv.err.1 = matrix(NA, n.candidate, fold)
    }else{
        cv.err.1 = NULL
    }
    
    if(2 %in% type){
        cv.err.2 = matrix(NA, n.candidate, fold)
    }else{
        cv.err.2 = NULL
    }

    if(3 %in% type){
        cv.err.3 = matrix(NA, n.candidate, fold)
    }else{
        cv.err.3 = NULL
    }
    
    ntest = floor(n/fold)
    ntrain = n - ntest
    
    cv.id = sample(1:n, n)
    
    for(k in 1:fold){
        test.id = cv.id[(1:ntest) + (k-1)*ntest]
        train.id = setdiff(1:n, test.id)
        for(i in 1:length(v1.id)){
            if(1 %in% type | 2 %in% type){
                tmpBVS = BVS.logistic(x = x[train.id,], y = y[train.id], m = m, normalize = FALSE, v1 = v1.list[v1.id[i]], iter.max = iter.max, thres = thres, param.ini = param.ini, a0 = a0, b0 = b0)
            }
            if(3 %in% type){
                tmpfit = glm(y~., data=data.frame(y=y[train.id], x=x[train.id,which(candidate.model[,i])]), family = "binomial")
            }
            
            if(1 %in% type){
                py.hat = predict.BVS.logistic(tmpBVS, x[test.id,], type = 1)
                if(LogLoss){
                    cv.err.1[i,k] = log.loss(py.hat$p, y[test.id])
                }else{
                    cv.err.1[i,k] = class.err(py.hat$y, y[test.id])
                }
            }
            if(2 %in% type){
                py.hat = predict.BVS.logistic(tmpBVS, x[test.id,], type = 2)
                if(LogLoss){
                    cv.err.2[i,k] = log.loss(py.hat$p, y[test.id])
                }else{
                    cv.err.2[i,k] = class.err(py.hat$y, y[test.id])
                }
            }
            if(3 %in% type){
                p.hat = predict(tmpfit, data.frame(x=x[test.id,which(candidate.model[,i])]), type = "response")
                y.hat = ifelse(p.hat>0.5, 1, 0)
                if(LogLoss){
                    cv.err.3[i,k] = log.loss(p.hat, y[test.id])
                }else{
                    cv.err.3[i,k] = class.err(y.hat, y[test.id])
                }
            }
        }
    }

    best.model.id.1 = NULL; best.model.id.2 = NULL; best.model.id.3 = NULL
    var.id.1 = NULL; var.id.2 = NULL; var.id.3 = NULL
    v1.best.1 = NULL; v1.best.2 = NULL; v1.best.3 = NULL
    BVS.1 = NULL; BVS.2 = NULL; BVS.3 = NULL    
    
    # find the largest v1 (if there is a tie) that has the lowest cv error.
    if(use.mean){
        if(1 %in% type){
            best.model.id.1 = length(v1.id) - which.min(rev(rowMeans(cv.err.1))) + 1
        }
        if(2 %in% type){
            best.model.id.2 = length(v1.id) - which.min(rev(rowMeans(cv.err.2))) + 1
        }
        if(3 %in% type){
            best.model.id.3 = length(v1.id) - which.min(rev(rowMeans(cv.err.3))) + 1
        }
    }else{
        if(1 %in% type){
            best.model.id.1 = length(v1.id) - which.min(rev(apply(cv.err.1, 1, median))) + 1
        }
        if(2 %in% type){
            best.model.id.2 = length(v1.id) - which.min(rev(apply(cv.err.2, 1, median))) + 1
        }
        if(3 %in% type){
            best.model.id.3 = length(v1.id) - which.min(rev(apply(cv.err.3, 1, median))) + 1
        }    
    }
    
    if(1 %in% type){
        var.id.1 = candidate.model[, best.model.id.1]
        v1.best.1 = v1.list[v1.id[best.model.id.1]]
        BVS.1 = BVS.all[[v1.id[best.model.id.1]]]
        BVS.1$x.mean = x.mean
        BVS.1$x.sd = x.sd
        BVS.1$remove.var = remove.var
    }
    
    if(2 %in% type){
        var.id.2 = candidate.model[, best.model.id.2]
        v1.best.2 = v1.list[v1.id[best.model.id.2]]
        BVS.2 = BVS.all[[v1.id[best.model.id.2]]]
        BVS.2$x.mean = x.mean
        BVS.2$x.sd = x.sd
        BVS.2$remove.var = remove.var
    }
        
    if(3 %in% type){
        var.id.3 = candidate.model[, best.model.id.3]
        v1.best.3 = v1.list[v1.id[best.model.id.3]]
        BVS.3 = BVS.all[[v1.id[best.model.id.3]]]
        BVS.3$x.mean = x.mean
        BVS.3$x.sd = x.sd
        BVS.3$remove.var = remove.var
    }
    
    sparse = list(var.id=var.id.1, v1.best = v1.best.1, BVS = BVS.1, cv.err = cv.err.1)
    
    dense = list(var.id=var.id.2, v1.best = v1.best.2, BVS = BVS.2, cv.err = cv.err.2)
    
    refit = list(var.id=var.id.3, v1.best = v1.best.3, BVS = BVS.3, cv.err = cv.err.3)
    
    list(sparse = sparse, dense = dense, refit = refit)
}
