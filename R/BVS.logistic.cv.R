BVS.logistic.cv <-
function(x, y, m=NULL, v1.list=NULL, fold=10, a0=1, b0=1, iter.max=10^3, thres=10^-3, normalize=TRUE, param.ini = list(alpha=NULL, mu=NULL, phi=NULL, sigma2=NULL, theta=NULL, w.E = NULL), my.seed = NULL, LogLoss = FALSE, type = 1){
    
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
  
    ntest = floor(n/fold)
    ntrain = n - ntest 
  
    if(!is.null(my.seed)){
        set.seed(my.seed)
    }
  
    cv.id = sample(1:n, n)
    cv.err = matrix(NA, length(v1.list), fold)
  
    for(k in 1:fold){
        test.id = cv.id[(1:ntest) + (k-1)*ntest]
        for(i in 1:length(v1.list)){
            bvs_temp = BVS.logistic(x=x[-test.id,], y=y[-test.id], v1=v1.list[i], normalize=FALSE)
            py.hat = predict.BVS.logistic(bvs_temp, x[test.id,], type = type)
            if(LogLoss){
                cv.err[i,k] = log.loss(py.hat$p, y[test.id])
            }else{
                cv.err[i,k] = class.err(py.hat$y, y[test.id])
            }
        }
    }
  
    v1.best = v1.list[which.min(rowMeans(cv.err))]
  
    BVS = BVS.logistic(x, y, v1 = v1.best, normalize = FALSE)
    BVS$x.mean = x.mean
    BVS$x.sd = x.sd
    BVS$remove.var = remove.var
  
    list(v1.best = v1.best, cv.err = cv.err, BVS = BVS)
  
}
