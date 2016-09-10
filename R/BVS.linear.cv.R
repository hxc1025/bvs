BVS.linear.cv <-
function(x, y, v1.all = NULL, v1.length = 100, FUN = BVS.linear.1, fold = 5, L = 100, a0 = 2, b0 = 2, lambda = 1, nu = 1, iter.max = 100, thres = 10^-3, normalize = TRUE, trace = FALSE, param.ini = list(mu=NULL, phi=NULL, sigma_j_sq=NULL, sigma2=NULL, theta=NULL), my.seed = NULL){

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

    # K-fold id. Shuffling the data. 
    if(!is.null(my.seed)){
        set.seed(my.seed)
    }
    sample.id = sample(1:n, n)
    test.id.matrix = matrix(head(sample.id, (n%/%fold)*fold), n%/%fold, fold)
    
    # initialize 
    param.ini.temp = param.ini
    cv.error = matrix(NA, fold, length(v1.all))

    # main loop
    for(k in 1:fold){
        test.id = drop(test.id.matrix[,k])
        x.norm.cv = Norm(x[-test.id,])
        for(i in 1:length(v1.all)){            
            BVS.temp = FUN(x.norm.cv$x, y[-test.id], v1 = v1.all[i], normalize = FALSE, trace = FALSE, a0 = a0, b0 = b0, lambda = lambda, nu = nu, iter.max = iter.max, thres = thres, param.ini = param.ini.temp)
            BVS.temp$x.mean = x.norm.cv$x.mean
            BVS.temp$x.sd = x.norm.cv$x.sd
            BVS.temp$remove.var = x.norm.cv$remove.var
            y.hat = predict.BVS.linear(BVS.temp, x.new = x[test.id,], type = 1)
            cv.error[k,i] = sum((y.hat - y[test.id])^2)
            param.ini.temp = list(
                theta = BVS.temp$theta,
                sigma2 = BVS.temp$sigma2
            )
        }
    }

    # The mean cross-validated error - a vector of length length(v1.all).
    cvm = colMeans(cv.error)
    cvsd = apply(cv.error, 2, sd)/sqrt(n)
    cvup = cvm + cvsd
    cvlo = cvm - cvsd
    
    v1.min = v1.all[tail(which(cvm == min(cvm)), 1)]
    v1.1se = v1.all[tail(which(cvm <= min(cvup)), 1)]
    
    BVS.cv = list(v1 = v1.all, cvm = cvm, cvsd = cvsd, cvup = cvup, cvlo = cvlo, v1.min = v1.min, v1.1se = v1.1se, x.mean = x.mean, x.sd = x.sd, remove.var = remove.var, y.mean = y.mean)
    
    BVS.cv
    
}
