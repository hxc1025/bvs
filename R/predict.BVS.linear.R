predict.BVS.linear <-
function(BVS, x.new, type = 1, phi.thres = 0.5){

    x.new = as.matrix(x.new)
    if(length(BVS$remove.var)>0){
        x.new = x.new[,-BVS$remove.var]
    }
    if((!is.null(BVS$x.mean)) && (!is.null(BVS$x.sd))){
        x.new = t((t(x.new) - BVS$x.mean)/BVS$x.sd)
    }    
    
    if(type == 1){
        beta.hat = ifelse(BVS$phi>phi.thres, BVS$mu, 0)
        y.hat = x.new %*% beta.hat
    }else if(type == 2){
        beta.hat = with(BVS, phi * mu)
        y.hat = x.new %*% beta.hat
    }else if(type == 3){
        used.var = which(BVS$phi > phi.thres)
        BVS.refit = lm(y~.+0, data = data.frame(x = BVS$x[,used.var], y = BVS$y))
        y.hat = predict(BVS.refit, newdata = data.frame(x = x.new[,used.var]), type = c("response"))
    }else{
        stop("Incorrect type of prediction!")
    }
    
    y.hat = y.hat + BVS$y.mean
    
    y.hat
    
}
