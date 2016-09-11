predict.BVS.logistic <-
function(BVS, x.new, m.new = NULL, type = 1, phi.thres = 0.5){

    x.new = as.matrix(x.new)
    if(length(BVS$remove.var)!=0){
        x.new = x.new[,-BVS$remove.var]
    }
    if((!is.null(BVS$x.mean)) && (!is.null(BVS$x.sd))){
        x.new = t((t(x.new) - BVS$x.mean)/BVS$x.sd)
    }    
    
    if(type == 1){
        beta.hat = ifelse(BVS$phi>phi.thres, BVS$mu, 0)
        p.hat = x.new %*% beta.hat + BVS$alpha
    }else if(type == 2){
        beta.hat = with(BVS, phi * mu)
        p.hat = x.new %*% beta.hat + BVS$alpha
    # }else if(type == 3){
        # if(any(c(is.null(x), is.null(y)))){
            # if(any(c(is.null(BVS$x), is.null(BVS$y)))){
                # stop("Please provide training data for refitting.")
            # }else{
                # x = BVS$x
                # y = BVS$y
            # }
        # }
        # if(length(BVS$remove.var)!=0){
            # x = x[,-BVS$remove.var]
        # }
        # var.id = which(BVS$phi > phi.thres)
        # BVS.refit = glm(y~., data = data.frame(x = x[,var.id], y = as.factor(y)), family="binomial")
        # p.hat = predict(BVS.refit, newdata = data.frame(x = x.new[,var.id]), type = c("response"))
    # }else if(type == 4){
        # p.hat = predict(BVS$refit, newdata = data.frame(x=x.new[,BVS$var.id]), type = "response")
    }else{
        stop("Incorrect type of prediction!")
    }
    
    # if(! (type %in% c(3,4))){
        p.hat = 1/(1+exp(-p.hat))
    # }
    
    y.hat = ifelse(p.hat>0.5, 1, 0)
    
    if(!is.null(m.new)){
        y.hat = round(m.new*p.hat)
    }else{
        y.hat = ifelse(p.hat>0.5, 1, 0)
    }
    
    list(y = as.vector(y.hat), p = as.vector(p.hat))
    
}
