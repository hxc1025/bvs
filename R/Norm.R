Norm <-
function(x){

    x.sd = apply(x, 2, sd)
    remove.var = which(x.sd==0)
    if(length(remove.var)>0){
        # print(cat("Varible ", remove.var, " removed."))    
        x = x[, -remove.var]
        x.sd = x.sd[-remove.var]
    }
    x.mean = colMeans(x)
    x = t((t(x)-x.mean)/x.sd)

    list(x = x, x.mean = x.mean, x.sd = x.sd, remove.var = remove.var)
    
}
