beta.convert.back <-
function(BVS, phi.cut = 0.5){

    beta.bar.temp = with(BVS, ifelse(phi>phi.cut, mu, 0))
    beta.bar.temp = with(BVS, beta.bar.temp/x.sd)
    p2 = length(BVS$remove.var)
    if(p2!=0){
        beta.bar = rep(NA, BVS$p + p2)
        beta.bar[-BVS$remove.var] = beta.bar.temp
    }else{
        beta.bar = beta.bar.temp
    }
    
    beta.bar
    
}
