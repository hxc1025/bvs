ent <-
function(p){
    if(any(p>1 | p<0)){stop("Invalid probabilities.")}
    p = as.vector(p)
    -min(log(p^p) + log((1-p)^(1-p)))
}
