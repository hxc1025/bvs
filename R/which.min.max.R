which.min.max <-
function(x){
    length(x) - which.min(rev(x)) + 1
}
