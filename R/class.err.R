class.err <-
function(y.hat, y.test){
    mean(y.hat != y.test)
}
