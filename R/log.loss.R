log.loss <-
function(p.hat, y.test, cut=0.01){
  for(i in 1:length(p.hat)){
    if(p.hat[i]<cut){
      p.hat[i] = cut
    }else if(p.hat[i]>1-cut){
      p.hat[i] = 1-cut
    }
  }
    -mean(ifelse(y.test, log(p.hat), log(1-p.hat)))
}
