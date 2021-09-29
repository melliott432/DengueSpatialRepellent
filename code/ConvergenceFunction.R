convergence_function = function(samples1, samples2, column="alpha"){
  
  min_rows = min(nrow(samples1), nrow(samples2))
  
  samples1 = head(samples1, min_rows)
  samples2 = head(samples2, min_rows)
  
  theta = cbind(as.matrix(samples1[,column]), as.matrix(samples2[,column]))
  res = GelmanRubin(theta)
  return(res$R)
  
}

convergence_all = function(samples1, samples2){
  min_rows = min(nrow(samples1), nrow(samples2))
  
  samples1 = head(samples1, min_rows)
  samples2 = head(samples2, min_rows)
  
  conv = rep(NA, ncol(samples1))
  for(i in 1:ncol(samples1)){
    
    theta = cbind(as.matrix(samples1[,i]), as.matrix(samples2[,i]))
    res = GelmanRubin(theta)
    conv[i] = res$R
  }
  
  return(mean(conv < 1.2) == 0)
}