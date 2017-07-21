OGA <- function(dataY, dataX, const){
  n = length(dataY)
  pn = ncol(dataX)
  y = dataY-mean(dataY)
  x = sweep(dataX, 2, apply(dataX, 2, mean))
  kn = floor(const*((n/log(pn))^(1/2)))
  x_rest_pos = 1:pn
  oga_path = NULL
  mse = NULL

  u = y
  tx = x
  for(i in 1:min(kn,pn)){
    if(i<pn){
      cor_all = apply(tx, 2, function(v) abs(cor(u,v)))
      cor_max_pos = which.max(cor_all)[1]
      oga_path = c(oga_path,x_rest_pos[cor_max_pos])
      x_rest_pos = x_rest_pos[-cor_max_pos]
      oga_x = as.data.frame(x[,oga_path])
      fit = lm(u~.-1, data=oga_x)   
      u = as.vector(fit$residuals)
      mse[i] = mean(u^2)
      tx = tx[,-cor_max_pos]
    }else if(i==pn){
      oga_path[i] = x_rest_pos[1]
      fit = lm(u~x-1) 
      u = as.vector(fit$residuals)
      mse[i] = mean(u^2)
    }
  }
  return(list("oga_path" = oga_path, "mse" = mse))
}