HDIC <- function(dataY, dataX, criterion, OGA_result, const){
  n = length(dataY)
  pn = ncol(dataX)

  oga_path = OGA_result$oga_path
  mse = OGA_result$mse
  kn = length(mse)
  
  if(criterion==1){
    omega = const*log(n)
  }else if(criterion==2){
    omega = const*log(log(n))
  }

  hdic = (n*log(mse))+((1:kn)*omega*(log(pn)))
  hat_kn = which.min(hdic)[1]

  return(list("hat_kn" = hat_kn, "hdic_value" = hdic[hat_kn]))
}