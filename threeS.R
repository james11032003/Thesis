threeS <- function(dataY, dataX, criterion1, criterion2, const1, const2, const3, const4){
  n=length(dataY)

  oga_result = OGA(dataY, dataX, const1)
  hdic_result = HDIC(dataY, dataX, criterion1, oga_result, const2)
  trim_result = Trim(dataY, dataX, criterion1, oga_result, hdic_result, const2)

  ffit1 = lm(dataY~dataX[,trim_result$relevent_variable])
  rrrooo=ffit1$residuals
  ppooss_res=which((rrrooo^2)<0.00000001)
  if(length(ppooss_res)>0){
    rrrooo[ppooss_res]=sign(rrrooo[ppooss_res])*0.0001
  }
  new_yy = log(as.vector(rrrooo)^2)
  
  oga_result1 = OGA(new_yy, dataX, const3)
  hdic_result1 = HDIC(new_yy, dataX, criterion2, oga_result1, const4)
  trim_result1 = Trim(new_yy, dataX, criterion2, oga_result1, hdic_result1, const4)

  return(list("location_var"=trim_result$relevent_variable,"dispersion_var"=trim_result1$relevent_variable))
}