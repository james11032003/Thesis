threeS_loc <- function(dataY, dataX, criterion, const1, const2){
  n=length(dataY)

  oga_result=OGA(dataY,dataX,const1)
  hdic_result=HDIC(dataY,dataX,criterion,oga_result,const2)
  trim_result=Trim(dataY,dataX,criterion,oga_result,hdic_result,const2)

  return(list("oga_path"=oga_result$oga_path,"oga_mse"=oga_result$mse,"hdic_hat_kn"=hdic_result$hat_kn,"hdic_value"=hdic_result$hdic_value,"trim"=trim_result$relevent_variable))
}