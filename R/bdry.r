bdry.ft=function(x){
  x[x<0]=0; x[x>1]=1
  return(x)
}

NA.ft=function(x){
  x2=sum(x); n2=length(x)
  if(is.na(x2)+is.infinite(x2)+is.nan(x2))
    x=rep(NA,n2)
  return(x)
}
