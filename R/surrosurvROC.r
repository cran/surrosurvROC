surrosurvROC=function(DATA, method, pred.time, wt=NULL, span=NULL, b.rep=200){

  Marker=DATA$Marker
  Time=DATA$Time
  Status=DATA$Status
  STime=DATA$STime
  SStatus=DATA$SStatus

  method=toupper(method)
  V=as.numeric(!is.na(Time) & !is.na(Status))
  n=length(V)
  if(is.null(wt)) wt=rep(1,n)
  wt0=rep(1,length(V))
  alldata=data.frame(V=V,Y=Marker,X2=Time,D2=Status,X1=STime,D1=SStatus,IPW0=wt0,IPW1=wt)

  RES=updtdROC(alldata=alldata, method=method, pred.time=pred.time, span=span, b.rep=b.rep)
  RES3=RES$res3 #update

  plot(RES3$tpf~RES3$fpf, type='l', ylim=c(0,1), xlim=c(0,1), ylab='TP', xlab='FP', main=paste("ROC at t=",pred.time,sep=''))
  points(RES3$tpf~RES3$fpf, type='p', pch=19, cex=0.5)

  points(RES3$tpf.bup~RES3$fpf, type='l')
  points(RES3$tpf.blw~RES3$fpf, type='l')

  return(list(TP=RES3$tpf, FP=RES3$fpf, AUC=RES3$auc,
              pred.time=pred.time, cut.values=RES3$y.a))
}
