KNN.ft=function(data0, eval.time, cpt, n.cpt, span){

  ###
  #1. KNN
  ###
  #1.1 sorting
  data0=data0[order(data0$X),]

  x=data0$X
  delta=data0$D
  y=data0$Y
  ipw=data0$IPW

  n=length(x)
  x.unq=unique(x)
  n.x=length(x.unq)
  n.ipw=sum(ipw)

  #1.2 comp couting and at-risk processes
  dN=Y=matrix(0,n,n.x)
  for(i in 1:n){
    j=which(x.unq==x[i])
    dN[i,j]=1*delta[i]
    Y[i,1:j]=rep(1,j)
  }

  #1.3 compute F.y
  F.y=rep(NA,n)
  for(i in 1:n)
    F.y[i]=sum(ipw*(y<=y[i]))/n.ipw

  y.unq=sort(unique(y))

  #1.4 compute span & KNN surv at each y_i, up to eval time, i.e. S.eval.pt: S_y(t=eval.pt)
  eval.pt=tail(which(x.unq<=eval.time),1)
  if(sum(eval.pt)==0) eval.pt=1

  #S=KNNsurv(n, eval.pt, F.y, ipw, span, Y, dN) #S is computed using Rcpp
  S=matrix(1,n,eval.pt)
  for(i in 1:n){
    wt=rep(0,n)
    F.bar=F.y[i]-F.y

    wt.idx=which( (-span<=F.bar) & (F.bar<=span) )
    wt[wt.idx]=1

    for(k in 1:eval.pt){
      if(k==1){
        S.prev=1
      }else{
        S.prev=S[i,k-1]
      }

      num=sum(wt*ipw*dN[,k])
      den=sum(wt*ipw*Y[,k])

      if(den==0){ #0/0=0
        S[i,k]=S.prev
      }else{
        S[i,k]=S.prev*(1-num/den)
      }
    }
  }

  ###
  #2. mzr
  #2.1. tpf, fpf, ppv, npv at eval.time
  S.eval.pt=S[,eval.pt]
  S.t=sum(ipw*S.eval.pt)/n.ipw

  tpf=fpf=rep(NA,n.cpt)
  for(j in 1:n.cpt){
    S.ct=sum(S.eval.pt*ipw*(y>cpt[j]))/n.ipw
    F.c=sum(ipw*(y<=cpt[j]))/n.ipw
    tpf[j]=(1-F.c-S.ct)/(1-S.t)
    fpf[j]=S.ct/S.t
  }

  #2.2 auc at eval.time
  #(1) tpf&fpf at all y.unq
  auc=NA
  n.unq=length(y.unq)
  tpf.a=fpf.a=rep(NA,n.unq)
  for(i in 1:n.unq){
    S.ct=sum(S.eval.pt*ipw*(y>y.unq[i]))/n.ipw
    F.c=sum(ipw*(y<=y.unq[i]))/n.ipw
    tpf.a[i]=(1-F.c-S.ct)/(1-S.t)
    fpf.a[i]=S.ct/S.t
  }

  #(2) sorted by y
  data0.a=data.frame(tpf.a,fpf.a,y.unq)
  data0.a=data0.a[order(data0.a$y.unq),]
  y.a=data0.a$y.unq
  tpf.a=data0.a$tpf.a
  fpf.a=data0.a$fpf.a

  #(3) including the end points
  y.a=c(-Inf,y.a,Inf)
  tpf.a=c(1,tpf.a,0)
  fpf.a=c(1,fpf.a,0)
  n.a=n.unq+2

  #(4) mid pt integration, like the same approach as SurvivalROC package, Heagerty et al.
  d.fpf=fpf.a[-n.a]-fpf.a[-1]
  mid.tpf=(tpf.a[-n.a] + tpf.a[-1])/2
  auc=sum(d.fpf*mid.tpf)

  #2.4. NA, NaN, Infinite due to S.t=1 or 0; F.c=1 or 0
  tpf=NA.ft(tpf.a)
  fpf=NA.ft(fpf.a)
  auc=NA.ft(auc)

  return(list(tpf=tpf,fpf=fpf,auc=auc,y.a=y.a))
}
