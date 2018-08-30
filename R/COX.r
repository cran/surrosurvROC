COX.ft=function(data0, eval.time, cpt,n.cpt){

  ###
  #1. cox
  ###
  #1.1 sorting
  data0=data0[order(data0$X),]

  x=data0$X
  delta=data0$D
  y=data0$Y
  ipw=data0$IPW

  n=length(x)
  n.ipw=sum(ipw)

  y.unq=sort(unique(y))

  #1.4 est beta using NR
  options(warn=2)
  cox.warning<-try(cox.fit <- coxph(Surv(x,delta)~y, weights=ipw),silent=TRUE)
  options(warn=0)
  if(class(cox.warning)=="try-error"){
    tpf=fpf=auc=y.cpt=y.a=beta=NA
  }else{
    beta<-unname(cox.fit$coefficient)
    bh.cox=basehaz(cox.fit, centered=FALSE)

    bh=bh.cox$hazar # it includes ipw, i.e the same as bh above.
    bt=bh.cox$time  # sorted by time
    eval.pt=tail(which(bt<=eval.time),1)
    if(sum(eval.pt)==0) eval.pt=1

    bh.eval.pt=bh[eval.pt]
    S.hat=function(y, bh, beta) return(exp(-bh*exp(beta*y))) #at eval.time

    ###
    #2. mzr
    ###
    #2.1. tpf, fpf at eval.time
    S.eval.pt=S.hat(y, bh=bh.eval.pt, beta=beta)
    S.t=sum(ipw*S.eval.pt)/n.ipw

    tpf=fpf=ppv=npv=rep(NA,n.cpt)
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
  }

  return(list(tpf=tpf,fpf=fpf,auc=auc,y.a=y.a))
}



