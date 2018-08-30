updtdROC=function(alldata, method, pred.time, span, b.rep, bootstrap='yes'){

  ###
  #1. initial
  ###
  n=nrow(alldata)
  cpt=sort(unique(alldata$Y))
  n.cpt=length(cpt)
  eval.time=pred.time

  res0=list();res0$tpf=res0$fpf=res0$auc=res0$y.a=NA;
  res1=list();res1$tpf=res1$fpf=res1$auc=res1$y.a=NA;
  res2=list();res2$tpf=res2$fpf=res2$auc=res2$y.a=NA;
  res3=list();res3$tpf=res3$fpf=res3$auc=res3$y.a=NA;

  ###
  #2. compute st, sp, auc using GK
  ###
  # data.P1A, phase 1, All
  # data.P1V, phase 1, V
  # data.P2V, phase 2, V
  data.P1A=data.frame(X=alldata$X1,D=alldata$D1,Y=alldata$Y,IPW=alldata$IPW0)
  idx.V=which(alldata$V==1);
  data.P1V=(data.frame(X=alldata$X1,D=alldata$D1,Y=alldata$Y,IPW=alldata$IPW1))[idx.V,]
  data.P2V=(data.frame(X=alldata$X2,D=alldata$D2,Y=alldata$Y,IPW=alldata$IPW1))[idx.V,]

  if(method=="KNN"){
    res0=KNN.ft(data.P1A, eval.time,cpt,n.cpt, span)
    res1=KNN.ft(data.P1V, eval.time,cpt,n.cpt, span)
    res2=KNN.ft(data.P2V, eval.time,cpt,n.cpt, span)
  }else if(method=="COX"){
    res0=COX.ft(data.P1A, eval.time,cpt,n.cpt)
    res1=COX.ft(data.P1V, eval.time,cpt,n.cpt)
    res2=COX.ft(data.P2V, eval.time,cpt,n.cpt)
  }

  ###
  #3.bootstrap & update
  ###
  n.cpt0=length(res0$y.a)
  n.cpt1=length(res1$y.a) #n.cpt1=n.cpt2

  res3$tpf=res3$fpf=rep(NA,n.cpt1)
  res3$y.a=res2$y.a

  if(bootstrap=='yes'){

    bootstrap='no'

    #a. single mzr
    b.tpf0=b.fpf0=matrix(NA,b.rep,n.cpt0) #ROC using all y
    b.tpf1=b.fpf1=matrix(NA,b.rep,n.cpt1)
    b.tpf2=b.fpf2=matrix(NA,b.rep,n.cpt1)
    b.tpf3=b.fpf3=matrix(NA,b.rep,n.cpt1) #Update for the Y, observed on P2V only

    #b. global mzr
    b.auc0=b.auc1=b.auc2=b.auc3=NA

    for(rr in 1:b.rep){
      samp.id=sample(1:n, n, replace = TRUE)
      alldata.B=alldata[samp.id,]
      b.res=updtdROC(alldata.B, method, pred.time, span, b.rep, bootstrap)
      b.res$res2

      #a. single mzr
      idx0=which(res0$y.a %in% b.res$res0$y.a); b.tpf0[rr,idx0]=b.res$res0$tpf; b.fpf0[rr,idx0]=b.res$res0$fpf;
      idx1=which(res1$y.a %in% b.res$res1$y.a); b.tpf1[rr,idx1]=b.res$res1$tpf; b.fpf1[rr,idx1]=b.res$res1$fpf;
      idx2=which(res2$y.a %in% b.res$res2$y.a); b.tpf2[rr,idx2]=b.res$res2$tpf; b.fpf2[rr,idx2]=b.res$res2$fpf;

      #b. global mzr
      b.auc0[rr]=b.res$res0$auc;     b.auc1[rr]=b.res$res1$auc;     b.auc2[rr]=b.res$res2$auc;
    }

    #var
    #a. single mzr
    v.tpf0=apply(b.tpf0,2,var,na.rm=TRUE); v.tpf1=apply(b.tpf1,2,var,na.rm=TRUE); v.tpf2=apply(b.tpf2,2,var,na.rm=TRUE);
    v.fpf0=apply(b.fpf0,2,var,na.rm=TRUE); v.fpf1=apply(b.fpf1,2,var,na.rm=TRUE); v.fpf2=apply(b.fpf2,2,var,na.rm=TRUE);
    res0$v.tpf=v.tpf0; res1$v.tpf=v.tpf1; res2$v.tpf=v.tpf2;
    res0$v.fpf=v.fpf0; res1$v.fpf=v.fpf1; res2$v.fpf=v.fpf2;

    #b. global mzr
    v.auc0=var(b.auc0,na.rm=TRUE); v.auc1=var(b.auc1,na.rm=TRUE); v.auc2=var(b.auc2,na.rm=TRUE);
    res0$v.auc=v.auc0;     res1$v.auc=v.auc1;     res2$v.auc=v.auc2;

    #3.com for res3 (updated)
    #3.a. single mz
    idx01=which(res0$y.a%in%res1$y.a)
    b.tpf01=b.tpf0[,idx01]-b.tpf1
    b.fpf01=b.fpf0[,idx01]-b.fpf1

    v.tpf01=apply(b.tpf01,2,var,na.rm=TRUE)
    v.fpf01=apply(b.fpf01,2,var,na.rm=TRUE)

    for(j in 1:n.cpt1){
      na.cov.tpf=sum(!is.na(b.tpf01[,j]+b.tpf2[,j]))
      na.cov.fpf=sum(!is.na(b.fpf01[,j]+b.fpf2[,j]))

      if(na.cov.tpf>1){
        cov.tpf=cov(b.tpf01[,j],b.tpf2[,j],use='complete')
        res3$tpf[j]=res2$tpf[j]-cov.tpf/v.tpf01[j]*((res0$tpf[idx01])[j]-res1$tpf[j])
        b.tpf3[,j] = b.tpf2[,j]-cov.tpf/v.tpf01[j]*((b.tpf0[,idx01])[,j]-  b.tpf1[,j])
      }
      if(na.cov.fpf>1){
        cov.fpf=cov(b.fpf01[,j],b.fpf2[,j],use='complete')
        res3$fpf[j]=res2$fpf[j]-cov.fpf/v.fpf01[j]*((res0$fpf[idx01])[j]-res1$fpf[j])
        b.fpf3[,j] = b.fpf2[,j]-cov.fpf/v.fpf01[j]*((b.fpf0[,idx01])[,j]-  b.fpf1[,j])
      }
    }

    v.tpf3=apply(b.tpf3,2,var,na.rm=TRUE);
    v.fpf3=apply(b.fpf3,2,var,na.rm=TRUE);

    res3$v.tpf=v.tpf3;
    res3$v.fpf=v.fpf3;

    #3.b. global mzr
    b.auc01=b.auc0-b.auc1
    v.auc01=var(b.auc01,na.rm=TRUE)

    na.cov.auc=sum(!is.na(b.auc01+b.auc2))

    if(na.cov.auc>1){
      cov.auc=cov(b.auc01,b.auc2,use='complete')
      res3$auc=res2$auc-cov.auc/v.auc1*(res0$auc-res1$auc)
      b.auc3=b.auc2-cov.auc/v.auc1*(b.auc0-b.auc1)
    }

    v.auc3=var(b.auc3,na.rm=TRUE);
    res3$v.auc=v.auc3;

    #4.confidence interval b/o quantile & bootstrap
    res3$tpf.blw=apply(b.tpf3,2,quantile, probs=0.025, na.rm=TRUE); res3$tpf.bup=apply(b.tpf3,2,quantile, probs=0.975, na.rm=TRUE);
    res3$fpf.blw=apply(b.fpf3,2,quantile, probs=0.025, na.rm=TRUE); res3$fpf.bup=apply(b.fpf3,2,quantile, probs=0.975, na.rm=TRUE);
    res3$auc.blw=quantile(b.auc3,2,probs=0.025); res3$auc.bup=quantile(b.auc3,2,probs=0.975);

    #5.bdry
    res3$tpf.blw=bdry.ft(res3$tpf.blw); res3$tpf.bup=bdry.ft(res3$tpf.bup)
    res3$fpf.blw=bdry.ft(res3$fpf.blw); res3$fpf.bup=bdry.ft(res3$fpf.bup)
    res3$auc.blw=bdry.ft(res3$auc.blw); res3$auc.bup=bdry.ft(res3$auc.bup)

    #6. sorting
    roc3.data=data.frame(tpf=res3$tpf,tpf.bup=res3$tpf.bup,tpf.blw=res3$tpf.blw, fpf=res3$fpf,fpf.bup=res3$fpf.bup,fpf.blw=res3$fpf.blw)
    roc3.rdata=roc3.data[order(roc3.data$fpf),]
    res3$tpf=roc3.rdata$tpf; res3$tpf.bup=roc3.rdata$tpf.bup; res3$tpf.blw=roc3.rdata$tpf.blw;
    res3$fpf=roc3.rdata$fpf; res3$fpf.bup=roc3.rdata$fpf.bup; res3$fpf.blw=roc3.rdata$fpf.blw;

    #7. NA pt are removed
    idx3=!is.na(res3$tpf);
    idx3=!is.na(res3$tpf.bup);
    idx3=!is.na(res3$tpf.blw);

    res3$tpf=res3$tpf[idx3]
    res3$fpf=res3$fpf[idx3]
    res3$tpf.bup=res3$tpf.bup[idx3]
    res3$tpf.blw=res3$tpf.blw[idx3]
    res3$y.a=res3$y.a[idx3]

    #8. pava
    res3$tpf=isoreg(res3$tpf)$yf
    res3$tpf.bup=isoreg(res3$tpf.bup)$yf
    res3$tpf.blw=isoreg(res3$tpf.blw)$yf

    #9 . dj end point
    if(min(res3$tpf)!=0 || min(res3$fpf)!=0){
      res3$tpf=c(0,res3$tpf)
      res3$fpf=c(0,res3$fpf)
      res3$tpf.bup=c(0,res3$tpf.bup)
      res3$tpf.blw=c(0,res3$tpf.blw)
      res3$y.a=c(-Inf,res3$y.a)
    }

    if(max(res3$tpf)!=1 || max(res3$fpf)!=1){
      res3$tpf=c(res3$tpf,1)
      res3$fpf=c(res3$fpf,1)
      res3$tpf.bup=c(res3$tpf.bup,1)
      res3$tpf.blw=c(res3$tpf.blw,1)
      res3$y.a=c(res3$y.a,Inf)
    }
  }

  return(list(res0=res0, res1=res1, res2=res2, res3=res3))
}
