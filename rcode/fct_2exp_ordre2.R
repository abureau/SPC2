# 1.11 La fonction qui resume pour un SNPs et 2 expositions
   
   fct_Rsum_2exp_ordre2<-function(dat,outc,vexp,vexp2,gm,gc,vvaraj,seu=NULL,seu2=NULL){
                                            # dat la base de donnée deja juméllé au snp
                                            # snp le snp 
                                            # vexp: la variable d'exposition
                                            # vstrat: la variable de stratification
                                            # vvaraj : le vecteur des facteur d'ajustement
                                            # qrt la variable qui permet de connaitre si nous sommes dans les seuille en cartile pour binaire 
                                            # seu : seuil pour la decouper la variable 
                                            # dat0<-dat[dat[,outc]==0,];vexp0<-dat0[vexp];vexp1<-dat[vexp]
                                            # le cas binaire où on precise le seuille de de la variable de dichotomisation
                                               nvar<-cbind(cathe_bi2(dat,outc,vexp,seu),cathe_bi2(dat,outc,vexp2,seu2));nch<-paste(vexp,".","dch",sep="");nch2<-paste(vexp2,".","dch",sep="");colnames(nvar)<-c(nch,nch2) 
                                                                   dat1<-data.frame(dat,nvar)
                                                                   fl0=formula(paste(outc,"~",paste(c(paste(c(nch,nch2,gm,gc),collapse="+"),paste(nch,":",gm,sep=""),paste(nch,":",gc,sep=""),paste(nch2,":",gm,sep=""),paste(nch2,":",gc,sep=""),paste(vvaraj,collapse="+")),collapse="+"),sep=""))   
                                                                   mod0<-glm(fl0,data=dat1,family=binomial)
                                                                   cof<-coef(mod0)
                                                                   Ic1<-confint(mod0)
                                                                   # or et intervalle de confiance pour l'effet de l'exposition 1
                                                                   OR<-round(c(1,exp(cof[2])),1)
                                                                   Icor<-rbind(c(0,0),round(exp(Ic1)[2,],1))
                                                                   ty_w<-cbind(OR,Icor)
                                                                   
                                                                   # la covariance
                                                                   mat1<-vcov(mod0);dg<-diag(mat1);
                                                                   
                                                                   # or pour mère 
                                                                   n<-length(cof)
                                                                   orjm<-c(1,exp(cof[2]+cof[n-3]));orjme<-round(orjm,1)
                                                                   Intm_inf=c(0,round(orjm[2]*exp(-1.96*sqrt(dg[2]+dg[n-3]+2*mat1[n-3,2])),1))
                                                                   Intm_sup=c(0,round(orjm[2]*exp(1.96*sqrt(dg[2]+dg[n-3]+2*mat1[n-3,2])),1))
                                                                   ty_me<-cbind(orjme,Intm_inf,Intm_sup)
                                                                   
                                                                   # enfant
                                                                   orje<-c(1,exp(cof[2]+cof[n-2]));orjen<-round(orje,1)
                                                                   Inte_inf=c(0,round(orje[2]*exp(-1.96*sqrt(dg[2]+dg[n-2]+2*mat1[n-2,2])),1))
                                                                   Inte_sup=c(0,round(orje[2]*exp(1.96*sqrt(dg[2]+dg[n-2]+2*mat1[n-2,2])),1))
                                                                   ty_en<-cbind(orjen,Inte_inf,Inte_sup)
                                                                   
                                                                   # or des interaction genotype mere expo
                                                                   ORme<-round(c(1,exp(cof[n-3])),1)
                                                                   IORme<-rbind(c(0,0),round(exp(Ic1)[n-3,],1))
                                                                   ty_ORme<-cbind(ORme,IORme)
                                                                   
                                                                   # or des interaction genotype enfant expo
                                                                   ORee<-round(c(1,exp(cof[n-2])),1)
                                                                   IORee<-rbind(c(0,0),round(exp(Ic1)[n-2,],1))
                                                                   ty_ORee<-cbind(ORee,IORee)
                                                                   
                                                                   # or et intervalle de confiance pour l'effet de l'exposition 2
                                                                   OR2<-round(c(1,exp(cof[3])),1)
                                                                   Icor2<-rbind(c(0,0),round(exp(Ic1)[3,],1))
                                                                   ty_w2<-cbind(OR,Icor)
                                                                                                                                   
                                                                   # or pour mère 
                                                                   orjm2<-exp(cof[3]+cof[n-1]);orjme2<-c(1,round(orjm2,1))
                                                                   Intm2_inf=c(0,round(orjm2*exp(-1.96*sqrt(dg[3]+dg[n-1]+2*mat1[n-1,3])),1))
                                                                   Intm2_sup=c(0,round(orjm2*exp(1.96*sqrt(dg[3]+dg[n-1]+2*mat1[n-1,3])),1))
                                                                   ty_me2<-cbind(orjme2,Intm2_inf,Intm2_sup)
                                                                   
                                                                   # enfant
                                                                   orje2<-exp(cof[3]+cof[n]);orjen2<-c(1,round(orje2,1))
                                                                   Inte2_inf=c(0,round(orje2*exp(-1.96*sqrt(dg[3]+dg[n]+2*mat1[n,3])),1))
                                                                   Inte2_sup=c(0,round(orje2*exp(1.96*sqrt(dg[3]+dg[n]+2*mat1[n,3])),1))
                                                                   ty_en2<-cbind(orjen2,Inte2_inf,Inte2_sup)
                                                                   
                                                                   # or des interaction genotype mere expo
                                                                   ORme2<-round(c(1,exp(cof[n-1])),1)
                                                                   IORme2<-rbind(c(0,0),round(exp(Ic1)[n-1,],1))
                                                                   ty_ORme2<-cbind(ORme2,IORme2)
                                                                   
                                                                   # or des interaction genotype enfant expo
                                                                   ORee2<-round(c(1,exp(cof[n])),1)
                                                                   IORee2<-rbind(c(0,0),round(exp(Ic1)[n,],1))
                                                                   ty_ORee2<-cbind(ORee2,IORee2)
                                                                   
                                                                   
                                                                   # la repartition 
                                                                    date1<-dat1[is.na(dat1[,vexp])!=TRUE,];nbr.na=dim(dat1[is.na(dat1[,vexp])==TRUE,])[1]
                                                                    dat10<-date1[date1$outc==0,];dat11<-date1[date1$outc==1,]
                                                                   #type wild
                                                                   tw0<-c(dim(dat11[dat11[,dim(dat11)[2]]==0,])[1],dim(dat10[dat10[,dim(dat10)[2]]==0,])[1])
                                                                   tw1<-c(dim(dat11[dat11[,dim(dat11)[2]]==1,])[1],dim(dat10[dat10[,dim(dat10)[2]]==1,])[1])
                                                                   tw<-rbind(tw0,tw1);
                                                                   
                                                                   tagm1<-dat11[dat11[,gm]>=1,];tagm0<-dat10[dat10[,gm]>=1,]
                                                                   tage1<-dat11[dat11[,gc]>=1,];tage0<-dat10[dat10[,gc]>=1,]
                                                                  
                                                                   twm0<-c(dim(tagm1[tagm1[,dim(tagm1)[2]]==0,])[1],dim(tagm0[tagm0[,dim(tagm0)[2]]==0,])[1])
                                                                   twm1<-c(dim(tagm1[tagm1[,dim(tagm1)[2]]==1,])[1],dim(tagm0[tagm0[,dim(tagm0)[2]]==1,])[1])
                                                                   twe0<-c(dim(tage1[tage1[,dim(tage1)[2]]==0,])[1],dim(tage0[tage0[,dim(tage0)[2]]==0,])[1])
                                                                   twe1<-c(dim(tage1[tage1[,dim(tage1)[2]]==1,])[1],dim(tage0[tage0[,dim(tage0)[2]]==1,])[1])
                                                                   twm<-rbind(twm0,twm1);twe<-rbind(twe0,twe1)
                                                                   # recap
                                                                   vgm1<-cbind(tw,ty_w,twm,ty_me,twe,ty_en,ty_ORme,ty_ORee);
                                                                   vgm2<-cbind(tw,ty_w,twm,ty_me2,twe,ty_en2,ty_ORme2,ty_ORee2);
                                                                   vgm = rbind(vgm1,vgm2)
                                                                   colnames(vgm)<-NULL;rownames(vgm)<-NULL
                                                                   nnc<-c("Cas1","Cont1","OR1","L1","U1","Cas2","Cont2","OR2","L2","U2","Cas3","Cont3","OR3","L3","U3","OR4","L4","U4","OR5","L5","U5")
                                                                   nnr<-c(paste(substr(gm,5,nchar(gm)),"_","Q1",nch,sep=""),paste(substr(gm,5,nchar(gm)),"_","Q2",nch,sep=""),paste(substr(gm,5,nchar(gm)),"_","Q1",nch2,sep=""),paste(substr(gm,5,nchar(gm)),"_","Q2",nch2,sep=""))
                                                                   colnames(vgm)<-nnc
                                                                   rownames(vgm)<-nnr
                                                                   
                                           return(list(matR=vgm,model1=mod0))
                                           }
