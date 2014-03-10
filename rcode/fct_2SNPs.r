# 1.11 La fonction qui resume pour deux SNPs et une seule exposition

   fct_2SNPs<-function(dat,outc,vexp,gm1,gc1,gm2,gc2,vvaraj,seu=NULL){
                                            # dat la base de donnÃ©e deja jumÃ©llÃ© au snp
                                            # snp le snp
                                            # la variable d'exposition
                                            # on suppose que la vecteur
                                            # vvaraj : le vecteur des facteur d'ajustement
                                            # qrt la variable qui permet de connaitre si nous sommes dans les seuille en cartile pour binaire
                                            # seu : seuil pour la decouper la variable
                                            # dat0<-dat[dat[,outc]==0,];vexp0<-dat0[vexp];vexp1<-dat[vexp]
                                            # le cas binaire oÃ¹ on precise le seuille de de la variable de dichotomisation
                                               nvar<-cathe_bi2(dat,outc,vexp,seu);nvar<-t(t(nvar));nch<-paste(vexp,".","dch",sep="");colnames(nvar)<-nch
                                                                   dat1<-data.frame(dat,nvar)
                                                                   fl0=formula(paste(outc,"~",paste(c(paste(c(nch,gm1,gc1,gm2,gc2),collapse="+"),paste(nch,":",gm1,sep=""),paste(nch,":",gc1,sep=""),paste(nch,":",gm2,sep=""),paste(nch,":",gc2,sep=""),paste(vvaraj,collapse="+")),collapse="+"),sep=""))
                                                                   mod0<-glm(fl0,data=dat1,family=binomial)
                                                                   cof<-coef(mod0)
                                                                   Ic1<-confint(mod0)
                                                                   # or et intervalle de confiance pour le model
                                                                   OR<-round(c(1,exp(cof[2])),1)
                                                                   Icor<-rbind(c(0,0),round(exp(Ic1)[2,],1))
                                                                   ty_w<-cbind(OR,Icor)

                                                                   # la covariance
                                                                   mat1<-vcov(mod0);dg<-diag(mat1);

                                                                   n<-length(cof)
                                                                   
                                                                   # or pour mÃ¨re gène 1
                                                                   orjm1<-c(1,exp(cof[2]+cof[n-3]));orjme1<-round(orjm1,1)
                                                                   Intm1_inf=c(0,round(orjm1[2]*exp(-1.96*sqrt(dg[2]+dg[n-3]+2*mat1[n-3,2])),1))
                                                                   Intm1_sup=c(0,round(orjm1[2]*exp(1.96*sqrt(dg[2]+dg[n-3]+2*mat1[n-3,2])),1))
                                                                   ty_me1<-cbind(orjme1,Intm1_inf,Intm1_sup)

                                                                   # or pour enfant gène 1
                                                                   orje1<-c(1,exp(cof[2]+cof[n-2]));orjen1<-round(orje1,1)
                                                                   Inte1_inf=c(0,round(orje1[2]*exp(-1.96*sqrt(dg[2]+dg[n-2]+2*mat1[n-2,2])),1))
                                                                   Inte1_sup=c(0,round(orje1[2]*exp(1.96*sqrt(dg[2]+dg[n-2]+2*mat1[n-2,2])),1))
                                                                   ty_en1<-cbind(orjen1,Inte1_inf,Inte1_sup)
                                                                   
                                                                   # or pour mÃ¨re gène 2
                                                                   orjm2<-c(1,exp(cof[2]+cof[n-1]));orjme2<-round(orjm2,1)
                                                                   Intm2_inf=c(0,round(orjm2[2]*exp(-1.96*sqrt(dg[2]+dg[n-1]+2*mat1[n-1,2])),1))
                                                                   Intm2_sup=c(0,round(orjm2[2]*exp(1.96*sqrt(dg[2]+dg[n-1]+2*mat1[n-1,2])),1))
                                                                   ty_me2<-cbind(orjme2,Intm2_inf,Intm2_sup)

                                                                   # or pour enfant gène 2
                                                                   orje2<-c(1,exp(cof[2]+cof[n]));orjen2<-round(orje2,1)
                                                                   Inte2_inf=c(0,round(orje2[2]*exp(-1.96*sqrt(dg[2]+dg[n]+2*mat1[n,2])),1))
                                                                   Inte2_sup=c(0,round(orje2[2]*exp(1.96*sqrt(dg[2]+dg[n]+2*mat1[n,2])),1))
                                                                   ty_en2<-cbind(orjen2,Inte2_inf,Inte2_sup)

                                                                   # or des interaction genotype 1 mere expo
                                                                   ORme1<-round(c(1,exp(cof[n-3])),1)
                                                                   IORme1<-rbind(c(0,0),round(exp(Ic1)[n-3,],1))
                                                                   ty_ORme1<-cbind(ORme1,IORme1)

                                                                   # or des interaction genotype 1 enfant expo
                                                                   ORee1<-round(c(1,exp(cof[n-2])),1)
                                                                   IORee1<-rbind(c(0,0),round(exp(Ic1)[n-2,],1))
                                                                   ty_ORee1<-cbind(ORee1,IORee1)

                                                                   # or des interaction genotype 2 mere expo
                                                                   ORme2<-round(c(1,exp(cof[n-1])),1)
                                                                   IORme2<-rbind(c(0,0),round(exp(Ic1)[n-1,],1))
                                                                   ty_ORme2<-cbind(ORme2,IORme2)

                                                                   # or des interaction genotype 2 enfant expo
                                                                   ORee2<-round(c(1,exp(cof[n])),1)
                                                                   IORee2<-rbind(c(0,0),round(exp(Ic1)[n,],1))
                                                                   ty_ORee2<-cbind(ORee2,IORee2)

                                                                   # recap
                                                                   vgm1<-cbind(ty_w,ty_me1,ty_en1,ty_ORme1,ty_ORee1);
                                                                   vgm2<-cbind(ty_w,ty_me2,ty_en2,ty_ORme2,ty_ORee2);
                                                                   vgm<-rbind(vgm1,vgm2)
                                                                   colnames(vgm)<-NULL;rownames(vgm)<-NULL
                                                                   nnc<-c("OR1","L1","U1","OR2","L2","U2","OR3","L3","U3","OR4","L4","U4","OR5","L5","U5")
                                                                   nnr<-c(paste(substr(gm1,5,nchar(gm1)),"_","Q1",sep=""),paste(substr(gm1,5,nchar(gm1)),"_","Q2",sep=""),
                                                                     paste(substr(gm2,5,nchar(gm2)),"_","Q1",sep=""),paste(substr(gm2,5,nchar(gm2)),"_","Q2",sep=""))
                                                                   colnames(vgm)<-nnc
                                                                   rownames(vgm)<-nnr

                                           return(list(matR=vgm,model1=mod0))
                                           }
