library(SPmlficmcm)

 # 1.9 la fonction qui transforme une variable continue en binaire (bis)
    cathe_bi2<-function(dat,outc=NULL,var1,seu=NULL,quant=0.75){
                                    # outcome : la variable d'interêt 
                                    dat0<-dat[dat[,outc]==0,];var0<-dat0[,var1]
                                    if(is.null(seu)){
                                                           var<-dat[,var1]
                                                           seu=quantile(var0,probs=quant,na.rm =TRUE)
                                                           qrt1<-ifelse(var>=seu,1,0)+ifelse(is.na(var)==TRUE,NA,0)
                                                           }else{
                                                           var<-dat[,var1]
                                                           qrt1<-ifelse(var>=seu,1,0)+ifelse(is.na(var)==TRUE,NA,0)
                                                           }
                                    return(qrt1)
                                    }    
# 1.11 La fonction qui resume pour un SNPs une seule exposition
   
   fct_Rsum<-function(dat,outc,vexp,gm,gc,vvaraj,seu=NULL){
                                            # dat la base de donn?e deja jum?ll? au snp
                                            # snp le snp 
                                            # la variable d'exposition 
                                            # on suppose que la vecteur 
                                            # vvaraj : le vecteur des facteur d'ajustement
                                            # qrt la variable qui permet de connaitre si nous sommes dans les seuille en cartile pour binaire 
                                            # seu : seuil pour la decouper la variable 
                                            # dat0<-dat[dat[,outc]==0,];vexp0<-dat0[vexp];vexp1<-dat[vexp]
                                            # le cas binaire o? on precise le seuille de de la variable de dichotomisation
                                               nvar<-cathe_bi2(dat,outc,vexp,seu);nvar<-t(t(nvar));nch<-paste(vexp,".","dch",sep="");colnames(nvar)<-nch 
                                                                   dat1<-data.frame(dat,nvar)                         
if (missing(vvaraj)) fl0=formula(paste(outc,"~",paste(c(paste(c(nch,gm,gc),collapse="+"),paste(nch,":",gm,sep=""),paste(nch,":",gc,sep="")),collapse="+"),sep=""))   
else                                                                  fl0=formula(paste(outc,"~",paste(c(paste(c(nch,gm,gc),collapse="+"),paste(nch,":",gm,sep=""),paste(nch,":",gc,sep=""),paste(vvaraj,collapse="+")),collapse="+"),sep=""))   
                                                                   mod0<-glm(fl0,data=dat1,family=binomial)
                                                                   smod<-summary(mod0)
                                                                   cof<-coef(mod0)
                                                                   Ic1<-confint(mod0)
                                                                   # or et intervalle de confiance pour le model
                                                                   OR<-round(c(1,exp(cof[2])),1)
                                                                   Icor<-rbind(c(0,0),round(exp(Ic1)[2,],1))
                                                                   ty_w<-cbind(OR,Icor)
                                                                   
                                                                   # la covariance
                                                                   mat1<-vcov(mod0);dg<-diag(mat1);
                                                                   
                                                                   # or pour m?re 
                                                                   n<-length(cof)
                                                                   orjm<-c(1,exp(cof[2]+cof[n-1]));orjme<-round(orjm,1)
                                                                   Intm_inf=c(0,round(orjm[2]*exp(-1.96*sqrt(dg[2]+dg[n-1]+2*mat1[n-1,2])),1))
                                                                   Intm_sup=c(0,round(orjm[2]*exp(1.96*sqrt(dg[2]+dg[n-1]+2*mat1[n-1,2])),1))
                                                                   ty_me<-cbind(orjme,Intm_inf,Intm_sup)
                                                                   
                                                                   # enfant
                                                                   orje<-c(1,exp(cof[2]+cof[n]));orjen<-round(orje,1)
                                                                   Inte_inf=c(0,round(orje[2]*exp(-1.96*sqrt(dg[2]+dg[n]+2*mat1[n,2])),1))
                                                                   Inte_sup=c(0,round(orje[2]*exp(1.96*sqrt(dg[2]+dg[n]+2*mat1[n,2])),1))
                                                                   ty_en<-cbind(orjen,Inte_inf,Inte_sup)
                                                                   
                                                                   # or de l'interaction genotype mere expo
                                                                   ORme<-round(c(1,exp(cof[n-1])),1)
                                                                   IORme<-rbind(c(0,0),round(exp(Ic1)[n-1,],1))
                                                                   ty_ORme<-cbind(ORme,IORme)
                                                                   
                                                                   # p de l'interaction genotype mere expo
                                                                   pinter_me = c(1,smod$coefficients[n-1,4])
                                                                                                                                      
                                                                   # or des interaction genotype enfant expo
                                                                   ORee<-round(c(1,exp(cof[n])),1)
                                                                   IORee<-rbind(c(0,0),round(exp(Ic1)[n,],1))
                                                                   ty_ORee<-cbind(ORee,IORee)
                                                                   
                                                                   # p de l'interaction genotype mere expo
                                                                   pinter_ee = c(1,smod$coefficients[n,4])
                                                                   
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
                                                                   vgm1<-cbind(tw,ty_w,twm,ty_me,pinter_me,twe,ty_en,pinter_ee,ty_ORme,ty_ORee);
                                                                   colnames(vgm1)<-NULL;rownames(vgm1)<-NULL
                                                                   #nnc<-c("Cas1","Cont1","OR1","Ic1","Ic2","Cas2","Cont2","OR2","Ic1","Ic2","Cas3","Cont3","OR3","Ic1","Ic1","OR4","Ic1","Ic2","OR5","Ic1","Ic2")
                                                                   nnc<-c("Cas1","Cont1","OR1","L1","U1","Cas2","Cont2","OR2","L2","U2","pI2v1","Cas3","Cont3","OR3","L3","U3","pI3v1","OR4","L4","U4","OR5","L5","U5")
                                                                   nnr<-c(paste(substr(gm,5,nchar(gm)),"_","Q1",sep=""),paste(substr(gm,5,nchar(gm)),"_","Q2",sep=""))
                                                                   colnames(vgm1)<-nnc
                                                                   rownames(vgm1)<-nnr
                                                                   
                                           return(list(matR=vgm1,model1=mod0))
                                           }
                                    
                                                                 
# 1.11 La fonction qui resume pour un SNPs une seule exposition avec m?thode de Chen

   fct_Rsum_CMCM<-function(dat,N,outc,vexp,gm,gc,vvaraj=NULL,gma=gm,gca=gc,seu=NULL,quant=0.75,typ=2,alpha=0.05,...){
      # dat la base de donnée deja juméllé au snp
      # snp le snp
      # la variable d'exposition
      # on suppose que la vecteur
      # vvaraj : le vecteur des facteur d'ajustement
      # qrt la variable qui permet de connaitre si nous sommes dans les seuille en cartile pour binaire
      # gma: variable de g?notype de la m?re ? utiliser dans le mod?le de r?gression
      # gca: variable de g?notype de l'enfant ? utiliser dans le mod?le de r?gression
      # seu : seuil pour la decouper la variable
      # le cas binaire o? on precise le seuil de la variable de dichotomisation
      nvar<-cathe_bi2(dat,outc,vexp,seu,quant);nvar<-t(t(nvar));nch<-paste(vexp,".","dch",sep="");colnames(nvar)<-nch
      dat1<-data.frame(dat,nvar)
      # code de d?buggage
      #print(names(dat1))
      #print(gma)
      #print(gca)
      if (!is.null(vvaraj))
         fl0=formula(paste(outc,"~",paste(c(paste(c(nch,gma,gca),collapse="+"),paste(nch,":",gma,sep=""),paste(nch,":",gca,sep=""),paste(vvaraj,collapse="+")),collapse="+"),sep=""))
      else
         fl0=formula(paste(outc,"~",paste(c(paste(c(nch,gma,gca),collapse="+"),paste(nch,":",gma,sep=""),paste(nch,":",gca,sep="")),collapse="+"),sep=""))  
      mod0<-try(Spmlficmcm(fl0,N,gma,gca,DatfE=dat1,typ=typ,...))
      if (inherits(mod0,"try-error"))  vgm1=NA
      else
      {
         z = qnorm(1-alpha/2)
         cof<-mod0$MatR[,colnames(mod0$MatR)=="Estimate"]
         et<-mod0$MatR[,colnames(mod0$MatR)=="Std.Error"]
         Ic1<-cbind(cof-z*et, cof+z*et)
         # or et intervalle de confiance pour le model
         OR<-round(c(1,exp(cof[2])),1)
         Icor<-rbind(c(0,0),round(exp(Ic1)[2,],1))
         ty_w<-cbind(OR,Icor)
         
         # la covariance
         mat1<-mod0$Matv;dg<-diag(mat1);
         
         # or pour mère
         n<-length(cof)
         orjm<-c(1,exp(cof[2]+cof[n-2]));orjme<-round(orjm,1)
         Intm_inf=c(0,round(orjm[2]*exp(-z*sqrt(dg[2]+dg[n-2]+2*mat1[n-2,2])),1))
         Intm_sup=c(0,round(orjm[2]*exp(z*sqrt(dg[2]+dg[n-2]+2*mat1[n-2,2])),1))
         ty_me<-cbind(orjme,Intm_inf,Intm_sup)
         
         # enfant
         orje<-c(1,exp(cof[2]+cof[n-1]));orjen<-round(orje,1)
         Inte_inf=c(0,round(orje[2]*exp(-z*sqrt(dg[2]+dg[n-1]+2*mat1[n-1,2])),1))
         Inte_sup=c(0,round(orje[2]*exp(z*sqrt(dg[2]+dg[n-1]+2*mat1[n-1,2])),1))
         ty_en<-cbind(orjen,Inte_inf,Inte_sup)
         
         # or des interaction genotype mere expo
         ORme<-round(c(1,exp(cof[n-2])),1)
         IORme<-rbind(c(0,0),round(exp(Ic1)[n-2,],1))
         ty_ORme<-cbind(ORme,IORme)
         
         # or des interaction genotype enfant expo
         ORee<-round(c(1,exp(cof[n-1])),1)
         IORee<-rbind(c(0,0),round(exp(Ic1)[n-1,],1))
         ty_ORee<-cbind(ORee,IORee)
         
         
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
         colnames(vgm1)<-NULL;rownames(vgm1)<-NULL
         nnc<-c("Cas1","Cont1","OR1","L1","U1","Cas2","Cont2","OR2","L2","U2","Cas3","Cont3","OR3","L3","U3","OR4","L4","U4","OR5","L5","U5")
         nnr<-c(paste(substr(gm,5,nchar(gm)),"_","Q1",sep=""),paste(substr(gm,5,nchar(gm)),"_","Q2",sep=""))
         colnames(vgm1)<-nnc
         rownames(vgm1)<-nnr
      }
      return(list(matR=vgm1,model1=mod0))
   }
   
# Fonction prenant les sorties de CCMO
   fct_Rsum_CCMO<-function(obj,SNP,alpha=0.05){
                                            # obj objet retourné par CCMO
                                                                   z = qnorm(1-alpha/2)
                                                                   cof<-obj$est
                                                                   cof=cof[-(length(cof)-0:1)]
                                                                   et<-obj$sd
                                                                   Ic1<-cbind(cof-z*et, cof+z*et)
                                                                   # or et intervalle de confiance pour le model
                                                                   OR<-round(c(1,exp(cof[4])),1)
                                                                   Icor<-rbind(c(0,0),round(exp(Ic1)[4,],1))
                                                                   ty_w<-cbind(OR,Icor)

                                                                   # la covariance
                                                                   mat1<-obj$Matv;dg<-diag(mat1);

                                                                   # or pour mère
                                                                   n<-length(cof)
                                                                   orjm<-c(1,exp(cof[4]+cof[n-1]));orjme<-round(orjm,1)
                                                                   Intm_inf=c(0,round(orjm[2]*exp(-z*sqrt(dg[4]+dg[n-1]+2*mat1[n-1,4])),1))
                                                                   Intm_sup=c(0,round(orjm[2]*exp(z*sqrt(dg[4]+dg[n-1]+2*mat1[n-1,4])),1))
                                                                   ty_me<-cbind(orjme,Intm_inf,Intm_sup)

                                                                   # enfant
                                                                   orje<-c(1,exp(cof[4]+cof[n]));orjen<-round(orje,1)
                                                                   Inte_inf=c(0,round(orje[2]*exp(-z*sqrt(dg[4]+dg[n]+2*mat1[n,2])),1))
                                                                   Inte_sup=c(0,round(orje[2]*exp(z*sqrt(dg[4]+dg[n]+2*mat1[n,2])),1))
                                                                   ty_en<-cbind(orjen,Inte_inf,Inte_sup)

                                                                   # or des interaction genotype mere expo
                                                                   ORme<-round(c(1,exp(cof[n-1])),1)
                                                                   IORme<-rbind(c(0,0),round(exp(Ic1)[n-1,],1))
                                                                   ty_ORme<-cbind(ORme,IORme)

                                                                   # or des interaction genotype enfant expo
                                                                   ORee<-round(c(1,exp(cof[n])),1)
                                                                   IORee<-rbind(c(0,0),round(exp(Ic1)[n,],1))
                                                                   ty_ORee<-cbind(ORee,IORee)

                                                                   # recap
                                                                   vgm1<-cbind(ty_w,ty_me,ty_en,ty_ORme,ty_ORee);
                                                                   colnames(vgm1)<-NULL;rownames(vgm1)<-NULL
                                                                   nnc<-c("OR1","L1","U1","OR2","L2","U2","OR3","L3","U3","OR4","L4","U4","OR5","L5","U5")
                                                                   nnr<-c(paste(SNP,"_","Q1",sep=""),paste(SNP,"_","Q2",sep=""))
                                                                   colnames(vgm1)<-nnc
                                                                   rownames(vgm1)<-nnr
                                           return(list(matR=vgm1))
                                           }
                                           
