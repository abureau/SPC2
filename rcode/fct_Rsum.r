  # 1.9 la fonction des quartiles 
    cathe_bi<-function(dat,var1,seu=NULL){
                                    if(is.null(seu)==TRUE){
                                                           var<-dat[,var1]
                                                           seu=quantile(var,probs=c(0.75),na.rm =TRUE)
                                                           qrt1<-ifelse(var>=seu,1,0)+ifelse(is.na(var)==TRUE,NA,0)
                                                           }else{
                                                           var<-dat[,var1]
                                                           qrt1<-ifelse(var>=seu,1,0)+ifelse(is.na(var)==TRUE,NA,0)
                                                           }
                                    return(qrt1)
                                    }                                   
# 1.11 La fonction qui resume pour un SNPs une seule exposition

   fct_Rsum<-function(dat,outc,vexp,gm,gc,vvaraj,gma=gm,gca=gc,seu=NULL){
                                            # dat la base de donnée deja juméllé au snp
                                            # snp le snp
                                            # la variable d'exposition
                                            # on suppose que la vecteur
                                            # vvaraj : le vecteur des facteur d'ajustement
                                            # qrt la variable qui permet de connaitre si nous sommes dans les seuille en cartile pour binaire
                                            # gma: variable de g�notype de la m�re � utiliser dans le mod�le de r�gression
                                            # gca: variable de g�notype de l'enfant � utiliser dans le mod�le de r�gression
                                            # seu : seuil pour la decouper la variable
                                            dat0<-dat[dat[,outc]==0,];vexp0<-dat0[vexp];vexp1<-dat[vexp]
                                            # le cas binaire o� on precise le seuil de la variable de dichotomisation
                                               nvar<-cathe_bi(dat,vexp,seu);nvar<-t(t(nvar));nch<-paste(vexp,".","dch",sep="");colnames(nvar)<-nch
                                                                   dat1<-data.frame(dat,nvar)
                                                                   fl0=formula(paste(outc,"~",paste(c(paste(c(nch,gma,gca),collapse="+"),paste(nch,":",gma,sep=""),paste(nch,":",gca,sep=""),paste(vvaraj,collapse="+")),collapse="+"),sep=""))
                                                                   mod0<-glm(fl0,data=dat1,family=binomial)
                                                                   cof<-coef(mod0)
                                                                   Ic1<-confint(mod0)
                                                                   # or et intervalle de confiance pour le model
                                                                   OR<-round(c(1,exp(cof[2])),1)
                                                                   Icor<-rbind(c(0,0),round(exp(Ic1)[2,],1))
                                                                   ty_w<-cbind(OR,Icor)

                                                                   # la covariance
                                                                   mat1<-vcov(mod0);dg<-diag(mat1);

                                                                   # or pour mère
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

                                                                   # or des interaction genotype mere expo
                                                                   ORme<-round(c(1,exp(cof[n-1])),1)
                                                                   IORme<-rbind(c(0,0),round(exp(Ic1)[n-1,],1))
                                                                   ty_ORme<-cbind(ORme,IORme)

                                                                   # or des interaction genotype enfant expo
                                                                   ORee<-round(c(1,exp(cof[n])),1)
                                                                   IORee<-rbind(c(0,0),round(exp(Ic1)[n,],1))
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
                                                                   nnc<-c("Cas1","Cont1","OR1","Ic1","Ic2","Cas2","Cont2","OR2","Ic1","Ic2","Cas3","Cont3","OR3","Ic1","Ic1","OR4","Ic1","Ic2","OR5","Ic1","Ic2")
                                                                   nnr<-c(paste(substr(gm,5,nchar(gm)),"_","Q1",sep=""),paste(substr(gm,5,nchar(gm)),"_","Q2",sep=""))
                                                                   colnames(vgm1)<-nnc
                                                                   rownames(vgm1)<-nnr

                                           return(list(matR=vgm1,model1=mod0))
                                           }