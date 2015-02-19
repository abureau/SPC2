# Boucle d'analyse de tous les SNPs d'un gène

     fct_gene_CMCM<-function(lis_dat,N,outc,vexp,vvaraj,vec.snp,seu=NULL,minp=0.05,varexclu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               	# varexclu: variable déterminant les observations qu'on veut exclure
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 if (is.null(varexclu)) Datd = Dat[!is.na(Dat[vexp]),]    
                 else
                 {                
                 if (varexclu %in% names(Dat))
                   Datd = Dat[!(is.na(Dat[vexp])|is.na(Dat[varexclu])),]
                 else Datd = Dat[!is.na(Dat[vexp]),]  
                 }               
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                                                    # calcul de fréquence d'allèle  
                                                     gm_tt<-unlist(Datd[gm]);np<-nrow(Datd)
                                                     n0<-length(gm_tt[gm_tt==0])
                                                     n1<-length(gm_tt[gm_tt==1]);n2<-length(gm_tt[gm_tt==2])
                                                     p<-(n2+(1/2)*n1)/np;
                 if (p > minp) llR<-fct_Rsum_CMCM(Datd,N,outc,vexp,gm,gc,vvaraj,seu=seu)
                 else {
                   gma = paste("Cme_","rs",substr(uu,3,nchar(uu)),sep="")
                   gca = paste("Cen_","rs",substr(uu,3,nchar(uu)),sep="")
                   Cme = pmin(gm_tt,1)
                   Cen = pmin(unlist(Datd[gc]),1)
                   Datt = data.frame(Datd,Cme,Cen)
                   names(Datt) = c(names(Datd),gma,gca)
                   # Code de débuggage
                   #print(names(Datt))
                   llR<-fct_Rsum_CMCM(Datt,N,outc,vexp,gm,gc,vvaraj,gma,gca,seu=seu);
                   }
                 if (!all(is.na(llR$matR)))
                   {
                   Tab_Reg<-rbind(Tab_Reg,llR$matR);
                   }
                 modl[[i]]<-llR$model1; 
                 i<-i+1 
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }

     fct_gene<-function(lis_dat,outc,vexp,vvaraj,vec.snp,seu=NULL,minp=0.05,varexclu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               	# varexclu: variable déterminant les observations qu'on veut exclure
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 if (is.null(varexclu)) Datd = Dat[!is.na(Dat[vexp]),]    
                 else
                 {                
                 if (varexclu %in% names(Dat))
                   Datd = Dat[!(is.na(Dat[vexp])|is.na(Dat[varexclu])),]
                 else Datd = Dat[!is.na(Dat[vexp]),]  
                 }               
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                                                    # calcul de fréquence d'allèle  
                                                     gm_tt<-unlist(Datd[gm]);np<-nrow(Datd)
                                                     n0<-length(gm_tt[gm_tt==0])
                                                     n1<-length(gm_tt[gm_tt==1]);n2<-length(gm_tt[gm_tt==2])
                                                     p<-(n2+(1/2)*n1)/np;
                 if (p > minp) llR<-fct_Rsum(Datd,outc,vexp,gm,gc,vvaraj,seu)
                 else {
                   gma = paste("Cme_","rs",substr(uu,3,nchar(uu)),sep="")
                   gca = paste("Cen_","rs",substr(uu,3,nchar(uu)),sep="")
                   Cme = pmin(gm_tt,1)
                   Cen = pmin(unlist(Datd[gc]),1)
                   Datt = data.frame(Datd,Cme,Cen)
                   names(Datt) = c(names(Datd),gma,gca)
                   llR<-fct_Rsum(Datt,outc,vexp,gma,gca,vvaraj,seu);
                   }
                 Tab_Reg<-rbind(Tab_Reg,llR$matR);modl[[i]]<-llR$model1;
                 i<-i+1
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }

     fct_gene_cond<-function(lis_dat,outc,vexp,vvaraj,strate,vec.snp,seu=NULL,minp=0.05,varexclu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               	# varexclu: variable déterminant les observations qu'on veut exclure
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
               	cat(uu,"\n")
                 Dat=lis_dat[[uu]]
                 if (is.null(varexclu)) Datd = Dat[!is.na(Dat[vexp]),]    
                 else
                 {                
                 if (varexclu %in% names(Dat))
                   Datd = Dat[!(is.na(Dat[vexp])|is.na(Dat[varexclu])),]
                 else Datd = Dat[!is.na(Dat[vexp]),]  
                 }               
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                                                    # calcul de fréquence d'allèle  
                                                     gm_tt<-unlist(Datd[gm]);np<-nrow(Datd)
                                                     n0<-length(gm_tt[gm_tt==0])
                                                     n1<-length(gm_tt[gm_tt==1]);n2<-length(gm_tt[gm_tt==2])
                                                     p<-(n2+(1/2)*n1)/np;
                 if (p > minp) llR<-fct_Rcond(Datd,outc,vexp,gm,gc,vvaraj,strate,seu)
                 else {
                   gma = paste("Cme_","rs",substr(uu,3,nchar(uu)),sep="")
                   gca = paste("Cen_","rs",substr(uu,3,nchar(uu)),sep="")
                   Cme = pmin(gm_tt,1)
                   Cen = pmin(unlist(Datd[gc]),1)
                   Datt = data.frame(Datd,Cme,Cen)
                   names(Datt) = c(names(Datd),gma,gca)
                   llR<-fct_Rcond(Datt,outc,vexp,gma,gca,vvaraj,strate,seu);
                   }
                 Tab_Reg<-rbind(Tab_Reg,llR$matR);modl[[i]]<-llR$model1;
                 i<-i+1
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }
                
     fct_gene_2exp_ordre2 <-function(lis_dat,outc,vexp,vexp2,vvaraj,vec.snp,seu=NULL,seu2=NULL,minp=0.05,varexclu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               	# varexclu: variable déterminant les observations qu'on veut exclure
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 if (is.null(varexclu)) Datd = Dat[!is.na(Dat[vexp]),]    
                 else
                 {                
                 if (varexclu %in% names(Dat))
                   Datd = Dat[!(is.na(Dat[vexp])|is.na(Dat[varexclu])),]
                 else Datd = Dat[!is.na(Dat[vexp]),]  
                 }               
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                                                    # calcul de fréquence d'allèle  
                                                     gm_tt<-unlist(Datd[gm]);np<-nrow(Datd)
                                                     n0<-length(gm_tt[gm_tt==0])
                                                     n1<-length(gm_tt[gm_tt==1]);n2<-length(gm_tt[gm_tt==2])
                                                     p<-(n2+(1/2)*n1)/np;
                 if (p > minp) llR<-fct_Rsum_2exp_ordre2(Datd,outc,vexp,vexp2,gm,gc,vvaraj,seu,seu2)
                 else {
                   gma = paste("Cme_","rs",substr(uu,3,nchar(uu)),sep="")
                   gca = paste("Cen_","rs",substr(uu,3,nchar(uu)),sep="")
                   Cme = pmin(gm_tt,1)
                   Cen = pmin(unlist(Datd[gc]),1)
                   Datt = data.frame(Datd,Cme,Cen)
                   names(Datt) = c(names(Datd),gma,gca)
                   llR<-fct_Rsum_2exp_ordre2(Datt,outc,vexp,vexp2,gma,gca,vvaraj,seu,seu2);
                   }
                 Tab_Reg<-rbind(Tab_Reg,llR$matR);modl[[i]]<-llR$model1;
                 i<-i+1
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }
                
     fct_gene_2exp_ordre2_cond <-function(lis_dat,outc,vexp,vexp2,vvaraj,strate,vec.snp,seu=NULL,seu2=NULL,minp=0.05,varexclu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               	# varexclu: variable déterminant les observations qu'on veut exclure
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 if (is.null(varexclu)) Datd = Dat[!is.na(Dat[vexp]),]    
                 else
                 {                
                 if (varexclu %in% names(Dat))
                   Datd = Dat[!(is.na(Dat[vexp])|is.na(Dat[varexclu])),]
                 else Datd = Dat[!is.na(Dat[vexp]),]  
                 }               
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                                                    # calcul de fréquence d'allèle  
                                                     gm_tt<-unlist(Datd[gm]);np<-nrow(Datd)
                                                     n0<-length(gm_tt[gm_tt==0])
                                                     n1<-length(gm_tt[gm_tt==1]);n2<-length(gm_tt[gm_tt==2])
                                                     p<-(n2+(1/2)*n1)/np;
                 if (p > minp) llR<-fct_Rsum_2exp_ordre2_cond(Datd,outc,vexp,vexp2,gm,gc,vvaraj,strate,seu,seu2)
                 else {
                   gma = paste("Cme_","rs",substr(uu,3,nchar(uu)),sep="")
                   gca = paste("Cen_","rs",substr(uu,3,nchar(uu)),sep="")
                   Cme = pmin(gm_tt,1)
                   Cen = pmin(unlist(Datd[gc]),1)
                   Datt = data.frame(Datd,Cme,Cen)
                   names(Datt) = c(names(Datd),gma,gca)
                   llR<-fct_Rsum_2exp_ordre2_cond(Datt,outc,vexp,vexp2,gma,gca,vvaraj,strate,seu,seu2);
                   }
                 Tab_Reg<-rbind(Tab_Reg,llR$matR);modl[[i]]<-llR$model1;
                 i<-i+1
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }
                                  