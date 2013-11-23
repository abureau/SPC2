# Boucle d'analyse de tous les SNPs d'un gène

     fct_gene_CMCM<-function(lis_dat,N,outc,vexp,vvaraj,vec.snp,seu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 Datd = Dat[!is.na(Dat[vexp]),]
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                 llR<-fct_Rsum_CMCM(Datd,N,outc,vexp,gm,gc,vvaraj,seu=seu);
                 if (!all(is.na(llR$matR)))
                   {
                   Tab_Reg<-rbind(Tab_Reg,llR$matR);
                   }
                 modl[[i]]<-llR$model1; 
                 i<-i+1 
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }

     fct_gene<-function(lis_dat,outc,vexp,vvaraj,vec.snp,seu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               #nsnp<-names(datCpx2)[-1]
               # v_snp<-nsnp[regexpr("rs",nsnp)==-1]
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 Datd = Dat[!is.na(Dat[vexp]),]
                 gm<-paste("gme_","rs",substr(uu,3,nchar(uu)),sep="")
                 gc<-paste("gen_","rs",substr(uu,3,nchar(uu)),sep="")
                 llR<-fct_Rsum(Datd,outc,vexp,gm,gc,vvaraj,seu);
                 Tab_Reg<-rbind(Tab_Reg,llR$matR);modl[[i]]<-llR$model1;
                 i<-i+1
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }
                 