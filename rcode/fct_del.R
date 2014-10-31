     fct_del<-function(lis_dat,outc,vexp,vvaraj,vec.snp,seu=NULL)
               {#datCpx2 : base de cimplex; datId : la basee des id; datexp la base des exposants vec.snp vecteur de snp
               i=1; Tab_Reg<-NULL;modl<-list();
               
               for(uu in vec.snp){
                 Dat=lis_dat[[uu]]
                 Datd = Dat[!is.na(Dat[vexp]),]
                 gm<-paste("gme_",uu,sep="")
                 gc<-paste("gen_",uu,sep="")
					llR<-fct_Rsum(Datd,outc,vexp,gm,gc,vvaraj,seu)
                 Tab_Reg<-rbind(Tab_Reg,llR$matR);modl[[i]]<-llR$model1;
                 i<-i+1
                 }
                return(list(Tab_Reg=Tab_Reg,modl=modl))
                }
