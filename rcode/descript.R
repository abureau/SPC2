fucomp_descC = function(dat,ve.vxp,vout){
                                         tabdes<-NULL;tabna<-NULL
                                         for(cc in ve.vxp){
                                                           ll<-desc_pvc(dat,cc,vout)
                                                           tabdes<-rbind(tabdes,ll$vcfr1)
                                                           tabna<-rbind(tabna,ll$vc.na)
                                                           }
                                        ncc<-c("Mean.cas","(sd)","nbr.na","Mean.controls","(sd)","3rd quartile","nbr.na","p")
                                        colnames(tabdes)<-ncc
                                        nam<-paste(ve.vxp,"(","z=1",")",sep="")
                                        rownames(tabdes)<-nam
                                        return(list(tabdes=tabdes,tabna=tabna))
                                        }
desc_pvc = function(dat,vxp,vout){
                                    
                                   efn<-dim(dat)[1] 
                                   dat0<-dat[dat[vout]==0,];dat1<-dat[dat[vout]==1,];
                                   # nombre na
                                   datcsna<-dat1[is.na(dat1[,vxp])!=TRUE,]
                                   nbrc<-dim(dat1[is.na(dat1[,vxp])==TRUE,])[1]
                                   datcona<-dat0[is.na(dat0[,vxp])!=TRUE,]
                                   nbrt<-dim(dat0[is.na(dat0[,vxp])==TRUE,])[1]
#                                   ttp = min(c(datcsna[,vxp]+0,datcona[,vxp]))
                                   ttp = wilcox.test(log(datcsna[,vxp]),log(datcona[,vxp]))$p.value
                                   vcfr1<-c(mean(datcsna[,vxp]),sd(datcsna[,vxp]),nbrc,mean(datcona[,vxp]),sd(datcona[,vxp]),quantile(datcona[,vxp],0.75),nbrt,ttp)
                                   nbnac<-dim(dat1)[1]-dim(datcsna)[1]
                                   nbnco<-dim(dat0)[1]-dim(datcona)[1]
                                   vc.na<-c(nbnac,nbnco)
                                   return(list(vcfr1=vcfr1,vc.na=vc.na))
                                   }
