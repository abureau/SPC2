#==============================================================================#
#           Programme de management des données                                #
#==============================================================================#

# PAckage 
           library(gdata)
           library(vcd)
           library(genetics)
  
 
 #======== les fonction necessaire pour le management des données =============#
 
 # 1.1 la fonction qui transorme le snp de la base g..t ou c..t , etc en une variable qui vaut o,1/2,1 NA(valeure manquante) (compte le nombre l'allèle mineur et divise par deux ) 
 dec_snp<-function(rsn,dat){
                            # rsn c'est le num rs
                            # dat : la base de donnée
                            bse<-dat[,rsn];vv<-table(bse)
                            if(("N"%in%bse)==TRUE){nn1<-names(vv);vv1<-vv[-which(nn1=="N")]}else{vv1<-vv}
                            bs_min<-substr(names(which.min(vv1)),1,1)
                            vec<-NULL
                            for(u in bse){vec<-rbind(vec,c(substr(u,1,1),substr(u,3,3)))}
                            
                            #traitement all�le 1
                            vid_vec1<-vec[,1]%in%""
                            na_vec1<-vec[,1]%in%"N"
                            maf_vec1<-vec[,1]%in%bs_min
                            #traitement all�le 2
                            vid_vec2<-vec[,2]%in%""
                            na_vec2<-vec[,2]%in%"N"
                            maf_vec2<-vec[,2]%in%bs_min
                       
                            Al1<-ifelse(is.na(vec[,1])|vid_vec1==TRUE|na_vec1==TRUE,NA,ifelse(maf_vec1==TRUE,1,0))
                            Al2<-ifelse(is.na(vec[,1])|vid_vec2==TRUE|na_vec2==TRUE,NA,ifelse(maf_vec2==TRUE,1,0))
                            sap<-(Al1+Al2);snp<-t(t(sap))
                            rsn1<-paste("rs",substr(rsn,3,nchar(rsn)),sep="")
                            tff<-data.frame(snp,bse);nams<-c(rsn1,rsn)
                            names(tff)<-nams
                            return(list(snp=tff,ma=bs_min))
                            }

# 1.2 la fonction qui transforme toute la table 
  TRansf_tab<-function(dat){
                           # dat la base 
                           vec<-names(dat)[-1]
                           SAMPLE<-dat[,1]
                           Ntsnp<-data.frame(SAMPLE);
                           for(u in vec){Ntsnp<-cbind(Ntsnp,dec_snp(u,dat)$snp)}
                           return(Ntsnp)
                           }
# 1.3 la fonction qui retrouve la mère et l'enfant 
  mere_enf<-function(dat1,rsn,dat2){
                                    # dat1 est la data.frame qui provient de la fonction TRansf_tab et rsn le nom du snp
                                    # dat2 la table qui contient les mère et les enfant 
                                    # rsn le snp qu'on regarde
                                    rsn1<-paste("rs",substr(rsn,3,nchar(rsn)),sep="")
                                    xsd<-dat1[[rsn1]] 
                                    xs<-dat1[[rsn]]
                                    zs<-dat1[["SAMPLE"]]
                                    dat3<-data.frame(zs,xs,xsd)
                                    names(dat3)<-c("SAMPLE",rsn,rsn1)
                                    codm<<-dat2[["CODMER"]]
                                    code<<-dat2[["CODENF"]]
                                    nosp1<-dat2[["NOSP1"]]
                                    dat4<-data.frame(nosp1,codm,code)
                                    Nvtm<-NULL;
                                    for(u in codm){
                                                  if(dim(dat3[dat3$SAMPLE==u,])[1]==0){
                                                                                       nam<-names(dat3[dat3$SAMPLE==u,]);tt<-data.frame(u,"",NA);names(tt)<-nam;
                                                                                       Nvtm<-rbind(Nvtm,tt)
                                                                                       }else{
                                                                                             Nvtm<-rbind(Nvtm,dat3[dat3$SAMPLE==u,])
                                                                                             }
                                                  }
                                    Nvte<-NULL;
                                    for(v in code){
                                                  if(dim(dat3[dat3$SAMPLE==v,])[1]==0){
                                                                                       nam<-names(dat3[dat3$SAMPLE==v,]);tt1<-data.frame(v,"",NA);names(tt1)<-nam;
                                                                                       Nvte<-rbind(Nvte,tt1)
                                                                                       }else{
                                                                                             Nvte<-rbind(Nvte,dat3[dat3$SAMPLE==v,])
                                                                                             }
                                                  }
                                  # table des génotype 
                                  ngm1<-paste("gme_",rsn,sep="")
                                  ngm2<-paste("gme_",rsn1,sep="")
                                  nge1<-paste("gen_",rsn,sep="")
                                  nge2<-paste("gen_",rsn1,sep="")
                                  gnen<-Nvte[,2:3];gnme<-Nvtm[,2:3];
                                  prenam<-names(dat4);nomc<-prenam[-1]
                                  datt<-data.frame(Nvtm,Nvte)
                                  names(datt)<-c(nomc[1],ngm1,ngm2,nomc[2],nge1,nge2)
                                  Tabge<-merge(dat4,datt,by=c(nomc[1],nomc[2]))
                                  return(Tabge)
                                  }
                                  
# modfication de la fonction pour les deletions (fonction qui ne marche pas pour l'insta)
 mere_enfD<-function(dat1,rsn,dat2){
                                    # dat1 est la data.frame qui provient de la fonction TRansf_tab et rsn le nom du snp
                                    # dat2 la table qui contient les mère et les enfant 
                                    # rsn le snp qu'on regarde
                                    #xsd<-dat1[rsn] 
                                    xs<-dat1[rsn]
                                    zs<-dat1["X"]
                                    dat3<-data.frame(zs,xs)
                                    names(dat3)<-c("SAMPLE",rsn)
                                    codm<-as.character(dat2[["CODMER"]])
                                    code<-as.character(dat2[["CODENF"]])
                                    nosp1<-as.character(dat2[["NOSP1"]])
                                    dat4<-data.frame(nosp1,codm,code)
                                    Nvtm<-NULL;
                                    for(u in codm){
                                                  if(dim(dat3[dat3$SAMPLE==u,])[1]==0){
                                                                                       nam<-names(dat3[dat3$SAMPLE==u,]);tt<-data.frame(u,NA);names(tt)<-nam;
                                                                                       Nvtm<-rbind(Nvtm,tt)
                                                                                       }else{
                                                                                             if(dim(dat3[dat3$SAMPLE==u,])[1]>1){Nvtm<-rbind(Nvtm,dat3[dat3$SAMPLE==u,][1,])
                                                                                                                                 }else{ Nvtm<-rbind(Nvtm,dat3[dat3$SAMPLE==u,])
                                                                                                                                      }
                                                                                             }
                                                  }
                                    Nvte<-NULL;
                                    for(v in code){
                                                  if(dim(dat3[dat3$SAMPLE==v,])[1]==0){
                                                                                       nam<-names(dat3[dat3$SAMPLE==v,]);tt1<-data.frame(v,NA);names(tt1)<-nam;
                                                                                       Nvte<-rbind(Nvte,tt1)
                                                                                       }else{
                                                                                             if(dim(dat3[dat3$SAMPLE==v,])>1){
                                                                                                                              Nvte<-rbind(Nvte,dat3[dat3$SAMPLE==v,][1,])
                                                                                                                              }else{
                                                                                                                                    Nvte<-rbind(Nvte,dat3[dat3$SAMPLE==v,])
                                                                                                                                    }
                                                                                             }
                                                  }
                                  # table des génotype 
                                  ngm1<-paste("gme_",rsn,sep="")
                                  #ngm2<-paste("gme_",rsn1,sep="")
                                  nge1<-paste("gen_",rsn,sep="")
                                  #nge2<-paste("gen_",rsn1,sep="")
                                  Nvtm1<-Nvtm[,2];Nvte1<-Nvte[,2]
                                  prenam<-names(dat4);naptou<-c(prenam,ngm1,nge1)
                                  Tabge<-data.frame(dat4,Nvtm1,Nvte1)
                                  names(Tabge)<-naptou
                                  return(Tabge)
                                  }

# 1.4 fonction du test de génotypage (elle teste si le génotype de l'enfant et la mère est cpmpatible. valeur d'incompatibilité au moins 90  )
  
  Test_geno<-function(datf){
                            #datf : l'objet de la fonction "mere_enf" 
                            fct_test<-function(v){
                                                  va<-1
                                                 # v : vecteur ne doit pas contenir de donn�es manquantes
                                                 if (all(!is.na(v)))
                                                 { 
                                                 gm<-v[1];ge<-v[2]; 
                                                 if(gm==0 & ge==2){va<-90}
                                                 if(gm==2 & ge==0){va<-90}
                                                 }
                                                 return(va)
                                                 }
                            v1<-datf[,5];v2<-datf[,7];
                            tgm<-cbind(v1,v2);test<-apply(tgm,1,fct_test)
                            datf1<-data.frame(datf,test)
                            return(datf1)
                            }
# 1.4 la fonction qui effectue le merge sur les données d'exposition

  fct_mergtt<-function(dat1,dat2,rsn,datz){
                                           # dat1 est la data.frame qui provient de la fonction TRansf_tab et rsn le nom du snp
                                           # dat2 la table qui contient les mère et les enfant 
                                           # rsn le snp qu'on regarde
                                           # datz table d'exposition 
                                           datgm<-mere_enf(dat1,rsn,dat2)
                                           # le test des compatibilité génétique 
                                           datg<-Test_geno(datgm)
                                           tabp<-merge(datg,datz,by=c("nosp1")) 
                                           return(tabp)
                                           }
# petite modification pour les délétion 
  fct_mergttD<-function(dat1,dat2,rsn,datz){
                                           # dat1 est la data.frame qui provient de la fonction TRansf_tab et rsn le nom du snp
                                           # dat2 la table qui contient les mère et les enfant 
                                           # rsn le snp qu'on regarde
                                           # datz table d'exposition 
                                           datgm<-mere_enfD(dat1,rsn,dat2)
                                           # le test des compatibilité génétique 
                                           # datg<-Test_geno(datgm)
                                           tabp<-merge(datgm,datz,by=c("nosp1")) 
                                           return(tabp)
                                           }

# 1.7 La fonction qui prepare toutes les tables et les concerve dans une liste 
 # 1.7.1 conserve les tables avec avec les données mananquantes sur le génotype de l'enfant  
        latSnpWMV<-function(dat1,dat2,vec.rsn,datz){
                                           #dat1 : data.frame obtenue à partir de transfer
                                           #dat2 : data.frame des id mere enfant
                                           #vec.rsn : vecteur des snp
                                           # datz : data.frame des covariable 
                                           listSNP<-list()
                                           i=1;
                                           for(ss in vec.rsn){ # debugging code
                                                              print(ss)
                                                              tab<-fct_mergtt(dat1,dat2,ss,datz)
                                                              vargme<-paste("gme_rs",substr(ss,3,nchar(ss)),sep="")
                                                              tab1<-tab[!is.na(tab[vargme]),]
                                                              listSNP[[i]]<-tab1[tab1["test"]!=90,];i=i+1}
                                           names(listSNP)<-vec.rsn
                                           return(listSNP)
                                           }
  #1.7.2 conserve les données sans les données manquantes 
        latSnpNMV<-function(dat1,dat2,vec.rsn,datz){
                                           #dat1 : data.frame obtenue à partir de transfer
                                           #dat2 : data.frame des id mere enfant
                                           #vec.rsn : vecteur des snp
                                           # datz : data.frame des covariable 
                                           listSNP<-list()
                                           i=1;
                                           for(ss in vec.rsn){tab<-fct_mergtt(dat1,dat2,ss,datz)
                                                              vargme<-paste("gme_rs",substr(ss,3,nchar(ss)),sep="")
                                                              vargen<-paste("gen_rs",substr(ss,3,nchar(ss)),sep="")
                                                              tab1<-tab[!is.na(tab[vargme])&!is.na(tab[vargen]),]
                                                              listSNP[[i]]<-tab1[tab1["test"]!=90,];i=i+1}
                                           names(listSNP)<-vec.rsn
                                           return(listSNP)
                                           }
                                           
# 1.7 La fonction qui prepare toutes les tables et les concerve dans une liste pour les delession
# 1.7.1 conserve les tables avec avec les données mananquantes sur le génotype de l'enfant  
        latSnpWMVD<-function(dat1,dat2,vec.rsn,datz){
                                           #dat1 : data.frame obtenue à partir de transfer
                                           #dat2 : data.frame des id mere enfant
                                           #vec.rsn : vecteur des snp
                                           #datz : data.frame des covariable 
                                           listSNP<-list()
                                           i=1;
                                           for(ss in vec.rsn){tab<-fct_mergttD(dat1,dat2,ss,datz)
                                                              vargme<-paste("gme_",ss,sep="")
                                                              tab1<-tab[tab[vargme]<2,]
                                                              listSNP[[i]]<-tab1[is.na(tab1[vargme])!=TRUE,];i=i+1}
                                           names(listSNP)<-vec.rsn
                                           return(listSNP)
                                           }
  #1.7.2 conserve les données sans les données manquantes 
        latSnpNMVD<-function(dat1,dat2,vec.rsn,datz){
                                           #dat1 : data.frame obtenue à partir de transfer
                                           #dat2 : data.frame des id mere enfant
                                           #vec.rsn : vecteur des snp
                                           # datz : data.frame des covariable 
                                           listSNP<-list()
                                           i=1;
                                           for(ss in vec.rsn){tab<-fct_mergttD(dat1,dat2,ss,datz)
                                                              vargme<-paste("gme_",ss,sep="")
                                                              vargen<-paste("gen_",ss,sep="")
                                                              tab1<-tab[tab[vargme]<2,] 
                                                              tab2<-tab1[is.na(tab1[vargme])!=TRUE,]
                                                              tab3<-tab2[tab2[vargen]<2,]
                                                              listSNP[[i]]<-tab3[is.na(tab3[vargen])!=TRUE,];i=i+1}
                                           names(listSNP)<-vec.rsn
                                           return(listSNP)
                                           }

# =================================== fin ======================================
