if(!require(sequoia)) install.packages("sequoia")
if(!require(dplyr)) install.packages("dplyr")
if(!require(openxlsx)) install.packages("openxlsx")
if(!require(dplyr)) install.packages("dplyr")


GPEDI.genealogy <- function(
                   data, #SNP型文件，第一列为个体名称，每一列为每一个个体的SNP位点
				   pedoutput = 0,
				   Tine = 0,
				   filename
				   ){
file = data
pairsss=rownames(file)
pairs = data.frame(ID1 = character(),
			       ID2 = character(),
			       stringsAsFactors = FALSE
			       )
for (i in 1:(length(pairsss)-1)){
     for (j in (i+1):length(pairsss)){
		 pair<- data.frame(ID1 = pairsss[i],
		 ID2 = pairsss[j]) 
		 pairs<- rbind(pairs, pair)
		}
}
PairL = CalcPairLL(Pairs=pairs,GenoM = as.matrix(file), Err = 1e-04, Plot=FALSE)
ped =NULL
ID = unique(rownames(file))				 
PO=PairL[-c(which(PairL$TopRel!="PO")),]	
##################################################			 
if(nrow(PO)!=0){
	 PO_1 = PO
	 familytree = 0
	 H_all = NULL
	 while(nrow(PO_1) != 0){
	 familytree = familytree+1
	 H = PO_1[familytree,]
	 H1 = H[[1]]
	 H2 = H[[2]]
	 H_c = H
	 #PO_1 = PO_1[!(PO_1$ID1==H1 & PO_1$ID2==H2),]
	 #nrow(PO_1)
	 h = NULL
	 H_all[[familytree]] = H_c
	 #familytree
	 #H_all[[familytree]]
	 while(nrow(H)!= 0){
		h = NULL
		for(i in 1:nrow(H)){
		H1 = H[i,1]
		H2 = H[i,2]
		#H1
		#H2
		PO_1 = PO_1[!(PO_1$ID1==H1 & PO_1$ID2==H2),]
		PO_1 = PO_1[!(PO_1$ID2==H1 & PO_1$ID1==H2),]
		H_3=PO_1[c(which(PO_1$ID1==H1)),]
		H_4=PO_1[c(which(PO_1$HID2==H1)),]
		H_5=PO_1[c(which(PO_1$ID1==H2)),]
		H_6=PO_1[c(which(PO_1$ID2==H2)),]
		h[[i]]= rbind(H_3,H_4,H_5,H_6)
		h[[i]] = unique(h[[i]])
		#nrow(PO_1)
		#nrow(H)
							}
		 H= do.call(rbind, h)
		 H=unique(H)
		 #nrow(H)
		 #nrow(H_all[[1]])
		 #nrow(PO_1)
		 H_all[[familytree]]= rbind(H,H_all[[familytree]])
		 H_all[[familytree]]=unique(H_all[[familytree]])
		 #nrow(H_all[[1]])
		}
		PO_1 = PO_1 [!is.na(PO_1 $ID1),]
	}				 
	#############################################################
	for(v in 1:length(H_all)){
		 PO_H = H_all[[v]]
		 PO_M = c(PO_H$ID1,PO_H$ID2)
		 PO_M = unique(PO_M)
		 PairL_1 = PairL[PairL$ID1 %in% PO_M, ]
		 PairL_1 = PairL_1[PairL_1$ID2 %in% PO_M, ]
		 FS_1 = PairL_1$TopRel
		 PO_T = PO_H
		 FS=PairL_1[-c(which(PairL_1$TopRel!="FS")),]
		 while(nrow(FS)!=0){
			 FS = FS[,1:2]
			 sson_1=FS[1,]$ID1
			 sson_3=c(FS[c(which(FS$ID1==sson_1)),])
			 sson_4=c(FS[c(which(FS$ID2==sson_1)),])
			 sson_2 = c(unlist(sson_3),unlist(sson_4))
			 sson_2= unique(sson_2)
			 sson = setdiff(sson_2,sson_1)
			 PO_sson_1 = PO_T[PO_T$ID1 == sson_1,]
			 PO_sson_1 = c(PO_sson_1$ID1,PO_sson_1$ID2)
			 PO_sson_1 = setdiff(PO_sson_1,sson_1)
			 #sson_2 = [which(sson_2!=sson_1)]
			 for(j in 1:length(sson)){
			 PO_sson = rbind(PO_T[PO_T$ID1 == sson[j],],PO_T[PO_T$ID2 == sson[j],])
			 PO_T = anti_join(PO_T, PO_sson)
			 PO_sson = PO_sson %>%
						mutate(across(everything(), ~ ifelse(. %in% PO_sson_1, "A1", .)))
			 PO_T = rbind(PO_T,PO_sson)
			  }
			  FS = FS[!(FS$ID1 %in% sson), ]
			  FS = FS[!(FS$ID2 %in% sson), ]	
			}
		###############################################
		 HS=PairL_1[-c(which(PairL_1$TopRel!="HS")),]
		 while(nrow(HS)!=0){
			 HS = HS[,1:2]
			 sson_1=HS[1,]$ID1
			 sson_3=c(HS[c(which(HS$ID1==sson_1)),])
			 sson_4=c(HS[c(which(HS$ID2==sson_1)),])
			 sson_2 = c(unlist(sson_3),unlist(sson_4))
			 sson_2= unique(sson_2)
			 sson = setdiff(sson_2,sson_1)
			 PO_sson_1 = PO_T[PO_T$ID1 == sson_1,]
			 PO_sson_1 = c(PO_sson_1$ID1,PO_sson_1$ID2)
			 PO_sson_1 = setdiff(PO_sson_1,sson_1)
			 #sson_2 = [which(sson_2!=sson_1)]
			 for(j in 1:length(sson)){
			 PO_sson = rbind(PO_T[PO_T$ID1 == sson[j],],PO_T[PO_T$ID2 == sson[j],])
			 PO_T = anti_join(PO_T, PO_sson)
			 PO_sson = PO_sson %>%
						mutate(across(everything(), ~ ifelse(. %in% PO_sson_1, "A1", .)))
			 PO_T = rbind(PO_T,PO_sson)
			  }
			  HS = HS[!(HS$ID1 %in% sson), ]
			  HS = HS[!(HS$ID2 %in% sson), ]	
			}
		####################################################
		ID_1 = c(PO_T$ID1)
		 ID_2 = c(PO_T$ID2)
		 ID_0=c(ID_1,ID_2)
		 name=NULL
		 min_ID=table(ID_0)
		 min_I=1
		 F1=names(min_ID[min_ID == min_I])
		 ID = F1
		 k=1
		 name [[k]]=ID
		 while(length(ID)!=0){
			 k=k+1
			 ID_1= PO_H[PO_H$ID1 %in% ID, ]
			 ID_1= c(ID_1$ID2)
			 ID_2= PO_H[PO_H$ID2 %in% ID, ]
			 ID_2= c(ID_2$ID1)
			 ID_0=c(ID_1,ID_2)
			 ID_0=unique(ID_0)
			 #ID_C=name[name %in% ID_0]
			 ID_0=ID_0[!(ID_0 %in% ID)]
			 ID_all=unlist(name)
			 ID = ID_0[!(ID_0 %in% ID_all)]
			 if(length(ID)!=0){
				 name[[k]]=ID
				 }else {
						}

			} 
		#########################################第一次确定代数
		 #ped_1 = matrix(,long,4)
		 #colnames(ped_1)=c("ID","sire","dam","sex")
		 #j=6
		 HS=PairL_1[-c(which(PairL_1$TopRel!="HS")),]
		 name_p=NULL
		 ped_m=NULL
		 j = length(name)
		 B = name[[j]]
		 p=1
		 name_p[[p]]=B
		 while (!all(is.na(B))){						  
			 B = na.omit(B)
			 B = B[B != "NA"]
			 name_p2=NULL
			 
			 for (k in 1:length(B)){
			  C = B[k]
			  ID_1= PO[which(PO$ID1==C), ]
			  ID_1=unique(ID_1)
			  ID_1= c(ID_1$ID2)
			  ID_2= PO[which(PO$ID2==C), ]
			  ID_2=unique(ID_2)
			  ID_2= c(ID_2$ID1)
			 # ID_father= (PO[which(PO$Sex1==2), ])[1,1]
			 # ID_mather= (PO[which(PO$Sex1==1), ])[1,1]												 
			 ID_all2=unlist(name_p)
			 ID_1=ID_1[!(ID_1 %in% ID_all2)]
			 ID_2=ID_2[!(ID_2 %in% ID_all2)]												 						                         
			 m=c(C,ID_1,ID_2)
			 if(length(m)!=3){
				 if(length(m)<3){																					 
					 m[(length(m)+1):3]=NA
					 m=matrix(m, nrow = 1,ncol=3 )
					 ped_m=rbind(ped_m,m)
					 while(length(ID_1)==0 & length(ID_2)==0){
													 ID_1=NA
													 ID_2=NA
													 }
					 name_p2[[k]]=c(ID_1,ID_2)
					}else{
						  m=c(C,NA,NA)
						  m=matrix(m, nrow = 1,ncol=3 )
						  ped_m=rbind(ped_m,m)
						  name_p2[[k]]=NA
						  }																
	  
					}else{
					 #if(length(m)=1)
						 m=matrix(m, nrow = 1,ncol=3 )
						 ped_m=rbind(ped_m,m)
						 while(length(ID_1)==0 & length(ID_2)==0){
																 ID_1=NA
																 ID_2=NA
																 }
						 name_p2[[k]]=c(ID_1,ID_2)					
						}							   
				}
			 p=p+1
			 name_p[[p]]=c(unlist(name_p2))
			 B=name_p[[p]]		 
			}	
			 ped_m=unique(ped_m)
			 sex_ID = matrix(0,length(ped_m[,1]),1)
			 S = c(PO$Sex1,PO$Sex2)
		#############################################################		 
			 if (!is.null(ped_m)){
			 for(i in 1:nrow(ped_m)){
				ID_sex = ped_m[i,1]
				ID_sex_row1 = PO[which(PO$ID1==ID_sex), ]
				ID_sex_row2 = PO[which(PO$ID2==ID_sex), ]
				ID_sex1 = unique(c((ID_sex_row1$Sex1),(ID_sex_row2$Sex2)))		  
				if(length(ID_sex1)==1){
				sex_ID[i,1] = ID_sex1
				}else{
				sex_ID[i,1] = NA
				}
			 }		
			 sex_ID[sex_ID == 3] = NA
			 sex_ID[sex_ID == 1] = 0
			 sex_ID[sex_ID == 2] = 1
			 ped_m = cbind(ped_m,sex_ID)
			 if(3 %in% S){
				 colnames(ped_m)=c("ID","parents1","parents2","gender")
				 }else{
					colnames(ped_m)=c("ID","sire","dam","gender")
					for(i in 1:length(ped_m[,2])){
					 m = ped_m[i,2]
					 sex_m1 = PO[which(PO$ID1==m), ] 
					 sex_m2 = PO[which(PO$ID2==m), ]
					 ID_sex1 = unique(c(unlist(ID_sex_row1$Sex1),unlist(ID_sex_row2$Sex2)))
					 if(ID_sex1==2){
						 ped_m[i,2] = ped_m[i,4]
						 ped_m[i,3] = m
						 }else{
							  }
					}

					 }
			}else{ 
				 ID_H = unique(c(PO_H$ID1,PO_H$ID2))
				 ped_m = matrix(NA,length(ID_H),4)
				 colnames(ped_m)=c("ID","parents","parents","gender")
				 ped_m[,1] = ID_H
				   }
		############################################################
		 HSFS_ID = matrix(0,length(ped_m[,1]),4)
		 HS=PairL[c(which(PairL$TopRel=="HS")),]
		 FS=PairL[c(which(PairL$TopRel=="FS")),]
		 if (nrow(FS)!=0){
			 for(i in 1:nrow(ped_m)){
				 ID_FS = ped_m[i,1]
				 ID_FS_row1 = FS[which(FS$ID1==ID_FS), ]
				 ID_FS_row2 = FS[which(FS$ID2==ID_FS), ]
				 ID_FS_1 = unique(c((ID_FS_row1$ID2),(ID_FS_row2$ID1)))		  
					 if(length(ID_FS_1)!=0){
					 HSFS_ID[i,1] = length(ID_FS_1)
					 HSFS_ID[i,2] = paste(ID_FS_1, collapse = ",")
					 }else{
						 HSFS_ID[i,2] = NA
						 }
			}
		 }
		 #########################################################
		 if (nrow(HS)!=0){
			 for(i in 1:nrow(ped_m)){
			 ID_HS = ped_m[i,1]
			 ID_HS_row1 = FS[which(HS$ID1==ID_HS), ]
			 ID_HS_row2 = FS[which(HS$ID2==ID_HS), ]
			 ID_HS_1 = unique(c((ID_HS_row1$ID2),(ID_HS_row2$ID1)))		  
			 if(length(ID_HS_1!=0)){
			 HSFS_ID[i,3] = length(ID_HS_1)
			 HSFS_ID[i,4] = paste(ID_HS_1, collapse = ",")
			 }else{
			 HSFS_ID[i,4] = NA
			}
		 }
		 }
		 HSFS_ID[HSFS_ID[, 2] == 0, 2] = NA
		 HSFS_ID[HSFS_ID[, 4] == 0, 4] = NA
		 colnames(HSFS_ID) = c("FSnum","FS","HSnum","HS")
		 ped_m = cbind(ped_m,HSFS_ID)
		############################################################
		 ped_m=unique(ped_m)
		 familynames=matrix(v, nrow =nrow(ped_m),ncol=1)
		 colnames(familynames) = "familytrees"
		 ped_m = cbind(ped_m,familynames)
		 ped[[v]] = ped_m
	 }
	 ped_V = do.call(rbind, ped)
	 ped_V = as.data.frame(ped_V)
	 ped_V = ped_V[!is.na(ped_V$ID), ]
	 ID_gen = rownames(file)
	 ID_ped = ped_V[,1]
	 result_ID=setdiff(ID_gen, ID_ped)
	 ped_B = matrix(0,length(result_ID),ncol(ped_V))
	 ped_B[,1] = result_ID
	 ped_B[, 2:3] = NA
	###########################################
	 HS=PairL[c(which(PairL$TopRel=="HS")),]
	 FS=PairL[c(which(PairL$TopRel=="FS")),]
	 if (nrow(FS)!=0){
		 for(i in 1:nrow(ped_B)){
			 ID_FS = ped_B[i,1]
			 ID_FS_row1 = FS[which(FS$ID1==ID_FS), ]
			 ID_FS_row2 = FS[which(FS$ID2==ID_FS), ]
			 ID_FS_1 = unique(c((ID_FS_row1$ID2),(ID_FS_row2$ID1)))		  
				 if(length(ID_FS_1)!=0){
				 ped_B[i,5] = length(ID_FS_1)
				 ped_B[i,6] = paste(ID_FS_1, collapse = ",")
				 }else{
					 ped_B[i,6] = NA
					 }
		}
	 }
	 #########################################################
	 if (nrow(HS)!=0){
		 for(i in 1:nrow(ped_B)){
		 ID_HS = ped_B[i,1]
		 ID_HS_row1 = FS[which(HS$ID1==ID_HS), ]
		 ID_HS_row2 = FS[which(HS$ID2==ID_HS), ]
		 ID_HS_1 = unique(c((ID_HS_row1$ID2),(ID_HS_row2$ID1)))		  
		 if(length(ID_HS_1!=0)){
		 ped_B[i,7] = length(ID_HS_1)
		 ped_B[i,8] = paste(ID_HS_1, collapse = ",")
		 }else{
		 ped_B[i,8] = NA
		}
	 }
	 }
	 ped_B[ped_B[, 5] == 0, 2] = NA
	 ped_B[ped_B[, 7] == 0, 4] = NA	
	 colnames(ped_B) = colnames(ped_V)	 
	 ped_result = rbind(ped_V,ped_B)

	
	###########################################
	}else{
		 ID_H = unique(rownames(file))
		 ped_m = matrix(NA,length(ID_H),4)
		 colnames(ped_m)=c("ID","parents","parents","gender")
		 ped_m[,1] = ID_H
		 HSFS_ID = matrix(0,length(ped_m[,1]),4)
		 HS=PairL[c(which(PairL$TopRel=="HS")),]
		 FS=PairL[c(which(PairL$TopRel=="FS")),]
		 if (nrow(FS)!=0){
			 for(i in 1:nrow(ped_m)){
				 ID_FS = ped_m[i,1]
				 ID_FS_row1 = FS[which(FS$ID1==ID_FS), ]
				 ID_FS_row2 = FS[which(FS$ID2==ID_FS), ]
				 ID_FS_1 = unique(c((ID_FS_row1$ID2),(ID_FS_row2$ID1)))		  
					 if(length(ID_FS_1)!=0){
					 HSFS_ID[i,1] = length(ID_FS_1)
					 HSFS_ID[i,2] = paste(ID_FS_1, collapse = ",")
					 }else{
						 HSFS_ID[i,2] = NA
						 }
			}
		 }
		 #########################################################
		 if (nrow(HS)!=0){
			 for(i in 1:nrow(ped_m)){
			 ID_HS = ped_m[i,1]
			 ID_HS_row1 = FS[which(HS$ID1==ID_HS), ]
			 ID_HS_row2 = FS[which(HS$ID2==ID_HS), ]
			 ID_HS_1 = unique(c((ID_HS_row1$ID2),(ID_HS_row2$ID1)))		  
			 if(length(ID_HS_1!=0)){
			 HSFS_ID[i,3] = length(ID_HS_1)
			 HSFS_ID[i,4] = paste(ID_HS_1, collapse = ",")
			 }else{
			 HSFS_ID[i,4] = NA
			}
		 }
		 }
		 HSFS_ID[HSFS_ID[, 2] == 0, 2] = NA
		 HSFS_ID[HSFS_ID[, 4] == 0, 4] = NA
		 ped_m = cbind(ped_m,HSFS_ID)
		############################################################
		 ped_m=unique(ped_m)
		 familynames=matrix(NA, nrow =nrow(ped_m),ncol=1 )
		 colnames(familynames) = "familytrees"
		 ped_m = cbind(ped_m,familynames)
		 ped_result = ped_m
	 }
if(Tine = 1){
ped_TU = rbind(ped_result[,1:2],ped_result[,c(1,3)])
ped_TU = na.omit(ped_TU)
ped_TU = unique(ped_TU)
colnames(ped_TU) = c("ID","parents")
write.csv(ped_TU,paste(filename,"Tline.csv",collapse = "_"))
}
if(pedoutput = 1){

write.csv(ped_result,paste(filename,"pedigree.csv",collapse = "_"))
}
return(ped_result)
}	
	

GPEDI.pedOUT <- function(
                          data，
						  filenames=0
						  ){
file = data
target=file[,1]
target = na.omit(target)
coi=matrix(NA,length(target),3)
for(m in 1:length(target))
	{
	 target0=target[m]#要计算的目标个体
	 j=0#目标个体E的代数
	 E=target0
	 while(!all(is.na(E)))
		{
		 search.result=file[file[,1]%in%E,]
		 j=j+1 
		 E=c(search.result[3],search.result[4])#
		 E = unlist(E)
		 E = unique(E)
			}#如果E在文件file第一列并且不为NA，则代数加一，否则代数不变
		 coi[m,1]=target0 
		 coi[m,3]=j 
	}
	daishu = paste(min(coi[,3]),"-",max(coi[,3]))
	gender =table(file$gender)
	sire =as.numeric(gender[as.character("0")])
	if(is.na(sire)){sire = 0}
	dam =as.numeric(gender[as.character("1")])
	if(is.na(dam)){dam = 0}
	unknow = nrow(file)-sire-dam
	########################################################
	 p1 = file[,c(1,3)]
	 p2 = file[,c(1,4)]
	 colnames(p1) = colnames(p2)
	 PO = rbind(p1,p2)
	 PO = na.omit(PO)
	 colnames(PO) = c("ID1","ID2")
	 PO_1 = PO
	 familytree = 0
	 H_all = NULL
	 while(nrow(PO_1) != 0){
	 familytree = familytree+1
	 H = PO_1[familytree,]
	 H1 = H[[1]]
	 H2 = H[[2]]
	 H_c = H
	 #PO_1 = PO_1[!(PO_1$ID1==H1 & PO_1$ID2==H2),]
	 #nrow(PO_1)
	 #h = NULL
	 H_all[[familytree]] = H_c
	 #familytree
	 #H_all[[familytree]]
	 while(nrow(H)!= 0){
		 h = NULL
		 for(i in 1:nrow(H)){
		 H1 = H[i,1]
		 H2 = H[i,2]
		 #H1
		 #H2
		 PO_1 = PO_1[!(PO_1$ID1==H1 & PO_1$ID2==H2),]
		 PO_1 = PO_1[!(PO_1$ID2==H1 & PO_1$ID1==H2),]
		 H_3=PO_1[c(which(PO_1$ID1==H1)),]
		 H_4=PO_1[c(which(PO_1$HID2==H1)),]
		 H_5=PO_1[c(which(PO_1$ID1==H2)),]
		 H_6=PO_1[c(which(PO_1$ID2==H2)),]
		 h[[i]]= rbind(H_3,H_4,H_5,H_6)
		 h[[i]] = unique(h[[i]])
		 #nrow(PO_1)
		 #nrow(H)
		 					}
		 H= do.call(rbind, h)
		 H=unique(H)
		 #nrow(H)
		 #nrow(H_all[[1]])
		 #nrow(PO_1)
		 H_all[[familytree]]= rbind(H,H_all[[familytree]])
		 H_all[[familytree]]=unique(H_all[[familytree]])
		 #nrow(H_all[[1]])
		}
		PO_1 = PO_1 [!is.na(PO_1 $ID1),]
	}
	
	
	 H_all = Filter(function(x) !any(is.na(x)), H_all)
	 jiaxi = length(H_all)
	 jiaxinei = length(unique(c(PO$ID1,PO$ID2)))
	 zong = nrow(file)
	 jiaxiwai = zong-jiaxinei
	############################################################
	 freq_table = table(PO$ID2)
	 values_to_keep = names(freq_table[freq_table > 1])
	 filtered_df = PO[PO$ID2 %in% values_to_keep, ]
	 HS_FS = length(unique(filtered_df$ID1))
	 information=matrix(NA,4,4)
	 information[1,] = c("家系","家系内个体","家系外个体","总个体")
	 information[3,] = c("代数","公","母","未知性别")
	 information[2,] = c(jiaxi,jiaxinei,jiaxiwai,zong)
	 information[4,] = c(daishu,dam,sire,unknow)
	 information = as.data.frame(information)
	 if(filenames != 0){
	 write.xlsx(information, file = "output.xlsx")
	 write.xlsx(coi, file = "output_daishu.xlsx")
	 }
	 
	return(information)
				
				}	
GPEDI.inbreeding <-function(data#Kinship file
	)
	{
file = data
target3=file[,1]
target3 = na.omit(target3)
inb = NULL
coi=matrix(NA,length(target3),3)
####################################
for(m in 1:length(target3))
{
 target0=target3[m]#要计算的目标个体
 j=0#目标个体E的代数
 E=target0
 while(!all(is.na(E)))
	{
	 search.result=file[file[,1]%in%E,]
	 j=j+1 
	 E=c(search.result[,3],search.result[,4])#
	 E = unlist(E)
	 E = unique(E)
		}#如果E在文件file第一列并且不为NA，则代数加一，否则代数不变
	 coi[m,1]=target0 
	 coi[m,2]=j-1 
}
coi = coi[order(coi[, 2]), ]
coi[coi[, 2] == 1, 3] = 0
coi_1 = as.data.frame(coi)
colnames(coi_1) = c("ID","gener","inb")
MAX_gender = max(as.numeric(coi[,2]))
for(gener in 2:MAX_gender){
	 coi_3 = coi_1[c(which(coi_1$gener== gener)),]
	for(o in 1:nrow(coi_3)){
	 target0 = coi_3[o,1]
	 ######################################
	 tar.s=file[file[,1]%in%target0,3]#寻找目标个体的父本
	 tar.d=file[file[,1]%in%target0,4]#寻找目标个体的母本
	 cross.generation=FALSE
	 tar.all = c(tar.s,tar.d)
	 # build a matrix for pedigree
	 inbreed=matrix(NA,gener,2^(gener-1))#制定谱系矩阵
	 inbreed[1,1]=target0#将谱系矩阵第一排第一个指定为目标个体
	 # fill with ID by pedigree
		for(i in 1:(gener-1)){
		 search.store=inbreed[i,!is.na(inbreed[i,])]#指定search.store为下一代
		 search.num=length(search.store)#指定search.num为下一代的个数
			for(j in 1:search.num){
			 if(search.store[j]%in%file[,1])
				{	
				   search.result=file[file[,1]%in%search.store[j],]
				}else{
				   search.result=matrix(NA,1,4)
				}
				inbreed[i+1,2*j-1]=search.result[,3]
				inbreed[i+1,2*j]=search.result[,4]
			
			}#寻找父母本并放置到指定的位置
		}
		 ##########################################
		line=list()
		mm=0
		if(cross.generation)
		{
		 print("This function is pending ...")
		}else{# !cross.generation
		   ## search duplicated ID
		   for(i in 2:gener)
		   {
			  search.store=as.character(inbreed[i,!is.na(inbreed[i,])])
			  dup.id=unique(search.store[duplicated(search.store)])#将一代的个体查找重复
			  if(length(dup.id)>0)
			  {
				 for(dup in 1:length(dup.id))
				 {
					dup.index=grep(dup.id[dup],search.store)
					dup.num=length(dup.index)
					for(m in 1:dup.num)
					{
						mm=mm+1
						line[[mm]]=inbreed[i,dup.index[m]]
					}
					for(j in (i-1):2)
					{
					   dup.index.2=ceiling(dup.index/2)
					   for(m in 1:dup.num)
					   {
						line[[mm-m+1]]=append(line[[mm-m+1]],inbreed[j,dup.index.2[m]])
					   }
					   dup.index=dup.index.2
					}#将重复的个体的通路全部列出来
				 }
				
				 
			  }else{# !length(dup.id)>0
				 next
			  } # end length(dup.id)>0

		    }# end i in 2:gener
		 }# end cross.generation
	 #################################
	n=length(line)
	if(n!=0){
	 first.id=NULL
	 for(i in 1:n)
	 {
	  first.id=append(first.id,line[[i]][1])
	 }
	 last.id=NULL
	 for(i in 1:n)
	 {
	  last.id=append(last.id,line[[i]][length(line[[i]])])
	 }#分别将每条通路的第一个ID和最后一个ID组成两条list
	 table.matrix=as.matrix(table(last.id))
	 if(length(table.matrix)!=1)
	{
		as.matrix(table.matrix[order(table.matrix[,1]),])    
		min.id=rownames(table.matrix)[1]
		max.id=rownames(table.matrix)[2]
		min.num=as.numeric(table.matrix[min.id,])
		max.num=as.numeric(table.matrix[max.id,])#这四行通过目标个体父母本在通路中重复的次数多少分别名为min.id和max.id。然后给出重复的次数
		whole.line=list()
		line.inb = NULL
		mm=0
		for(i in 1:min.num)
		{
		  min.index=grep(min.id,last.id)[i]
		  min.first=line[[min.index]][1]
		  #for(j in 1:max.num)
		  #{
			last.id.2=last.id
			last.id.2[!first.id%in%min.first]=NA#指定前面找到通路的first.ID也就是重复个体（每次循环的重复个体不一样）的半条通路赛选出来
			last.id.2[min.index]=NA#将重复少的父本或者母本中的那半条通路剔除。得到重复多的父本或者母本的通路
			max.index=grep(max.id,last.id.2)#重复多的通路中通路所在的位置
			for(k in max.index)
			{
			  max.line=line[[k]]
			  max.line2=NULL
			  for(f in length(max.line):1)
			  {
				max.line2=append(max.line2,max.line[f])
			  }#这个循环是将重复多的通路中通路路倒过来
			  
			  the.line = append(max.line2,line[[min.index]][-1])
			  table.line=as.matrix(table(the.line))
			  Anumber = max(table.line[,1])
			  if(the.line[1]%in%tar.all & the.line[length(the.line)]%in%tar.all)
			  {	  
				 if(Anumber == 1)
				 {
					  mm=mm+1
					  whole.line[[mm]]=the.line
					  ancestors = min.index[1]
					  ancestors.inb = coi_1[c(which(coi_1$ID== ancestors)),3]
					  line.inb = c(line.inb,ancestors.inb)
				  }
			  } 
			  
			}
		}
			if(length(whole.line)!=0)
			{
				 x=0
				 i=1
				 for(i in (1:length(whole.line)))
				 {
					  T1=whole.line[[i]]
					  if(!is.na(line.inb[i])){
					  x1=((1/2)^length(T1))*(1+as.numeric(line.inb[i]))
					  }else{x1=((1/2)^length(T1))}
					  i=1+1
					  x=x+x1
				  }
			  
			}else{x=0}
	}else{x=0}#end length(table.matrix)!=1
	}else{x=0}#end n!=0
	 coi_1[which(coi_1$ID== target0),3] = x
	}#end nrow(coi_3)
}#end MAX_gender
  
return(coi_1)
}




























				