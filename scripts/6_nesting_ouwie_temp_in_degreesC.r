require(phytools)
require(phangorn)
require(diversitree)
require(circlize)

########################################
########################################
####PRE-STUFF FINISHED - NOW ON To OUwie
########################################
#this is for temp vars/10 so it is now in degrees C instead of degrees C*10

require("OUwie")
#filepath<-"/home/mpnelsen/ant_ecogeo_ouwie_analyses/"
filepath<-"/home/mpnelsen/ant_ecogeo_ouwie_analyses/analyses_for_fewer_variables/temp_in_C/"

#broke some of these up to...traits[1:3,4:8,9:11]
#read continuous data and restrict to env
traits<-c(paste("bio_",1:11,sep=""))
conts<-read.csv(file=paste(filepath,"ind_var_median_values_for_ouwie.csv",sep=""),stringsAsFactors=FALSE,row.names=1)
colnames(conts)<-traits
traits<-c(paste("bio_",1:11,sep=""))
#divide by 10
conts<-conts[,traits]/10

#vector of codings to explore
codings<-"NestingModallupper"

#read discrete data
disc<-read.csv(file=paste(filepath,"discrete_tip_data_for_ouwie.csv",sep=""),stringsAsFactors=FALSE,row.names=1)

#read tree
tree<-read.tree(file=paste(filepath,"tree_w_nesting_nodes_for_ouwie.tre",sep=""))
tree$node.label<-as.numeric(tree$node.label)
tree$node.label

#vector of ouwie trait evolution models to fit
models<-c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")

#vector of column names for results
column.names<-c("CodingScheme","Trait","Model","lnL","AIC","AICc","ntax","alpha1","alpha2","alpha3","sigmassq1","sigmassq2","sigmassq3","optimum1","optimum2","optimum3","se1","se2","se3","eigval1","eigval2","eigval3","eigval4","eigval5","eigval6","EigvalIssue")

#make empty data frame where results of ouwie analyses will go
df.out<-data.frame(matrix(nrow=length(codings)*length(traits)*length(models),ncol=length(column.names)))

#assign column names
colnames(df.out)<-column.names

#iter.no will be used to assign results to a row number in df.out.  after each ouwie analysis, 1 is added (see below)
iter.no<-0

#for each coding...
for(c in 1:length(codings)){
	#for each trait...
	for(t in 9:11){
		#add continuous data
		df.clean<-disc
		df.clean[,traits[t]]<-NA
		for(xz in 1:nrow(conts)){
			df.clean[df.clean$Taxon %in% rownames(conts)[xz],traits[t]]<-conts[xz,traits[t]]
		}
		df.clean<-df.clean[order(match(rownames(df.clean),tree$tip.label)),]
		#fit each model
		for(m in 1:length(models)){
			#add 1 to iteration number (starts at 0), and this will be the row in the results data frame that we add info to
			iter.no<-iter.no+1
			#for bms, ouwie documentation recommends have root.station=FALSE
			if(models[m] %in% "BMS"){
				fit<-OUwie(tree,df.clean,model=models[m],simmap.tree=FALSE,root.station=FALSE,diagn=TRUE)		
			}
			if(!models[m] %in% "BMS"){
				fit<-OUwie(tree,df.clean,model=models[m],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)		
			}
			#save the tree used, the data used, and the R results
			save(tree,df.clean,fit,file=paste(filepath,codings[c],"_",traits[t],"_",models[m],".Rsave",sep=""))
			#add relevant info from results to a data frame
			df.out[iter.no,column.names[1]]<-codings[c]
			df.out[iter.no,column.names[2]]<-traits[t]
			df.out[iter.no,column.names[3]]<-models[m]
			df.out[iter.no,column.names[4]]<-fit$loglik
			df.out[iter.no,column.names[5]]<-fit$AIC
			df.out[iter.no,column.names[6]]<-fit$AICc
			df.out[iter.no,column.names[7]]<-nrow(fit$res)
			df.out[iter.no,column.names[8]]<-fit$solution[1,1]
			df.out[iter.no,column.names[9]]<-fit$solution[1,2]
			df.out[iter.no,column.names[11]]<-fit$solution[2,1]
			df.out[iter.no,column.names[12]]<-fit$solution[2,2]
			df.out[iter.no,column.names[14]]<-fit$theta[1,1]
			df.out[iter.no,column.names[15]]<-fit$theta[2,1]
			df.out[iter.no,column.names[17]]<-fit$theta[1,2]
			df.out[iter.no,column.names[18]]<-fit$theta[2,2]
			#if 3 states (instead of 2), add info for third state (this equals 4 in 2 state situations (2 alpha and 2 sigma sq))
			if(length(fit$solution)==6){
				df.out[iter.no,column.names[10]]<-fit$solution[1,3]	
				df.out[iter.no,column.names[13]]<-fit$solution[2,3]	
				df.out[iter.no,column.names[16]]<-fit$theta[3,1]
				df.out[iter.no,column.names[19]]<-fit$theta[3,2]
			}
			df.out[iter.no,column.names[20:25]]<-fit$eigval[1:6]
			pos<-sum(fit$eigval>0)
			neg<-sum(fit$eigval<0)
			if(sum(c(pos,neg)>0)==2){
				df.out[iter.no,column.names[26]]<-"Yes"
			}
			if(sum(c(pos,neg)>0)!=2){
				df.out[iter.no,column.names[26]]<-"No"
			}
			#save data frame after each model fit
			write.csv(df.out,file=paste(filepath,"summary_ouwie_analyses_median_temp_in_C_9to11.csv",sep=""),row.names=FALSE)
		}
	}
}

