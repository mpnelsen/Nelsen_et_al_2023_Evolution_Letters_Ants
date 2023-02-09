require("OUwie")
#filepath<-"/home/mpnelsen/ant_ecogeo_ouwie_analyses/"
filepath<-"/home/mpnelsen/ant_ecogeo_ouwie_analyses/analyses_for_fewer_variables/temp_in_C/"

traits<-c(paste("bio_",1:19,sep=""))
traits<-c(traits,"NPP","PET","ELEV")

#vector of codings to explore
codings<-"NestingModallupper"

#vector of ouwie trait evolution models to fit
models<-c("BM1","BMS","OU1","OUM","OUMV","OUMA","OUMVA")

#bio1 OUMVA
load(paste(filepath,codings,"_",traits[1],"_",models[7],".Rsave",sep=""))

for(s in 1:10){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}

for(s in 11:20){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}

for(s in 21:30){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}

for(s in 31:40){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}

for(s in 41:50){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}


for(s in 51:60){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}

for(s in 61:70){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}
for(s in 71:80){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}
for(s in 81:90){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}
for(s in 91:100){
	df.sim<-OUwie.sim(phy=tree,data=df.clean,simmap.tree=FALSE,alpha=fit$solution[1,],sigma.sq=fit$solution[2,],theta=fit$theta[,1],theta0=fit$theta[1,1],mserr="none")
	sim.fit<-OUwie(tree,df.sim,model=models[7],simmap.tree=FALSE,root.station=TRUE,diagn=TRUE)
	#save the tree used, the data used, and the R results
	save(tree,df.sim,sim.fit,file=paste(filepath,codings,"_",traits[1],"_",models[7],"_","sim_",s,".Rsave",sep=""))
}