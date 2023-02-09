################################
#get discrete tip data in order
################################
tree<-read.tree(file="~/Desktop/ants_ecogeography/tree_w_nesting_nodes_for_ouwie.tre")
tree$node.label<-as.numeric(tree$node.label)
tree$node.label

spec<-read.csv(file="~/Desktop/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/ant_plant_diversification1_manuscript/manuscript_tables_etc/combined.updated.updated.pruned.wouts.badsremoved.badsremovedagain.nopheidrectispina_w_bammgen10nov2016.cleaned27april2017.csv",stringsAsFactors=FALSE)
rownames(spec)<-spec$Taxon
spec<-spec[spec$Family %in% "Formicidae",]
spec.mod<-spec[,c("Taxon","NestingModallupper")]
unique(spec.mod$NestingModallupper)
nrow(spec.mod)
#remove ambiguous taxa
spec.mod<-spec.mod[spec.mod$NestingModallupper %in% c(0,1),]
unique(spec.mod$NestingModallupper)
nrow(spec.mod)

#change to 1/2 instead of 0/1
spec.mod$NestingModallupper[spec.mod$NestingModallupper==1]<-2
spec.mod$NestingModallupper[spec.mod$NestingModallupper==0]<-1

#reduce spec.mod to those in tree
spec.mod<-spec.mod[spec.mod$Taxon %in% tree$tip.label,]
tree
nrow(spec.mod)

#reorder data to be in order of tree
spec.mod<-spec.mod[order(match(spec.mod$Taxon,tree$tip.label)),]
write.csv(spec.mod,file="~/Desktop/ants_ecogeography/discrete_tip_data_for_ouwie.csv",row.names=TRUE)
################################


################################
#get median raw trait value in order
################################
tree<-read.tree(file="~/Desktop/ants_ecogeography/tree_w_nesting_nodes_for_ouwie.tre")
dat.in.tre<-read.csv(file="~/Desktop/ants_ecogeography/ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.csv",stringsAsFactors=FALSE)
cols.to.get<-paste("bio_",1:19,"_orig_Median",sep="")
cols.to.get<-c(cols.to.get,"NPP_Median","PET_Median","ELEV_Median",paste(c("T_GRAVEL","T_SAND","T_SILT","T_CLAY","T_BULK_DENSITY","T_OC","T_PH_H2O","T_CEC_SOIL","T_BS","T_TEB","T_CACO3","T_CASO4","T_ESP","T_ECE"),"_Median",sep=""))
clean.dat<-dat.in.tre[cols.to.get]
rownames(clean.dat)<-dat.in.tre$Taxon
clean.dat<-clean.dat[rownames(clean.dat) %in% tree$tip.label,]
colnames(clean.dat)<-gsub("_orig_Median","",colnames(clean.dat))
colnames(clean.dat)<-gsub("_Median","",colnames(clean.dat))
clean.dat<-clean.dat[match(tree$tip.label,rownames(clean.dat)),]
nrow(clean.dat)
head(clean.dat)
write.csv(clean.dat,file="~/Desktop/ants_ecogeography/ind_var_median_values_for_ouwie.csv",row.names=TRUE)

