
#note that the habcod column read in is coded incorrectly (as 0.33, not 33.3), but that does not matter for this analysis here since that column is not used at all.
sum.table<-read.csv(file="~/Documents/papers_reviews/papers/ants_ecogeography/ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.temps_in_C.csv",stringsAsFactors=FALSE)


open<-c("GPcnt","HPcnt","IPcnt","JPcnt","KPcnt","MPcnt","LPcnt")
closed<-c("APcnt","BPcnt","CPcnt","DPcnt","EPcnt","NPcnt","FPcnt")

sum.table$OpenMedOpen<-0
sum.table$ClosedMedOpen<-0

for(x in 1:nrow(sum.table)){
	if(sum(sum.table[x,open])>33.3){
		sum.table$OpenMedOpen[x]<-1
	}
	if(sum(sum.table[x,closed])>33.3){
		sum.table$ClosedMedOpen[x]<-1
	}
}

sum.table$habcodMedOpen<-NA

for(x in 1:nrow(sum.table)){
	if(sum.table$OpenMedOpen[x]==0 & sum.table$ClosedMedOpen[x]==1){
		sum.table$habcodMedOpen[x]<-0
	}
	if(sum.table$OpenMedOpen[x]==1 & sum.table$ClosedMedOpen[x]==0){
		sum.table$habcodMedOpen[x]<-1
	}
	if(sum.table$OpenMedOpen[x]==1 & sum.table$ClosedMedOpen[x]==1){
		sum.table$habcodMedOpen[x]<-"0&1"
	}
}

write.csv(sum.table,file="~/Documents/papers_reviews/papers/ants_ecogeography/ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.temps_in_C_mediterropen_33.3.csv",row.names=FALSE)


#ASR of Open/Closed
require(corHMM)
new.output<-read.csv(file="~/Documents/papers_reviews/papers/ants_ecogeography/ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.temps_in_C_mediterropen_33.3.csv",stringsAsFactors=FALSE)
rownames(new.output)<-new.output$Taxon
tree<-read.tree(file="~/Documents/papers_reviews/papers/ants_ecogeography/hab.data.plates.min.thresh.29may2018.wsoils.nona.tre")
new.output <-new.output[match(tree$tip.label,rownames(new.output)),]
new.output$habcodMedOpenthree<-NA
for(x in 1:nrow(new.output)){
	if(new.output$OpenMedOpen[x]==0 & new.output$ClosedMedOpen[x]==1){
		new.output$habcodMedOpenthree[x]<-0
	}
	if(new.output$OpenMedOpen[x]==1 & new.output$ClosedMedOpen[x]==0){
		new.output$habcodMedOpenthree[x]<-2
	}
	if(new.output$OpenMedOpen[x]==1 & new.output$ClosedMedOpen[x]==1){
		new.output$habcodMedOpenthree[x]<-1
	}
}
cleaned.tri<-new.output[!is.na(new.output$habcodMedOpenthree),c("Taxon","habcodMedOpenthree")]
my.mat<-rate.mat.maker(rate.cat=1,hrm=FALSE,ntraits=1,nstates=3,model="ARD")
my.mat.test<-rayDISC(phy=tree,data=cleaned.tri,rate.mat=my.mat,model="ARD",node.states="marginal",root.p="maddfitz")
save(my.mat.test,file="~/Documents/papers_reviews/papers/ants_ecogeography/my.mat.test.raydisc.MEDOPEN_33.3.Rsave")
#plotRECON(tree, my.mat.test$states,show.tip.label=FALSE,edge.width=0.2,pie.cex=0.75,piecolors=c("dark green", "light green", "yellow"))

my.mat.con<-rate.par.drop(rate.mat.index=my.mat,drop.par=c(2,5))
my.mat.test.con<-rayDISC(phy=tree,data=cleaned.tri,rate.mat=my.mat.con,model="ARD",node.states="marginal",root.p="maddfitz")
save(my.mat.test.con,file="~/Documents/papers_reviews/papers/ants_ecogeography/my.mat.test.con.raydisc.MEDOPEN_33.3.Rsave")
#plotRECON(tree, my.mat.test.con$states,show.tip.label=FALSE,edge.width=0.2,pie.cex=0.75,piecolors=c("dark green", "light green", "yellow"))

cleaned<-new.output[!is.na(new.output$habcodMedOpen),c("Taxon","habcodMedOpen")]
my.test<-rayDISC(phy=tree,data=cleaned,model="ARD",node.states="marginal",root.p="maddfitz")
save(my.test,file="~/Documents/papers_reviews/papers/ants_ecogeography/my.test.raydisc.MEDOPEN_33.3.Rsave")
#plotRECON(new.tree,my.test$states,show.tip.label=FALSE,edge.width=0.2,pie.cex=0.75,piecolors=c("dark green", "yellow"))





#Plot biome type at tips and pie-charts for open/closed
require(BAMMtools)
library(circlize)
require(RColorBrewer)
require(diversitree)
require(phytools)
mcmcout <- read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ant_plant_interactions/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/Hymenoptera/drop_bads_again/drop_bads_again_partfinder/proper_raxml_ml500searches_Pristo_Out/formicidae_ml_treepl_147/bamm_analyses_185_10nov2016_w_correct_priors_background1/bamm_analyses_185_10nov2016_w_correct_priors_background1_medusa/mcmc_out_medusa.txt", header=T)

burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#library(coda)
#effectiveSize(postburn$N_shifts)
#effectiveSize(postburn$logLik)

#tree<-read.tree("~/Desktop/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/Hymenoptera/drop_bads_again/drop_bads_again_partfinder/proper_raxml_ml500searches_Pristo_Out/formicidae_ml_treepl_147/bamm_analyses_185_10nov2016_w_correct_priors_background1/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")

#load my.mat.test.con
load("~/Documents/papers_reviews/papers/ants_ecogeography/my.mat.test.con.raydisc.MEDOPEN_33.3.Rsave")
seqs<-read.csv(file="/Users/matthewnelsen/Documents/papers_reviews/papers/ant_plant_interactions/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/Hymenoptera/drop_bads_again/drop_bads_again_partfinder/proper_raxml_ml500searches_Pristo_Out/formicidae_ml_treepl_147/hisse_analyses_ants_28april2017/combined.updated.updated.pruned.wouts.badsremoved.badsremovedagain.nopheidrectispina_w_bammgen10nov2016.cleaned27april2017.csv",stringsAsFactors=FALSE)
in.tree<-seqs[seqs$Taxon %in% my.mat.test.con$phy$tip.label,]
in.tree<-in.tree[order(match(in.tree$Taxon,my.mat.test.con$phy$tip.label)),]
rownames(in.tree)<-in.tree$Taxon


dat.in.tre<-read.csv(file="~/Documents/papers_reviews/papers/ants_ecogeography/ant_summaries_w_latlong_merged_30may2018_onlywsoilsandenv.CLEANED.temps_in_C_mediterropen_33.3.csv",stringsAsFactors=FALSE)
dat.in.tre<-dat.in.tre[order(match(dat.in.tre$Taxon,my.mat.test.con$phy$tip.label)),]
rownames(dat.in.tre)<-dat.in.tre$Taxon




bcs<-cbind(dat.in.tre,in.tree)
bcs<-bcs[order(match(bcs$Taxon,my.mat.test.con$phy$tip.label)),]

tree<-read.tree("/Users/matthewnelsen/Documents/papers_reviews/papers/ant_plant_interactions/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/Hymenoptera/drop_bads_again/drop_bads_again_partfinder/proper_raxml_ml500searches_Pristo_Out/formicidae_ml_treepl_147/bamm_analyses_185_10nov2016_w_correct_priors_background1/ML_TREE_treepl_185_outsdroppeds.ladderized.increasing.tre")
edata <- getEventData(tree, eventdata = "/Users/matthewnelsen/Documents/papers_reviews/papers/ant_plant_interactions/14july2016_antcat_aligns_fullalign_einsi_nogappy_wadds_heavilymodified_correctedmacse/Hymenoptera/drop_bads_again/drop_bads_again_partfinder/proper_raxml_ml500searches_Pristo_Out/formicidae_ml_treepl_147/bamm_analyses_185_10nov2016_w_correct_priors_background1/bamm_analyses_185_10nov2016_w_correct_priors_background1_medusa/event_data_medusa.txt", burnin=0.1,nsamples=2000)

edata.sub<-subtreeBAMM(edata,tips=my.mat.test.con$phy$tip.label)
edata<-edata.sub

shift_probs <- summary(edata)
postfile<-mcmcout
bfmat <- computeBayesFactors(postfile, expectedNumberOfShifts=50, burnin=0.1)

#plots single best shift configuration
msc.set <- maximumShiftCredibility(edata, maximize='product')
msc.config <- subsetEventData(edata, index = msc.set$sampleindex)


#plot any occurrences in a biome, or shade according to abundance?  This is because biomes data includes any present.  need to use code for plotting position on ppca at tip and ramp from white (0) to color...very few if any will have 100% color as likely spread across multiple biomes


#cats<-c("APcnt","BPcnt","CPcnt","DPcnt","EPcnt","FPcnt","GPcnt","HPcnt","IPcnt","JPcnt","KPcnt","LPcnt","MPcnt","NPcnt")
cats<-c("APcnt","BPcnt","CPcnt","DPcnt","EPcnt","FPcnt","NPcnt","GPcnt","HPcnt","IPcnt","JPcnt","KPcnt","MPcnt","LPcnt")

bcs[,cats]<-round(bcs[,cats])
#mycols<-c("forestgreen","darkolivegreen","palegreen","darkgreen","blue","darkseagreen","orange","lemonchiffon","lightblue","palegoldenrod","cyan","red","salmon","purple")
mycols<-c("forestgreen","darkolivegreen","palegreen","darkgreen","blue","darkseagreen","purple","orange","lemonchiffon","lightblue","palegoldenrod","cyan","salmon","red")


for(x in 1:length(cats)){
	colfunc <- colorRampPalette(c("white",mycols[x]))
	colstouse<-colfunc(101)
	for(p in 1:nrow(bcs)){
		if(is.na(bcs[p,cats[x]])){
			bcs[p,cats[x]]<-"#FFFFFF"
		}
		if(!bcs[p,cats[x]] %in% "#FFFFFF"){
			bcs[p,cats[x]]<-colstouse[as.numeric(bcs[p,cats[x]])+1]
		}
	}
}


p2p<-diversitree:::plot2.phylo
t <- max(branching.times(my.mat.test.con$phy))
obj<-p2p(my.mat.test.con$phy,type="fan",label.offset = t * 1/4, show.tip.label = FALSE)
vz<-pb(msc.config,breaksmethod="jenks",pal="terrain",spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.75,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE)
vz$coords[,1][row.names(vz$coords)==nrow(bcs)]
max(vz$coords[,4])

ht<-max(nodeHeights(my.mat.test.con$phy))
#since bamm does not put root right at center, need to make height greater than root height
ht<-ht*max(vz$coords[,4])

##############
##############
##############
png("~/Documents/papers_reviews/papers/ants_ecogeography/biomes_canopy_type_medopen_33.3_norates.png",width=11,height=11,units="in",res=600)
circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1.524,1.524),canvas.ylim=c(-1.524,1.524))
circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 


timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
require(phytools)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<ht,]

for(x in 1:nrow(timey)){
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-timey$End[x])/ht,rou2=(ht-timey$Start[x])/ht,col=timey$RGB[x],border=NA)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-50)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-100)/ht,col=NA,border="gray77",lwd=0.75,lty=3)
}
circos.clear()
par(new = T)


bammplot<-pb(msc.config,breaksmethod="jenks",pal="black",spex="s",method="polar",vtheta=0,xlim=c(-5,5),lwd=0.5,ofs=0.627,mar=c(0.55,0.55,0.55,0.55),par.reset=FALSE,legend=FALSE)

xy<-obj$xy
theta<-xy$theta[seq_along(my.mat.test.con$phy$tip.label)]
dt<-diff(sort(theta))[1]/2

comb.mat<-bcs
edata$tip.label->tiplabs
rownames(bcs)<-bcs$Taxon
tax2<-bcs[match(tiplabs,rownames(bcs)),]
rownames(comb.mat)<-rownames(bcs)
dat<-comb.mat[match(tiplabs,rownames(comb.mat)),]
dat<-dat[,cats]

div.fa<-diversitree:::filled.arcs
w = 1/50
for (i in seq_along(cats)){
	div.fa(theta-dt,theta+dt,rep(1.14+(i*0.02),length(theta)),w,dat[i])
}

nodelabels(pie=my.mat.test.con$states[,1:3],cex=0.2,piecol=c("dark green", "light green", "yellow"))

keepers<-names(table(bcs$Subfamily))[table(bcs$Subfamily)>20]
keepers<-keepers[!keepers %in% c("Amblyoponinae", "Myrmeciinae")]
bcs$ClassPlot<-NA
for(x in 1:nrow(bcs)){
	if(bcs$Subfamily[x] %in% keepers){
		bcs$ClassPlot[x]<-bcs$Subfamily[x]
	}
}
div.arcs<-diversitree:::arcs
n.taxa <- obj$n.taxa
col.bar<-"black"
n <- obj$n.spp
offset.bar=w * (n + 3)
offset.lab=w * (n + 4)
dy <-1/6
dy <- dy/n * 2 * pi
yy <- obj$xy$theta[seq_len(obj$Ntip)]
y0 <- tapply(yy - dy, bcs$ClassPlot, min)
y1 <- tapply(yy + dy, bcs$ClassPlot, max)
x.bar <- rep(max(obj$xx) + offset.bar, length(y0))
x.lab <- rep(max(obj$xx) + offset.lab, length(y0))
div.arcs(y0, y1, x.bar-max(x.bar)+1.47, col = col.bar, lwd = 1)

classification<-bcs$ClassPlot[order(match(bcs$Taxon, my.mat.test.con$phy$tip.label))]
clades.to.plot<-unique(classification[!is.na(classification)])
nodes.for.clades<-NULL
for(x in 1:length(clades.to.plot)){
	nodes.for.clades <-c(nodes.for.clades,findMRCA(my.mat.test.con$phy,tips=bcs$Taxon[bcs$Subfamily %in% clades.to.plot[x]]))
}

clades.to.plot[clades.to.plot %in% "Myrmicinae"]<-"Myrmicinae                                          "
for(i in 1:length(nodes.for.clades)){
	arc.cladelabels.mod(text=clades.to.plot[i],node=nodes.for.clades[i],dr.arc=FALSE,lab.offset=1.34,cex=1.5)
}

#legend(x=1,y=-1,title="Tip States (Biomes)",c("Tropical and Subtropical Moist Broadleaf Forests","Tropical and Subtropical Dry Broadleaf Forests","Tropical and Subtropical Coniferous Forests","Temperate Broadleaf and Mixed Forests","Temperate Coniferous Forests","Boreal Forests and Taiga","Mediterranean Forests, Woodlands and Scrubs","Mangrove","Tropical and Subtropical Grasslands, Savannas and Shrublands","Temperate Grasslands, Savannas and Shrublands","Flooded Grasslands and Savannas","Montane Grasslands and Shrublands","Tundra","Deserts and Xeric Shrublands"),col=mycols,pch=19,bty="n",cex=0.8,pt.cex=1)

text("Tip States (Biomes)",font=2,x=-1.48,y=-1.44,cex=1.5)

text("Closed Canopy",font=3,x=-1.68,y=-1.5,cex=1)
segments(x0=-1.9,y0=-1.52,x1=-0.75,y1=-1.52,lwd=1.5)
legend(x=-1.9,y=-1.5,c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Coniferous Forests","Temperate Broadleaf & Mixed Forests","Temperate Coniferous Forests","Boreal Forests & Taiga","Mangrove"),col=mycols[1:7],pch=19,bty="n",cex=0.75,pt.cex=1.5)

text("Open Canopy",font=3,x=-0.55,y=-1.50,cex=1)
segments(x0=-0.75,y0=-1.52,x1=0.5,y1=-1.52,lwd=1.5)
legend(x=-0.75,y=-1.50,c("Tropical & Subtropical Grasslands, Savannas & Shrublands","Temperate Grasslands, Savannas & Shrublands","Flooded Grasslands & Savannas","Montane Grasslands & Shrublands","Tundra","Deserts & Xeric Shrublands","Mediterranean Forests, Woodlands & Scrubs"),col=mycols[8:14],pch=19,bty="n",cex=0.75,pt.cex=1.5)


text("Node States (Canopy Type)",font=2,x=1.325,y=-1.44,cex=1.5)
segments(x0=0.75,y0=-1.52,x1=1.9,y1=-1.52,lwd=1.5)
legend(x=0.75,y=-1.5,c("Open Only", "Open & Closed", "Closed Only"),col=c("yellow", "light green", "dark green"),pch=19,bty="n",cex=0.75,pt.cex=1.5)
#addBAMMlegend(bammplot,direction="horizontal",location=c(0.75,1.85,-1.6,-1.55),cex.axis=1.2)
#text(x=1.28,y=-1.35,expression(bold("Speciation Rate")),cex=1.5)

dev.off()







#from http://blog.phytools.org/2017/03/clade-labels-on-circular-fan-tree.html
arc.cladelabels.mod<-function(tree=NULL,text,node,ln.offset=1.02,lab.offset=1.25,cex=1,dr.arc=TRUE,orientation="curved",...){
    require(plotrix)
    obj<-get("last_plot.phylo",envir=.PlotPhyloEnv)
    if(obj$type!="fan") stop("method works only for type=\"fan\"")
    h<-max(sqrt(obj$xx^2+obj$yy^2))
    if(hasArg(mark.node)) mark.node<-list(...)$mark.node
    else mark.node<-TRUE
    #if(mark.node) points(obj$xx[node],obj$yy[node],pch=21,bg="red")
    if(is.null(tree)){
        tree<-list(edge=obj$edge,tip.label=1:obj$Ntip,
            Nnode=obj$Nnode)
        class(tree)<-"phylo"
    }
    d<-getDescendants(tree,node)
    d<-sort(d[d<=Ntip(tree)])
    deg<-atan(obj$yy[d]/obj$xx[d])*180/pi
    ii<-intersect(which(obj$yy[d]>=0),which(obj$xx[d]<0))
    deg[ii]<-180+deg[ii]
    ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]<0))
    deg[ii]<-180+deg[ii]
    ii<-intersect(which(obj$yy[d]<0),which(obj$xx[d]>=0))
    deg[ii]<-360+deg[ii]
    if(isTRUE(dr.arc)){
    		draw.arc(x=0,y=0,radius=ln.offset*h,deg1=min(deg),deg2=max(deg))
    }
    if(orientation=="curved")
        arctext(text,radius=lab.offset*h,
            middle=mean(range(deg*pi/180)),cex=cex,font=3)
    else if(orientation=="horizontal"){
        x0<-lab.offset*cos(median(deg)*pi/180)*h
        y0<-lab.offset*sin(median(deg)*pi/180)*h
        text(x=x0,y=y0,label=text,
        adj=c(if(x0>=0) 0 else 1,if(y0>=0) 0 else 1),
        offset=0)
    }
}


#took out ofs=0 and moved to options in plot.bammdata in BAMMtools 2.1.4
pb<-function (x, tau = 0.01, method = "phylogram", ofs=0, xlim = NULL, ylim = NULL, 
    vtheta = 5, rbf = 0.001, show = TRUE, labels = FALSE, legend = FALSE, 
    spex = "s", lwd = 1, cex = 1, pal = "RdYlBu", mask = integer(0), 
    mask.color = gray(0.5), colorbreaks = NULL, logcolor = FALSE, 
    breaksmethod = "linear", color.interval = NULL, JenksSubset = 20000, 
    par.reset = FALSE, direction = "rightwards", ...) 
{
    as.phylo.bammdata<-BAMMtools:::as.phylo.bammdata
    colorMap<-BAMMtools:::colorMap
    setPolarTreeCoords<-BAMMtools:::setPolarTreeCoords
    mkdtsegsPolar<-BAMMtools::: mkdtsegsPolar
    div.arcs<-diversitree:::arcs
    arc<-BAMMtools:::arc
	if ("bammdata" %in% class(x)) {
        if (attributes(x)$order != "cladewise") {
            stop("Function requires tree in 'cladewise' order")
        }
        phy <- as.phylo.bammdata(x)
    }
    else stop("Object ephy must be of class bammdata")
    if (!spex %in% c("s", "e", "netdiv")) {
        stop("spex must be 's', 'e' or 'netdiv'.")
    }
    if (length(pal) == 1 && !pal %in% names(get("palettes", envir = .colorEnv)) && 
        pal != "temperature" && pal != "terrain") 
        pal <- rep(pal, 3)
    else if (length(pal) == 2) 
        pal <- c(pal, pal[2])
    if (breaksmethod == "linear" & !is.null(color.interval)) {
        if (length(color.interval) != 2) {
            stop("color.interval must be a vector of 2 numeric values.")
        }
    }
    if (!is.binary.tree(phy)) {
        stop("Function requires fully bifurcating tree")
    }
    if (any(phy$edge.length == 0)) {
        warning("Tree contains zero length branches. Rates for these will be NA and coerced to zero")
    }
    if (!("dtrates" %in% names(x))) {
        x <- dtRates(x, tau)
    }
    if (is.null(colorbreaks)) {
        colorbreaks <- assignColorBreaks(x$dtrates$rates, 64, 
            spex, logcolor, breaksmethod, JenksSubset)
    }
    if (x$type == "trait") {
        colorobj <- colorMap(x$dtrates$rates, pal, colorbreaks, 
            logcolor, color.interval)
    }
    else if (x$type == "diversification") {
        if (tolower(spex) == "s") {
            colorobj <- colorMap(x$dtrates$rates[[1]], pal, colorbreaks, 
                logcolor, color.interval)
        }
        else if (tolower(spex) == "e") {
            colorobj <- colorMap(x$dtrates$rates[[2]], pal, colorbreaks, 
                logcolor, color.interval)
        }
        else if (tolower(spex) == "netdiv") {
            colorobj <- colorMap(x$dtrates$rates[[1]] - x$dtrates$rates[[2]], 
                pal, colorbreaks, logcolor, color.interval)
        }
    }
    else {
        stop("Unrecognized/corrupt bammdata class. Type does not equal 'trait' or 'diversification'")
    }
    edge.color <- colorobj$cols
    tH <- max(x$end)
    phy$begin <- x$begin
    phy$end <- x$end
    tau <- x$dtrates$tau
    if (method == "polar") {
        ret <- setPolarTreeCoords(phy, vtheta, rbf)
        rb <- tH * rbf
        p <- mkdtsegsPolar(ret$segs[-1, ], tau, x$edge)
    }
    else if (method == "phylogram") {
        ret <- setPhyloTreeCoords(phy)
        p <- mkdtsegsPhylo(ret$segs[-1, ], tau, x$edge)
    }
    else {
        stop("Unimplemented method")
    }
    x0 <- c(ret$segs[1, 1], p[, 1])
    x1 <- c(ret$segs[1, 3], p[, 2])
    y0 <- c(ret$segs[1, 2], p[, 3])
    y1 <- c(ret$segs[1, 4], p[, 4])
    offset <- table(p[, 5])[as.character(unique(p[, 5]))]
    if (length(mask)) {
        edge.color[p[, 5] %in% mask] <- mask.color
    }
    arc.color <- c(edge.color[1], edge.color[match(unique(p[, 
        5]), p[, 5]) + offset])
    edge.color <- c(edge.color[1], edge.color)
    if (show) {
        op <- par(no.readonly = TRUE)
        if (length(list(...))) {
            par(...)
        }
        if (legend) {
            par(mar = c(5, 4, 4, 5))
        }
        plot.new()
        #ofs <- 0     #####THIS
        if (labels) {
            if (method == "phylogram") 
                ofs <- max(nchar(phy$tip.label) * 0.03 * cex * 
                  tH)
            else ofs <- max(nchar(phy$tip.label) * 0.03 * cex)
        }
        if (method == "polar") {
            plot.window(xlim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                ofs), ylim = c(-1, 1) + c(-rb, rb) + c(-ofs, 
                ofs), asp = 1)
            segments(x0, y0, x1, y1, col = edge.color, lwd = lwd, 
                lend = 2)
            arc(0, 0, ret$arcs[, 1], ret$arcs[, 2], c(rb, rb + 
                phy$end/tH), border = arc.color, lwd = lwd)
            if (labels) {
                for (k in 1:length(phy$tip.label)) {
                  text(ret$segs[-1, ][phy$edge[, 2] == k, 3], 
                    ret$segs[-1, ][phy$edge[, 2] == k, 4], phy$tip.label[k], 
                    cex = cex, srt = (180/pi) * ret$arcs[-1, 
                      ][phy$edge[, 2] == k, 1], adj = c(0, NA))
                }
            }
        }
        if (method == "phylogram") {
            direction <- match.arg(direction, c("rightwards", 
                "leftwards", "downwards", "upwards"))
            if (direction == "rightwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), 0)
                arcs <- redirect(ret$arcs, 0)
                bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
                arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
                ret$segs[-1, c(1, 3)] <- tH * ret$segs[-1, c(1, 
                  3)]
            }
            else if (direction == "leftwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), pi)
                bars[, c(2, 4)] <- abs(bars[, c(2, 4)])
                arcs <- redirect(ret$arcs, pi)
                arcs[, c(2, 4)] <- abs(arcs[, c(2, 4)])
                bars[, c(1, 3)] <- tH * bars[, c(1, 3)]
                arcs[, c(1, 3)] <- tH * arcs[, c(1, 3)]
                ret$segs[-1, c(1, 3)] <- -tH * ret$segs[-1, c(1, 
                  3)]
            }
            else if (direction == "downwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), -pi/2)
                arcs <- redirect(ret$arcs, -pi/2)
                bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
                arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
                ret$segs <- redirect(ret$segs, -pi/2)
                ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2, 4)]
            }
            else if (direction == "upwards") {
                bars <- redirect(cbind(x0, y0, x1, y1), pi/2)
                bars[, c(1, 3)] <- abs(bars[, c(1, 3)])
                arcs <- redirect(ret$arcs, pi/2)
                arcs[, c(1, 3)] <- abs(arcs[, c(1, 3)])
                bars[, c(2, 4)] <- tH * bars[, c(2, 4)]
                arcs[, c(2, 4)] <- tH * arcs[, c(2, 4)]
                ret$segs <- redirect(ret$segs, pi/2)
                ret$segs[, c(1, 3)] <- abs(ret$segs[, c(1, 3)])
                ret$segs[, c(2, 4)] <- tH * ret$segs[, c(2, 4)]
            }
            if (is.null(xlim) && direction == "rightwards") 
                xlim <- c(0, tH + ofs)
            if (is.null(xlim) && direction == "leftwards") 
                xlim <- c(-(tH + ofs), 0)
            if (is.null(ylim) && (direction == "rightwards" || 
                direction == "leftwards")) 
                ylim <- c(0, phy$Nnode)
            if (is.null(xlim) && (direction == "upwards" || direction == 
                "downwards")) 
                xlim <- c(0, phy$Nnode)
            if (is.null(ylim) && direction == "upwards") 
                ylim <- c(0, tH + ofs)
            if (is.null(ylim) && direction == "downwards") 
                ylim <- c(-(tH + ofs), 0)
            plot.window(xlim = xlim, ylim = ylim)
            segments(bars[-1, 1], bars[-1, 2], bars[-1, 3], bars[-1, 
                4], col = edge.color[-1], lwd = lwd, lend = 2)
            isTip <- phy$edge[, 2] <= phy$Nnode + 1
            isTip <- c(FALSE, isTip)
            segments(arcs[!isTip, 1], arcs[!isTip, 2], arcs[!isTip, 
                3], arcs[!isTip, 4], col = arc.color[!isTip], 
                lwd = lwd, lend = 2)
            if (labels) {
                if (direction == "rightwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 4, offset = 0.25)
                else if (direction == "leftwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 2, offset = 0.25)
                else if (direction == "upwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 4, srt = 90, offset = 0)
                else if (direction == "downwards") 
                  text(ret$segs[isTip, 3], ret$segs[isTip, 4], 
                    phy$tip.label[phy$edge[isTip[-1], 2]], cex = cex, 
                    pos = 2, srt = 90, offset = 0)
            }
        }
    }
    index <- order(as.numeric(rownames(ret$segs)))
    if (show) {
        if (method == "phylogram") {
            assign("last_plot.phylo", list(type = "phylogram", 
                direction = direction, Ntip = phy$Nnode + 1, 
                Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                  3], yy = ret$segs[index, 4], pp = par(no.readonly = TRUE)), 
                envir = .PlotPhyloEnv)
        }
        else if (method == "polar") {
            assign("last_plot.phylo", list(type = "fan", Ntip = phy$Nnode + 
                1, Nnode = phy$Nnode, edge = phy$edge, xx = ret$segs[index, 
                3], yy = ret$segs[index, 4], theta = ret$segs[index, 
                5], rb = rb, pp = par(no.readonly = TRUE)), envir = .PlotPhyloEnv)
        }
        if (legend) {
            addBAMMlegend(x = list(coords = ret$segs[-1, ], colorbreaks = colorbreaks, 
                palette = colorobj$colpalette, colordens = colorobj$colsdensity), 
                location = "right")
        }
    }
    if (par.reset) {
        par(op)
    }
    invisible(list(coords = ret$segs[-1, ], colorbreaks = colorbreaks, 
        palette = colorobj$colpalette, colordens = colorobj$colsdensity))
}










#from https://github.com/macroevolution/bammtools/blob/587117a2a79159cc9809e346d7dedf58a4aa2add/BAMMtools/R/bammcolors.R
palettes <- list(

BrBG =     rev(c("#543005","#8c510a","#bf812d","#dfc27d","#f6e8c3",

                 "#f5f5f5","#c7eae5","#80cdc1","#35978f","#01665e",

                 "#003c30")),

PiYG =     rev(c("#8e0152","#c51b7d","#de77ae","#f1b6da","#fde0ef",

                 "#f7f7f7","#e6f5d0","#b8e186","#7fbc41","#4d9221",

                "#276419")),

PRGn =     rev(c("#40004b","#762a83","#9970ab","#c2a5cf","#e7d4e8",

                 "#f7f7f7","#d9f0d3","#a6dba0","#5aae61","#1b7837",

                "#00441b")),

PuOr =     rev(c("#7f3b08","#b35806","#e08214","#fdb863","#fee0b6",

                 "#f7f7f7","#d8daeb","#b2abd2","#8073ac","#542788",

                 "#2d004b")),

RdBu =     rev(c("#67001f","#b2182b","#d6604d","#f4a582","#fddbc7",

                 "#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac",

                 "#053061")),

RdYlBu =   rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee090",

                 "#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4",

                 "#313695")),

BuOr =     c("#002bff","#1a66ff","#3399ff","#66CCff","#99eeff",

            "#ccffff","#ffffcc","#ffee99","#ffee66","#ff9933",

            "#ff661a","#ff2b00"),

BuOrRd =   c("#085aff","#3377ff","#5991ff","#8cb2ff","#bfd4FF",

             "#e6eeff","#f7faff","#ffffcc","#ffff99","#ffff00",

             "#ffcc00","#ff9900","#ff6600","#ff0000"),

DkRdBu =   c("#2a0bd9","#264eff","#40a1ff","#73daff","#abf8ff",

             "#e0ffff","#ffffbf","#ffe099","#ffad73","#f76e5e",

             "#d92632","#a60021"),

BuDkOr =   c("#1f8f99","#52c4cc","#99faff","#b2fcff","#ccfeff",

             "#e6ffff","#ffe6cc","#ffca99","#ffad66","#ff8f33",

             "#cc5800","#994000"),

GnPu =     c("#005100","#008600","#00bc00","#00f100","#51ff51",

             "#86ff86","#bcffbc","#ffffff","#fff1ff","#ffbcff",

             "#ff86ff","#ff51ff","#f100f1","#bc00bc","#860086",

             "#510051"),

RdYlGn =   rev(c("#a50026","#d73027","#f46d43","#fdae61","#fee08b",

                 "#ffffbf","#d9ef8b","#a6d96a","#66bd63","#1a9850",

                 "#006837")),

Spectral = rev(c("#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b",

             "#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd",

             "#5e4fa2")),

grayscale = c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000"),

revgray = rev(c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000")),

greyscale = c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000"),

revgrey = rev(c("#ffffff","#f0f0f0","#d9d9d9","#bdbdbd","#969696","#737373","#525252",

        "#252525",

        "#000000"))

);

.colorEnv <- new.env();

assign("palettes", palettes, env = .colorEnv);


