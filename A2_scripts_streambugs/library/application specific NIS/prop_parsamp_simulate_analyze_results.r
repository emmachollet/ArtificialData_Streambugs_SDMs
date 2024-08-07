# analyze results

# calculate results for parameter sample:
# =======================================

res.mat <- simulation.parsamp(par.samp=par.samp,
                              par.prior.delta=par.prior.delta,
                              f.survivals=simulate.streambugs.return.likeli,
                              name.run=name.run,tout=tout, y.names=y.names,M.taxa=M.taxa,
                              p.obs=p.obs,p.abs=p.abs,D.drift=D.drift,outputfolder=outputfolder)


# delete NA rows (if calculation is interrupted, res.mat is not returned! but can be read from file)
# res.mat <- read.table(paste(outputfolder,"/res_ss_sample_",name.run,".dat",sep=""),sep="\t",header=TRUE) 
# res.mat <- res.mat[,-ncol(res.mat)]

ind.na <- NULL

for(i in 1:nrow(res.mat))
{
  if( sum(is.na(res.mat[i,])) >0 ) ind.na <- c(ind.na,i)
}

if(length(ind.na)>0)
{
  par.samp.nona <- par.samp[,-ind.na]
  res.mat <- res.mat[-ind.na,]
} else par.samp.nona <- par.samp


par.prior.delta.matrix <- matrix(rep(par.prior.delta,each=ncol(par.samp.nona)),nrow=ncol(par.samp.nona))
colnames(par.prior.delta.matrix) <- names(par.prior.delta)
par.prior.delta.matrix <- t(par.prior.delta.matrix)

par.samp.nona <- rbind(par.samp.nona, par.prior.delta.matrix)

#Number of simulations
cat("No of Simulations: ",nrow(res.mat),"\n")

save.image(file = paste(outputfolder,"/streambugs_",name.run,".RData",sep=""))

##############################################################

Offset <- -max(res.mat[,"loglikeli"])
log.mean.likeli <- log(mean(exp(res.mat[,"loglikeli"]+Offset)))-Offset

cat("log mean likelihood: ",signif(log.mean.likeli,digits=4))


p.y.mat <- read.table(paste(outputfolder,"/res_py_sample_",name.run,".dat",sep=""),sep="\t",header=TRUE,blank.lines.skip=TRUE,fill=TRUE)
p.y.mat <-  p.y.mat[-ind.na,-ncol(p.y.mat)]
colnames(p.y.mat) <- gsub("X","",colnames(p.y.mat))

p.1.mat <- read.table(paste(outputfolder,"/res_p1_sample_",name.run,".dat",sep=""),sep="\t",header=TRUE,blank.lines.skip=TRUE,fill=TRUE)
p.1.mat <- p.1.mat[-ind.na,-ncol(p.1.mat)]
colnames(p.1.mat) <- gsub("X","",colnames(p.1.mat))

# p.0.mat <- read.table(paste(outputfolder,"/res_p0_sample_",name.run,".dat",sep=""),sep="\t",header=TRUE,blank.lines.skip=TRUE,fill=TRUE)
# p.0.mat <- p.0.mat[-ind.na,-ncol(p.0.mat)]
# colnames(p.0.mat) <- gsub("X","",colnames(p.0.mat))

comp.obs.pred.samp <- compare.observed.predicted.sample(p.y.mat,observed.abund,Invertebrates)

freq.obs <- calc.freq.observation(observed.abund,Invertebrates)

# pdf(paste("output/boxplot_freqobs_pobs",name.run,".pdf",sep=""),height=6,width=6.5)
# boxplot.freqobs.pobs(freq.obs,p.1.mat,cex=1.5,cex.axis=1.5)
# dev.off()

pdf(paste(outputfolder,"/boxplot_freqobs_pobs_binary",name.run,".pdf",sep=""),height=6,width=6.5)
boxplot.freqobs.pobs.binary(freq.obs,p.1.mat,cex=1.5,cex.axis=1.5)
dev.off()