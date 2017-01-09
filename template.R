#################################################################################
# Some Notes:
# Only analyses with the same set of initial entities are allowed (i.e. they must
# have the same IMap file).

# To Do: 
# - get a summary of the ESS values
# - get priors and implement a function to plot the prior distributions
# - implement detection of analyses that were stuck in the one species model
# - a lot of testing


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# This script can be used to summarize the results of a BPPmulti run.
# 
# Load ape and the functions of BPPmultitool.
source("subs.R")
#
# Set the working directory to the one that was specified in BPPmulti.pl
# 
setwd("/home/jonas/my_BPPmulti_run")


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# The next step will collect information on the BPPmulti run based on the folders
# and control files in the working directory.
# 
settings <- getSettings()


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# plot the guide trees
# 
plotGuideTrees(filename="BPP_guide_trees.pdf", settings=settings, label.offset=.1, cex=.6)


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# parse all BPP output - this may take a while...
# 
BPPresults <- parseBPPmulti(settings)


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# a look at the ESS values
#
BPPresults$ESS$repeat0 # all values of one repeat
lapply(BPPresults$ESS$repeat0, min) # minimum values


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# some first checks of single analyses
# - for "mol-NNI"-sets, the BPPresultsMatrix is for instance BPPresults$unguided_mean_supports or BPPresults$unguided$repeat0
# - for other sets, the BPPresultsMatrix is for instance BPPresults$mean_supports or BPPresults$supports$repeat0
# 
extractSet(BPPresults$supports$repeat1, "Tree1", "intgr")
(check <- extractSet(BPPresults$supports$repeat2, "Tree1", "intgr"))
plot.pp(settings$trees$Tree1, check, nrow=3, ncol=3, sep=c(.17,.16), exp=.5, add.order=F)


#=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~
# write data tables for all analyses and plot
# 
write.BPP.tables("tables", settings, BPPresults)

pdf("BPPsummary.pdf", width=8, height=3.9, useDingbats=F)
plot.all.pp(settings, BPPresults, nrow=3, ncol=3, sep=c(.2,.18), exp=.5, cex=.75, edge.width=.5, add.order=F)
dev.off()
