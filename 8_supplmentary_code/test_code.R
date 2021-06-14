MEturquoise = sample(1:100,50)
MEblue      = sample(1:100,50)
MEbrown     = sample(1:100,50)
MEyellow    = sample(1:100,50) 
MEgreen     = c(MEyellow[1:30], sample(1:100,20))
MEred	    = c(MEbrown [1:20], sample(1:100,30))
MEblack	    = c(MEblue  [1:25], sample(1:100,25))
ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred, MEblack)
dat1   = simulateDatExpr(ME,300,c(0.2,0.1,0.08,0.051,0.05,0.042,0.041,0.3), signed=TRUE)
colorh = labels2colors(dat1$allLabels)
hubs    = chooseTopHubInEachModule(dat1$datExpr, colorh)
hubs
require(WGCNA,quietly=T)
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengene

