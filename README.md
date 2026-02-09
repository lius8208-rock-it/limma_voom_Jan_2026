# This is differential gene expression analysis for population as random effect using limma-voom
#In this script, I first compare the model with GBL population placement with coastR or coastL
# The results showed that regional and interaction of region and treatment effect better explained total variance when placed GBL with coastR 
# Second analysis was PCA with and without introduced model. Should do PCA without introduced model as people would like to see, without any manipulation, transcript global expression pattern
#Third analysis, DE analysis, as GBL assigned to coastR, there is only one population in one region under one treatment, there's a confounding effect. Population was blocked.
# GO enrichment and KEGG enrichment analysis for significantly expression genes across region, treatment and interaction of these two effects. 
