# TMI identification
lib.loc =  "/specific/a/home/cc/cs/tomhait/R/x86_64-pc-linux-gnu-library/3.5"
.libPaths(c(lib.loc,.libPaths()))

require(parallel)
require(igraph)

folder <- 'ct-focs/'
data_folder <- 'data/'
script_directory <- 'scripts/'
data_type <- 'encode'
tmp_directory <- 'tmp/'


# Load V11 physical interactions with score >=500:
head(string_human <- read.table('data/human_ppi.direct.500.txt', header=TRUE))
statistics.mat["stats","total string interactions"] <- dim(string_human)[1]

# Create igraph object:
ppi.gp <- graph.data.frame(string_human[,1:2], directed=FALSE)
ppi.gp <- simplify(ppi.gp)

# Load proteins:
proteins <- read.delim('data/9606.protein.info.v11.0.txt', header=TRUE, quote="", fill=FALSE)

# Load 402 HOCOMOCOv11 TFs:
df.annot.hoc <- readRDS(file="rds/df.annot.hoc.rds")

# Load the normalized betweeness matrix
betweenness.in.tfs.sp.norm <- read.csv(file="data/TF-TF_involved_proteins_table.csv", header=T, row.names=1)
# Filter out entries with 0 betweeness value
betweenness.in.tfs.sp.norm <- betweenness.in.tfs.sp.norm[betweenness.in.tfs.sp.norm$betweeness_value>0,]
# log2 transformation
betweenness.in.tfs.sp.norm$log2_bw <- log2(betweenness.in.tfs.sp.norm$betweeness_value)
hist(betweenness.in.tfs.sp.norm$log2_bw, breaks=50) #should be normal-like dist.
#hist(betweenness.in.tfs.sp.norm$betweeness_value, breaks=50) #should be normal-like dist.

# calculate right-sided 95% confidence interval threshold
a <- mean(betweenness.in.tfs.sp.norm$log2_bw,na.rm=T)
s <- sd(betweenness.in.tfs.sp.norm$log2_bw,na.rm=T)
n <- length(betweenness.in.tfs.sp.norm$log2_bw)
margin_of_error <- qt(.95,df=n-1)*s/sqrt(n) # qt from t-distribution, qnorm(.95) from z distribution
right_ci_th <- a+margin_of_error

hist(betweenness.in.tfs.sp.norm$log2_bw, breaks=50); abline(v=right_ci_th,col="red")

# Load TF-TF motif count matrix
tf_pair.list <- readRDS(file="rds/TF_pairs.filt_by_cv.rds")
cells <- names(tf_pair.list)

# load betweenness values per cell type
list_bw <- readRDS(file="rds/TFiTFj_betweenness.rds")

# Protein names in ppi.gp (should be 376 TFs for String ppi above 500 combined score):
prot_in_net <- proteins$preferred_name[match(names(V(ppi.gp)),as.character(proteins$protein_external_id))]


background_set <- betweenness.in.tfs.sp.norm$log2_bw
a <- mean(background_set)


names(background_set) <- betweenness.in.tfs.sp.norm$gene_sym

list_modules <- rep(list(NULL),length(cells))
names(list_modules) <- cells

pvals <- rep(1,length(cells))
names(pvals) <- cells

for(i in 1:length(cells)){
	print(i)
	s.cell <- cells[i]
	tg_set <- list_bw[[s.cell]]
	tg_set <- tg_set[tg_set$betweeness_value>0,]
	tg_set$log2_bw <- log2(tg_set$betweeness_value)
	tmp <- tg_set$log2_bw
	names(tmp) <- tg_set$gene_sym
	tg_set <- tmp
	tg_set <- sort(tg_set)
	pnames <- names(tg_set)
	min_pval <- t.test(tg_set[pnames[1:length(pnames)]],background_set[pnames[1:length(pnames)]],paired=F,alternative="greater")$p.value
	min_j <- 1
	for(j in 2:(length(pnames)-1)){
		new_pval <- t.test(tg_set[pnames[j:length(pnames)]],background_set[pnames[j:length(pnames)]],paired=F,alternative="greater")$p.value
		if(new_pval<min_pval){
			min_pval <- new_pval
			min_j <- j
		}
	}
	list_modules[[s.cell]] <- pnames[min_j:length(pnames)]
	pvals[s.cell] <- min_pval
}


saveRDS(list_modules,file='rds/cell_specific_modules.TMI.rds')
saveRDS(pvals,file='rds/TMI.pvals.rds')
