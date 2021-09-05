### A script for creating the normalized betweeness of every protein in the shortest paths between TF pairs ###

#lib.loc =  "/specific/elkon/amitlevon/R/x86_64-pc-linux-gnu-library/3.5"
lib.loc =  "/specific/a/home/cc/cs/tomhait/R/x86_64-pc-linux-gnu-library/3.5"
.libPaths(c(lib.loc,.libPaths()))

# require libraries
require(parallel)
require(igraph)

folder <- 'ct-focs/'
data_folder <- 'data/'
script_directory <- 'scripts/'
data_type <- 'encode'
tmp_directory <- 'tmp/'

summary.list <- c("Min", "1st Q", "Median", "Mean", "3rd Q", "Max")
med.summary.mat <- matrix(0,nrow=1,ncol=length(summary.list))
rownames(med.summary.mat) <- c("stats")
colnames(med.summary.mat) <- summary.list

metrics.list <- c("total string interactions", "ppi vertices", "ppi edges", "num of proteins", "proteins in network", "protein names in ppi",
					"avg edges between pairs of proteins", "avg edges between pairs of TFs", "TFs with path to different TFs",
					"TFs that each TF can get to", "MED betweenness sum")
statistics.mat <- matrix(0,nrow=1,ncol=length(metrics.list))
rownames(statistics.mat) <- c("stats")
colnames(statistics.mat) <- metrics.list

# Load V11 physical interactions with score >=500:
head(string_human <- read.table('data/human_ppi.direct.500.txt', header=TRUE))
statistics.mat["stats","total string interactions"] <- dim(string_human)[1]

# Create igraph object:
ppi.gp <- graph.data.frame(string_human[,1:2], directed=FALSE)
ppi.gp <- simplify(ppi.gp)
statistics.mat["stats","ppi vertices"] <- vcount(ppi.gp)
statistics.mat["stats","ppi edges"] <- ecount(ppi.gp)

# Load proteins:
proteins <- read.delim('data/9606.protein.info.v11.0.txt', header=TRUE, quote="", fill=FALSE)
statistics.mat["stats","num of proteins"] <- dim(proteins)[1]
head(proteins[match(names(V(ppi.gp)),proteins$protein_external_id),])
tail(proteins[match(names(V(ppi.gp)),proteins$protein_external_id),])

# Load 402 HOCOMOCOv11 TFs:
df.annot.hoc <- readRDS(file="rds/df.annot.hoc.rds")

# Protein names in ppi.gp (should be 376 TFs for String ppi above 500 combined score):
prot_in_net <- proteins$preferred_name[match(names(V(ppi.gp)),as.character(proteins$protein_external_id))]
statistics.mat["stats","proteins in network"] <- length(prot_in_net)
statistics.mat["stats","protein names in ppi"] <- length(intersect(df.annot.hoc$string_name,prot_in_net))
sort(setdiff(df.annot.hoc$string_name,prot_in_net)) # 26 missing TFs

# Go over pairwise proteins and find the shortest paths between each pair:
dist_tb <- distance_table(ppi.gp)
vals <- dist_tb$res
names(vals) <- as.character(1:length(vals))

png(file="shortest_paths[between pairs of proteins].png")

barplot(vals,xlab="Number of edges in the shortest path", ylab="Number of shortest paths", main="Number of shortest paths \n 
		between pairs of proteins in the PPI")

dev.off()
statistics.mat["stats","avg edges between pairs of proteins"] <- sum(vals*as.numeric(names(vals)))/sum(vals)

# Go over pair-wise TFs and count how many direct/indirect interactions there are:
vert.tfs <- intersect(df.annot.hoc$string_name,proteins$preferred_name)
vert.tfs <- proteins$protein_external_id[match(vert.tfs,proteins$preferred_name)]
v = match(as.character(vert.tfs),names(V(ppi.gp)))
v <- v[!is.na(v)]
dist_tb.tfs <- distances(ppi.gp,v=v,to=v)
x <- upper.tri(dist_tb.tfs, diag = FALSE)
vals.tf <- dist_tb.tfs[x]
tb <- table(dist_tb.tfs[x])

png(file="shortest_paths[between pairs of TFs].png")

barplot(tb[1:length(tb)],xlab="Number of edges in the shortest path", ylab="Number of shortest paths", main="Number of shortest paths \n 
		between pairs of TFs (n=376) in the PPI")

dev.off()

# Average num of edges:
statistics.mat["stats","avg edges between pairs of TFs"] <- sum(tb[1:length(tb)]*as.numeric(names(tb[1:length(tb)])))/sum(tb[1:length(tb)])

# How many TFs have a path to a different TF?
flags <- rep(FALSE,dim(dist_tb.tfs)[1])
for(i in 1:(dim(dist_tb.tfs)[1]-1)){
	if(flags[i]) next
	for(j in (i+1):dim(dist_tb.tfs)[1]){
		if(dist_tb.tfs[i,j]==Inf) next
		flags[i] = TRUE
		flags[j] = TRUE
	}
}

statistics.mat["stats","TFs with path to different TFs"] <- sum(flags)

# How many TFs can each TF get to?
sum_reachable_tfs <- rep(0,dim(dist_tb.tfs)[1])
for(i in 1:(dim(dist_tb.tfs)[1]-1)){
	for(j in (i+1):dim(dist_tb.tfs)[1]){
		if(dist_tb.tfs[i,j]==Inf) next
		sum_reachable_tfs[i]=sum_reachable_tfs[i]+1
		sum_reachable_tfs[j]=sum_reachable_tfs[j]+1
	}
}

statistics.mat["stats","TFs that each TF can get to"] <- as.numeric(names(table(sum_reachable_tfs)))

proteins[match(rownames(dist_tb.tfs)[which(sum_reachable_tfs==0)],proteins$protein_external_id),]


### Check the MED## proteins - part of the mediator complex ###

ind <- grep("^MED\\d+|^CCNC|CDK8|CDK19",proteins$preferred_name)
sort(med_names <- proteins$preferred_name[ind])
ppi.ind.med <- match(proteins$protein_external_id[ind],names(V(ppi.gp)))
med_names <- med_names[!is.na(ppi.ind.med)]
ppi.ind.med <- ppi.ind.med[!is.na(ppi.ind.med)]
ppi.med.ens <- names(V(ppi.gp))[ppi.ind.med]
med_betweenness <- betweenness(ppi.gp, v = V(ppi.gp)[ppi.ind.med], directed = FALSE, nobigint = TRUE, normalized = TRUE)

names(med_betweenness) <- med_names
med_table <- sort(med_betweenness)
write.csv(med_table, "MED_table.betweeness.all_protein_pairs.csv")

statistics.mat["stats","MED betweenness sum"] <- sum(med_betweenness)*100

med_summ <- summary(med_betweenness)
med.summary.mat["stats",] <- t(matrix(med_summ))[1,]
write.csv(med.summary.mat, "MED_summary.betweeness.all_protein_pairs.csv")
write.csv(statistics.mat, "network_statistics.csv")

### Measure betweennes of TFi TFj pairs for each cell type ###

sort(ppi_names <- proteins$preferred_name)
ppi.ind <- match(proteins$protein_external_id,names(V(ppi.gp)))
ppi_names <- ppi_names[!is.na(ppi.ind)]
ppi.ind <- ppi.ind[!is.na(ppi.ind)]
ppi.ens <- names(V(ppi.gp))[ppi.ind]

tf_to_ens_list <- rep(list(NULL), length(unique(df.annot.hoc$TF)))
names(tf_to_ens_list) <- unique(df.annot.hoc$TF)
for(tf in names(tf_to_ens_list)){
	m <- match(tf,df.annot.hoc$TF)
	tf_string <- df.annot.hoc$string_name[m]
	m <- match(proteins$preferred_name,tf_string)
	m <- which(!is.na(m))
	tfs_ens <- proteins$protein_external_id[m]
	tfs_ens <- intersect(tfs_ens, names(V(ppi.gp)))
	tf_to_ens_list[[tf]] <- tfs_ens
}

tf_pair.list <- readRDS(file="rds/TF_pairs.filt_by_cv.rds") # Load TF-TF motif count matrix.
cell_types <- names(tf_pair.list)
list_bw <- rep(list(NULL),length(cell_types))
names(list_bw) <- cell_types

celltype.metrics.list <- c("Cell types TFi-TFj involved proteins betweenness sum")
celltype.statistics.mat <- matrix(0,nrow=length(cell_types),ncol=length(celltype.metrics.list))
rownames(celltype.statistics.mat) <- cell_types
colnames(celltype.statistics.mat) <- celltype.metrics.list

celltype.summary.mat <- matrix(0,nrow=length(cell_types),ncol=length(summary.list))
rownames(celltype.summary.mat) <- cell_types
colnames(celltype.summary.mat) <- summary.list

for(s.cell in cell_types){

	print(s.cell)
	cell.pairs <- tf_pair.list[[s.cell]]
	betweenness.in.tfs.sp <- rep(0, length(ppi.ens))
	names(betweenness.in.tfs.sp) <- ppi.ens

	tot_sp <- 0
	for(pair in cell.pairs){
		#print(pair)
		pair.split <- strsplit(pair,"#")[[1]]
		pro.1.ens <- tf_to_ens_list[pair.split[1]]
		pro.1.ind <- match(as.character(pro.1.ens),names(V(ppi.gp)))
		if(is.na(pro.1.ind)){ # Protein not in network #
			next}
		pro.2.ens <- tf_to_ens_list[pair.split[2]]
		pro.2.ind <- match(as.character(pro.2.ens),names(V(ppi.gp)))
		if(is.na(pro.2.ind)){ # Protein not in network #
			next}
		sh_paths <- all_shortest_paths(ppi.gp, from = pro.1.ind, to = pro.2.ind)
		sh_paths <- sh_paths$res
		if(length(sh_paths)==0) next
		tb <- sapply(sh_paths,function(x) names(x)[length(x)])
		tb <- table(tb)
		tot_paths_per_tg_tf <- tb
		names(tot_paths_per_tg_tf) <- names(tb)
		inds <- cumsum(tot_paths_per_tg_tf)
		res <- mclapply(sh_paths, function(x) names(x),mc.cores=20) # Taking all proteins in the ppi net and intersecting with every set of proteins x will give us x.
		tbl <- mclapply(1:length(inds),function(j) {if(j==1) {return(table(unlist(res[1:inds[j]])))}; return(table(unlist(res[(inds[j-1]+1):inds[j]])))},mc.cores=20)
		tbl <- lapply(1:length(tbl), function(j) {x <- tbl[[j]]; if(length(x)==0) {return(NULL)};  return(x/tot_paths_per_tg_tf[j])})
		names(tbl) <- names(inds)
		tbl[sapply(tbl,is.null)] <- NULL
		for (x in tbl){ betweenness.in.tfs.sp[names(x)] <- betweenness.in.tfs.sp[names(x)] + x}
		tot_sp <- tot_sp + length(sh_paths)
	}

	### Make a data.frame with gene symbols and betweeness values ###
	betweenness.in.tfs.sp.norm <- data.frame(gene_sym = as.character(ppi_names), betweeness_value = betweenness.in.tfs.sp/tot_sp)
	rownames(betweenness.in.tfs.sp.norm) <- ppi.ens

	### Reorder from high to low betweeness ###
	betweenness.in.tfs.sp.norm <- betweenness.in.tfs.sp.norm[rev(order(betweenness.in.tfs.sp.norm$betweeness_value)),,drop=F]
	
	list_bw[[s.cell]] <- betweenness.in.tfs.sp.norm
	summ <- summary(betweenness.in.tfs.sp.norm$betweeness_value)
	celltype.summary.mat[s.cell,] <- t(matrix(summ))[1,]
	celltype.statistics.mat[s.cell,"Cell types TFi-TFj involved proteins betweenness sum"] <- sum(betweenness.in.tfs.sp.norm$betweeness_value)*100

}

saveRDS(list_bw,file="TFiTFj_betweenness.rds")

write.csv(celltype.summary.mat, "celltype_TFi-TFj_involved_proteins_summary.csv")

write.csv(celltype.statistics.mat, "celltype_network_statistics.csv")

