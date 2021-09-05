# A script for preparing GLADIATOR's files
lib.loc =  "/specific/a/home/cc/cs/tomhait/R/x86_64-pc-linux-gnu-library/3.5"
.libPaths(c(lib.loc,.libPaths()))

require(parallel)
require(igraph)

folder <- 'ct-focs/'
data_folder <- 'data/'
script_directory <- 'scripts/'
data_type <- 'encode'
tmp_directory <- 'tmp/'

# Load ENCODE gene expression data
exp.data <- readRDS(file="/home/gaga/tomhait/projects/ePPI/rds/Sheffield_GR_2013.exp.data.rds")
#exp.data <- log2(exp.data+1)

# Load V11 physical interactions with score >=500:
head(string_human <- read.table('data/human_ppi.direct.500.txt', header=TRUE))

# Create igraph object:
ppi.gp <- graph.data.frame(string_human[,1:2], directed=FALSE)
ppi.gp <- simplify(ppi.gp)

# Load proteins:
proteins <- read.delim('data/9606.protein.info.v11.0.txt', header=TRUE, quote="", fill=FALSE)

# Load 402 HOCOMOCOv11 TFs:
df.annot.hoc <- readRDS(file="rds/df.annot.hoc.rds")

# Load the normalized betweeness matrix
betweenness.in.tfs.sp.norm <- read.csv(file="data/TF-TF_involved_proteins_table.csv", header=T, row.names=1)

# Load TF-TF motif count matrix
tf_pair.list <- readRDS(file="rds/TF_pairs.filt_by_cv.rds")


# Protein names in ppi.gp (should be 376 TFs for String ppi above 500 combined score):
prot_in_net <- proteins$preferred_name[match(names(V(ppi.gp)),as.character(proteins$protein_external_id))]

# create Interactome.tsv for GLADIATOR
df <- string_human
df$combined_score <- NULL
df$comment <- rep("STRING_physical_int_V11",dim(df)[1])

write.table(df, file="GLADIATOR/data/Interactome.tsv", quote=F, sep="\t", row.names=F, col.names=F)

# create PhenSimMat.tsv for GLADIATOR
cell_types <- names(tf_pair.list)
cell_types <- gsub("-",".",cell_types)
tmp <- c('NT2.D1', 'NH.A', 'HMVEC.dAd', 'HRCEpiC', 'Monocytes.CD14+', 'SK.N.MC')
tmp.rep <- c('Ntera2', 'NHA', 'HMVECdAd', 'HRCE', 'CD14', 'SKNMC')
idx <- match(tmp.rep,colnames(exp.data))
colnames(exp.data)[idx] <- tmp
names(tf_pair.list) <- cell_types

filt_cell_types <- sort(intersect(cell_types,colnames(exp.data)))

tf_pair.list.filt <- tf_pair.list[filt_cell_types]
exp.data.filt <- exp.data[,filt_cell_types]
exp.data.filt <- exp.data.filt[which(apply(exp.data.filt,1,function(row) all(row!=0))),]


# plot coefficient of variation
sd_val <- apply(exp.data.filt,1,sd)
mean_val <- apply(exp.data.filt,1,mean)
cv <- sd_val/mean_val
mean_val <- log2(mean_val)
cv <- log2(cv)
plot(mean_val,cv,type='p')

# binning by log2(mean)
bins <- 20
th <- quantile(mean_val,probs=seq(0,1,by=1/bins))
inds <- rep(list(NULL), bins)
inds_above_med_cv <- rep(list(NULL), bins)
cv_med <- rep(0,bins)

for(i in 2:(bins+1)){
	tmp_idx <- which(mean_val>th[i-1] & mean_val<=th[i])
	inds[[i-1]] <- tmp_idx
	cv_med[i-1] <- quantile(cv[tmp_idx],probs=0.99)
	inds_above_med_cv[[i-1]] <- tmp_idx[cv[tmp_idx]>=cv_med[i-1]]
}

col <- rep('black',length(mean_val))
col[unlist(inds_above_med_cv)] <- 'red'
plot(mean_val,cv,type='p',col=col)

exp.data.cv.filt <- exp.data.filt[unlist(inds_above_med_cv),]
length(unlist(inds_above_med_cv))

# create PhenSimMat.tsv file for GLADIATOR
raw_exp <- exp.data.cv.filt
mat <- matrix(0,nrow=length(filt_cell_types),ncol=length(filt_cell_types))
rownames(mat) <- colnames(mat) <- filt_cell_types

for(i in 1:(length(filt_cell_types)-1)){
	cell_1 <- filt_cell_types[i]
	sqrt_cell_1 <- sqrt(sum(raw_exp[,cell_1]*raw_exp[,cell_1]))
	for(j in (i+1):length(filt_cell_types)){
		cell_2 <- filt_cell_types[j]
		numer <- sum(raw_exp[,cell_1]*raw_exp[,cell_2])
		sqrt_cell_2 <- sqrt(sum(raw_exp[,cell_2]*raw_exp[,cell_2]))
		sim_cell <- numer/(sqrt_cell_1*sqrt_cell_2)
		mat[cell_1,cell_2] <- sim_cell
		mat[cell_2,cell_1] <- sim_cell
	}
}


summary(apply(mat,1,sd))
summary(apply(mat,1,mean))

df <- as.dataframe(mat)

df <- cbind(cell_type=rownames(df),df)

write.table(df, file='GLADIATOR/data/PhenSimMat.tsv', quote=F, sep='\t', row.names=F)

## create initial SeedPSfName.tsv containing for each cell type the set of TFs in tf_pair.list.filt
cell_types <- names(tf_pair.list.filt)
df <- data.frame(cell_type=cell_types, num_genes = rep(0,length(cell_types)), num_prot=rep(0,length(cell_types)), 
		num_dumb=rep(0,length(cell_types)), TFs=rep("",length(cell_types)))

rownames(df) <- cell_types
df$TFs <- as.character(df$TFs)

for(s.cell in cell_types){
	tf_pairs <- names(tf_pair.list.filt[[s.cell]])
	first_tf <- unique(sapply(tf_pairs,function(x) unlist(strsplit(x,"#"))[1]))
	second_tf <- unique(sapply(tf_pairs,function(x) unlist(strsplit(x,"#"))[2]))
	tfs <- union(first_tf,second_tf)
	m <- match(tfs,df.annot.hoc$TF)
	tfs_string <- df.annot.hoc$string_name[m]
	m <- match(proteins$preferred_name,tfs_string)
	m <- which(!is.na(m))
	tfs_ens <- proteins$protein_external_id[m]
	tfs_ens <- intersect(tfs_ens, names(V(ppi.gp)))
	df[s.cell,'num_prot'] <- length(tfs_ens)
	df[s.cell,'TFs'] <- paste(tfs_ens,collapse=";")
}

write.table(df, file='GLADIATOR/data/KnownDisPS.tsv', quote=F, sep='\t', row.names=F, col.names=F)

df_new <- df

df_new$TFs <- sapply(df$TFs, function(x) gsub(';',',',x))


write.table(df_new[,c('cell_type','TFs')], file='GLADIATOR/data/SeedPS.tsv', quote=F, sep='\t', row.names=F, col.names=F)

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



# create Interactome.tsv for GLADIATOR with added edges between pairs of TFs
df <- string_human
df$combined_score <- NULL
df$comment <- rep("STRING_physical_int_V11",dim(df)[1])

cell_types <- names(tf_pair.list.filt)

virt_edges <- data.frame()

for(s.cell in cell_types){
	tf_pairs <- names(tf_pair.list.filt[[s.cell]])
	first_tf <- sapply(tf_pairs,function(x) unlist(strsplit(x,"#"))[1])
	second_tf <- sapply(tf_pairs,function(x) unlist(strsplit(x,"#"))[2])
	tfs <- union(first_tf,second_tf)
	tmp_list <- tf_to_ens_list[tfs]
	tmp_list[which(sapply(tmp_list,length)==0)] <- NULL
	tfs <- names(tmp_list)
	m1 <- match(first_tf, tfs)
	m2 <- match(second_tf, tfs)
	inds <- which(!is.na(m1) & !is.na(m2))
	first_tf <- first_tf[inds]
	second_tf <- second_tf[inds]
	res <- lapply(1:length(first_tf), function(i) {f1 <- tf_to_ens_list[[first_tf[i]]]; f2 <- tf_to_ens_list[[second_tf[i]]];
					n1 <- length(f1); n2 <- length(f2); f1rep <- rep(f1, n2); f2rep <- rep(f2, n1);
					return(data.frame(protein1=f1rep,protein2=f2rep))})
	res <- do.call(rbind,res)
	res$comment <- rep('virtual_edge',dim(res)[1])
	virt_edges <- rbind(virt_edges,res)

}

virt_edges <- unique(virt_edges)

df <- rbind(df,virt_edges)

write.table(df, file="GLADIATOR/data/Interactome.edgesBetweenTFpairs.tsv", quote=F, sep="\t", row.names=F, col.names=F)


