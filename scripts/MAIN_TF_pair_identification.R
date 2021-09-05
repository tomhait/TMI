# A script for identifying TF-TF pairs in Enhancer-Promoter links

#lib.loc =  "/specific/a/home/cc/cs/tomhait/R/x86_64-pc-linux-gnu-library/3.5"
#.libPaths(c(lib.loc,.libPaths()))

# required libraries
require(BSgenome.Hsapiens.UCSC.hg19);require(seqinr); require(topGO)
require(parallel);require(plyr)

genome <- BSgenome.Hsapiens.UCSC.hg19
source('/scripts/FUNCTIONS_motif_finding.R')
# location of meme suite - should include the FIMO tool
Sys.setenv(PATH = paste(Sys.getenv("PATH"),"/home/elkon/tools/meme/bin",sep = ":"))

### Variables:
qValTh <- 0.1 # Q-value.
k <- 10 # k closest enhancers to each gene.
win.size <- 500 # upstream/downstream bases from each center position.
mc.cores <- 40

# raw data - EP links
folder <- '/home/gaga/html/ct-focs/'
data_folder <- 'data/'
script_directory <- 'scripts/'
data_type <- 'encode'
tmp_directory <- 'tmp/'

# digitial genomic footprints
dgf.bs <- readRDS(file="/home/gaga/tomhait/projects/ePPI/rds/dgf.bs.rds") 
g_to_e = readRDS(paste(folder,data_folder,data_type,'/',data_type,'.cand.enh.rds',sep=''))

enh.bs <- readRDS(file=paste0(folder,data_folder,data_type,'/',data_type,'.enh.pos.rds')) # Enhancer genomic positions.
prom.bs <- readRDS(file=paste0(folder,data_folder,data_type,'/',data_type,'.prom.pos.rds')) # Genes positions.

sample_annot <- readRDS(file=paste0(folder,data_folder,data_type,'/',data_type,'.sample.annot.rds')) # Sample annotations.
ep_links.ct <- readRDS(file=paste0(folder,data_folder,data_type,'/',data_type,'.EP_links.cell.rds')) # CT-FOCS EP links.
df.annot.hoc <- readRDS(file="/home/gaga/tomhait/projects/ePPI/rds/df.annot.hoc.rds") # TF hocomoco.

names(dgf.bs) <- paste0(seqnames(dgf.bs),":",start(dgf.bs),"-",end(dgf.bs))

num_links_per_cell <- sort(sapply(ep_links.ct,function(x) length(unlist(x))))
num_links_per_cell <- num_links_per_cell[num_links_per_cell>=50] # include only cell types with at least 50 predicted EP links
cells <- names(num_links_per_cell)
cells <- sort(cells)

tf_pair.list <- rep(list(NULL),length(cells))
names(tf_pair.list) <- cells

metrics.list <- c("num of EP", "num of E", "num of P", "raw width E", "window width E", "raw width P", "window width P", "DGF hits union",
					"DGF hits E", "DGF hits P", "FIMO results E", "Passed qval E", "FIMO results P", "Passed qval P", "TFs for E", "TFs for P",
					"TF pairs", "TRUE", "FALSE", "Recieved pairs", "Final pairs (after 0.1 filtration)", "Final E TFs", "Final P TFs", "UNION")
statistics.mat <- matrix(0,nrow=length(cells),ncol=length(metrics.list))
rownames(statistics.mat) <- cells
colnames(statistics.mat) <- metrics.list

summary.list <- c("Min", "1st Q", "Median", "Mean", "3rd Q", "Max")
summary.mat <- matrix(0,nrow=length(cells),ncol=length(summary.list))
rownames(summary.mat) <- cells
colnames(summary.mat) <- summary.list

for(i in 1:length(cells)){

cell_type <- cells[i]
print(paste0(i,' ',cell_type))
ep_links <- ep_links.ct[[cell_type]]

tmp <- ep_links
(statistics.mat[cell_type,"num of EP"] <- length(unlist(tmp)))
(statistics.mat[cell_type,"num of E"] <- length(unique(unlist(tmp))))
(statistics.mat[cell_type,"num of P"] <- length(tmp))

prom.bs.all <- prom.bs[names(ep_links)]
enh.bs.all <- enh.bs[unique(unlist(ep_links))]

statistics.mat[cell_type,"raw width P"] <- width(prom.bs.all)[1]
tmp <- prom.bs.all
shift.size <- abs(start(tmp)-end(tmp))/2
tmp <- shift(tmp,shift=shift.size)
tmp <- promoters(tmp,upstream=win.size,downstream=win.size)
prom.bs.all <- tmp
statistics.mat[cell_type,"window width P"] <- width(prom.bs.all)[1]

statistics.mat[cell_type,"raw width E"] <- width(enh.bs.all)[1]
tmp <- enh.bs.all
shift.size <- abs(start(tmp)-end(tmp))/2
tmp <- shift(tmp,shift=shift.size)
tmp <- promoters(tmp,upstream=win.size,downstream=win.size)
enh.bs.all <- tmp
statistics.mat[cell_type,"window width E"] <- width(enh.bs.all)[1]

TF_count.p <- matrix(0, nrow=length(prom.bs.all), ncol=dim(df.annot.hoc)[1])
rownames(TF_count.p) <- names(prom.bs.all)
colnames(TF_count.p) <- df.annot.hoc$TF

TF_count.e <- matrix(0, nrow=length(enh.bs.all), ncol=dim(df.annot.hoc)[1])
rownames(TF_count.e) <- names(enh.bs.all)
colnames(TF_count.e) <- df.annot.hoc$TF

### Find overlaps with dgf:
hit.p <- findOverlaps(dgf.bs,prom.bs.all,ignore.strand=TRUE,type="within")
hit.e <- findOverlaps(dgf.bs,enh.bs.all,ignore.strand=TRUE,type="within")
dgf.id <- union(queryHits(hit.p),queryHits(hit.e))
statistics.mat[cell_type,"DGF hits union"] <- length(dgf.id)

dgf.bs.reduced.p <- dgf.bs[sort(queryHits(hit.p))]
dgf.bs.reduced.e <- dgf.bs[sort(queryHits(hit.e))]
statistics.mat[cell_type,"DGF hits E"] <- length(dgf.bs.reduced.e)
statistics.mat[cell_type,"DGF hits P"] <- length(dgf.bs.reduced.p)

### Run FIMO on target sequences:
fimo.gr.e <- runFimo(dgf.bs.reduced.e,genome=genome, type="encode")
fimo.gr.p <- runFimo(dgf.bs.reduced.p,genome=genome, type="encode")
statistics.mat[cell_type,"FIMO results E"] <- length(fimo.gr.e)
statistics.mat[cell_type,"FIMO results P"] <- length(fimo.gr.p)

fimo.gr.p.q <- fimo.gr.p[as.numeric(fimo.gr.p$qvalue)<=qValTh]
fimo.gr.e.q <- fimo.gr.e[as.numeric(fimo.gr.e$qvalue)<=qValTh]
statistics.mat[cell_type,"Passed qval E"] <- length(fimo.gr.p.q)
statistics.mat[cell_type,"Passed qval P"] <- length(fimo.gr.e.q)

hit <- hit.p
fimo.gr <- fimo.gr.p.q
TF_count <- TF_count.p
tmp.bs <- prom.bs.all
shit <- unique(subjectHits(hit))
for(s in shit){
	a <- queryHits(hit)[subjectHits(hit)==s]
	a <- names(dgf.bs)[a]
	m <- match(fimo.gr$orig_seqnames,a)
	TFs <- fimo.gr$TF[!is.na(m)]
	TFs <- table(TFs)
	hoc_names = names(TFs)
	names(TFs) = as.character(df.annot.hoc[hoc_names,"TF"])
	TF_count[names(tmp.bs)[s],match(names(TFs),colnames(TF_count))] <- TFs
}
TF_count.p <- TF_count

hit <- hit.e
fimo.gr <- fimo.gr.e.q
TF_count <- TF_count.e
tmp.bs <- enh.bs.all
shit <- unique(subjectHits(hit))
for(s in shit){
	a <- queryHits(hit)[subjectHits(hit)==s]
	a <- names(dgf.bs)[a]
	m <- match(fimo.gr$orig_seqnames,a)
	TFs <- fimo.gr$TF[!is.na(m)]
	TFs <- table(TFs)
	hoc_names = names(TFs)
	names(TFs) = as.character(df.annot.hoc[hoc_names,"TF"])
	TF_count[names(tmp.bs)[s],match(names(TFs),colnames(TF_count))] <- TFs
}
TF_count.e <- TF_count

tb <- colSums(TF_count.p)
TFs.p <- names(tb[tb!=0])
p_rows <- rownames(TF_count.p)[apply(TF_count.p,1,function(row) any(row>0))]
statistics.mat[cell_type,"TFs for P"] <- length(TFs.p)

tb <- colSums(TF_count.e)
TFs.e <- names(tb[tb!=0])
e_rows <- rownames(TF_count.e)[apply(TF_count.e,1,function(row) any(row>0))]
statistics.mat[cell_type,"TFs for E"] <- length(TFs.e)

pairwiseTF.mat <- matrix(0,nrow=length(unlist(ep_links)),ncol=length(TFs.e)*length(TFs.p))
pairNames <- as.character(unlist(sapply(TFs.e,function(x) paste0(x,"#",TFs.p))))
ep_names <- unlist(sapply(names(ep_links), function(x) paste0(ep_links[[x]],"#",x)))
names(ep_names) <- NULL
rownames(pairwiseTF.mat) <- ep_names
colnames(pairwiseTF.mat) <- pairNames
statistics.mat[cell_type,"TF pairs"] <- dim(pairwiseTF.mat)[2]

e_tfs <- unlist(sapply(pairNames,function(x) unlist(strsplit(x,"#"))[1]))
p_tfs <- unlist(sapply(pairNames,function(x) unlist(strsplit(x,"#"))[2]))
e_tfs_ids <- match(e_tfs,colnames(TF_count.e))
p_tfs_ids <- match(p_tfs,colnames(TF_count.p))

for(ep in rownames(pairwiseTF.mat)){
	e_chr <- unlist(strsplit(ep,"#"))[1]
	p_chr <- unlist(strsplit(ep,"#"))[2]
	e_count <- TF_count.e[e_chr,e_tfs_ids]
	p_count <- TF_count.p[p_chr,p_tfs_ids]
	joint_tfs <- e_count * p_count
	names(joint_tfs) <- paste0(names(e_count),'#',names(p_count))
	joint_tfs[joint_tfs!=0] <- 1
	pairwiseTF.mat[ep,] <- joint_tfs
}

tb <- colSums(pairwiseTF.mat)
pair_names <- names(tb)
self_names <- sapply(pair_names,function(x) {tokens<-unlist(strsplit(x,"#")); return(tokens[1]!=tokens[2])})
statistics.mat[cell_type,"TRUE"] <- sum(self_names==TRUE)
statistics.mat[cell_type,"FALSE"] <- sum(self_names==FALSE)

tb <- tb[self_names]
statistics.mat[cell_type,"Recieved pairs"] <- length(tb[tb>0])

tb.f <- tb/length(unlist(ep_links))#tb[tb>=0.1*length(unlist(ep_links))]
statistics.mat[cell_type,"Final pairs (after 0.1 filtration)"] <- length(tb.f)
summ <- summary(tb.f)
summary.mat[cell_type,] <- t(matrix(summ))[1,]

tfs.e.f <- unique(unlist(sapply(names(tb.f),function(x) unlist(strsplit(x,"#"))[1])))
statistics.mat[cell_type,"Final E TFs"] <- length(tfs.e.f)
tfs.p.f <- unique(unlist(sapply(names(tb.f),function(x) unlist(strsplit(x,"#"))[2])))
statistics.mat[cell_type,"Final P TFs"] <- length(tfs.p.f)
statistics.mat[cell_type,"UNION"] <- length(union(tfs.e.f,tfs.p.f))

tf_pair.list[[cell_type]] <- tb.f
}

saveRDS(tf_pair.list,file="rds/TF_pairs.rds")
#write.csv(statistics.mat, "statistics.csv")
#write.csv(summary.mat, "summary.csv")

# Filter out TF pairs with Coefficient of variation (CV) below the 75% percentile of all CVs
vals <- unlist(tf_pair.list)
vals <- vals[vals>0]
log2_hits <- log2(vals)
summary(log2_hits)
sd(log2_hits)

tf_pairs <- unique(unlist(sapply(tf_pair.list, function(x) unique(names(x)))))

mat <- mclapply(tf_pairs, function(x){vals <- sapply(tf_pair.list,function(y){if(is.na(y[x])) return(0); return(y[x])}); if(sd(vals)==0) return(NULL); return(vals)},mc.cores=20)
names(mat) <- tf_pairs
mat[sapply(mat,is.null)] <- NULL

vals <- unlist(mat)
vals <- vals[vals>0]
log2_hits <- log2(vals)
summary(log2_hits)
sd(log2_hits)

sd_vals <- sapply(mat,sd)
mean_vals <- sapply(mat,mean)
cv <- log2(sd_vals/mean_vals)

table(cv>=quantile(cv,probs=0.75))

idx <- which(cv>=quantile(cv,probs=0.75))

mat_filt <- mat[idx]

tf_filt.list <- lapply(names(mat_filt),function(x) names(tf_pair.list)[which(mat_filt[[x]]!=0)])
names(tf_filt.list) <- names(mat_filt)

tf_pair.list.filt <- inverseList(tf_filt.list)

saveRDS(tf_pair.list.filt,file="rds/TF_pairs.filt_by_cv.rds")



