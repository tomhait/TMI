# A script for enrichment analyses on detected cell type modules by each algorithm

# required libraries
library(httr)
library(stringi)
library(dplyr)
library(clusterProfiler)

# Load proteins:
proteins <- read.delim('/home/gaga/tomhait/projects/ePPI/data/9606.protein.info.v11.0.txt', header=TRUE, quote="", fill=FALSE)

# Load the normalized betweeness matrix - for background set
betweenness.in.tfs.sp.norm <- read.csv(file="/home/gaga/tomhait/projects/ePPI/data/TF-TF_involved_proteins_table.csv", header=T, row.names=1)
betweenness.in.tfs.sp.norm <- betweenness.in.tfs.sp.norm[betweenness.in.tfs.sp.norm$betweeness_value>0,]

gaf_file <- read.table("/home/gaga/tomhait/projects/ePPI/data/goa_human.gaf",sep = "\t",skip = 40,quote=NULL) # Get GO annotations.

## Choose one of the algorithms' output below
## 1. TMI Enrichment analysis
seedPs <- readRDS(file="/home/gaga/tomhait/projects/ePPI/rds/TF_pairs.filt_by_cv.rds")
modules_sym <- readRDS(file='/home/gaga/tomhait/projects/ePPI/rds/cell_specific_modules.TMI.rds')

summary(sapply(seedPs,length))

res <- lapply(names(seedPs), function(x) unlist(strsplit(as.character(seedPs[[x]]),"#")))
names(res) <- names(seedPs)
seedPs <- res
summary(sapply(seedPs,length))

summary(sapply(modules_sym,length))

summary(sapply(names(modules_sym),function(x) length(intersect(modules_sym[[x]],seedPs[[x]]))))
summary(sapply(names(modules_sym),function(x) length(setdiff(modules_sym[[x]],seedPs[[x]]))))
summary(sapply(names(modules_sym),function(x) length(setdiff(seedPs[[x]],modules_sym[[x]]))))

############################### GLADIATOR modules ################################################
## 2. GLADIATOR original Enrichment analysis
seedPs <- read.table(file='/home/gaga/tomhait/projects/ePPI/GLADIATOR/data/SeedPS.tsv', sep='\t', row.names=1, header=F)
modules <- read.table(file='/home/gaga/tomhait/projects/ePPI/GLADIATOR/data/GLADIATOR_modules.txt', sep='\t', row.names=1, header=T)


## 3. GLADIATOR virtual edges Enrichment analysis
seedPs <- read.table(file='/home/gaga/tomhait/projects/ePPI/GLADIATOR/data/SeedPS.tsv', sep='\t', row.names=1, header=F)
modules <- read.table(file='/home/gaga/tomhait/projects/ePPI/GLADIATOR/data/GLADIATOR_modules.withVirtEdges.txt', sep='\t', row.names=1, header=T)

# convert to lists
res <- lapply(rownames(seedPs), function(x) unlist(strsplit(as.character(seedPs[x,1]),",")))
names(res) <- rownames(seedPs)
seedPs <- res
summary(sapply(seedPs,length))

res <- lapply(rownames(modules), function(x) unlist(strsplit(as.character(modules[x,1]),", ")))
names(res) <- rownames(modules)
modules <- res
summary(sapply(modules,length))

summary(sapply(names(modules),function(x) length(intersect(modules[[x]],seedPs[[x]]))))
summary(sapply(names(modules),function(x) length(setdiff(modules[[x]],seedPs[[x]]))))
summary(sapply(names(modules),function(x) length(setdiff(seedPs[[x]],modules[[x]]))))

res <- lapply(1:length(modules),function(i) sort(sapply(names(modules[-i]), function(x) length(intersect(modules[[i]],unlist(modules[[x]]))))))
names(res) <- names(modules)
sapply(res,summary)

res <- lapply(1:length(modules),function(i) sort(sapply(names(modules[-i]), function(x) length(setdiff(modules[[i]],unlist(modules[[x]]))))))
names(res) <- names(modules)
sapply(res,summary)

res <- sapply(1:length(modules),function(i) length(setdiff(modules[[i]],unlist(modules[-i]))))
names(res) <- names(modules)
summary(res)

cell_types <- names(modules)

## convert to symbols
modules_sym <- rep(list(NULL), length(modules))
names(modules_sym) <- cell_types

for(s.cell in cell_types){
	tmp <- modules[[s.cell]]
	m <- match(tmp,proteins$protein_external_id)
	sym <- proteins$preferred_name[m]
	modules_sym[[s.cell]] <- unique(as.character(sym))
}

seedPs_sym <- rep(list(NULL), length(seedPs))
names(seedPs_sym) <- cell_types

for(s.cell in cell_types){
	tmp <- seedPs[[s.cell]]
	m <- match(tmp,proteins$protein_external_id)
	sym <- proteins$preferred_name[m]
	seedPs_sym[[s.cell]] <- unique(as.character(sym))
}

######################################End of GLADIATOR modules############################################

#
length(background_set <- unique(as.character(betweenness.in.tfs.sp.norm$gene_sym)))

length(unique(unlist(modules_sym)))
go_classes <- c("F","P","C")

res_list <- rep(list(NULL),3)
names(res_list) <- go_classes

for(j in 1:3){
	print(j)
	go_class <- go_classes[j]
	gaf_P <- gaf_file[gaf_file$V9 == go_class,]
	term2gene <- data.frame(gaf_P$V5,gaf_P$V3)
	all_results <- data.frame()
	for(i in 1:length(modules_sym)){
		print(i)
		s.cell <- names(modules_sym)[i]
		cluster_C <- modules_sym[[s.cell]]
		results <- enricher(gene=cluster_C, pvalueCutoff=1, pAdjustMethod="BH", background_set, minGSSize=10, maxGSSize=1000, qvalueCutoff=1, term2gene) # Compute functional enrichment of cluster 'C' using a hyper-geometric test.
		
  		results_table <- as.data.frame(results) # Get enrichment results for current cluster 'C'.
  		results_table <- results_table[results_table$Count>3,] # Consider only annotations that cover at least 3 proteins from cluster.
		results_table <- results_table[results_table$p.adjust<=0.05,]
		if(dim(results_table)[1]==0) next
		results_table <- cbind(cell_type=rep(s.cell,dim(results_table)[1]),results_table)
		desc <- go2term(results_table$ID)
		(results_table$Description <- desc$Term)
		head(results_table)
		all_results <- rbind(all_results,results_table)
	}
	rownames(all_results) <- 1:dim(all_results)[1]
	res_list[[go_class]] <- all_results
}

saveRDS(res_list,file='rds/enrichment_results.TMI.rds')
#saveRDS(res_list,file='rds/enrichment_results.GLADIATOR.orig.rds')
#saveRDS(res_list,file='rds/enrichment_results.GLADIATOR.virt.rds')

################################################# REVIGO analysis#######################################################

results_files = c('enrichment_results.GLADIATOR.orig','enrichment_results.GLADIATOR.virt','enrichment_results.TMI')

for(f in results_files){
	print(f)
	
	res_list <- readRDS(file=paste('rds/',f,".rds",sep=""))
	cells <- unique(res_list$F$cell_type)
	REVIGOdir <- "REVIGO_results"
	dir.create(REVIGOdir)
	resultsdir <- file.path(REVIGOdir,f)
	dir.create(resultsdir)
	
	REVIGO.metrics.list <- c("BP","CC","MF")
	REVIGO.statistics.mat <- matrix(0,nrow=length(cells),ncol=length(REVIGO.metrics.list))
	rownames(REVIGO.statistics.mat) <- cells
	colnames(REVIGO.statistics.mat) <- REVIGO.metrics.list
	
	summary.list <- c("Min", "1st Q", "Median", "Mean", "3rd Q", "Max")
	REVIGO.summary.mat <- matrix(0,nrow=length(REVIGO.metrics.list),ncol=length(summary.list))
	rownames(REVIGO.summary.mat) <- REVIGO.metrics.list
	colnames(REVIGO.summary.mat) <- summary.list
	
	for(i in 1:length(cells)){
		print(i)
		
		s.cell <- cells[i]
		outputpath <- file.path(resultsdir,s.cell)
		dir.create(outputpath)
		fileNameInput <- paste(s.cell,"input.txt",sep="_")
		input_filepath <- file.path(outputpath,fileNameInput)
		sapply(res_list,function(x) summary(c(table(x$cell_type))))
		res_table <- lapply(res_list,function(x) x[which(!is.na(match(x$cell_type,s.cell))),])
		res_table_edited <- do.call(rbind,res_table)[,c('ID','p.adjust')]
		res_table_final <- filter(res_table_edited,p.adjust<0.01)
		write.table(res_table_final,file=input_filepath,sep='\t',col.names=F,row.names=F,quote=F)

		# Read user data from a file:
		fileName <- input_filepath
		userData <- readChar(fileName,file.info(fileName)$size)

		# Submit job to Revigo [speciesTaxon = Homo sapiens NCBI:txid9606]:
		httr::POST(
		  url = "http://revigo.irb.hr/StartJob.aspx",
		  body = list(
			cutoff = "0.4",
			valueType = "pvalue",
			speciesTaxon = "9606",
			measure = "SIMREL",
			goList = userData
		  ),
		  # application/x-www-form-urlencoded
		  encode = "form"
		) -> res
		dat <- httr::content(res, encoding = "UTF-8")
		jobid <- jsonlite::fromJSON(dat,bigint_as_char=TRUE)$jobid

		# Check job status:
		running <- "1"
		while (running != "0" ) {
			httr::POST(
			  url = "http://revigo.irb.hr/QueryJobStatus.aspx",
			  query = list( jobid = jobid )
			) -> res2
			dat2 <- httr::content(res2, encoding = "UTF-8")
			running <- jsonlite::fromJSON(dat2)$running
			Sys.sleep(1)
		}

		for(j in 1:3){
			# Fetch results:
			httr::POST(
				url = "http://revigo.irb.hr/ExportJob.aspx",
				query = list(
				jobid = jobid, 
				namespace = j,
				type = "csvtable"
				)
			) -> res3
			dat3 <- httr::content(res3, encoding = "UTF-8")

			# Write results to a file:
			dat3 <- stri_replace_all_fixed(dat3, "\r", "")
			if(j==1){
				GO = "_BP" # BIOLOGICAL_PROCESS
				column = "BP"
			}
			if(j==2){
				GO = "_CC" # CELLULAR_COMPONENT
				column = "CC"
			}
			if(j==3){
				GO = "_MF" # MOLECULAR_FUNCTION
				column = "MF"
			}
			fileNameOutput <- paste(s.cell,GO,".csv",sep="")
			output_filepath <- file.path(outputpath,fileNameOutput)
			cat(dat3, file=output_filepath, fill = FALSE)
			revigo_results <- read.csv(output_filepath)
			revigo_results_filtered <- filter(revigo_results,Eliminated==" False"&Dispensability<0.1)
			term_num <- nrow(revigo_results_filtered)
			REVIGO.statistics.mat[s.cell,column] <- term_num
		}
	}
	
	statfile_name <- paste("REVIGOstats_",f,".csv",sep="")
	write.csv(REVIGO.statistics.mat, statfile_name)
	for(ann in REVIGO.metrics.list){
		cur_results <- REVIGO.statistics.mat[,ann]
		curr_summ <- summary(cur_results)
		REVIGO.summary.mat[ann,] <- t(matrix(curr_summ))[1,]
	}
	summfile_name <- paste("REVIGOsummary_",f,".csv",sep="")
	write.csv(REVIGO.summary.mat, summfile_name)
	
}
