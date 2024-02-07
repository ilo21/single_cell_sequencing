library(Seurat)
# Define the output file path 
output_file1 <- "C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample1\\sample1_R.rds"
output_file2 <- "C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample2\\sample2_R.rds"
output_file3 <- "C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample3\\sample3_R.rds"
output_file4 <- "C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample4\\sample4_R.rds"
output_file5 <- "C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample5\\sample5_R.rds"

###########################################
# SAMPLE1
received_data1 <- readRDS("E:\\SCS\\Plate1_February_23\\SS3_22_291\\zUMIs_output\\expression\\SS3_22_291.dgecounts.rds")
new_data1 <- received_data1["umicount"]
exon1 <- new_data1$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix1 <- as.data.frame(exon1)
# Save data frame to an RDS file
saveRDS(my_matrix1, file = output_file1)
# read attached txt file with matched gene ids and gene names (55488 rows)
genes1 <- read.table("E:\\SCS\\Plate1_February_23\\SS3_22_291\\zUMIs_output\\expression\\SS3_22_291.gene_names.txt")
gene_id_row_names1 <- rownames(my_matrix1)
all_gene_ids1 <- genes1[[1]]
gene_names_row_names1 <- list() # new list with matched gene names (if found in the txt file)
cnt <- 0 # check how many ids didn't have gene names
for (i in 1:length(gene_id_row_names1)) {
  id <- gene_id_row_names1[i]
  if (!(id %in% all_gene_ids1)) {
    #cat(paste(x, "is not in the list.\n"))
    cnt <- cnt+1
    gene_names_row_names1[i] <- id
  }
  else {
    name <- genes1[genes1$V1==id,][2]
    gene_names_row_names1[i] <- name
  }
}
# Convert the result to a simple character vector without list structures
gene_names_row_names1 <- unlist(gene_names_row_names1)


# Generate unique row names for gene_names_row_names
unique_gene_names_row_names1 <- make.unique(gene_names_row_names1, sep = "-")

# Set the unique character vector as new row names for my_matrix1
rownames(my_matrix1) <- unique_gene_names_row_names1


# READY TO CREATE SEURAT OBJECT
sample1 <- CreateSeuratObject(counts = my_matrix1)

###########################################
# SAMPLE2
received_data2 <- readRDS("E:\\SCS\\Plate2_mAY_23\\SS3_23_049\\zUMIs_output\\expression\\SS3_23_049.dgecounts.rds")
new_data2 <- received_data2["umicount"]
exon2 <- new_data2$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix2 <- as.data.frame(exon2)

# Save data frame to an RDS file
saveRDS(my_matrix2, file = output_file2)

# read attached txt file with matched gene ids and gene names (55488 rows)
genes2 <- read.table("E:\\SCS\\Plate2_mAY_23\\SS3_23_049\\zUMIs_output\\expression\\SS3_23_049.gene_names.txt")
gene_id_row_names2 <- rownames(my_matrix2)
all_gene_ids2 <- genes2[[1]]
gene_names_row_names2 <- list() # new list with matched gene names (if found in the txt file)
cnt <- 0 # check how many ids didn't have gene names
for (i in 1:length(gene_id_row_names2)) {
  id <- gene_id_row_names2[i]
  if (!(id %in% all_gene_ids2)) {
    #cat(paste(x, "is not in the list.\n"))
    cnt <- cnt+1
    gene_names_row_names2[i] <- id
  }
  else {
    name <- genes2[genes2$V1==id,][2]
    gene_names_row_names2[i] <- name
  }
}
# Convert the result to a simple character vector without list structures
gene_names_row_names2 <- unlist(gene_names_row_names2)


# Generate unique row names for gene_names_row_names
unique_gene_names_row_names2 <- make.unique(gene_names_row_names2, sep = "-")

# Set the unique character vector as new row names for my_matrix1
rownames(my_matrix2) <- unique_gene_names_row_names2


# READY TO CREATE SEURAT OBJECT
sample2 <- CreateSeuratObject(counts = my_matrix2)

############################################################################
# SAMPLE3
received_data3 <- readRDS("E:\\SCS\\Plate3_August_23\\SS3_23_127\\zUMIs_output\\expression\\SS3_23_127.dgecounts.rds")
new_data3 <- received_data3["umicount"]
exon3 <- new_data3$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix3 <- as.data.frame(exon3)

# Save data frame to an RDS file
saveRDS(my_matrix3, file = output_file3)

# read attached txt file with matched gene ids and gene names (55488 rows)
genes3 <- read.table("E:\\SCS\\Plate3_August_23\\SS3_23_127\\zUMIs_output\\expression\\SS3_23_127.gene_names.txt")
gene_id_row_names3 <- rownames(my_matrix3)
all_gene_ids3 <- genes3[[1]]
gene_names_row_names3 <- list() # new list with matched gene names (if found in the txt file)
cnt <- 0 # check how many ids didn't have gene names
for (i in 1:length(gene_id_row_names3)) {
  id <- gene_id_row_names3[i]
  if (!(id %in% all_gene_ids3)) {
    #cat(paste(x, "is not in the list.\n"))
    cnt <- cnt+1
    gene_names_row_names3[i] <- id
  }
  else {
    name <- genes3[genes3$V1==id,][2]
    gene_names_row_names3[i] <- name
  }
}
# Convert the result to a simple character vector without list structures
gene_names_row_names3 <- unlist(gene_names_row_names3)


# Generate unique row names for gene_names_row_names
unique_gene_names_row_names3 <- make.unique(gene_names_row_names3, sep = "-")

# Set the unique character vector as new row names for my_matrix1
rownames(my_matrix3) <- unique_gene_names_row_names3


# READY TO CREATE SEURAT OBJECT
sample3 <- CreateSeuratObject(counts = my_matrix3)

###########################################
# SAMPLE4

received_data4 <- readRDS("C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample4\\SS3_23_193.dgecounts.rds")
new_data4 <- received_data4["umicount"]
exon4 <- new_data4$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix4 <- as.data.frame(exon4)
# Save data frame to an RDS file
saveRDS(my_matrix4, file = output_file4)
# read attached txt file with matched gene ids and gene names (55488 rows)
genes4 <- read.table("C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample4\\SS3_23_193.gene_names.txt")
gene_id_row_names4 <- rownames(my_matrix4)
all_gene_ids4 <- genes4[[1]]
gene_names_row_names4 <- list() # new list with matched gene names (if found in the txt file)
cnt <- 0 # check how many ids didn't have gene names
for (i in 1:length(gene_id_row_names4)) {
  id <- gene_id_row_names4[i]
  if (!(id %in% all_gene_ids4)) {
    #cat(paste(x, "is not in the list.\n"))
    cnt <- cnt+1
    gene_names_row_names4[i] <- id
  }
  else {
    name <- genes4[genes4$V1==id,][2]
    gene_names_row_names4[i] <- name
  }
}
# Convert the result to a simple character vector without list structures
gene_names_row_names4 <- unlist(gene_names_row_names4)


# Generate unique row names for gene_names_row_names
unique_gene_names_row_names4 <- make.unique(gene_names_row_names4, sep = "-")

# Set the unique character vector as new row names for my_matrix1
rownames(my_matrix4) <- unique_gene_names_row_names4


# READY TO CREATE SEURAT OBJECT
sample4 <- CreateSeuratObject(counts = my_matrix4)

###########################################
# SAMPLE5
received_data5 <- readRDS("C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample5\\SS3_23_195.dgecounts.rds")
new_data5 <- received_data5["umicount"]
exon5 <- new_data5$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix5 <- as.data.frame(exon5)
# Save data frame to an RDS file
saveRDS(my_matrix5, file = output_file5)
# read attached txt file with matched gene ids and gene names (55488 rows)
genes5 <- read.table("C:\\Users\\ilosz01\\OneDrive - Linköpings universitet\\MarcinLab\\SingleCellSequencing\\scs_analysis\\my_data\\sample5\\SS3_23_195.gene_names.txt")
gene_id_row_names5 <- rownames(my_matrix5)
all_gene_ids5 <- genes5[[1]]
gene_names_row_names5 <- list() # new list with matched gene names (if found in the txt file)
cnt <- 0 # check how many ids didn't have gene names
for (i in 1:length(gene_id_row_names5)) {
  id <- gene_id_row_names5[i]
  if (!(id %in% all_gene_ids5)) {
    #cat(paste(x, "is not in the list.\n"))
    cnt <- cnt+1
    gene_names_row_names5[i] <- id
  }
  else {
    name <- genes5[genes5$V1==id,][2]
    gene_names_row_names5[i] <- name
  }
}
# Convert the result to a simple character vector without list structures
gene_names_row_names5 <- unlist(gene_names_row_names5)


# Generate unique row names for gene_names_row_names
unique_gene_names_row_names5 <- make.unique(gene_names_row_names5, sep = "-")

# Set the unique character vector as new row names for my_matrix1
rownames(my_matrix5) <- unique_gene_names_row_names5


# READY TO CREATE SEURAT OBJECT
sample5 <- CreateSeuratObject(counts = my_matrix5)



#################################################################################################

# MERGE DATA FROM 3 SAMPLES
adata <- merge(sample1, y = c(sample2, sample3), add.cell.ids = c("sample1", "sample2", "sample3"), project = "sample123")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
adata[["percent.mt"]] <- PercentageFeatureSet(adata, pattern = "^mt-")
# Show QC metrics for the first 5 cells
# head(sample1@meta.data, 5)
# Visualize QC metrics as a violin plot
VlnPlot(adata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


################################################################################################################
cnt <- 0 # check how many genes are missing
marker_genes = c('Trpm3','Piezo2','Trpm2','Smr2','Sstr2','Bmpr1b','Trpm8','Trpv1','Piezo2','Piezo1','Nppb',
                'Sst','Pvalb','Prokr2','Mrgprd','Mrgpra3','Cd34',
                'Th','Trpa1','Ntrk3','Ntrk2','Ntrk1','Ret','Tac1','Calca','Calcb','Nefh',
                'S100b','Scn10a','Slc17a8','Atf3','Pou4f3','Calb1','Calb2','Avil','Asic3',
                'Asic2','Asic1','Pou6f2','Avpr1a','Pou4f2','Sox10','Casq2','Chrna7','Chrna3',
                'P2rx3','Gfra2','Ldhb','Necab2','Spp1','Adm','Hpse')

for (gene in marker_genes) {
  if (!(gene %in% gene_names_row_names)) {
    #cat(paste(x, "is not in the list.\n"))
    cnt <- cnt+1
  }
}
# all are there :)

