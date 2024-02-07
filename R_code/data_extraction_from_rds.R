# Define the original files paths
# umi counts
received_data1_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample1/SS3_22_291.dgecounts.rds"
received_data2_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample2/SS3_23_049.dgecounts.rds"
received_data3_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample3/SS3_23_127.dgecounts.rds"
received_data4_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample4/SS3_23_193.dgecounts.rds"
received_data5_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample5/SS3_23_195.dgecounts.rds"
received_data6_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample6/SS3_23_325.dgecounts.rds"
received_data7_path <- "/home/ilosz01/Documents/single_cell_sequencing/data/sample7/SS3_23_327.dgecounts.rds"

##############################
# Define the output file path 
output_file1 <- "/home/ilosz01/Documents/single_cell_sequencing/sample1/counts_umi/sample1_R.rds"
output_file2 <- "/home/ilosz01/Documents/single_cell_sequencing/sample2/counts_umi/sample2_R.rds"
output_file3 <- "/home/ilosz01/Documents/single_cell_sequencing/sample3/counts_umi/sample3_R.rds"
output_file4 <- "/home/ilosz01/Documents/single_cell_sequencing/sample4/counts_umi/sample4_R.rds"
output_file5 <- "/home/ilosz01/Documents/single_cell_sequencing/sample5/counts_umi/sample5_R.rds"
output_file6 <- "/home/ilosz01/Documents/single_cell_sequencing/sample6/counts_umi/sample6_R.rds"
output_file7 <- "/home/ilosz01/Documents/single_cell_sequencing/sample7/counts_umi/sample7_R.rds"

###########################################
# SAMPLE1
received_data1 <- readRDS(received_data1_path)
new_data1 <- received_data1["umicount"]
exon1 <- new_data1$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix1 <- as.data.frame(as.matrix(exon1))
# Save data frame to an RDS file
saveRDS(my_matrix1, file = output_file1)
###########################################
# SAMPLE2
received_data2 <- readRDS(received_data2_path)
new_data2 <- received_data2["umicount"]
exon2 <- new_data2$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix2 <- as.data.frame(as.matrix(exon2))

# Save data frame to an RDS file
saveRDS(my_matrix2, file = output_file2)
############################################################################
# SAMPLE3
received_data3 <- readRDS(received_data3_path)
new_data3 <- received_data3["umicount"]
exon3 <- new_data3$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix3 <- as.data.frame(as.matrix(exon3))

# Save data frame to an RDS file
saveRDS(my_matrix3, file = output_file3)
###########################################
# SAMPLE4
received_data4 <- readRDS(received_data4_path)
new_data4 <- received_data4["umicount"]
exon4 <- new_data4$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix4 <- as.data.frame(as.matrix(exon4))
# Save data frame to an RDS file
saveRDS(my_matrix4, file = output_file4)
###########################################
# SAMPLE5
received_data5 <- readRDS(received_data5_path)
new_data5 <- received_data5["umicount"]
exon5 <- new_data5$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix5 <- as.data.frame(as.matrix(exon5))
# Save data frame to an RDS file
saveRDS(my_matrix5, file = output_file5)
###########################################
# SAMPLE6
received_data6 <- readRDS(received_data6_path)
new_data6 <- received_data6["umicount"]
exon6 <- new_data6$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix6 <- as.data.frame(as.matrix(exon6))
# Save data frame to an RDS file
saveRDS(my_matrix6, file = output_file6)
###########################################
# SAMPLE7
received_data7 <- readRDS(received_data7_path)
new_data7 <- received_data7["umicount"]
exon7 <- new_data7$umicount$exon$all
# my matrix with gene ids as rows and barcodes as columns (29247 rows)
my_matrix7 <- as.data.frame(as.matrix(exon7))
# Save data frame to an RDS file
saveRDS(my_matrix7, file = output_file7)
