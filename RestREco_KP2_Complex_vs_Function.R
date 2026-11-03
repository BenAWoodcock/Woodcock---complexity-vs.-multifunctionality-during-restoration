# Restoration of ecological complexity drives system multifunctionality

# Assessing the relationships between complexity and function.
# Developed under R 4.5.1

library(pvclust)
library(DHARMa)
library(dplyr)


# Specify working directory 
setwd("P:\\07583 RestREco (Restoring Resilient Ecosystems)\\PAPERS\\Ecosystem complexity -  processes paper")
getwd()
list.files()

# files provided for raw data
site_info<-read.csv("Site_info.csv")
site_info_grass <- site_info[1:18, ]
site_info_wood <- site_info[19:36, ]

process_metrics<-read.csv("Process_metrics.csv")
process_metrics_grass <- process_metrics[1:18, ]
process_metrics_wood <- process_metrics[19:36, ]
# Note some NA values are present

complexity_metrics<-read.csv("Complexity_metrics.csv")
complexity_metrics_grass <- complexity_metrics[1:18, ]
complexity_metrics_wood <- complexity_metrics[19:36, ]

par(mfrow = c(2, 2))

#=============================================================================================================
#  Approaches for defining multifunctionlaity follow those described in 
# Manning, P., van der Plas, F., Soliveres, S., Allan, E., Maestre, F.T., Mace, G., Whittingham, M.J. & Fischer, M. (2018) Redefining ecosystem multifunctionality. Nature Ecology & Evolution, 2, 427-436.
# Byrnes, J.E.K., Gamfeldt, L., Isbell, F., Lefcheck, J.S., Griffin, J.N., Hector, A., Cardinale, B.J., Hooper, D.U., Dee, L.E. & Emmett Duffy, J. (2014) Investigating the relationship between biodiversity and ecosystem multifunctionality: challenges and solutions. Methods in Ecology and Evolution, 5, 111-124.

#Increasing measures of all function metrics are considered to indicate higher process rates

#----------------------------------------------------------------------------------------------------
# grassland derivation of multifunctionality index based on thresholds of 0.25, 0.5 and 0.75
working_df<-process_metrics_grass
working_df <- working_df [, !apply(is.na(working_df ), 2, all)] #  remove columns with just NA values
head(working_df,2)
working_df[is.na(working_df)] <- 0   #  missing values are assumed to have a function of 0
sapply(working_df, is.numeric)  #  check what columns are numeric
# Convert to matrix
rownames_for_matrix <- as.character(working_df[,1])
working_df_no_site_names <- working_df[, -1]
working_matrix <- as.matrix(working_df_no_site_names)
rownames(working_matrix ) <- rownames_for_matrix

# Wards clustering to identify weighting for the multifunctionality metric
working_matrix_scaled <- scale(working_matrix)  #  standardize matrix values for each metric across the sites 
# identify statistically significant clusters
library(pvclust)
set.seed(123)
pv <- pvclust(working_matrix_scaled, method.hclust="ward.D2", method.dist="euclidean", nboot=1000)
plot(pv, main = "Grassland function metrics cluster dendrogram")
pvrect(pv, alpha=0.95)  # Highlight clusters with approx. unbiased p-value >= 0.95
clusters_95 <- pvpick(pv, alpha=0.95)
# identify clusters that are significant and non-clustered elements
all_columns <- colnames(working_matrix) # names of each column
strong_clusters <- clusters_95$clusters# strong_clusters from pvpick()
clustered_columns <- unique(unlist(strong_clusters))# Flatten list of clustered columns
unclustered_columns <- setdiff(all_columns, clustered_columns)# Identify columns NOT clustered

# Create a data frame listing clusters and their members
cluster_df <- do.call(rbind, lapply(seq_along(strong_clusters), function(i) {
  data.frame(
    Cluster = paste0("Cluster ", i),
    Member = strong_clusters[[i]],
    stringsAsFactors = FALSE
  )
}))
# Add unclustered columns with unique sequential names
if(length(unclustered_columns) > 0) {
  unclustered_df <- data.frame(
    Cluster = paste0("Unclustered ", seq_along(unclustered_columns)),
    Member = unclustered_columns,
    stringsAsFactors = FALSE
  )
  
  # Combine clustered and unclustered
  cluster_df <- rbind(cluster_df, unclustered_df)
}

#------------------------------------------------------------
#  Deriving the threshold multifunctionality metric for 75th percentile
# In the dataframe identify any value greater than the 75th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_75thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.75, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of multifunctionality
site_info_grass$MF_un_75th <- rowSums(working_df_binary_75thentile)

#the following produces a wght  measure of Multifunctionality, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_75thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold multifunctionality metric
site_info_grass$MF_wght_75th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold multifunctionality metric for 50th percentile
# In the dataframe identify any value greater than the 50th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_50thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.50, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of multifunctionality
site_info_grass$MF_un_50th <- rowSums(working_df_binary_50thentile)

#the following produces a wght  measure of Multifunctionality, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_50thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold multifunctionality metric
site_info_grass$MF_wght_50th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold multifunctionality metric for 25th percentile
# In the dataframe identify any value greater than the 25th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_25thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.25, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of multifunctionality
site_info_grass$MF_un_25th <- rowSums(working_df_binary_25thentile)

#the following produces a wght  measure of Multifunctionality, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_25thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold multifunctionality metric
site_info_grass$MF_wght_25th<- rowSums(cluster_means_df)


#========================================================================================
# woodland multifunctionality index based on thresholds of 0.25, 0.5 and 0.75
working_df<-process_metrics_wood
working_df <- working_df [, !apply(is.na(working_df ), 2, all)] #  remove columns with just NA values
head(working_df,2)
working_df[is.na(working_df)] <- 0   #  missing values are assumed to have a function of 0
sapply(working_df, is.numeric)  #  check what columns are numeric
# Convert to matrix
rownames_for_matrix <- as.character(working_df[,1])
working_df_no_site_names <- working_df[, -1]
working_matrix <- as.matrix(working_df_no_site_names)
rownames(working_matrix ) <- rownames_for_matrix

# Wards clustering to identify weighting for the multifunctionality metric
working_matrix_scaled <- scale(working_matrix)  #  standardize matrix values for each metric across the sites 
# identify statistically significant clusters
library(pvclust)
set.seed(123)
pv <- pvclust(working_matrix_scaled, method.hclust="ward.D2", method.dist="euclidean", nboot=1000)
plot(pv,main = "Woodland function metrics cluster dendrogram")
pvrect(pv, alpha=0.95)  # Highlight clusters with approx. unbiased p-value >= 0.95
clusters_95 <- pvpick(pv, alpha=0.95)
# identify clusters that are significant and non-clustered elements
all_columns <- colnames(working_matrix) # names of each column
strong_clusters <- clusters_95$clusters# strong_clusters from pvpick()
clustered_columns <- unique(unlist(strong_clusters))# Flatten list of clustered columns
unclustered_columns <- setdiff(all_columns, clustered_columns)# Identify columns NOT clustered

# Create a data frame listing clusters and their members
cluster_df <- do.call(rbind, lapply(seq_along(strong_clusters), function(i) {
  data.frame(
    Cluster = paste0("Cluster ", i),
    Member = strong_clusters[[i]],
    stringsAsFactors = FALSE
  )
}))
# Add unclustered columns with unique sequential names
if(length(unclustered_columns) > 0) {
  unclustered_df <- data.frame(
    Cluster = paste0("Unclustered ", seq_along(unclustered_columns)),
    Member = unclustered_columns,
    stringsAsFactors = FALSE
  )
  
  # Combine clustered and unclustered
  cluster_df <- rbind(cluster_df, unclustered_df)
}

#------------------------------------------------------------
#  Deriving the threshold multifunctionality metric for 75th percentile
# In the dataframe identify any value greater than the 75th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_75thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.75, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of multifunctionality
site_info_wood$MF_un_75th <- rowSums(working_df_binary_75thentile)

#the following produces a wght  measure of Multifunctionality, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_75thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold multifunctionality metric
site_info_wood$MF_wght_75th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold multifunctionality metric for 50th percentile
# In the dataframe identify any value greater than the 50th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_50thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.50, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of multifunctionality
site_info_wood$MF_un_50th <- rowSums(working_df_binary_50thentile)

#the following produces a wght  measure of Multifunctionality, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_50thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold multifunctionality metric
site_info_wood$MF_wght_50th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold multifunctionality metric for 25th percentile
# In the dataframe identify any value greater than the 25th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_25thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.25, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of multifunctionality
site_info_wood$MF_un_25th <- rowSums(working_df_binary_25thentile)

#the following produces a wght  measure of Multifunctionality, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_25thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold multifunctionality metric
site_info_wood$MF_wght_25th<- rowSums(cluster_means_df)




#===================================================================================================================
# deriving metacomplexity metrics

#=============================================================================================================
#  Approaches for defining complexity will  follow two different approaches.
#  IN both cases they will broadly apply the rpinciples applied in 
# Manning, P., van der Plas, F., Soliveres, S., Allan, E., Maestre, F.T., Mace, G., Whittingham, M.J. & Fischer, M. (2018) Redefining ecosystem multifunctionality. Nature Ecology & Evolution, 2, 427-436.
#  This is to say the complexity metrics will be agglomeration clustered.  We will then apply the threshold based approach
# with wght averaging as described in 
# Byrnes, J.E.K., Gamfeldt, L., Isbell, F., Lefcheck, J.S., Griffin, J.N., Hector, A., Cardinale, B.J., Hooper, D.U., Dee, L.E. & Emmett Duffy, J. (2014) Investigating the relationship between biodiversity and ecosystem multifunctionality: challenges and solutions. Methods in Ecology and Evolution, 5, 111-124.
#  we will then apply the wght averaging to standardised values (not thresholds) as well for the clustered
# complexity metrics.

#Increasing values of all complexity  metrics are considered to indicate higher process rates

#----------------------------------------------------------------------------------------------------
#  Threshold based metacomplexity indexes for 25th, 50th and 75th percentiles
# grassland derivation of complexity index based on thresholds of 0.25, 0.5 and 0.75
working_df<-complexity_metrics_grass
working_df <- working_df [, !apply(is.na(working_df ), 2, all)] #  remove columns with just NA values
head(working_df,2)
working_df[is.na(working_df)] <- 0   #  missing values are assumed to have a function of 0
sapply(working_df, is.numeric)  #  check what columns are numeric
# Convert to matrix
rownames_for_matrix <- as.character(working_df[,1])
working_df_no_site_names <- working_df[, -1]
working_matrix <- as.matrix(working_df_no_site_names)
rownames(working_matrix ) <- rownames_for_matrix

# Wards clustering to identify weighting for the metacomplexity metric
working_matrix_scaled <- scale(working_matrix)  #  standardize matrix values for each metric across the sites 
# identify statistically significant clusters
library(pvclust)
set.seed(123)
pv <- pvclust(working_matrix_scaled, method.hclust="ward.D2", method.dist="euclidean", nboot=1000)
plot(pv,main = "Grassland complexity cluster dendrogram")
pvrect(pv, alpha=0.95)  # Highlight clusters with approx. unbiased p-value >= 0.95
clusters_95 <- pvpick(pv, alpha=0.95)
# identify clusters that are significant and non-clustered elements
all_columns <- colnames(working_matrix) # names of each column
strong_clusters <- clusters_95$clusters# strong_clusters from pvpick()
clustered_columns <- unique(unlist(strong_clusters))# Flatten list of clustered columns
unclustered_columns <- setdiff(all_columns, clustered_columns)# Identify columns NOT clustered

# Create a data frame listing clusters and their members
cluster_df <- do.call(rbind, lapply(seq_along(strong_clusters), function(i) {
  data.frame(
    Cluster = paste0("Cluster ", i),
    Member = strong_clusters[[i]],
    stringsAsFactors = FALSE
  )
}))
# Add unclustered columns with unique sequential names
if(length(unclustered_columns) > 0) {
  unclustered_df <- data.frame(
    Cluster = paste0("Unclustered ", seq_along(unclustered_columns)),
    Member = unclustered_columns,
    stringsAsFactors = FALSE
  )
  
  # Combine clustered and unclustered
  cluster_df <- rbind(cluster_df, unclustered_df)
}

#------------------------------------------------------------
#  Deriving the threshold metacomplexity metric for 75th percentile
# In the dataframe identify any value greater than the 75th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_75thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.75, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of metacomplexity
site_info_grass$CX_Th_un_75th <- rowSums(working_df_binary_75thentile)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_75thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_grass$CX_Th_wght_75th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold metacomplexity metric for 50th percentile
# In the dataframe identify any value greater than the 50th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_50thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.50, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of metacomplexity
site_info_grass$CX_Th_un_50th <- rowSums(working_df_binary_50thentile)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_50thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_grass$CX_Th_wght_50th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold metacomplexity metric for 25th percentile
# In the dataframe identify any value greater than the 25th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_25thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.25, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of metacomplexity
site_info_grass$CX_Th_un_25th <- rowSums(working_df_binary_25thentile)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_25thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_grass$CX_Th_wght_25th<- rowSums(cluster_means_df)

#----------------------------------------------------------------------------------------------------

#  Deriving the averaged metacomplexity metric for grassland standardises values

site_info_grass$CX_StdV_un<- rowSums(working_matrix_scaled)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75th percentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_matrix_scaled[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_grass$CX_StdV_wght<- rowSums(cluster_means_df)



#######################################################################################
#========================================================================================
# woodland metacomplexity index based on thresholds of 0.25, 0.5 and 0.75
working_df<-complexity_metrics_wood
working_df <- working_df [, !apply(is.na(working_df ), 2, all)] #  remove columns with just NA values
head(working_df,2)
working_df[is.na(working_df)] <- 0   #  missing values are assumed to have a function of 0
sapply(working_df, is.numeric)  #  check what columns are numeric
# Convert to matrix
rownames_for_matrix <- as.character(working_df[,1])
working_df_no_site_names <- working_df[, -1]
working_matrix <- as.matrix(working_df_no_site_names)
rownames(working_matrix ) <- rownames_for_matrix

# Wards clustering to identify weighting for the metacomplexity metric
working_matrix_scaled <- scale(working_matrix)  #  standardize matrix values for each metric across the sites 
# identify statistically significant clusters
library(pvclust)
set.seed(123)
pv <- pvclust(working_matrix_scaled, method.hclust="ward.D2", method.dist="euclidean", nboot=1000)
plot(pv, main = "Woodland complexity cluster dendrogram")
pvrect(pv, alpha=0.95)  # Highlight clusters with approx. unbiased p-value >= 0.95
clusters_95 <- pvpick(pv, alpha=0.95)
# identify clusters that are significant and non-clustered elements
all_columns <- colnames(working_matrix) # names of each column
strong_clusters <- clusters_95$clusters# strong_clusters from pvpick()
clustered_columns <- unique(unlist(strong_clusters))# Flatten list of clustered columns
unclustered_columns <- setdiff(all_columns, clustered_columns)# Identify columns NOT clustered

# Create a data frame listing clusters and their members
cluster_df <- do.call(rbind, lapply(seq_along(strong_clusters), function(i) {
  data.frame(
    Cluster = paste0("Cluster ", i),
    Member = strong_clusters[[i]],
    stringsAsFactors = FALSE
  )
}))
# Add unclustered columns with unique sequential names
if(length(unclustered_columns) > 0) {
  unclustered_df <- data.frame(
    Cluster = paste0("Unclustered ", seq_along(unclustered_columns)),
    Member = unclustered_columns,
    stringsAsFactors = FALSE
  )
  
  # Combine clustered and unclustered
  cluster_df <- rbind(cluster_df, unclustered_df)
}

#------------------------------------------------------------
#  Deriving the threshold metacomplexity metric for 75th percentile
# In the dataframe identify any value greater than the 75th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_75thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.75, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of metacomplexity
site_info_wood$CX_Th_un_75th <- rowSums(working_df_binary_75thentile)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_75thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_wood$CX_Th_wght_75th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold metacomplexity metric for 50th percentile
# In the dataframe identify any value greater than the 50th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_50thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.50, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of metacomplexity
site_info_wood$CX_Th_un_50th <- rowSums(working_df_binary_50thentile)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_50thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_wood$CX_Th_wght_50th<- rowSums(cluster_means_df)

#------------------------------------------------------------
#  Deriving the threshold metacomplexity metric for 25th percentile
# In the dataframe identify any value greater than the 25th percentile as 1 (less than as zero)

rownames<- as.character(working_df[,1])
working_df_binary_25thentile <- as.data.frame(lapply(working_df[,-1], function(col) {     # note excluded col1 site names
  thresh <- quantile(col, 0.25, na.rm = TRUE)
  as.integer(col > thresh)
}))
# now combine these values across the rows to get an un-wght measure of metacomplexity
site_info_wood$CX_Th_un_25th <- rowSums(working_df_binary_25thentile)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75thentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_df_binary_25thentile[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_wood$CX_Th_wght_25th<- rowSums(cluster_means_df)

#----------------------------------------------------------------------------------------------------

#  Deriving the averaged metacomplexity metric for woodland standardises values

site_info_wood$CX_StdV_un<- rowSums(working_matrix_scaled)

#the following produces a wght  measure of metacomplexity, so that each cluster has a total weight
# of 1  - i.e. if there are 4 elements in a cluster these are divided by 4.
clusters <- unique(cluster_df$Cluster)
#For each cluster, calculate the row mean of corresponding columns in working_df_binary_75th percentile
cluster_means <- lapply(clusters, function(clust) {
  members <- cluster_df$Member[cluster_df$Cluster == clust]  # Get members of this cluster
  rowMeans(working_matrix_scaled[, members, drop = FALSE])# Calculate row sums of members columns, then divide by number of members
})
# Combine cluster means into a data frame
cluster_means_df <- do.call(cbind, cluster_means)
colnames(cluster_means_df) <- clusters
# Sum these cluster means row-wise to get the wght sum representing the threshold metacomplexity metric
site_info_wood$CX_StdV_wght<- rowSums(cluster_means_df)



###################################################################################################################
#  Finalizing data set for analysis
# we can look for relationships between 

#-----------------------------------------------------------------
#  Standardize each of the new metrics so that we can combine the woodland and grassland sites within a single analysis
# for select derived multifunctionality (MF_) and metacomplexity (CX_) metrics.

# Grasslands
selected <- site_info_grass [, which(names(site_info_grass) == "MF_un_75th") : which(names(site_info_grass) == "CX_StdV_wght")]
# Min-max scale (standardize to 0–1)
scaled <- as.data.frame(lapply(selected, function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}))
# Rename columns with "Std_" prefix
names(scaled) <- paste0("Std_", names(scaled))
# Combine with original data frame
site_info_grass<- cbind(site_info_grass, scaled)

#Woodlands
selected <- site_info_wood [, which(names(site_info_wood) == "MF_un_75th") : which(names(site_info_wood) == "CX_StdV_wght")]
# Min-max scale (standardize to 0–1)
scaled <- as.data.frame(lapply(selected, function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}))
# Rename columns with "Std_" prefix
names(scaled) <- paste0("Std_", names(scaled))
# Combine with original data frame
site_info_wood <- cbind(site_info_wood, scaled)

#Combine woodland and grassland data, noting for analysis we need to use the standardized (0-1) multifunctionality and 
# metacomplexity metrics
site_info<-rbind(site_info_grass, site_info_wood)


##==============================================================================================================
#  Main paper analyses 
#===============================================================================================================
library(broom)

# Effect of plant species richness on multifunctionality

# Main paper fig 1


#----------------------------------------------------------------------------
# Combined site responses
# The response of low multifunctionality (Std_MF_wght_25th)

mod_ful <- lm(Std_MF_wght_25th ~ Plant.species.incl.sp + Habitat +Habitat*Plant.species.incl.sp , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_25th ~ Plant.species.incl.sp + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_25th ~ Plant.species.incl.sp , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)

ggplot(site_info, aes(x = Plant.species.incl.sp,  y = Std_MF_wght_25th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Plant species richness",
       y = "Multifunctionality (>25th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Plant.species.incl.sp),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 15,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )

# The response of medium multifunctionality (Std_MF_wght_50th)


mod_ful <- lm(Std_MF_wght_50th ~ Plant.species.incl.sp + Habitat +Habitat*Plant.species.incl.sp , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_50th ~ Plant.species.incl.sp + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_50th ~ Plant.species.incl.sp , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)


ggplot(site_info, aes(x = Plant.species.incl.sp,  y = Std_MF_wght_50th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Plant species richness",
       y = "Multifunctionality (>50th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Plant.species.incl.sp),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² < 0.01"),
    hjust = 0, vjust = 1.5,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )


# The response of high multifunctionality (Std_MF_wght_75th)

mod_ful <- lm(Std_MF_wght_75th ~ Plant.species.incl.sp + Habitat +Habitat*Plant.species.incl.sp , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_75th ~ Plant.species.incl.sp + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_75th ~ Plant.species.incl.sp , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)

ggplot(site_info, aes(x = Plant.species.incl.sp,  y = Std_MF_wght_75th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Plant species richness",
       y = "Multifunctionality (>75th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Plant.species.incl.sp),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )

#------------------------------------------------
# Grass only sites
# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Plant.species.incl.sp, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)

ggplot(site_info_grass, aes(x = Plant.species.incl.sp, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Plant species richness", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Plant.species.incl.sp),
    y = max(site_info_grass$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 2, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# The response of medium multifunctionality (Std_MF_wght_50th)

mod <- lm(Std_MF_wght_50th ~ Plant.species.incl.sp, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Plant.species.incl.sp, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Plant species richness", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Plant.species.incl.sp),
    y = max(site_info_grass$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of high multifunctionality (Std_MF_wght_75th)

mod <- lm(Std_MF_wght_75th ~ Plant.species.incl.sp, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Plant.species.incl.sp, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Plant species richness", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Plant.species.incl.sp),
    y = max(site_info_grass$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

#-----------------------------------------------------------------------------
#   Woodland only sites

# The response of low multifunctionality (Std_MF_wght_25th)


mod <- lm(Std_MF_wght_25th ~ Plant.species.incl.sp, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)

ggplot(site_info_wood, aes(x = Plant.species.incl.sp, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Plant species richness", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Plant.species.incl.sp),
    y = max(site_info_wood$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 0, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# The response of medium multifunctionality (Std_MF_wght_50th)

mod <- lm(Std_MF_wght_50th ~ Plant.species.incl.sp, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Plant.species.incl.sp, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Plant species richness", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Plant.species.incl.sp),
    y = max(site_info_wood$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 2, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of high multifunctionality (Std_MF_wght_75th)

mod <- lm(Std_MF_wght_75th ~ Plant.species.incl.sp, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Plant.species.incl.sp, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Plant species richness", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Plant.species.incl.sp),
    y = max(site_info_wood$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


#====================================================================================================
#  Effect fo system complexity on multifunctionality
#====================================================================================================

#--------------------------------------------------------------------------------
# Complexity at 75th percentile threshold (high complexity)
# Main paper fig 2
#------------------------------------------------------------------------------
# combined sites
# The response of low multifunctionality (Std_MF_wght_25th)

mod_ful <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_75th + Habitat +Habitat*Std_CX_Th_wght_75th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_75th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_75th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info)

anova(mod_2, mod_null)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)

#DHARMAa
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


ggplot(site_info, aes(x = Std_CX_Th_wght_75th,  y = Std_MF_wght_25th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>75th percentile)",
       y = "Multifunctionality (>25th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_75th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 0,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
 #   legend.title = element_text(size = 20),
 #   legend.text = element_text(size = 18)
legend.position = "none"
  )


# The response of medium multifunctionality (Std_MF_wght_50th)
mod_ful <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_75th + Habitat +Habitat*Std_CX_Th_wght_75th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_75th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_75th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

anova(mod_2, mod_null)

AIC(mod_2, mod_null)


ggplot(site_info, aes(x = Std_CX_Th_wght_75th,  y = Std_MF_wght_50th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>75th percentile)",
       y = "Multifunctionality (>50th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_75th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )


# The response of high multifunctionality (Std_MF_wght_75th)

mod_ful <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_75th + Habitat +Habitat*Std_CX_Th_wght_75th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_75th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_75th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2,  mod_null)




ggplot(site_info, aes(x = Std_CX_Th_wght_75th,  y = Std_MF_wght_75th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>75th percentile)",
       y = "Multifunctionality (>75th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_75th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )


#-----------------------------------------------------------------------------
# Grass only habitats
# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_75th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)

ggplot(site_info_grass, aes(x = Std_CX_Th_wght_75th, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>75th percentile)", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_75th),
    y = max(site_info_grass$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 0, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of medium multifunctionality (Std_MF_wght_50th)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_75th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Std_CX_Th_wght_75th, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>75th percentile)", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_75th),
    y = max(site_info_grass$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of high multifunctionality (Std_MF_wght_75th)



mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_75th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Std_CX_Th_wght_75th, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>75th percentile)", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_75th),
    y = max(site_info_grass$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


#-----------------------------------------------------------------------------
#   Woodland only sites

# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_75th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)

ggplot(site_info_wood, aes(x = Std_CX_Th_wght_75th, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>75th percentile)", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_75th),
    y = max(site_info_wood$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 0, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# The response of medium multifunctionality (Std_MF_wght_50th)


mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_75th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Std_CX_Th_wght_75th, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>75th percentile)", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_75th),
    y = max(site_info_wood$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of high multifunctionality (Std_MF_wght_75th)

mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_75th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions



mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Std_CX_Th_wght_75th, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>75th percentile)", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_75th),
    y = max(site_info_wood$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )




#################################################################################
#  the response of multifunctionality to sites achieving minimum levels of complexity (25th percentile)
# Fig S3 in the supplementary material
############################################################################

#------------------------------------------------------------------------------
# combined site responses

# The response of low multifunctionality (Std_MF_wght_25th)

mod_ful <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_25th + Habitat +Habitat*Std_CX_Th_wght_25th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_25th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_25th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)

ggplot(site_info, aes(x = Std_CX_Th_wght_25th,  y = Std_MF_wght_25th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>25th percentile)",
       y = "Multifunctionality (>25th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_25th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )

# The response of medium multifunctionality (Std_MF_wght_50th)


mod_ful <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_25th + Habitat +Habitat*Std_CX_Th_wght_25th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_25th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_25th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)


ggplot(site_info, aes(x = Std_CX_Th_wght_25th,  y = Std_MF_wght_50th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>25th percentile)",
       y = "Multifunctionality (>50th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_25th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )


# The response of high multifunctionality (Std_MF_wght_75th)

mod_ful <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_25th + Habitat +Habitat*Std_CX_Th_wght_25th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_25th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_25th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)




ggplot(site_info, aes(x = Std_CX_Th_wght_25th,  y = Std_MF_wght_75th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>25th percentile)",
       y = "Multifunctionality (>75th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_25th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )
#------------------------------------------------
# Grassland only sites
# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_25th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)

ggplot(site_info_grass, aes(x = Std_CX_Th_wght_25th, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>25th percentile)", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_25th),
    y = max(site_info_grass$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 0, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )



mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_25th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Std_CX_Th_wght_25th, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>25th percentile)", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_25th),
    y = max(site_info_grass$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )



# The response of medium multifunctionality (Std_MF_wght_50th)


mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_25th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Std_CX_Th_wght_25th, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>25th percentile)", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_25th),
    y = max(site_info_grass$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

#-----------------------------------------------------------------------------
#   woodland only sites
# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_25th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)

ggplot(site_info_wood, aes(x = Std_CX_Th_wght_25th, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>25th percentile)", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_25th),
    y = max(site_info_wood$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# The response of medium multifunctionality (Std_MF_wght_50th)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_25th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Std_CX_Th_wght_25th, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>25th percentile)", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_25th),
    y = max(site_info_wood$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )



# The response of high multifunctionality (Std_MF_wght_75th)


mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_25th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Std_CX_Th_wght_25th, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>25th percentile)", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_25th),
    y = max(site_info_wood$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )





##############################################################################

# The effect of sites achieving moderate Complexity (50th percentile) on multifunctionality
# Fig. S4 in the supplementary material

############################################################################

#------------------------------------------------------------------------------
# combined
# The response of low multifunctionality (Std_MF_wght_25th)

mod_ful <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_50th + Habitat +Habitat*Std_CX_Th_wght_50th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_50th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_50th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)





ggplot(site_info, aes(x = Std_CX_Th_wght_50th,  y = Std_MF_wght_25th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>50th percentile)",
       y = "Multifunctionality (>25th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_50th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 2,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )

# The response of medium multifunctionality (Std_MF_wght_50th)


mod_ful <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_50th + Habitat +Habitat*Std_CX_Th_wght_50th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_50th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_50th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)


ggplot(site_info, aes(x = Std_CX_Th_wght_50th,  y = Std_MF_wght_50th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>50th percentile)",
       y = "Multifunctionality (>50th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_50th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )

# The response of low multifunctionality (Std_MF_wght_75th)

mod_ful <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_50th + Habitat +Habitat*Std_CX_Th_wght_50th , data = site_info)
r2_val <- broom::glance(mod_ful)$r.squared
summary(mod_ful)

mod_2 <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_50th + Habitat , data = site_info)
r2_val <- broom::glance(mod_2)$r.squared
summary(mod_2)

mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_50th , data = site_info)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info)

anova(mod_2, mod_null)

AIC(mod_2, mod_null)




ggplot(site_info, aes(x = Std_CX_Th_wght_50th,  y = Std_MF_wght_75th)) +
  # points by habitat
  geom_point(aes(shape = Habitat), size = 6, stroke = 1.5) +
  # single global trend line
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  # manually assign shapes: 1 = open circle, 19 = solid circle
  scale_shape_manual(values = c("Grass" = 1,
                                "Wood" = 19)) +
  labs(x = "Complexity (>50th percentile)",
       y = "Multifunctionality (>75th percentile)",
       shape = "Habitat") +
  annotate(
    "text",
    x = min(site_info$Std_CX_Th_wght_50th),
    y = max(site_info$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3,
    size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    #   legend.title = element_text(size = 20),
    #   legend.text = element_text(size = 18)
    legend.position = "none"
  )

#-----------------------------------------------------------------------------
#   grass only sites
# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_50th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)

ggplot(site_info_grass, aes(x = Std_CX_Th_wght_50th, y = Std_MF_wght_25th)) +
  geom_point(size = 7, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>50th percentile)", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_50th),
    y = max(site_info_grass$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 2, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# The response of medium multifunctionality (Std_MF_wght_50th)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_50th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Std_CX_Th_wght_50th, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>50th percentile)", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_50th),
    y = max(site_info_grass$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of low multifunctionality (Std_MF_wght_75th)

mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_50th, data = site_info_grass)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_grass)
AIC(mod_null)
AIC(mod)


ggplot(site_info_grass, aes(x = Std_CX_Th_wght_50th, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 1) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>50th percentile)", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_grass$Std_CX_Th_wght_50th),
    y = max(site_info_grass$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


#-----------------------------------------------------------------------------
#   woodland only sites
# The response of low multifunctionality (Std_MF_wght_25th)

mod <- lm(Std_MF_wght_25th ~ Std_CX_Th_wght_50th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_25th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)

ggplot(site_info_wood, aes(x = Std_CX_Th_wght_50th, y = Std_MF_wght_25th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>50th percentile)", y = "Multifunctionality (>25th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_50th),
    y = max(site_info_wood$Std_MF_wght_25th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 0, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )

# The response of medium multifunctionality (Std_MF_wght_50th)

mod <- lm(Std_MF_wght_50th ~ Std_CX_Th_wght_50th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions

mod_null <- lm(Std_MF_wght_50th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Std_CX_Th_wght_50th, y = Std_MF_wght_50th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>50th percentile)", y = "Multifunctionality (>50th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_50th),
    y = max(site_info_wood$Std_MF_wght_50th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 1, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


# The response of high multifunctionality (Std_MF_wght_75th)



mod <- lm(Std_MF_wght_75th ~ Std_CX_Th_wght_50th, data = site_info_wood)
r2_val <- broom::glance(mod)$r.squared
summary(mod)
simulationOutput<-simulateResiduals(fittedModel=mod, plot=F)
plot(simulationOutput)   # meets model assumptions


mod_null <- lm(Std_MF_wght_75th ~ 1, data = site_info_wood)
AIC(mod_null)
AIC(mod)


ggplot(site_info_wood, aes(x = Std_CX_Th_wght_50th, y = Std_MF_wght_75th)) +
  geom_point(size = 6, shape = 19) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", fill = "lightblue") +
  labs(x = "Complexity (>50th percentile)", y = "Multifunctionality (>75th percentile)") +
  annotate(
    "text",
    x = min(site_info_wood$Std_CX_Th_wght_50th),
    y = max(site_info_wood$Std_MF_wght_75th),
    label = paste0("R² = ", round(r2_val, 3)),
    hjust = 0, vjust = 3, size = 15
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 30, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 30, face = "bold", margin = margin(r = 10)),
    axis.text = element_text(size = 25),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14)
  )


