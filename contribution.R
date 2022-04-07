load("fusion_v3.0.RData")
contribution=matrix(data = 0,nrow = 4,ncol = 7)
rownames(contribution)=c("Inside Clusters","Between Clusters","Total","Inside Ratio")
colnames(contribution)=c(1:7)
contribution[3,]=c(15120,6213,812,2335,194,39,40 )

#inner edges
for (i in 1:length(contribution_sta_group_list))
{
  for (j in 1:7)
  {
    num=which(names(contribution_sta_group_list[[i]])==j)
    if (length(num)>0)
    {
      contribution[1,j]=contribution[1,j]+contribution_sta_group_list[[i]][num]
    }
  }
}
#ratios
for (i in 1:7) 
{
  contribution[2,i]=contribution[3,i]-contribution[1,i]
  contribution[4,i]=contribution[1,i]/contribution[3,i]
}
write.table(contribution,file = "insider_cluster_contribution.txt",quote = F,sep = "\t")

#### contribution heatmap####
library(pheatmap) 
contribution_heatmap=matrix(data = 0,nrow = 7,ncol = length(contribution_sta_group_list))
for (i in 1:length(contribution_sta_group_list))
{
  for (j in 1:7)
  {
    num=which(names(contribution_sta_group_list[[i]])==j)
    if (length(num)>0)
    {
      contribution_heatmap[j,i]=contribution_sta_group_list[[i]][num]
    }
  }
}
temp=colSums(contribution_heatmap)
for(i in 1:nrow(contribution_heatmap))
{
  for (j in 1:ncol(contribution_heatmap))
  {
    contribution_heatmap[i,j]=contribution_heatmap[i,j]/temp[j]
  }
}
rownames(contribution_heatmap)=paste("Contribution",1:7,sep = "")
colnames(contribution_heatmap)=paste("Cluster",1:24,sep = "")
pheatmap(contribution_heatmap,cluster_cols=F,color = colorRampPalette(c("white", "firebrick3"))(50),cluster_rows=F,border_color="grey60",fontsize=10,show_rownames = T,show_colnames = T)
write.table(contribution_heatmap,file = "cluster_contribution.txt",quote = F,sep = "\t")

