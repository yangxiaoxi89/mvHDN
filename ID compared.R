load("edge.RData")
library(xlsx)

x=as.numeric(edge_fusion[,3])
edge_fusion_temp = edge_fusion[which(x>=quantile(x,probs = 0.9)),]
TreeNumbers_disease = read.xlsx(file = "Truelabel_TreeNumbers.xlsx",sheetIndex = 1)
TreeNumbers_disease = apply(TreeNumbers_disease, 2, as.character)
temp = matrix(data = NA, nrow = nrow(edge_fusion_temp), ncol = 3)
edge_fusion_temp = cbind(edge_fusion_temp, temp)
rm(temp)

#
for (i in 1:nrow(edge_fusion_temp))
{
  for (j in 1:223)
  {
    if(edge_fusion_temp[i,1]==TreeNumbers_disease[j,1])
    {
      edge_fusion_temp[i,5]=TreeNumbers_disease[j,3]
    }
  }
}
#
for (m in 1:nrow(edge_fusion_temp)) 
{
  for (n in 1:223)
  {
    if(edge_fusion_temp[m,2]==TreeNumbers_disease[n,1])
    {
      edge_fusion_temp[m,6]=TreeNumbers_disease[n,3]
    }
  }
}
rm(i,j,m,n)


### Hierarchical analysis
for (k in 1:nrow(edge_fusion_temp)) 
{
  if(substr(edge_fusion_temp[k,5],1,3)!=substr(edge_fusion_temp[k,6],1,3))
  {
    edge_fusion_temp[k,7]=1
  }else
  {
    if(substr(edge_fusion_temp[k,5],5,7)!=substr(edge_fusion_temp[k,6],5,7))
    {
      edge_fusion_temp[k,7]=2
    }else
    {
      if(substr(edge_fusion_temp[k,5],9,11)!=substr(edge_fusion_temp[k,6],9,11))
      {
        edge_fusion_temp[k,7]=3
      }else
      {
        if(substr(edge_fusion_temp[k,5],13,15)!=substr(edge_fusion_temp[k,6],13,15))
        {
          edge_fusion_temp[k,7]=4
        }else
        {
          if(substr(edge_fusion_temp[k,5],17,19)!=substr(edge_fusion_temp[k,6],17,19))
          {
            edge_fusion_temp[k,7]=5
          }else
          {
            edge_fusion_temp[k,7]=6
          }
        }
      }
    }
  }
}

# 
edge_fusion_temp = read.xlsx(file = "edge_fusion_temp.xlsx",sheetIndex = 1)
edge_fusion_temp = apply(edge_fusion_temp, 2, as.character)
cluster_group = as.matrix(read.table(file = "cluster_group.txt",header=F,sep="\t",row.names = NULL,check.names = F))

temp1 = matrix(data = NA, nrow = nrow(edge_fusion_temp), ncol = 3)
edge_fusion_temp = cbind(edge_fusion_temp, temp1)
rm(temp1)

#
for (i in 1:nrow(edge_fusion_temp))
{
  for (j in 1:223)
  {
    if(edge_fusion_temp[i,1]==cluster_group[j,1])
    {
      edge_fusion_temp[i,8]=cluster_group[j,2]
    }
  }
}
#
for (m in 1:nrow(edge_fusion_temp))
{
  for (n in 1:223)
  {
    if(edge_fusion_temp[m,2]==cluster_group[n,1])
    {
      edge_fusion_temp[m,9]=cluster_group[n,2]
    }
  }
}
rm(i,j,m,n)

edge_fusion_incluster = edge_fusion_temp[which(edge_fusion_temp[,8]==edge_fusion_temp[,9]),] 
edge_fusion_outcluster = edge_fusion_temp[which(edge_fusion_temp[,8]!=edge_fusion_temp[,9]),] 
edge_fusion_incluster[,c(8,9)] = as.numeric(edge_fusion_incluster[,c(8,9)])
edge_fusion_outcluster[,c(8,9)] = as.numeric(edge_fusion_outcluster[,c(8,9)])
# statistics
edge_fusion_incluster_order = edge_fusion_incluster[order(edge_fusion_incluster[,8]),]
edge_fusion_incluster_order1 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==1),]
edge_fusion_incluster_order2 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==2),]
edge_fusion_incluster_order3 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==3),]
edge_fusion_incluster_order4 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==4),]
edge_fusion_incluster_order5 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==5),]
edge_fusion_incluster_order6 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==6),]
edge_fusion_incluster_order7 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==7),]
edge_fusion_incluster_order8 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==8),]
edge_fusion_incluster_order9 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==9),]
edge_fusion_incluster_order10 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==10),]
edge_fusion_incluster_order11 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==11),]
edge_fusion_incluster_order12 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==12),]
edge_fusion_incluster_order13 = as.matrix(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==13),])
edge_fusion_incluster_order13 = t(edge_fusion_incluster_order13)
edge_fusion_incluster_order14 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==14),]
edge_fusion_incluster_order15 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==15),]
edge_fusion_incluster_order16 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==16),]
edge_fusion_incluster_order17 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==17),]
edge_fusion_incluster_order18 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==18),]
edge_fusion_incluster_order19 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==19),]
edge_fusion_incluster_order20 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==20),]
edge_fusion_incluster_order21 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==21),]
edge_fusion_incluster_order22 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==22),]
edge_fusion_incluster_order23 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==23),]
edge_fusion_incluster_order24 = edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==24),]

# 
edge_fusion_incluster_order[,c(8,9)] = as.numeric(edge_fusion_incluster_order[,c(8,9)])
cluster_hierarchy = matrix(data = NA, nrow = 24, ncol = 6)
rownames(cluster_hierarchy) = paste("cluster", 1:24, sep = "")
colnames(cluster_hierarchy) = paste("hierarchy", 1:6, sep = "")

for (i in 1:24)
{
  for (j in 1:6)
  {
    cluster_hierarchy[i,j]=length(which(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==i),7]=="1"))
    cluster_hierarchy[i,2]=length(which(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==i),7]=="2"))
    cluster_hierarchy[i,3]=length(which(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==i),7]=="3"))
    cluster_hierarchy[i,4]=length(which(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==i),7]=="4"))
    cluster_hierarchy[i,5]=length(which(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==i),7]=="5"))
    cluster_hierarchy[i,6]=length(which(edge_fusion_incluster_order[which(edge_fusion_incluster_order[,8]==i),7]=="6"))
  }
}
rm(i,j)
#
temp2 = matrix(data = NA, nrow = 24, ncol = 1)
cluster_hierarchy_percent = cbind(cluster_hierarchy, temp2)
rm(temp2)
cluster_hierarchy_percent[,7] = rowSums(cluster_hierarchy)
for (i in 1:24)
{
  for (j in 1:6)
  {
    cluster_hierarchy_percent[i,j]=round(cluster_hierarchy_percent[i,j]/cluster_hierarchy_percent[i,7]*100,2)
  }
}
cluster_hierarchy_percent = cluster_hierarchy_percent[,-7]
write.table(cluster_hierarchy_percent,file = "disease_hierarchy_percentage.txt",quote = F,sep = "\t",row.names = T,col.names = T)



###### Figure 5 ######
### 
cluster_hierarchy_percent = as.data.frame(read.table(file = "disease_hierarchy_percentage.txt",header = T, sep = "\t",row.names = 1,check.names = F,quote = "",stringsAsFactors = F))
library(pheatmap)
library(ggplot2)
library(reshape2)
cluster_hierarchy_percent = t(cluster_hierarchy_percent)
cluster_hierarchy_percent = cluster_hierarchy_percent*0.01
cluster_hierarchy_percent_temp = melt(cluster_hierarchy_percent)
colnames(cluster_hierarchy_percent_temp)=c("hierarchy","cluster","value")
# 
cluster_hierarchy_percent_temp1 = cluster_hierarchy_percent_temp[which(cluster_hierarchy_percent_temp[,3]>="0.0001"),]
# 
# mid = mean(cluster_hierarchy_percent_temp$value)
p = ggplot(data = cluster_hierarchy_percent_temp, aes(x = cluster, y = hierarchy, size = value)) + 
  xlab("Cluster") + ylab("Hierarchy") + theme_bw() +
  theme(panel.grid.major = element_blank(),panel.border = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) +
  geom_point(aes(color = value)) + scale_color_gradientn(colours = c('#6699CC','#FFFF99','#CC3333')) + scale_size(range = c(3,10)) +
  theme(aspect.ratio = 3/9)

p

# pheatmap(cluster_hierarchy_percent,cluster_cols=F,color = colorRampPalette(c("white", "firebrick3"))(50),cluster_rows=F,border_color="grey60",fontsize=10,show_rownames = T,show_colnames = T)

##### 
p = ggplot(data = cluster_hierarchy_percent_temp1, aes(x = cluster, y = hierarchy, size = value)) + 
  xlab("Cluster") + ylab("Hierarchy") + theme_bw() +
  theme(panel.grid.major = element_blank(),panel.border = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) +
  geom_point(aes(color = value)) + scale_color_gradientn(colours = c('#6699CC','#FFFF99','#CC3333')) + scale_size(range = c(3,10)) +
  theme(aspect.ratio = 3/9)

p








