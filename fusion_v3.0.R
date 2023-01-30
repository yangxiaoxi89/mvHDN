load("disease_matrixs_to_Hesong—v3.rdata")
library(SNFtool)
library(clValid)
library(pheatmap)
library(plotrix)
library(xlsx)
library(RColorBrewer)
library(reshape2)
library(ggplot2)

neighbor_num=7  # number of neighbors, usually (10~30),20 default.
alpha=0.5  # hyperparameter, usually (0.3~0.8),0.5 default.
iteration_num=20  # Number of Iterations, usually (10~20), 20 default.
rownames(GOsim_sci2015)=rownames(sim_natcom2014)
colnames(GOsim_sci2015)=colnames(sim_natcom2014)
rownames(dist_sci2015)=rownames(sim_natcom2014)
colnames(dist_sci2015)=colnames(sim_natcom2014)

#  dist_sci2015: M-HDN, GOsim_sci2015: B-HDN, sim_natcom2014: S-HDN

#### Compute the affinity matrix ####
affini_list=NULL
affini_list=c(affini_list,list(affinityMatrix(dist_sci2015,neighbor_num, alpha)))
affini_list=c(affini_list,list(affinityMatrix(1-GOsim_sci2015,neighbor_num, alpha)))
affini_list=c(affini_list,list(affinityMatrix(1-sim_natcom2014,neighbor_num, alpha)))

#### SNF fusion ####
fusion_result = SNF(affini_list, neighbor_num, iteration_num)  # similarity matrix of data_result
rownames(fusion_result)=rownames(dist_sci2015)
colnames(fusion_result)=rownames(dist_sci2015)
write.table(fusion_result,file = "fusion_result.txt",quote = F,sep="\t")

cor_matrix=matrix(data = 1,nrow = 4,ncol = 4)
rownames(cor_matrix)=c("dist_sci2015","GOsim_sci2015","sim_natcom2014","fusion")
colnames(cor_matrix)=rownames(cor_matrix)

# the correlation between the matrices
for (i in 1:3) 
{
  for (j in 1:3)
  {
    cor_matrix[i,j]=cor(as.vector(affini_list[[i]]),as.vector(affini_list[[j]]))
    
  }
}
cor_matrix[1,4]=cor(as.vector(affini_list[[1]]),as.vector(fusion_result))
cor_matrix[2,4]=cor(as.vector(affini_list[[2]]),as.vector(fusion_result))
cor_matrix[3,4]=cor(as.vector(affini_list[[3]]),as.vector(fusion_result))
cor_matrix[4,1]=cor_matrix[1,4]
cor_matrix[4,2]=cor_matrix[2,4]
cor_matrix[4,3]=cor_matrix[3,4]

write.table(cor_matrix,file = "cor_matrix.txt",quote = F,sep = "\t")


#### the cluster number with spectral clustering ####
connectivity_num=NULL
dunn_num=NULL
for (i in 4:34)
{
  cluster_group =spectralClustering(fusion_result,i)
  connectivity_num=c(connectivity_num,connectivity(distance = 1-fusion_result,clusters = cluster_group))
  dunn_num=c(dunn_num,dunn(distance = 1-fusion_result,clusters = cluster_group))
}
twoord.plot(lx=c(3:33),ly=connectivity_num,rx=c(3:33),ry=dunn_num
            ,xlim = c(3,33),lylim = c(40,220),rylim = c(0.74,0.94),lcol = "red",rcol = "blue"
            ,xlab = "Cluster Number",ylab = "Connectivity Index",rylab = "Dunn Index"
            ,main = "Cluster Validity",type = c("b","b"),lwd=2)

######################
cluster_group =spectralClustering(fusion_result,24)
order_fusion_result=fusion_result[order(cluster_group),order(cluster_group)] 
temp=order_fusion_result-diag(diag(order_fusion_result)) 
temp=temp/max(temp)+diag(nrow(order_fusion_result)) 

#### heatmap ####
Mesh_disease=read.xlsx(file = "MeshDiseaseName_total223.xlsx",sheetIndex = 1)

###2017年6月2日加，为了将热图上的Meshcode颜色条调整得更好看###
disease_slim_temp=as.matrix(read.xlsx(file = "Mesh_disease_slim.xlsx",sheetIndex = 1))
disease_slim=NULL
for (i in 1:nrow(disease_slim_temp)) 
{
  disease_slim=c(disease_slim,list(unique(disease_slim_temp[i,which(!is.na(disease_slim_temp[i,]))])))
}
rm(disease_slim_temp)

for (i in 1:24) 
{
  num=which(cluster_group==i)
  cluster_disease_temp=NULL
  for (j in 1:length(num))
  {
    cluster_disease_temp=c(cluster_disease_temp,disease_slim[[num[j]]])
  }
  cluster_sta=table(cluster_disease_temp)
  max_nam=names(cluster_sta)[which(cluster_sta==max(cluster_sta))] #最多的那个mesh ID作为max_nam
  for (j in 1:length(num))
  {
    if(length(disease_slim[[num[j]]])>1)
    {
      num1=which(disease_slim[[num[j]]]==max_nam)
      if (length(num1)==1)
      {
        disease_slim[[num[j]]]=max_nam
      }else
      {
        disease_slim[[num[j]]]=disease_slim[[num[j]]][1]
      }
      rm(num1)
    }
  }
  rm(num)
  rm(max_nam)
  rm(cluster_disease_temp)
  rm(cluster_sta)
}

disease_slim[[23]]="C17"
disease_slim[[32]]="C10"
disease_slim[[93]]="C13"
disease_slim[[179]]="C13"
disease_slim[[205]]="C10"

truelabel=unlist(disease_slim)
rm(disease_slim)
### annotation ###
#truelabel=as.vector(Mesh_disease[,8])
#for (i in 1:length(truelabel))
#{
#  truelabel[i]=substr(truelabel[i],start = 1,stop = 3)
#}

truelabel_mapping=matrix(data = NA,nrow = length(unique(truelabel)),ncol = 2)
truelabel_mapping[,1]=unique(sort(truelabel))
truelabel_mapping[,2]=1:length(unique(truelabel))
truelabel_temp=matrix(data = 0,nrow = length(truelabel),ncol = 1)

for (i in 1:length(truelabel)) 
{
  truelabel_temp[i,1]=as.numeric(truelabel_mapping[which(truelabel_mapping[,1]==truelabel[i]),2])
}
truelabel=truelabel_temp 
rm(truelabel_temp)

##
###### Figure 2 ######
truelabel_colors=c("#6C6B74","#80ADD7","#C05640","#1D65A6","#A882C1","#0ABDA0","#D6618F","#F3D4A0","#D4DCA9","#5C868D","#9199BE","#D9AC2A","#F46A4E","#BF9D7A","#00743F","#192E5B","#74593D")
truelabel_mapping=cbind(truelabel_mapping,truelabel_colors)
cluster_group_colors=c("bisque","steelblue4","burlywood","cadetblue","palevioletred4","chocolate","coral","cornsilk","firebrick","dimgrey","gold","maroon3","hotpink","khaki","lightblue","skyblue","seagreen","pink","slateblue","orange","red","papayawhip","violet","mediumaquamarine")
truelabel_mapping=cbind(truelabel_mapping,as.matrix(read.xlsx(file = "Mesh_ID_map.xlsx",sheetIndex = 1))[,2])
colnames(truelabel_mapping)[4]="Name"

M_label=cbind(cluster_group,truelabel)
colnames(M_label)=c("spectralClustering","TrueLabel")
#M_label_colors=t(apply(M_label,1,getColorsForGroups))
#cl <- colors()
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],colors=cluster_group_colors),
                     "TrueLabel"=getColorsForGroups(M_label[,"TrueLabel"],colors=truelabel_colors))

displayClustersWithHeatmap(fusion_result, cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.5,cexCol = 0.5)

displayClustersWithHeatmap(affini_list[[1]], cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.5,cexCol = 0.5)
displayClustersWithHeatmap(affini_list[[2]], cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.5,cexCol = 0.5)
displayClustersWithHeatmap(affini_list[[3]], cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.5,cexCol = 0.5)


#### Contribution of data types ####
contribution_matrix=fusion_result
contribution_matrix[which(contribution_matrix!=0)]=0
for (i in 1:nrow(contribution_matrix))
{
  for (j in i:nrow(contribution_matrix))
  {
    if (j!=i)
    {
      temp=matrix(data = 0,nrow = 3,ncol = 1)
      temp[1]=affini_list[[1]][i,j]
      temp[2]=affini_list[[2]][i,j]
      temp[3]=affini_list[[3]][i,j]
      rank_temp=order(-temp)
      temp=temp[rank_temp]
      chazhi=(temp[1]-temp[2])/temp[2]
      chazhi=c(chazhi,(temp[2]-temp[3])/temp[3])
      num=length(which(chazhi<=0.1))
      if (chazhi[1]>0.1)
      {
        contribution_matrix[i,j]=rank_temp[1]
      }else
      {
        if (chazhi[2]>0.1)
        {
          a=c(rank_temp[1],rank_temp[2])
          if (length(intersect(a,c(1,2)))==2)
          {
            contribution_matrix[i,j]=4
          }
          if (length(intersect(a,c(1,3)))==2)
          {
            contribution_matrix[i,j]=5
          }
          if (length(intersect(a,c(2,3)))==2)
          {
            contribution_matrix[i,j]=6
          }
        }else
        {
          contribution_matrix[i,j]=7
        }
      }
      #rm(chazhi)
      #rm(temp)
      #rm(rank_temp)
    }
  }
}
contribution_sta=table(contribution_matrix[which(contribution_matrix!=0)])
##
contribution_matrix_group_list=NULL
contribution_sta_group_list=NULL
for (i in 1:24)
{
  contribution_matrix_group_temp=contribution_matrix[which(cluster_group==i),which(cluster_group==i)]
  contribution_sta_group_list=c(contribution_sta_group_list,list(table(contribution_matrix_group_temp[which(contribution_matrix_group_temp!=0)])))
  contribution_matrix_group_list=c(contribution_matrix_group_list,list(contribution_matrix_group_temp))
  
  rm(contribution_matrix_group_temp)
}


##########
node=matrix(data = NA,nrow = nrow(fusion_result),ncol = 5)
node[,1:3]=as.matrix(Mesh_disease[,c(1,3,10)])
node[,4]=truelabel
node[,5]=cluster_group
colnames(node)=c("Name","DiseaseID","SlimMappings","Truelabel","Cluster_group")

edge_fusion=matrix(data = NA,nrow = nrow(fusion_result)*(nrow(fusion_result)-1)/2,ncol = 4) #223*222/2
edge_sci2015=matrix(data = NA,nrow = nrow(fusion_result)*(nrow(fusion_result)-1)/2,ncol = 3)
edge_sci2015GO=edge_sci2015
edge_nc2014=edge_sci2015
k=1
for (i in 1:nrow(fusion_result))
{
  for (j in i:nrow(fusion_result))
  {
    if (j!=i)
    {
    edge_fusion[k,1]=rownames(fusion_result)[i]
    edge_fusion[k,2]=rownames(fusion_result)[j]
    edge_fusion[k,3]=fusion_result[i,j]
    edge_fusion[k,4]=contribution_matrix[i,j]
    edge_sci2015[k,1]=rownames(fusion_result)[i]
    edge_sci2015[k,2]=rownames(fusion_result)[j]
    edge_sci2015[k,3]=affini_list[[1]][i,j]
    edge_sci2015GO[k,1]=rownames(fusion_result)[i]
    edge_sci2015GO[k,2]=rownames(fusion_result)[j]
    edge_sci2015GO[k,3]=affini_list[[2]][i,j]
    edge_nc2014[k,1]=rownames(fusion_result)[i]
    edge_nc2014[k,2]=rownames(fusion_result)[j]
    edge_nc2014[k,3]=affini_list[[3]][i,j]
    k=k+1
    }
  }
}
x=as.numeric(edge_fusion[,3])
x1=as.numeric(edge_sci2015[,3])
x2=as.numeric(edge_sci2015GO[,3])
x3=as.numeric(edge_nc2014[,3])

###### Statistic the relationship between disease tree numbers and clusters ######
save(edge_fusion,edge_nc2014,edge_sci2015,edge_sci2015GO,file = "edge.RData")
######以下为原始代码######
write.table(node,file = "node.txt",quote = F,sep = "\t",row.names = F)
write.table(edge_fusion[which(x>=quantile(x,probs = 0.8)),],file = "edge_fusion.txt",quote = F,sep = "\t",row.names = F,col.names = F) #quantile百分位数，20%
write.table(edge_sci2015[which(x1>=quantile(x1,probs = 0.8)),],file = "edge_sci2015.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(edge_sci2015GO[which(x2>=quantile(x2,probs = 0.8)),],file = "edge_sci2015GO.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(edge_nc2014[which(x3>=quantile(x3,probs = 0.8)),],file = "edge_nc2014.txt",quote = F,sep = "\t",row.names = F,col.names = F)
rm(x)
rm(x1)
rm(x2)
rm(x3)
save(contribution_matrix,cor_matrix,dist_sci2015,fusion_result,GOsim_sci2015,Mesh_disease,sim_natcom2014,affini_list,cluster_group,contribution_matrix_group_list,contribution_sta_group_list,truelabel,truelabel_mapping,file = "fusion_v3.0.RData")


#### heat of disease mesh and cluster ####
mesh_cluster_matrix=matrix(data = 0,nrow = nrow(truelabel_mapping),ncol = 24)
for (i in 1:nrow(mesh_cluster_matrix))
{
  num1=which(truelabel==truelabel_mapping[i,2])
  for (j in 1:ncol(mesh_cluster_matrix))
  {
    num2=which(cluster_group==j)
    mesh_cluster_matrix[i,j]=length(intersect(num1,num2))
  }
}
rownames(mesh_cluster_matrix)=truelabel_mapping[,4]
colnames(mesh_cluster_matrix)=paste("Cluster",1:24,sep = "")
# pheatmap(mesh_cluster_matrix,cluster_cols=F,color = colorRampPalette(c("white", "firebrick3"))(50),cluster_rows=F,border_color="grey60",fontsize=10,show_rownames = T,show_colnames = T)
write.table(mesh_cluster_matrix,file = "mesh_cluster_matrix.txt",quote = F,sep = "\t")



##### Figure 4 #####
node_temp = as.data.frame(read.table(file = "node.txt",header = T, sep = "\t",row.names = NULL,check.names = F,quote = "",stringsAsFactors = F))
label_mapping = as.matrix(read.csv(file = "Mesh_ID_map.csv",header=T,sep=",",row.names = NULL,check.names = F))
#
label_mapping = cbind(label_mapping,"NA")
label_mapping[,3] = c(1:17)
colnames(label_mapping)[3] = "TrueNumber" 
#
node_temp = cbind(node_temp,"NA")
colnames(node_temp)[6] = "DescripName"
node_temp[,4] = as.numeric(node_temp[,4])
node_temp[,5] = as.numeric(node_temp[,5])
for (i in 1:nrow(node_temp)) 
{
  node_temp[i,6] = label_mapping[which(label_mapping[,3]==node_temp[i,4]),2]
  print(i)
}
node_temp = node_temp[,c(1,5,6)]

##  Function
Ratio = function(Data)
{
  #
  Data[,2]=as.numeric(Data[,2]) 
  ratio_matrix=matrix(0, nrow = length(unique(Data[,2])), ncol = length(unique(Data[,3]))) 
  rownames(ratio_matrix)=paste("cluster",c(1:length(unique(Data[,2]))),sep = "")
  colnames(ratio_matrix)=sort(unique(Data[,3]))
  
  for (i in 1:length(unique(Data[,2])))
  {
    for (j in 1:length(unique(Data[,3]))) 
    {
      ratio_matrix[i,j]=length(which(Data[which(Data[,2]==i),3]==colnames(ratio_matrix)[j]))/length(which(Data[,2]==i))
    }
  }
  rm(i,j)
  return(ratio_matrix)
}
# 
ratio_node = Ratio(node_temp)

node_temp = ratio_node
node_temp = t(node_temp)
node_temp = melt(node_temp)
colnames(node_temp)=c("true_label","cluster","value")

## 
color = c("#6C6B74","#80ADD7","#C05640","#1D65A6","#A882C1","#0ABDA0","#D6618F","#F3D4A0","#D4DCA9","#5C868D","#9199BE","#D9AC2A","#F46A4E","#BF9D7A","#00743F","#192E5B","#74593D")
color_label = c("Bacterial infection or mycosis","Neoplasms","Musculoskeletal disease",
                "Digestive system disease","Respiratory tract disease","Nervous system disease",
                "Eye diseases","Urologic and male genital diseases","Female genital diseases and pregnancy complications",
                "Cardiovascular disease","Hemic and lymphatic diseases","Congenital abnormality",
                "Skin and connective tissue diseases","Nutritional and metabolic diseases","Endocrine system disease",
                "Immune system diseases","Pathological Conditions, signs and symptoms")
names(color)=color_label
#
p = ggplot(data = node_temp, aes(cluster, 100*value, fill = true_label)) +
  geom_col(position = "stack", width = 0.8,size = 0.5) +
  scale_fill_manual(values = color) + 
  labs(x = "", y = "Percentage (%)") + scale_y_continuous(expand = c(0,0)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'transparent', fill = 'transparent'), strip.text = element_text(size = 15)) +
  theme(axis.text = element_text(angle = 45, vjust = 1, hjust = 1,size = 20), axis.title = element_text(size = 15), legend.title = element_blank(), legend.text = element_text(size = 15))
p

write.table(mesh_cluster_matrix,file = "mesh_cluster_matrix.txt",quote = F,sep = "\t")



####
write.table(affini_list[[1]],file = "affinity_matrix_sci2015.txt",quote = F,sep="\t")
write.table(affini_list[[2]],file = "affinity_matrix_sci2015GO.txt",quote = F,sep="\t")
write.table(affini_list[[3]],file = "affinity_matrix_nc2014.txt",quote = F,sep="\t")
write.table(cbind(rownames(fusion_result),cluster_group),file = "cluster_group.txt",quote = F,sep = "\t",row.names = F,col.names = F)




#######################################################################################################
### Performance Comparison (ANF、COCA、SNF-CC) ####

###### SNF ######
load("disease_matrixs_to_Hesong—v3.rdata")
library(SNFtool)
library(clValid) 
library(cluster)

rownames(GOsim_sci2015)=rownames(sim_natcom2014)
colnames(GOsim_sci2015)=colnames(sim_natcom2014)
rownames(dist_sci2015)=rownames(sim_natcom2014)
colnames(dist_sci2015)=colnames(sim_natcom2014)

##
dist_list = NULL
dist_list = c(dist_list,list(dist_sci2015))
dist_list = c(dist_list,list(1-GOsim_sci2015))
dist_list = c(dist_list,list(1-sim_natcom2014))
affini_list = lapply(dist_list, function(x) affinityMatrix(x,20,0.5))

##
fusion_result = SNF(affini_list, 10, 20)
cluster_group = spectralClustering(fusion_result,24)

# connectivity
connectivity(distance = 1-fusion_result,clusters = cluster_group)    
# dunn
dunn(distance = 1-fusion_result,clusters = cluster_group)   



###### ANF ######
load("disease_matrixs_to_Hesong—v3.rdata")
library(ANF)

rownames(GOsim_sci2015)=rownames(sim_natcom2014)
colnames(GOsim_sci2015)=colnames(sim_natcom2014)
rownames(dist_sci2015)=rownames(sim_natcom2014)
colnames(dist_sci2015)=colnames(sim_natcom2014)

##
dist_list = NULL
dist_list = c(dist_list,list(dist_sci2015))
dist_list = c(dist_list,list(1-GOsim_sci2015))
dist_list = c(dist_list,list(1-sim_natcom2014))
affini_list = lapply(dist_list, function(x) affinity_matrix(x,20))

##
fusion_result = ANF(affini_list, K=10)
cluster_group = spectral_clustering(fusion_result,k=24)

# connectivity
connectivity(distance = 1-fusion_result,clusters = cluster_group)   
# dunn
dunn(distance = 1-fusion_result,clusters = cluster_group)   



###### COCA ######
load("disease_matrixs_to_Hesong—v3.rdata")
library(ConsensusClusterPlus)
library(kernlab)

rownames(GOsim_sci2015)=rownames(sim_natcom2014)
colnames(GOsim_sci2015)=colnames(sim_natcom2014)
rownames(dist_sci2015)=rownames(sim_natcom2014)
colnames(dist_sci2015)=colnames(sim_natcom2014)


set.seed(798)
#
temp = specc(1-dist_sci2015,centers=24)
clusternum = 24
inputdataconsen1 = matrix(data = NA, nrow = clusternum, ncol = nrow(dist_sci2015))
for (i in 1:ncol(inputdataconsen1)) 
{
  if(temp@.Data[i] == 1)
  {
    inputdataconsen1[1,i] = 1
  }else if(temp@.Data[i] == 2)
  {
    inputdataconsen1[2,i] = 1
  }else if(temp@.Data[i] == 3)
  {
    inputdataconsen1[3,i] = 1
  }else if(temp@.Data[i] == 4)
  {
    inputdataconsen1[4,i] = 1
  }else if(temp@.Data[i] == 5)
  {
    inputdataconsen1[5,i] = 1
  }else if(temp@.Data[i] == 6)
  {
    inputdataconsen1[6,i] = 1
  }else if(temp@.Data[i] == 7)
  {
    inputdataconsen1[7,i] = 1
  }else if(temp@.Data[i] == 8)
  {
    inputdataconsen1[8,i] = 1
  }else if(temp@.Data[i] == 9)
  {
    inputdataconsen1[9,i] = 1
  }else if(temp@.Data[i] == 10)
  {
    inputdataconsen1[10,i] = 1
  }else if(temp@.Data[i] == 11)
  {
    inputdataconsen1[11,i] = 1
  }else if(temp@.Data[i] == 12)
  {
    inputdataconsen1[12,i] = 1
  }else if(temp@.Data[i] == 13)
  {
    inputdataconsen1[13,i] = 1
  }else if(temp@.Data[i] == 14)
  {
    inputdataconsen1[14,i] = 1
  }else if(temp@.Data[i] == 15)
  {
    inputdataconsen1[15,i] = 1
  }else if(temp@.Data[i] == 16)
  {
    inputdataconsen1[16,i] = 1
  }else if(temp@.Data[i] == 17)
  {
    inputdataconsen1[17,i] = 1
  }else if(temp@.Data[i] == 18)
  {
    inputdataconsen1[18,i] = 1
  }else if(temp@.Data[i] == 19)
  {
    inputdataconsen1[19,i] = 1
  }else if(temp@.Data[i] == 20)
  {
    inputdataconsen1[20,i] = 1
  }else if(temp@.Data[i] == 21)
  {
    inputdataconsen1[21,i] = 1
  }else if(temp@.Data[i] == 22)
  {
    inputdataconsen1[22,i] = 1
  }else if(temp@.Data[i] == 23)
  {
    inputdataconsen1[23,i] = 1
  }else
  {
    inputdataconsen1[24,i] = 1
  }
  print(i)
}
inputdataconsen1[is.na(inputdataconsen1)] = 0
#
temp = specc(GOsim_sci2015,centers=24)
clusternum = 24
inputdataconsen2 = matrix(data = NA, nrow = clusternum, ncol = nrow(GOsim_sci2015))
for (i in 1:ncol(inputdataconsen2)) 
{
  if(temp@.Data[i] == 1)
  {
    inputdataconsen2[1,i] = 1
  }else if(temp@.Data[i] == 2)
  {
    inputdataconsen2[2,i] = 1
  }else if(temp@.Data[i] == 3)
  {
    inputdataconsen1[3,i] = 1
  }else if(temp@.Data[i] == 4)
  {
    inputdataconsen2[4,i] = 1
  }else if(temp@.Data[i] == 5)
  {
    inputdataconsen2[5,i] = 1
  }else if(temp@.Data[i] == 6)
  {
    inputdataconsen2[6,i] = 1
  }else if(temp@.Data[i] == 7)
  {
    inputdataconsen2[7,i] = 1
  }else if(temp@.Data[i] == 8)
  {
    inputdataconsen2[8,i] = 1
  }else if(temp@.Data[i] == 9)
  {
    inputdataconsen2[9,i] = 1
  }else if(temp@.Data[i] == 10)
  {
    inputdataconsen2[10,i] = 1
  }else if(temp@.Data[i] == 11)
  {
    inputdataconsen2[11,i] = 1
  }else if(temp@.Data[i] == 12)
  {
    inputdataconsen2[12,i] = 1
  }else if(temp@.Data[i] == 13)
  {
    inputdataconsen2[13,i] = 1
  }else if(temp@.Data[i] == 14)
  {
    inputdataconsen2[14,i] = 1
  }else if(temp@.Data[i] == 15)
  {
    inputdataconsen2[15,i] = 1
  }else if(temp@.Data[i] == 16)
  {
    inputdataconsen2[16,i] = 1
  }else if(temp@.Data[i] == 17)
  {
    inputdataconsen2[17,i] = 1
  }else if(temp@.Data[i] == 18)
  {
    inputdataconsen2[18,i] = 1
  }else if(temp@.Data[i] == 19)
  {
    inputdataconsen2[19,i] = 1
  }else if(temp@.Data[i] == 20)
  {
    inputdataconsen2[20,i] = 1
  }else if(temp@.Data[i] == 21)
  {
    inputdataconsen2[21,i] = 1
  }else if(temp@.Data[i] == 22)
  {
    inputdataconsen2[22,i] = 1
  }else if(temp@.Data[i] == 23)
  {
    inputdataconsen2[23,i] = 1
  }else
  {
    inputdataconsen2[24,i] = 1
  }
  print(i)
}
inputdataconsen2[is.na(inputdataconsen2)] = 0
#
temp = specc(sim_natcom2014,centers=24)
clusternum = 24
inputdataconsen3 = matrix(data = NA, nrow = clusternum, ncol = nrow(sim_natcom2014))
for (i in 1:ncol(inputdataconsen3)) 
{
  if(temp@.Data[i] == 1)
  {
    inputdataconsen3[1,i] = 1
  }else if(temp@.Data[i] == 2)
  {
    inputdataconsen3[2,i] = 1
  }else if(temp@.Data[i] == 3)
  {
    inputdataconsen3[3,i] = 1
  }else if(temp@.Data[i] == 4)
  {
    inputdataconsen3[4,i] = 1
  }else if(temp@.Data[i] == 5)
  {
    inputdataconsen3[5,i] = 1
  }else if(temp@.Data[i] == 6)
  {
    inputdataconsen3[6,i] = 1
  }else if(temp@.Data[i] == 7)
  {
    inputdataconsen3[7,i] = 1
  }else if(temp@.Data[i] == 8)
  {
    inputdataconsen3[8,i] = 1
  }else if(temp@.Data[i] == 9)
  {
    inputdataconsen3[9,i] = 1
  }else if(temp@.Data[i] == 10)
  {
    inputdataconsen3[10,i] = 1
  }else if(temp@.Data[i] == 11)
  {
    inputdataconsen3[11,i] = 1
  }else if(temp@.Data[i] == 12)
  {
    inputdataconsen3[12,i] = 1
  }else if(temp@.Data[i] == 13)
  {
    inputdataconsen3[13,i] = 1
  }else if(temp@.Data[i] == 14)
  {
    inputdataconsen3[14,i] = 1
  }else if(temp@.Data[i] == 15)
  {
    inputdataconsen3[15,i] = 1
  }else if(temp@.Data[i] == 16)
  {
    inputdataconsen3[16,i] = 1
  }else if(temp@.Data[i] == 17)
  {
    inputdataconsen3[17,i] = 1
  }else if(temp@.Data[i] == 18)
  {
    inputdataconsen3[18,i] = 1
  }else if(temp@.Data[i] == 19)
  {
    inputdataconsen3[19,i] = 1
  }else if(temp@.Data[i] == 20)
  {
    inputdataconsen3[20,i] = 1
  }else if(temp@.Data[i] == 21)
  {
    inputdataconsen3[21,i] = 1
  }else if(temp@.Data[i] == 22)
  {
    inputdataconsen3[22,i] = 1
  }else if(temp@.Data[i] == 23)
  {
    inputdataconsen3[23,i] = 1
  }else
  {
    inputdataconsen3[24,i] = 1
  }
  print(i)
}
inputdataconsen3[is.na(inputdataconsen3)] = 0

inputdataconsen = rbind(inputdataconsen1, inputdataconsen2, inputdataconsen3)

## clustering
rcc = ConsensusClusterPlus(inputdataconsen,maxK=24,reps=10,pItem=0.8,pFeature=1,title="example",distance="pearson",clusterAlg="hc",seed = 798)
a = rcc[[24]][["consensusClass"]] 
b = rcc[[24]][["consensusMatrix"]]

# connectivity
connectivity(-log(rcc[[24]][["consensusMatrix"]]),rcc[[24]][["consensusClass"]])    
# dunn
dunn(1-rcc[[24]][["consensusMatrix"]],rcc[[24]][["consensusClass"]])  



###### SNF-CC ######
load("disease_matrixs_to_Hesong—v3.rdata")
library(CancerSubtypes)


rownames(GOsim_sci2015)=rownames(sim_natcom2014)
colnames(GOsim_sci2015)=colnames(sim_natcom2014)
rownames(dist_sci2015)=rownames(sim_natcom2014)
colnames(dist_sci2015)=colnames(sim_natcom2014)
dist_sci2015 = 1-dist_sci2015
dist_sci2015 = dist_sci2015-diag(diag(dist_sci2015)) 
GOsim_sci2015 = GOsim_sci2015-diag(diag(GOsim_sci2015)) 
sim_natcom2014 = sim_natcom2014-diag(diag(sim_natcom2014)) 

list_temp = list(dist_sci2015,GOsim_sci2015,sim_natcom2014)
result = ExecuteSNF.CC(list_temp, clusterNum=24, K=10, alpha=0.5, t=20,
                       maxK = 24, pItem = 0.8,reps=10, 
                       finalLinkage ="average")


# connectivity
connectivity(1-result$distanceMatrix,result$group)  
# dunn
dunn(1-result$distanceMatrix,result$group)    




######### Figure 6 #############
#
compar_result = read.xlsx(file = "Compare.xlsx",sheetIndex = 1)
compar_result$Model = factor(compar_result$Model,levels = c("SNF","ANF","COCA","SNFCC"))


# 
ggplot(compar_result,aes(Performance, Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.6) +
  scale_fill_manual(values = c("indianred3","chocolate1","goldenrod1","mediumturquoise")) +
 # theme(panel.grid = element_blank(), panel.background = element_rect(color = 'transparent', fill = 'transparent')) +
  scale_y_continuous(name = "Conn",
        sec.axis = sec_axis(~./100,name = "Denn"))


   










