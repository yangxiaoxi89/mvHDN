load("disease_matrixs¡ªv3.rdata")
library(SNFtool)
library(clValid) 
library(pheatmap)
library(plotrix) 
library(xlsx)
neighbor_num=7#number of neighbors, usually (10~30),20 default.
alpha=0.5#hyperparameter, usually (0.3~0.8),0.5 default.
iteration_num=20#Number of Iterations, usually (10~20), 20 default.
rownames(GOsim_sci2015)=rownames(sim_natcom2014)
colnames(GOsim_sci2015)=colnames(sim_natcom2014)
rownames(dist_sci2015)=rownames(sim_natcom2014)
colnames(dist_sci2015)=colnames(sim_natcom2014)

#dist_sci2015:M-HDN,GOsim_sci2015:B-HDN,sim_natcom2014:S-HDN

####Compute the affinity matrix####
affini_list=NULL
affini_list=c(affini_list,list(affinityMatrix(dist_sci2015,neighbor_num, alpha))) 
affini_list=c(affini_list,list(affinityMatrix(1-GOsim_sci2015,neighbor_num, alpha)))
affini_list=c(affini_list,list(affinityMatrix(1-sim_natcom2014,neighbor_num, alpha)))

####SNF fusion####
fusion_result = SNF(affini_list, neighbor_num, iteration_num)#similarity matrix of data_result(after the SNF calculation) for heatmap
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

####the cluster number with spectral clustering###
connectivity_num=NULL
dunn_num=NULL
for (i in 2:40)
{
  cluster_group =spectralClustering(fusion_result,i)
  connectivity_num=c(connectivity_num,connectivity(distance = 1-fusion_result,clusters = cluster_group))
  dunn_num=c(dunn_num,dunn(distance = 1-fusion_result,clusters = cluster_group))
}
twoord.plot(lx=c(2:40),ly=connectivity_num,rx=c(2:40),ry=dunn_num
            ,xlim = c(1,40),lylim = c(40,220),rylim = c(0.74,0.94),lcol = "red",rcol = "blue"
            ,xlab = "Cluster Number",ylab = "Connectivity Index",rylab = "Dunn Index"
            ,main = "Cluster Validity",type = c("line","line"),lwd=2)
            #,do.first = 'plot_bg(col = \'gray\'); grid(col = \'white\', lty = 2)')

#the cluster number = 24£¬connectivity = 160.06270,dunn = 0.9232009
cluster_group =spectralClustering(fusion_result,24)
order_fusion_result=fusion_result[order(cluster_group),order(cluster_group)] 
temp=order_fusion_result-diag(diag(order_fusion_result)) 
temp=temp/max(temp)+diag(nrow(order_fusion_result)) 

####heatmap####
Mesh_disease=read.xlsx(file = "MeshDiseaseName_total223.xlsx",sheetIndex = 1)

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
  max_nam=names(cluster_sta)[which(cluster_sta==max(cluster_sta))] 
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
###annotation##
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

truelabel_colors=c("bisque","springgreen2","burlywood","cadetblue","chartreuse","chocolate","coral","cornsilk","firebrick","dimgrey","gold","maroon3","hotpink","khaki","lightblue","skyblue","green")
truelabel_mapping=cbind(truelabel_mapping,truelabel_colors)
cluster_group_colors=c("bisque","springgreen2","burlywood","cadetblue","chartreuse","chocolate","coral","cornsilk","firebrick","dimgrey","gold","maroon3","hotpink","khaki","lightblue","skyblue","green","pink","yellow","orange","red","papayawhip","violet","purple")
truelabel_mapping=cbind(truelabel_mapping,as.matrix(read.xlsx(file = "Mesh_ID_map.xlsx",sheetIndex = 1))[,2])
colnames(truelabel_mapping)[4]="Name"

#
M_label=cbind(cluster_group,truelabel) 
colnames(M_label)=c("spectralClustering","TrueLabel")
#M_label_colors=t(apply(M_label,1,getColorsForGroups))
#cl <- colors()
M_label_colors=cbind("spectralClustering"=getColorsForGroups(M_label[,"spectralClustering"],colors=cluster_group_colors),
                     "TrueLabel"=getColorsForGroups(M_label[,"TrueLabel"],colors=truelabel_colors))

displayClustersWithHeatmap(fusion_result, cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.3,cexCol = 0.3)

displayClustersWithHeatmap(affini_list[[1]], cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.3,cexCol = 0.3)
displayClustersWithHeatmap(affini_list[[2]], cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.3,cexCol = 0.3)
displayClustersWithHeatmap(affini_list[[3]], cluster_group, M_label_colors,col=colorRampPalette(c("aliceblue","navy"))(1000),margins = c(10,10),cexRow=0.3,cexCol = 0.3)


####Contribution of data####
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

#
contribution_matrix_group_list=NULL
contribution_sta_group_list=NULL
for (i in 1:24)
{
  contribution_matrix_group_temp=contribution_matrix[which(cluster_group==i),which(cluster_group==i)]
  contribution_sta_group_list=c(contribution_sta_group_list,list(table(contribution_matrix_group_temp[which(contribution_matrix_group_temp!=0)])))
  contribution_matrix_group_list=c(contribution_matrix_group_list,list(contribution_matrix_group_temp))

  rm(contribution_matrix_group_temp)
}


########
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

######Statistic the relationship between disease treenumbers and clusters######
save(edge_fusion,edge_nc2014,edge_sci2015,edge_sci2015GO,file = "edge.RData")
############
write.table(node,file = "node.txt",quote = F,sep = "\t",row.names = F)
write.table(edge_fusion[which(x>=quantile(x,probs = 0.8)),],file = "edge_fusion.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(edge_sci2015[which(x1>=quantile(x1,probs = 0.8)),],file = "edge_sci2015.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(edge_sci2015GO[which(x2>=quantile(x2,probs = 0.8)),],file = "edge_sci2015GO.txt",quote = F,sep = "\t",row.names = F,col.names = F)
write.table(edge_nc2014[which(x3>=quantile(x3,probs = 0.8)),],file = "edge_nc2014.txt",quote = F,sep = "\t",row.names = F,col.names = F)
rm(x)
rm(x1)
rm(x2)
rm(x3)
save(contribution_matrix,cor_matrix,dist_sci2015,fusion_result,GOsim_sci2015,Mesh_disease,sim_natcom2014,affini_list,cluster_group,contribution_matrix_group_list,contribution_sta_group_list,truelabel,truelabel_mapping,file = "fusion_v3.0.RData")


####heatmap of disease_mesh and cluster####
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
pheatmap(mesh_cluster_matrix,cluster_cols=F,color = colorRampPalette(c("white", "firebrick3"))(50),cluster_rows=F,border_color="grey60",fontsize=10,show_rownames = T,show_colnames = T)
write.table(mesh_cluster_matrix,file = "mesh_cluster_matrix.txt",quote = F,sep = "\t")


########
write.table(affini_list[[1]],file = "affinity_matrix_sci2015.txt",quote = F,sep="\t")
write.table(affini_list[[2]],file = "affinity_matrix_sci2015GO.txt",quote = F,sep="\t")
write.table(affini_list[[3]],file = "affinity_matrix_nc2014.txt",quote = F,sep="\t")
write.table(cbind(rownames(fusion_result),cluster_group),file = "cluster_group.txt",quote = F,sep = "\t",row.names = F,col.names = F)


