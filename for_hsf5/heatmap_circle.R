library(tidyr)

library(readxl)
library(ComplexHeatmap)
mmc2_2_ <- read_excel("mmc2 (2).xlsx")
chosed_gene_set <- spermatid_gene_set
#############################
undif <- mmc2_2_[(mmc2_2_$gene)%in%chosed_gene_set$spermatid_development,]
undif  <- undif[,c(8,7,3)]
undif_wide<-spread(undif,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值undif_c1 <- undif[(undif$cluster)%in% "6",]
undif_wide <- as.data.frame(undif_wide)
rownames(undif_wide) <- undif_wide$gene
undif_wide <- undif_wide[,-1]
undif_mat<- as.matrix(undif_wide)
undif_mat_scale = t(scale(t(undif_mat)))
column_order = c("10","8","6","9","12","2","7","13","11","5","3","1","4","14")
svg(file="sperm_assembly.svg", width=6, height=30)
Heatmap(undif_mat_scale, 
        name = 'mat3', 
        column_order = column_order,
        cluster_rows = F, 
        cluster_columns = F)
dev.off()
#####
undif_mat<- as.matrix(undif_wide)
undif_mat_scale = t(scale(t(undif_mat)))
p1<-Heatmap(undif_mat_scale,column_names_side = "top",column_title = "sperm_assembly",name = "normalized_avg_logFC")

svg(file="sperm_assembly.svg", width=6, height=30)
p1 
dev.off()
svg(file="sperm_assembly.svg", width=10, height=30)
Heatmap(undif_wide, 
        name = 'mat3', 
        
        cluster_rows = T, 
        cluster_columns = T)

dev.off()

undif_mat<- as.matrix(rownames(undif_wide))
write.csv(undif_mat,file = "sp.csv")
