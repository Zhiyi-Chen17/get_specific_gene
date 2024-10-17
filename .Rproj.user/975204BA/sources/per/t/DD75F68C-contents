

library(tidyr)

library(readxl)
mmc2_2_ <- read_excel("mmc2 (2).xlsx")
nitaharaDatasetS2 <- read_excel("nitaharaDatasetS2.xlsx", 
                                col_names = FALSE)
nitaharaDatasetS2 <- nitaharaDatasetS2[-1,]
colnames(nitaharaDatasetS2) <- nitaharaDatasetS2[1,]
nitaharaDatasetS2 <- nitaharaDatasetS2[-1,]
#############################
undif <- mmc2_2_[(mmc2_2_$gene)%in%nitaharaDatasetS2$`SSC self-renewal ~ initiation of differentiation`,]
undif  <- undif[,c(8,7,3)]
undif_wide<-spread(undif,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值undif_c1 <- undif[(undif$cluster)%in% "6",]
undif_wide <- as.data.frame(undif_wide)
rownames(undif_wide) <- undif_wide$gene
undif_wide <- undif_wide[,-1]
#######
dif <- mmc2_2_[(mmc2_2_$gene)%in%nitaharaDatasetS2$`Spermatogonial differentiation ~ meiotic entry`,]
dif  <- dif[,c(8,7,3)]
dif_wide<-spread(dif,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值dif_c1 <- dif[(dif$cluster)%in% "6",]
dif_wide <- as.data.frame(dif_wide)
rownames(dif_wide) <- dif_wide$gene
dif_wide <- dif_wide[,-1]
###########
MEIOi <- mmc2_2_[(mmc2_2_$gene)%in%nitaharaDatasetS2$`meiosis prophase I`,]
MEIOi  <- MEIOi[,c(8,7,3)]
MEIOi_wide<-spread(MEIOi,cluster,avg_logFC,fill = 0) #year为需要分解的变量，gdp为分解后的列的取值MEIOi_c1 <- MEIOi[(MEIOi$cluster)%in% "6",]
MEIOi_wide <- as.data.frame(MEIOi_wide)
rownames(MEIOi_wide) <- MEIOi_wide$gene
MEIOi_wide <- MEIOi_wide[,-1]




###############################
p1<-Heatmap(undif_wide,column_names_side = "top",column_title = "SSC self-renewal ~ initiation of differentiation",name = "avg_logFC")

p2<-Heatmap(dif_wide,column_names_side = "top",column_title = "Spermatogonial differentiation ~ meiotic entry",name = "avg_logFC")

p3<-Heatmap(MEIOi_wide,column_names_side = "top",column_title = "meiosis prophase I",name = "avg_logFC")


svg(file="SSC self-renewal ~ initiation of differentiation.svg", width=6, height=30)
p1 
dev.off()
svg(file="Spermatogonial differentiation ~ meiotic entry.svg", width=6, height=30)
p2
dev.off()
svg(file="meiosis prophase I.svg", width=6, height=30)
p3
dev.off()











################

undif_mat<- as.matrix(undif_wide)
undif_mat_scale = t(scale(t(undif_mat)))
p1<-Heatmap(undif_mat_scale,column_names_side = "top",column_title = "SSC self-renewal ~ initiation of differentiation",name = "normalized_avg_logFC")


dif_mat<- as.matrix(dif_wide)
dif_mat_scale = t(scale(t(dif_mat)))
p2<-Heatmap(dif_mat_scale,column_names_side = "top",column_title = "Spermatogonial differentiation ~ meiotic entry",name = "normalized_avg_logFC")



MEIOi_mat<- as.matrix(MEIOi_wide)
MEIOi_mat_scale = t(scale(t(MEIOi_mat)))
p3<-Heatmap(MEIOi_mat_scale,column_names_side = "top",column_title = "meiosis prophase I",name = "normalized_avg_logFC")



###############################

svg(file="SSC self-renewal ~ initiation of differentiation.svg", width=6, height=30)
p1 
dev.off()
svg(file="Spermatogonial differentiation ~ meiotic entry.svg", width=6, height=30)
p2
dev.off()
svg(file="meiosis prophase I.svg", width=6, height=30)
p3
dev.off()
p1+p2+p3
