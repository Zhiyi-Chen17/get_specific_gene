library(tidyverse)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
options(timeout = 4000000)









?function()
# 编写函数  
getGeneInfoByGO <- function(go_id, mart, custom_colname = "spermatid_differentiation") {  
  # 从GO ID获取基因Entrez ID  
  gene_ids <- get(go_id, org.Mm.egGO2ALLEGS)  
  
  # 检查是否获取到基因ID  
  if (length(gene_ids) == 0) {  
    stop("No gene IDs found for GO ID: ", go_id)  
  }  
  
  # 使用biomaRt获取基因信息  
  option_info <- getBM(  
    attributes = c("mgi_symbol"),  # 这里可以根据需要修改或添加更多属性  
    filters = "entrezgene_id",  
    values = gene_ids,  
    mart = mart  
  )  
  
  # 设置列名  
  colnames(option_info) <- custom_colname  
  
  # 返回结果  
  return(option_info)  
}  



# 调用函数  

# 查看结果  



gene_info_spermatid_differentiation                                          <- getGeneInfoByGO("GO:0048515", mart, custom_colname = "spermatid_differentiation                                            ")                                                                                                         
gene_info_spermatid_cytoplasm_removal_during_spermiation_of_flagellated_sperm<- getGeneInfoByGO("GO:0160087", mart, custom_colname = "spermatid_cytoplasm_removal_during_spermiation_of_flagellated_sperm  ")                                                                                                             
gene_info_spermatid_development                                              <- getGeneInfoByGO("GO:0007286", mart, custom_colname = "spermatid_development                                                ")                                                                                                                        
gene_info_spermatid_nucleus_differentiation                                  <- getGeneInfoByGO("GO:0007289", mart, custom_colname = "spermatid_nucleus_differentiation                                    ")                                                                                                            
gene_info_spermatid_nucleus_elongation                                       <- getGeneInfoByGO("GO:0007290", mart, custom_colname = "spermatid_nucleus_elongation                                         ")                                                                                                      
gene_info_regulation_of_spermatid_nuclear_differentiation                    <- getGeneInfoByGO("GO:0045700", mart, custom_colname = "regulation_of_spermatid_nuclear_differentiation                      ")                                                                            
gene_info_negative_regulation_of_spermatid_nuclear_differentiation           <- getGeneInfoByGO("GO:0045701", mart, custom_colname = "negative_regulation_of_spermatid_nuclear_differentiation             ")                                                                                     
gene_info_positive_regulation_of_spermatid_nuclear_differentiation           <- getGeneInfoByGO("GO:0045702", mart, custom_colname = "positive_regulation_of_spermatid_nuclear_differentiation             ")                                                                                     
gene_info_spermatocyte_division                                              <- getGeneInfoByGO("GO:0048137", mart, custom_colname = "spermatocyte_division                                                ")                                                  
gene_info_sperm_DNA_condensation                                             <- getGeneInfoByGO("GO:0035092", mart, custom_colname = "sperm_DNA_condensation                                               ")                                                   
gene_info_concave_side_of_sperm_head                                         <- getGeneInfoByGO("GO:0061830", mart, custom_colname = "concave_side_of_sperm_head                                           ")                                                       
gene_info_apical_ectoplasmic_specialization                                  <- getGeneInfoByGO("GO:0061831", mart, custom_colname = "apical_ectoplasmic_specialization                                    ")                                                              
gene_info_sperm_head                                                         <- getGeneInfoByGO("GO:0061827", mart, custom_colname = "sperm_head                                                           ")                                          
gene_info_apical_tubulobulbar_complex                                        <- getGeneInfoByGO("GO:0061828", mart, custom_colname = "apical_tubulobulbar_complex                                          ")                                                        
gene_info_ciliary_cap                                                        <- getGeneInfoByGO("GO:0061822", mart, custom_colname = "ciliary_cap                                                          ")                                          
gene_info_ring_centriole                                                     <- getGeneInfoByGO("GO:0061823", mart, custom_colname = "ring_centriole                                                       ")                                           
gene_info_tubulobulbar_complex                                               <- getGeneInfoByGO("GO:0036284", mart, custom_colname = "tubulobulbar_complex                                                 ")                                                 
gene_info_sperm_flagellum                                                    <- getGeneInfoByGO("GO:0036126", mart, custom_colname = "sperm_flagellum                                                      ")                                            
gene_info_acroblast                                                          <- getGeneInfoByGO("GO:0036063", mart, custom_colname = "acroblast                                                            ")                                          
gene_info_sperm_cytoplasmic_droplet                                          <- getGeneInfoByGO("GO:0097598", mart, custom_colname = "sperm_cytoplasmic_droplet                                            ")                                                      
gene_info_Nebenkern                                                          <- getGeneInfoByGO("GO:0016006", mart, custom_colname = "Nebenkern                                                            ")                                          
gene_info_manchette                                                          <- getGeneInfoByGO("GO:0002177", mart, custom_colname = "manchette                                                            ")                                          
gene_info_Nebenkern_assembly                                                 <- getGeneInfoByGO("GO:0007287", mart, custom_colname = "Nebenkern_assembly                                                   ")                                               
gene_info_sperm_individualization                                            <- getGeneInfoByGO("GO:0007291", mart, custom_colname = "sperm_individualization                                              ")                                                    
gene_info_acrosome_assembly                                                  <- getGeneInfoByGO("GO:0001675", mart, custom_colname = "acrosome_assembly                                                    ")    



x<-plyr::rbind.fill(gene_info_spermatid_differentiation                              
                    ,gene_info_spermatid_nucleus_elongation
                    ,gene_info_spermatid_development  
                            
                    ,gene_info_spermatocyte_division                                              
                    ,gene_info_sperm_DNA_condensation                                             
                    ,gene_info_concave_side_of_sperm_head                                         
                    ,gene_info_apical_ectoplasmic_specialization                                  
                    ,gene_info_sperm_head                                                         
                    ,gene_info_apical_tubulobulbar_complex                                        
                                                                         
                                                                     
                    ,gene_info_tubulobulbar_complex                                               
                    ,gene_info_sperm_flagellum                                                    
                                                                         
                    ,gene_info_sperm_cytoplasmic_droplet                                          
                                                                        
                    ,gene_info_manchette                                                          
                                                                  
                    ,gene_info_sperm_individualization                                            
                    ,gene_info_acrosome_assembly    )

write.csv(x,file="spermatid,gene_set.csv")







GO:0048515	spermatid differentiation
GO:0160087	spermatid cytoplasm removal during spermiation of flagellated sperm
GO:0007286	spermatid development
GO:0007289	spermatid nucleus differentiation
GO:0007290	spermatid nucleus elongation
GO:0045700	regulation of spermatid nuclear differentiation
GO:0045701	negative regulation of spermatid nuclear differentiation
GO:0045702	positive regulation of spermatid nuclear differentiation
GO:0048137	spermatocyte division
GO:0035092	sperm DNA condensation
GO:0061830	concave side of sperm head
GO:0061831	apical ectoplasmic specialization
GO:0061827	sperm head
GO:0061828	apical tubulobulbar complex
GO:0061822	ciliary cap
GO:0061823	ring centriole
GO:0036284	tubulobulbar complex
GO:0036126	sperm flagellum
GO:0036063	acroblast
GO:0097598	sperm cytoplasmic droplet
GO:0016006	Nebenkern
GO:0002177	manchette
GO:0007287	Nebenkern assembly
GO:0007291	sperm individualization
GO:0001675	acrosome assembly

bind_cols(data1, data3, id = NULL)   
