library(tidyverse)
library(Seurat)
library(gtools)
library(foreach)
library(doParallel)


#* This function is to do DEG calling between any one of 
#* treatment group over wildtype/untreat group for all clusters in  
#* input seurat object. 
#* The grouping factor is a column in seuart object
#* meta.data and more than one groups/labels in that column can be used 
#* as wildtype group. 
#* 
#* 


DEGcallingPerCluster<-function(seuratObj,
                               control_group=c("wt"),
                               group_name="condition",
                               idents_name,n_cores=1){
  registerDoParallel(n_cores)
  
  # get all cluster
  Idents(seuratObj)<- idents_name
  clusters <- mixedsort(unique(Idents(seuratObj)))
  
  # get all groups 
  groups <- unique(seuratObj$sample)
  no_wt_groups <- setdiff(groups, control_group)
  
  # split the job for one group over control group
  deg_all <- NULL
  for (grp in no_wt_groups){
    # select cells that are in target and control group
    select_cells <- which(seuratObj@meta.data[,group_name]%in%c(control_group,grp))
    s <- seuratObj[,select_cells]
    
    # DEG calling per cluster
    Idents(s) <- idents_name
    deg <- foreach(cl=clusters, .combine = rbind)%dopar%{
      tmp<-NULL
      try({
        tmp <- FindMarkers(object = s,
                           ident.1 = grp,
                           group.by = group_name,
                           subset.ident = cl)
        tmp$gene <- row.names(tmp)
        tmp$cluster <- cl
        tmp$avg_log2FC <- round(tmp$avg_log2FC, 3)
        tmp$p_val<- signif(tmp$p_val, 4)
        tmp$p_val_adj <- signif(tmp$p_val_adj, 4)
      })
      return(tmp)
    }
    
    if (!is.null(deg)) {
      deg$comparison <- paste0(grp, "_vs_", 
                               paste(control_group,collapse=";"))
      deg_all <- rbind(deg_all, deg)
    }
  }
  stopImplicitCluster()
  return(deg_all)
}
