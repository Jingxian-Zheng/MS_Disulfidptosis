Sys.setlocale('LC_CTYPE','English_United States.utf8')
options(stringsAsFactors=FALSE)
.libPaths(c(Sys.getenv('R_LIBS_USER'), .libPaths()))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(jsonlite))
root <- Sys.getenv('LEGACY_STUDY_ROOT',unset=file.path(Sys.getenv('REVISION_ROOT',unset=getwd()),'legacy_inputs'))
out <- file.path(root,'manuscript-Revised/revision_20261026')
logfile <- file.path(out,'audit/source_audit.txt')
sink(logfile,split=TRUE)
cat('RUNTIME\n'); print(sessionInfo())
write.csv(read.delim(file.path(root,'Data/\u4e8c\u786b\u4e2d\u6bd2\u76f8\u5173\u57fa\u56e0\u96c6.txt')), file.path(out,'tables/S01_original_gene_set_complete.csv'),row.names=FALSE)
gs <- read.delim(file.path(root,'Data/\u4e8c\u786b\u4e2d\u6bd2\u76f8\u5173\u57fa\u56e0\u96c6.txt')); print(table(gs[[1]]));cat('Unique background genes',length(unique(gs[[2]])),'\n')
markers<-scan(file.path(root,'Result1/disulfidptosis_markers.txt'),what='',quiet=TRUE)
write.csv(data.frame(gene=markers,original_set_member=markers %in% gs[[2]]),file.path(out,'tables/S01_selected_gene_provenance.csv'),row.names=FALSE)
for(i in c(3,4,5,6)){
 e<-new.env(); load(file.path(root,paste0('Result',i),'.RData'),envir=e)
 keys<-switch(as.character(i),'3'=c('auc_matrix','consensus_df','gbm_summary','GSE17048_sample_type','GSE32915_sample_type'), '4'=c('ssgsea_res','cibersort_res','valid_checkpoints','target_checkpoints','GSE17048_sample_type'), '5'=c('global_screening_res','predict_ic50'), '6'=c('sig_data','go_res','kegg_res'))
 for(n in keys){
  if(!exists(n,envir=e))next
  o<-e[[n]]
  if(is.matrix(o)||is.data.frame(o))write.csv(o,file.path(out,'audit',paste0('original_R',i,'_',n,'.csv')),row.names=TRUE)
  else if(is.character(o))writeLines(o,file.path(out,'audit',paste0('original_R',i,'_',n,'.txt')))
 }
 if(i==3){write.csv(data.frame(sample=e$train_groups$sample,y=e$y_train,old_DRSS=e$all_train_prob),file.path(out,'audit/original_DRSS.csv'),row.names=FALSE);saveRDS(e$fit_gbm,file.path(out,'machine_learning/original_GBM_model.rds'))}
 if(i==4)saveRDS(list(ssgsea=e$ssgsea_res,cibersort=e$cibersort_res,checkpoints=e$valid_checkpoints),file.path(out,'immune/original_immune_scores.rds'))
 if(i==6){cat('SIG DATA\n');print(head(e$sig_data));print(e$myPaths)}
 rm(e);invisible(gc())
}
cat('Reading referenced atlas\n')
atlasfile<-Sys.getenv('ATLAS_RDS')
suppressPackageStartupMessages(library(Seurat))
sce<-readRDS(atlasfile)
cat('ATLAS\n');print(sce);print(names(sce@meta.data));print(head(sce@meta.data));print(Assays(sce))
write.csv(sce@meta.data,file.path(out,'single_nucleus/nuclei_metadata.csv'),row.names=TRUE)
for(an in Assays(sce)){
 aa<-sce[[an]];cat('ASSAY',an,'class',class(aa),'\n'); print(slotNames(aa))
 for(sl in intersect(c('counts','data','scale.data'),slotNames(aa))){
  mm<-slot(aa,sl); cat('SLOT',sl,'DIM',dim(mm),'\n')
  if(inherits(mm,'sparseMatrix') && length(mm@x))cat('RANGE',range(mm@x),'fraction_integer',mean(abs(mm@x[seq_len(min(100000,length(mm@x)))]-round(mm@x[seq_len(min(100000,length(mm@x)))]))<1e-6),'\n')
 }
}
meta<-sce@meta.data
for(cn in names(meta)){if(length(unique(meta[[cn]]))<=100){cat('META',cn,'\n');print(table(meta[[cn]],useNA='always'))}}
sc<-which(meta$celltype=='SC');saveRDS(sce[,sc],file.path(out,'single_nucleus/stromal_original.rds'))
# Leave original atlas unchanged. Its exact read path is recorded for reproducibility.
write_json(list(atlas=atlasfile,atlas_size=file.info(atlasfile)$size,n_nuclei=ncol(sce),n_genes=nrow(sce),metadata_columns=names(meta)),file.path(out,'audit/atlas_source.json'),pretty=TRUE)
cat('AUDIT_COMPLETE\n');sink()
