
addMetaDataToSpatialData <- function(SeuratObject, newMetaData ){
  md = SeuratObject@meta.data
  if(length(grep(pattern = "SpotID", x = colnames(md))) == 0){
    md %>% rownames_to_column(var = "SpotID") %>% 
      separate(col = SpotID, into = c("Sample","NTtag"), sep = "_", remove = FALSE )
  }
  left_join(md,  newMetaData) %>%as.data.frame() -> new_md
  rownames(new_md) = new_md$SpotID
  new_md = new_md[rownames(md),]
  SeuratObject@meta.data = new_md
  
  return(SeuratObject)
} 


integrateSamples <-  function(SeuratObject, splitColumn = "orig.ident" ){
  ifnb.list <- SplitObject(SeuratObject, split.by = splitColumn)

  for(name in names(ifnb.list)){
    ifnb.list[[name]] %>%
      SCTransform(assay = "RNA",  
                  method = "glmGamPoi",
                  verbose = FALSE)   %>%
      FindVariableFeatures(xselection.method = "vst", nfeatures = 3000)->
      test
    ifnb.list[[name]] = test
  }
  
  features <- SelectIntegrationFeatures(object.list = ifnb.list)
  
  ifnb.list = PrepSCTIntegration(ifnb.list)
  major.anchors <- FindIntegrationAnchors(object.list = ifnb.list, 
                                          anchor.features = features, 
                                          normalization.method = "SCT",
                                          scale = TRUE)
  
  integrated.seuratObject <- IntegrateData(anchorset = major.anchors, normalization.method = "SCT")

  DefaultAssay(integrated.seuratObject) <- "integrated"
  
  return(integrated.seuratObject)
}

addMetaDataToSingleCellData <- function(SeuratObject, newMetaData){
  md = SeuratObject@meta.data
  if(length(grep(pattern = "cellID", x = colnames(md))) == 0){
    md %>% rownames_to_column(var = "cellID") %>% as.data.frame()->
      md
    rownames(md) = md$cellID 
  }
  left_join(md,  newMetaData) %>% as.data.frame() -> new_md
  rownames(new_md) = new_md$cellID
  new_md = new_md[rownames(md),]
  SeuratObject@meta.data = new_md
  
  return(SeuratObject)
} 


wssplot <- function(data, nc=15, seed=1234){
  wss <- (nrow(data)-1)*sum(apply(data,2,var))
  for (i in 2:nc){
    set.seed(seed)
    wss[i] <- sum(kmeans(data, centers=i)$withinss)}
  plot(1:nc, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
  wss
}



writeUMAPtoLoupe <- function(seuratObject, 
                             dir= "~/Box/bridgeSpider/spatialRNA/loupeFiles/UMAP",
                             preFile =  "intergrated"){
  md = seuratObject@meta.data
  
  md = as.data.frame(Embeddings(object = seuratObject[["umap"]])) %>%
    rownames_to_column(var = "SpotID") %>%
    inner_join(md)
  
  md %>% select(NTtag, Sample, UMAP_1,UMAP_2) %>%
    rename(X.Coordinate = "UMAP_1", Y.Coordinate = "UMAP_2", Barcode = "NTtag") ->
    md.umap
  
  for( slide in unique(md.umap$Sample)){
    UMAPfile = paste(dir,paste(preFile,slide,"UMAP.csv", sep = "."), sep = "/")
    md.umap %>%
      filter(Sample == slide) %>%
      select(Barcode,X.Coordinate ,Y.Coordinate)->
      md.umap.slide
    write_csv(x = md.umap.slide, file = UMAPfile)
  }
  
}

writeMetaDataToLoupe <- function(seuratObject,classes,preFile , 
                                 dir= "~/Box/bridgeSpider/spatialRNA/loupeFiles/UMAP"){
  md = seuratObject@meta.data
  
  
  for( slide in unique(md$Sample)){
    GBfile = paste(dir,paste(preFile,slide,"Graph-Based.csv", sep = "."), sep = "/")
    md %>%
      filter(Sample == slide) %>%
      select(NTtag, all_of(classes)) %>%
      rename(Barcode = "NTtag")->
      md.class
    write_csv(x = md.class, file = GBfile)
  }
  
}

addMetaDataToSeruatObject <-function(seuratObject, newMetaData){
  md.seurat = seuratObject@meta.data
  md.seurat = addMetaDataToSeruatMetaData(md.seurat,newMetaData )
  seuratObject@meta.data = md.seurat
  return(seuratObject)
}

addMetaDataToSeruatMetaData = function(md, newMetaData){
  if(length(intersect("cellID", colnames(md)) )== 0){
    md %>%rownames_to_column("cellID") ->md.temp  
  }else{
    md->md.temp  
  }
  if(length(intersect(colnames(md.temp), colnames(newMetaData)) )>0){
    md.temp2 = left_join(md.temp, newMetaData) %>% as.data.frame()
    rownames(md.temp2) = md.temp2$cellID
    md.temp2 = md.temp2[rownames(md), ]
    return(md.temp2)
  }
  else{
    print("No common columns to merge")
    return(md)
  }
  
  
}


GeneSetEnrichment <- function(gene2function, background, enrichedGenes){
  library(clusterProfiler)
  
  
  
}



getVSTexpression = function(count.df, 
                            geneColumn = "gene_id",
                            sampleColumn = "sample",
                            countColumn = "readCount",
                            md ,
                            gi, 
                            design = "Tissue"){
  
  
  count.df %>%
    dplyr::select(all_of(geneColumn), all_of(sampleColumn),all_of(countColumn)) %>%
    pivot_wider(names_from = all_of(sampleColumn),
                values_from = all_of(countColumn) ) %>%
    as.data.frame() %>%
    column_to_rownames(var = "gene_id")->
    CM
  
  
  CM = CM[gi[[all_of(geneColumn)]],
          md[[ all_of(sampleColumn)]]
  ]
  
  md$design = as.factor(md[[all_of(design)]])
  
  dds.tissue <- DESeqDataSetFromMatrix(countData = CM,
                                       colData = md,
                                       design = ~design)
  
  if(nrow(gi) >=1000){
    vst <- vst(dds.tissue)
  }else{
    vst = varianceStabilizingTransformation(dds.tissue)
  }
  normExpression = as.data.frame(assay(vst))
  return(normExpression)
}


getNormalisedExpression = function(count.df, 
                                   geneColumn = "gene_id",
                                   sampleColumn = "sample",
                                   countColumn = "readCount",
                                   md ,
                                   gi, 
                                   design = "Tissue"){
  
  
  count.df %>%
    dplyr::select(all_of(geneColumn), all_of(sampleColumn),all_of(countColumn)) %>%
    pivot_wider(names_from = all_of(sampleColumn),
                values_from = all_of(countColumn) ) %>%
    as.data.frame() %>%
    column_to_rownames(var = "gene_id")->
    CM
  
  
  CM = CM[gi[[all_of(geneColumn)]],
          md[[ all_of(sampleColumn)]]
  ]
  
  md$design = as.factor(md[[all_of(design)]])
  
  dds.tissue <- DESeqDataSetFromMatrix(countData = CM,
                                       colData = md,
                                       design = ~design)
  dds.tissue <- DESeq(dds.tissue)
  
  
  normCounts <- DESeq2::counts(dds.tissue, normalized=TRUE)
  
  normExpression = as.data.frame(assay(normCounts))
  return(normExpression)
}


getTPM = function(count.df, 
                  geneColumn = "gene_id",
                  sampleColumn = "sample",
                  countColumn = "readCount",
                  lengthColumn = "length",
                  md ,
                  gi, 
                  design = "Tissue"){
  
  
  count.df %>% group_by(all_of(sampleColumn)) %>%
    summarise(sampleSum = sum(all_of(countColumn))) %>%
    inner_join(count.df) 
    
    dplyr::select(all_of(geneColumn), all_of(sampleColumn),all_of(countColumn)) %>%
    pivot_wider(names_from = all_of(sampleColumn),
                values_from = all_of(countColumn) ) %>%
    as.data.frame() %>%
    column_to_rownames(var = "gene_id")->
    CM
  
  
  CM = CM[gi[[all_of(geneColumn)]],
          md[[ all_of(sampleColumn)]]
  ]
  
  md$design = as.factor(md[[all_of(design)]])
  
  dds.tissue <- DESeqDataSetFromMatrix(countData = CM,
                                       colData = md,
                                       design = ~design)
  dds.tissue <- DESeq(dds.tissue)
  
  
  normCounts <- DESeq2::counts(dds.tissue, normalized=TRUE)
  
  normExpression = as.data.frame(assay(normCounts))
  return(normExpression)
}





getDEseq2Object = function(count.df, 
                           geneColumn = "gene_id",
                           sampleColumn = "sample",
                           countColumn = "readCount",
                           md ,
                           gi, 
                           design = "Tissue"){
  
  
  count.df %>%
    dplyr::select(all_of(geneColumn), all_of(sampleColumn),all_of(countColumn)) %>%
    pivot_wider(names_from = all_of(sampleColumn),
                values_from = all_of(countColumn) ) %>%
    as.data.frame() %>%
    column_to_rownames(var = "gene_id")->
    CM
  
  
  CM = CM[gi[[all_of(geneColumn)]],
          md[[ all_of(sampleColumn)]]
  ]
  
  md$design = as.factor(md[[all_of(design)]])
  
  dds.tissue <- DESeqDataSetFromMatrix(countData = CM,
                                       colData = md,
                                       design = ~design)
  
  
  return(dds.tissue)
}




getOPLS = function(normExpression, 
                   nrOfPredictive = 1,
                   nrOfOrtho = NA,
                   md,
                   gi, 
                   design = "Tissue"){
  
  Parts.pls <- opls(t(normExpression),md[[design]] ,
                    predI = nrOfPredictive, orthoI = nrOfOrtho)
}

getPLS = function(normExpression, 
                  md,
                  gi, 
                  design = "Tissue"){
  Parts.pls <- opls(t(normExpression),md[[design]])
}




getMetaDataFromOPLS <- function(OPLS =  OPLS.parts,
                                md =sampleInfo.QC.merged.parts,
                                extension,sampleColumn = "sample")  {
  
  
  OPLS@scoreMN %>% as.data.frame() ->md.opls
  cn = paste(extension, colnames(md.opls), sep = "_")
  colnames(md.opls) = cn
  md.opls %>%
    rownames_to_column(var = all_of(sampleColumn)) ->
    md.opls 
  md = md %>% left_join(md.opls)
  return(md)
}  


getGeneInfoFromOPLS <- function(OPLS,md ,gi, 
                                extension ,
                                sampleColumn = "sample",
                                geneColumn = "gene_id", 
                                design = "Tissue")  {
  
  
  OPLS@scoreMN %>% as.data.frame() ->md.opls
  md.opls %>%
    rownames_to_column(var = all_of(sampleColumn)) ->
    md.opls 
  md %>% left_join(md.opls)%>%
    rename(design = all_of(design)) %>%
    group_by(design) %>%
    summarise(mean = mean(p1)) %>%
    mutate(dir = "neg" ) %>%
    mutate(dir = replace(x = dir, mean >=0 , values = "pos") ) %>%
    select(design,dir)  ->
    md.summary 
  
  OPLS@loadingMN %>% as.data.frame() -> loading
  
  
  loading %>% rownames_to_column(var = all_of(geneColumn)) %>%
    mutate(dir = "neg" ) %>%
    mutate(dir = replace(x = dir, p1 >=0 , values = "pos") ) %>%
    inner_join(md.summary) %>%
    select(all_of(geneColumn),design,p1) ->
    loading
  
  colnames(loading)[3] = paste(extension,colnames(loading)[3], sep = "_" )
  colnames(loading)[2] = all_of(extension)
  
  gi %>% left_join(loading) -> gi
  
  
  
  return(gi)
}

getGeneInfoFrom2DPLS <- function(OPLS,md ,gi, 
                                 extension ,
                                 sampleColumn = "sample",
                                 geneColumn = "gene_id", 
                                 design = "Tissue")  {
  
  
  OPLS@scoreMN %>% as.data.frame() ->md.opls
  md.opls %>%
    rownames_to_column(var = all_of(sampleColumn)) ->
    md.opls 
  md %>% left_join(md.opls)%>%
    rename(design = all_of(design)) %>%
    group_by(design) %>%
    summarise(mean_p1 = mean(p1),mean_p2 = mean(p2)) %>%
    mutate(dir_p1 = "neg" ) %>%
    mutate(dir_p1 = replace(x = dir_p1, mean_p1 >=0 , values = "pos") ) %>%
    mutate(dir_p2 = "neg" ) %>%
    mutate(dir_p2 = replace(x = dir_p2, mean_p2 >=0 , values = "pos") ) %>%
    select(design,dir_p1,dir_p2)  ->
    md.summary 
  
  OPLS@loadingMN %>% as.data.frame() -> loading
  
  
  loading %>% rownames_to_column(var = all_of(geneColumn)) %>%
    mutate(dir_p1 = "neg" ) %>%
    mutate(dir_p1 = replace(x = dir_p1, p1 >=0 , values = "pos") ) %>%
    mutate(dir_p2 = "neg" ) %>%
    mutate(dir_p2 = replace(x = dir_p2, p2 >=0 , values = "pos") ) %>%
    left_join(md.summary) %>%
    select(all_of(geneColumn),dir_p1,dir_p2,design,p1,p2) ->
    loading
  
  
  
  colnames(loading)[5] = paste(extension,colnames(loading)[5], sep = "_" )
  colnames(loading)[6] = paste(extension,colnames(loading)[6], sep = "_" )
  colnames(loading)[4] = all_of(extension)
  
  return(loading)
}



