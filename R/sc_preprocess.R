#' Preprocess a single-cell dataset for Scissor
#'
#' Moves TCR/BCR genes into a separate \code{AIR} assay, filters cells by QC
#' metrics, selects highly variable genes, optionally corrects batch effects
#' with Harmony, and builds the cell-cell graph, clusters and UMAP.
#'
#' @param object A Seurat (v5) object.
#' @param assay Assay to process. Default \code{"RNA"}.
#' @param filter_TCR_BCR Logical. Move TCR/BCR genes (Ensembl biotypes \code{IG_*} / \code{TR_*}, plus TRBV11-2) into an \code{AIR} assay.
#' @param keep_c_genes Logical. Keep constant-region genes (e.g. TRAC, IGHM) in \code{assay}.
#' @param plot_QC_metrics Logical. Draw violin plots of the QC metrics.
#' @param threshold_min_nCount,threshold_max_nCount Range of UMI counts per cell to keep.
#' @param threshold_min_nFeature,threshold_max_nFeature Range of detected genes per cell to keep.
#' @param threshold_min_percentMT,threshold_max_percentMT Range of mitochondrial percentage to keep.
#' @param correct_batch Logical. Run Harmony batch correction.
#' @param group.by.vars Metadata column(s) giving the batch.
#' @param threshold_batch Drop batches with fewer cells than this when selecting HVGs.
#' @param dims.use,harmony_theta,harmony_lambda Passed to \code{RunHarmony}.
#' @param nfeatures_batch Number of HVGs per batch.
#' @param nfeatures Number of HVGs in total.
#' @param hvg_custom Extra genes to add to the HVGs.
#' @param vars.to.regress Variables to regress out in \code{ScaleData}.
#' @param rm_cycling Logical. Regress out the cell-cycle scores.
#' @param npcs Number of PCs to compute.
#' @param k.param Number of neighbours for \code{FindNeighbors}.
#' @param dims Dimensions used for the neighbour graph and UMAP.
#' @param resolution Clustering resolution.
#'
#' @return The processed Seurat object.
#' @export
sc_preprocess = function(object, assay='RNA',
                      filter_TCR_BCR=T,keep_c_genes=T,
                      plot_QC_metrics=T, threshold_min_nCount=0, threshold_max_nCount=Inf, threshold_min_nFeature=0, threshold_max_nFeature=Inf, threshold_min_percentMT=0, threshold_max_percentMT=100,
                      correct_batch=T, group.by.vars = 'batch', threshold_batch=0, dims.use=1:50, harmony_theta=2, harmony_lambda=1,
                      nfeatures_batch=5000, nfeatures=4000, hvg_custom=NULL,
                      vars.to.regress=c(), rm_cycling = F,
                      npcs=50,
                      k.param=20,
                      dims=1:20, resolution=.5){

    DefaultAssay(object) = assay

    if(inherits(object[[assay]], "Assay5")){
        object[[assay]] = JoinLayers(object[[assay]])
    }

    ### Filter TCR & BCR genes ###
    if(filter_TCR_BCR){

        remove_genes = .remove_genes # internal TCR/BCR gene list, defined at the bottom of this file

        if(keep_c_genes){
            # Keep common region genes, not associated with clonotypes
            c_genes = c("TRAC","TRDC","IGKC","IGHM","IGHD","IGHE","TRBC1","TRBC2","TRGC2","TRGC1","IGHA1","IGHA2","IGHG1","IGHG2","IGHG3","IGHG4","IGLC1","IGLC2","IGLC3","IGLC7")
            remove_genes = remove_genes[-which(remove_genes %in% c_genes)]

        }

        # Only touch genes still present in the target assay (handles re-running
        # on an object that already had some/all TCR/BCR genes moved to "air")
        genes_to_move = intersect(remove_genes, rownames(object[[assay]]))

        if(length(genes_to_move) > 0){
            new_air = GetAssayData(object, assay = assay, layer = "counts")[genes_to_move, , drop = FALSE]

            if("AIR" %in% Assays(object)){
                message(length(genes_to_move), " new TCR/BCR gene(s) found in '", assay, "' assay; merging into existing **air** assay!")
                old_air = GetAssayData(object, assay = "AIR", layer = "counts")
                combined_air = rbind(old_air, new_air[, colnames(old_air), drop = FALSE])
                object[["AIR"]] = CreateAssayObject(counts = combined_air)
            } else {
                message("The TCR & BCR genes are stored into AIR assay!")
                object[["AIR"]] = CreateAssayObject(counts = new_air)
            }
            object[[assay]] = subset(object[[assay]], features = setdiff(rownames(object[[assay]]), genes_to_move))
        } else {
            message("No new TCR/BCR genes found in '", assay, "' assay; AIR assay is already up to date.")
        }

    }

    ### QC metrics for the target assay ###
    nCount_col = paste0("nCount_", assay)
    nFeature_col = paste0("nFeature_", assay)
    counts = GetAssayData(object, assay = assay, layer = "counts")
    object[[nCount_col]] = colSums(counts)
    object[[nFeature_col]] = colSums(counts > 0)
    rm(counts)

    object[["percent.mt"]] = PercentageFeatureSet(object, pattern = "^MT-", assay = assay)
    print(dim(object))

    ### Plot QC metrics ###
    if(plot_QC_metrics){
        options(repr.plot.width=15,repr.plot.height=10)
        object$QC_check = 'SeuratObject'
        Idents(object) = 'QC_check'
        p1 = VlnPlot(object,features = c(nCount_col,nFeature_col,'percent.mt'),pt.size = 0,ncol = 3)
        print(p1)
    }

    meta = object@meta.data
    keep_cells = rownames(meta)[which((meta$percent.mt<=threshold_max_percentMT) & (meta$percent.mt>=threshold_min_percentMT) &
                                (meta[[nFeature_col]]<=threshold_max_nFeature) & (meta[[nFeature_col]]>=threshold_min_nFeature) &
                                (meta[[nCount_col]]<=threshold_max_nCount) & (meta[[nCount_col]]>=threshold_min_nCount))]
    object = subset(object, cells = keep_cells)

    object = NormalizeData(object)
    object = FindVariableFeatures(object, nfeatures = nfeatures)
    
    ### Collect HVGs ###
    if(correct_batch){
        obj.list = SplitObject(object, split.by = group.by.vars[1])
        obj.list = obj.list[sapply(obj.list, ncol) >= threshold_batch]

        obj.list = lapply(obj.list,FUN = function(x){
            x = NormalizeData(x) %>% FindVariableFeatures(nfeatures = nfeatures_batch)
            hvg_x = VariableFeatures(x)
            gene.pattern = c("MALAT1", "^MT-", "^RPL", "^RPS","^MTRNR")
            hvg_x = hvg_x[!hvg_x %in% grep(paste0(gene.pattern, collapse = "|"), hvg_x, value = T)]
            if(filter_TCR_BCR){
                hvg_x = hvg_x[!hvg_x %in% remove_genes]
            }
            VariableFeatures(x) = hvg_x
            return(x)})
        hvg = SelectIntegrationFeatures(obj.list,nfeatures = nfeatures)
    } 
    else{
        hvg = VariableFeatures(object)
        gene.pattern = c("MALAT1", "^MT-", "^RPL", "^RPS","^MTRNR")
        hvg = hvg[!hvg %in% grep(paste0(gene.pattern, collapse = "|"), hvg, value = T)]
        if(filter_TCR_BCR){
                hvg = hvg[!hvg %in% remove_genes]
            }
    }


    if(!is.null(hvg_custom)){
        hvg = c(hvg,hvg_custom[hvg_custom %in% rownames(object)]) %>% unique
    }
    VariableFeatures(object) = hvg
    object = NormalizeData(object) 
    
    ### Cycling genes ###
    s.genes = Seurat::cc.genes$s.genes
    g2m.genes = Seurat::cc.genes$g2m.genes
    object = CellCycleScoring(object, s.features = s.genes, g2m.features = g2m.genes, set.ident = T)
    if(rm_cycling){
        vars.to.regress = c(vars.to.regress,"S.Score","G2M.Score")
    }

    object = ScaleData(object,features = hvg, vars.to.regress=vars.to.regress)
    object = RunPCA(object,npcs = npcs,features = hvg)

    ### Select nPC ###
    options(repr.plot.width=7,repr.plot.height=6)
    p2 = ElbowPlot(object, ndims = npcs)
    print(p2)

    ### Batch correction with harmony ###
    if(correct_batch){
        object = RunHarmony(object,dims.use = dims.use,reduction.use = 'pca',group.by.vars = group.by.vars,theta = harmony_theta,lambda=harmony_lambda)
        
        options(repr.plot.width=7,repr.plot.height=6)
        p3 = ElbowPlot(object,ndims = npcs,reduction = 'harmony')
        print(p3)
        
        object = FindNeighbors(object,reduction = 'harmony',dims = dims,k.param=k.param) %>% FindClusters(resolution = resolution)
        object = RunUMAP(object,reduction = 'harmony',dims = dims)
    }
    else{
        object = FindNeighbors(object,reduction = 'pca',dims = dims,k.param=k.param) %>% FindClusters(resolution = resolution)
        object = RunUMAP(object,reduction = 'pca',dims = dims)
    }

    return(object)
}

# TCR/BCR genes
.remove_genes = c(
    "AC134981.1","AC135068.2","AC135068.8","AC135068.9","AC244452.3","IGHA1","IGHA2","IGHD",
    "IGHD1-1","IGHD1-14","IGHD1-20","IGHD1-26","IGHD1-7","IGHD1OR15-1A","IGHD1OR15-1B","IGHD2-15",
    "IGHD2-2","IGHD2-21","IGHD2-8","IGHD2OR15-2A","IGHD2OR15-2B","IGHD3-10","IGHD3-16","IGHD3-22",
    "IGHD3-3","IGHD3-9","IGHD3OR15-3A","IGHD3OR15-3B","IGHD4-11","IGHD4-17","IGHD4-23","IGHD4-4",
    "IGHD4OR15-4A","IGHD4OR15-4B","IGHD5-12","IGHD5-18","IGHD5-24","IGHD5-5","IGHD5OR15-5A","IGHD5OR15-5B",
    "IGHD6-13","IGHD6-19","IGHD6-25","IGHD6-6","IGHD7-27","IGHE","IGHEP1","IGHEP2",
    "IGHG1","IGHG2","IGHG3","IGHG4","IGHGP","IGHJ1","IGHJ1P","IGHJ2",
    "IGHJ2P","IGHJ3","IGHJ3P","IGHJ4","IGHJ5","IGHJ6","IGHM","IGHV1-12",
    "IGHV1-14","IGHV1-17","IGHV1-18","IGHV1-2","IGHV1-24","IGHV1-3","IGHV1-45","IGHV1-46",
    "IGHV1-58","IGHV1-67","IGHV1-68","IGHV1-69","IGHV1-69-2","IGHV1-69D","IGHV1OR15-1","IGHV1OR15-2",
    "IGHV1OR15-3","IGHV1OR15-4","IGHV1OR15-6","IGHV1OR15-9","IGHV1OR16-1","IGHV1OR16-2","IGHV1OR16-3","IGHV1OR16-4",
    "IGHV1OR21-1","IGHV2-26","IGHV2-5","IGHV2-70","IGHV2-70D","IGHV2OR16-5","IGHV3-11","IGHV3-13",
    "IGHV3-15","IGHV3-16","IGHV3-19","IGHV3-20","IGHV3-21","IGHV3-22","IGHV3-23","IGHV3-25",
    "IGHV3-29","IGHV3-30","IGHV3-30-2","IGHV3-32","IGHV3-33","IGHV3-33-2","IGHV3-35","IGHV3-36",
    "IGHV3-37","IGHV3-38","IGHV3-41","IGHV3-42","IGHV3-43","IGHV3-47","IGHV3-48","IGHV3-49",
    "IGHV3-50","IGHV3-52","IGHV3-53","IGHV3-54","IGHV3-57","IGHV3-6","IGHV3-60","IGHV3-62",
    "IGHV3-63","IGHV3-64","IGHV3-64D","IGHV3-65","IGHV3-66","IGHV3-69-1","IGHV3-7","IGHV3-71",
    "IGHV3-72","IGHV3-73","IGHV3-74","IGHV3-75","IGHV3-76","IGHV3-79","IGHV3OR15-7","IGHV3OR16-10",
    "IGHV3OR16-11","IGHV3OR16-12","IGHV3OR16-13","IGHV3OR16-15","IGHV3OR16-6","IGHV3OR16-7","IGHV3OR16-8","IGHV3OR16-9",
    "IGHV4-28","IGHV4-31","IGHV4-34","IGHV4-39","IGHV4-4","IGHV4-55","IGHV4-59","IGHV4-61",
    "IGHV4-80","IGHV4OR15-8","IGHV5-10-1","IGHV5-51","IGHV5-78","IGHV6-1","IGHV7-27","IGHV7-34-1",
    "IGHV7-4-1","IGHV7-40","IGHV7-56","IGHV7-81","IGHVII-1-1","IGHVII-15-1","IGHVII-22-1","IGHVII-26-2",
    "IGHVII-28-1","IGHVII-30-1","IGHVII-30-21","IGHVII-33-1","IGHVII-40-1","IGHVII-43-1","IGHVII-44-2","IGHVII-46-1",
    "IGHVII-49-1","IGHVII-51-2","IGHVII-53-1","IGHVII-60-1","IGHVII-62-1","IGHVII-65-1","IGHVII-67-1","IGHVII-74-1",
    "IGHVII-78-1","IGHVIII-11-1","IGHVIII-13-1","IGHVIII-16-1","IGHVIII-2-1","IGHVIII-22-2","IGHVIII-25-1","IGHVIII-26-1",
    "IGHVIII-38-1","IGHVIII-44","IGHVIII-47-1","IGHVIII-5-1","IGHVIII-5-2","IGHVIII-51-1","IGHVIII-67-2","IGHVIII-67-3",
    "IGHVIII-67-4","IGHVIII-76-1","IGHVIII-82","IGHVIV-44-1","IGKC","IGKJ1","IGKJ2","IGKJ3",
    "IGKJ4","IGKJ5","IGKV1-12","IGKV1-13","IGKV1-16","IGKV1-17","IGKV1-22","IGKV1-27",
    "IGKV1-32","IGKV1-33","IGKV1-35","IGKV1-37","IGKV1-39","IGKV1-5","IGKV1-6","IGKV1-8",
    "IGKV1-9","IGKV1D-12","IGKV1D-13","IGKV1D-16","IGKV1D-17","IGKV1D-22","IGKV1D-27","IGKV1D-32",
    "IGKV1D-33","IGKV1D-35","IGKV1D-37","IGKV1D-39","IGKV1D-42","IGKV1D-43","IGKV1D-8","IGKV1OR-2",
    "IGKV1OR-3","IGKV1OR1-1","IGKV1OR10-1","IGKV1OR2-1","IGKV1OR2-108","IGKV1OR2-11","IGKV1OR2-118","IGKV1OR2-3",
    "IGKV1OR2-6","IGKV1OR2-9","IGKV1OR22-1","IGKV1OR22-5","IGKV1OR9-1","IGKV1OR9-2","IGKV2-10","IGKV2-14",
    "IGKV2-18","IGKV2-19","IGKV2-23","IGKV2-24","IGKV2-26","IGKV2-28","IGKV2-29","IGKV2-30",
    "IGKV2-36","IGKV2-38","IGKV2-4","IGKV2-40","IGKV2D-10","IGKV2D-14","IGKV2D-18","IGKV2D-19",
    "IGKV2D-23","IGKV2D-24","IGKV2D-26","IGKV2D-28","IGKV2D-29","IGKV2D-30","IGKV2D-36","IGKV2D-38",
    "IGKV2D-40","IGKV2OR2-1","IGKV2OR2-10","IGKV2OR2-2","IGKV2OR2-7","IGKV2OR2-7D","IGKV2OR2-8","IGKV2OR22-3",
    "IGKV2OR22-4","IGKV3-11","IGKV3-15","IGKV3-20","IGKV3-25","IGKV3-31","IGKV3-34","IGKV3-7",
    "IGKV3D-11","IGKV3D-15","IGKV3D-20","IGKV3D-25","IGKV3D-31","IGKV3D-34","IGKV3D-7","IGKV3OR2-268",
    "IGKV3OR2-5","IGKV3OR22-2","IGKV4-1","IGKV5-2","IGKV6-21","IGKV6D-21","IGKV6D-41","IGKV7-3",
    "IGLC1","IGLC2","IGLC3","IGLC4","IGLC5","IGLC6","IGLC7","IGLCOR22-1",
    "IGLCOR22-2","IGLJ1","IGLJ2","IGLJ3","IGLJ4","IGLJ5","IGLJ6","IGLJ7",
    "IGLJCOR18","IGLV1-36","IGLV1-40","IGLV1-41","IGLV1-44","IGLV1-47","IGLV1-50","IGLV1-51",
    "IGLV1-62","IGLV10-54","IGLV10-67","IGLV11-55","IGLV2-11","IGLV2-14","IGLV2-18","IGLV2-23",
    "IGLV2-28","IGLV2-33","IGLV2-34","IGLV2-5","IGLV2-8","IGLV3-1","IGLV3-10","IGLV3-12",
    "IGLV3-13","IGLV3-15","IGLV3-16","IGLV3-17","IGLV3-19","IGLV3-2","IGLV3-21","IGLV3-22",
    "IGLV3-24","IGLV3-25","IGLV3-26","IGLV3-27","IGLV3-29","IGLV3-30","IGLV3-31","IGLV3-32",
    "IGLV3-4","IGLV3-6","IGLV3-7","IGLV3-9","IGLV4-3","IGLV4-60","IGLV4-69","IGLV5-37",
    "IGLV5-45","IGLV5-48","IGLV5-52","IGLV6-57","IGLV7-35","IGLV7-43","IGLV7-46","IGLV8-61",
    "IGLV8OR8-1","IGLV9-49","IGLVI-20","IGLVI-38","IGLVI-42","IGLVI-56","IGLVI-63","IGLVI-68",
    "IGLVI-70","IGLVIV-53","IGLVIV-59","IGLVIV-64","IGLVIV-65","IGLVIV-66-1","IGLVIVOR22-1","IGLVIVOR22-2",
    "IGLVV-58","IGLVV-66","IGLVVI-22-1","IGLVVI-25-1","IGLVVII-41-1","TRAC","TRAJ1","TRAJ10",
    "TRAJ11","TRAJ12","TRAJ13","TRAJ14","TRAJ16","TRAJ17","TRAJ18","TRAJ19",
    "TRAJ2","TRAJ20","TRAJ21","TRAJ22","TRAJ23","TRAJ24","TRAJ25","TRAJ26",
    "TRAJ27","TRAJ28","TRAJ29","TRAJ3","TRAJ30","TRAJ31","TRAJ32","TRAJ33",
    "TRAJ34","TRAJ35","TRAJ36","TRAJ37","TRAJ38","TRAJ39","TRAJ4","TRAJ40",
    "TRAJ41","TRAJ42","TRAJ43","TRAJ44","TRAJ45","TRAJ46","TRAJ47","TRAJ48",
    "TRAJ49","TRAJ5","TRAJ50","TRAJ51","TRAJ52","TRAJ53","TRAJ54","TRAJ55",
    "TRAJ56","TRAJ57","TRAJ58","TRAJ59","TRAJ6","TRAJ60","TRAJ61","TRAJ7",
    "TRAJ8","TRAJ9","TRAV1-1","TRAV1-2","TRAV10","TRAV11","TRAV12-1","TRAV12-2",
    "TRAV12-3","TRAV13-1","TRAV13-2","TRAV14DV4","TRAV15","TRAV16","TRAV17","TRAV18",
    "TRAV19","TRAV2","TRAV20","TRAV21","TRAV22","TRAV23DV6","TRAV24","TRAV25",
    "TRAV26-1","TRAV26-2","TRAV27","TRAV28","TRAV29DV5","TRAV3","TRAV30","TRAV31",
    "TRAV32","TRAV33","TRAV34","TRAV35","TRAV36DV7","TRAV37","TRAV38-1","TRAV38-2DV8",
    "TRAV39","TRAV4","TRAV40","TRAV41","TRAV5","TRAV6","TRAV7","TRAV8-1",
    "TRAV8-2","TRAV8-3","TRAV8-4","TRAV8-5","TRAV8-6","TRAV8-7","TRAV9-1","TRAV9-2",
    "TRBC1","TRBC2","TRBD1","TRBJ1-1","TRBJ1-2","TRBJ1-3","TRBJ1-4","TRBJ1-5",
    "TRBJ1-6","TRBJ2-1","TRBJ2-2","TRBJ2-2P","TRBJ2-3","TRBJ2-4","TRBJ2-5","TRBJ2-6",
    "TRBJ2-7","TRBV1","TRBV10-1","TRBV10-2","TRBV10-3","TRBV11-1","TRBV11-2","TRBV11-3",
    "TRBV12-1","TRBV12-2","TRBV12-3","TRBV12-4","TRBV12-5","TRBV13","TRBV14","TRBV15",
    "TRBV16","TRBV17","TRBV18","TRBV19","TRBV2","TRBV20-1","TRBV20OR9-2","TRBV21-1",
    "TRBV21OR9-2","TRBV22-1","TRBV22OR9-2","TRBV23-1","TRBV23OR9-2","TRBV24-1","TRBV24OR9-2","TRBV25-1",
    "TRBV25OR9-2","TRBV26","TRBV26OR9-2","TRBV27","TRBV28","TRBV29-1","TRBV29OR9-2","TRBV3-1",
    "TRBV30","TRBV4-1","TRBV4-2","TRBV5-1","TRBV5-2","TRBV5-3","TRBV5-4","TRBV5-5",
    "TRBV5-6","TRBV5-7","TRBV6-1","TRBV6-2","TRBV6-4","TRBV6-5","TRBV6-6","TRBV6-7",
    "TRBV6-8","TRBV7-1","TRBV7-2","TRBV7-3","TRBV7-4","TRBV7-5","TRBV7-6","TRBV7-7",
    "TRBV7-9","TRBV8-1","TRBV8-2","TRBV9","TRBVA","TRBVB","TRDC","TRDD1",
    "TRDD2","TRDD3","TRDJ1","TRDJ2","TRDJ3","TRDJ4","TRDV1","TRDV2",
    "TRDV3","TRGC1","TRGC2","TRGJ1","TRGJ2","TRGJP","TRGJP1","TRGJP2",
    "TRGV1","TRGV10","TRGV11","TRGV2","TRGV3","TRGV4","TRGV5","TRGV5P",
    "TRGV6","TRGV7","TRGV8","TRGV9","TRGVA","TRGVB"
)

# Internal QC for both single-cell and Bulk expression matrices
.qc_on_mtx = function(mtx,
                     row_min_cnt = 1, row_expr_cutoff = 0,
                     col_min_cnt = 200, col_expr_cutoff = 0){
    dims = c(0,0)
    counter = 0
    print(paste0('The original matrix is in shape of (',dim(mtx)[1],', ',dim(mtx)[2],').'))
    while(!(dims[1] == dim(mtx)[1] & dims[2] == dim(mtx)[2])){
        print(paste0('Round',counter))
        dims = dim(mtx)
        mtx = mtx[apply(mtx,1,function(x){sum(x>row_expr_cutoff)}) >= row_min_cnt,]
        mtx = mtx[,apply(mtx,2,function(x){sum(x>col_expr_cutoff)}) >= col_min_cnt]
        counter = counter + 1
    }
    print(paste0('The final matrix is in shape of (',dims[1],', ',dims[2],').'))
    return(mtx)
}