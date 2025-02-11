#### set of maftools function to use in the app.
### #WRITTEN BY FARID AHADLI (farid.ahadli@outlook.com)


##### subsetting MAF
#load(file = "turkish_subtype_added_pcgr_added.RData")



all_pathways = unique(read.delim(file = system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools"))$Pathway)
load(file = "tcga_maf_data_subtype_added.RData")
hugo.to.ms = system.file("extdata", "hugo_to_mutSigSymbol.txt.gz", 
                         package = "maftools")
hugo.to.ms = as.data.frame(data.table::fread(file = hugo.to.ms, sep = "\t", 
                                             stringsAsFactors = FALSE))
#load(file = "turkish_IHC_added_MAF.RData") do not load
#below will come from input
#############
#turkish_maf_genes = c(unique(c(turkish_ihc_added@data$Hugo_Symbol,turkish_ihc_added@maf.silent$Hugo_Symbol)),"all")
#turkish_maf_patients = c(unique(turkish_ihc_added@clinical.data$Tumor_Sample_Barcode),"all")

####################
#no eqtl analysis

#significant_eqtls = read.csv(header = T,row.names = 1,file = "significant_eqtls.csv")
#significant_eqtls_anova = read.csv(header = T,row.names = 1,file = "ANOVA_significant_eqtls.csv")
#eqtl_expression_plot = read.csv(file = "standadized_filtered_rna-seq_input_to_peer.csv",row.names = 1)  
#eqtl_mut_plot = read.csv(file = "matrixEqtl_muts.csv",row.names = 1)
#eqtl_mut_plot_anova  =read.csv(file = "matrixEqtl_muts_ANOVA.csv",row.names = 1)
#eqtl_cis_linear = read.csv(file = "cis_eqtls.csv",row.names = 1)
#eqtl_cis_anova=read.csv(file = "ANOVA_cis_eqtls.csv",row.names = 1) 
#load(file = "peer_out.RData")
tcga_sig_estimate_num_nsig_added=readRDS("updated_sigs.rds")
#######################
# NO TMB ANALYSIS already carried out

#tmb_mat=read.csv("turk_tcga_tmb_per_VAF_cutoff.csv",row.names = 1)

################


onco_plotter2= function(x,fdr_cutoff){
  x$`-log10(FDR)`=-log10(x$fdr)
  x=pivot_wider(data = x,id_cols = c(Hugo_Symbol,cohort),names_from = cohort,names_prefix = "-log10(FDR)",values_from = `-log10(FDR)`,values_fill = 0)
  return(ggplot(data = x)+geom_point(aes(x =.data[[paste0("-log10(FDR)",unique(x$cohort)[1])]],x =.data[[paste0("-log10(FDR)",unique(x$cohort)[2])]] )))
}

onco_plotter= function(x,fdr_cutoff){
  x$`-log10(FDR)`=-log10(x$fdr)
  
  x$repel = ""
  x$repel[x$fdr<fdr_cutoff] = paste0(x$Hugo_Symbol[x$fdr<fdr_cutoff],"[",x$clusters[x$fdr<fdr_cutoff],"]")
  x$`if significant`="No"
  x$`if significant`[x$fdr<fdr_cutoff] = "Yes"
  return(ggplot(data = x)+aes(x=muts_in_clusters,y=`-log10(FDR)`)+geom_point(aes(col=`if significant`))+geom_text_repel(aes(label=repel),max.overlaps = 100)+facet_wrap(facets = vars(cohort)))
}

subset_MAF = function(maf_obj,
                      subtype_method = "Subtype_PAM50",
                      subtype_subset = NULL,
                      sample_subset = NULL,
                      vaf_subset = 0,
                      gene_subset = NULL,
                      tier_subset = NULL){
  if(!is.null(vaf_subset)){
    vaf_subset = paste0("VAF >= ",vaf_subset)
  }
  if(!is.null(subtype_subset)){
    if(subtype_method== "Subtype_PAM50"){
      subtype_subset = paste0("Subtype_PAM50 %in% c(",paste0(subtype_subset,collapse = ","),")")
    }
    else if(subtype_method== "Subtype_IHC"){
      subtype_subset = paste0("Subtype_IHC %in% c(",paste0(subtype_subset,collapse = ","),")")
    }
    
    else if(subtype_method== "Subtype_SCMOD1"){
      subtype_subset = paste0("Subtype_SCMOD1 %in% c(",paste0(subtype_subset,collapse = ","),")")
    }
  }
  subsetted_MAF = subsetMaf(maf = maf_obj,
                            tsb = sample_subset,
                            genes = gene_subset,
                            query = vaf_subset,
                            clinQuery = subtype_subset
                            )
  ### now subset by tiers
  if(!is.null(tier_subset)){
  tier_subset = paste0("pcgr_TIER %in% c(",paste0(tier_subset,collapse = ","),")")
  subsetted_MAF = subsetMaf(maf = subsetted_MAF,
                            query = tier_subset)
  }
  return(subsetted_MAF)
}

do_cophnetic = function(maf,vaf_cutoff,num_try,bootstrap){
  print(vaf_cutoff)
  maf_new = subsetMaf(maf = maf,query = paste0("VAF>=",vaf_cutoff))
  num_pat= length(unique(maf_new@data$Tumor_Sample_Barcode))
  if(num_pat<=3){
    return(NA)
  }else if(num_try>=num_pat){
    num_try = num_pat-1
  }
  tnm = tryCatch({
    trinucleotideMatrix(maf = maf_new, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
  },
    error = function(e){
      return(trinucleotideMatrix(maf = maf_new, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19",prefix = "chr"))
    },
  warning=function(w){
    return(trinucleotideMatrix(maf = maf_new, add = TRUE, ref_genome = "BSgenome.Hsapiens.UCSC.hg19",prefix = "chr"))
  }
  )
  estimated_num = estimateSignatures(mat = tnm, nTry = num_try,nrun = bootstrap)
  return(list(num_pat=num_pat,tnm=tnm,estimated_num=estimated_num))
}


###### Oncoplot fuction
onco_plot = function(  maf_obj, 
                       if_display_titv = F,
                       display_subtype = NULL,
                       sort_by_subtype=F,
                       top_genes = 2,
                       gene_fontSize = 0.8,
                       SampleNamefontSize = 1,
                       titleFontSize = 1.5,
                       legendFontSize=1.2,
                       annotationFontSize = 1.2
                       )
{
  
  oncoplot(
    maf = maf_obj,
    top = top_genes,
    draw_titv = if_display_titv,
    clinicalFeatures = display_subtype,
    sortByAnnotation = sort_by_subtype,
    fontSize = gene_fontSize,
    SampleNamefontSize = SampleNamefontSize,
    titleFontSize = titleFontSize,
    legendFontSize=legendFontSize,
    annotationFontSize =annotationFontSize
  )
}

##### Plotting summary VAF distribution


plot_vaf = function(maf_obj,
                    order_by_median = F,
                    top_genes_to_display = 20
                    )
  {
  plotVaf(maf = maf_obj,orderByMedian = order_by_median,top = top_genes_to_display)
}


do_kegg_maf=function(eqs,genes){
  input_eqs = eqs[eqs$Hugo_Symbol%in%genes,]
  return(compareCluster(entrez_expression~Hugo_Symbol+direction,fun="enrichKEGG",data = input_eqs))
}


#### Some modified maftools pathway analysis codes modified

get_pw_summary2 = function(maf, pathways = NULL){
  #maf = turkish_ihc_added
  if(is.null(pathways)){
    pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
    pathdb = data.table::fread(input = pathdb)
  }else{
    pathdb = data.table::copy(x = pathways)
    colnames(pathdb) = c("Gene", "Pathway")
    data.table::setDT(x = pathdb)
  }
  pathdb_size = pathdb[,.N,Pathway]
  pathdb = split(pathdb, as.factor(pathdb$Pathway))
  
  
  altered_pws = lapply(pathdb, function(pw){
    x = suppressMessages(try(genesToBarcodes(maf = maf, genes = pw$Gene)))
    if(class(x) != "NULL"){
      pw_genes = names(genesToBarcodes(maf = maf, genes = pw$Gene, justNames = TRUE, verbose = FALSE))
    }else{
      pw_genes = NULL
    }
    pw_genes
  })
  
  mut_load = lapply(altered_pws, function(x){
    if(is.null(x)){
      nsamps =  0
    }else{
      nsamps = length(unique(as.character(unlist(
        genesToBarcodes(
          maf = maf,
          genes = x,
          justNames = TRUE
        )
      ))))
    }
    nsamps
  })
  
  altered_pws = as.data.frame(t(data.frame(lapply(altered_pws, length))))
  data.table::setDT(x = altered_pws, keep.rownames = TRUE)
  colnames(altered_pws) = c("Pathway", "n_affected_genes")
  altered_pws$Pathway = names(pathdb)
  
  altered_pws = merge(pathdb_size, altered_pws, all.x = TRUE,)
  # all(pathdb_size$Pathway==altered_pws$Pathway)
  altered_pws[, fraction_affected := n_affected_genes/N]
  Mutated_samples = unlist(mut_load)
  Mutated_samples = Mutated_samples[order(names(Mutated_samples))]
  altered_pws = as.data.frame(altered_pws)
  altered_pws = altered_pws[order(altered_pws$Pathway),]
  all(names(Mutated_samples)==altered_pws$Pathway)
  altered_pws$Mutated_samples = Mutated_samples
  altered_pws = as.data.table(altered_pws)
  nsamps = as.numeric(maf@summary[ID == "Samples", summary])
  altered_pws[,Fraction_mutated_samples := Mutated_samples/nsamps]
  altered_pws = altered_pws[order(n_affected_genes, fraction_affected, decreasing = FALSE)]
  altered_pws
}

mafCompare_pathway= function (m1, m2, m1Name = NULL, m2Name = NULL, minMut = 5, useCNV = TRUE, 
                              pathways = FALSE, custom_pw = NULL, pseudoCount = FALSE) 
{
  m1.gs <- getGeneSummary(x = m1)
  m2.gs <- getGeneSummary(x = m2)
  if (is.null(m1Name)) {
    m1Name = "M1"
  }
  if (is.null(m2Name)) {
    m2Name = "M2"
  }
  m1.sampleSize = as.numeric(m1@summary[3, summary])
  m2.sampleSize = as.numeric(m2@summary[3, summary])
  sampleSummary = data.table::data.table(Cohort = c(m1Name, 
                                                    m2Name), SampleSize = c(m1.sampleSize, m2.sampleSize))
  if (pathways) {
    m1_pw = get_pw_summary2(maf = m1, pathways = custom_pw)[, 
                                                            .(Pathway, Mutated_samples)]
    m2_pw = get_pw_summary2(maf = m2, pathways = custom_pw)[, 
                                                            .(Pathway, Mutated_samples)]
    m.gs.meged = merge(m1_pw, m2_pw, by = "Pathway", all = TRUE)
    m.gs.meged = as.data.frame(m.gs.meged)
  }
  else {
    if (useCNV) {
      m1.genes = as.character(m1.gs[AlteredSamples >= minMut, 
                                    Hugo_Symbol])
      m2.genes = as.character(m2.gs[AlteredSamples >= minMut, 
                                    Hugo_Symbol])
      uniqueGenes = unique(c(m1.genes, m2.genes))
    }
    else {
      m1.genes = as.character(m1.gs[MutatedSamples >= minMut, 
                                    Hugo_Symbol])
      m2.genes = as.character(m2.gs[MutatedSamples >= minMut, 
                                    Hugo_Symbol])
      uniqueGenes = unique(c(m1.genes, m2.genes))
    }
    m1.gs.comGenes = m1.gs[Hugo_Symbol %in% uniqueGenes]
    m2.gs.comGenes = m2.gs[Hugo_Symbol %in% uniqueGenes]
    if (useCNV) {
      m.gs.meged = merge(m1.gs.comGenes[, .(Hugo_Symbol, 
                                            AlteredSamples)], m2.gs.comGenes[, .(Hugo_Symbol, 
                                                                                 AlteredSamples)], by = "Hugo_Symbol", all = TRUE)
    }
    else {
      m.gs.meged = merge(m1.gs.comGenes[, .(Hugo_Symbol, 
                                            MutatedSamples)], m2.gs.comGenes[, .(Hugo_Symbol, 
                                                                                 MutatedSamples)], by = "Hugo_Symbol", all = TRUE)
    }
    m.gs.meged[is.na(m.gs.meged)] = 0
    m.gs.meged = as.data.frame(m.gs.meged)
  }
  if (nrow(m.gs.meged) == 0) {
    stop("No genes pass the minMut threshold. Try decreasing the value..")
  }
  fisherTable = lapply(seq_len(nrow(m.gs.meged)), function(i) {
    gene = m.gs.meged[i, 1]
    m1Mut = m.gs.meged[i, 2]
    m2Mut = m.gs.meged[i, 3]
    ft_mat = matrix(c(m1Mut, m1.sampleSize - m1Mut, m2Mut, 
                      m2.sampleSize - m2Mut), byrow = TRUE, nrow = 2)
    if (length(which(x = ft_mat == 0)) > 0) {
      if (pseudoCount) {
        ft_mat = ft_mat + 1
      }
    }
    xf = fisher.test(ft_mat, conf.int = TRUE, conf.level = 0.95)
    pval = xf$p.value
    or = xf$estimate
    ci.up = xf$conf.int[2]
    ci.low = xf$conf.int[1]
    tdat = data.table::data.table(Hugo_Symbol = gene, m1Mut, 
                                  m2Mut, pval = pval, or = or, ci.up = ci.up, ci.low = ci.low)
    tdat
  })
  fisherTable = data.table::rbindlist(l = fisherTable, use.names = TRUE, 
                                      fill = TRUE)
  fisherTable = fisherTable[order(pval)]
  fisherTable[, `:=`(adjPval, p.adjust(p = pval, method = "fdr"))]
  colnames(fisherTable)[2:3] = c(m1Name, m2Name)
  return(list(results = fisherTable, SampleSummary = sampleSummary))
}

OncogenicPathways2 = function (maf, pathways = NULL, fontSize = 1, panelWidths = c(2,4, 4)) 
{
  #maf = turkish_ihc_added
  altered_pws = get_pw_summary2(maf = maf, pathways = pathways)
  altered_pws = altered_pws[!n_affected_genes %in% 0]
  if (nrow(altered_pws) == 0) {
    stop("None of the provided pathways are altered!")
  }
  nsamps = as.numeric(maf@summary[ID == "Samples", summary])
  lo = graphics::layout(mat = matrix(c(1, 2, 3), ncol = 3, nrow = 1), 
              widths = panelWidths)
  par(mar = c(2, 2, 2, 0))
  plot(NA, xlim = c(0, 1), ylim = c(0, nrow(altered_pws)), 
       axes = FALSE)
  text(x = 1, y = seq(0.5, nrow(altered_pws), by = 1), labels = altered_pws$Pathway, 
       adj = 1, xpd = TRUE, cex = fontSize)
  title(main = "Pathway", adj = 0)
  par(mar = c(2, 0, 2, 1), xpd = TRUE)
  plot(NA, xlim = c(0, 1.2), ylim = c(0, nrow(altered_pws)), 
       axes = FALSE)
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), 
       xright = 1, ytop = seq(0.2, nrow(altered_pws), by = 1) + 
         0.7, col = "#ecf0f1", border = "white")
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), 
       xright = altered_pws$fraction_affected, ytop = seq(0.2, 
                                                          nrow(altered_pws), by = 1) + 0.7, col = "#c0392b", 
       border = "white")
  text(x = 1.05, y = seq(0.5, nrow(altered_pws), by = 1), labels = paste0(altered_pws$n_affected_genes, 
                                                                          "/", altered_pws$N), adj = 0, cex = fontSize)
  axis(side = 1, at = seq(0, 1, 0.25), line = -1, cex.axis = fontSize)
  title(main = "Fraction of pathway affected", adj = 0)
  par(mar = c(2, 0, 2, 1), xpd = TRUE)
  plot(NA, xlim = c(0, 1.2), ylim = c(0, nrow(altered_pws)), 
       axes = FALSE)
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), 
       xright = 1, ytop = seq(0.2, nrow(altered_pws), by = 1) + 
         0.7, col = "#ecf0f1", border = "white")
  rect(xleft = 0, ybottom = seq(0.1, nrow(altered_pws), by = 1), 
       xright = altered_pws$Fraction_mutated_samples, ytop = seq(0.2, 
                                                                 nrow(altered_pws), by = 1) + 0.7, col = "#c0392b", 
       border = "white")
  text(x = 1.05, y = seq(0.5, nrow(altered_pws), by = 1), labels = paste0(altered_pws$Mutated_samples, 
                                                                          "/", nsamps), adj = 0, cex = fontSize)
  axis(side = 1, at = seq(0, 1, 0.25), line = -1, cex.axis = fontSize)
  title(main = "Fraction of samples affected")
  altered_pws
}


##### deconstruct the clonal structure
##### same func with infer hetero but with two clusters by default
# 
# inferHeterogeneity_2clust = function(maf, tsb = NULL, top = 5, vafCol = NULL, segFile = NULL, 
#                                       ignChr = NULL, minVaf = 0, maxVaf = 1, useSyn = FALSE, dirichlet = FALSE,num_clust = 2) 
# {
#   if (is.null(tsb)) {
#     tsb = as.character(getSampleSummary(x = maf)[1:top, Tumor_Sample_Barcode])
#   }
#   dat.tsb = subsetMaf(maf = maf, tsb = tsb, includeSyn = useSyn, 
#                       mafObj = FALSE)
#   if (nrow(dat.tsb) == 0) {
#     stop(paste(tsb, "not found in MAF"))
#   }
#   onlyContigs = as.character(seq(1:22))
#   dirichlet = FALSE
#   if (!"t_vaf" %in% colnames(dat.tsb)) {
#     if (is.null(vafCol)) {
#       if (all(c("t_ref_count", "t_alt_count") %in% colnames(dat.tsb))) {
#         message("t_vaf field is missing, but found t_ref_count & t_alt_count columns. Estimating vaf..")
#         dat.tsb[, `:=`(t_vaf, t_alt_count/(t_ref_count + 
#                                              t_alt_count))]
#       }
#       else {
#         print(colnames(dat.tsb))
#         stop("t_vaf field is missing. Use vafCol to manually specify vaf column name.")
#       }
#     }
#     else {
#       colnames(dat.tsb)[which(colnames(dat.tsb) == vafCol)] = "t_vaf"
#     }
#   }
#   if (!is.null(segFile)) {
#     seg.dat = readSegs(segFile)
#     seg.dat = seg.dat[Chromosome %in% onlyContigs]
#     seg.dat = seg.dat[order(as.numeric(Chromosome))]
#     setkey(x = seg.dat, Chromosome, Start_Position, End_Position)
#     seg.tsbs = unique(seg.dat[, Sample])
#     if (sum(!tsb %in% seg.tsbs) > 0) {
#       message("CN data for following samples not found. Ignoring them ..")
#       print(tsb[!tsb %in% seg.tsbs])
#       seg.tsbs = tsb[tsb %in% seg.tsbs]
#     }
#     else {
#       seg.tsbs = tsb
#     }
#     if (length(seg.tsbs) > 0) {
#       seg.dat = seg.dat[Sample %in% seg.tsbs]
#     }
#     else {
#       stop("None of the provided samples have copy number data.")
#     }
#   }
#   clust.dat = c()
#   dat.tsb = dat.tsb[, .(Hugo_Symbol, Chromosome, Start_Position, 
#                         End_Position, Tumor_Sample_Barcode, t_vaf)]
#   dat.tsb$Chromosome = as.character(dat.tsb$Chromosome)
#   dat.tsb$t_vaf = as.numeric(as.character(dat.tsb$t_vaf))
#   dat.tsb = dat.tsb[order(Chromosome)]
#   if (max(dat.tsb$t_vaf, na.rm = TRUE) > 1) {
#     dat.tsb$t_vaf = dat.tsb$t_vaf/100
#   }
#   if (!is.null(ignChr)) {
#     dat.tsb = dat.tsb[!Chromosome %in% ignChr]
#   }
#   dat.tsb = dat.tsb[t_vaf > minVaf][t_vaf < maxVaf]
#   dat.tsb$Chromosome = gsub(pattern = "chr", replacement = "", 
#                             x = dat.tsb$Chromosome, fixed = TRUE)
#   dat.tsb = suppressWarnings(dat.tsb[order(as.numeric(Chromosome))])
#   for (i in 1:length(tsb)) {
#     message("Processing ", tsb[i], "..")
#     tsb.dat = dat.tsb[Tumor_Sample_Barcode %in% tsb[i]]
#     tsb.dat = tsb.dat[!is.na(tsb.dat$t_vaf), ]
#     tempCheck = 0
#     if (!is.null(segFile)) {
#       if (tsb[i] %in% seg.tsbs) {
#         seg = seg.dat[Sample %in% tsb[i]]
#         seg.res = filterCopyNumber(seg, tsb.dat, tempCheck, 
#                                    tsb[i])
#         tsb.dat = seg.res[[1]]
#         tsb.dat.cn.vars = seg.res[[2]]
#         tempCheck = seg.res[[3]]
#       }
#     }
#     if (nrow(tsb.dat) < 3) {
#       message("Too few mutations for clustering. Skipping..")
#     }
#     else {
#       if (dirichlet) {
#         tsb.dat = dirichletClusters(tsb.dat)
#       }
#       else {
#         tsb.cluster = mclust::densityMclust(tsb.dat[, 
#                                                     t_vaf], G = num_clust, verbose = FALSE)
#         tsb.dat$cluster = as.character(tsb.cluster$classification)
#         abs.med.dev = abs(tsb.dat[, t_vaf] - median(tsb.dat[, 
#                                                             t_vaf]))
#         pat.mad = median(abs.med.dev) * 100
#         pat.math = pat.mad * 1.4826/median(tsb.dat[, 
#                                                    t_vaf])
#         tsb.dat$MATH = pat.math
#         tsb.dat$MedianAbsoluteDeviation = pat.mad
#       }
#       tsb.dat = refineClusters(clusters = tsb.dat)
#       if (tempCheck == 1) {
#         tsb.dat = rbind(tsb.dat, tsb.dat.cn.vars, fill = TRUE)
#       }
#       clust.dat = rbind(clust.dat, tsb.dat, fill = TRUE)
#     }
#   }
#   if (is.null(clust.dat)) {
#     message("No result, this is basically caused by no copy number neutral variants,\n you may re-run this without copy number data.")
#   }
#   else {
#     clust.dat.mean = clust.dat[, mean(t_vaf), by = .(Tumor_Sample_Barcode, 
#                                                      cluster)]
#     colnames(clust.dat.mean)[ncol(clust.dat.mean)] = "meanVaf"
#     return(list(clusterData = clust.dat, clusterMeans = clust.dat.mean))
#   }
# }
# 
# 
# ### lets just make a test runc with infer.
# 
# refineClusters = function(clusters){
#   totclusts = unique(clusters[,cluster])
#   refinedClusters = c()
#   
#   for(i in 1:length(totclusts)){
#     clust.dat = clusters[cluster == totclusts[i]]
#     if(nrow(clust.dat) > 1){
#       clust.dat = clusters[cluster == totclusts[i]]
#       
#       #Based on boxplot stats
#       outl = boxplot.stats(x = clust.dat[,t_vaf])$out
#       clust.dat$cluster = ifelse(test = clust.dat$t_vaf %in% outl, yes = 'outlier', no = clust.dat$cluster)
#       
#       #Based on zscore
#       #clust.dat$z = scale(clust.dat[,t_vaf], center = TRUE, scale = TRUE)
#       #clust.dat$cluster = ifelse(test = abs(clust.dat$z) >= 2.5 , no = clust.dat$cluster, yes = 'outlier')
#       #clust.dat$cluster = ifelse(test = abs(clust.dat$z) >= 2.5 , no = clust.dat$cluster, yes = 'outlier')
#       #clust.dat[,z := NULL]
#     }else{
#       clust.dat$cluster = 'outlier'
#     }
#     refinedClusters = rbind(refinedClusters, clust.dat)
#   }
#   
#   temp.dat = refinedClusters[cluster %in% c('CN_altered', 'outlier')]
#   refinedClusters = refinedClusters[!cluster %in% c('CN_altered', 'outlier')]
#   
#   clust.lvls = levels(factor(refinedClusters$cluster))
#   refinedClusters$cluster = as.character(factor(refinedClusters$cluster, levels = clust.lvls, labels = 1:length(clust.lvls)))
#   refinedClusters = rbind(refinedClusters, temp.dat)
#   
#   return(refinedClusters)
# }
# 
# test_infer2 = inferHeterogeneity_2clust(maf = turkish_subtype_added_pcgr_added_MAF,vafCol = "VAF",useSyn = F )
# 
# plotClusters(clusters = test_infer2)
# plot(1,1)
# 
# test_infer2 = inferHeterogeneity_2clust(maf = turkish_subtype_added_pcgr_added_MAF,vafCol = "VAF",useSyn = T,tsb = turkish_subtype_added_pcgr_added_MAF@clinical.data$Tumor_Sample_Barcode)
# clusters =  test_infer2$clusterData
# 
# max_vaf_cluster1=clusters%>%group_by(Tumor_Sample_Barcode,cluster)%>%summarise(max(t_vaf))%>%filter(cluster==1)
# 
# library(ggplot2)
# plot_clust2_syn=ggplot(max_vaf_cluster1)+geom_boxplot(aes(x = cluster,y = `max(t_vaf)`))
# 
# test_infer2_non_syn = inferHeterogeneity_2clust(maf = turkish_subtype_added_pcgr_added_MAF,vafCol = "VAF",useSyn = F,tsb = turkish_subtype_added_pcgr_added_MAF@clinical.data$Tumor_Sample_Barcode)
# 
# clusters_ns =  test_infer2_non_syn$clusterData
# 
# max_vaf_cluster1_ns=clusters_ns%>%group_by(Tumor_Sample_Barcode,cluster)%>%summarise(max(t_vaf))%>%filter(cluster==1)
# 
# library(ggplot2)
# plot_clust2_ns = ggplot(max_vaf_cluster1_ns)+geom_boxplot(aes(x = cluster,y = `max(t_vaf)`))
# 
# test_infer2 = inferHeterogeneity_2clust(maf = turkish_subtype_added_pcgr_added_MAF,vafCol = "VAF",useSyn = F )
# 
# plotClusters(clusters = test_infer2)
# plot(1,1)
# 
# test_infer2 = inferHeterogeneity_2clust(maf = turkish_subtype_added_pcgr_added_MAF,vafCol = "VAF",useSyn = T,tsb = turkish_subtype_added_pcgr_added_MAF@clinical.data$Tumor_Sample_Barcode,num_clust = 3)
# clusters =  test_infer2$clusterData
# 
# max_vaf_cluster1=clusters%>%group_by(Tumor_Sample_Barcode,cluster)%>%summarise(max(t_vaf))%>%filter(cluster==1)
# 
# library(ggplot2)
# plot_clust3_syn=ggplot(max_vaf_cluster1)+geom_boxplot(aes(x = cluster,y = `max(t_vaf)`))
# plot_clust3_syn
# test_infer2_non_syn = inferHeterogeneity_2clust(maf = turkish_subtype_added_pcgr_added_MAF,vafCol = "VAF",useSyn = F,tsb = turkish_subtype_added_pcgr_added_MAF@clinical.data$Tumor_Sample_Barcode,num_clust = 3)
# 
# clusters_ns =  test_infer2_non_syn$clusterData
# 
# max_vaf_cluster1_ns=clusters_ns%>%group_by(Tumor_Sample_Barcode,cluster)%>%summarise(max(t_vaf))%>%filter(cluster==1)
# 
# library(ggplot2)
# plot_clust3_ns = ggplot(max_vaf_cluster1_ns)+geom_boxplot(aes(x = cluster,y = `max(t_vaf)`))
# plot_clust3_ns
# 
# combine_vaf = c(turkish_subtype_added_pcgr_added_MAF@data$VAF,tcga_subtype_added@data$VAF)
# combine_vaf = data.frame(VAF = combine_vaf,center = rep(c("TBCP","TCGA"),c(length(turkish_subtype_added_pcgr_added_MAF@data$VAF),length(tcga_subtype_added@data$VAF))))
# 
# library(ggsignif) 
# ggplot(data = combine_vaf,aes(x=VAF,colour = center))+stat_ecdf(geom = "step") #+geom_signif(aes(x = center,y= VAF),test = ks.test)
# 
# combine_vaf2 = combine_vaf[combine_vaf$VAF>=0.07,]
# 
# ggplot(data = combine_vaf2,aes(x=VAF,colour = center))+stat_ecdf(geom = "step") #+geom_signif(aes(x = center,y= VAF),test = ks.test)
# 
# ggplot(data = combine_vaf2,aes(x=VAF,colour = center))+geom_density()
