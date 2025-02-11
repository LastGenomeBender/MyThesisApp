
#WRITTEN BY FARID AHADLI (farid.ahadli@outlook.com)

#SERVER PART OF THE APP!!!


server <- function(input, output, session) {
  
  ############################################################################################################################################################     
  ##### server side of the Mutation analysis
  
  ## server of  "Mutation and Clinical Data" 
  ## turkish_subtype_added_pcgr_added_MAF and tcga_subtype_added are fixed.
  
  
  uploaded_maf = reactive({
    req(input$upload_files)
    maf =  read.delim(file = input$maf_file$datapath,header = T,sep = "\t",comment.char = '#',na.strings = c("NA","N/A"))
    maf$VAF = maf$t_alt_count/( maf$t_alt_count+ maf$t_ref_count)
    mtch = match(maf$Hugo_Symbol,hugo.to.ms$MutSig_Synonym)
    maf$Hugo_Symbol[which(mtch>0)] = hugo.to.ms$Hugo_Symbol[na.exclude(mtch)]
    maf = read.maf(maf =maf,
                   clinicalData = read.delim(file = input$clinical_file$datapath,header = T,sep = "\t",comment.char = '#'),
                   isTCGA = F)
    return(maf)
  }
  )
  tcga_in_use = reactive({
    req(input$upload_files)
    tcga_mut = rbind.data.frame(tcga_subtype_added@data,tcga_subtype_added@maf.silent)
    tcga_clin = tcga_subtype_added@clinical.data
    mtch = match(tcga_mut$Hugo_Symbol,hugo.to.ms$MutSig_Synonym)
    tcga_mut$Hugo_Symbol[which(mtch>0)] = hugo.to.ms$Hugo_Symbol[na.exclude(mtch)]
    tcga_subtype_added = read.maf(maf =tcga_mut,
                                  clinicalData = tcga_clin,
                                  isTCGA = T)
    if(as.character(input$if_targeted_seq)=="FALSE"){
      
      return(tcga_subtype_added)
    }
    else{
      capture = read.delim(file = input$capture_file$datapath,header = T,sep = "\t",comment.char = '#')  
      genes = as.character(capture[,1])
      mtch = match(genes,hugo.to.ms$MutSig_Synonym)
      genes[which(mtch>0)] = hugo.to.ms$Hugo_Symbol[na.exclude(mtch)]
      subsetted_tcga = subsetMaf(maf = tcga_subtype_added,genes = genes)
      
      return(subsetted_tcga)
    }
  })  
  
  ######################################################  
  
  maf_for_filter = reactiveValues(maf = 2)
  
  observeEvent(input$upload_files,{
    
    req(input$upload_files)
    req(uploaded_maf())
    maf_for_filter$maf = uploaded_maf()
  })
  
  observeEvent(input$maf_reset_filt,{
    req(input$upload_files)
    req(uploaded_maf())
    maf_for_filter$maf = uploaded_maf()
  })
  
  observeEvent(input$maf_clin_or_mut,{
    req(input$upload_files)
    req(uploaded_maf())
    maf =maf_for_filter$maf
    if(input$maf_clin_or_mut == "Mutation"){
      updateSelectizeInput(session = session,
                           inputId = "maf_column",
                           label = "Select the column you would like to apply filter",
                           choices = colnames(maf@data),
                           selected =colnames(maf@data)[1],server = T )
      
    }else{
      updateSelectizeInput(session = session,
                           inputId = "maf_column",
                           label = "Select the column you would like to apply filter",
                           choices = colnames(maf@clinical.data),
                           selected =colnames(maf@clinical.data)[1],server = T)
    }
  })
  
  output$filter_controls = renderUI({
    req(input$upload_files)
    req(uploaded_maf())
    maf =maf_for_filter$maf
    if(input$maf_clin_or_mut == "Mutation"){
      data = as.data.frame(maf@data)
      data = data %>% mutate_if(all.is.numeric, as.numeric)
    }else{
      data = as.data.frame(maf@clinical.data)
      data = data %>%mutate_if(all.is.numeric, as.numeric)
      print(class(data[,input$maf_column]))
    }
    if((class(data[1,input$maf_column])=="numeric")|(class(data[1,input$maf_column])=="integer") ){
      data_rng = range(data[,input$maf_column])
      out = paste0("your data range is ",data_rng )
      tagList(
        h3(cat(out)),
        numericInput(inputId = "maf_number",label = "Enter a number between the range", value = 0,min = data_rng[1],max = data_rng[2]),
        selectInput(inputId = "maf_num_compare",label = "Select the operation",choices =c(`Equals`="==",
                                                                                          `Greater_than`=">",
                                                                                          `Smaller_than`="<",
                                                                                          `Greater_than_or_equal_to`=">=",
                                                                                          `Smaller_than_or_equal_to`="<=",
                                                                                          `Not_equals` = "!="
        ),multiple = F )
      )
    }else{
      data_rng = unique(data[,input$maf_column])
      selectizeInput(inputId = "maf_select_char",label = "Select the values to include",choices =data_rng,multiple = T )
    }
  })
  
  observeEvent(input$maf_apply_filt,{
    req(input$upload_files)
    req(uploaded_maf())
    maf = maf_for_filter$maf
    if(input$maf_clin_or_mut == "Mutation"){
      data = as.data.frame(maf@data)
    }else{
      data = as.data.frame(maf@clinical.data)
    }
    if(class(data[,input$maf_column])=="numeric"){
      filter = paste0(input$maf_column,input$maf_num_compare,input$maf_number)
      if(input$maf_clin_or_mut == "Mutation"){
        maf_for_filter$maf = subsetMaf(maf = maf,query = filter)
      }else{
        maf_for_filter$maf = subsetMaf(maf = maf,clinQuery = filter)
      }
    }else{
      deprsed = deparse(input$maf_select_char)
      filter = paste0(input$maf_column,"%in%",deprsed)
      if(input$maf_clin_or_mut == "Mutation"){
        maf_for_filter$maf = subsetMaf(maf = maf,query = filter)
      }else{
        maf_for_filter$maf = subsetMaf(maf = maf,clinQuery = filter)
      }
    }
    
  })
  
  
  output$maf_coding_muts = renderDT({
    data = maf_for_filter$maf@data
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
    
  output$download_nonsyn <- downloadHandler(
    filename = function() {
    "nonsyn.csv"
    },
    content = function(file) {
      write.csv(maf_for_filter$maf@data, file)
    }
  )
    
  output$maf_silent_muts = renderDT({
    data = maf_for_filter$maf@maf.silent
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
    
    
  output$download_syn <- downloadHandler(
    filename = function() {
      "syn.csv"
    },
    content = function(file) {
      write.csv(maf_for_filter$maf@maf.silent, file)
    }
  )  
  
  
  output$maf_clinic = renderDT({
    data = maf_for_filter$maf@clinical.data
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
   
  output$download_clinical <- downloadHandler(
    filename = function() {
      "clinical.csv"
    },
    content = function(file) {
      write.csv(maf_for_filter$maf@clinical.data, file)
    }
  )  
  
    
  ########################Server oncoplot#############  
  observe({
    req(input$upload_files)
    req(uploaded_maf())
    maf = maf_for_filter$maf
    updateSelectizeInput(session = session,inputId = "oncoplot_if_annotate_subtype",label = "Choose Clinical Data Column to Add",choices = colnames(maf@clinical.data),server = T)
  })
  observeEvent(input$oncoplot_draw,{
    req(input$upload_files)
    req(uploaded_maf())
    maf = maf_for_filter$maf
    updateSelectizeInput(session = session,inputId = "oncoplot_if_annotate_subtype",label = "Choose Clinical Data Column to Add",choices = colnames(maf@clinical.data),server = T)
    output$oncoplot = renderPlot({
      add_subtype = isolate(input$oncoplot_if_annotate_subtype)
      # if(nchar(add_subtype)==5){
      #   add_subtype = NULL
      # }
      onco_plot(maf_obj = maf,
                if_display_titv = isolate(input$oncoplot_if_annotate_titv),
                display_subtype = add_subtype,
                sort_by_subtype = isolate(input$oncoplot_group_by_subtype),
                top_genes = isolate(input$oncoplot_num_gene),
                gene_fontSize = isolate(input$oncoplot_gene_font),
                SampleNamefontSize = isolate(input$oncoplot_sample_font),
                titleFontSize = isolate(input$oncoplot_title_font),
                legendFontSize = isolate(input$oncoplot_legend_font),
                annotationFontSize = isolate(input$oncoplot_legend_font)
                
                  
      )
    },width = 1000,height = 900)
  })
  output$download_oncoplot = downloadHandler(
    filename = function() {
      "oncoplot.png"
    },
    content = function(file) {
     png(filename = file,width = 9,height = 9,units = "in",res = 300)
      maf = maf_for_filter$maf
      onco_plot(maf_obj = maf,
                if_display_titv = isolate(input$oncoplot_if_annotate_titv),
                display_subtype = isolate(input$oncoplot_if_annotate_subtype),
                sort_by_subtype = isolate(input$oncoplot_group_by_subtype),
                top_genes = isolate(input$oncoplot_num_gene),
                gene_fontSize = isolate(input$oncoplot_gene_font),
                SampleNamefontSize = isolate(input$oncoplot_sample_font),
                titleFontSize = isolate(input$oncoplot_title_font),
                legendFontSize = isolate(input$oncoplot_legend_font),
                annotationFontSize = isolate(input$oncoplot_legend_font)
                
                
      )
      dev.off()
    }
  )
  ########################### VAF Filter Cutoff Server #####################
  tmb_mat = eventReactive(input$tmb_calculate,{
    req(uploaded_maf())
    tmb = lapply(seq(0,0.85,0.05),function(x){
      #x= 0.1
      tmb_mat_turk = data.frame(tmb = getSampleSummary(subsetMaf(maf_for_filter$maf,query = paste0("VAF>=",x)))$total,cohort = input$cohort_name,cutoff = x)
      tmb_mat_tcga = data.frame(tmb = getSampleSummary(subsetMaf(tcga_in_use(),query = paste0("VAF>=",x)))$total,cohort = "TCGA",cutoff = x)
      mat = rbind.data.frame(tmb_mat_turk,tmb_mat_tcga)
    })
    tmb_mat=do.call(rbind.data.frame,tmb)
    return(tmb_mat)
  })
  
  tmb_plot = reactive({
    req(tmb_mat())
    input$tmb_test
    input$tmb_axis_label_font
    input$tmb_axis_tick_font
    tmb_mat = tmb_mat()
    tmb_mat$tmb = log2(tmb_mat$tmb)
    tmb_mat_list = split(x = tmb_mat,tmb_mat$cutoff)
    plot_tmb_list = lapply(tmb_mat_list,function(x){
        ggplot(data = x,aes(x=cohort,y=tmb))+geom_boxplot(outlier.shape = NA)+stat_pvalue_manual(data = compare_means(formula = tmb~cohort,data = x,method =input$tmb_test,paired = F),y.position = 14)+geom_quasirandom(alpha=0.5)+ylab(label = "log2(tmb)")+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$tmb_axis_label_font),axis.text =element_text(face = "bold",size = input$tmb_axis_tick_font),text = element_text(face = "bold"))+ggtitle(paste0("VAF>=",x$cutoff[1]))+ylim(c(0,17))
      
    })
    figure_4 = cowplot::plot_grid(plotlist = plot_tmb_list,labels = "AUTO",ncol = 3)
    return(figure_4)
  })
  output$tmb_plot = renderPlot({
    tmb_plot()
  },height = 2000)
  output$download_tmb_plot <- downloadHandler(
    filename = function() {
     "tmb_plot.png"
    },
    content = function(file) {
      png(filename = file,width = 20,height = 20,units = "in",res = 300)
      plot(tmb_plot())
      dev.off()
    }
  )
  output$tmb_table = renderDT({
    data = tmb_mat()
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_tmb_table <- downloadHandler(
    filename = function() {
     "tmb_mat.csv"
    },
    content = function(file) {
      write.csv(tmb_mat(), file)
    }
  )  
    
  
  ########################### Signatures#####################
  
  tur_selected_vaf_sigs = eventReactive(input$sig_run,{
    do_cophnetic(maf = maf_for_filter$maf,vaf_cutoff = input$sig_vaf_cutoff,num_try = input$sig_nmax,bootstrap = input$sig_bootstrap)
    #turkish_sig_estimate_num_nsig_added[[as.character(input$sig_vaf_cutoff)]]
  })
  
  # turkish_sig_estimate_num_nsig_added[[x]][["sigs"]]=extractSignatures(mat = turkish_sig_estimate_num_nsig_added[[x]]$tnm,
  #                                                                      n = turkish_sig_estimate_num_nsig_added[[x]]$nsig)
  # turkish_sig_estimate_num_nsig_added[[x]][["cosm"]]=compareSignatures(nmfRes = turkish_sig_estimate_num_nsig_added[[x]][["sigs"]], sig_db = "SBS")
  
 extracted_sigs = eventReactive(input$extract_sigs,{
   extractSignatures(mat = tur_selected_vaf_sigs()$tnm,
                     n = input$num_signatures)
 })
  
 # cossim_sig = eventReactive(input$extract_sigs,{
 #   req(extracted_sigs())
 #   compareSignatures(nmfRes = extracted_sigs(), sig_db = "SBS")
 # })
  output$sig_coph = renderPlot({
    req(tur_selected_vaf_sigs())
    plotCophenetic(res =tur_selected_vaf_sigs()$estimated_num)
  },width = 900,height = 1000)
  
  output$download_sig_coph <- downloadHandler(
    filename = function() {
      "cophnetic_plot.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plotCophenetic(res =tur_selected_vaf_sigs()$estimated_num)
      dev.off()
    }
  )
  
  output$sig_sigs = renderPlot({
    req(extracted_sigs())
    plotSignatures(nmfRes = extracted_sigs(), sig_db = "SBS",contributions = F,show_barcodes = F,font_size = input$sig_sig_plot_font,title_size = input$sig_sig_plot_font_title,yaxisLim = NA)
  },width = 1000,height = 900)
  
  output$download_sig_sigs = downloadHandler(
    filename = function() {
      "signature_plot.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plotSignatures(nmfRes = extracted_sigs(), sig_db = "SBS",contributions = F,show_barcodes = F,font_size = input$sig_sig_plot_font,title_size = input$sig_sig_plot_font_title,yaxisLim = NA)
      dev.off()
    }
  )
  
  
  output$sig_conts = renderPlot({
    req(extracted_sigs())
    maftools::plotSignatures(nmfRes = extracted_sigs(), sig_db = "SBS",contributions = T,show_barcodes = T,font_size = input$sig_cont_plot_font)
  },width = 1000,height = 900)
  
  output$download_sig_conts = downloadHandler(
    filename = function() {
      "signature_contribution_plot.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      maftools::plotSignatures(nmfRes = extracted_sigs(), sig_db = "SBS",contributions = T,show_barcodes = T,font_size = input$sig_cont_plot_font)
      dev.off()
    }
  )
  
  
  output$sig_conts_table = renderDT({
    data = t(extracted_sigs()$contributions)
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
    
   output$download_sig_conts_table = downloadHandler(
     filename = function() {
      "signature_contriburion_table.csv"
     },
     content = function(file) {
       write.csv( t(extracted_sigs()$contributions), file)
     }
   )
    
    
    
  
  ####
  
  sig_correlations = reactive({
    cont_try=extracted_sigs()$contributions
    tmb_try = maf_for_filter$maf@data%>%group_by(Tumor_Sample_Barcode)%>%summarise(tmb=n())
    tmb_try = tmb_try[order(tmb_try$Tumor_Sample_Barcode),]
    print(nrow(tmb_try))
    #clin=maf_for_filter$maf@clinical.data
    #clin = clin[match(tmb_try$Tumor_Sample_Barcode,clin$Tumor_Sample_Barcode),]
    
    cont_try = as.data.frame(t(cont_try))
    cont_try = cont_try[order(rownames(cont_try)),]
    cont_try = cont_try[tmb_try$Tumor_Sample_Barcode,]
    tmb_try = dplyr::bind_cols(tmb_try,cont_try)
    tmb_try_long = pivot_longer(data = tmb_try,cols = 3:ncol(tmb_try),names_to = "signature",values_to = "contribution")
    #tmb_try_long$A = as.numeric(tmb_try_long$Age)
    return(tmb_try_long)
  })
  
  
  
  clin_updated = reactive({
    req(uploaded_maf())
    clin = getClinicalData(maf_for_filter$maf)
    clin = clin %>%mutate(across(where(~all.is.numeric(.x,extras = c("NA",NA,"n/a","N/A","NOT KNOWN",tolower("NOT KNOWN"),"NOT DETECTED"))),as.numeric))
   
    if_cont = unlist(lapply(clin,FUN = all.is.numeric,extras = c("NA",NA,"n/a","N/A","NOT KNOWN",tolower("NOT KNOWN"),"NOT DETECTED")))
    clin = as.data.frame(clin)
    tsb_col = which(colnames(clin)=="Tumor_Sample_Barcode")
    return(list(continous=clin[,c(tsb_col,which(if_cont))],categoric=clin[,c(tsb_col,which(!if_cont))]))
  })
  
  observe({
    req(clin_updated())
    updateSelectizeInput(session = session,inputId = "clin_cont_col",choices = colnames(clin_updated()$continous),server = T)
  })
  observe({
    req(clin_updated())
    updateSelectizeInput(session = session,inputId = "clin_char_col",choices = colnames(clin_updated()$categoric),server = T)
  })
  
  sig_correlations_continous = eventReactive(input$clin_cont_run,{
    cont_try=extracted_sigs()$contributions
    tmb_try = clin_updated()$continous[,c("Tumor_Sample_Barcode",input$clin_cont_col)]
    tmb_try = tmb_try[order(tmb_try$Tumor_Sample_Barcode),]
    #clin=maf_for_filter$maf@clinical.data
    #clin = clin[match(tmb_try$Tumor_Sample_Barcode,clin$Tumor_Sample_Barcode),]
    cont_try = as.data.frame(t(cont_try))
    cont_try = cont_try[order(rownames(cont_try)),]
    tmb_try = cbind.data.frame(tmb_try,cont_try)
    tmb_try_long = pivot_longer(data = tmb_try,cols = 3:ncol(tmb_try),names_to = "signature",values_to = "contribution")
    #tmb_try_long$A = as.numeric(tmb_try_long$Age)
    return(tmb_try_long)
  })
  sig_correlations_char = eventReactive(input$clin_char_run,{
    cont_try=extracted_sigs()$contributions
    tmb_try = clin_updated()$categoric[,c("Tumor_Sample_Barcode",input$clin_char_col)]
    tmb_try = tmb_try[order(tmb_try$Tumor_Sample_Barcode),]
    #clin=maf_for_filter$maf@clinical.data
    #clin = clin[match(tmb_try$Tumor_Sample_Barcode,clin$Tumor_Sample_Barcode),]
    cont_try = as.data.frame(t(cont_try))
    cont_try = cont_try[order(rownames(cont_try)),]
    tmb_try = cbind.data.frame(tmb_try,cont_try)
    tmb_try_long = pivot_longer(data = tmb_try,cols = 3:ncol(tmb_try),names_to = "signature",values_to = "contribution")
    #tmb_try_long$A = as.numeric(tmb_try_long$Age)
    return(tmb_try_long)
  })
  
  output$sig_tmb_cor = renderPlot({
    ggplot(data = sig_correlations(),aes(x = contribution,y=tmb))+geom_point()+stat_poly_line()+stat_poly_eq(aes(label=paste(after_stat(eq.label),after_stat(p.value.label), sep = "*\", \"*")))+facet_wrap(facets = vars(signature))+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_cor_label_font),axis.text =element_text(face = "bold",size = input$clin_cor_tick_font),text = element_text(face = "bold"))+ylim(c(0,input$clin_cor_ylim))
  },width = 900,height = 1000)
  
  output$download_sig_tmb_cor =  downloadHandler(
    filename = function() {
      "sig_tmb_cor.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = sig_correlations(),aes(x = contribution,y=tmb))+geom_point()+stat_poly_line()+stat_poly_eq(aes(label=paste(after_stat(eq.label),after_stat(p.value.label), sep = "*\", \"*")))+facet_wrap(facets = vars(signature))+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_cor_label_font),axis.text =element_text(face = "bold",size = input$clin_cor_tick_font),text = element_text(face = "bold"))+ylim(c(0,input$clin_cor_ylim)))
      dev.off()
    }
  )
 
  output$sig_cont_cor = renderPlot({
    ggplot(data = sig_correlations_continous(),aes_string(x = "contribution",y=input$clin_cont_col))+geom_point()+stat_poly_line()+stat_poly_eq(aes(label=paste(after_stat(eq.label),after_stat(p.value.label), sep = "*\", \"*")))+facet_wrap(facets = vars(signature))+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_cont_label_font),axis.text =element_text(face = "bold",size = input$clin_cont_tick_font),text = element_text(face = "bold"))+ylim(c(0,input$clin_cont_ylim))
  },width = 900,height = 1000)
  
  output$download_sig_cont_cor =  downloadHandler(
    filename = function() {
      "sig_cont_cor.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = sig_correlations_continous(),aes_string(x = "contribution",y=input$clin_cont_col))+geom_point()+stat_poly_line()+stat_poly_eq(aes(label=paste(after_stat(eq.label),after_stat(p.value.label), sep = "*\", \"*")))+facet_wrap(facets = vars(signature))+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_cont_label_font),axis.text =element_text(face = "bold",size = input$clin_cont_tick_font),text = element_text(face = "bold"))+ylim(c(0,input$clin_cont_ylim)))
      dev.off()
    }
  )
  
  output$sig_anova = renderPlot({
    ggplot(data = sig_correlations(),aes(x = signature,y=contribution))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "anova",label.y = 1.1)+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_anova_label_font),axis.text =element_text(face = "bold",size = input$clin_anova_tick_font),text = element_text(face = "bold"))
  },width = 900,height = 1000)
  
  output$download_sig_anova =  downloadHandler(
    filename = function() {
      "sig_anova.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = sig_correlations(),aes(x = signature,y=contribution))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "anova",label.y = 1.1)+theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_anova_label_font),axis.text =element_text(face = "bold",size = input$clin_anova_tick_font),text = element_text(face = "bold")))
      dev.off()
    }
  )
  
  
  
  output$sig_char_cor = renderPlot({
    ggplot(data = sig_correlations_char(),aes(x = signature,y=contribution))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "anova",label.y = 1.1)+facet_wrap(as.formula(paste("~", input$clin_char_col)))+
      theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_char_label_font),axis.text =element_text(face = "bold",size = input$clin_char_tick_font),text = element_text(face = "bold"),axis.text.x = element_text(angle = 90))
  },width = 900,height = 1000)
  
  output$download_sig_char_cor =  downloadHandler(
    filename = function() {
      "sig_char_cor.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = sig_correlations_char(),aes(x = signature,y=contribution))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "anova",label.y = 1.1)+facet_wrap(as.formula(paste("~", input$clin_char_col)))+
        theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_char_label_font),axis.text =element_text(face = "bold",size = input$clin_char_tick_font),text = element_text(face = "bold"),axis.text.x = element_text(angle = 90)))
      dev.off()
    }
  )
  
  ##########
  
  observe({
    req(clin_updated())
    updateSelectizeInput(session = session,inputId = "clin_pair_col",choices = colnames(clin_updated()$categoric),server = T)
  })
  
  observe({
    req(clin_updated())
    req(input$clin_pair_col)
    updateSelectizeInput(session = session,inputId = "clin_pair_f",choices = unique(clin_updated()$categoric[,input$clin_pair_col]),server = T)
    
  })
  observe({
    req(input$clin_pair_col)
    req(clin_updated())
    updateSelectizeInput(session = session,inputId = "clin_pair_s",choices = unique(clin_updated()$categoric[,input$clin_pair_col]),server = T)
    
  })
  
  
  sig_pair_char = eventReactive(input$clin_pair_run,{
    cont_try=extracted_sigs()$contributions
    tmb_try = clin_updated()$categoric[,c("Tumor_Sample_Barcode",input$clin_pair_col)]
    bool = tmb_try[,input$clin_pair_col] %in% c(input$clin_pair_f,input$clin_pair_s)
    
    tmb_try = tmb_try[bool,]
    print(nrow(tmb_try))
    tmb_try = tmb_try[order(tmb_try$Tumor_Sample_Barcode),]
    #clin=maf_for_filter$maf@clinical.data
    #clin = clin[match(tmb_try$Tumor_Sample_Barcode,clin$Tumor_Sample_Barcode),]
    cont_try = as.data.frame(t(cont_try))
    cont_try = cont_try[tmb_try$Tumor_Sample_Barcode,]
    print(nrow(cont_try))
    cont_try = cont_try[order(rownames(cont_try)),]
    tmb_try = cbind.data.frame(tmb_try,cont_try)
    
    tmb_try_long = pivot_longer(data = tmb_try,cols = 3:ncol(tmb_try),names_to = "signature",values_to = "contribution")
    print(colnames(tmb_try_long))
    #tmb_try_long$A = as.numeric(tmb_try_long$Age)
    return(tmb_try_long)
  })
  
  output$sig_pair_cor = renderPlot({
    ggplot(data = sig_pair_char(),aes_string(x = input$clin_pair_col,y="contribution"))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "t.test",label.y = 1.1)+facet_wrap(as.formula(paste("~", "signature")))+
      theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_pair_label_font),axis.text =element_text(face = "bold",size = input$clin_pair_tick_font),text = element_text(face = "bold"),axis.text.x = element_text(angle = 90))
  },width = 900,height = 1000)
  
  
  output$download_sig_pair_cor =  downloadHandler(
    filename = function() {
      "sig_pair_cor.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = sig_pair_char(),aes_string(x = input$clin_pair_col,y="contribution"))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "t.test",label.y = 1.1)+facet_wrap(as.formula(paste("~", "signature")))+
        theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_pair_label_font),axis.text =element_text(face = "bold",size = input$clin_pair_tick_font),text = element_text(face = "bold"),axis.text.x = element_text(angle = 90)))
      dev.off()
    }
  )
  
  
  ################ within signature analysis
  
  observe({
    req(clin_updated())
    updateSelectizeInput(session = session,inputId = "clin_pair_col_sig",choices = colnames(clin_updated()$categoric),server = T)
  })
  
  observe({
    req(extracted_sigs())
    cont_try=extracted_sigs()$contributions
    cont_try = as.data.frame(t(cont_try))
    updateSelectizeInput(session = session,inputId = "clin_pair_f_sig",choices = colnames(cont_try),server = T)
    
  })
  observe({
    req(extracted_sigs())
    cont_try=extracted_sigs()$contributions
    cont_try = as.data.frame(t(cont_try))
    updateSelectizeInput(session = session,inputId = "clin_pair_s_sig",choices = colnames(cont_try),server = T)
    
  })
  
  
  sig_pair_char_sig = eventReactive(input$clin_pair_run_sig,{
    cont_try=extracted_sigs()$contributions
    tmb_try = clin_updated()$categoric[,c("Tumor_Sample_Barcode",input$clin_pair_col_sig)]
    #bool = tmb_try[,input$clin_pair_col] %in% c(input$clin_pair_f,input$clin_pair_s)
    
    #tmb_try = tmb_try[bool,]
    #print(nrow(tmb_try))
    tmb_try = tmb_try[order(tmb_try$Tumor_Sample_Barcode),]
    #clin=maf_for_filter$maf@clinical.data
    #clin = clin[match(tmb_try$Tumor_Sample_Barcode,clin$Tumor_Sample_Barcode),]
    cont_try = as.data.frame(t(cont_try))
    cont_try = cont_try[tmb_try$Tumor_Sample_Barcode,]
    print(nrow(cont_try))
    cont_try = cont_try[order(rownames(cont_try)),]
    tmb_try = cbind.data.frame(tmb_try,cont_try)
    
    tmb_try_long = pivot_longer(data = tmb_try,cols = 3:ncol(tmb_try),names_to = "signature",values_to = "contribution")
    tmb_try_long = tmb_try_long[tmb_try_long$signature %in% c(input$clin_pair_f_sig,input$clin_pair_s_sig),]
    #print(colnames(tmb_try_long))
    #tmb_try_long$A = as.numeric(tmb_try_long$Age)
    return(tmb_try_long)
  })
  
  output$sig_pair_cor_sig = renderPlot({
    ggplot(data = sig_pair_char_sig(),aes_string(x = "signature",y="contribution"))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "t.test",label.y = 1.1)+facet_wrap(as.formula(paste("~", input$clin_pair_col_sig)))+
      theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_pair_label_font_sig),axis.text =element_text(face = "bold",size = input$clin_pair_tick_font_sig),text = element_text(face = "bold"),axis.text.x = element_text(angle = 90))
  },width = 900,height = 1000)
  
  
  output$download_sig_pair_cor_sig =  downloadHandler(
    filename = function() {
      "sig_pair_cor_sig.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = sig_pair_char_sig(),aes_string(x = "signature",y="contribution"))+geom_boxplot()+geom_quasirandom()+stat_compare_means(method = "t.test",label.y = 1.1)+facet_wrap(as.formula(paste("~", input$clin_pair_col_sig)))+
             theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$clin_pair_label_font),axis.text =element_text(face = "bold",size = input$clin_pair_tick_font),text = element_text(face = "bold"),axis.text.x = element_text(angle = 90)))
      dev.off()
    }
  )
  
  
  
  
#########clinical enrichment###############
  observe({
    req(uploaded_maf())
    maf = maf_for_filter$maf
    updateSelectizeInput(session = session,
                         inputId = "clin_enrich_column",
                         label = "Select the column you would like to apply filter",
                         choices = colnames(maf@clinical.data),
                         selected =colnames(maf@clinical.data)[1],server = T)
  })
  
 enrich_result = eventReactive(input$clin_enrich_run,{
   clinicalEnrichment(maf = maf_for_filter$maf,clinicalFeature = input$clin_enrich_column)
 })
    
 output$clin_enrich_pairwise = renderDT({
   data = enrich_result()$pairwise_comparision
   datatable(data = data,class = 'nowrap',filter = "top", options = list(  
     scrollX = TRUE,
     pageLength=5))
 })
 output$clin_enrich_groupwise  = renderDT({
   data = enrich_result()$groupwise_comparision
   datatable(data = data,class = 'nowrap',filter = "top", options = list(  
     scrollX = TRUE,
     pageLength=5))
 })
 output$enrich_plot = renderPlot({
   plotEnrichmentResults(enrich_res = enrich_result(),pVal = input$clin_enrich_pval,ORthr = input$clin_enrich_OR,annoFontSize = input$clin_enrich_annotation_font,geneFontSize = input$clin_enrich_gene_font,legendFontSize = input$clin_enrich_legend_font)
 },width = 900,height = 1000)
 ###
 #download buttons
    output$download_pairwise_enrich = downloadHandler(
      filename = function() {
        "pairwise_enrichmemt_results.csv"
      },
      content = function(file) {
        write.csv(enrich_result()$pairwise_comparision, file)
      }
    )
    output$download_groupwise_enrich = downloadHandler(
      filename = function() {
        "group_enrichmemt_results.csv"
      },
      content = function(file) {
        write.csv(enrich_result()$groupwise_comparision, file)
      }
    )
    output$download_groupwise_enrich = downloadHandler(
      filename = function() {
        "enrichmemt_results.png"
      },
      content = function(file) {
        png(filename = file,width = 9,height = 9,units = "in",res = 300)
        plotEnrichmentResults(enrich_res = enrich_result(),pVal = input$clin_enrich_pval,ORthr = input$clin_enrich_OR,annoFontSize = input$clin_enrich_annotation_font,geneFontSize = input$clin_enrich_gene_font,legendFontSize = input$clin_enrich_legend_font)
        dev.off()
      }
    )
######################################################### Comparisons##################
  tcga_maf_compare = eventReactive(input$vaf_cutoff_run,{
    subsetMaf(maf = tcga_in_use(),query = paste0("VAF>=",input$vaf_cutoff_compare))
  })
  turk_maf_compare = reactive({
    subsetMaf(maf = maf_for_filter$maf,query = paste0("VAF>=",input$vaf_cutoff_compare))
    
  })
  tcga_vs_turk =  eventReactive(input$vaf_cutoff_run,{
    req(turk_maf_compare)
    req(tcga_maf_compare)
    maf1= turk_maf_compare()
    maf2= tcga_maf_compare()
    mafCompare(m1 = maf1, m2 = maf2, m1Name = input$cohort_name, m2Name = 'TCGA', minMut = 5)
  })
  
  output$diffmut_table = renderDT({
    data =  tcga_vs_turk()$results
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  output$download_diffmut_table = downloadHandler(
    filename = function() {
      "diffrentially_mutated_genes.csv"
    },
    content = function(file) {
      write.csv(tcga_vs_turk()$results, file)
    }
  )
  
    
  
  
  output$diffmut_forest = renderPlot({
    req(tcga_vs_turk)
    comp_maf = tcga_vs_turk()
    forestPlot(mafCompareRes = comp_maf,fdr =input$vaf_cutoff_fdr ,geneFontSize = input$vaf_cutoff_gene_size,titleSize = input$vaf_cutoff_title_size )
  },height = 900)
  
  
  output$download_diffmut_forest = downloadHandler(
    filename = function() {
      "forest_plot.png"
    },
    content = function(file) {
      comp_maf = tcga_vs_turk()
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      
      forestPlot(mafCompareRes = comp_maf,fdr =input$vaf_cutoff_fdr ,geneFontSize = input$vaf_cutoff_gene_size,titleSize = input$vaf_cutoff_title_size )
      dev.off()
    }
  )
  
  
  ###################################################
  #### oncoCLUST
  tcga_maf_clust= eventReactive(input$oncoclust_run,{
    subsetMaf(maf = tcga_in_use(),query = paste0("VAF>=",input$oncoclust_vaf_cutoff))
  })
  turk_maf_clust = eventReactive(input$oncoclust_run,{
    subsetMaf(maf = maf_for_filter$maf,query = paste0("VAF>=",input$oncoclust_vaf_cutoff))
  })
  
  ## perform the oncoclust analysis
  
  tcga_oncoclust = eventReactive(input$oncoclust_run,{
    req(tcga_maf_clust)
    oncodrive(maf = tcga_maf_clust(),AACol ="HGVSp_Short",bgEstimate = F)
  })
  
  turk_oncoclust = eventReactive(input$oncoclust_run,{
    req(turk_maf_clust)
    oncodrive(maf = turk_maf_clust(),AACol ="HGVSp_Short",bgEstimate = F)
  })
  
  output$oncoclust_tcga = renderDT({
    data = as.data.frame(tcga_oncoclust())
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_oncoclust_tcga =  downloadHandler(
    filename = function() {
      "TCGA_oncoclust_res.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(tcga_oncoclust()), file)
    }
  )
  
    
  
  output$oncoclust_turk = renderDT({
    data = as.data.frame(turk_oncoclust())
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
    
  output$download_oncoclust_turk =  downloadHandler(
    filename = function() {
      "Uploaded_oncoclust_res.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(turk_oncoclust()), file)
    }
  )
  
  output$oncoclust_plot = renderPlotly({
    req(tcga_oncoclust())
    req(turk_oncoclust())
    test_df = dplyr::bind_rows(tcga_oncoclust(), turk_oncoclust())
    test_df$cohort = rep(c("TCGA",input$cohort_name),c(nrow(tcga_oncoclust()),nrow(turk_oncoclust())))
    x =test_df
    x$`-log10(FDR)`=-log10(x$fdr)
    dta=pivot_wider(data = x,id_cols = c(Hugo_Symbol,cohort),names_from = cohort,names_prefix = "-log10(FDR) ",values_from = `-log10(FDR)`,values_fill = 0)
    ggplotly(ggplot(data = dta,aes(x=.data[[colnames(dta[2])]],y=.data[[colnames(dta[3])]],label = Hugo_Symbol))+geom_point(position = position_jitter(width = .1, height = .1),alpha=.1) + geom_vline(xintercept = -log10(input$oncoclust_fdr_cutoff),color="red")+geom_hline(yintercept = -log10(input$oncoclust_fdr_cutoff),color="red")+
      theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$onco_axis_label_font),axis.text =element_text(face = "bold",size = input$onco_axis_tick_font),text = element_text(face = "bold")))
  })
      observe({
        req(uploaded_maf())
        updateSelectizeInput(session = session,"oncoclust_lolipop_genes","Select Genes",choices = c(unique(maf_for_filter$maf@data$Hugo_Symbol),unique(tcga_in_use()@data$Hugo_Symbol)))
      })
      
  
  output$download_oncoclust_plot = downloadHandler(
    filename = function() {
      "onclust.png"
    },
    content = function(file) {
      test_df = dplyr::bind_rows(tcga_oncoclust(), turk_oncoclust())
      test_df$cohort = rep(c("TCGA",input$cohort_name),c(nrow(tcga_oncoclust()),nrow(turk_oncoclust())))
      x =test_df
      x$`-log10(FDR)`=-log10(x$fdr)
      dta=pivot_wider(data = x,id_cols = c(Hugo_Symbol,cohort),names_from = cohort,names_prefix = "-log10(FDR) ",values_from = `-log10(FDR)`,values_fill = 0)
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(ggplot(data = dta)+geom_point(aes(x=.data[[colnames(dta[2])]],y=.data[[colnames(dta[3])]]),position = position_jitter(width = .1, height = .1),alpha=.1) + geom_vline(xintercept = -log10(input$oncoclust_fdr_cutoff),color="red")+geom_hline(yintercept = -log10(input$oncoclust_fdr_cutoff),color="red")+
             theme_classic(base_size = 12)+theme(axis.title =  element_text(face = "bold",size = input$onco_axis_label_font),axis.text =element_text(face = "bold",size = input$onco_axis_tick_font),text = element_text(face = "bold")))
      dev.off()
    }
  )
  
  
  output$oncoclust_lollipop2 = renderPlot({
    #oncoclust_lolipop_genes
    lollipopPlot2(m1 = turk_maf_clust(), m2 = tcga_maf_clust(), gene = input$oncoclust_lolipop_genes, AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = input$cohort_name, m2_name = "TCGA",m1_label ="all",m2_label ="all", showDomainLabel = F,labPosAngle =input$oncoclust_lollipop_angle,labPosSize = input$oncoclust_lollipop_label_text_size,axisTextSize = rep(input$oncoclust_lollipop_axis_text_size,2),pointSize = input$oncoclust_lollipop_head_size,legendTxtSize =input$oncoclust_lollipop_legend_text_size )
  })
  
  output$download_oncoclust_lollipop2 = downloadHandler(
    filename = function() {
      "lollipop.png"
    },
    content = function(file) {
      
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      lollipopPlot2(m1 = turk_maf_clust(), m2 = tcga_maf_clust(), gene = input$oncoclust_lolipop_genes, AACol1 = "HGVSp_Short", AACol2 = "HGVSp_Short", m1_name = input$cohort_name, m2_name = "TCGA",m1_label ="all",m2_label ="all", showDomainLabel = F,labPosAngle =input$oncoclust_lollipop_angle,labPosSize = input$oncoclust_lollipop_label_text_size,axisTextSize = rep(input$oncoclust_lollipop_axis_text_size,2),pointSize = input$oncoclust_lollipop_head_size,legendTxtSize =input$oncoclust_lollipop_legend_text_size )
      dev.off()
    }
  )
  #############################################
  ##### somtatic interactions
  tcga_maf_somint= eventReactive(input$somint_run,{
    subsetMaf(maf = tcga_in_use(),query = paste0("VAF>=",input$somint_vaf_cutoff))
  })
  turk_maf_somint = eventReactive(input$somint_run,{
    subsetMaf(maf = maf_for_filter$maf,query = paste0("VAF>=",input$somint_vaf_cutoff))
  })
  
  
  
  
  output$somint_turk_plot = renderPlot({
    somaticInteractions(maf = turk_maf_somint(),top = 25,  fontSize = input$somint_fontsize,
                        sigSymbolsSize = input$somint_asteriks_size,
                        sigSymbolsFontSize = input$somint_legend_fontsize,
                        showSum = input$somint_if_show_sum)
  },width = 900,height = 1000)
  
  output$download_somint_turk_plot = downloadHandler(
    filename = function() {
      "uploaded_somint.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      somaticInteractions(maf = turk_maf_somint(),top = 25,fontSize = input$somint_fontsize,
                          sigSymbolsSize = input$somint_asteriks_size,
                          sigSymbolsFontSize = input$somint_legend_fontsize,
                          showSum = input$somint_if_show_sum)
      dev.off()
    }
  )
  
  output$somint_tcga_plot = renderPlot({
    somaticInteractions(maf = tcga_maf_somint(),top = 25,fontSize = input$somint_fontsize,
                        sigSymbolsSize = input$somint_asteriks_size,
                        sigSymbolsFontSize = input$somint_legend_fontsize,
                        showSum = input$somint_if_show_sum)
  },width = 900,height = 1000)
  
  output$download_somint_tcga_plot = downloadHandler(
    filename = function() {
      "tcga_somint.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      somaticInteractions(maf = tcga_maf_somint(),top = 25,fontSize = input$somint_fontsize,
                          sigSymbolsSize = input$somint_asteriks_size,
                          sigSymbolsFontSize = input$somint_legend_fontsize,
                          showSum = input$somint_if_show_sum)
      dev.off()
    }
  )
  
  output$somint_turk_table = renderDT({
    data = as.data.frame(somaticInteractions(maf = turk_maf_somint(),top = 25))
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
    
   output$download_somint_turk_table = downloadHandler(
     filename = function() {
      "uploaded_somatic_interaction_table.csv"
     },
     content = function(file) {
       write.csv(as.data.frame(somaticInteractions(maf = turk_maf_somint(),top = 25)), file)
     }
   )
    
    
  output$somint_tcga_table = renderDT({
    data = as.data.frame(somaticInteractions(maf = tcga_maf_somint(),top = 25))
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  output$download_somint_tcga_table  = downloadHandler(
    filename = function() {
      "tcga_somatic_interaction_table.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(somaticInteractions(maf = tcga_maf_somint(),top = 25)), file)
    }
  )
    

  tcga_somint = reactive({
    data = as.data.frame(somaticInteractions(maf = tcga_maf_somint(),top = 25))
  })
  turk_somint = reactive({
    data = as.data.frame(somaticInteractions(maf = turk_maf_somint(),top = 25))
    
  })
  comparison_somint = reactive({
    data_turk = turk_somint()[turk_somint()$pValue<0.05,]
    data_tcga = tcga_somint()[tcga_somint()$pValue<0.05,]
    colnames(data_turk) = paste0("Uploaded_",colnames(data_turk))
    colnames(data_tcga) = paste0("TCGA_",colnames(data_tcga))
    data_tcga = data_tcga[order(data_tcga$TCGA_pair),]
    data_turk = data_turk[order(data_turk$Uploaded_pair),]
    unique_in_turk = setdiff(data_turk$Uploaded_pair,data_tcga$TCGA_pair)
    unique_in_tcga = setdiff(data_tcga$TCGA_pair,data_turk$Uploaded_pair)
    common_overall = intersect(data_turk$Uploaded_pair,data_tcga$TCGA_pair)
    common_in_same= intersect(paste0(data_turk$Uploaded_pair,"|",data_turk$Uploaded_Event),paste0(data_tcga$TCGA_pair,"|",data_tcga$TCGA_Event))
    common_in_same_only_pairs = unlist(strsplit(x = common_in_same,"|",fixed=T))
    common_in_same_only_pairs = common_in_same_only_pairs[seq(1,length(common_in_same_only_pairs),2)]
    common_but_opposite = setdiff(common_overall,common_in_same_only_pairs)
    #### now the dataframe
    df_unique_in_turk = data_turk[data_turk$Uploaded_pair %in% unique_in_turk,]
    df_unique_in_tcga = data_tcga[data_tcga$TCGA_pair %in% unique_in_tcga,]
    df_common_in_both = cbind.data.frame(data_turk[data_turk$Uploaded_pair %in% common_in_same_only_pairs,],data_tcga[data_tcga$TCGA_pair %in% common_in_same_only_pairs,])
    df_opposite =cbind.data.frame(data_turk[data_turk$Uploaded_pair %in% common_but_opposite,],data_tcga[data_tcga$TCGA_pair %in% common_but_opposite,])
    return(list(df_unique_in_turk=df_unique_in_turk,
                df_unique_in_tcga=df_unique_in_tcga,
                df_common_in_both = df_common_in_both,
                df_opposite = df_opposite
    ))
  })
  
  output$somint_uniq_tcga = renderDT({
    data = comparison_somint()$df_unique_in_tcga
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_somint_uniq_tcga  = downloadHandler(
    filename = function() {
      "tcga_uniq_int.csv"
    },
    content = function(file) {
      write.csv(comparison_somint()$df_unique_in_tcga, file)
    }
  )
  
  output$somint_uniq_turk = renderDT({
    data = comparison_somint()$df_unique_in_turk
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_somint_uniq_turk  = downloadHandler(
    filename = function() {
      "turk_uniq_int.csv"
    },
    content = function(file) {
      write.csv(comparison_somint()$df_unique_in_turk, file)
    }
  )
  
  output$somint_both = renderDT({
    data = comparison_somint()$df_common_in_both
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_somint_both  = downloadHandler(
    filename = function() {
      "significant_in_both_int.csv"
    },
    content = function(file) {
      write.csv(comparison_somint()$df_common_in_both, file)
    }
  )
  
  output$somint_opopsite = renderDT({
    data = comparison_somint()$df_opposite
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_somint_opopsite  = downloadHandler(
    filename = function() {
      "opposite_directions_int.csv"
    },
    content = function(file) {
      write.csv(comparison_somint()$df_opposite, file)
    }
  )
  
  
  ### make the tables
  
  #######################################
  ########### signature comparison.
  
  comp_tur_selected_vaf_sigs = reactive({
    turkish_sig_estimate_num_nsig_added[[as.character(input$comp_sig_vaf_cutoff)]]
  })
  
  comp_tcga_selected_vaf_sigs = reactive({
    tcga_sig_estimate_num_nsig_added[[as.character(input$sig_vaf_cutoff)]]
  })
  
  output$comp_tcga_coph = renderPlot({
    plotCophenetic(res =comp_tcga_selected_vaf_sigs()$estimated_num,bestFit =  comp_tcga_selected_vaf_sigs()$nsig)
  },width = 900,height = 1000)
  output$download_comp_turk_coph = downloadHandler(
    filename = function() {
      "uploaded_cophnetic.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plotCophenetic(res =tur_selected_vaf_sigs()$estimated_num)
      dev.off()
    }
  )
  output$download_comp_tcga_coph = downloadHandler(
    filename = function() {
      "tcga_cophnetic.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plotCophenetic(res =comp_tcga_selected_vaf_sigs()$estimated_num,bestFit =  comp_tcga_selected_vaf_sigs()$nsig)
      dev.off()
    }
  )
  
  
  output$comp_turk_coph = renderPlot({
    plotCophenetic(res =tur_selected_vaf_sigs()$estimated_num)
  },width = 900,height = 1000)
  
  output$comp_turk_sigs = renderPlot({
    plotSignatures(nmfRes = extracted_sigs(), sig_db = "SBS",contributions = F,show_barcodes = F,font_size = 2,title_size = 1.5)
  },width = 900,height = 1000)
  
  output$download_comp_turk_sigs = downloadHandler(
    filename = function() {
      "uploaded_sigs.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plotSignatures(nmfRes = extracted_sigs(), sig_db = "SBS",contributions = F,show_barcodes = F,font_size = 2,title_size = 1.5)
      dev.off()
    }
  )
  
  
  output$comp_tcga_sigs = renderPlot({
    plotSignatures(nmfRes = comp_tcga_selected_vaf_sigs()$sigs, sig_db = "SBS",contributions = F,show_barcodes = F,font_size = 2,title_size = 1.5)
  },width = 900,height = 1000)
  
  output$download_comp_tcga_sigs = downloadHandler(
    filename = function() {
      "tcga_sigs.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plotSignatures(nmfRes = comp_tcga_selected_vaf_sigs()$sigs, sig_db = "SBS",contributions = F,show_barcodes = F,font_size = 2,title_size = 1.5)
      dev.off()
    }
  )
  
  compare_cosine_similarities = reactive({
    cos_sims=cos_sim_matrix(mut_matrix1 = extracted_sigs()$signatures,comp_tcga_selected_vaf_sigs()$sigs$signatures)
    rownames(cos_sims) = paste0(input$cohort_name,"",rownames(cos_sims))
    colnames(cos_sims) = paste0("TCGA_",colnames(cos_sims)) 
    return(cos_sims)
  })
  output$comp_cosine = renderPlot({
    plot_cosine_heatmap(compare_cosine_similarities(), 
                        cluster_rows = TRUE, cluster_cols = TRUE,plot_values = T)
  },width = 900,height = 1000)
  
  output$download_comp_cosine = downloadHandler(
    filename = function() {
      "cosine_similarities.png"
    },
    content = function(file) {
      cos_sims=cos_sim_matrix(mut_matrix1 = extracted_sigs()$signatures,comp_tcga_selected_vaf_sigs()$sigs$signatures)
      rownames(cos_sims) = paste0(input$cohort_name,"",rownames(cos_sims))
      colnames(cos_sims) = paste0("TCGA_",colnames(cos_sims))
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      plot(plot_cosine_heatmap(cos_sim_matrix = cos_sims,cluster_rows = TRUE, cluster_cols = TRUE,plot_values = T))
      dev.off()
    }
  )
  
  
  
  #######################################
  ######Pathway Comparisons
  
  turk_maf_pathway = reactive({
    subsetMaf(maf = maf_for_filter$maf,query = paste0("VAF>=",input$path_summary_vaf_cutoff))
  })
  
  tcga_maf_pathway= reactive({
    subsetMaf(maf = tcga_in_use(),query = paste0("VAF>=",input$path_summary_vaf_cutoff))
  })
  
  pathway_compare = reactive({
    mafCompare_pathway(m1 = turk_maf_pathway(), m2 = tcga_maf_pathway(), m1Name = input$cohort_name, m2Name = 'TCGA', minMut = 5,pathways = T)
  })
  
  ### sumary tab
  
  summary_path_chosen = reactive({
    if(input$path_summary_choose=="Uploaded Dataset"){
      return(turk_maf_pathway())
    }
    else{
      return(tcga_maf_pathway())
    }
  }) 
  
  output$path_summary_table = renderDT({
    data = as.data.frame(OncogenicPathways2(maf = summary_path_chosen()))
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
  
  output$download_path_summary_table = downloadHandler(
    filename = function() {
      "pathway_summary.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(OncogenicPathways2(maf = summary_path_chosen())), file)
    }
  )
  
  output$path_summary_plot  =  renderPlot({
    OncogenicPathways2(maf = summary_path_chosen())
  },width = 1000,height = 900)
  
  output$download_path_summary_plot = downloadHandler(
    filename = function() {
      "pathway_summary.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      OncogenicPathways2(maf = summary_path_chosen())
      
      dev.off()
    }
  )
  
  ##comparison tab
  
  #pathway_compare
  
  output$path_comparison_table = renderDT({
    #print(pathway_compare()$results)
    data = as.data.frame(pathway_compare()$results)
    datatable(data = data,class = 'nowrap',filter = "top", options = list(  
      scrollX = TRUE,
      pageLength=5))
  })
    
  output$download_path_comparison_table = downloadHandler(
    filename = function() {
      "fisher_pathways.csv"
    },
    content = function(file) {
      write.csv(as.data.frame(pathway_compare()$results), file)
    }
  )  
  
  
  output$path_comparison_plot = renderPlot({
    forestPlot(mafCompareRes = pathway_compare(),fdr =input$path_comparison_fdr ,geneFontSize = input$path_comparison_gene_size,titleSize = input$path_comparison_title_size )
  },width = 1000,height = 900)
  
  output$download_path_comparison_plot = downloadHandler(
    filename = function() {
      "pathway_comparison.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      forestPlot(mafCompareRes = pathway_compare(),fdr =input$path_comparison_fdr ,geneFontSize = input$path_comparison_gene_size,titleSize = input$path_comparison_title_size )
      
      dev.off()
    }
  )
  
  ## pathway oncoplots
  
  onco_path_chosen = reactive({
    if(input$path_display_choose_dataset=="Uploaded Dataset"){
      return(turk_maf_pathway())
    }
    else{
      return(tcga_maf_pathway())
    }
  }) 
  
  output$path_display_plot = renderPlot({
    PlotOncogenicPathways(maf = onco_path_chosen(),pathways = input$path_display_choose_pathway,removeNonMutated = T,showTumorSampleBarcodes = input$path_display_show_tsb,fontSize =input$path_display_font_size,SampleNamefontSize =  input$path_display_sample_font_size)
  },width = 1000,height = 900)
  
  output$download_path_display_plot = downloadHandler(
    filename = function() {
      "pathway.png"
    },
    content = function(file) {
      png(filename = file,width = 9,height = 9,units = "in",res = 300)
      PlotOncogenicPathways(maf = onco_path_chosen(),pathways = input$path_display_choose_pathway,removeNonMutated = T,showTumorSampleBarcodes = input$path_display_show_tsb,fontSize =input$path_display_font_size,SampleNamefontSize =  input$path_display_sample_font_size)
      
      dev.off()
    }
  )
  
}





