
##WRITTEN BY FARID AHADLI (farid.ahadli@outlook.com)

#UI PART OF THE APP.
options(shiny.maxRequestSize = 300*1024^2)
#load packages
library(data.table)
library(ggpmisc)
library(ggpubr)
library(MutationalPatterns)
library(ggfortify)
library(ggrepel)
library(NMF)
library(ggsignif)
library(ggbreak)
library(ggbeeswarm)
library("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(clusterProfiler)
library(shiny)
library(shinydashboard)
library(DT)
library(rmarkdown)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)
library(shinycssloaders)
library(InteractiveComplexHeatmap)
library(DESeq2)
library(edgeR)
library(gplots)
#library(genefu)
library(ggplot2)
library(dplyr)
library(plotly)
library(shinythemes)
require(graphics)
library(tidyverse)
library(RColorBrewer)
library(heatmaply)
library(maftools)
library(shiny)
library(shinyjs)
library(Hmisc)
library(shinyWidgets)
##for mutation_analysis
#required files
#source("ClinicalData.R")
#source("GeneLists.R")
#source("TCGA_PCA.R")
source("ahadli_mutation_function.R")

#XQuartz app is needed for heatmaply.

#set repositories
options(repos = BiocManager::repositories())

#gene_list <- rownames(numReads)

navbarPage("Ahadli",
           
           useShinyjs(),
           
           ######################################## 
           ## Ui of upload
           tabPanel("Data Upload",
                    sidebarPanel(width=3,
                                 fluidPage(
                                   
                                   fileInput(inputId = "maf_file",
                                             
                                             label = "Upload your MAF file as a .txt document",
                                             multiple = F,
                                             accept = ".txt"),
                                   fileInput(inputId = "clinical_file",
                                             label = "Upload your clinical data file as a .txt document",
                                             multiple = F,
                                             accept = ".txt"),
                                   checkboxInput(inputId = "if_targeted_seq",label = "Is your data generated with a targeted sequencing protocol?",value = F),
                                   
                                   conditionalPanel(condition = "input.if_targeted_seq == true",
                                                    fileInput(inputId = "capture_file",
                                                              label = "Upload the name of captrued genes",multiple = F,
                                                              accept = ".txt")
                                   ), ## MIGHT GIVE ERROR
                                   actionButton(inputId = "upload_files",
                                                label = "Click to upload",
                                                icon = icon("fa-fire")),
                                   br(),
                                   br(),
                                   dropMenu(
                                     actionButton("go0", "info",icon = icon("info")),
                                     "You can read more about MAF specification from this ",
                                     tags$a(href = "https://docs.gdc.cancer.gov/Encyclopedia/pages/Mutation_Annotation_Format_TCGAv2/", "link."),
                                     "If you have a vcf file, you can use the ",
                                     tags$a("vcf2maf",href="https://github.com/mskcc/vcf2maf"),
                                     "utility to convert it to the MAF format . Clinical file has to contain a column named ",
                                     strong('Tumor_Sample_Barcode'),
                                     ".",
                                     
                                     
                                     theme = "light-border",
                                     placement = "bottom",
                                     arrow = T
                                   ),
                                   br(),
                                   
                                 )
                    )
           ),
           ###################################################################################################################################################
           ######## Start of the MUTATION analysis
           tabPanel("Mutation and Clinical Data",
                    sidebarPanel(width = 3,
                                 fluidPage(
                                   selectizeInput(inputId = "maf_clin_or_mut","Select the attributes you wolud like to filter on",choices = c("Clinical","Mutation"),multiple=F, selected= "Mutation"), 
                                   selectizeInput(inputId = "maf_column","Select the column you would like to apply filter",choices =NULL,multiple=F),
                                   uiOutput(outputId = "filter_controls"),
                                   actionButton(inputId = "maf_apply_filt",label = "Apply the filters"),
                                   actionButton(inputId = "maf_reset_filt",label = "Reset the filters"),
                                   br(),
                                   br(),
                                   dropMenu(
                                     actionButton("go2", "info",icon = icon("info")),
                                     "The filters will be applied sequentially, e.g. if you first filter samples based on ER status, let's say ER+, by clicking",
                                     strong('Apply the filters'),
                                     "button, then any new filter you apply will work on the subsetted dataset. In essence you use an ",
                                     strong('AND'),
                                     "operator with each click. All the filters will be reset as soon as you click the ",
                                     strong("Reset the filters"),
                                     "Button",
                                     
                                     
                                     theme = "light-border",
                                     placement = "bottom",
                                     arrow = T
                                   ),
                                   
                                 )),
                    mainPanel(width = 9,
                              h3("Table of the coding mutations"),
                              dataTableOutput(outputId = "maf_coding_muts") %>% withSpinner(color="black"),
                              downloadButton(outputId ="download_nonsyn" ),
                              h3("Table of the silent mutations"),
                              dataTableOutput(outputId = "maf_silent_muts") %>% withSpinner(color="black"),
                              downloadButton(outputId = "download_syn"),
                              h3("Table of the associated clinical data"),
                              dataTableOutput(outputId = "maf_clinic") %>% withSpinner(color="black"),
                              downloadButton(outputId = "download_clinical")
                    )),
           ##################################################################################
           tabPanel("Oncoplot",
                    sidebarPanel(width = 3,
                                 fluidPage(
    
                                   selectizeInput(inputId = "oncoplot_if_annotate_subtype","Choose Clinical Data Column(s) to Add",choices = NULL,multiple=T),
                                   numericInput(inputId = "oncoplot_num_gene","Select number of genes to display",value = 2,min = 1,max = 100,step = 1),
                                   selectInput(inputId = "oncoplot_group_by_subtype","Group samples by subtype by the first selected Column",choices = list("Yes"=T,"No"=F),selected = "Yes",multiple = F),
                                   selectInput(inputId = "oncoplot_if_annotate_titv","Add transition,transversion annotation to oncoplot",choices = list("Yes"=T,"No"=F),selected = "Yes",multiple = F),
                                   numericInput(inputId = "oncoplot_gene_font", "Choose the gene fontsize",value = 0.8,min = 0.1,max = 20,step = 0.05),
                                   numericInput(inputId = "oncoplot_sample_font", "Choose the sample fontsize",value = 1,min = 0.1,max = 20,step = 0.05),
                                   numericInput(inputId = "oncoplot_title_font", "Choose the title fontsize",value = 1.5,min = 0.1,max = 20,step = 0.05),
                                   numericInput(inputId = "oncoplot_legend_font", "Choose the legend fontsize",value = 1.2,min = 0.1,max = 20,step = 0.05),
                                   numericInput(inputId = "oncoplot_annotation_font", "Choose the annotation fontsize",value = 1.2,min = 0.1,max = 20,step = 0.05),
                                   actionButton(inputId = "oncoplot_draw",label = "Click to Draw The Oncoplot")
                                 )),
                    mainPanel(width = 9,
                              downloadButton(outputId = "download_oncoplot"),
                              plotOutput(outputId = "oncoplot") %>% withSpinner(color="black")
                             
                    )
                    
           ),
           #################################################################
           tabPanel("Comparison of TMB",
                    sidebarPanel(width = 3,
                                 fluidPage(textInput(inputId = "cohort_name",label = "Enter the cohort name",value = "",placeholder ="metastasis" ),
                                           numericInput("tmb_axis_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                           numericInput("tmb_axis_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                           selectInput(inputId = "tmb_test",label = "Choose the test to carry out",choices = c("t.test","wilcox.test"),selected = "t.test"),
                                           actionButton(inputId = "tmb_calculate",label = "Click to carry out TMB analysis"),
                                           br(),
                                           br(),
                                           dropMenu(
                                             actionButton("go3", "info",icon = icon("info")),
                                             "If the samples are sequenced at a higher depth compared to TCGA, use this tab to decide on a optimal",
                                             strong("VAF"),
                                             "cutoff so that the TMBs are comparable",
                                             theme = "light-border",
                                             placement = "bottom",
                                             arrow = T
                                           ))
                                 
                                 ),
                    mainPanel(width = 9,
                                h3("Data associated with the plot"),
                                downloadButton(outputId = "download_tmb_table"),
                                dataTableOutput(outputId = "tmb_table") %>% withSpinner(color="black"),
                                h3("TMB difference at several VAF cutoffs"),
                                downloadButton(outputId = "download_tmb_plot"),
                                plotOutput(outputId = "tmb_plot") %>% withSpinner(color="black")
                              )
                  
           ),
           ################################################################
           tabPanel("Mutational Signatures",
                    mainPanel(width = 12,
                              tabsetPanel(
                                tabPanel("Signature Detection",
                                         
                                         sidebarPanel(width = 3,
                                                     
                                                      fluidPage(numericInput(inputId = "sig_vaf_cutoff",label = "Select the VAF cutoff for Signature Analysis",min = 0,max = 0.45,value = 0,step = 0.05),
                                                                numericInput(inputId = "sig_nmax",label = "Select the maximum number of signatures to try",min = 3,max = 10,value = 3,step = 1),
                                                                numericInput(inputId = "sig_bootstrap",label = "Select the bootstrap value",min = 2,max = 10,value = 2,step = 1),
                                                                actionButton(inputId = "sig_run",label = "Click to run Cophnetic correlation analysis"),
                                                                numericInput(inputId = "num_signatures",label = "Choose the optimal number signatures, after observing the Cophnetic plot",value = 2,min = 2,max = 10,step = 1),
                                                                actionButton(inputId = "extract_sigs",label = "Click to run analyses based on chosen number of signatures"),
                                                                h3("Signature Plot  Aesthetics"),
                                                                numericInput(inputId = "sig_sig_plot_font",label = "Select the font size",value = 1.2,min = 0.1,max = 10,step = 0.05),
                                                                numericInput(inputId = "sig_sig_plot_font_title",label = "Select the font size for title",value = 1.3,min = 0.1,max = 10,step = 0.05),
                                                                
                                                                br(),
                                                                br(),
                                                                dropMenu(
                                                                  actionButton("go4", "info",icon = icon("info")),
                                                                  "Mutational signature detection is based on the non-negative matrix factorization. Stability of the signatures are assessed via",
                                                                  strong("cophnetic correlation coefficient"),
                                                                  ". Detailed information for each method can be found in the",
                                                                  tags$a(href = "https://pubmed.ncbi.nlm.nih.gov/15016911/", "article."),
                                                                  theme = "light-border",
                                                                  placement = "bottom",
                                                                  arrow = T
                                                                ))
                                                      
                                                      ),
                                         mainPanel(
                                           tabsetPanel(
                                             tabPanel("Cophnetic Correlation vs Signature Number",
                                                      downloadButton("download_sig_coph"),
                                                      plotOutput(outputId = "sig_coph") %>% withSpinner(color="black")
                                                      ),
                                             tabPanel("Identified signatures and Their Cosine similarity to SBS signatures",
                                                      downloadButton("download_sig_sigs"),
                                                      plotOutput(outputId = "sig_sigs") %>% withSpinner(color="black")
                                                      )
                                           )
                                           
                                  
                                         )
                                        
                                         ),
                                tabPanel("Signature Contribution",
                                         sidebarPanel(
                                           fluidPage(h3("Contribution Plot  Aesthetics"),
                                                     numericInput(inputId = "sig_cont_plot_font",label = "Select the font size",value = 1.2,min = 0.1,max = 10,step = 0.05),
                                                     
                                         )),
                                         mainPanel(h3("Table of Contributions"),
                                                   downloadButton("download_sig_conts_table"),
                                                   dataTableOutput(outputId = "sig_conts_table")%>% withSpinner(color="black"),
                                                   h3("Contribution of Signatures"),
                                                   downloadButton("download_sig_conts"),
                                                   plotOutput(outputId = "sig_conts") %>% withSpinner(color="black"))
                                         
                                         ),
                                tabPanel("Clinical Correlates",
                                         tabsetPanel(
                                           tabPanel(title = "TMB vs Contributions",
                                                    sidebarPanel(
                                                      numericInput("clin_cor_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_cor_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_cor_ylim","Choose the Y-axis limit",value = 1000,min = 2,max = 1000,step = 1),
                                                      
                                                    ),
                                                    mainPanel(
                                                      downloadButton("download_sig_tmb_cor"),
                                                      plotOutput(outputId = "sig_tmb_cor") %>% withSpinner(color="black"),
                                                    )
                                                    
                                           ),
                                           tabPanel(title = "ANOVA of contributions",
                                                    sidebarPanel(
                                                      numericInput("clin_anova_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_anova_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                    ),
                                                    mainPanel(
                                                      downloadButton("download_sig_anova"),
                                                      plotOutput(outputId = "sig_anova") %>% withSpinner(color="black")
                                                    )
                                                   
                                           ),
                                           tabPanel(title = "Correlation with a Continous Variable",
                                                    sidebarPanel(
                                                      selectizeInput(inputId = "clin_cont_col","Select the column for analysis",choices =NULL),
                                                      numericInput("clin_cont_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_cont_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_cont_ylim","Choose the Y-axis limit",value = 1000,min = 2,max = 1000,step = 1),
                                                      actionButton(inputId = "clin_cont_run",label = "Click to run")
                                                    ),
                                                    mainPanel(
                                                      downloadButton("download_sig_cont_cor"),
                                                      plotOutput(outputId = "sig_cont_cor") %>% withSpinner(color="black")
                                                    )
                                                    
                                           ),
                                           
                                             tabPanel(title = "ANOVA by a character column",
                                                      sidebarPanel(
                                                        selectizeInput(inputId = "clin_char_col","Select the column for analysis",choices =NULL),
                                                        numericInput("clin_char_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                        numericInput("clin_char_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                        actionButton(inputId = "clin_char_run",label = "Click to run")
                                                      ),
                                                      mainPanel(
                                                        downloadButton(outputId = "download_sig_char_cor"),
                                                        plotOutput(outputId = "sig_char_cor") %>% withSpinner(color="black")
                                                      )
                                             
                                           ),
                                           tabPanel(title = "Pairwise Comparison",
                                                    sidebarPanel(
                                                      selectizeInput(inputId = "clin_pair_col","Select the column for analysis",choices =NULL),
                                                      selectizeInput(inputId = "clin_pair_f","Select the first group",choices =NULL),
                                                      selectizeInput(inputId = "clin_pair_s","Select the second group",choices =NULL),
                                                      numericInput("clin_pair_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_pair_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                      actionButton(inputId = "clin_pair_run",label = "Click to run")
                                                    ),
                                                    mainPanel(
                                                      downloadButton("download_sig_pair_cor"),
                                                      plotOutput(outputId = "sig_pair_cor") %>% withSpinner(color="black")
                                                    )
                                                    
                                           ),
                                           tabPanel(title = "Pairwise Comparison Within Signatures",
                                                    sidebarPanel(
                                                      selectizeInput(inputId = "clin_pair_col_sig","Select the column for analysis",choices =NULL),
                                                      selectizeInput(inputId = "clin_pair_f_sig","Select the first signature",choices =NULL),
                                                      selectizeInput(inputId = "clin_pair_s_sig","Select the second signature",choices =NULL),
                                                      numericInput("clin_pair_label_font_sig","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                      numericInput("clin_pair_tick_font_sig","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                      actionButton(inputId = "clin_pair_run_sig",label = "Click to run")
                                                    ),
                                                    mainPanel(
                                                      downloadButton("download_sig_pair_cor_sig"),
                                                      plotOutput(outputId = "sig_pair_cor_sig") %>% withSpinner(color="black")
                                                    )
                                                    
                                           )
                                           
                                         )
                                        
                                         
                                
                                
                              
                              
                    )))),
           tabPanel("Clinical enrichment",
                    sidebarPanel(width = 3,
                      fluidPage(
                        selectizeInput(inputId = "clin_enrich_column","Select the column you would like to carry the enrichment analysis on",choices =NULL,multiple=F),
                        h3("Enrichment Plot Aesthetics"),
                        numericInput(inputId = "clin_enrich_pval","Eneter the p-value threshold",value = 0.001,min = 0,max = 1,step = 0.01),
                        numericInput(inputId = "clin_enrich_OR","Eneter the odds ratio threshold",value = 0.001,min = 0,max = 1,step = 0.01),
                        numericInput(inputId = "clin_enrich_annotation_font","Eneter the annotation font size",value = 0.8,min = 0,max = 20,step = 0.1),
                        numericInput(inputId = "clin_enrich_gene_font","Eneter the gene font size",value = 0.8,min = 0,max = 20,step = 0.1),
                        numericInput(inputId = "clin_enrich_legend_font","Eneter the legend font size",value = 0.8,min = 0,max = 20,step = 0.1),
                        actionButton("clin_enrich_run",label = "click to run")

                      )
                    ),
                    mainPanel(width = 9,
                        tabsetPanel(
                          tabPanel(title = "Pairwise Enrichment Results Table",
                                    downloadButton(outputId = "download_pairwise_enrich"),
                                              dataTableOutput(outputId = "clin_enrich_pairwise")%>% withSpinner(color="black")

                                   ),
                          tabPanel(title = "Groupwise Enrichment results",
                                   downloadButton(outputId = "download_groupwise_enrich"),
                                             dataTableOutput(outputId = "clin_enrich_groupwise")%>% withSpinner(color="black")

                                   ),
                          tabPanel(title = "Clinical Enrichment Plot",
                                   mainPanel(downloadButton(outputId = "download_enrich_plot"),
                                             plotOutput(outputId = "enrich_plot",)%>% withSpinner(color="black"))

                                   )
                        )

                    )
                    
                    ),
           tabPanel("Comparison with TCGA",
                    mainPanel(width = 12,
                              tabsetPanel(
                                tabPanel(title = "Differentially Mutated Genes",
                                         sidebarPanel(width = 3,
                                                      fluidPage(numericInput(inputId = "vaf_cutoff_compare",label = "Select the VAF cutoff for comparison",min = 0,max = 1,value = 0.1,step = 0.05),
                                                                actionButton(inputId = "vaf_cutoff_run",label = "Click to run"),
                                                                h3("Forest Plot Aesthetics"),
                                                                numericInput(inputId = "vaf_cutoff_fdr",label = "Select the FDR cutoff",value = 0.05,min = 0,max = 1,step = 0.05),
                                                                numericInput(inputId = "vaf_cutoff_gene_size",label = "Select the gene font size",value = 0.8,min = 0.05,max = 20,step = 0.05),
                                                                numericInput(inputId = "vaf_cutoff_title_size",label = "Select the title font size",value = 1.2,min = 0.05,max = 20,step = 0.05),
                                                                )

                                                      ),
                                         mainPanel(width = 9,
                                                   tabsetPanel(
                                                     tabPanel(title = "Fisher's exact test results",
                                                              downloadButton("download_diffmut_table"),
                                                              dataTableOutput(outputId = "diffmut_table")%>% withSpinner(color="black"),
                                                     ),
                                                     tabPanel("Forest Plot",
                                                              downloadButton("download_diffmut_forest"),
                                                              plotOutput(outputId = "diffmut_forest") %>% withSpinner(color="black")
                                                              )
                                                   )
                                                   )




                                ),
                                tabPanel("Cancer Driver Genes Based on Positional Clustering",
                                         sidebarPanel(width = 3,
                                                      fluidPage( numericInput(inputId = "oncoclust_vaf_cutoff",label = "Select the VAF cutoff for comparison",min = 0,max = 1,value = 0.1,step = 0.05),
                                                                 numericInput(inputId = "oncoclust_fdr_cutoff",label = "Select the FDR cutoff for OncoCLUST",min = 0,max = 1,value = 0.001,step = 0.05),
                                                                 numericInput("onco_axis_label_font","Choose axis label font size",value = 12,min = 2,max = 30,step = 1),
                                                                 numericInput("onco_axis_tick_font","Choose axis tick font size",value = 12,min = 2,max = 30,step = 1),
                                                                 actionButton(inputId = "oncoclust_run",label = "Click to run"),
                                                                 h3("Co-lollipop plot aesthetics"),
                                                                 numericInput("oncoclust_lollipop_angle","Select the label angle",min = 0,max=90,step = 1,value = 0),
                                                                 numericInput("oncoclust_lollipop_label_text_size","Select the label font size",min = 0,max=20,step = 0.1,value = 3),
                                                                 numericInput("oncoclust_lollipop_axis_text_size","Select the axis text size",min = 0,max=20,step = 0.1,value = 1),
                                                                 numericInput("oncoclust_lollipop_head_size","Select the lollipop head size",min = 0,max=20,step = 0.1,value = 1.2),
                                                                 numericInput("oncoclust_lollipop_legend_text_size","Select the legend text",min = 0,max=20,step = 0.1,value = 1.2),
                                                                 )

                                                      ),
                                         mainPanel(width = 9,
                                                   tabsetPanel(
                                                     tabPanel("Uploaded cohort OncoCLUST results",
                                                              downloadButton(outputId = "download_oncoclust_turk"),
                                                              dataTableOutput(outputId = "oncoclust_turk")%>% withSpinner(color="black")
                                                              ),
                                                     tabPanel("TCGA cohort OncoCLUST results",
                                                                downloadButton(outputId = "download_oncoclust_tcga"),
                                                                dataTableOutput(outputId = "oncoclust_tcga")%>% withSpinner(color="black")
                                                              ),
                                                     tabPanel("OncoCLUST plot",
                                                              downloadButton(outputId = "download_oncoclust_plot"),
                                                              plotlyOutput(outputId = "oncoclust_plot") %>% withSpinner(color="black")
                                                              ),
                                                     tabPanel("Co-lollipop plot",
                                                              selectizeInput(inputId = "oncoclust_lolipop_genes","Select Genes",choices = NULL,multiple=F),
                                                              downloadButton(outputId = "download_oncoclust_lollipop2"),
                                                              plotOutput(outputId = "oncoclust_lollipop2") %>% withSpinner(color="black")
                                                              )
                                                   )
                                                   ),


                                ),

                                tabPanel("Somatic Interactions",
                                         sidebarPanel(width = 3,
                                                      fluidPage( numericInput(inputId = "somint_vaf_cutoff",label = "Select the VAF cutoff",min = 0,max = 1,value = 0.1,step = 0.05),
                                                                 actionButton(inputId = "somint_run",label = "Click to run"),
                                                                 h3("Somatic interaction plot Aesthetics"),
                                                                 numericInput(inputId = "somint_fontsize",label = "Select the fontsize for somatic interaction plot",min = 0,max = 20,value = 0.8,step = 0.1),
                                                                 numericInput(inputId = "somint_asteriks_size",label = "Select the size of the significance asteriks",min = 0,max = 20,value = 2,step = 0.1),
                                                                 numericInput(inputId = "somint_legend_fontsize",label = "Select the fontsize for the legend",min = 0,max = 20,value = 0.9,step = 0.1),
                                                                 selectInput(inputId = "somint_if_show_sum",label = "Show the sum of events for each gene",list("Yes"=T,"No"=F),selected = "No")
                                                                 )
                                         
                                                      ),
                                         mainPanel(width = 9,
                                                    tabsetPanel(
                                                      tabPanel(title = "Uploaded cohort somatic interaction table",
                                                               downloadButton("download_somint_turk_table"),
                                                               dataTableOutput(outputId = "somint_turk_table")%>% withSpinner(color="black")
                                                      ),
                                                      tabPanel(title = "Uploaded cohort somatic interaction plot",
                                                               downloadButton("download_somint_turk_plot"),
                                                               plotOutput(outputId = "somint_turk_plot") %>% withSpinner(color="black")
                                                               ),
                                                      tabPanel(title = "TCGA cohort somatic interaction table",
                                                               downloadButton("download_somint_tcga_table"),
                                                               dataTableOutput(outputId = "somint_tcga_table")%>% withSpinner(color="black")
                                                               ),
                                                      tabPanel("TCGA cohort somatic interaction plot",
                                                               downloadButton("download_somint_tcga_plot"),
                                                               plotOutput(outputId = "somint_tcga_plot") %>% withSpinner(color="black")
                                                               ),
                                                      tabPanel("Interactions unique to TCGA",
                                                               downloadButton("download_somint_uniq_tcga"),
                                                               dataTableOutput(outputId = "somint_uniq_tcga")%>% withSpinner(color="black")
                                                               ),
                                                      tabPanel("Interactions unique to uploaded cohort",
                                                               downloadButton("download_somint_uniq_turk"),
                                                               dataTableOutput(outputId = "somint_uniq_turk")%>% withSpinner(color="black")
                                                               ),
                                                      tabPanel("Interactions common to both cohorts",
                                                               downloadButton("download_somint_both"),
                                                               dataTableOutput(outputId = "somint_both")%>% withSpinner(color="black")
                                                               ),
                                                      tabPanel("Interactions in opposite direction",
                                                               downloadButton("download_somint_opopsite"),
                                                               dataTableOutput(outputId = "somint_opopsite")%>% withSpinner(color="black")
                                                               )
                                                    )
                                                   )




                                ),
                                tabPanel("Comparison of Mutational Signatures",

                                         tabsetPanel(
                                           tabPanel("Cophnetic Correlation vs Signature Number, Uploaded cohort",
                                                    downloadButton("download_comp_turk_coph"),
                                                    plotOutput(outputId = "comp_turk_coph") %>% withSpinner(color="black")
                                           ),
                                           tabPanel("Cophnetic Correlation vs Signature Number, TCGA cohort",
                                                    downloadButton("download_comp_tcga_coph"),
                                                    plotOutput(outputId = "comp_tcga_coph") %>% withSpinner(color="black")
                                                    ),
                                           tabPanel("Identified signatures and Their Cosine similarity to SBS signatures, Uploaded cohort",
                                                    downloadButton("download_comp_turk_sigs"),
                                                    plotOutput(outputId = "comp_turk_sigs") %>% withSpinner(color="black")
                                                    ),
                                           tabPanel("Identified signatures and Their Cosine similarity to SBS signatures, TCGA cohort cohort",
                                                    downloadButton("download_comp_tcga_sigs"),
                                                    plotOutput(outputId = "comp_tcga_sigs") %>% withSpinner(color="black")
                                                    ),
                                           tabPanel("Cosine similarity between TCGA and Uploaded cohort signatures",
                                                    downloadButton("download_comp_cosine"),
                                                    plotOutput(outputId = "comp_cosine") %>% withSpinner(color="black")
                                                    )
                                         )

                                ),
                                
                                                    
                                                                tabPanel("Comaprison of Pathways",
                                                                         tabsetPanel(
                                                                           tabPanel("Summary",
                                                                                    sidebarPanel(width = 3,
                                                                                                 fluidPage(numericInput(inputId = "path_summary_vaf_cutoff",label = "Select the VAF cutoff for comparison",min = 0,max = 0.45,value = 0.1,step = 0.05),
                                                                                                           selectInput(inputId = "path_summary_choose",label = "Choose the dataset to display",choices = c("Uploaded Dataset","TCGA"),selected = "Uploaded Dataset",multiple = F))

                                                                                    ),

                                                                                    mainPanel(width = 9,
                                                                                              tabsetPanel(
                                                                                                tabPanel(title = "Summary Table",
                                                                                                         downloadButton("download_path_summary_table"),
                                                                                                         dataTableOutput(outputId = "path_summary_table")%>% withSpinner(color="black")
                                                                                                ),
                                                                                                tabPanel(title = "Summary Plot",
                                                                                                         downloadButton(outputId = "download_path_summary_plot"),
                                                                                                         plotOutput(outputId = "path_summary_plot") %>% withSpinner(color="black")
                                                                                                )
                                                                                              )
                                                                                    )


                                                                           ),
                                                                           tabPanel(title = "Comparison",
                                                                                    tabsetPanel(
                                                                                      tabPanel(title = "Fisher's exact test table",
                                                                                               downloadButton("download_path_comparison_table"),
                                                                                               dataTableOutput(outputId = "path_comparison_table")%>% withSpinner(color="black")
                                                                                      ),
                                                                                      tabPanel(title = "Forest plot",
                                                                                               sidebarPanel(
                                                                                                 numericInput(inputId = "path_comparison_fdr",label = "Select the FDR cutoff",value = 0.05,min = 0,max = 1,step = 0.05),
                                                                                                 numericInput(inputId = "path_comparison_gene_size",label = "Select the gene font size",value = 0.8,min = 0.05,max = 20,step = 0.05),
                                                                                                 numericInput(inputId = "path_comparison_title_size",label = "Select the title font size",value = 1.2,min = 0.05,max = 20,step = 0.05),
                                                                                               ),
                                                                                              mainPanel(
                                                                                                downloadButton("download_path_comparison_plot"),
                                                                                                plotOutput(outputId = "path_comparison_plot") %>% withSpinner(color="black")
                                                                                              )
                                                                                               
                                                                                      )
                                                                                    )

                                                                           ),
                                                                           tabPanel("Pathways",
                                                                                    sidebarPanel(
                                                                                      selectInput(inputId = "path_display_choose_dataset",label = "Choose the dataset to display",choices = c("Uploaded Dataset","TCGA"),selected = "Uploaded Dataset",multiple = F),
                                                                                      selectInput(inputId = "path_display_choose_pathway",label = "Choose the dataset to display",choices = all_pathways,multiple = F),
                                                                                      selectInput(inputId = "path_display_show_tsb",label = "Show the tumor sample barcode",choices = list("Yes"=T,"No"=F),selected = "Yes",multiple = F),
                                                                                      numericInput(inputId = "path_display_font_size",label = "Select the fontsize",value = 0.6,min = 0,max = 20,step = 0.1),
                                                                                      numericInput(inputId = "path_display_sample_font_size",label = "Select the tumor sample barcode  font size",value = 0.6,min = 0,max = 20,step = 0.1),
                                                                                      
                                                                                    ),
                                                                                   mainPanel(
                                                                                     downloadButton("download_path_display_plot"),
                                                                                     plotOutput(outputId = "path_display_plot") %>% withSpinner(color="black")
                                                                                   )
                                                                                    
                                                                           )
                                                                         )
                                                                )

                                                              )
                                                  


                                           )
                                          
                                           
                                         )
                                )

                              

                    
                    
           
           ###############################
           
