#LoadPackage
library(tidyverse)
library(shiny)
library(rsconnect)
library(plotly)
library(shinythemes)
library(ggpubr)
library(DT)
library(flextable)
library(shinydashboard)
library(shinyWidgets)
library(data.table)
library(ggsci)
library(rmarkdown)
library(RColorBrewer)
library(randomcoloR)
library(ComplexHeatmap)
library(genekitr)
library(PhenoExam)
library(ggrepel)
library(ggh4x)
library(VennDiagram)
library(rstatix)
library(Hmisc)
library(pheatmap)
library(AnnotationDbi)
library(clusterProfiler)
library("org.Hs.eg.db")
library(shinycssloaders)


#Load data
#Data
#AstraZeneca PheWAS Portal
az_gene_binary_raw <- fread("Data/AZ_Gene-level_Binary_1e-08.csv")
az_gene_continuous_raw <- fread("Data/AZ_Gene-level_Continuous_1e-08.csv")
az_var_binary_raw <- fread("Data/AZ_Variant-level_Binary_1e-08_gnomAD_clinvar.csv")
az_var_continuous_raw <- fread("Data/AZ_Variant-level_Continuous_1e-08_gnomAD_clinvar.csv")
#FinnGen 
finngen_raw <- fread("Data/Finngen11_Result_5e-8_SmallestPvalueGene_gnomAD_clinvar.csv")
#GWAS Catalog
gwas_raw <- fread("Data/gwas_catalog_v1.0.2-associations_e112_r2024-07-08.tsv")
gwas_raw <- gwas_raw %>% filter(PVALUE_MLOG > abs(log10(5e-8)))

gtex_raw <- fread("Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct")
#EFO MAP
efo_raw <- fread("Data/trait_mappings.txt")
#gnomAD GWAS Catalog
gnomad_raw <- fread("Data/gnomAD_GWAS.csv")
#Clinvar
clinvar_raw <- fread("Data/clinvar_extracted_data.txt.gz")
example <- fread("Data/ExampleGeneList.csv")

#############################
#function
plot_exception <-function(
    ...,
    sep=" ",
    type=c("message","warning","cat","print"),
    color="auto",
    console=TRUE,
    size = 6){
  type=match.arg(type)
  txt = paste(...,collapse=sep)
  if(console){
    if(type == "message") message(txt)
    if(type == "warning") warning(txt)
    if(type == "cat") cat(txt)
    if(type == "print") print(txt)
  }
  if(color =="auto") color <- if(type == "cat") "black" else "red"
  if(txt == "warning") txt <- paste("warning:",txt)
  print(ggplot2::ggplot() +
          ggplot2::geom_text(ggplot2::aes(x=0,y=0,label=txt),color=color,size=size) +
          ggplot2::theme_void())
  invisible(NULL)
}


plotly_exception <- function(
    ...,
    sep = " ",
    type = c("message", "warning", "cat", "print"),
    color = "auto",
    console = TRUE,
    size = 6
) {
  # Match argument for type
  type <- match.arg(type)
  
  # Combine the inputs into a single text string
  txt <- paste(..., collapse = sep)
  
  # Print the message to the console if console = TRUE
  if (console) {
    switch(type,
           message = message(txt),
           warning = warning(txt),
           cat = cat(txt, "\n"),
           print = print(txt)
    )
  }
  
  # Set color automatically based on type if color is "auto"
  if (color == "auto") {
    color <- if (type == "cat") "black" else "blue"
  }
  
  # Modify the warning text if necessary
  if (type == "warning") txt <- paste("warning:", txt)
  
  # Create the ggplot object
  p <- ggplot2::ggplot() +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label = txt), color = color, size = size) +
    ggplot2::theme_void()
  
  # Convert ggplot to plotly and print the plot
  p_interactive <- plotly::ggplotly(p)
  print(p_interactive)
  #invisible(NULL)  # Return NULL invisibly
}



infoBtn <- function(id) {
  actionButton(id,
               label = "",
               icon = icon("question"),
               style = "info",
               size = "extra-small",
               class='btn action-button btn-info btn-xs shiny-bound-input'
  )
}


#Define ui
ui <- navbarPage(inverse = TRUE,
                 h4(strong("GeneSetPheno")), 
                 #theme = shinytheme("spacelab"),
                 theme = shinytheme("cosmo"),
                 id = "panels_main",
                 tabPanel(h4(icon("home")), 
                          fluidPage(
                            br(),
                            br(),
                            fluidRow(
                              column(12,
                                        
                                     verbatimTextOutput("sessioninfo"),
                                     fileInput("file1", strong("Upload Gene Set of Interest CSV File"), accept = ".csv",width = "1000px"),
                                     actionButton("RunAnalysis",strong("Run Gene Set Analysis"), icon = icon("play")),
                                     # br(),
                                     # br(),
                                     downloadButton(outputId = "report",label = strong("Generate report")),
                                     downloadButton("downloadData", strong("Download Example Input File"))
                                     
                              ), 
                              
                              align = "center"
                              
                              
                            ),
                            
                            br(),
                            br(),
                            
                            p(strong("About"),style = "font-size: 17pt"),
                            p(style="text-align: justify;","TThis web application integrates data from public databases such as the ", 
                              a(strong("AstraZeneca PheWAS Portal, "),style = "color:#337ab7", href = "https://azphewas.com"),
                              a(strong("FinnGen, "),style = "color:#337ab7", href = "https://r11.finngen.fi"),
                              a(strong("GWAS Catalog, "),style = "color:#337ab7", href = "https://www.ebi.ac.uk/gwas"),
                              a(strong("Human Phenotype Ontology, "),style = "color:#337ab7", href = "https://hpo.jax.org"),
                              a(strong("gnomAD, "),style = "color:#337ab7", href = "https://gnomad.broadinstitute.org"),
                              a(strong("ClinVar, "),style = "color:#337ab7", href = "https://www.ncbi.nlm.nih.gov/clinvar/"),
                              
                              "and ",
                              a(strong("GTEx, "),style = "color:#337ab7", href = "https://gtexportal.org/home/"),
                           "to generate visual summaries of gene, genetic variant, phenotype, and association information.",
                            style = "font-size: 13pt"),
                         
                            br(),
                            br(),
                            img(src='IntroMindMap.png', align = "center", height="80%", width="80%"),
                            br(),
                            br(),
                            
                            p(strong("Data Sources"),style = "font-size: 13pt"),
                            p(style="text-align: justify;",
                              a(strong("AstraZeneca PheWAS Portal: "),style = "color:#337ab7", href = "https://azphewas.com"), "Gene-phenotype associations, the largest and most comprehensive exome-wide genotype-phenotype dataset, rare-variant genetic association data. Provides both gene and variant-level phenome-wide association statistics (PheWAS) using the exome sequences of the UK Biobank participants and considered ~17K binary and ~1.4K quantitative phenotypes.
",
                              br(),
                              a(strong("FinnGen: "),style = "color:#337ab7", href = "https://r11.finngen.fi"), "The FinnGen research project is an academic industrial collaboration aiming to identify genotype-phenotype correlations in the Finnish founder population designed to develop the potential of these resources to serve medicine initiate and enrich drug discovery programs.",
                              br(),
                              a(strong("GWAS Catalog: "),style = "color:#337ab7", href = "https://www.ebi.ac.uk/gwas"), "The NHGRI-EBI Catalog of human genome-wide association studies.",
                              br(),
                              a(strong("Human Phenotype Ontology (HPO): "),style = "color:#337ab7", href = "https://hpo.jax.org"), "A standardized vocabulary of phenotypic abnormalities encountered in human disease.",
                              br(),
                              a(strong("gnomAD: "),style = "color:#337ab7", href = "https://gnomad.broadinstitute.org"), "Provides information on human genetic variation in healthy individuals across a diverse range of genetic ancestry groups.",
                              br(),
                              a(strong("ClinVar: "),style = "color:#337ab7", href = "https://www.ncbi.nlm.nih.gov/clinvar/"), "Report the relationships among human variations and phenotypes, with supporting evidence.",
                              br(),
                              a(strong("Genotype-Tissue Expression (GTEx): "),style = "color:#337ab7", href = "https://gtexportal.org/home/"), "A data resource and tissue bank established to study the relationships between variants and gene expression across several tissue types and different people.",
                              br(),
                              
                              style = "font-size: 10pt"),
                            br(),
                            br(),
                            div(p(strong("This was developed by Jiru Han with input from Melanie Bahlo and other "), a("Bahlo Lab members", href = "https://www.wehi.edu.au/people/melanie-bahlo/372/melanie-bahlo-lab-team")), 
                                p(strong("For any queries or suggestions, please contact Jiru Han: han.ji@wehi.edu.au")), 
                                style="text-align: right;")
                          )
                 ),
                 
                
                 #Gene Summary
                 tabPanel(h4(strong("Gene Summary")),
                          fluidPage(
                            
                            br(),
                            p(strong("Gene Summary"),style = "font-size: 15pt"),
                            p("This module enables the extraction of comprehensive gene set information, including Entrez, Ensembl, and Uniprot IDs, genomic locations, gene function summaries, gene sequences, transcript counts, gene biotypes, and additional details. Additionally, it summarizes significant gene-phenotype associations from databases such as AZPheWAS, GWAS Catalog, and FinnGen, identifying the number of associated genes and presenting the results in both summary tables and figures.",style = "font-size: 14pt"),
                            br(),
                            actionButton("ViewGeneInfo",strong("Generate Gene Summary")),
                            br(),
                            hr(),
                            
                            mainPanel(tabsetPanel(
                              tabPanel(strong("Gene Info"),
                                       br(),
                                       p("Detailed gene set information",style = "font-size: 14pt"),
                                       
                                       style = "background: white",
                                       fluidRow(
                                       
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(DT::dataTableOutput("p1_dt1")),
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(plotOutput("p1_plot1"))
                                       )))),
                              
                              tabPanel(strong("Summary Table of Significant Gene–Phenotype Associations"),
                                       br(),
                                       p("A table summarizing significant gene-phenotype associations from multiple databases",style = "font-size: 14pt"),
                                       style = "background: white",
                                       fluidRow(
                                         
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(uiOutput("p1_dt2"))
                                       )))),
                              
                              
                              tabPanel(strong("Summary Plot of Significant Gene–Phenotype Associations"),
                                       br(),
                                       p("A summary plot visualizing these associations for the gene sets",style = "font-size: 14pt"),
                                       fluidRow(
                                                
                                                fluidRow(column(12,
                                                                br(),
                                                                hr(),
                                                                shinycssloaders::withSpinner(plotOutput("p1_plot2",width = "100%"))
                                                               
                                                ))))
                              
                            ))
                            
                            
                          )
                 ),
                 
                 
                 
                 
                 #Gene-Phenotype Associations
                 tabPanel(h4(strong("Gene-Phenotype Association")),
                          fluidPage(
                            br(),
                            p(strong("Gene-Phenotype Association"),style = "font-size: 15pt"),
                            p("This module in GeneSetPheno provides a detailed summary of gene-phenotype associations. Significant associations from AZPheWAS, GWAS Catalog, and FinnGen databases are summarized for the gene set. This information is presented in an interactive table and figure, offering comprehensive details for each gene, including all associated phenotypes, phenotype categories, and variants.",style = "font-size: 14pt"),
                            
                            actionButton("ViewGenePhenotypeAsso",strong("Generate Gene-Phenotype Associations")),
                            br(),
                            hr(),

                            
                            mainPanel(tabsetPanel(
                              tabPanel(strong("Summary Plot of Significant Gene–Phenotype Associations"),
                                       br(),
                                       p("A summary plot that displays significant gene–phenotype associations for the gene sets. Users can easily access detailed information about each gene’s phenotype associations across databases, such as phenotype categories, by hovering over a gene or database.",style = "font-size: 14pt"),
                                       
                                       style = "background: white",
                                       fluidRow(
                                       
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(plotlyOutput("p2_plot1"))
                                       )))),

                              
                              tabPanel(strong("Summary Table of Significant Gene–Phenotype Associations"),
                                       br(),
                                       p("A table that displays significant gene–phenotype associations for the gene sets. Users can access all relevant phenotype categories, phenotypes, and variants across multiple databases for each gene, offering a clear and comprehensive overview of significant gene–phenotype associations.",style = "font-size: 14pt"),
                                       
                                       fluidRow(
                                         
                                                fluidRow(column(12,
                                                                br(),
                                                                hr(),
                                                                shinycssloaders::withSpinner(DT::dataTableOutput("p2_dt1"))
                                                ))))

                            ))


                          )
                 ),


                 
                 #Variant-Phenotype Associations
                 tabPanel(h4(strong("Variant-Phenotype Association")),
                          fluidPage(
                            br(),
                            p(strong("Variant-Phenotype Association"),style = "font-size: 15pt"),
                            p("This module focuses on four key components: a summary table of variant-phenotype associations, AZPheWAS, GWAS Catalog, and Finngen, with the aim of displaying significant associations between genetic variants and various phenotypes from different databases.",style = "font-size: 14pt"),
                            br(),
                            hr(),
                            
                            
                            
                            mainPanel(tabsetPanel(
                              
                              
                              
                              tabPanel(strong("Summary Table of Significant Variant–Phenotype Associations"),
                                       br(),
                                       p("The summary table provides a comprehensive overview of variant associations from various databases, including details like the gene list group, gene symbol, variant, rsID, allele frequency in gnomAD, link-outs to the gnomAD browser, and clinical significance from ClinVar. It also consolidates significant phenotype data from AZPheWAS, GWAS Catalog, and Finngen, along with additional information on phenotype categories.",style = "font-size: 14pt"),
                                       
                                       style = "background: white",
                                       fluidRow(column(12,
                                                       br(),
                                                       actionButton("ViewVarTable",strong("Generate Variant-Phenotype Associations Table"))
                                       ),
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                      shinycssloaders::withSpinner(DT::dataTableOutput("p3_dt1"))
                                                       
                                       )))
                                       
                              ),
                              
                              tabPanel(strong("AstraZeneca PheWAS Portal"),
                                       style = "background: white",
                                       tabsetPanel(
                                         tabPanel(strong("Phenotypic Profile Clustering"),
                                                  br(),
                                                  p("Phenotypic profile clustering is performed by grouping genes according to their associations with upper-level phenotype categories. For each category, significant variant-phenotype associations for each gene are assessed, followed by hierarchical clustering. This analysis highlights genes with similar phenotype associations. The x-axis represents genes, and the y-axis represents phenotype categories. The color indicates whether there is a significant association of each gene across phenotype categories: red represents a significant association, while white represents no significant association.",style = "font-size: 14pt"),
                                                  
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewAZCluster",strong("View Phenotypic Profile Clustering"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotOutput("p3_plot1"))
                                                                  
                                                  )))
                                                  
                                                  
                                                 
                                         ),
                                         
                                         tabPanel(strong("Phenotype Distribution Overview"),
                                                  br(),
                                                  p("The phenotype distribution overview section visualizes phenotype association patterns across different gene set groups. It aggregates all phenotypes associated with each gene list within each phenotype category, either by counting unique phenotypes or calculating the proportions of genes.",style = "font-size: 14pt"),
                                                  
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  actionButton("ViewAZVar",strong("View Phenotype Distribution"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  shinycssloaders::withSpinner(plotlyOutput("p3_plot2")),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotlyOutput("p3_plot3"))
                                                                  
                                                  )))
                                                  
                                                  
                                                  
                                         ),
                                         
                                         
                                         
                                         tabPanel(strong("Variant-Phenotype Gene Effect"),
                                                  br(),
                                                  p("The mean effect of each gene across various phenotype categories is calculated by averaging the odds ratios or absolute effect sizes from all significant variant-phenotype associations within each gene, for both binary and continuous phenotypes. This reflects the estimated overall effect of each gene within each phenotype category.",style = "font-size: 14pt"),
                                                  
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewAZEffect",strong("View Variant-Phenotype Gene Effect"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotOutput("p3_plot_add1")),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("p3_dt2"))
                                                                  
                                                  )))
                                                  
                                                  
                                                  
                                         )
                                         
                                         
                                         
                                       ) #Close inner tabsetPanel
                              ),

                              
                              
                              tabPanel(strong("GWAS Catalog"),
                                       style = "background: white",
                                       tabsetPanel(
                                         tabPanel(strong("Phenotypic Profile Clustering"),
                                                  br(),
                                                  p("Phenotypic profile clustering is performed by grouping genes according to their associations with upper-level phenotype categories. For each category, significant variant-phenotype associations for each gene are assessed, followed by hierarchical clustering. This analysis highlights genes with similar phenotype associations. The x-axis represents genes, and the y-axis represents phenotype categories. The color indicates whether there is a significant association of each gene across phenotype categories: red represents a significant association, while white represents no significant association.",style = "font-size: 14pt"),
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewGWASCatalogCluster",strong("View Phenotypic Profile Clustering"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotOutput("p3_plot6"))

                                                  )))


                                                  
                                         ),

                                        
                                        


                                         tabPanel(strong("Phenotype Distribution Overview"),
                                                  br(),
                                                  p("The phenotype distribution overview section visualizes phenotype association patterns across different gene set groups. It aggregates all phenotypes associated with each gene list within each phenotype category, either by counting unique phenotypes or calculating the proportions of genes.",style = "font-size: 14pt"),
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewGWASCatalogPheno",strong("View Phenotype Distribution"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotlyOutput("p3_plot7")),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotlyOutput("p3_plot8"))

                                                  )))


                                         )







                                       ) #Close inner tabsetPanel
                              ),

                              
                              
                              
                              
                              tabPanel(strong("FinnGen"),
                                       style = "background: white",
                                       tabsetPanel(
                                         tabPanel(strong("Phenotypic Profile Clustering"),
                                                  br(),
                                                  p("Phenotypic profile clustering is performed by grouping genes according to their associations with upper-level phenotype categories. For each category, significant variant-phenotype associations for each gene are assessed, followed by hierarchical clustering. This analysis highlights genes with similar phenotype associations. The x-axis represents genes, and the y-axis represents phenotype categories. The color indicates whether there is a significant association of each gene across phenotype categories: red represents a significant association, while white represents no significant association.",style = "font-size: 14pt"),
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewFinnGenCluster",strong("View Phenotypic Profile Clustering"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotOutput("p3_plot9")),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotOutput("p3_plot10"))
                                                                  
                                                  )))
                                                  
                                                  
                                                  
                                         ),
                                         
                                         
                                         
                                         
                                         
                                         tabPanel(strong("Phenotype Distribution Overview"),
                                                  br(),
                                                  p("The phenotype distribution overview section visualizes phenotype association patterns across different gene set groups. It aggregates all phenotypes associated with each gene list within each phenotype category, either by counting unique phenotypes or calculating the proportions of genes.",style = "font-size: 14pt"),
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewFinnGenPheno",strong("View Phenotype Distribution"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotlyOutput("p3_plot11")),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotlyOutput("p3_plot12"))
                                                                  
                                                  )))
                                                  
                                                  
                                         ),
                                         
                                         tabPanel(strong("Variant-Phenotype Gene Effect"),
                                                  br(),
                                                  p("The mean effect of each gene across various phenotype categories is calculated by averaging the odds ratios or absolute effect sizes from all significant variant-phenotype associations within each gene, for both binary and continuous phenotypes. This reflects the estimated overall effect of each gene within each phenotype category.",style = "font-size: 14pt"),
                                                  
                                                  style = "background: white",
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  actionButton("ViewFinnGenEffect",strong("View Variant-Phenotype Gene Effect"), icon = icon("play"))
                                                  ),
                                                  fluidRow(column(12,
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(plotOutput("p3_plot_add2")),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  br(),
                                                                  hr(),
                                                                  shinycssloaders::withSpinner(DT::dataTableOutput("p3_dt3"))
                                                                  
                                                  )))
                                                  
                                                  
                                                  
                                         )
                                         
                                         
                                         
                                         
                                         
                                         
                                         
                                       ) #Close inner tabsetPanel
                              )
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                             
                             
                              
                            ))
                            
                            
                          )
                 ),
                 
                 
                 
                 
                 #HPO Phenotype
                 tabPanel(h4(strong("HPO Phenotype")),
                          fluidPage(
                            br(),
                            p(strong("HPO Phenotype"),style = "font-size: 15pt"),
                            p("This module enables visualization of phenotype enrichment results and comparative phenotype analysis of different gene sets using the HPO resource.",style = "font-size: 14pt"),
                            
                            br(),
                            hr(),
                            
                            
                            mainPanel(tabsetPanel(
                              
                              
                              tabPanel(strong("Phenotype Enrichment Visualization"),
                                       br(),
                                       p("Phenotype enrichment analysis on the gene set. The results display the top enriched terms for the gene set in a plot.",style = "font-size: 14pt"),
                                       
                                       style = "background: white",
                                       fluidRow(column(12,
                                                       br(),
                                                       selectInput('GeneListGroup', label = 'Select Group', choices = 'No choices here yet'),
                                                       actionButton("ViewHPOPlot",strong("Phenotype Enrichment Visualization"), icon = icon("play"))
                                       ),
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(plotlyOutput("p4_plot1"))
                                       )))),
                              
                              
                              tabPanel(strong("Phenotype Enrichment Summary Table"),
                                       br(),
                                       p("Phenotype enrichment analysis on the gene set. The results display the top enriched terms for the gene set in a summary table.",style = "font-size: 14pt"),
                                       
                                       style = "background: white",
                                       fluidRow(column(12,
                                                       br(),
                                                       selectInput('GeneListGroupTab', label = 'Select Group', choices = 'No choices here yet'),
                                                       actionButton("ViewHPOTab",strong("View Phenotype Enrichment Summary Table"), icon = icon("play"))
                                       ),
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(DT::dataTableOutput("p4_dt1"))
                                       )))),
                              
                              
                              tabPanel(strong("Comparator Phenotype analysis"),
                                       br(),
                                       p("This generates an interactive graph displaying the relevant phenotypic terms for each gene set, highlighting unique and shared phenotypes, differentiated by color.",style = "font-size: 14pt"),
                                       
                                       style = "background: white",
                                       fluidRow(column(12,
                                                       br(),
                                                       selectInput('Group1', label = 'Select Group1', choices = 'No choices here yet'),
                                                       selectInput('Group2', label = 'Select Group2', choices = 'No choices here yet'),
                                                       actionButton("ViewHPOCompar",strong("Comparator Phenotype Analysis - Gene Set"), icon = icon("play"))
                                       ),
                                       fluidRow(column(12,
                                                       br(),
                                                       hr(),
                                                       shinycssloaders::withSpinner(plotlyOutput("p4_plot2"))
                                       ))))
                              
                              
                              
                              
                              
                            ))
                            
                            
                          )
                 ),
                 
                 
                 
                 #GTEx Analysis
                 tabPanel(h4(strong("GTEx Analysis")),
                          
                          fluidPage(
                            
                            br(),
                            p(strong("Genotype-Tissue Expression (GTEx) Analysis"),style = "font-size: 15pt"),
                            p("This module enables the exploration of tissue-specific gene expression within a gene set using GTEx resources. For this analysis, tissue-specific gene expression data from the GTEx V8 release was downloaded, focusing on the median transcripts per million (TPM) for each gene across various tissues.",style = "font-size: 14pt"),
                            
                            br(),
                            hr(),
                            mainPanel(tabsetPanel(
                            
                              
                              tabPanel(strong("Gene Set Tissue Expression Distribution"),
                                       br(),
                                       p("This feature generates an interactive boxplot that visualizes the distribution of median tissue TPM across GTEx tissues, categorized by gene set groups and distinguished by tissue class through color-coding.",style = "font-size: 14pt"),
                                       
                                                             style = "background: white",
                                                             fluidRow(column(12,
                                                                             br(),
                                                                             actionButton("ViewGTExBoxplot",strong("View Expression Distribution"), icon = icon("play"))
                                                             ),
                                                             fluidRow(column(12,
                                                                             br(),
                                                                             hr(),
                                                                             shinycssloaders::withSpinner(plotlyOutput("p5_plot1"))
                                                             )))),


                                  tabPanel(strong("Gene Set Tissue-Specific Expression Comparison"),
                                           br(),
                                           p("The comparison of gene sets to detect significant differences in gene expression across tissue types. For each tissue type, t-tests were conducted using the rstatix R package to compare gene expression between gene set groups.",style = "font-size: 14pt"),
                                           
                                                             style = "background: white",
                                                             fluidRow(column(12,
                                                                             br(),
                                                                             actionButton("ViewGTExCompar",strong("View Gene Set Expression Comparison"), icon = icon("play"))
                                                             ),
                                                             fluidRow(column(12,
                                                                             br(),
                                                                             hr(),
                                                                             shinycssloaders::withSpinner(plotOutput("p5_plot2"))
                                                             )))),
                              
                              tabPanel(strong("Gene Expression Clustering"),
                                       br(),
                                       p("Clustering of gene expression is conducted by grouping genes based on their expression levels across various tissues. This analysis allows for the identification of genes with similar expression patterns, as well as those exhibiting unique tissue-specific patterns.",style = "font-size: 14pt"),
                                       
                                                             fluidRow(column(12,
                                                                             br(),
                                                                             actionButton("ViewGTExCluster",strong("View Gene Expression Clustering"), icon = icon("play"))),

                                                                      fluidRow(column(12,
                                                                                      br(),
                                                                                      hr(),
                                                                                      shinycssloaders::withSpinner(plotOutput("p5_plot3"))
                                                                      ))))
                                       
                              
                              
                              
                            ))
                            
                            
                          )
                 )
                 

                 
)

# Define server logic to summarize and view selected dataset ----
server <- function(input, output,session) {
  
  #Gene raw data
  gene_raw <- eventReactive(input$RunAnalysis,{
    file <- input$file1
    ext <- tools::file_ext(file$datapath)
    
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    
    data <- read.csv(file$datapath)
    data <- data %>% 
      arrange(Group,Gene)
    data
  }) 
  
  
  #Filter the AZ, GWAS Catalog, and FinnGen data
  p1_dt1_tab <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    gene_info_id <- gene_raw$Gene
    gene_info <- genInfo(gene_info_id, org = "hs", unique = TRUE)
    colnames(gene_info)[1] <- "Gene"
    gene_info <- gene_info %>% 
      full_join(gene_raw) %>% 
      dplyr::select(Group,Gene,gene_name,chr,start,end, width,strand,symbol,hgnc_id, entrezid,ensembl,uniprot,summary,
                    gc_content,gene_biotype,transcript_count,omim)
    gene_info
  })
  
  
  #AstraZeneca PheWAS Portal
  az_gene_binary_flt <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    az_gene_binary_flt <- az_gene_binary_raw %>% 
      filter(Gene %in% gene_raw$Gene)
    az_gene_binary_flt
  })
  
  az_gene_continuous_flt <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    az_gene_continuous_flt <- az_gene_continuous_raw %>% 
      filter(Gene %in% gene_raw$Gene)
    az_gene_continuous_flt
  })

  
  az_var_binary_flt <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    az_var_binary_flt <- az_var_binary_raw %>% 
      filter(Gene %in% gene_raw$Gene)
    az_var_binary_flt
  })
  
  
  az_var_continuous_flt <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    az_var_continuous_flt <-  az_var_continuous_raw %>% 
      filter(Gene %in% gene_raw$Gene)
    az_var_continuous_flt
  })
  

  
  
  #GWAS Catalog
  
  gwas_raw_flt <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    gwas_raw_flt <- gwas_raw %>% 
      mutate(Gene=MAPPED_GENE) %>% 
      separate_rows(Gene, sep = " - ",convert = TRUE) %>% 
      filter(Gene %in% gene_raw$Gene)
    gwas_raw_flt
  })
  

  #FinnGen 
  
  finngen_raw_flt <- eventReactive(input$RunAnalysis,{
    gene_raw <- gene_raw()
    finngen_raw_flt <- finngen_raw %>% 
      mutate(Gene_multi=Gene) %>% 
      separate_rows(Gene, sep = ",",convert = TRUE) %>% 
      filter(Gene %in% gene_raw$Gene)
    finngen_raw_flt
  })
  
  
########################################Home Page#####################################################
  
  
  
  output$downloadData <- downloadHandler(
    
    
    filename = function() {
      paste("GeneSetList-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(example, file, row.names = FALSE)
    }
    
    
  )
  
  
#########################################Page 1 Gene Summary #########################################
  
  p1_dt1_tab <- eventReactive(input$ViewGeneInfo,{
    gene_raw <- gene_raw()
    gene_info_id <- gene_raw$Gene
    gene_info <- genInfo(gene_info_id, org = "hs", unique = TRUE)
    colnames(gene_info)[1] <- "Gene"
    gene_info <- gene_info %>% 
      full_join(gene_raw) %>% 
      dplyr::select(Group,Gene,gene_name,chr,start,end, width,strand,symbol,hgnc_id, entrezid,ensembl,uniprot,summary,
                    gc_content,gene_biotype,transcript_count,omim)
    
    gene_info <- gene_info %>% 
      group_by(Group,Gene) %>% 
      dplyr::slice(1) %>% 
      ungroup() %>% 
      arrange(Group)
    gene_info
  })
  
  # Gene Info
  output$p1_dt1 <- DT::renderDataTable(
    datatable(
      p1_dt1_tab(), 
      extensions = 'Buttons',
      rownames = FALSE,
      options = list(
      dom = 'Blfrtip',
      buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
      columnDefs = list(list(
        targets = c(2,11,12,13),
        render = JS(
          "function(data, type, row, meta) {",
          "return type === 'display' && data.length > 6 ?",
          "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
          "}")
      ))

      )
      , callback = JS('table.page(3).draw(false);')
      )
    
    
  )


  #Summary of Genes with Significant Associations Across Different Databases
  gene_list_summary_tab_num  <- eventReactive(input$ViewGeneInfo,{
    gene_raw <- gene_raw()
    
    az_gene <- unique(c(unique(az_gene_binary_raw$Gene),unique(az_gene_continuous_raw$Gene),unique(az_var_binary_raw$Gene),unique(az_var_continuous_raw$Gene)))
    az_gene_list <- gene_raw %>% 
      filter(Gene %in% az_gene) %>% 
      mutate(AZ="AstraZeneca PheWAS Portal")
    
    if (nrow(az_gene_list)!=0) {
      az_gene_list
    } else {
      az_gene_list <- data.frame(Group=gene_raw$Group,Gene=gene_raw$Gene,AZ=NA)
    }
    
    
    #GWAS Catalog
    gwas_gene <- unlist(strsplit(gwas_raw$MAPPED_GENE, split= " - ", fixed=TRUE))
    gwas_gene <- unique(gwas_gene)
    gwas_gene_list <- gene_raw %>% 
      filter(Gene %in% gwas_gene) %>% 
      mutate(GWASCatalog="GWAS Catalog")
    
    if (nrow(gwas_gene_list)!=0) {
      gwas_gene_list
    } else {
      gwas_gene_list <- data.frame(Group=gene_raw$Group,Gene=gene_raw$Gene,GWASCatalog=NA)
    }
    
    
    #FinnGen 
    finngen_gene <- unlist(strsplit(finngen_raw$Gene, split= ",", fixed=TRUE))
    finngen_gene <- unique(finngen_gene)
    finngen_gene_list <- gene_raw %>% 
      filter(Gene %in% finngen_gene) %>% 
      mutate(FinnGen="FinnGen")
    
    if (nrow(finngen_gene_list)!=0) {
      finngen_gene_list
    } else {
      finngen_gene_list <- data.frame(Group=gene_raw$Group,Gene=gene_raw$Gene,FinnGen=NA)
    }
    
    #Combine Databases
    gene_list_summary_tab <- Reduce(full_join,list(gene_raw, az_gene_list, finngen_gene_list,gwas_gene_list))
    
    #gene_list_summary_pheat <- gene_list_summary_tab
    gene_num <- as.data.frame(table(gene_list_summary_tab$Group))
    colnames(gene_num) <- c("Group", "Gene Number")
    
    if (all(is.na(gene_list_summary_tab$AZ))) {
      az_num <- data.frame(Group=gene_raw$Group,N=0)
      colnames(az_num) <- c("Group", "AstraZeneca PheWAS Portal Significant Association (#N Gene)")
      
    } else {
      az_num<- as.data.frame(table(gene_list_summary_tab$Group,gene_list_summary_tab$AZ)) %>% 
        dplyr::select(-Var2)
      colnames(az_num) <- c("Group", "AstraZeneca PheWAS Portal Significant Association (#N Gene)") 
      
    }
    
    
    
    if (all(is.na(gene_list_summary_tab$GWASCatalog))) {
      gwas_num <- data.frame(Group=gene_raw$Group,N=0)
      colnames(gwas_num) <- c("Group", "GWAS Catalog Significant Association (#N Gene)")
      
    } else {
      gwas_num<- as.data.frame(table(gene_list_summary_tab$Group,gene_list_summary_tab$GWASCatalog)) %>% 
        dplyr::select(-Var2)
      colnames(gwas_num) <- c("Group", "GWAS Catalog Significant Association (#N Gene)")
    }
    
    
    if (all(is.na(gene_list_summary_tab$FinnGen))) {
      finngen_num <- data.frame(Group=gene_raw$Group,N=0)
      colnames(finngen_num) <- c("Group", "Finngen Significant Association (#N Gene)")
      
    } else {
      finngen_num<- as.data.frame(table(gene_list_summary_tab$Group,gene_list_summary_tab$FinnGen)) %>% 
        dplyr::select(-Var2)
      colnames(finngen_num) <- c("Group", "Finngen Significant Association (#N Gene)")
    }
    
    gene_list_summary_tab_num <- Reduce(full_join,list(gene_num, az_num, gwas_num,finngen_num))
    total <- data.frame(Group="Total",
                        `Gene Number`=length(unique(gene_list_summary_tab$Gene)),
                        `AstraZeneca PheWAS Portal Significant Association (#N Gene)` = as.list(gene_list_summary_tab %>% filter(is.na(AZ)==FALSE) %>% summarise(n=length(unique(Gene))))$n,
                        `GWAS Catalog Significant Association (#N Gene)`= as.list(gene_list_summary_tab %>% filter(is.na(GWASCatalog)==FALSE) %>% summarise(n=length(unique(Gene))))$n,
                        `Finngen Significant Association (#N Gene)` = as.list(gene_list_summary_tab %>% filter(is.na(FinnGen)==FALSE) %>% summarise(n=length(unique(Gene))))$n
    )
    colnames(total) <- colnames(gene_list_summary_tab_num)
    gene_list_summary_tab_num <- do.call(rbind,list(gene_list_summary_tab_num,total))
    gene_list_summary_tab_num
    
  })
  
  
  

    output$p1_dt2 <- renderUI({
      ft <- flextable(gene_list_summary_tab_num(),cwidth = 1.5)
      ft <- theme_vanilla(ft)
      ft <- color(ft, part = "footer", color = "#666666")
      ft <- set_caption(ft, caption = "Overview of Genes with Significant Associations Across Different Databases")
      ft %>% 
        autofit() %>%
        htmltools_value()
    })

  
    
    #Gene List Venn diagram
    gene_venn  <- eventReactive(input$ViewGeneInfo,{
      gene_raw <- gene_raw()
      gene_venn <- gene_raw()
      gene_venn
      
    })
    
    
    
    output$p1_plot1 <- renderPlot({
      gene_venn <- gene_venn()
      # Identify the unique groups
      unique_groups <- unique(gene_venn$Group)
      # Create a list of sets based on the Group column
      sets <- split(gene_venn$Gene, gene_venn$Group)
      # Generate a palette of colors based on the number of unique groups
      set.seed(1111)
      n1 <- length(unique(unique_groups))
      colors <- distinctColorPalette(k = n1, altCol = FALSE, runTsne = FALSE)
      
      # Generate the Venn diagram with custom colors
      venn.plot <- venn.diagram(
        x = sets,
        category.names = unique_groups,
        filename = NULL,
        output = TRUE,
        main.cex = 2,
        main="Venn Diagram of Gene List Sets",
        fill = colors,  # Dynamic color assignment
        alpha = 0.5,  # Transparency level
        cat.col = colors,  # Colors for category names
        cat.cex = 1.5,  # Size of category names
        cat.fontface = "bold",  # Bold category names
        cex = 1.5,  # Size of numbers in the diagram
        fontface = "bold",  # Bold numbers
        lwd = 2  # Line width of the circles
      )
      
      # To plot the Venn diagram
      grid.draw(venn.plot)
      
    }, height = 500, width = 500
    )
    
    
    # Summary pheatmap
    gene_list_summary_pheat  <- eventReactive(input$ViewGeneInfo,{
      gene_raw <- gene_raw()
      #AZPheWAS
      az_gene <- unique(c(unique(az_gene_binary_raw$Gene),unique(az_gene_continuous_raw$Gene),unique(az_var_binary_raw$Gene),unique(az_var_continuous_raw$Gene)))
      az_gene_list <- gene_raw %>% 
        filter(Gene %in% az_gene) %>% 
        mutate(AZ="AstraZeneca PheWAS Portal")
      
      if (nrow(az_gene_list)!=0) {
        az_gene_list
      } else {
        az_gene_list <- data.frame(Group=gene_raw$Group,Gene=gene_raw$Gene,AZ=NA)
      }
      
      
      #GWAS Catalog
      gwas_gene <- unlist(strsplit(gwas_raw$MAPPED_GENE, split= " - ", fixed=TRUE))
      gwas_gene <- unique(gwas_gene)
      gwas_gene_list <- gene_raw %>% 
        filter(Gene %in% gwas_gene) %>% 
        mutate(GWASCatalog="GWAS Catalog")
      
      if (nrow(gwas_gene_list)!=0) {
        gwas_gene_list
      } else {
        gwas_gene_list <- data.frame(Group=gene_raw$Group,Gene=gene_raw$Gene,GWASCatalog=NA)
      }
      
      
      
      #FinnGen 
      finngen_gene <- unlist(strsplit(finngen_raw$Gene, split= ",", fixed=TRUE))
      finngen_gene <- unique(finngen_gene)
      finngen_gene_list <- gene_raw %>% 
        filter(Gene %in% finngen_gene) %>% 
        mutate(FinnGen="FinnGen")
      
      if (nrow(finngen_gene_list)!=0) {
        finngen_gene_list
      } else {
        finngen_gene_list <- data.frame(Group=gene_raw$Group,Gene=gene_raw$Gene,FinnGen=NA)
      }
      
      
      gene_list_summary_tab <- Reduce(full_join,list(gene_raw, az_gene_list, finngen_gene_list,gwas_gene_list))
      gene_list_summary_pheat <- as.data.frame(gene_list_summary_tab)
      for (i in 3:ncol(gene_list_summary_pheat)) {
        gene_list_summary_pheat[,i][is.na(gene_list_summary_pheat[,i])==FALSE] <- "Significant Association"
        gene_list_summary_pheat[,i][is.na(gene_list_summary_pheat[,i])] <- "No Significant Association"
      }
      
      gene_list_summary_pheat <- gene_list_summary_pheat %>% 
        pivot_longer(cols = -c(Group,Gene))
      colnames(gene_list_summary_pheat) <- c("Group","Gene","Databases", "Type")
      gene_list_summary_pheat$Type <- factor(gene_list_summary_pheat$Type,
                                             levels = c("Significant Association","No Significant Association"))
      
      
      gene_list_summary_pheat <- gene_list_summary_pheat %>% 
        group_by(Gene,Databases) %>% 
        mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
        ungroup() %>% 
        group_by(Group,Gene,Databases) %>% 
        dplyr::slice(1)
      
      gene_list_summary_pheat
      
    })
    
    
    
    output$p1_plot2 <- renderPlot({
      gene_list_summary_pheat() %>% 
        ggplot(aes(x=Databases,y=Gene,fill=Type)) + 
        geom_tile(size=0.5, color="black",alpha=0.5)+
        scale_fill_manual(values = c("#89ABE3FF", "#b8b8b8")) +
        ggpubr::theme_classic2()+
        facet_grid(Group~., scales = "free",space = "free")+
        theme(strip.text.y = element_text(angle = 0))
      
      
    }, height = 1000, width = 800
    )
    
    

    #########################################Page 2 Gene-Phenotype Associations#########################################
    
    ## Summary pheatmap
      summary_gene_pheat  <- eventReactive(input$ViewGenePhenotypeAsso,{
        gene_raw <- gene_raw()
        #AstraZeneca PheWAS Portal
        az_gene_binary_flt <- az_gene_binary_flt()
        az_gene_continuous_flt <- az_gene_continuous_flt()
        az_var_binary_flt <- az_var_binary_flt()
        az_var_continuous_flt <- az_var_continuous_flt()
        
        if (nrow(az_gene_binary_flt) != 0) {
          az_gene_binary_flt_asso <- az_gene_binary_flt %>% 
            mutate(Variant=NA) %>% 
            dplyr::select(Gene,Variant,Phenotype,Category)
        } else {
          az_gene_binary_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                                Variant=NA,
                                                Phenotype=NA,
                                                Category=NA
          )
        }
        
        #az_gene_continuous_flt
        if (nrow(az_gene_continuous_flt) != 0) {
          az_gene_continuous_flt_asso <- az_gene_continuous_flt %>% 
            mutate(Variant=NA) %>% 
            dplyr::select(Gene,Variant,Phenotype,Category) 
        } else {
          az_gene_continuous_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                                    Variant=NA,
                                                    Phenotype=NA,
                                                    Category=NA
          )
        }
        
        #az_var_binary_flt
        
        if (nrow(az_var_binary_flt) != 0) {
          az_var_binary_flt_asso <- az_var_binary_flt %>% 
            dplyr::select(Gene,Variant,Phenotype,Category) 
          
        } else {
          az_var_binary_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                               Variant=NA,
                                               Phenotype=NA,
                                               Category=NA
          )
        }
        
        #az_var_continuous_flt       
        if (nrow(az_var_continuous_flt) != 0) {
          az_var_continuous_flt_asso <- az_var_continuous_flt %>% 
            dplyr::select(Gene,Variant,Phenotype,Category)
          
        } else {
          az_var_continuous_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                                   Variant=NA,
                                                   Phenotype=NA,
                                                   Category=NA
          )
        }
        
        
        
        az_gene_asso <- do.call(rbind,list(az_gene_binary_flt_asso,az_gene_continuous_flt_asso,az_var_binary_flt_asso,az_var_continuous_flt_asso))
        
        az_gene_asso <- az_gene_asso %>% 
          group_by(Gene) %>%
          summarise(Phenotype=paste0(na.omit(unique(Phenotype)),collapse = "; "),
                    Category=paste0(na.omit(unique(Category)),collapse = "; "), 
                    Variant=paste0(na.omit(unique(Variant)),collapse = "; ")
          ) %>%
          ungroup()
        
        colnames(az_gene_asso) <- c("Gene","AZPhenotype","AZCategory","AZVariant")
        
        #GWAS Catalog
        #Mapping efo id to trait category
        efo <- efo_raw %>% 
          dplyr::select(`Parent term`,`EFO URI`) %>% 
          group_by(`EFO URI`) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        gwas_raw_flt <- gwas_raw_flt()
        if (nrow(gwas_raw_flt) !=0) {
          gwas_raw_flt_asso <- gwas_raw_flt %>% 
            mutate(`EFO URI`=MAPPED_TRAIT_URI) %>% 
            separate_rows(`EFO URI`, sep = ", ",convert = TRUE) %>% 
            left_join(efo) %>% 
            group_by(Gene) %>% 
            summarise(Phenotype=paste0(na.omit(unique(MAPPED_TRAIT)),collapse = "; "),
                      Category=paste0(na.omit(unique(`Parent term`)),collapse = "; "),
                      Variant=paste0(na.omit(unique(SNPS)),collapse = "; ")
            ) %>%
            ungroup()
          colnames(gwas_raw_flt_asso) <- c("Gene","GWASCatalogPhenotype","GWASCatalogCategory","GWASCatalogVariant")
          
        } else {
          gwas_raw_flt_asso <- data.frame(
            Gene=unique(gene_raw$Gene),
            GWASCatalogPhenotype=NA,
            GWASCatalogCategory=NA,
            GWASCatalogVariant=NA
            
          )
        }
        
        
        
        #FinnGen 
        finngen_raw_flt <- finngen_raw_flt()
        if (nrow(finngen_raw_flt) != 0) {
          finngen_raw_flt_asso <- finngen_raw_flt %>% 
            dplyr::select(Gene,phenotype,category,rsids) %>% 
            group_by(Gene) %>% 
            summarise(Phenotype=paste0(na.omit(unique(phenotype)),collapse = "; "),
                      Category=paste0(na.omit(unique(category)),collapse = "; "),
                      Variant=paste0(na.omit(unique(rsids)),collapse = "; ")
            ) %>%
            ungroup()
          colnames(finngen_raw_flt_asso) <- c("Gene","FinnGenPhenotype","FinnGenCategory","FinnGenVariant")
          
        } else {
          finngen_raw_flt_asso  <- data.frame(
            Gene=unique(gene_raw$Gene),
            FinnGenPhenotype=NA,
            FinnGenCategory=NA,
            FinnGenVariant=NA
          )
        }
        
    
        
        #overlapping gene
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          summarise(Group=paste0(na.omit(unique(Group)),collapse = "; ")
          ) %>%
          ungroup() %>% 
          dplyr::select(Group,Gene)
        
        
        summary_gene <- gene_change %>% 
          left_join(az_gene_asso) %>% 
          left_join(gwas_raw_flt_asso) %>% 
          left_join(finngen_raw_flt_asso)
        summary_gene[summary_gene==''] <- NA
        
        summary_gene_pheat <- summary_gene %>% 
          dplyr::select(Group,Gene,AZCategory,GWASCatalogCategory,FinnGenCategory)
        summary_gene_pheat <- summary_gene_pheat %>% 
          pivot_longer(cols = -c(Group,Gene))
        colnames(summary_gene_pheat) <- c("Group","Gene","Databases", "Category")
        summary_gene_pheat <- summary_gene_pheat %>% 
          mutate(Type = ifelse(is.na(Category)==FALSE, "Significant Association","No Significant Association"))
        
        summary_gene_pheat$Type <- factor(summary_gene_pheat$Type,
                                          levels = c("Significant Association","No Significant Association"))
        
        summary_gene_pheat
      
      })
    
    
      output$p2_plot1 = renderPlotly({
        
        p <- summary_gene_pheat() %>% 
          ggplot(aes(x=Databases,y=Gene,fill=Type,
                     label=Category
          )) + 
          geom_tile(size=0.5, color="black",alpha=0.5)+
          scale_fill_manual(values = c("#89ABE3FF", "#b8b8b8")) +
          ggpubr::theme_classic2()+
          facet_grid(Group~Databases, scales = "free",space = "free")+
          theme(strip.text.y = element_text(angle = 0))
        
        ggplotly(p,tooltip=c("Category")) %>%
          layout(autosize = T, width = 1200, height = 1500)
        
        
      })
    
    
    ## Summary of Genes with Significant Associations Across Different Databases
      summary_gene_tab  <- eventReactive(input$ViewGenePhenotypeAsso,{
        gene_raw <- gene_raw()
        #AstraZeneca PheWAS Portal
        az_gene_binary_flt <- az_gene_binary_flt()
        az_gene_continuous_flt <- az_gene_continuous_flt()
        az_var_binary_flt <- az_var_binary_flt()
        az_var_continuous_flt <- az_var_continuous_flt()
        
        if (nrow(az_gene_binary_flt) != 0) {
          az_gene_binary_flt_asso <- az_gene_binary_flt %>% 
            mutate(Variant=NA) %>% 
            dplyr::select(Gene,Variant,Phenotype,Category)
        } else {
          az_gene_binary_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                                Variant=NA,
                                                Phenotype=NA,
                                                Category=NA
          )
        }
        
        #az_gene_continuous_flt
        if (nrow(az_gene_continuous_flt) != 0) {
          az_gene_continuous_flt_asso <- az_gene_continuous_flt %>% 
            mutate(Variant=NA) %>% 
            dplyr::select(Gene,Variant,Phenotype,Category) 
        } else {
          az_gene_continuous_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                                    Variant=NA,
                                                    Phenotype=NA,
                                                    Category=NA
          )
        }
        
        #az_var_binary_flt
        
        if (nrow(az_var_binary_flt) != 0) {
          az_var_binary_flt_asso <- az_var_binary_flt %>% 
            dplyr::select(Gene,Variant,Phenotype,Category) 
          
        } else {
          az_var_binary_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                               Variant=NA,
                                               Phenotype=NA,
                                               Category=NA
          )
        }
        
        #az_var_continuous_flt       
        if (nrow(az_var_continuous_flt) != 0) {
          az_var_continuous_flt_asso <- az_var_continuous_flt %>% 
            dplyr::select(Gene,Variant,Phenotype,Category)
          
        } else {
          az_var_continuous_flt_asso <- data.frame(Gene=unique(gene_raw$Gene),
                                                   Variant=NA,
                                                   Phenotype=NA,
                                                   Category=NA
          )
        }
        
        
        
        az_gene_asso <- do.call(rbind,list(az_gene_binary_flt_asso,az_gene_continuous_flt_asso,az_var_binary_flt_asso,az_var_continuous_flt_asso))
        
        az_gene_asso <- az_gene_asso %>% 
          group_by(Gene) %>%
          summarise(Phenotype=paste0(na.omit(unique(Phenotype)),collapse = "; "),
                    Category=paste0(na.omit(unique(Category)),collapse = "; "), 
                    Variant=paste0(na.omit(unique(Variant)),collapse = "; ")
          ) %>%
          ungroup()
        
        colnames(az_gene_asso) <- c("Gene","AZPhenotype","AZCategory","AZVariant")
        
        #GWAS Catalog
        #Mapping efo id to trait category
        efo <- efo_raw %>% 
          dplyr::select(`Parent term`,`EFO URI`) %>% 
          group_by(`EFO URI`) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        gwas_raw_flt <- gwas_raw_flt()
        if (nrow(gwas_raw_flt) !=0) {
          gwas_raw_flt_asso <- gwas_raw_flt %>% 
            mutate(`EFO URI`=MAPPED_TRAIT_URI) %>% 
            separate_rows(`EFO URI`, sep = ", ",convert = TRUE) %>% 
            left_join(efo) %>% 
            group_by(Gene) %>% 
            summarise(Phenotype=paste0(na.omit(unique(MAPPED_TRAIT)),collapse = "; "),
                      Category=paste0(na.omit(unique(`Parent term`)),collapse = "; "),
                      Variant=paste0(na.omit(unique(SNPS)),collapse = "; ")
            ) %>%
            ungroup()
          colnames(gwas_raw_flt_asso) <- c("Gene","GWASCatalogPhenotype","GWASCatalogCategory","GWASCatalogVariant")
          
        } else {
          gwas_raw_flt_asso <- data.frame(
            Gene=unique(gene_raw$Gene),
            GWASCatalogPhenotype=NA,
            GWASCatalogCategory=NA,
            GWASCatalogVariant=NA
            
          )
        }
        
        
        
        #FinnGen 
        finngen_raw_flt <- finngen_raw_flt()
        if (nrow(finngen_raw_flt) != 0) {
          finngen_raw_flt_asso <- finngen_raw_flt %>% 
            dplyr::select(Gene,phenotype,category,rsids) %>% 
            group_by(Gene) %>% 
            summarise(Phenotype=paste0(na.omit(unique(phenotype)),collapse = "; "),
                      Category=paste0(na.omit(unique(category)),collapse = "; "),
                      Variant=paste0(na.omit(unique(rsids)),collapse = "; ")
            ) %>%
            ungroup()
          colnames(finngen_raw_flt_asso) <- c("Gene","FinnGenPhenotype","FinnGenCategory","FinnGenVariant")
          
        } else {
          finngen_raw_flt_asso  <- data.frame(
            Gene=unique(gene_raw$Gene),
            FinnGenPhenotype=NA,
            FinnGenCategory=NA,
            FinnGenVariant=NA
          )
        }
        
        
        
        #overlapping gene
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          summarise(Group=paste0(na.omit(unique(Group)),collapse = "; ")
          ) %>%
          ungroup() %>% 
          dplyr::select(Group,Gene)
        
        
        summary_gene <- gene_change %>% 
          left_join(az_gene_asso) %>% 
          left_join(gwas_raw_flt_asso) %>% 
          left_join(finngen_raw_flt_asso)
        summary_gene[summary_gene==''] <- NA
        summary_gene
        
        
        
      })
      
   
      
      output$p2_dt1 <- DT::renderDataTable(
        
        summary_gene_tab() %>%
          DT::datatable(escape = FALSE,
                        rownames = FALSE,
                        plugins = "ellipsis",
                        extensions = 'Buttons',

                        options = list(
                          dom = 'Bfrtip',
                          buttons = list(
                            list(extend = "csv", text = "Download Full Results", filename = "Full_data",
                                 exportOptions = list(
                                   modifier = list(page = "all"),
                                   orthogonal = "export"
                                 )
                            )
                          ),
                          # paging = TRUE,    # Enable pagination
                          # searching = TRUE, # Enable search functionality
                          # ordering = TRUE,
                          #
                          #buttons = c('copy', 'csv'),
                          columnDefs = list(list(
                            targets = c(2:10),
                            render = JS("$.fn.dataTable.render.ellipsis(17, false )")
                          ))


                        ))

        
      )
      
    
    
      
      #########################################Page 3 Variant-Phenotype Associations#########################################
      #Summary Table
      summary_var <- eventReactive(input$ViewVarTable,{
        gene_raw <- gene_raw()
        #AstraZeneca PheWAS Portal
        az_var_binary_flt_asso <- az_var_binary_flt() %>% 
          dplyr::select(Gene,Variant,Phenotype,Category,rsID,AF, ClinVar) 
        az_var_continuous_flt_asso <- az_var_continuous_flt() %>% 
          dplyr::select(Gene,Variant,Phenotype,Category,rsID,AF, ClinVar)
        az_var_asso <- do.call(rbind,list(az_var_binary_flt_asso,az_var_continuous_flt_asso))
        az_var_asso <- az_var_asso %>% 
          mutate(gnomADLink=paste0("<a target = '_blank' href= ", "https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4",">","https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4","</a>"))
        
        az_var_asso <- az_var_asso %>% 
          group_by(Gene,Variant,rsID,AF,gnomADLink,ClinVar) %>%
          summarise(
            Category=paste0(na.omit(unique(Category)),collapse = "; "),
            Phenotype=paste0(na.omit(unique(Phenotype)),collapse = "; ")
          ) %>%
          ungroup()
        colnames(az_var_asso) <- c("Gene","Variant","rsID","AF","gnomADLink","ClinVar","AZCategory","AZPhenotype")
        
        #GWAS Catalog
        #Mapping efo id to trait category
        efo <- efo_raw %>% 
          dplyr::select(`Parent term`,`EFO URI`) %>% 
          group_by(`EFO URI`) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        gwas_raw_flt <- gwas_raw_flt()
        if (nrow(gwas_raw_flt) > 0) {
          gwas_raw_flt_asso <- gwas_raw_flt %>% 
            mutate(`EFO URI`=MAPPED_TRAIT_URI) %>% 
            separate_rows(`EFO URI`, sep = ", ",convert = TRUE) %>% 
            left_join(efo) %>% 
            mutate(Label=row_number())
          
          gnomad <- gnomad_raw
          colnames(gnomad) <- c("Variant","rsID","AF")
          gnomad <- gnomad %>% 
            #mutate(gnomADLink=paste0("https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4")) %>% 
            mutate(gnomADLink=paste0("<a target = '_blank' href= ", "https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4",">","https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4","</a>"))
          
          gwas_raw_flt_asso <- gwas_raw_flt_asso %>% 
            left_join(gnomad,by=c("SNPS"="rsID"))
          clinvar <- clinvar_raw
          clinvar <- clinvar %>% 
            mutate(Variant=paste0(V1,"-",V2,"-",V4,"-",V5)) %>% 
            dplyr::select(Variant,V6)
          colnames(clinvar) <- c("Variant","ClinVar")
          gwas_raw_flt_asso <- gwas_raw_flt_asso %>% 
            left_join(clinvar)
          
          gwas_raw_flt_asso <- gwas_raw_flt_asso %>% 
            group_by(Label) %>% 
            mutate(Variant=paste0(Variant,collapse = "; "),
                   AF=paste0(AF,collapse = "; "),
                   gnomADLink=paste0(gnomADLink,collapse = "; "),
                   ClinVar=paste0(ClinVar,collapse = "; ")
            ) %>% 
            dplyr::slice(1) %>% 
            ungroup()
          
          
          gwas_raw_flt_asso <- gwas_raw_flt_asso %>% 
            group_by(Gene,Variant,SNPS,AF,gnomADLink,ClinVar) %>%
            summarise(
              Category=paste0(na.omit(unique(`Parent term`)),collapse = "; "),
              Phenotype=paste0(na.omit(unique(MAPPED_TRAIT)),collapse = "; ")
            ) %>%
            ungroup() 
          colnames(gwas_raw_flt_asso) <- c("Gene","GWASCatalogVariant","rsID","GWASCatalogAF","GWASCataloggnomADLink","GWASCatalogClinVar","GWASCatalogCategory","GWASCatalogPhenotype")
          
        } else {
          gwas_raw_flt_asso <- data.frame(
            Gene=NA,
            GWASCatalogVariant=NA,
            rsID=NA,
            GWASCatalogAF=NA,
            GWASCataloggnomADLink=NA,
            GWASCatalogClinVar=NA,
            GWASCatalogCategory=NA,
            GWASCatalogPhenotype=NA
          ) %>% 
            na.omit() %>% 
            as.tibble()
          gwas_raw_flt_asso$Gene <- as.character(gwas_raw_flt_asso$Gene)
          gwas_raw_flt_asso$rsID <- as.character(gwas_raw_flt_asso$rsID)
        }
        
        
        
        
        
        #FinnGen 
        finngen_raw_flt_asso <- finngen_raw_flt() %>% 
          mutate(gnomADLink=paste0("<a target = '_blank' href= ", "https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4",">","https://gnomad.broadinstitute.org/variant/",Variant,"?dataset=gnomad_r4","</a>")) %>% 
          dplyr::select(Gene,Variant,rsID,AF,gnomADLink, ClinVar,phenotype,category) %>% 
          group_by(Gene,Variant,rsID,AF,gnomADLink, ClinVar) %>% 
          summarise(
            Category=paste0(na.omit(unique(category)),collapse = "; "),
            Phenotype=paste0(na.omit(unique(phenotype)),collapse = "; ")
          ) %>%
          ungroup()
        colnames(finngen_raw_flt_asso) <- c("Gene","Variant","rsID","AF","gnomADLink","ClinVar","FinnGenCategory","FinnGenPhenotype")
        
        summary_var <- az_var_asso %>% 
          full_join(gwas_raw_flt_asso) %>% 
          full_join(finngen_raw_flt_asso) 
        
        
        summary_var$Variant[is.na(summary_var$Variant)] <- summary_var$GWASCatalogVariant[is.na(summary_var$Variant)]
        gene_change <- gene_raw() %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        summary_var <- summary_var %>% 
          left_join(gene_change) %>% 
          dplyr::select(Group,everything())
        
        summary_var
        
      })
      
      
      # Summary Table Output
      output$p3_dt1 <- DT::renderDataTable(

        summary_var() %>% 
          DT::datatable(escape = FALSE,
                        rownames = FALSE,
                        plugins = "ellipsis",
                        extensions = 'Buttons',
                        options = list(
                          
                          buttons = list(
                            list(extend = "csv", text = "Download Full Results", filename = "Full_data",
                                 exportOptions = list(
                                   modifier = list(page = "all"),
                                   orthogonal = "export"
                                 )
                            )
                          ),
                          dom = 'Bfrtip',
                          buttons = c('copy', 'csv'),
                          columnDefs = list(list(
                            targets = c(7:10,12:16),
                            render = JS("$.fn.dataTable.render.ellipsis(17, false )")
                          ))
                        )) 

        
        
      )
      
      #AstraZeneca PheWAS Portal
      
      # Phenotype-GeneSet Cluster
      
      az_sum_pheat <- eventReactive(input$ViewAZCluster,{
        gene_raw <- gene_raw()
        #AstraZeneca PheWAS Portal
        az_gene_binary_flt_sum <- az_gene_binary_flt() %>% 
          dplyr::select(Gene,Category,Phenotype)
        az_gene_continuous_flt_sum <- az_gene_continuous_flt() %>% 
          dplyr::select(Gene,Category,Phenotype)
        az_var_binary_flt_sum  <- az_var_binary_flt() %>% 
          dplyr::select(Gene,Category,Phenotype)
        az_var_continuous_flt_sum <- az_var_continuous_flt() %>% 
          dplyr::select(Gene,Category,Phenotype)
        
        az_sum <- do.call(rbind,list(az_gene_binary_flt_sum,az_gene_continuous_flt_sum,az_var_binary_flt_sum,az_var_continuous_flt_sum))
        
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        az_sum <-az_sum %>% 
          left_join(gene_change) %>% 
          group_by(Group,Gene,Category) %>% 
          summarise(
            n=length(unique(Phenotype)))  %>% 
          ungroup()  %>% 
          mutate(n=ifelse(n>1,1,n))
        
        #Summary table and figure
        multi_gene <- az_sum %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          group_by(Gene) %>% 
          mutate(n=n()) %>% 
          filter(n>1) %>% 
          ungroup()
        
        
        az_sum_pheat <- az_sum %>% 
          pivot_wider(
            names_from = "Category",
            values_from = "n",
            values_fill = 0
          )
        
        az_sum_pheat <- az_sum_pheat %>% 
          mutate(Gene=ifelse((Gene %in% multi_gene$Gene),paste0(Gene,"_",Group),Gene)) %>% 
          as.data.frame()
        rownames(az_sum_pheat) <- az_sum_pheat$Gene
        az_sum_pheat
      })
      
      
      output$p3_plot1 <- renderPlot({
        az_sum_pheat <- az_sum_pheat()
        if (nrow(az_sum_pheat)>0) {
          az_sum_pheat <- az_sum_pheat
          ann <- data.frame(az_sum_pheat$Group)
          rownames(ann) <- az_sum_pheat$Gene
          colnames(ann) <- "Group"
          az_sum_pheat <- az_sum_pheat[,-c(1,2)]
          
          set.seed(1111)
          n1 <- length(unique(ann$Group))
          col1=distinctColorPalette(k = n1, altCol = FALSE, runTsne = FALSE)
          names(col1) <- unique(ann$Group)
          
          ann_colors = list(
            Group = col1
          )
          
          ComplexHeatmap::pheatmap(t(az_sum_pheat),
                                   angle_col = "45",
                                   name=str_wrap("Significant association within the gene",20),
                                   row_labels = str_wrap(rownames(t(az_sum_pheat)),30),
                                   annotation_col = ann,
                                   color = c("white", "red"),
                                   legend = F,
                                   main = "AstraZeneca PheWAS Portal",
                                   annotation_colors = ann_colors)
        } else {
          plot_exception("No significant association within the gene list")
        }
        
        
        
      }, height = 800, width = 1000
      )
      
      
      
      
      # Variant-level
      
     
      az_var_flt_pheno <- eventReactive(input$ViewAZVar,{
        gene_raw <- gene_raw()
        #AstraZeneca PheWAS Portal
        # az_var_binary_flt
        # az_var_continuous_flt
        # Number of unique phenotype within each gene list
        az_var_binary_flt_pheno <- az_var_binary_flt() %>% 
          dplyr::select(Category,Phenotype,Gene,rsID)
        az_var_continuous_flt_pheno <- az_var_continuous_flt() %>% 
          dplyr::select(Category,Phenotype,Gene,rsID)
        az_var_flt_pheno <- do.call(rbind,list(az_var_binary_flt_pheno,az_var_continuous_flt_pheno))
        
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        az_var_flt_pheno <- az_var_flt_pheno %>% 
          left_join(gene_change) %>% 
          group_by(Group,Category) %>% 
          summarise(n_pheno = length(unique(Phenotype)),
                    Gene=paste0(unique(Gene),collapse = "; "),
                    Phenotype=paste0(unique(Phenotype),collapse = "; "),
                    rsID = paste0(unique(rsID),collapse = "; ")
          ) %>% 
          ungroup()
        
        az_var_flt_pheno
      })
      
      
      az_var_flt_pheno_pro  <- eventReactive(input$ViewAZVar,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        gene_number <- gene_change %>% 
          group_by(Group) %>% 
          mutate(GeneNumber=length(unique(Gene))) %>% 
          ungroup()
        az_var_binary_flt_pheno_pro <- az_var_binary_flt() %>% 
          dplyr::select(Category,Phenotype,Gene,rsID)
        az_var_continuous_flt_pheno_pro <- az_var_continuous_flt() %>% 
          dplyr::select(Category,Phenotype,Gene,rsID)
        az_var_flt_pheno_pro <- do.call(rbind,list(az_var_binary_flt_pheno_pro,az_var_continuous_flt_pheno_pro))
        
        az_var_flt_pheno_pro <- az_var_flt_pheno_pro %>% 
          left_join(gene_number) %>% 
          group_by(Group,GeneNumber,Category) %>% 
          summarise(n_gene = length(unique(Gene)),
                    Gene=paste0(unique(Gene),collapse = "; "),
                    Phenotype=paste0(unique(Phenotype),collapse = "; "),
                    rsID = paste0(unique(rsID),collapse = "; ")
          ) %>% 
          ungroup() %>% 
          mutate(GeneFraction=n_gene/GeneNumber*100)
        az_var_flt_pheno_pro
      })
      
      
      output$p3_plot2 = renderPlotly({
        az_var_flt_pheno <- az_var_flt_pheno()
        if (nrow(az_var_flt_pheno)>0) {
          p <- az_var_flt_pheno %>% 
            ggplot(aes(x=n_pheno,y=Category,
                       fill=Group,
                       label=Gene,label2=rsID,label3=Phenotype
            ))+
            geom_col()+
            facet_grid(~Group) +
            theme(
              legend.position = "none"
            )+
            # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            #    theme(strip.text.y = element_text(size = 1),
            #          legend.position = "none"
            #          ) +
            labs(x="The unique phenotype number",y="Phenotype Category")+
            theme_bw()
          
          #ggplotly format
          ggplotly(p,tooltip=c("Gene", "rsID","Phenotype")) %>%
            layout(autosize = T, width = 1200, height = 600)
        } else {
          plotly_exception("No significant association within the gene list")

        }
        
        
        
      })
      
      
      output$p3_plot3 = renderPlotly({
        az_var_flt_pheno_pro <- az_var_flt_pheno_pro()
        if (nrow(az_var_flt_pheno_pro)>0) {
          p <- az_var_flt_pheno_pro %>% 
            ggplot(aes(x=GeneFraction,y=Category,
                       fill=Group,
                       label=Gene,label2=rsID,label3=Phenotype
            ))+
            geom_col()+
            facet_grid(~Group) +
            theme(
              legend.position = "none"
            )+
            # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            #    theme(strip.text.y = element_text(size = 1),
            #          legend.position = "none"
            #          ) +
            labs(x="Percentage of genes in each gene set group",y="Phenotype Category") +
            theme_bw()
          
          #ggplotly format
          ggplotly(p,tooltip=c("Gene", "rsID","Phenotype")) %>%
            layout(autosize = T, width = 1200, height = 600)
        } else {
          plotly_exception("No significant association within the gene list")
        }
        

      })
      
      
      
      
      # Variant-Phenotype Gene Effect
      
      az_var_flt_pheno_effect <- eventReactive(input$ViewAZEffect,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        #AstraZeneca PheWAS Portal
        az_var_binary_flt_pheno_effect <- az_var_binary_flt() %>% 
          dplyr::select(Gene,Variant,Phenotype,Category,`Odds ratio`,Model) %>% 
          left_join(gene_change) %>% 
          group_by(Group,Gene,Category) %>% 
          mutate(MeanEffect=mean(`Odds ratio`)) %>% 
          ungroup() %>% 
          group_by(Group,Gene,Category,MeanEffect) %>% 
          summarise(
            Phenotype=paste0(unique(Phenotype),collapse = "; "),
            Variant=paste0(unique(Phenotype),collapse = "; "),
            Model=paste0(unique(Model),collapse = "; ")
          ) %>% 
          ungroup() %>% 
          mutate(Trait="Binary")

        az_var_continuous_flt_pheno_effect <- az_var_continuous_flt() %>% 
          dplyr::select(Gene,Variant,Phenotype,Category,`Effect size`,Model) %>% 
          left_join(gene_change) %>% 
          group_by(Group,Gene,Category) %>% 
          mutate(MeanEffect=mean(abs(`Effect size`))) %>% 
          ungroup() %>% 
          group_by(Group,Gene,Category,MeanEffect) %>% 
          summarise(
            Phenotype=paste0(unique(Phenotype),collapse = "; "),
            Variant=paste0(unique(Phenotype),collapse = "; "),
            Model=paste0(unique(Model),collapse = "; ")
          ) %>% 
          ungroup() %>% 
          mutate(Trait="Continuous")
        
        az_var_flt_pheno_effect <- do.call(rbind,list(az_var_binary_flt_pheno_effect,az_var_continuous_flt_pheno_effect))
        
        az_var_flt_pheno_effect
        
      
      })
      
      
      output$p3_plot_add1 <- renderPlot({
        az_var_flt_pheno_effect <- az_var_flt_pheno_effect()
        if (nrow(az_var_flt_pheno_effect)>0) {
          az_var_flt_pheno_effect <- az_var_flt_pheno_effect
          az_var_flt_pheno_effect_size <- az_var_flt_pheno_effect %>% 
            group_by(Category) %>% 
            summarise(n=n()) %>% 
            ungroup() %>% 
            mutate(size=(n/10+0.4))
          
          az_var_flt_pheno_effect %>%
            mutate(MeanEffect=ifelse(MeanEffect>10,10,MeanEffect)) %>% 
            #filter(Trait!="Binary") %>% 
            #filter(MeanEffect<10) %>% 
            ggplot(aes(x=log10(MeanEffect),y=Category,
                       fill=Group,
                       color=Group,
                       label=Gene
            ))+
            geom_jitter(position = position_dodge2(0.9))+
            #geom_col(position = position_dodge2(preserve = "single", padding = 0.1))+
            facet_grid(Category~Trait,scales = "free",space = "free") +
            theme_bw()+
            theme(
              #legend.position = "none",
              strip.text.y = element_blank()
            )+
            # # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            # #    theme(strip.text.y = element_text(size = 1),
            # #          legend.position = "none"
            # #          ) +
            labs(x="Log10(MeanEffect) (Odds ratio or Effect size)",y="Phenotype Category")+
            ggrepel::geom_text_repel(position = position_dodge2(0.9))+
            ggh4x::force_panelsizes(rows = az_var_flt_pheno_effect_size$size)
        } else {
          plot_exception("No significant association within the gene list")
        }
        

      }, height = 850, width = 1500
      )
      
      #options(digits=3)
      
      # Summary Table Output
      output$p3_dt2 <- DT::renderDataTable(
        
        az_var_flt_pheno_effect() %>% 
          dplyr::select(Group,Gene,Category,Phenotype,Variant,Model,Trait,MeanEffect) %>% 
          mutate_if(is.numeric, round, digits = 3)  %>% 
        DT::datatable(escape = FALSE,
                      rownames = FALSE,
                      plugins = "ellipsis",
                      extensions = 'Buttons',
                      options = list(
                        
                        buttons = list(
                          list(extend = "csv", text = "Download Full Results", filename = "Full_data",
                               exportOptions = list(
                                 modifier = list(page = "all"),
                                 orthogonal = "export"
                               )
                          )
                        ),
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv'),
                        columnDefs = list(list(
                          targets = c(3:4),
                          render = JS("$.fn.dataTable.render.ellipsis(17, false )")
                        ))
                      )) 
        
        
        
        
      )
      
      
      
      
      
      
      
      
      
      #GWAS Catalog
     #Phenotype cluster
      gwas_raw_flt_gene_sum_pheat <- eventReactive(input$ViewGWASCatalogCluster,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        gene_number <- gene_change %>%
          group_by(Group) %>%
          mutate(GeneNumber=length(unique(Gene))) %>%
          ungroup()
        efo <- efo_raw %>%
          dplyr::select(`Parent term`,`EFO URI`) %>%
          group_by(`EFO URI`) %>%
          dplyr::slice(1) %>%
          ungroup()

        #Each gene each trait only select one
        gwas_raw_flt <- gwas_raw_flt()
        if (nrow(gwas_raw_flt)>0) {
          gwas_raw_flt_gene_sum <- gwas_raw_flt %>%
            left_join(gene_number) %>%
            mutate(`EFO URI`=MAPPED_TRAIT_URI) %>%
            separate_rows(`EFO URI`, sep = ", ",convert = TRUE) %>%
            left_join(efo) %>%
            dplyr::select(Group,`Parent term`,MAPPED_TRAIT,Gene,SNPS) %>%
            group_by(Group,Gene,`Parent term`) %>%
            summarise(
              n=length(unique(SNPS)))  %>%
            ungroup()  %>%
            mutate(n=ifelse(n>1,1,n))
          
          #Summary table and figure
          multi_gene <- gwas_raw_flt_gene_sum %>%
            group_by(Group,Gene) %>%
            dplyr::slice(1) %>%
            ungroup() %>%
            group_by(Gene) %>%
            mutate(n=n()) %>%
            filter(n>1) %>%
            ungroup()
          
          gwas_raw_flt_gene_sum_pheat <- gwas_raw_flt_gene_sum %>%
            pivot_wider(
              names_from = `Parent term`,
              values_from = "n",
              values_fill = 0
            )
          
          gwas_raw_flt_gene_sum_pheat <- gwas_raw_flt_gene_sum_pheat %>%
            mutate(Gene=ifelse((Gene %in% multi_gene$Gene),paste0(Gene,"_",Group),Gene)) %>%
            as.data.frame()
          rownames(gwas_raw_flt_gene_sum_pheat) <- gwas_raw_flt_gene_sum_pheat$Gene
          
          gwas_raw_flt_gene_sum_pheat
        } else {
          gwas_raw_flt_gene_sum_pheat <- NULL
        }
        
      })
      


      output$p3_plot6 <- renderPlot({
        gwas_raw_flt_gene_sum_pheat <- gwas_raw_flt_gene_sum_pheat()
        if (is.null(gwas_raw_flt_gene_sum_pheat)==FALSE) {
          gwas_raw_flt_gene_sum_pheat <- gwas_raw_flt_gene_sum_pheat
          ann <- data.frame(gwas_raw_flt_gene_sum_pheat$Group)
          rownames(ann) <- gwas_raw_flt_gene_sum_pheat$Gene
          colnames(ann) <- "Group"
          gwas_raw_flt_gene_sum_pheat <- gwas_raw_flt_gene_sum_pheat[,-c(1,2)]
          
          
          set.seed(1111)
          n1 <- length(unique(ann$Group))
          col1=distinctColorPalette(k = n1, altCol = FALSE, runTsne = FALSE)
          names(col1) <- unique(ann$Group)
          
          ann_colors = list(
            Group = col1
          )
          
          ComplexHeatmap::pheatmap(t(gwas_raw_flt_gene_sum_pheat),
                                   angle_col = "45",
                                   name=str_wrap("Significant association within the gene",20),
                                   row_labels = str_wrap(rownames(t(gwas_raw_flt_gene_sum_pheat)),30),
                                   annotation_col = ann,
                                   color = c("white", "red"),
                                   legend = F,
                                   main = "GWAS Catalog",
                                   annotation_colors = ann_colors)
        } else {
          plot_exception("No significant association within the gene list")
        }
        


      }, height = 800, width = 1200
      ) 


      
      
      
      #gwas_raw_flt_gene
      
      gwas_raw_flt_gene <- eventReactive(input$ViewGWASCatalogPheno,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        gene_number <- gene_change %>% 
          group_by(Group) %>% 
          mutate(GeneNumber=length(unique(Gene))) %>% 
          ungroup()
        efo <- efo_raw %>% 
          dplyr::select(`Parent term`,`EFO URI`) %>% 
          group_by(`EFO URI`) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        #Each gene each trait only select one
        if (nrow(gwas_raw_flt())>0) {
          gwas_raw_flt_gene <- gwas_raw_flt() %>% 
            left_join(gene_number) %>% 
            mutate(`EFO URI`=MAPPED_TRAIT_URI) %>% 
            separate_rows(`EFO URI`, sep = ", ",convert = TRUE) %>% 
            left_join(efo) %>% 
            dplyr::select(Group,`Parent term`,MAPPED_TRAIT,Gene,SNPS) %>% 
            group_by(Group,`Parent term`) %>% 
            mutate(
              n=length(unique(MAPPED_TRAIT)),
              Phenotype=paste0(unique(MAPPED_TRAIT),collapse = "; "),
              rsids = paste0(unique(SNPS),collapse = "; "),
              Gene = paste0(unique(Gene),collapse = "; ")) %>% 
            dplyr::slice(1) %>% 
            ungroup()
          
          
          gwas_raw_flt_gene
        } else {
          gwas_raw_flt_gene <- NULL
        }
        
      })
      
    
      
      #gwas_raw_flt_gene_pro
      gwas_raw_flt_gene_pro  <- eventReactive(input$ViewGWASCatalogPheno,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        gene_number <- gene_change %>% 
          group_by(Group) %>% 
          mutate(GeneNumber=length(unique(Gene))) %>% 
          ungroup()
        
        efo <- efo_raw %>% 
          dplyr::select(`Parent term`,`EFO URI`) %>% 
          group_by(`EFO URI`) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        #Each gene each trait only select one
        if (nrow(gwas_raw_flt())>0) {
          gwas_raw_flt_gene_pro <- gwas_raw_flt() %>% 
            left_join(gene_number) %>% 
            mutate(`EFO URI`=MAPPED_TRAIT_URI) %>% 
            separate_rows(`EFO URI`, sep = ", ",convert = TRUE) %>% 
            left_join(efo) %>% 
            dplyr::select(Group,GeneNumber,`Parent term`,MAPPED_TRAIT,Gene,SNPS) %>% 
            group_by(Group,GeneNumber,`Parent term`) %>% 
            mutate(
              n=length(unique(Gene)),
              Phenotype=paste0(unique(MAPPED_TRAIT),collapse = "; "),
              rsids = paste0(unique(SNPS),collapse = "; "),
              Gene = paste0(unique(Gene),collapse = "; ")) %>% 
            dplyr::slice(1) %>% 
            ungroup()  %>% 
            mutate(GeneFraction=n/GeneNumber*100)
          
          gwas_raw_flt_gene_pro
        } else {
          gwas_raw_flt_gene_pro <- NULL
        }
        
      })
      
      
      output$p3_plot7 = renderPlotly({
        gwas_raw_flt_gene <- gwas_raw_flt_gene()
        if (is.null(gwas_raw_flt_gene)==FALSE) {
          p <- gwas_raw_flt_gene %>% 
            ggplot(aes(x=n,y=`Parent term`,
                       #color=MultiGene,
                       fill=Group,
                       label=Gene,label2=SNPS,label3=Phenotype
            ))+
            geom_col()+
            facet_grid(~Group,scales = "free_y",space = "free") +
            theme(
              legend.position = "none"
            )+
            # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            #    theme(strip.text.y = element_text(size = 1),
            #          legend.position = "none"
            #          ) +
            labs(x="The unique phenotype number",y="MAPPED_TRAIT") +
            theme_bw()
          
          #ggplotly format
          ggplotly(p,tooltip=c("Gene", "SNPS","Phenotype")) %>%
            layout(autosize = T, width = 1200, height = 600)
        } else {
          plotly_exception("No significant association within the gene list")
        }
        
        
        
      })
      
      
      output$p3_plot8 = renderPlotly({
        
        gwas_raw_flt_gene_pro <- gwas_raw_flt_gene_pro()
        if (is.null(gwas_raw_flt_gene_pro)==FALSE) {
          p <- gwas_raw_flt_gene_pro %>% 
            ggplot(aes(x=GeneFraction,y=`Parent term`,
                       #color=MultiGene,
                       fill=Group,
                       label=Gene,label2=SNPS,label3=Phenotype
            ))+
            geom_col()+
            facet_grid(~Group,scales = "free_y",space = "free") +
            theme(
              legend.position = "none"
            )+
            # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            #    theme(strip.text.y = element_text(size = 1),
            #          legend.position = "none"
            #          ) +
            labs(x="Percentage of genes in each gene set group",y="MAPPED_TRAIT") +
            theme_bw()
          
          #ggplotly format
          ggplotly(p,tooltip=c("Gene", "SNPS","Phenotype")) %>%
            layout(autosize = T, width = 1200, height = 600)
        } else {
          plotly_exception("No significant association within the gene list")
        }
        
        
        
        
      })
      
      
      #FinnGen
      
      
      #Phenotype cluster
      #Cluster ICD-classification
      finngen_raw_flt_gene_sum_icd_pheat <- eventReactive(input$ViewFinnGenCluster,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        finngen_raw_flt_gene <- finngen_raw_flt() %>%
          left_join(gene_change)
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>%
          filter(category != "Quantitative endpoints")
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>%
          mutate(Endpoint=ifelse(category %in%
                                   c("I Certain infectious and parasitic diseases (AB1_)",
                                     "XI Diseases of the digestive system (K11_)",
                                     "IV Endocrine, nutritional and metabolic diseases (E4_)",
                                     "II Neoplasms, from cancer register (ICD-O-3)",
                                     "II Neoplasms from hospital discharges (CD2_)",
                                     "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
                                     "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                                     "VII Diseases of the eye and adnexa (H7_)",
                                     "V Mental and behavioural disorders (F5_)",
                                     "IX Diseases of the circulatory system (I9_)",
                                     "VI Diseases of the nervous system (G6_)",
                                     "VIII Diseases of the ear and mastoid process (H8_)",
                                     "X Diseases of the respiratory system (J10_)",
                                     "XII Diseases of the skin and subcutaneous tissue (L12_)",
                                     "XIV Diseases of the genitourinary system (N14_)",
                                     "XV Pregnancy, childbirth and the puerperium (O15_)",
                                     "XVI Certain conditions originating in the perinatal period (P16_)",
                                     "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)",
                                     "XXI Factors influencing health status and contact with health services (Z21_)"
                                     
                                   ), "ICD-classification", "Specific requests")
                 
          )
        finngen_raw_flt_gene_sum_icd <- finngen_raw_flt_gene %>% 
          filter(Endpoint=="ICD-classification") %>% 
          group_by(Group,Gene,category) %>% 
          summarise(n=length(unique(Variant))) %>% 
          ungroup() %>% 
          mutate(n=ifelse(n>1,1,n))
        
        #Summary table and figure
        finngen_raw_flt_gene_sum_icd <- finngen_raw_flt_gene_sum_icd
        multi_gene <- finngen_raw_flt_gene_sum_icd %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          group_by(Gene) %>% 
          mutate(n=n()) %>% 
          filter(n>1) %>% 
          ungroup()
        
        
        finngen_raw_flt_gene_sum_icd_pheat <- finngen_raw_flt_gene_sum_icd %>% 
          pivot_wider(
            names_from = "category",
            values_from = "n",
            values_fill = 0
          )
        
        finngen_raw_flt_gene_sum_icd_pheat <- finngen_raw_flt_gene_sum_icd_pheat %>% 
          mutate(Gene=ifelse((Gene %in% multi_gene$Gene),paste0(Gene,"_",Group),Gene)) %>% 
          as.data.frame()
        rownames(finngen_raw_flt_gene_sum_icd_pheat) <- finngen_raw_flt_gene_sum_icd_pheat$Gene
        
        
        
        finngen_raw_flt_gene_sum_icd_pheat
      })
      
      
      
      output$p3_plot9 <- renderPlot({
        finngen_raw_flt_gene_sum_icd_pheat <- finngen_raw_flt_gene_sum_icd_pheat()
        if (nrow(finngen_raw_flt_gene_sum_icd_pheat)>0) {
          finngen_raw_flt_gene_sum_icd_pheat <- finngen_raw_flt_gene_sum_icd_pheat
          ann <- data.frame(finngen_raw_flt_gene_sum_icd_pheat$Group)
          rownames(ann) <- finngen_raw_flt_gene_sum_icd_pheat$Gene
          colnames(ann) <- "Group"
          finngen_raw_flt_gene_sum_icd_pheat <- finngen_raw_flt_gene_sum_icd_pheat[,-c(1,2)]
          
          
          set.seed(1111)
          n1 <- length(unique(ann$Group))
          col1=distinctColorPalette(k = n1, altCol = FALSE, runTsne = FALSE)
          names(col1) <- unique(ann$Group)
          
          ann_colors = list(
            Group = col1
          )
          
          
          ComplexHeatmap::pheatmap(t(finngen_raw_flt_gene_sum_icd_pheat),
                                   angle_col = "45",
                                   name=str_wrap("Significant association within the gene",20),
                                   row_labels = str_wrap(rownames(t(finngen_raw_flt_gene_sum_icd_pheat)),30),
                                   annotation_col = ann,
                                   color = c("white", "red"),
                                   legend = F,
                                   main = "Endpoints: ICD-classification",
                                   annotation_colors = ann_colors)
        } else {
          plot_exception("No significant association within the gene list")
        }
        
        
        
        
        
      }, height = 800, width = 1200
      )
      
      
      # Cluster Endpoints: specific requests
      
      
      finngen_raw_flt_gene_sum_spe_pheat <- eventReactive(input$ViewFinnGenCluster,{
        
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        finngen_raw_flt_gene <- finngen_raw_flt() %>%
          left_join(gene_change)
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>%
          filter(category != "Quantitative endpoints")
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>%
          mutate(Endpoint=ifelse(category %in%
                                   c("I Certain infectious and parasitic diseases (AB1_)",
                                     "XI Diseases of the digestive system (K11_)",
                                     "IV Endocrine, nutritional and metabolic diseases (E4_)",
                                     "II Neoplasms, from cancer register (ICD-O-3)",
                                     "II Neoplasms from hospital discharges (CD2_)",
                                     "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
                                     "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                                     "VII Diseases of the eye and adnexa (H7_)",
                                     "V Mental and behavioural disorders (F5_)",
                                     "IX Diseases of the circulatory system (I9_)",
                                     "VI Diseases of the nervous system (G6_)",
                                     "VIII Diseases of the ear and mastoid process (H8_)",
                                     "X Diseases of the respiratory system (J10_)",
                                     "XII Diseases of the skin and subcutaneous tissue (L12_)",
                                     "XIV Diseases of the genitourinary system (N14_)",
                                     "XV Pregnancy, childbirth and the puerperium (O15_)",
                                     "XVI Certain conditions originating in the perinatal period (P16_)",
                                     "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)",
                                     "XXI Factors influencing health status and contact with health services (Z21_)",
                                     "XIX Injury, poisoning and certain other consequences of external causes (ST19_)" 
                                     
                                   ), "ICD-classification", "Specific requests")
                 
          )
        finngen_raw_flt_gene_sum_spe <- finngen_raw_flt_gene %>% 
          filter(Endpoint=="Specific requests") %>% 
          group_by(Group,Gene,category) %>% 
          summarise(n=length(unique(Variant))) %>% 
          ungroup() %>% 
          mutate(n=ifelse(n>1,1,n))
        
        #Summary table and figure
        finngen_raw_flt_gene_sum_spe <- finngen_raw_flt_gene_sum_spe
        multi_gene <- finngen_raw_flt_gene_sum_spe %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          group_by(Gene) %>% 
          mutate(n=n()) %>% 
          filter(n>1) %>% 
          ungroup()
        
        
        finngen_raw_flt_gene_sum_spe_pheat <- finngen_raw_flt_gene_sum_spe %>% 
          pivot_wider(
            names_from = "category",
            values_from = "n",
            values_fill = 0
          )
        
        finngen_raw_flt_gene_sum_spe_pheat <- finngen_raw_flt_gene_sum_spe_pheat %>% 
          mutate(Gene=ifelse((Gene %in% multi_gene$Gene),paste0(Gene,"_",Group),Gene)) %>% 
          as.data.frame()
        rownames(finngen_raw_flt_gene_sum_spe_pheat) <- finngen_raw_flt_gene_sum_spe_pheat$Gene
    
        finngen_raw_flt_gene_sum_spe_pheat
      })
      
      
      
      output$p3_plot10 <- renderPlot({
        
        finngen_raw_flt_gene_sum_spe_pheat <- finngen_raw_flt_gene_sum_spe_pheat()
        if (nrow(finngen_raw_flt_gene_sum_spe_pheat)>0) {
          finngen_raw_flt_gene_sum_spe_pheat <- finngen_raw_flt_gene_sum_spe_pheat
          ann <- data.frame(finngen_raw_flt_gene_sum_spe_pheat$Group)
          rownames(ann) <- finngen_raw_flt_gene_sum_spe_pheat$Gene
          colnames(ann) <- "Group"
          finngen_raw_flt_gene_sum_spe_pheat <- finngen_raw_flt_gene_sum_spe_pheat[,-c(1,2)]
          
          
          set.seed(1111)
          n1 <- length(unique(ann$Group))
          col1=distinctColorPalette(k = n1, altCol = FALSE, runTsne = FALSE)
          names(col1) <- unique(ann$Group)
          
          ann_colors = list(
            Group = col1
          )
          
          ComplexHeatmap::pheatmap(t(finngen_raw_flt_gene_sum_spe_pheat),
                                   angle_col = "45",
                                   name=str_wrap("Significant association within the gene",20),
                                   row_labels = str_wrap(rownames(t(finngen_raw_flt_gene_sum_spe_pheat)),30),
                                   annotation_col = ann,
                                   color = c("white", "red"),
                                   legend = F,
                                   main = "Endpoints: Specific requests",
                                   annotation_colors = ann_colors)
        } else {
          plot_exception("No significant association within the gene list")
        }
        
        
        
        
        
      }, height = 800, width = 1200
      )
      
      
      
      
      
      #finngen_raw_flt_gene

      finngen_raw_flt_gene <- eventReactive(input$ViewFinnGenPheno,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        finngen_raw_flt_gene <- finngen_raw_flt() %>% 
          left_join(gene_change)
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>% 
          filter(category != "Quantitative endpoints")
        #select the lowest p value SNP for each gene each phenotype
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>% 
          group_by(Group,Gene,category,phenotype) %>% 
          arrange(pval) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          dplyr::select(Group,category,phenotype,Gene,rsids) %>% 
          group_by(Group,category) %>% 
          mutate(
            n=length(unique(phenotype)),
            Phenotype=paste0(unique(phenotype),collapse = "; "),
            rsids = paste0(unique(rsids),collapse = "; "),
            Gene = paste0(unique(Gene),collapse = "; ")) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        finngen_raw_flt_gene <- finngen_raw_flt_gene %>% 
          mutate(Endpoint=ifelse(category %in% 
                                   c("I Certain infectious and parasitic diseases (AB1_)",
                                     "XI Diseases of the digestive system (K11_)",
                                     "IV Endocrine, nutritional and metabolic diseases (E4_)",
                                     "II Neoplasms, from cancer register (ICD-O-3)",
                                     "II Neoplasms from hospital discharges (CD2_)",
                                     "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
                                     "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                                     "VII Diseases of the eye and adnexa (H7_)",
                                     "V Mental and behavioural disorders (F5_)",
                                     "IX Diseases of the circulatory system (I9_)",
                                     "VI Diseases of the nervous system (G6_)",
                                     "VIII Diseases of the ear and mastoid process (H8_)",
                                     "X Diseases of the respiratory system (J10_)",
                                     "XII Diseases of the skin and subcutaneous tissue (L12_)",
                                     "XIV Diseases of the genitourinary system (N14_)",
                                     "XV Pregnancy, childbirth and the puerperium (O15_)",
                                     "XVI Certain conditions originating in the perinatal period (P16_)",
                                     "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)",
                                     "XXI Factors influencing health status and contact with health services (Z21_)"
                                     
                                   ), "ICD-classification", "Specific requests")
                 
          )
        

        finngen_raw_flt_gene
      })



      #finngen_raw_flt_gene_pro
      finngen_raw_flt_gene_pro  <- eventReactive(input$ViewFinnGenPheno,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        gene_number <- gene_change %>% 
          group_by(Group) %>% 
          mutate(GeneNumber=length(unique(Gene))) %>% 
          ungroup()
        
        finngen_raw_flt_gene_pro <- finngen_raw_flt() %>% 
          left_join(gene_number)
        
        
        finngen_raw_flt_gene_pro <- finngen_raw_flt_gene_pro %>% 
          filter(category != "Quantitative endpoints")
        
        #select the lowest p value SNP for each gene each phenotype
        finngen_raw_flt_gene_pro <- finngen_raw_flt_gene_pro %>% 
          group_by(Group,Gene,category,phenotype) %>% 
          arrange(pval) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          dplyr::select(Group,GeneNumber,category,phenotype,Gene,rsids) %>% 
          group_by(Group,GeneNumber,category) %>% 
          mutate(
            n=length(unique(Gene)),
            Phenotype=paste0(unique(phenotype),collapse = "; "),
            rsids = paste0(unique(rsids),collapse = "; "),
            Gene = paste0(unique(Gene),collapse = "; ")) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          mutate(GeneFraction=n/GeneNumber*100)
        
        
        
        finngen_raw_flt_gene_pro <- finngen_raw_flt_gene_pro %>% 
          mutate(Endpoint=ifelse(category %in% 
                                   c("I Certain infectious and parasitic diseases (AB1_)",
                                     "XI Diseases of the digestive system (K11_)",
                                     "IV Endocrine, nutritional and metabolic diseases (E4_)",
                                     "II Neoplasms, from cancer register (ICD-O-3)",
                                     "II Neoplasms from hospital discharges (CD2_)",
                                     "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
                                     "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                                     "VII Diseases of the eye and adnexa (H7_)",
                                     "V Mental and behavioural disorders (F5_)",
                                     "IX Diseases of the circulatory system (I9_)",
                                     "VI Diseases of the nervous system (G6_)",
                                     "VIII Diseases of the ear and mastoid process (H8_)",
                                     "X Diseases of the respiratory system (J10_)",
                                     "XII Diseases of the skin and subcutaneous tissue (L12_)",
                                     "XIV Diseases of the genitourinary system (N14_)",
                                     "XV Pregnancy, childbirth and the puerperium (O15_)",
                                     "XVI Certain conditions originating in the perinatal period (P16_)",
                                     "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)",
                                     "XXI Factors influencing health status and contact with health services (Z21_)"
                                     
                                   ), "ICD-classification", "Specific requests")
                 
          )
        
        

        finngen_raw_flt_gene_pro
      })


      output$p3_plot11 = renderPlotly({
        
        finngen_raw_flt_gene <- finngen_raw_flt_gene()
        if (nrow(finngen_raw_flt_gene)>0) {
          p <- finngen_raw_flt_gene %>% 
            ggplot(aes(x=n,y=category,
                       #color=MultiGene,
                       fill=Group,
                       label=Gene,label2=rsids,label3=Phenotype
            ))+
            geom_col()+
            facet_grid(Endpoint~Group,scales = "free_y",space = "free") +
            theme(
              legend.position = "none"
            )+
            # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            #    theme(strip.text.y = element_text(size = 1),
            #          legend.position = "none"
            #          ) +
            labs(x="The unique phenotype number",y="Endpoints (Phenotype)") +
            theme_bw()
          
          #ggplotly format
          ggplotly(p,tooltip=c("Gene", "rsids","Phenotype")) %>%
            layout(autosize = T, width = 1200, height = 800)
        } else {
          plotly_exception("No significant association within the gene list")
        }
        
        
        
        



      })


      output$p3_plot12 = renderPlotly({
        finngen_raw_flt_gene_pro <- finngen_raw_flt_gene_pro()
        if (nrow(finngen_raw_flt_gene_pro)>0) {
          p <- finngen_raw_flt_gene_pro %>% 
            ggplot(aes(x=GeneFraction,y=category,
                       #color=MultiGene,
                       fill=Group,
                       label=Gene,label2=rsids,label3=Phenotype
            ))+
            geom_col()+
            facet_grid(Endpoint~Group,scales = "free_y",space = "free") +
            theme(
              legend.position = "none"
            )+
            # facet_grid(Category~Oddsratio,scales = "free",space = "free") +
            #    theme(strip.text.y = element_text(size = 1),
            #          legend.position = "none"
            #          ) + 
            labs(x="Percentage of genes in each gene set group",y="Endpoints (Phenotype)") +
            theme_bw()
          
          #ggplotly format
          ggplotly(p,tooltip=c("Gene", "rsids","Phenotype")) %>%
            layout(autosize = T, width = 1200, height = 800)
        } else {
          plotly_exception("No significant association within the gene list")
        }
        
      })

      
      
      
      
      
      # Variant-Phenotype Gene Effect
      
      finngen_raw_flt_gene_effect <- eventReactive(input$ViewFinnGenEffect,{
        gene_raw <- gene_raw()
        gene_change <- gene_raw %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        finngen_raw_flt_gene_effect <- finngen_raw_flt() %>% 
          left_join(gene_change)
        finngen_raw_flt_gene_effect <- finngen_raw_flt_gene_effect %>% 
          filter(category != "Quantitative endpoints")
        
        #select the lowest p value SNP for each gene each phenotype
        finngen_raw_flt_gene_effect <- finngen_raw_flt_gene_effect %>% 
          group_by(Group,Gene,category,phenotype) %>% 
          arrange(pval) %>% 
          dplyr::slice(1) %>% 
          ungroup() %>% 
          dplyr::select(Group,category,phenotype,Gene,rsids,beta) %>% 
          group_by(Group,category,Gene) %>%
          mutate(MeanEffect=mean(abs(beta))) %>% 
          mutate(
            Phenotype=paste0(unique(phenotype),collapse = "; "),
            rsids = paste0(unique(rsids),collapse = "; ")
          ) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        finngen_raw_flt_gene_effect <- finngen_raw_flt_gene_effect %>% 
          mutate(Endpoint=ifelse(category %in% 
                                   c("I Certain infectious and parasitic diseases (AB1_)",
                                     "XI Diseases of the digestive system (K11_)",
                                     "IV Endocrine, nutritional and metabolic diseases (E4_)",
                                     "II Neoplasms, from cancer register (ICD-O-3)",
                                     "II Neoplasms from hospital discharges (CD2_)",
                                     "III Diseases of the blood and blood-forming organs and certain disorders involving the immune mechanism (D3_)",
                                     "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
                                     "VII Diseases of the eye and adnexa (H7_)",
                                     "V Mental and behavioural disorders (F5_)",
                                     "IX Diseases of the circulatory system (I9_)",
                                     "VI Diseases of the nervous system (G6_)",
                                     "VIII Diseases of the ear and mastoid process (H8_)",
                                     "X Diseases of the respiratory system (J10_)",
                                     "XII Diseases of the skin and subcutaneous tissue (L12_)",
                                     "XIV Diseases of the genitourinary system (N14_)",
                                     "XV Pregnancy, childbirth and the puerperium (O15_)",
                                     "XVI Certain conditions originating in the perinatal period (P16_)",
                                     "XVII Congenital malformations, deformations and chromosomal abnormalities (Q17)",
                                     "XXI Factors influencing health status and contact with health services (Z21_)"
                                     
                                   ), "ICD-classification", "Specific requests")
                 
          )
        
        finngen_raw_flt_gene_effect
        
        
        
      })
      
      
      output$p3_plot_add2 <- renderPlot({
        finngen_raw_flt_gene_effect <- finngen_raw_flt_gene_effect()
        if (nrow(finngen_raw_flt_gene_effect)) {
          finngen_raw_flt_gene_effect <- finngen_raw_flt_gene_effect
          
          
          finngen_raw_flt_gene_effect %>%
            mutate(MeanEffect=ifelse(MeanEffect>10,10,MeanEffect)) %>% 
            ggplot(aes(x=log10(MeanEffect),y=category,
                       fill=Group,
                       color=Group,
                       label=Gene
            ))+
            geom_jitter(position = position_dodge2(0.9))+
            facet_grid(Endpoint~.,scales = "free",space = "free") +
            theme_bw()+
            labs(x="Log10(MeanEffect) (Odds ratio or Effect size)",y="Phenotype Category")+
            ggrepel::geom_text_repel(position = position_dodge2(0.9))
        } else {
          plot_exception("No significant association within the gene list")
        }
         
        
        
      }, height = 850, width = 1500
      )
      
      
      
      # Summary Table Output
      output$p3_dt3 <- DT::renderDataTable(
        
        finngen_raw_flt_gene_effect() %>% 
          dplyr::select(Group,Gene,category,Phenotype,rsids,Endpoint,MeanEffect) %>% 
          mutate_if(is.numeric, round, digits = 3)  %>% 
        DT::datatable(escape = FALSE,
                      rownames = FALSE,
                      plugins = "ellipsis",
                      extensions = 'Buttons',
                      options = list(
                        
                        buttons = list(
                          list(extend = "csv", text = "Download Full Results", filename = "Full_data",
                               exportOptions = list(
                                 modifier = list(page = "all"),
                                 orthogonal = "export"
                               )
                          )
                        ),
                        dom = 'Bfrtip',
                        buttons = c('copy', 'csv'),
                        columnDefs = list(list(
                          targets = c(3:4),
                          render = JS("$.fn.dataTable.render.ellipsis(17, false )")
                        ))
                      )) 
        
        
        
        
      )
      
      
      
      
  ####################################################Page 4 HPO Phenotypes####################################################
      
      
      
      observeEvent(input$RunAnalysis, {
        
        gene_raw <- gene_raw()
        group <- unique(gene_raw$Group)
        updateSelectInput(session, "GeneListGroup", label = "Select", choices = group)
      })
      
      
      
      hpo_gene  <- eventReactive(input$ViewHPOPlot,{
        gene_raw <- gene_raw()
        hpo_gene <- gene_raw %>% 
          filter(Group==input$GeneListGroup)
        hpo_gene <- hpo_gene$Gene
      
        hpo_gene
      })
      
      
      output$p4_plot1 = renderPlotly({
        hpo_gene <- hpo_gene()
        phedb <- c("HPO")
        
        # Run PhenoEnrichGenes to get the results
        hpo_gene_enrich <- PhenoEnrichGenes(genes= hpo_gene, database = phedb) 
        
        # Show the graph
        hpo_gene_enrich$graph
      })
      
      
      #Enrichment Table
      
      observeEvent(input$RunAnalysis, {
        gene_raw <- gene_raw()
        group <- unique(gene_raw$Group)
        updateSelectInput(session, "GeneListGroupTab", label = "Select", choices = group)
      })
      
      
    
      hpo_gene_tab  <- eventReactive(input$ViewHPOTab,{
        gene_raw <- gene_raw()
        hpo_gene <- gene_raw %>% 
          filter(Group==input$GeneListGroupTab)
        hpo_gene <- hpo_gene$Gene
        phedb <- c("HPO")
        # Run PhenoEnrichGenes to get the results
        hpo_gene_enrich <- PhenoEnrichGenes(genes= hpo_gene, database = phedb) 
        # Save the table
        hpo_gene_table <- hpo_gene_enrich$htmltabla
        hpo_gene_table
      })
      
      
      # Summary Table Output
      output$p4_dt1 <- DT::renderDataTable(

        hpo_gene_tab() %>% 
        DT::datatable(
                      caption = 'This table shows the gene associated with different traits.',
                      filter = 'top',
                      extensions = 'Buttons',
                      escape = FALSE,
                      options = list(dom = 'Blfrtip',
                                     buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
                                     lengthMenu = list(c(10,25,50,-1),
                                                       c(10,25,50,"All")))
        )
        
        
        
      )
      
      
      
      #Compare plot
      
      observeEvent(input$RunAnalysis, {
        gene_raw <- gene_raw()
        group <- unique(gene_raw$Group)
        updateSelectInput(session, "Group1", label = "Select", choices = group)
      })
      
      
      observeEvent(input$RunAnalysis, {
        gene_raw <- gene_raw()
        group <- unique(gene_raw$Group)
        updateSelectInput(session, "Group2", label = "Select", choices = group)
      })
      
      
      
      
      comparare_hpo_group  <- eventReactive(input$ViewHPOCompar,{
        gene_raw <- gene_raw()
        phedb <- c("HPO")
        
        hpo_gene_group1 <- gene_raw() %>% 
          filter(Group==input$Group1)
        hpo_gene_group1 <- hpo_gene_group1$Gene
        
        
        #Different group
        hpo_gene_group2 <- gene_raw() %>% 
          filter(Group==input$Group2)
        hpo_gene_group2 <- hpo_gene_group2$Gene
        
        # Run the anaysis
        comparare_hpo_group <- RandomComparePheno(geneset = hpo_gene_group1, genesetcompare = hpo_gene_group2, database = phedb, nulltestnumber = 50)

        comparare_hpo_group
      })
      
      
      
      
      output$p4_plot2 = renderPlotly({
        comparare_hpo_group <- comparare_hpo_group()
        comparare_hpo_group$plotdif
        
      })
      
      
      
      
  ####################################################Page 5 GTEx Analysis####################################################
      
      #ViewGTExBoxplot
      
      
      gtex_gene  <- eventReactive(input$ViewGTExBoxplot,{
      
        gtex_gene <- gene_raw() %>% 
          left_join(gtex_raw,by=c("Gene"="Description"))
        
        gtex_gene <- gtex_gene %>% 
          dplyr::select(-Name) %>% 
          #mutate(id = seq_len(n())) %>% 
          pivot_longer(-c("Group","Gene"),
                       names_to = 'tissue',
                       values_to = 'expression') %>% 
          mutate(tissue = str_replace(tissue, ' -', ':') %>% str_remove('\\s?\\(.+\\)'),
                 tissue_class = coalesce(str_extract(tissue, '^.+(?=:)'),
                                         tissue)) %>%
          mutate(expression=log10(expression+1))
        
        gtex_gene
      })
      
      
      output$p5_plot1 = renderPlotly({
        
        
        gtex_gene() %>%
          group_by(Group) %>%
          do(p=plot_ly(., y=~tissue,x=~expression,color=~tissue_class,type="box",text = ~paste(Group, Gene)))  %>%
          subplot(nrows = 1, shareX = TRUE, shareY = TRUE)  %>% 
          layout(xaxis = list(title = "Median gene expression levels log10(TPM+1)"), yaxis = list(title = "Tissue"),autosize = T, width = 1200, height = 1200) 
        
        
      })
      
      
      
      
      #ViewGTExCompar
    
      
      
      gtex_geneset  <- eventReactive(input$ViewGTExCompar,{
        
        gtex_geneset <- gene_raw() %>% 
          left_join(gtex_raw,by=c("Gene"="Description"))
        
        gtex_geneset <- gtex_geneset %>% 
          dplyr::select(-Name) %>% 
          #mutate(id = seq_len(n())) %>% 
          pivot_longer(-c("Group","Gene"),
                       names_to = 'tissue',
                       values_to = 'expression') %>% 
          mutate(tissue = str_replace(tissue, ' -', ':') %>% str_remove('\\s?\\(.+\\)'),
                 tissue_class = coalesce(str_extract(tissue, '^.+(?=:)'),
                                         tissue)) %>%
          mutate(expression=log10(expression+1))
        gtex_geneset
        
      })
      
      
      output$p5_plot2 = renderPlot({
        
        gtex_geneset <- gtex_geneset()
        
        # Step 1: Perform Statistical Testing
        stat.test <- tryCatch({
          gtex_geneset %>%
            group_by(tissue) %>%
            
            # Perform t-test between 'Group' and 'expression', adjusting p-values using Bonferroni
            t_test(expression ~ Group, p.adjust.method = "bonferroni") %>%
            
            # Add significance markers based on adjusted p-values
            add_significance()
        }, error = function(e) {
          plot_exception("Not enough observations for comparison")
          return(NULL)  # Return NULL if an error occurs
        })
        
        # Step 2: Check if stat.test was successful and proceed with plotting
        if (!is.null(stat.test)) {
          
          # Create Box Plot
          bxp <- ggboxplot(
            gtex_geneset, x = "Group", y = "expression", 
            fill = "Group",
            xlab = "Group",
            ylab = "Median gene expression levels log10(TPM+1)",
            facet.by = "tissue"
          )
          
          # Add P-values to the Plot
          stat.test <- stat.test %>% add_xy_position(x = "Group")
          p <- bxp + stat_pvalue_manual(stat.test)
          
          # Print the figure if it's generated
          print(p)
          
        } else {
          # If stat.test is NULL, the error message has already been printed
          print("Not enough observations for comparison")
        }
        
        
        
      } , height = 800, width = 1200)
      
      
      
      
      
      #ViewGTExCluster
      
      
      gtex_gene_pheat_flt  <- eventReactive(input$ViewGTExCluster,{
        gtex_gene_pheat_flt <- gene_raw() %>% 
          left_join(gtex_raw,by=c("Gene"="Description"))
        
        gtex_gene_pheat_flt <- gtex_gene_pheat_flt %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        gtex_gene_pheat_flt <- gtex_gene_pheat_flt %>% 
          dplyr::select(-c("Group","Name")) %>% 
          as.data.frame()
        rownames(gtex_gene_pheat_flt) <- gtex_gene_pheat_flt$Gene
        gtex_gene_pheat_flt <- gtex_gene_pheat_flt[,-1]
        gtex_gene_pheat_flt <- t(gtex_gene_pheat_flt)
        gtex_gene_pheat_flt <- log10(gtex_gene_pheat_flt+1)
        gtex_gene_pheat_flt[is.na(gtex_gene_pheat_flt)] <- 0
        gtex_gene_pheat_flt
        
        
      })
      
      
      output$p5_plot3 = renderPlot({
        
        gtex_gene_pheat_flt <- gtex_gene_pheat_flt()
        annotation_col = data.frame(
          Gene=colnames(gtex_gene_pheat_flt)
        )
        
        gene_change <- gene_raw() %>% 
          group_by(Gene) %>% 
          mutate(Group=paste0(unique(Group),collapse = "; ")) %>% 
          ungroup() %>% 
          group_by(Group,Gene) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        annotation_col <- annotation_col %>% 
          left_join(gene_change)
        rownames(annotation_col) <- annotation_col$Gene
        annotation_col <- annotation_col %>% 
          dplyr::select(-Gene)
        
        
        tissue_class <- gene_raw() %>% 
          left_join(gtex_raw,by=c("Gene"="Description"))
        
        tissue_class <- tissue_class %>% 
          dplyr::select(-Name) %>% 
          #mutate(id = seq_len(n())) %>% 
          pivot_longer(-c("Group","Gene"),
                       names_to = 'tissue',
                       values_to = 'expression') %>% 
          mutate(tissue = str_replace(tissue, ' -', ':') %>% str_remove('\\s?\\(.+\\)'),
                 tissue_class = coalesce(str_extract(tissue, '^.+(?=:)'),
                                         tissue)) %>% 
          dplyr::select(tissue,tissue_class) %>% 
          group_by(tissue) %>% 
          dplyr::slice(1) %>% 
          ungroup()
        
        
        annotation_row = data.frame(
          tissue = rownames(gtex_gene_pheat_flt)
        ) %>% 
          mutate(tissue = str_replace(tissue, ' -', ':') %>% str_remove('\\s?\\(.+\\)'))
        
        annotation_row <- annotation_row %>% 
          left_join(tissue_class)
        rownames(annotation_row) <- annotation_row$tissue
        annotation_row <- annotation_row %>% 
          dplyr::select(-tissue)
        
        
        
        set.seed(1111)
        n1 <- length(unique(annotation_col$Group))
        n2 <- length(unique(annotation_row$tissue_class))
        col1=distinctColorPalette(k = n1, altCol = FALSE, runTsne = FALSE)
        names(col1) <- unique(annotation_col$Group)
        col2=distinctColorPalette(k = n2, altCol = FALSE, runTsne = FALSE)
        names(col2) <- unique(annotation_row$tissue_class)
        
        
        ann_colors = list(
          Group = col1,
          tissue_class = col2
        )
        
        p <- pheatmap(gtex_gene_pheat_flt, annotation_col = annotation_col, annotation_row = annotation_row,name = "log10(TPM+1)",
                      annotation_colors = ann_colors)
        
        draw(p, legend_grouping = "original")
        
        
      } , height = 1000, width = 1400)
      
      
  

  ###########################Report###########################
  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    #filename = "report.html",
    filename = function() {
      paste('my-report.html')
    },
    
    content = function(file) {
      tempReport <- file.path("Data/report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      #params <- list(n = input$search)
      params <- list(n = input$file1)
     
      #params <- read.csv(file$datapath, header = input$header)
      
      rmarkdown::render(tempReport, 
                        output_file = file,
                        #HTML = html_document(),
                        # switch(
                        #   input$format,
                        #   PDF = pdf_document(), HTML = html_document(), Word = word_document()
                        # ),

                        params = params,
                        envir = new.env(parent = globalenv())
      )
    }
  )
  

  
}

# Create Shiny app ----
shinyApp(ui, server)


