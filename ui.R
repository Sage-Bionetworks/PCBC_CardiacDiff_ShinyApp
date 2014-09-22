library(shiny)
library(shinyIncubator)

shinyUI(fluidPage(
  
  sidebarLayout(
    #SIDE BAR PANEL FOR USER OPTIONS
    sidebarPanel(
      progressInit(),
      h4('1. Select a gene list'),
      
      tabsetPanel(
        id = 'genelist_type',
        
        #TAB PANEL 1 : custom gene list
        tabPanel('My gene list', h5('Search on a custom gene list:'),
                 tags$textarea(id="custom_gene_list",rows=8,cols=200,sample_gene_list),
                 helpText("Accepts HUGO/Ensembl/Entrez gene id's separated either by comma, space, line"),
                 actionButton("searchByGenes_button", "Search by Genes"),
                 value='custom_gene_list'), #END TAB PANEL 1
        
        
        # Tab Panel 2 : select a pathway
        tabPanel(
          'Pathways',
          selectInput("selected_pathways",
                      "Pathways",
                      choices = names(pathways_list), selectize=T, multiple=T, width='400px',
                      selected = names(pathways_list)[c(1:2)]),
          br(), br(), br(), br(), br(), br(),
          value='pathway' ), #END TAB PANEL 2
        
        
      
        #TAB PANEL 3: precomputed Enriched genelists
        tabPanel( 'Sig gene lists',
                  selectInput("selected_geneLists",
                              "Precomputed Significant gene lists", selectize=T, 
                              choices =  geneLists,
                              selected = geneLists[1],
                              multiple=F),
                  br(),br(),br(),br(),
                  br(),br(),br(),br(),
                  value='precomputed_geneList'
        ) #END Tab Panel 3
      ), #END tabsetPanel
      
      
      br(), br(),
      h4('2. Filter samples by:'),
      #1. heatmap annotation labels
      selectInput("selected_timepoints","Time Points", 
                  choices =  metadata$TimePoint, selected=c('D60','D03'),
                  selectize=T, multiple=T),
      br(), br(),
      
      h4('3. Heatmap Settings:'),
      #distance metric
      selectInput("clustering_distance","Distance Calculation",
                  choices=c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski"),
                  selectize=T, multiple=F, selected="correlation"),
      #linkage 
      selectInput("clustering_method","Clustering Method",
                  choices=c("ward", "single", "complete", "average", "mcquitty", "median", "centroid"),
                  selectize=T, multiple=F, selected="complete")
      
    ),#END sideBarPanel
    
    mainPanel(
      tabsetPanel(
        #tab panel 1
        tabPanel("expression HeatMap", 
                 plotOutput("heatmap",height="700px",width="auto",hoverId=NULL),
                 br(),
                 br(),
                 h4('Summary'),
                 tableOutput("summary")
        ), #END tabPanel 1
        tabPanel("Explore Data", 
                 htmlOutput("topgene_linkOut"),
                 downloadButton('downloadData','Download Expression Data'),
                 br(), br(), br(),
                 dataTableOutput("geneExpTable")
      )#END tabset panel
    ) #END tabsetPanel under mainPanel
    ) #END mainPanel
  ),# END SidebarLayout
  responsive=T,
  title="PCBC Cardiac Differentiation Study Data Explorer"
  )
)