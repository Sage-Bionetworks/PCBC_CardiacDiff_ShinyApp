library(shiny)
library(shinyIncubator)


shinyServer(function(input, output, session) {
  
  genes_in_selectedPathways <- reactive({
   
  })
  
  #get the list of user submitted genes
  user_submitted_geneList <- reactive({
    input$searchByGenes_button
    selected_genes <- isolate(input$custom_gene_list)
    geneList <- unlist(strsplit(selected_genes,split=c('[\\s+,\\n+\\r+)]'),perl=T))
    #conevert everything to upper case
    geneList <- toupper(geneList)
    geneList <- geneList[ !geneList == "" ] #remove the blank entries
    geneList
  })
  
  #list of genes : based on the tab opened by the user
  selected_genes <- reactive({
    if( input$genelist_type == 'custom_gene_list'  ){
      genes <- unique(user_submitted_geneList())
    }
    else if( input$genelist_type == 'precomputed_geneList'){
     rows_to_keep <- enriched_geneLists$geneListName  %in% input$selected_geneLists
     genes <- as.character(enriched_geneLists[rows_to_keep, 'ensemblId'])
    }
    #list present in the selected pathways
    else if( input$genelist_type == 'pathway'){
      genes <- as.character(unlist(pathways_list[input$selected_pathways]))
      }
    genes
  })
  
  #get the ensembl id's for the selected genes
  selected_ensemblIds <- reactive({
    filtered_df <- hg19_gene_annot %>%
      filter(SYMBOL %in% selected_genes()  | ENSEMBL %in% selected_genes() |
             ENSEMBLTRANS %in% selected_genes()  | ENTREZID %in% selected_genes())
    #na.omit(unique(c(filtered_df$ENSEMBL, filtered_df$ENSEMBLTRANS)))
    na.omit(unique(c(filtered_df$ENSEMBL)))
  })
  
  
  #get the HUGO ids for the selected genes
  selected_HUGOIds <- reactive({
    filtered_df <- hg19_gene_annot %>%
      filter(SYMBOL %in% selected_genes()  | ENSEMBL %in% selected_genes() |
               ENSEMBLTRANS %in% selected_genes()  | ENTREZID %in% selected_genes())
    #na.omit(unique(c(filtered_df$ENSEMBL, filtered_df$ENSEMBLTRANS)))
    na.omit(unique(c(filtered_df$SYMBOL)))
  })
  

  output$test <- renderPrint({
#     print(length(selected_ensemblIds()))
#      print(dim(filtered_expMatrix()))
#     print(input$selected_timepoints)
#     print(selected_samples())
  })
  

  selected_samples <- reactive({
    #filter based on timepoint
    if(length(input$selected_timepoints) > 0){
      selected_samples <- metadata %>% filter(TimePoint %in% input$selected_timepoints)
      selected_samples <- selected_samples[['Sample.Name']]
    }
  })

  
  filtered_expMatrix <- reactive({
    flt_mat <- expMat
    
    #fitler based on genes
     if( length(selected_ensemblIds()) > 0 ){
       rows_to_keep <- rownames(flt_mat) %in% selected_ensemblIds()
       flt_mat <- flt_mat[rows_to_keep,]
     }
    
    #filter based on samples
    if (! is.null(selected_samples()) ){
      cols_to_keep <- colnames(flt_mat) %in% selected_samples()
      flt_mat <- flt_mat[, cols_to_keep]
    }
    flt_mat
  })
  
  
  filtered_annotation <- reactive({
    flt_metadata <- metadata['TimePoint',drop=F]
    rownames(flt_metadata) <- metadata[['Sample.Name']]
    flt_metadata
  })
  
  
  #heatmap for pubData pathway enrichment 
  output$heatmap <- renderPlot({
    m <- filtered_expMatrix()
    #remove any rows with NA values
    m <- na.omit(m)
    #remove any rows with 0 variance
    m <- m[apply(m,1,var) != 0, ]
    
    validate( need(nrow(m) < 10000, "Filtered expression matrix contains > 1000 genes ") )
    fontsize_row=8
    if(nrow(m) > 100){ fontsize_row = 0 }
    m = t(scale(t(m)))
    withProgress(session, {
      setProgress(message = "clustering & rendering heatmap, please wait", 
                  detail = "This may take a few moments...")
      matrix_heatMap(m, filtered_annotation(),
                     fontsize_col=0, 
                     fontsize_row=fontsize_row,
                     clustering_method = "average"
                     )
    })
  })
  
  
  
  output$topgene_linkOut <- reactive({
    prefix =  '<form action="https://toppgene.cchmc.org/CheckInput.action" method="post" target="_blank" display="inline">
    <input type="hidden" name="query" value="TOPPFUN">
    <input type="hidden" id="type" name="type" value="HGNC">
    <input type="hidden" name="training_set" id="training_set" value="'
    suffix = '">
    <input type="Submit" class="btn shiny-download-link shiny-bound-output", value="Enrichment Analysis in ToppGene">
    </form>'
    genes <- paste(selected_HUGOIds(),collapse=" ")
    #generate the HTML content
    htmlContent <- paste(c(prefix,genes,suffix), collapse="")
  })
  
  
  #handle data download request
  output$downloadData <- downloadHandler(
    filename = function() { paste('PCBC_Cardiac_geneExpr_data.csv')},
    content  = function(file){
      write.csv(filtered_expMatrix(),file,row.names=T)
    }
  )
  
  output$geneExpTable <- renderDataTable({ 
    df <- filter(hg19_grpd, ENSEMBL %in% rownames(filtered_expMatrix()))
  })

})