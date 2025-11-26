#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(factoextra)


clean_data <- read.csv("./clean_stage2.csv", header = TRUE, sep = ",")

# visualization pvalue by Parameter Name for Gene
  # convert all gene_symbols to lower case and ensure distinct gene_symbols are stored 
clean_data$gene_symbol <- tolower(clean_data$gene_symbol)
genes <- unique(clean_data$gene_symbol)
genes <- as.vector(genes) 

  # transform pvalue using -log10(pvalue) for visualization
clean_data$pvalue_transformed <- -log10(clean_data$pvalue)

  # function to subset the data for a selected gene
subset_gene <- function(gene){
  df_subsetted_gene <- clean_data[clean_data$gene_symbol == gene, ]
  df_subsetted_gene
}

  # function to visualize the parameter_name pvalue for selected gene
plot_gene_knockout <- function(df_subsetted_gene, gene){
  df_subsetted_gene$parameter_name <- reorder(df_subsetted_gene$parameter_name, -df_subsetted_gene$pvalue_transformed) # reorder the parameter_name by pvalue
  ggplot(df_subsetted_gene,
         aes(x = parameter_name, y = pvalue_transformed)) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # adds dashed red line to indicate pvalue of 0.05 
    labs(
      title = sprintf("pvalue by Parameter Name for %s", gene),
      x = "Parameter Name",
      y = "Transformed pvalue -log10(pvalue)"
    )
}


# visualization pvalue by Gene Knockout for Phenotype
  # select all phenotype/parameter_name from dataset
phenotypes <- as.vector(clean_data$parameter_name) 

  # function to subset the data for a selected phenotype 
subset_phenotype <- function(phenotype){
  df_subsetted_phenotype <- clean_data[clean_data$parameter_name == phenotype, ]
  df_subsetted_phenotype
}

  # function to visualize the gene_symbol pvalue for a selected phenotype
plot_knockout_phenotype <- function(df_subsetted_phenotype, phenotype){
  df_subsetted_phenotype$gene_symbol <- reorder(df_subsetted_phenotype$gene_symbol, -df_subsetted_phenotype$pvalue_transformed) # reorder the gene_symbol by pvalue
  ggplot(df_subsetted_phenotype,
         aes(x = gene_symbol, y = pvalue_transformed)) + 
    geom_point() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") + # adds a dashed red line to indicate pvalue of 0.05 
    labs(
      title = sprintf("pvalue by Gene Knockout for %s", phenotype),
      x = "Gene Symbol",
      y = "Transformed pvalue -log10(pvalue)"
    )
}

# visualization PCA + Clustering of Genes with Similar Phenotype Scores
  # subset data and select gene_symbol, parameter_name, pvalue 
gp_data <- subset(clean_data, select = c(gene_symbol, parameter_name, pvalue))

gp_data_agg <- gp_data %>%
  group_by(gene_symbol, parameter_name) %>%   # group by both gene and phenotype
  summarize(pvalue = min(pvalue, na.rm = TRUE), .groups = "drop") %>% # takes the minimum pvalue for each gene-phenotype combination
  # converts the data to a wide format with rows = genes, columns = phenotypes and cells = minimum pvalue
  pivot_wider(
    names_from = parameter_name,
    values_from = pvalue,
    values_fill = 1   # fill missing p-values with 1 -> no significance 
  )

  # saves gene_names present in the gp_data_agg dataframe 
gene_names <- gp_data_agg$gene_symbol

  # convert to dataframe
gpmat <- as.data.frame(gp_data_agg)

  # sets gene_symbols as row names 
rownames(gpmat) <- gpmat$gene_symbol

gpmat <- subset(gpmat, select = -c(gene_symbol))

gpmat_pca <- gpmat %>% 
  mutate(across(everything(), ~ -log1p(.))) # transform pvalue with -log1p(pvalue) for before clustering so larger values indicate higher significance

  # applying principal component analysis (PCA) to reduce dimensionality before clustering using KMeans clustering 
pca_clustering <- function(gpmat_pca, k){
  pca <- prcomp(gpmat_pca, scale. = TRUE)
  numPCs <- min(50, ncol(pca$x))
  km <- kmeans(pca$x[, 1:numPCs], centers = k) # runs the kmeans clustering algorithm on the first 50 principal components
  list(pca = pca, km = km) 
}

  # visualizes the first two PCs | each point is a gene, and gene_symbol are added, clusters are indicated by color 
cluster_vis <- function(pca, km, gene_names){
  pca_df <- data.frame(
    gene = gene_names,
    PC1  = pca$x[,1],
    PC2  = pca$x[,2],
    cluster = factor(km$cluster)
  )
  
  ggplot(pca_df, aes(x = PC1, y = PC2, color = cluster)) +
    geom_point(alpha = 0.6) +
    theme_minimal() +
    geom_text(aes(label = gene), vjust = -0.4) +
    labs(title = "K-means clusters of Genes in PCA space")
}


# Define UI for application that draws a histogram
ui <- navbarPage(

    # Application title
    title = "Gene Phenotype Visualizations",
    
    # page 1
    tabPanel(
      "Knockout Gene - p-value by Parameter Name",
      
      titlePanel("Knockout Gene - p-value by Parameter Name"),
      
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "gene_name",
            label = "Select Knockout Gene",
            choices = genes,
          )
        ),
        
        mainPanel(
          plotOutput("gene_knockout_plot")
        )
      )
    ),
    
    # page 2
    tabPanel(
      "Phenotype - p-value by Gene",
      
      titlePanel("Phenotype - p-value by Gene"),
      
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "phenotype_name",
            label = "Select Phenotype",
            choices = phenotypes,
          )
        ),
        
        mainPanel(
          plotOutput("phenotype_knockout_plot")
        )
      )
    ),
    
    # page 3
    tabPanel(
      "Clustering of Genes with Similar Phenotype Scores",
      
      titlePanel("PCA K-Means Clusters of Genes with Similar Phenotype Scores"),
      
      # Sidebar with a slider input for number of bins 
      sidebarLayout(
        sidebarPanel(
          numericInput(
            inputId = "k_clusters",
            label = "Select Number of Clusters - k",
            min = 2,
            max = 12,
            step = 1,
            value = 2
          )
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
          plotOutput("clusterPlot")
        )
      )
    ),
    
    # page 4 
    tabPanel(
      "More about the Genes of Interest",
      
      titlePanel("p-values for Parameters of Selected Genes"),
      
      fluidRow(
        column(
          width = 6,
          plotOutput("plot_slc25a48")
        ),
        
        column(
          width = 6,
          plotOutput("plot_ppp2r2a")
        ),
        
        column(
          width = 6,
          plotOutput("plot_miga1")
        ),
        
        column(
          width = 6,
          plotOutput("plot_mbtd1")
        ),
      )
    )

)

# Define server logic 
server <- function(input, output) {
    # for pvalue by phenotype for selected gene
    output$gene_knockout_plot <- renderPlot({
      df_subsetted_gene <- subset_gene(input$gene_name)
      plot_gene_knockout(df_subsetted_gene, input$gene_name)
    })
    
    # for pvalue by gene for selected phenotype
    output$phenotype_knockout_plot <- renderPlot({
      df_subsetted_phenotype <- subset_phenotype(input$phenotype_name)
      plot_knockout_phenotype(df_subsetted_phenotype, input$phenotype_name)
    })
    
    # clustering visualization
    output$clusterPlot <- renderPlot({
      result <- pca_clustering(gpmat, input$k_clusters)
      cluster_vis(result$pca, result$km, gene_names)
    })
    
    # for the genes of interest defined in query_genes.csv 
    output$plot_slc25a48 <- renderPlot({
      gene <- "slc25a48"
      df <- subset_gene(gene)
      plot_gene_knockout(df, gene)
    })
    
    output$plot_ppp2r2a <- renderPlot({
      gene <- "ppp2r2a"
      df <- subset_gene(gene)
      plot_gene_knockout(df, gene)
    })
    
    output$plot_miga1 <- renderPlot({
      gene <- "miga1"
      df <- subset_gene(gene)
      plot_gene_knockout(df, gene)
    })
    
    output$plot_mbtd1 <- renderPlot({
      gene <- "mbtd1"
      df <- subset_gene(gene)
      plot_gene_knockout(df, gene)
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
