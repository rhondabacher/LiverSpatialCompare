options(shiny.maxRequestSize = 10000*1024^2) 

  
ui <- fluidPage(

    headerPanel("Vizualize UMI and Full-length Liver Data Correlations"),

    fluidRow(
        column(width = 7, 
            tags$div(tags$h4("Choose a Kegg or GO category from the drop-down list below. 
            Also select the adjusted p-value significant threshold for genes to consider. This is applied to both datasets.")),
            tags$br(),
            tags$div(tags$h4("Only 10 significantly zonated genes in the Kegg category are shown. The ten are the best correlated across datasets.")),
            tags$br()     
        )
    ),
    fluidRow(
        column(width = 7, 
          uiOutput("choose_category")),
        column(width = 3, 
           numericInput("sigThresh", "Significance Threshold:", .1, min = 0, max = 1, step= .1))
    ),

    fluidRow(
        column(10, align='right',
                downloadButton('downloadPlot', 'Download Plot'))
        ),
		fluidRow(	
        column(10, align="center",
            mainPanel(plotOutput('genePlot'), width = "100%"),
            tags$br(),
            tags$br()
        )
    )

)    

    