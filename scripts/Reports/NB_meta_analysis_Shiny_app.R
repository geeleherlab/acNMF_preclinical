library(shiny)
library(shinythemes)

source("navbarpage_srcdata.R")

# Define UI ----
ui <- navbarPage(tags$strong("Neuroblastoma GEP Network"),
           theme = shinytheme("cerulean"),
           
           tabPanel("Network",
                 
                    sidebarLayout(
                      
                      sidebarPanel(
                        helpText("Customize the network visualization"),
      
                        selectInput("var", 
                                    label = "Display nodes by:",
                                    choices = c("Dataset", "Community"),
                                    selected = "Community"
                        )
                      ),
                      
                      forceNetworkOutput("force", height = 850)
                      
                    ),
                    
                    #guided tutorial
                    actionButton("openTutorial", "Click Here For Guided Instructions")
           ),
           
           navbarMenu("More",
                      
                      tabPanel("Metadata",
                               
                               tabsetPanel(
                                 tabPanel(tags$strong("Human"), tabsetPanel(
                                   tabPanel(tags$i("Dong"), tableOutput("dong")), 
                                   tabPanel(tags$i("GOSH"), tableOutput("gosh")), 
                                   tabPanel(tags$i("Jansky"), tableOutput("jansky")),
                                   tabPanel(tags$i("PMC"), tableOutput("pmc")),
                                   tabPanel(tags$i("Verhoeven/Kameneva"), tableOutput("kameneva")) 
                                 )),
                                 tabPanel(tags$strong("Preclinical"), tabsetPanel(
                                   tabPanel(tags$i("PDX"), tableOutput("pdx")), 
                                   tabPanel(tags$i("Mouse"), tableOutput("mouse")), 
                                   tabPanel(tags$i("Cell Lines"), tableOutput("NB_cell_line")) 
                                 )),
                               )
                               
                      ),
                      
                      tabPanel("Usage Score Heatmaps",
                               
                               tabsetPanel(
                                 tabPanel(tags$strong("Human"), tabsetPanel(
                                   tabPanel(tags$i("Dong"), fluidRow(
                                     column(3, img(src = 'images/Dong_activityscores.png', align = 'left'))
                                   )), 
                                   tabPanel(tags$i("GOSH"), fluidRow(
                                     column(3, img(src = 'images/GOSH_activityscores.png', align = 'left'))
                                   )), 
                                   tabPanel(tags$i("Jansky"), fluidRow(
                                     column(3, img(src = 'images/Jansky_activityscores.png', align = 'left'))
                                   )),
                                   tabPanel(tags$i("PMC"), fluidRow(
                                     column(3, img(src = 'images/PMC_activityscores.png', align = 'left'))
                                   )),
                                   tabPanel(tags$i("Verhoeven/Kameneva"), fluidRow(
                                     column(3, img(src = 'images/Kameneva_activityscores.png', align = 'left'))
                                   ))
                                 )),
                                 tabPanel(tags$strong("Preclinical"), tabsetPanel(
                                   tabPanel(tags$i("PDX"), fluidRow(
                                     column(3, img(src = 'images/PDX_Human_activityscores.png', align = 'left'))
                                   )), 
                                   tabPanel(tags$i("Mouse"), fluidRow(
                                     column(3, img(src = 'images/mouse_activityscores.png', align = 'left'))
                                   )), 
                                   tabPanel(tags$i("Cell Lines"), fluidRow(
                                     column(3, img(src = 'images/cellline_activityscores.png', align = 'left'))
                                   )) 
                                 )),
                                 tabPanel(tags$strong("Fetal"), tabsetPanel(
                                   tabPanel(tags$i("Human"), fluidRow(
                                     column(3, img(src = 'images/Kameneva_human_fetal_activityscores.png', align = 'left'))
                                   )), 
                                   tabPanel(tags$i("Mouse"), fluidRow(
                                     column(3, img(src = 'images/Kameneva_mouse_fetal_activityscores.png', align = 'left'))
                                   )) 
                                 )),
                               )
                               
                                
                      ),
                      
                      tabPanel("Dataset Specific GEPs",

                               tabsetPanel(
                                 tabPanel(tags$strong("Human"), tabsetPanel(
                                   tabPanel(tags$i("Dong"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_10.html", "Program 10. Neuroblastoma #2"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_19.html", "Program 19. Neuroblastoma #5"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_23.html", "Program 23. Pro-B Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_27.html", "Program 27. Schwannian Stromal Cell: NB"), br(), br(),
                                                              h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_1.html", "Program 1. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_2.html", "Program 2. Neuroblastoma #1"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_3.html", "Program 3. Translation (Activity Program)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_4.html", "Program 4. Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_5.html", "Program 5. M2 Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_6.html", "Program 6. Cancer Associated Fibroblast: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_7.html", "Program 7. Cancer Associated Fibroblast: Myofibroblastic (POSTN+)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_8.html", "Program 8. Cell Cycle (G2-M)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_9.html", "Program 9. Adrenal Cortex Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_11.html", "Program 11. Neuroblastoma #3"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_12.html", "Program 12. Schwannian Stromal Cell: GNB"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_13.html", "Program 13. Erythroblasts"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_14.html", "Program 14. Plasmacytoid Dendritic Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_15.html", "Program 15. Pre-B Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_16.html", "Program 16. Plasma Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_17.html", "Program 17. Neuroblastoma #4"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_18.html", "Program 18. Endothelial: Lymphatic"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_19.html", "Program 19. Neuroblastoma: Adrenergic II (pre-neuronal like)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_20.html", "Program 20. Cytotoxic T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_21.html", "Program 21. Ganglioneuroblastoma Stromal"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_22.html", "Program 22. Intermediate Monocyte"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_24.html", "Program 24. Neuroblastoma: Adrenergic I (sympathoblast like)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_25.html", "Program 25. Cancer Associated Fibroblast: Intermediate [Myo:Inf]"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_26.html", "Program 26. Naive B Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_28.html", "Program 28. Cell Cycle (G1/S)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_29.html", "Program 29. Myelocyte"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_30.html", "Program 30. Cancer Associated Fibroblast: Intermediate [Inf:AP]"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_31.html", "Program 31. Neuroblastoma: Adrenergic III"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_32.html", "Program 32. Regulatory T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_33.html", "Program 33. Neuroblastoma: MYCN"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_34.html", "Program 34. Cancer Associated Fibroblast: Myofibroblastic (Contractile)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_35.html", "Program 35. Naive T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_36.html", "Program 36. Neuroblastoma #5"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_37.html", "Program 37. Neuroblastoma #6"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_38.html", "Program 38. M1 Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_39.html", "Program 39. M-MDSC"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_40.html", "Program 40. Cell Cycle (G1)"), br(),
                                                              tags$a(href="GEP_reports/human/GSE137804/gepSummaryKnitr_41.html", "Program 41. IFN Response (Activity Program)"), br()
                                            
                                            ),
                                   tabPanel(tags$i("GOSH"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_12.html", "Program 12. GMP"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_15.html", "Program 15. Megakaryocyte"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_20.html", "Program 20. Keratin Expression"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_23.html", "Program 23. Megakaryocyte-erythroid progenitors"), br(), br(),
                                                            h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_1.html", "Program 1. Naive B Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_2.html", "Program 2. Inflammatory:cytokine signalling I"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_3.html", "Program 3. Normal Adrenal Fibroblasts"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_4.html", "Program 4. T Cell Activation"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_5.html", "Program 5. Neuroblastoma: Adrenergic I (sympathoblast like)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_6.html", "Program 6. Hepatocyte-like Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_7.html", "Program 7. Endothelial: Vascular"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_8.html", "Program 8. Monocyte/M-MDSC"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_9.html", "Program 9. Adrenocortical Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_10.html", "Program 10. Erythroblasts"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_11.html", "Program 11. Plasma Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_13.html", "Program 13. Cell Cycle (G2-M)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_14.html", "Program 14. Mature NKT Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_16.html", "Program 16. Mononuclear Phagocyte"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_17.html", "Program 17. Undifferentiated lymphocyte"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_18.html", "Program 18. Antigen Presentation"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_19.html", "Program 19. Cancer Associated Fibroblast: Myofibroblastic (IGF1+)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_21.html", "Program 21. Mast Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_22.html", "Program 22. Neuroblastoma: MYCN"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_24.html", "Program 24. M2 Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_25.html", "Program 25. Cancer Associated Fibroblast: Myofibroblastic (POSTN+)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_26.html", "Program 26. Pro-B Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_27.html", "Program 27. Endothelial: Lymphatic"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_28.html", "Program 28. TNFa signalling via NFkB"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_29.html", "Program 29. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_30.html", "Program 30. Plasmacytoid Dendritic Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_31.html", "Program 31. Myelocyte"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_32.html", "Program 32. Regulatory T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_33.html", "Program 33. Monocyte: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_34.html", "Program 34. Natural killer cell (CD16+)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_35.html", "Program 35. Cancer Associated Fibroblast: Myofibroblastic (Contractile)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_36.html", "Program 36. Endothelial: Angiogenesis"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_37.html", "Program 37. Naive T Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_38.html", "Program 38. CD8+ T Cells"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_39.html", "Program 39. Neuroblastoma: Adrenergic II (pre-neuronal like)"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_40.html", "Program 40. Fibroblast: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/human/GOSH/gepSummaryKnitr_41.html", "Program 41. Inflammatory: cytokine signalling II"), br()
                                            ),
                                   tabPanel(tags$i("Jansky"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_1.html", "Program 1. Neuroblastoma #1"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_7.html", "Program 7. Neuroblastoma #3"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_13.html", "Program 13. Neuroblastoma #5"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_14.html", "Program 14. Unknown NFkB"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_16.html", "Program 16. Fetal Liver Hepatoblast-Like"), br(), br(),
                                                              h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_2.html", "Program 2. Neuroblastoma #2"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_3.html", "Program 3. Cancer Associated Fibroblast: Myofibroblastic (IGF1+)"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_4.html", "Program 4. Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_5.html", "Program 5. Cancer Associated Fibroblast: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_6.html", "Program 6. Cancer Associated Fibroblast: Intermediate [Myo:Inf]"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_8.html", "Program 8. Fibroblast/Hepatoblast-Like"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_9.html", "Program 9. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_10.html", "Program 10. Neuroblastoma: Adrenergic #1"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_11.html", "Program 11. Neuroblastoma #4"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_12.html", "Program 12. Neuroblastoma: Adrenergic #2"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_15.html", "Program 15. Neuroblastoma #6"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_17.html", "Program 17. Cancer Associated Fibroblast: Myofibroblastic (POSTN+)"), br(),
                                                              tags$a(href="GEP_reports/human/Jansky/gepSummaryKnitr_18.html", "Program 18. Neuroblastoma: MYCN"), br()
                                            ),
                                   tabPanel(tags$i("PMC"), h3("Dataset Specific GEPs"), br(), h5("These GEPs appear in the network"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_1.html", "Program 1. Neuroblastoma: MYCN"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_2.html", "Program 2. Myelocyte: Unknown"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_3.html", "Program 3. Schwannain Stromal Cell"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_4.html", "Program 4. M2 Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_5.html", "Program 5. Mast Cells"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_6.html", "Program 6. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_7.html", "Program 7. Cancer Associated Fibroblast: Intermediate [Myo:Inf] (IGF1+)"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_8.html", "Program 8. Erythroblasts"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_9.html", "Program 9. Neuroblastoma: Adrenergic I (sympathoblast like)"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_10.html", "Program 10. Cytotoxic T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_11.html", "Program 11. Cancer Associated Fibroblast: Myofibroblastic (POSTN+)"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_12.html", "Program 12. Cancer Associated Fibroblast: Myofibroblastic"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_13.html", "Program 13. Intermediate Monocyte"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_14.html", "Program 14. Neuroblastoma: Adrenergic II (pre-neuronal like)"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_15.html", "Program 15. Cancer Associated Fibroblast"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_16.html", "Program 16. M-MDSC"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_17.html", "Program 17. Cancer Associated Fibroblast: Myofibroblastic (Contractile)"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_18.html", "Program 18. Cancer Associated Fibroblast: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/human/PMC/gepSummaryKnitr_19.html", "Program 19. Regulatory T Cell"), br()
                                            ),
                                   #230817
                                   tabPanel(tags$i("Verhoeven/Kameneva"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_10.html", "Program 10. Neutrophil"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_14.html", "Program 14. Fetal Kidney II (Podocytes)"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_20.html", "Program 20. Hepatocyte-Like Cell II"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_24.html", "Program 24. Fetal Kidney III"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_29.html", "Program 29. Fetal Kidney IV"), br(), br(),
                                                              h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_1.html", "Program 1. Plasma Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_2.html", "Program 2. Neuroblastoma: Adrenergic I"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_3.html", "Program 3. Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_4.html", "Program 4. Fibroblast"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_5.html", "Program 5. T Cell: Naive Helper T"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_6.html", "Program 6. Fetal Kidney I"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_7.html", "Program 7. Cancer Associated Fibroblast: Intermediate [Myo:Inf]"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_8.html", "Program 8. Hepatocyte-Like Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_9.html", "Program 9. Endothelial I"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_11.html", "Program 11. T Cell: Effector T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_12.html", "Program 12. Adrenocortical Cells"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_13.html", "Program 13. Monocyte"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_15.html", "Program 15. Schwanninan Stromal Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_16.html", "Program 16. Macrophage II"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_17.html", "Program 17. Erythroblasts"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_18.html", "Program 18. T Cells: CD8+"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_19.html", "Program 19. Cancer Associated Fibroblast: Myofibroblastic"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_21.html", "Program 21. Stress Response"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_22.html", "Program 22. Mast Cells"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_23.html", "Program 23. M-MDSC"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_25.html", "Program 25. Neuroblastoma: Adrenergic II"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_26.html", "Program 26. Cell Cycle"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_27.html", "Program 27. Plasmacytoid Dendritic Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_28.html", "Program 28. Cancer Associated Fibroblast: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_30.html", "Program 30. Natural Killer T Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_31.html", "Program 31. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_32.html", "Program 32. Naive B Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_33.html", "Program 33. Neuroblastoma: Adrenergic III"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_34.html", "Program 34. Neuroblastoma: Unknown"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_35.html", "Program 35. Neuroblastoma: Adrenergic IV"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_36.html", "Program 36. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_37.html", "Program 37. Cancer Associated Fibroblast: Myofibroblastic (TAGLN+)"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_38.html", "Program 38. Neuroblastoma: Adrenergic V"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_39.html", "Program 39. Schwannian Stromal Cell II"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_40.html", "Program 40. T Cell: Cytotoxic"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_41.html", "Program 41. Inflammatory Macrophage"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_42.html", "Program 42. Dendritic Cell"), br(),
                                                              tags$a(href="GEP_reports/human/Kameneva/gepSummaryKnitr_43.html", "Program 43. T Cell: Regulatory"), br()
                                            
                                   )
                                 )),
                                 tabPanel(tags$strong("Preclinical"), tabsetPanel(
                                   tabPanel(tags$i("PDX"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_4.html", "Program 4. PDX Human #4"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_6.html", "Program 6. PDX Human #6"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_7.html", "Program 7. PDX Human #7"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_8.html", "Program 8. PDX Human #8"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_9.html", "Program 9. PDX Human #9"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_11.html", "Program 11. PDX Human #11"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_15.html", "Program 15. PDX Human #15"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_20.html", "Program 20. PDX Human #20"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_1.html", "Program 1. PDX Mouse #1"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_5.html", "Program 5. PDX Mouse #5"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_8.html", "Program 8. PDX Mouse #8"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_9.html", "Program 9. PDX Mouse #9"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_10.html", "Program 10. PDX Mouse #10"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_13.html", "Program 13. PDX Mouse #13"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_14.html", "Program 14. PDX Mouse #14"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_15.html", "Program 15. PDX Mouse #15"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_20.html", "Program 20. PDX Mouse #20"), br(), br(),
                                                              h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_1.html", "Program 1. PDX Human #1"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_2.html", "Program 2. PDX Human #2"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_3.html", "Program 3. PDX Human #3"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_5.html", "Program 5. PDX Human #5"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_10.html", "Program 10. PDX Human #10"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_12.html", "Program 12. PDX Human #12"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_13.html", "Program 13. PDX Human #13"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_14.html", "Program 14. PDX Human #14"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_16.html", "Program 16. PDX Human #16"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_17.html", "Program 17. PDX Human #17"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_18.html", "Program 18. PDX Human #18"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_19.html", "Program 19. PDX Human #19"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_21.html", "Program 21. PDX Human #21"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Human/gepSummaryKnitr_22.html", "Program 22. PDX Human #22"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_2.html", "Program 2. PDX Mouse #2"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_3.html", "Program 3. PDX Mouse #3"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_4.html", "Program 4. PDX Mouse #4"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_6.html", "Program 6. PDX Mouse #6"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_7.html", "Program 7. PDX Mouse #7"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_11.html", "Program 11. PDX Mouse #11"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_12.html", "Program 12. PDX Mouse #12"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_16.html", "Program 16. PDX Mouse #16"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_17.html", "Program 17. PDX Mouse #17"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_18.html", "Program 18. PDX Mouse #18"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_19.html", "Program 19. PDX Mouse #19"), br(),
                                                              tags$a(href="GEP_reports/preclinical/PDX_Mouse/gepSummaryKnitr_21.html", "Program 21. PDX Mouse #21"), br()
                                            ),
                                   tabPanel(tags$i("Mouse"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_2.html", "Program 2. Neuroblastoma drug-EMT (NB847)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_3.html", "Program 3. Neuroblastoma EMT (NB837:NB847)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_12.html", "Program 12. Neuroblastoma EMT II (NB837)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_13.html", "Program 13. Neuroblastoma drug-EMT II (NB847)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_15.html", "Program 15. Neuroblastoma (NB849)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_18.html", "Program 18. CD4 T Cell"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_19.html", "Program 19. Dendritic Cell"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_23.html", "Program 23. Monocyte"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_25.html", "Program 25. Unclear program"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_31.html", "Program 31. Neuroblastoma EMT III (NB837:NB847)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_33.html", "Program 33. Neuroblastoma EMT IV (NB847)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_34.html", "Program 34. Neuroblastoma Broad"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_36.html", "Program 36. Neuroblastoma drug-EMT III (NB847)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_39.html", "Program 39. Neuroblastoma drug-EMT IV (NB853)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_42.html", "Program 42. Neuroblastoma (NB831, NB855, NB864)"), br(), br(),
                                                              h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_1.html", "Program 1. Antigen Presentation"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_4.html", "Program 4. Cancer Associated Fibroblast: Intermediate [Myo:Inf]"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_5.html", "Program 5. Neuroblastoma drug (NB839)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_6.html", "Program 6. Immune unclear"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_7.html", "Program 7. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_8.html", "Program 8. Macrophage"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_9.html", "Program 9. Cancer Associated Fibroblast: Inflammatory"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_10.html", "Program 10. Schwannian Stromal Cell"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_11.html", "Program 11. Immune unclear II"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_14.html", "Program 14. Fetal Liver Hepatoblast-Like"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_16.html", "Program 16. NFkB: Myeloid"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_17.html", "Program 17. Neuroblastoma: MYCN"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_20.html", "Program 20. Natural Killer T Cell"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_21.html", "Program 21. Cancer Associated Fibroblast: Myofibroblastic (Contractile)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_22.html", "Program 22. NFkB: Activity Program"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_24.html", "Program 24. Cancer Associated Fibroblast: Intermediate [Inf:AP]"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_26.html", "Program 26. Neuroblastoma (NB849:NB887)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_27.html", "Program 27. Putative Neuroblastoma"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_28.html", "Program 28. M-MDSC"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_29.html", "Program 29. IFN Response"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_30.html", "Program 30. Neuroblastoma: Adrenergic II"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_32.html", "Program 32. NFkB: Fibrotic Neuroblastoma"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_35.html", "Program 35. CD8+ T Cells"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_37.html", "Program 37. Myofibroblastic Cancer Associated Fibroblast"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_38.html", "Program 38. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_40.html", "Program 40. Neuroblastoma #11"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_41.html", "Program 41. Endothelial: Lymphatic"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_43.html", "Program 43. Cell Cycle (G2M)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_44.html", "Program 44. Neuroblastoma: Adrenergic I"), br(),
                                                              tags$a(href="GEP_reports/preclinical/mouse/gepSummaryKnitr_45.html", "Program 45. Plasma Cells"), br()
                                            ),
                                   tabPanel(tags$i("Cell Lines"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_1.html", "Program 1. NB Cell Line #1 (SKNFI #1)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_3.html", "Program 3. NB Cell Line #3 (BE2C)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_4.html", "Program 4. NB Cell Line #4 (IMR5)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_6.html", "Program 6. NB Cell Line #6 (NGP)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_9.html", "Program 9. NB Cell Line #9 (SKNSH Mesenchymal)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_11.html", "Program 11. NB Cell Line #11 (CHP212)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_17.html", "Program 17. NB Cell Line #17 (SKNSH Mesenchymal)"), br(), br(),
                                                              h5("These GEPs appear in the network"), br(), 
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_2.html", "Program 2. NB Cell Line #2 (Mixed)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_5.html", "Program 5. NB Cell Line #5 (Mixed)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_7.html", "Program 7. NB Cell Line #7 (Mixed Mesenchymal)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_8.html", "Program 8. NB Cell Line #8 (SKNAS Mesenchymal)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_10.html", "Program 10. NB Cell Line #10 (Mixed)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_12.html", "Program 12. NB Cell Line #12 (Mixed Adrenergic)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_13.html", "Program 13. NB Cell Line #13 (SKNFI #2)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_14.html", "Program 14. NB Cell Line #14 (Mixed Adrenergic)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_15.html", "Program 15. NB Cell Line #15 (Mixed)"), br(),
                                                              tags$a(href="GEP_reports/preclinical/NB_cell_line/gepSummaryKnitr_16.html", "Program 16. NB Cell Line #16 (Mixed)"), br()
                                            )
                                 )),
                                 tabPanel(tags$strong("Fetal"), tabsetPanel(
                                   tabPanel(tags$i("Human"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_1.html", "Program 1. Unknown I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_2.html", "Program 2. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_3.html", "Program 3. Hepatocytes"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_4.html", "Program 4. Myeloid"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_5.html", "Program 5. Unknown II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_6.html", "Program 6. Erythroblasts"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_7.html", "Program 7. Peripheral Nervous System (PNS) I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_8.html", "Program 8. Unknown III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_9.html", "Program 9. Stellate Cell"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_10.html", "Program 10. Lymphocyte"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_11.html", "Program 11. Kidney"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_12.html", "Program 12. Melanocyte"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_13.html", "Program 13. Lymphocyte II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_14.html", "Program 14. Unknown IV"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_15.html", "Program 15. PNS II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_16.html", "Program 16. Unknown V"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_17.html", "Program 17. Stress Response"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_18.html", "Program 18. Megakaryocyte"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_19.html", "Program 19. Schwann Cell"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_20.html", "Program 20. Myeloid II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_21.html", "Program 21. Lymphocyte III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_22.html", "Program 22. Unknown VI"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_23.html", "Program 23. Lymphocyte IV"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_24.html", "Program 24. PNS III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_25.html", "Program 25. Stromal I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_26.html", "Program 26. Endothelial II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_27.html", "Program 27. Cell Cycle"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_28.html", "Program 28. Stromal II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_29.html", "Program 29. Endothelial III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_30.html", "Program 30. Lymphocyte V"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_31.html", "Program 31. Erythroblasts II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_32.html", "Program 32. Adrenal Cortex I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_33.html", "Program 33. Endothelial IV"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_34.html", "Program 34. Adrenal Cortex II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_35.html", "Program 35. Unknown VII"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_36.html", "Program 36. Lymphocyte VI"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_37.html", "Program 37. Unknown VII"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_38.html", "Program 38. Stromal III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_39.html", "Program 39. Endothelial V"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_40.html", "Program 40. Unknown IX"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_41.html", "Program 41. Myeloid III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_human/gepSummaryKnitr_42.html", "Program 42. Endothelial VI"), br()
                                        ),
                                   tabPanel(tags$i("Mouse"), h3("Dataset Specific GEPs"), br(), h5("These GEPs don't appear in the network"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_1.html", "Program 1. Cardiomyocytes"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_2.html", "Program 2. Fibroblasts"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_3.html", "Program 3. Peripheral Nervous System (PNS) I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_4.html", "Program 4. Endothelial"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_5.html", "Program 5. Stromal I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_6.html", "Program 6. Myeloid"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_7.html", "Program 7. PNS II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_8.html", "Program 8. Epithelial Cells I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_9.html", "Program 9. Adrenocortical Cells"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_10.html", "Program 10. Erythroblasts"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_11.html", "Program 11. Hepatocytes"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_12.html", "Program 12. Cell Cycle"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_13.html", "Program 13. Schwann Cell I"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_14.html", "Program 14. Myofibroblasts"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_15.html", "Program 15. PNS III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_16.html", "Program 16. Myocytes"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_17.html", "Program 17. PNS III"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_18.html", "Program 18. Schwann Cell II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_19.html", "Program 19. PNS IV"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_20.html", "Program 20. Stromal II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_21.html", "Program 21. Myocytes II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_22.html", "Program 22. Lymphocytes"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_23.html", "Program 23. Epithelial Cell II"), br(),
                                                              tags$a(href="GEP_reports/fetal/Kameneva_mouse/gepSummaryKnitr_24.html", "Program 24. PNS V"), br()
                                   )
                               )),
                            )   

                      )
                      
          )
)

#makes nodes interactive with GEP html reports
Myclickaction <- "window.open(d.hyperlink)"

# Define server logic ----
server <- function(input, output) {

  output$force <- renderForceNetwork({
    node.groups <- switch(input$var, 
                          "Dataset" = db.vertex.name,
                          "Community" = community.group)
  
    node.color <- switch(input$var,
                         "Dataset" = 'd3.scaleOrdinal()
                                          .domain(["GOSH", "PMC", "Dong", "Jansky", "Verhoeven/Kameneva", "NB_cellline", "Mouse", "PDX_Human", "PDX_Mouse"])
                                          .range(["#C6DBEF", "#9ECAE1", "#6BAED6", "#3182BD", "#08519C", "#CC0000", "#FDCC0D", "#336600", "#00CC99"]);',
                         "Community" = 'd3.scaleOrdinal()
                                          .domain(["Neuroblastoma: Adrenergic", "Neuroblastoma: MYCN", "Neuroblastoma: ALK", "Neuroblastoma: Unknown #1", "Neuroblastoma: Unknown #2", "Cell Cycle", "Schwannian Stromal Cells", "Adrenal Cortex Cells", "Plasma Cells", "T Cells", "B Cells", "Mast Cells", "Plasmacytoid Dendritic Cells", "Myeloid", "Endothelial", "Cancer Associated Fibroblast (Myofibroblastic)", "Cancer Associated Fibroblast (Other)",  "Erythroblast", "Hypoxia"])
                                          .range(["#FF0000B3", "#FF5100B3", "#FFA100B3", "#FFF200B3", "#BCFF00B3", "#6BFF00B3", "#1BFF00B3", "#00FF36B3", "#00FF86B3", "#00FFD7B3", "#00D7FFB3", "#0086FFB3", "#0036FFB3", "#1B00FFB3", "#6B00FFB3", "#BC00FFB3", "#FF00F2B3", "#FF00A1B3", "#FF0051B3"]);')
    
                                         
    meta_d3 <- igraph_to_networkD3(meta.graph, group = node.groups)
    meta_d3$nodes$degree <- 3 * meta.degree
    
    meta_d3$links$value <- meta_d3$links$value + 0.000001
    
    fn <- forceNetwork(Links = meta_d3$links, 
                       Nodes = meta_d3$nodes, 
                       Source = 'source',
                       Target = 'target',
                       NodeID = 'name',
                       Group = 'group',
                       Value = 'value', 
                       Nodesize = 'degree',
                       radiusCalculation = JS("Math.sqrt(d.nodesize)*3"),
                       opacity = 1,
                       zoom = T,
                       fontSize = 20,
                       fontFamily = "Helvetica",
                       charge = -12,
                       legend = T,
                       linkDistance = JS("function(d){return (1/d.value)/1e+06}"),
                       linkWidth = 0.7,
                       linkColour = "grey",
                       bounded = T,
                       colourScale = JS(node.color),
                       clickAction = Myclickaction)
    fn$x$nodes$hyperlink <- hyperlink
    fn
    
    })
  
    #TABSET PANEL 
    #Metadata Tables
    output$gosh <- renderTable({
      gosh.table
    }) 
  
    output$pmc <- renderTable({
      pmc.table
    })
    
    output$dong <- renderTable({
      dong.table
    }) 
    
    output$jansky <- renderTable({
      jansky.table
    })
    
    output$kameneva <- renderTable({
      kameneva.table
    })
    
    output$mouse <- renderTable({
      mouse.table
    }) 
  
    output$pdx <- renderTable({
      pdx.table
    })
    
    output$NB_cell_line <- renderTable({
      cellline.table
    })
    
    #Guided Tutorial
    observeEvent(input$openTutorial, {
      showModal(
        modalDialog(
          title = "Welcome to the NB GEP Network!",
          fade = T, easyClose = T,
          tags$ul(
            "This tutorial will walk you through how to use this app.", tags$p(""),
            "You can close any time by clicking outside the popup box."), tags$p(""),
          footer = actionButton("first", "Proceed")
          )
        )
    })
    
    observeEvent(input$first, {
      showModal(
        modalDialog(
          title = "The Home Screen",
          fade = T, easyClose = T,
          "Displayed on the home screen is the NB GEP Network described in Figure 4.", tags$p(""),
          HTML('<img src="images/Tutorial/homescreen.png" style="width:540px;height:300px;" />'), br(),
          tags$p(strong("Features include:")),
          tags$ul(
            tags$li("An interactive network of GEPs (nodes) that highlight the landscape of neuroblastoma gene expression signatures."),
            tags$li("Multiple displays of said network"),
            tags$li("Additional metadata that describes the patient/sample cohort for each dataset")
          ),
          "Each of these features will be demonstrated in detail in upcoming steps.", tags$p(""),
          footer = actionButton("second", "Next")
          )
        )
    })
    
    observeEvent(input$second, {
      showModal(
        modalDialog(
           title = "Changing the Display: Node Color", 
           fade = T, easyClose = T,
          "The Display Toggle Switch changes the color of the nodes", tags$p(""),
          "You can view nodes by community or by their corresponding dataset.",br(),
          HTML('<video width="540" height ="480" controls> <source src="images/Tutorial/dataset.mp4" type="video/mp4" </video>'),
          footer = actionButton("third", "Next")
          
        )
      )
    })
    
    observeEvent(input$third, {
      showModal(
        modalDialog(
          title = "Changing the Display: Layout", 
          fade = T, easyClose = T,
          "The layout of the network can be changed by the user.", tags$p(""),
          "Simply click and hold a node and drag it to new location.", tags$p(""),
          HTML('<video width="540" height ="480" controls> <source src="images/Tutorial/drag_drop.mp4" type="video/mp4" </video>'),
          footer = actionButton("fourth", "Next")
        )
      )
    })
    
    observeEvent(input$fourth, {
      showModal(
        modalDialog(
          title = "Interacting with the Nodes", 
          fade = T, easyClose = T,
          "Hovering over each node displays the name and number of the corresponding GEP.", tags$p(""),
          "Additionally, all nodes connected to that node are highlighted.", tags$p(""),
          "Clicking a node opens the report for that GEP.", tags$p(""),
          HTML('<video width="540" height ="480" controls> <source src="images/Tutorial/clickable.mp4" type="video/mp4" </video>'),
          footer = actionButton("fifth", "Next")
        )
      )
    })
    
    observeEvent(input$fifth, {
      showModal(
        modalDialog(
          title = "Interacting with the GEP Report", 
          fade = T, easyClose = T,
          "The GEP is arguably the most important aspect of this app.", tags$p(""),
          "The GEP report is the portal that contains all of the post-acNMF computational analyses.", tags$p(""),
          tags$p(strong("Features include:")),
          tags$ul(
            tags$li("Interactive plots."),
            tags$li("Links to additional pertinent information"),
            tags$li("Downloadable tables and graphs"),
            tags$li("The ability to start dialogue about the interpretation of the GEP")
          ),
          tags$p(strong("Analyses include (but are not limited to):")),
          tags$ul(
            tags$li("InferCNV displayed as interactive BioCircos plots."),
            tags$li("UMAPs of activity scores"),
            tags$li("QQ plots for inter- and intra-dataset similarity"),
            tags$li("Enrichment of gene sets from MSigDB, literature-curated, Descartes, etc.")
          ),
          HTML('<video width="540" height ="480" controls> <source src="images/Tutorial/GEP.mp4" type="video/mp4" </video>'),
          footer = actionButton("sixth", "Next")
        )
      )
    })
    
    observeEvent(input$sixth, {
      showModal(
        modalDialog(
          title = "Additional Metadata", 
          fade = T, easyClose = T,
          "The 'More' tab includes additional metadata from each study.", tags$p(""),
          tags$p(strong("The menu includes:")),
          tags$ul(
            tags$li("Metadata - Patient/Sample metrics."),
            tags$li("Usage Scores"),
            tags$li("Dataset Specific GEPs - Includes whether they are connected to other nodes in the network."),
          ),
          HTML('<video width="540" height ="480" controls> <source src="images/Tutorial/more_tab.mp4" type="video/mp4" </video>'),
          footer = actionButton("seventh", "Next")
        )
      )
    })
    
    observeEvent(input$seventh, {
      showModal(
        modalDialog(
          title = "Thanks for Watching!",
          fade = T, easyClose = T,
          tags$ul(
            "We hope you find this app useful.", tags$p(""),
            "Feel free to add comments."), tags$p(""),
        )
      )
    })
    
}


# Run the app ----
shinyApp(ui = ui, server = server)
