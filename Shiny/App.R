library(shiny)
library(tidyverse)
library(ggforce)
#library(shinycssloaders)
library(HiMC); data(nodes)

#setwd('~/Dropbox/src/MitoImpute/shiny')
#setwd('~/GitCode/MitoImpute/')

MT_haps <- readRDS("MT_haps.rds")

## Reference Panel with mtSNPS MAF == 1% 
imp.info01 <- readRDS("imp.info01.rds")
imp.dat01 <- readRDS("imp.dat.rds")

## Reference Panel with mtSNPS MAF == 0.5% 
imp.info005 <- readRDS("imp.info005.rds")
imp.dat005 <- readRDS("imp.dat005.rds")

## Reference Panel with mtSNPS MAF == 0.1% 
imp.info001 <- readRDS("imp.info001.rds")
imp.dat001 <- readRDS("imp.dat001.rds")

# Define UI for app that draws a histogram ----
ui <- fluidPage(
  
  # App title ----
  titlePanel("mitoImpute Validation"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      
      selectInput("select", label = h3("Select Microarray"), 
                  choices = list('BDCHP-1X10-HUMANHAP240S_11216501_A-b37' = 'BDCHP-1X10-HUMANHAP240S_11216501_A-b37',
                                 'BDCHP-1X10-HUMANHAP550_11218540_C-b37' = 'BDCHP-1X10-HUMANHAP550_11218540_C-b37',
                                 'Cardio-Metabo_Chip_11395247_A-b37' = 'Cardio-Metabo_Chip_11395247_A-b37',
                                 'cardio-metabo_chip_11395247_c-b37' = 'cardio-metabo_chip_11395247_c-b37',
                                 'Consortium-OncoArray_15047405_A-b37' = 'Consortium-OncoArray_15047405_A-b37',
                                 'GenomeWideSNP_6.na33.annot-b37' = 'GenomeWideSNP_6.na33.annot-b37',
                                 'GSA-24v1-0_A1-b37' = 'GSA-24v1-0_A1-b37',
                                 'GSA-24v1-0_A2-b37' = 'GSA-24v1-0_A2-b37',
                                 'GSA-24v1-0_A6-b37' = 'GSA-24v1-0_A6-b37',
                                 'GSA-24v2-0_A1-b37' = 'GSA-24v2-0_A1-b37',
                                 'GSAMD-24v1-0_20011747_A1-b37' = 'GSAMD-24v1-0_20011747_A1-b37',
                                 'GSAMD-24v2-0_20024620_A1-b37' = 'GSAMD-24v2-0_20024620_A1-b37',
                                 'GSCA-24v1-0_20010441_A1-b37' = 'GSCA-24v1-0_20010441_A1-b37',
                                 'Human1-2M-DuoCustom_v1_A-b37' = 'Human1-2M-DuoCustom_v1_A-b37',
                                 'Human1M-Duov3_B-b37' = 'Human1M-Duov3_B-b37',
                                 'Human1M-Duov3_C-b37' = 'Human1M-Duov3_C-b37',
                                 'human1m-duov3_h-b37' = 'human1m-duov3_h-b37',
                                 'Human1Mv1_C-b37' = 'Human1Mv1_C-b37',
                                 'Human610-Quadv1_B-b37' = 'Human610-Quadv1_B-b37',
                                 'human610-quadv1_h-b37' = 'human610-quadv1_h-b37',
                                 'Human660W-Quad_v1_A-b37' = 'Human660W-Quad_v1_A-b37',
                                 'Human660W-Quad_v1_C-b37' = 'Human660W-Quad_v1_C-b37',
                                 'human660w-quad_v1_h-b37' = 'human660w-quad_v1_h-b37',
                                 'Human670-QuadCustom_v1_A-b37' = 'Human670-QuadCustom_v1_A-b37',
                                 'HumanCore-12-v1-0-B-b37' = 'HumanCore-12-v1-0-B-b37',
                                 'humancore-12v1-0_a-b37' = 'humancore-12v1-0_a-b37',
                                 'humancore-24-v1-0-manifest-file-a-b37' = 'humancore-24-v1-0-manifest-file-a-b37',
                                 'HumanCoreExome-12-v1-0-D-b37' = 'HumanCoreExome-12-v1-0-D-b37',
                                 'HumanCoreExome-12-v1-1-C-b37' = 'HumanCoreExome-12-v1-1-C-b37',
                                 'HumanCoreExome-12v1-0_B-b37' = 'HumanCoreExome-12v1-0_B-b37',
                                 'humancoreexome-12v1-0_c-b37' = 'humancoreexome-12v1-0_c-b37',
                                 'humancoreexome-12v1-1_a-b37' = 'humancoreexome-12v1-1_a-b37',
                                 'HumanCoreExome-12v1-1_B-b37' = 'HumanCoreExome-12v1-1_B-b37',
                                 'HumanCoreExome-24v1-0_A-b37' = 'HumanCoreExome-24v1-0_A-b37',
                                 'HumanExome-12v1_A-b37' = 'HumanExome-12v1_A-b37',
                                 'HumanExome-12v1-1_A-b37' = 'HumanExome-12v1-1_A-b37',
                                 'HumanExome-12v1-2_A-b37' = 'HumanExome-12v1-2_A-b37',
                                 'HumanHap550-2v3_B-b37' = 'HumanHap550-2v3_B-b37',
                                 'HumanHap550v3_A-b37' = 'HumanHap550v3_A-b37',
                                 'HumanHap650Yv3_A-b37' = 'HumanHap650Yv3_A-b37',
                                 'HumanOmni1-Quad_v1-0_B-b37' = 'HumanOmni1-Quad_v1-0_B-b37',
                                 'HumanOmni1-Quad_v1-0_C-b37' = 'HumanOmni1-Quad_v1-0_C-b37',
                                 'HumanOmni1-Quad_v1-0_H-b37' = 'HumanOmni1-Quad_v1-0_H-b37',
                                 'humanomni1-quad_v1-0-multi_h-b37' = 'humanomni1-quad_v1-0-multi_h-b37',
                                 'HumanOmni1S-8v1_H-b37' = 'HumanOmni1S-8v1_H-b37',
                                 'HumanOmni2-5-8-v1-0-D-b37' = 'HumanOmni2-5-8-v1-0-D-b37',
                                 'HumanOmni2-5-8-v1-1-C-b37' = 'HumanOmni2-5-8-v1-1-C-b37',
                                 'HumanOmni2-5-8-v1-2-A-b37' = 'HumanOmni2-5-8-v1-2-A-b37',
                                 'HumanOmni2-5Exome-8-v1-1-A-b37' = 'HumanOmni2-5Exome-8-v1-1-A-b37',
                                 'HumanOmni2.5-4v1_B-b37' = 'HumanOmni2.5-4v1_B-b37',
                                 'HumanOmni2.5-4v1_D-b37' = 'HumanOmni2.5-4v1_D-b37',
                                 'HumanOmni2.5-4v1_H-b37' = 'HumanOmni2.5-4v1_H-b37',
                                 'HumanOmni2.5-8v1_A-b37' = 'HumanOmni2.5-8v1_A-b37',
                                 'HumanOmni2.5-8v1_C-b37' = 'HumanOmni2.5-8v1_C-b37',
                                 'HumanOmni2.5S-8v1_B-b37' = 'HumanOmni2.5S-8v1_B-b37',
                                 'HumanOmni25-8v1-1_A-b37' = 'HumanOmni25-8v1-1_A-b37',
                                 'HumanOmni25-8v1-2_A1-b37' = 'HumanOmni25-8v1-2_A1-b37',
                                 'HumanOmni25Exome-8v1_A-b37' = 'HumanOmni25Exome-8v1_A-b37',
                                 'humanomni25m-8v1-1_b-b37' = 'humanomni25m-8v1-1_b-b37',
                                 'HumanOmni5-4-v1-0-D-b37' = 'HumanOmni5-4-v1-0-D-b37',
                                 'HumanOmni5-4v1_B-b37' = 'HumanOmni5-4v1_B-b37',
                                 'humanomni5-4v1_c-b37' = 'humanomni5-4v1_c-b37',
                                 'HumanOmni5-4v1-1_A-b37' = 'HumanOmni5-4v1-1_A-b37',
                                 'HumanOmni5Exome-4v1_A-b37' = 'HumanOmni5Exome-4v1_A-b37',
                                 'HumanOmni5Exome-4v1-1_A-b37' = 'HumanOmni5Exome-4v1-1_A-b37',
                                 'HumanOmni5Exome-4v1-2_A-b37' = 'HumanOmni5Exome-4v1-2_A-b37',
                                 'HumanOmniExpressExome-8-v1-0-B-b37' = 'HumanOmniExpressExome-8-v1-0-B-b37',
                                 'HumanOmniExpressExome-8-v1-1-C-b37' = 'HumanOmniExpressExome-8-v1-1-C-b37',
                                 'HumanOmniExpressExome-8-v1-2-B-b37' = 'HumanOmniExpressExome-8-v1-2-B-b37',
                                 'HumanOmniExpressExome-8v1_A-b37' = 'HumanOmniExpressExome-8v1_A-b37',
                                 'humanomniexpressexome-8v1-0_a-b37' = 'humanomniexpressexome-8v1-0_a-b37',
                                 'HumanOmniExpressExome-8v1-2_A-b37' = 'HumanOmniExpressExome-8v1-2_A-b37',
                                 'HumanOmniZhongHua-8v1_B-b37' = 'HumanOmniZhongHua-8v1_B-b37',
                                 'humanomnizhonghua-8v1-1_a-b37' = 'humanomnizhonghua-8v1-1_a-b37',
                                 'HumanOmniZhongHua-8v1-2_A-b37' = 'HumanOmniZhongHua-8v1-2_A-b37',
                                 'InfiniumCore-24v1-1_A1-b37' = 'InfiniumCore-24v1-1_A1-b37',
                                 'InfiniumCoreExome-24v1-1_A-b37' = 'InfiniumCoreExome-24v1-1_A-b37',
                                 'InfiniumCoreExome-24v1-2_A1-b37' = 'InfiniumCoreExome-24v1-2_A1-b37',
                                 'InfiniumExome-24v1-1_A1-b37' = 'InfiniumExome-24v1-1_A1-b37',
                                 'InfiniumOmni2-5-8v1-3_A1-b37' = 'InfiniumOmni2-5-8v1-3_A1-b37',
                                 'InfiniumOmni2-5Exome-8v1-3_A1-b37' = 'InfiniumOmni2-5Exome-8v1-3_A1-b37',
                                 'InfiniumOmni2.5Exome-8v1-2_A-b37' = 'InfiniumOmni2.5Exome-8v1-2_A-b37',
                                 'InfiniumOmni5-4v1-2_A1-b37' = 'InfiniumOmni5-4v1-2_A1-b37',
                                 'InfiniumOmni5Exome-4v1-3_A1-b37' = 'InfiniumOmni5Exome-4v1-3_A1-b37',
                                 'InfiniumOmniExpressExome-8v1-3_A-b37' = 'InfiniumOmniExpressExome-8v1-3_A-b37',
                                 'InfiniumOmniExpressExome-8v1-4_A1-b37' = 'InfiniumOmniExpressExome-8v1-4_A1-b37',
                                 'InfiniumOmniZhongHua-8v1-3_A1-b37' = 'InfiniumOmniZhongHua-8v1-3_A1-b37',
                                 'InfiniumOmniZhongHua-8v1-3_A2-b37' = 'InfiniumOmniZhongHua-8v1-3_A2-b37',
                                 'InfiniumPsychArray-24v1-1_A1-b37' = 'InfiniumPsychArray-24v1-1_A1-b37',
                                 'InfiniumPsychArray-24v1-2_A1-b37' = 'InfiniumPsychArray-24v1-2_A1-b37',
                                 'InfiniumPsychArray-24v1-2_A2-b37' = 'InfiniumPsychArray-24v1-2_A2-b37',
                                 'Multi-EthnicGlobal_A1-b37' = 'Multi-EthnicGlobal_A1-b37',
                                 'NeuroX_15036164_A-b37' = 'NeuroX_15036164_A-b37',
                                 'OmniExpressExome-8v1-1_A-b37' = 'OmniExpressExome-8v1-1_A-b37',
                                 'OmniExpressExome-8v1-1_B-b37' = 'OmniExpressExome-8v1-1_B-b37',
                                 'OncoArray-500K_B-b37' = 'OncoArray-500K_B-b37',
                                 'OncoArray-500K-C-b37' = 'OncoArray-500K-C-b37',
                                 'PsychArray_A-b37' = 'PsychArray_A-b37',
                                 'PsychArray-B-b37' = 'PsychArray-B-b37',
                                 'PsychChip_15048346_B-b37' = 'PsychChip_15048346_B-b37',
                                 'PsychChip_v1-1_15073391_A1-b37' = 'PsychChip_v1-1_15073391_A1-b37'), 
                  selected = 'GSA-24v1-0_A2-b37'),
      
      br(),
      
      # Input: Slider for the number of bins ----
      sliderInput(inputId = "bins",
                  label = "Info Score Cut Off:",
                  min = 0,
                  max = 1,
                  value = 0.3), 
      
      br(),
      
      radioButtons("MAF", label = "Reference MAF",
                   choices = list("1%" = 1, "0.5%" = 2, "0.1%" = 3), 
                   selected = 1),
      
      br(), 
      
      checkboxInput("ShowHiMC", label = "Show Hi-MC", value = TRUE),
      
      br(),
      
      checkboxInput("SnpType", label = "SNP Type", value = TRUE)
      
      
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      
      # Output: Tabset w/ plot, summary, and table ----
      tabsetPanel(type = "tabs",
                  tabPanel("Summary", 
                           htmlOutput("Platform"),
                           htmlOutput("type0"),
                           textOutput("type2"),
                           textOutput("type3"),
                           htmlOutput("Proportion"),
                           htmlOutput("PreConc"), 
                           textOutput("PostConc")),
                  tabPanel("Info Score",  plotOutput(outputId = "distPlot"), DT::dataTableOutput("table")),
                  tabPanel("Haplogroups, Pre", plotOutput(outputId = "PrePlot"), DT::dataTableOutput("PreMatch")),
                  tabPanel("Haplogroups, Post", plotOutput(outputId = "PostPlot"), DT::dataTableOutput("PostMatch")) 
                  
      )
    )
  )
)

server <- function(input, output) {
  
  mt.haps <- reactive({
    ## Cut put to include snps
    info.cut <- input$bins
    select <- input$select
    MAF <- input$MAF
    
    imp.info <- if(MAF == 1){imp.info01}else if(MAF == 2){imp.info005}else{imp.info001}
    imp.dat <- if(MAF == 1){imp.dat01}else if(MAF == 2){imp.dat005}else{imp.dat001}
    
    ## Filter SNPs
    rm.info <- filter(imp.info[[select]], info > info.cut)
    imp.dat_filt <- imp.dat[[select]][ ,colnames(imp.dat[[select]]) %in% 
                                               c('Individual', rm.info$position)]
    
    ## Assign haplogroups
    MTimp.classifications <- HiMC::getClassifications(as.data.frame(imp.dat_filt))
    mt.haps <- MT_haps[[select]] %>% 
      left_join(MTimp.classifications, by = 'Individual') %>% 
      rename(haplogroup_imp = haplogroup, full_path_imp = full_path)
  })
  
    ## Summary Text
  output$Platform = renderUI({
    HTML(paste0(tags$br(), "Platform: ", input$select))
  })
  output$type0 = renderUI({
    imp.info <- if(input$MAF == 1){imp.info01}else if(input$MAF == 2){imp.info005}else{imp.info001}
    SnpType <- filter(imp.info[[input$select]], type == 0)
    HTML(paste0(tags$br(), "SNPs in Reference Panel only: ", nrow(SnpType)))
  })
  output$type2 = renderText({
    imp.info <- if(input$MAF == 1){imp.info01}else if(input$MAF == 2){imp.info005}else{imp.info001}
    SnpType <- filter(imp.info[[input$select]], type == 2)
    paste0("SNPs in Reference & Sample Panel only: ", nrow(SnpType))
  })
  output$type3 = renderText({
    imp.info <- if(input$MAF == 1){imp.info01}else if(input$MAF == 2){imp.info005}else{imp.info001}
    SnpType <- filter(imp.info[[input$select]], type == 3)
    paste0("SNPs in Sample Panel only: ", nrow(SnpType))
  })
  
  output$Proportion = renderUI({
    info.cut <- input$bins
    imp.info <- if(input$MAF == 1){imp.info01[[input$select]]}else if(input$MAF == 2){imp.info005[[input$select]]}else{imp.info001[[input$select]]}
    HTML(paste0(tags$br(), nrow(filter(imp.info, info > info.cut)), ' out of ', nrow(imp.info), ' mtSNPs (', 
                nrow(filter(imp.info, info > info.cut & himc == 'yes')), '/', nrow(filter(imp.info, himc == 'yes')), ' Hi-MC mtSNPs) are retained using an Info threshold of ', info.cut 
                ))
  })
  
  
  output$PreConc = renderUI({
    hap.concordance <- MT_haps[[input$select]] %>%
      count(haplogroup_typ, haplogroup_wgs) %>% 
      mutate(perc = (n/sum(n))*100) %>%
      mutate(match = haplogroup_typ == haplogroup_wgs) %>% 
      group_by(match) %>%
      summarise(sum = round(sum(perc), 2))
    
    HTML(paste0(tags$br(), "Typed vs WGS Haplogroup concordance: ", as.character(hap.concordance[2,2]), '%'))
  })
  output$PostConc = renderText({

    ## Count pairs of haplogroups of imputed and WGS assignments
    hap.concordance <- mt.haps() %>%
      count(haplogroup_imp, haplogroup_wgs) %>% 
      mutate(perc = (n/sum(n))*100) %>% 
      mutate(match = haplogroup_imp == haplogroup_wgs) %>% 
      group_by(match) %>%
      summarise(sum = round(sum(perc), 2))

    HTML(paste0("Typed + Imputed vs WGS Haplogroup concordance: ", as.character(hap.concordance[2,2]), '%'))
  })
  
  ##===============================================##  
  ## Plot displaying info score across mitochondrial genome
  ##===============================================##  
  output$distPlot <- renderPlot({
    
    info.cut <- input$bins
    select <- input$select
    MAF <- input$MAF
    ShowHiMC <- input$ShowHiMC
    
    imp.info <- if(MAF == 1){imp.info01}else if(MAF == 2){imp.info005}else{imp.info001}

    if(ShowHiMC == F){
      ggplot(imp.info[[select]], 
             aes(x = position, y = info_comb, size = exp_freq_a1,
                 colour = if(input$SnpType == F){info_comb > info.cut}
                 else{as.factor(type)})) + 
        geom_point() + theme_bw() + 
        labs(x = 'mtDNA position', y = 'Info Score', 
             title = 'Info Score of imputed mtSNPs') + 
        geom_hline(yintercept = info.cut, linetype = 2, colour = 'red') +
        guides(colour=guide_legend(title="SNP Type")) 
    }else{
      ggplot(imp.info[[select]], 
             aes(x = position, y = info_comb, size = exp_freq_a1, alpha = himc, 
                 colour = if(input$SnpType == F){info_comb > info.cut}
                 else{as.factor(type)})) + 
        geom_point() + theme_bw() + 
        labs(x = 'mtDNA position', y = 'Info Score', 
             title = 'Info Score of imputed mtSNPs') + 
        geom_hline(yintercept = info.cut, linetype = 2, colour = 'red') +
        guides(colour=guide_legend(title="SNP Type")) 
      
    }
    
  })
  
  output$table <- DT::renderDataTable(DT::datatable({
    MAF <- input$MAF
    imp.info <- if(MAF == 1){imp.info01}else if(MAF == 2){imp.info005}else{imp.info001}
    
    if(input$ShowHiMC == F){
      imp.info[[input$select]] %>% 
        select(position, a0, a1, exp_freq_a1, info, type, info_type0, himc, Haplogroup)
    }else{
      imp.info[[input$select]] %>% 
        select(position, a0, a1, exp_freq_a1, info, type, info_type0, himc, Haplogroup) %>%
        filter(himc == 'yes')
      
    }
  }))
  ##===============================================##  
  ## Haplogroup Concordance plot - Typed vs WGS
  ##===============================================##  
  output$PrePlot <- renderPlot({
    select <- input$select
    
    ## Count pairs of haplogroups of imputed and WGS assignments
    hap.match <- MT_haps[[select]] %>%
      count(haplogroup_typ, haplogroup_wgs) %>% 
      mutate(perc = round((n/sum(n))*100,2)) %>% 
      mutate(match = haplogroup_typ == haplogroup_wgs) 
    
    ## Use ggforce to tidy data for geom_parallel Sets 
    # Requires developmental version of ggforce
    dat_ggforce <- hap.match  %>%
      gather_set_data(1:2) %>%       
      arrange(x,haplogroup_wgs,desc(haplogroup_typ))
    
    ggplot(dat_ggforce, aes(x = x, id = id, split = y, value = n)) +
      geom_parallel_sets(aes(fill = match), alpha = 0.5, axis.width = 0.2) +
      geom_parallel_sets_labels(colour = 'black', angle = 0, size = 3) + 
      theme_classic() + theme(legend.position = 'bottom') + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            text = element_text(size=12)) + 
      scale_x_discrete(labels=c("Typed", "WGS")) + 
      labs(x = "Mitochondrial Haplogroups") + scale_fill_brewer(palette = 'Set1')
    
  })
  
  output$PreMatch <- DT::renderDataTable(DT::datatable({
    MT_haps[[input$select]] %>%
      count(haplogroup_typ, haplogroup_wgs) %>% 
      mutate(perc = round((n/sum(n))*100)) %>% 
      mutate(match = haplogroup_typ == haplogroup_wgs) 
  }))
  
  ##===============================================##  
  ## Haplogroup Concordance plot - Imputed vs WGS
  ##===============================================##  
  output$PostPlot <- renderPlot({
    ## Count pairs of haplogroups of imputed and WGS assignments
    hap.match <- mt.haps() %>%
      count(haplogroup_imp, haplogroup_wgs) %>% 
      mutate(perc = round((n/sum(n))*100)) %>% 
      mutate(match = haplogroup_imp == haplogroup_wgs) 
    
    ## Use ggforce to tidy data for geom_parallel Sets 
    # Requires developmental version of ggforce
    dat_ggforce <- hap.match  %>%
      gather_set_data(1:2) %>%       
      arrange(x,haplogroup_wgs,desc(haplogroup_imp))
    
    ggplot(dat_ggforce, aes(x = x, id = id, split = y, value = n)) +
      geom_parallel_sets(aes(fill = match), alpha = 0.5, axis.width = 0.2) +
      geom_parallel_sets_labels(colour = 'black', angle = 0, size = 3) + 
      theme_classic() + theme(legend.position = 'bottom') + 
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            text = element_text(size=12)) + 
      scale_x_discrete(labels=c("Imputed", "WGS")) + 
      labs(x = "Mitochondrial Haplogroups") + scale_fill_brewer(palette = 'Set1')
    
      })

  output$PostMatch <- DT::renderDataTable(DT::datatable({
     ## Count pairs of haplogroups of imputed and WGS assignments
    mt.haps() %>%
      count(haplogroup_imp, haplogroup_wgs) %>% 
      mutate(perc = round((n/sum(n))*100)) %>% 
      mutate(match = haplogroup_imp == haplogroup_wgs) 
  }))
  
}

shinyApp(ui = ui, server = server)
