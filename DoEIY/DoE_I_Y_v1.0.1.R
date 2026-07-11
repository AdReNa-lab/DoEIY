library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(rhandsontable)
library(DT)
library(ggplot2)
library(plotly)
library(broom)

# Source utility modules
source("formula_utils.R")
source("optimization_utils.R")

# Source designs
source("Box_Behnken_Designs.R")
source("Central_Composite_Designs.R")
source("D_Optimal_Design.R")
source("D_Optimal_Augment.R")
source("Fractional_Factorial_Designs.R")
source("Full_Factorial_Designs.R")
source("Latin_Hypercube_Designs.R")
source("Plackett_Burman_Designs.R")

# Source Shiny namespace modules
source("modules/make_design.R")
source("modules/enter_results.R")
source("modules/analyze_design.R")
source("modules/model_explorer.R")
source("modules/resources.R")

# ===== UI Layout =============================================================
ui <- dashboardPage(
  skin = "black",
  header = dashboardHeader(
    title = span(
      img(src = "DoEIYLogo.png", height = "40px", style = "margin-right: 15px; margin-bottom: 5px;"),
      strong("DoEIY", style = "color: white; font-size: 1.25em;")
    ),
    titleWidth = 350
  ),
  sidebar = dashboardSidebar(
    width = 350,
    sidebarMenu(
      id = "tabs",
      menuItem(text = "Make a Design", tabName = "Make_Design_id"),
      menuItem(text = "Enter/Edit Results", tabName = "Enter_edit_results_id"),
      menuItem(text = "Analyze the Design", tabName = "Analyse_id"),
      menuItem(text = "Model Explorer", tabName = "Explorer_id"),
      menuItem(text = "Resources & References", tabName = "Resources_id")
    )
  ),
  body = dashboardBody(
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "Style_sheet_v1.css")
    ),
    tabItems(
      tabItem(tabName = "Make_Design_id", make_design_ui("make_design")),
      tabItem(tabName = "Enter_edit_results_id", enter_results_ui("enter_results")),
      tabItem(tabName = "Analyse_id", analyze_design_ui("analyze_design")),
      tabItem(tabName = "Explorer_id", model_explorer_ui("model_explorer")),
      tabItem(tabName = "Resources_id", resources_ui("resources"))
    )
  ),
  footer = dashboardFooter(
    left = p(
      paste0("\U00A9", format(Sys.Date(), "%Y")),
      a(
        href = "https://www.ucl.ac.uk/responsive-nanomaterials/",
        strong("AdReNa"), target = "_blank"
      )
    ),
    right = "Version: 1.0.1"
  )
)

# ===== Server Shell ==========================================================
server <- function(input, output, session) {
  # Global shared reactive variables
  design_reactive <- reactiveVal(value = NULL)
  factor_types_reactive <- reactiveVal(value = NULL)
  factor_data_reactive <- reactiveVal(value = NULL)
  
  factor_names_reactive <- reactiveVal(value = NULL)
  factor_blocks_reactive <- reactiveVal(value = NULL)
  responses_reactive <- reactiveVal(value = NULL)
  
  design_matrix_reactive <- reactiveVal(value = NULL)
  model_formula_reactive <- reactiveVal(value = NULL)
  model_aov_reactive <- reactiveVal(value = NULL)
  
  # Temporary storage for D-Optimal model build
  d_optimal_model_components <- reactiveVal(value = NULL)
  
  # Call module servers
  make_design_server("make_design", 
                     design_reactive = design_reactive, 
                     factor_data_reactive = factor_data_reactive, 
                     factor_types_reactive = factor_types_reactive,
                     d_optimal_model_components = d_optimal_model_components)
  
  enter_results_server("enter_results", 
                       design_reactive = design_reactive, 
                       factor_types_reactive = factor_types_reactive, 
                       factor_data_reactive = factor_data_reactive)
  
  analyze_design_server("analyze_design", 
                        design_reactive = design_reactive, 
                        factor_types_reactive = factor_types_reactive, 
                        factor_names_reactive = factor_names_reactive, 
                        factor_blocks_reactive = factor_blocks_reactive, 
                        responses_reactive = responses_reactive,
                        design_matrix_reactive = design_matrix_reactive, 
                        model_formula_reactive = model_formula_reactive, 
                        model_aov_reactive = model_aov_reactive)
  
  model_explorer_server("model_explorer", 
                        design_matrix_reactive = design_matrix_reactive, 
                        factor_types_reactive = factor_types_reactive, 
                        model_formula_reactive = model_formula_reactive, 
                        model_aov_reactive = model_aov_reactive)
  
  resources_server("resources")
}

# ===== Start the App ==========================================================
shinyApp(ui, server)
