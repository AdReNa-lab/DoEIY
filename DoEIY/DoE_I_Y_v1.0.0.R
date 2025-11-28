library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinyWidgets)
library(gfonts)
library(waiter)
library(ggplot2)

library(rhandsontable)
library(data.table)
library(mltools)
library(broom)
library(plotly)
library(rlang)
library(FrF2)
library(rsm)
library(lhs)
library(DoE.base)
library(AlgDesign)
library(dplyr)


# ===== Load the Design Functions ============================================
source("Box_Behnken_Designs.R")
source("Full_Factorial_Designs.R")
source("Plackett_Burman_Designs.R")
source("Fractional_Factorial_Designs.R")
source("Central_Composite_Designs.R")
source("Latin_Hypercube_Designs.R")
source("D_Optimal_Design.R")
source("D_Optimal_Augment.R")


# ===== UI =====================================================================
ui <- dashboardPage(

  # ===== Header ===============================================================
  title = "DoEIY",
  header <- dashboardHeader(
    title = uiOutput(outputId = "app_title_id")
  ),

  # ===== Sidebar ==============================================================
  sidebar <- dashboardSidebar(
    width = 230,
    sidebarMenu(
      id = "tabs",
      menuItem("Make a Design", tabName = "Make_Design_id", icon = icon("gears")),
      menuItem("Enter/Edit Results", tabName = "Enter_edit_results_id", icon = icon("table")),
      menuItem("Analyze the Design", tabName = "Analyse_id", icon = icon("chart-line")),
      menuItem("Model Explorer", tabName = "Explorer_id", icon = icon("magnifying-glass-chart")),
      menuItem("Resources & References", tabName = "Resources_id", icon = icon("circle-info"))
    )
  ),

  # ===== Body =================================================================
  body <- dashboardBody(

    # ===== Style ==============================================================
    tags$style("@import url(https://use.fontawesome.com/releases/v6.1.1/css/all.css);"),
    tags$head(tags$style(HTML(
      "
      .skin-blue .main-header .navbar {
        background-color: #0097a9;
      }
      .skin-blue .main-header .logo:hover {
        background-color: #0097a9;
      }
      .skin-blue .main-header .logo {
        background-color: #0097a9;
      }
      .skin-blue .left-side, .skin-blue .main-sidebar, .skin-blue .wrapper {
          background-color: #fff;
      }
      .skin-blue .sidebar-menu>li.active>a {
          border-left-color: #0097a9;
      }
      .skin-blue .sidebar-menu>li.active>a, .skin-blue .sidebar-menu>li.menu-open>a, .skin-blue .sidebar-menu>li:hover>a {
          color: #333333;
          background: #fff;
      }
      .skin-blue .sidebar-menu>li.active>a, .skin-blue .sidebar-menu>li.menu-open>a {
        font-weight: bold;
      }
      .skin-blue .sidebar a {
          color: #333333;
      }
      a {
        color: #0097a9;
      }
      .main-footer {
          padding-bottom: 5px;
          padding-top: 8px;
          padding-left: 15px;
          padding-right: 15px;
      }
      hr {
        border-top: 2px solid #0097a9;
      }
      .box.box-primary {
        border-color: #0097a9;
      }
      .box-header {
        margin-top: -40px;
      }
      .box.box-solid.box-primary .box-header {
        border-radius: 3px;
      }
      .Grey-Button, .btn-default {
        height: 40px;
        border-radius: 20px;
        background-color: #ddd;
      }
      .modal-content  {
        border-radius: 10px !important;
      }
      .Blue-Button, .btn-primary {
        background-color: #0097a9;
        border-radius: 20px;
        color: #FFFFFF;
        text-align:  center;
        font-weight: bold;
        font-size: 14px;
        height: 40px;
        border: 2px solid #0097a9;
      }
      .Blue-Button:hover, .btn-primary:hover {
        background-color: #0097a9;
        color: #FFFFFF;
        border: 2px solid #0097a9;
        box-shadow: 0px 0px 3px #0097a9;
      }
      .Blue-Button:focus, .btn-primary:focus {
        background-color: #0097a9;
        color: #FFFFFF;
        border: 2px solid #0097a9;
      }
      .modal-dialog {
        width: 1200px !important;
      }
       .input-group-btn:first-child>.btn {
        height: 34px;
        background-color: #0097a9;
        color: #FFFFFF;
       }
       .handsontable-container {
        overflow: auto;
        padding-bottom: 100px;
       }
       .htDropdownMenu {
        overflow: visible !important;
       }
      "
    ))),

    # ===== Tab Items ==========================================================
    tabItems(
      tabItem(
        tabName = "Make_Design_id",
        box(
          id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
          fluidRow(
            column(
              4,
              selectizeInput(inputId = "type_of_design_select_id", label = "Select a Design:", choices = list("Screening Designs" = list("Plackett-Burman"), "Factorial" = list("Full Factorial", "Fractional Factorial"), "Response Surface Designs" = list("Box-Behnken", "Central Composite"), "Space Filling" = list("Latin Hypercube Sampling"), "Custom" = list("D-Optimal")), selected = NULL, width = "100%", options = list(
                placeholder = "Please select a suitable design",
                onInitialize = I('function() { this.setValue(""); }')
              )),
              br(),
              actionButton(inputId = "make_design_id", label = "Make Design", width = "100%", class = "Blue-Button", style = "height: 34px;")
            ),
            column(
              8,
              uiOutput(outputId = "selected_design_text_id")
            )
          )
        )
      ),
      tabItem(
        tabName = "Enter_edit_results_id",
        box(
          id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
          fluidRow(
            column(4, 
                uiOutput("augment_design_btn_ui_id", )     
            ),
            column(4,
              div(
                align = "center",
                downloadButton(outputId = "save_design_id", "Save", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px; width: 200px;")
              )
            ),
            column(
              4,
              div(
                align = "center",
                actionButton(inputId = "load_design_id", label = "Load", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;", icon = icon("upload"))
              )
            )
          ),
          hr(),
          uiOutput(outputId = "design_table_UI_id")
        )
      ),
      tabItem(
        tabName = "Analyse_id",
        box(
          id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
          tabBox(
            id = "Analysis_results_tabs_id", title = NULL, side = "left", width = "100%",
            tabPanel(
              "Create Model",
              p("Note:  Load a design for analysis on the Enter/Edit Results tab."),
              uiOutput("analysis_table_ui_id")
            ),
            tabPanel(
              "Model Fitting Results",
              uiOutput("ANOVA_output_ui_id")
            ),
            tabPanel(
              "Design Diagnostics",
              uiOutput("Design_diagnostics_output_ui_id")
            )
          )
        )
      ),
      
      tabItem(tabName = "Explorer_id", 
              box(id = NULL, title = NULL, status = "primary", solidHeader = FALSE, 
                  collapsible = FALSE, headerBorder = FALSE, width = 12,
                  tabBox(id = "Explorer_tabs_id", title = NULL, side = "left", 
                         width = "100%",
                         tabPanel(
                           "Visualizer",
                           uiOutput("model_visualizer_output_ui_id")
                         ),
                         tabPanel(
                           "Response Optimizer",
                           uiOutput("model_interaction_output_ui_id")
                         )
                  )
              )
      ),
      
      tabItem(tabName = "Resources_id", 
              box(id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
                  tabBox(id = "Resources_tabs_id", title = NULL, side = "left", 
                         width = "100%",
                         tabPanel("DoEIY Manual",
                                  tags$iframe(style="height:600px; width:100%", src="User_Guide_v1.pdf")
                         ),
                         tabPanel("References",
                                  p(HTML("
                                    <p style='margin: 0 0 16px 0;'>This app's design generators were build on foundational work by the researchers and package authors listed below.</p>
                                    <ul>
                                    <li><strong>Box-Behnken</strong><br />Box, G. E. P., &amp; Behnken, D. W. (1960). Some New Three Level Designs for the Study of Quantitative Variables. <em>Technometrics</em>, 2(4), 455-475. DOI: <a href='https://www.tandfonline.com/
                                                                               doi/abs/10.1080/00401706.1960.10489912'> 10.1080/00401706.1960.10489912</a>.</li>
                                    <li><strong>Central Composite</strong><br />Lenth, R. V. (2008). <code>rsm</code>: Response-Surface Analysis (R package). DOI: <a href='https://doi.org/10.32614/
                                                                               CRAN.package.rsm'>10.32614/CRAN.package. rsm</a>.</li>
                                    <li><strong>D-Optimal</strong><br />AlgDesign: Algorithmic Experimental Design (R package). DOI: <a href='https://doi.org/10.32614/
                                                                               CRAN.package.AlgDesign'>10.32614/ CRAN.package.AlgDesign</a>.</li>
                                    <li><strong>Fractional Factorial</strong> <br />Gr&ouml;mping, U. (2025). <code>FrF2</code>: Fractional Factorial Designs with 2-Level Factors (R package). DOI: <a href='https://doi.org/10.32614/
                                                                               CRAN.package.FrF2'>10.32614/ CRAN.package.FrF2</a>.</li>
                                    <li><strong>Latin Hypercube</strong> <br />Carnell, R. (2024). <code>lhs</code> : Latin Hypercube Samples (R package). DOI: <a href='https://doi.org/10.32614/
                                                                               CRAN.package.lhs'>10.32614/ CRAN.package.lhs</a>.</li>
                                    <li><strong>Plackett&ndash;Burman</strong> <br />Plackett, R. L., &amp; Burman, J. P. (1946). The Design of Optimum Multifactorial Experiments. <em>Biometrika</em>, 33(4), 305-325. DOI: <a href='https://
                                                                               academic.oup.com/biomet/article-abstract/33/
                                                                               4/305/225377?redirectedFrom=fulltext'>10.1093/ biomet/33.4.305</a>.</li>
                                    </ul>
                                  "
                                  ))
                         )
                  )
              )
      )
    )
  ),

  # ===== Footer ===============================================================
  footer = dashboardFooter(
    left = p(
      paste0("\U00A9", format(Sys.Date(), "%Y")),
      a(
        href = "",
        strong(""), target = "_blank"
      )
    ),
    right = "Version: 1.0.0"
  )
)

# ===== Server =================================================================
server <- function(input, output, session) {
  # ===== Header Formatting ====================================================
  observe({
    state <- input$sidebarCollapsed
    if (input$sidebarCollapsed == TRUE) {
      output$app_title_id <- renderUI({
        h4(strong(HTML("D<sub>oE</sub>IY")), align = "center", style = "color: #fff; margin-left:-13px; margin-top: 15px;")
      })
    } else {
      output$app_title_id <- renderUI({
        h3(strong("DoE-It-Yourself"), align = "center", style = "color: #fff; margin-top: 15px;")
      })
    }
  })


  # ===== Make a Design ========================================================
  # stores factors and limits/requirements for the selected design
  column_widths_reactive <- reactiveVal(value = NULL) # column widths for factor table
  max_factors_reactive <- reactiveVal(value = NULL) # max allowed factors for chosen design
  min_factors_reactive <- reactiveVal(value = NULL) # min allowed factors for chosen design
  set_num_levels_reactive <- reactiveVal(value = NULL) # required number of levels, if appropriate
  factors_string_reactive <- reactiveVal(value = NULL) # description of design
  

  factor_data_reactive <- reactiveVal(value = NULL) # Store the factor data
  d_optimal_model_components <- reactiveVal(value = NULL) # Stores the selected model components for a D-Optimal Design

  factor_types_reactive <- reactiveVal(value = NULL) # Saves the factor types (i.e. continuous, categorical, etc)

  design_reactive <- reactiveVal(value = NULL) # Stores the created design


  # Change the text describing the design based on the users' selection
  observeEvent(input$type_of_design_select_id, {
    if (input$type_of_design_select_id == "") {
      output$selected_design_text_id <- renderUI({
        p("")
      })
    }
    if (input$type_of_design_select_id == "Plackett-Burman") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            '<p style="text-align: justify;"><strong>Plackett-Burman Designs</strong> are efficient screening tools for estimating main effects in an experiment. These designs allow for independent estimation of main effects, though they may be confounded with higher-order interactions. This software offers 2-level Plackett-Burman designs suitable for experiments involving 4 to 23 factors.</p>'
          ),
          style = "margin-top: 20px;"
        )
      })
    }
    if (input$type_of_design_select_id == "Full Factorial") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            '<p style="text-align: justify;"><strong>Full Factorial Designs</strong> investigate every possible combination of factor levels, using two or more levels for each factor. This comprehensive approach enables the independent estimation of all main effects and higher-order interactions, providing a complete understanding of factor relationships.</p>'
          ),
          style = "margin-top: 20px;"
        )
      })
    }
    if (input$type_of_design_select_id == "Fractional Factorial") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            "<p style='text-align: justify;'><strong>Fractional Factorial Designs</strong> are efficient alternatives to full factorial designs, utilizing only a fraction of the total possible runs. They allow for the estimation of main effects and, depending on the design's resolution, certain interactions. This software provides 2-level fractional factorial designs with the following resolutions: </p>
          <ul>
          <li style='text-align: justify;'><strong>Resolution III</strong> - main effects are not confounded with each other but may be confounded with two-factor interactions.</li>
          <li style='text-align: justify;'><strong>Resolution IV</strong> - main effects are not confounded with other main effects or two-factor interactions; two-factor interactions may still be confounded with each other.</li>
          <li style='text-align: justify;'><strong>Resolution V</strong> - two-factor interactions are not confounded with main effects or other two-factor interactions, but may be confounded with three-factor interactions.</li>
          </ul>"
          ),
          style = "margin-top: 20px;"
        )
      })
    }
    if (input$type_of_design_select_id == "Box-Behnken") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            "<p><strong>Box-Behnken Designs</strong> are used to estimate quadratic models by employing three levels for each factor. This design enables the independent estimation of main effects, two-factor interactions, and quadratic effects. This software offers Box-Behnken designs for experiments with 3 to 7 factors.</p>"
          ),
          style = "margin-top: 20px;"
        )
      })
    }
    if (input$type_of_design_select_id == "Central Composite") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            "<p style='text-align: justify;'><strong>Central Composite Designs</strong> are used to fit quadratic models by incorporating three levels for each factor.&nbsp; The are composed of 3 parts.&nbsp;&nbsp;</p>
          <ol>
          <li style='text-align: justify;'><em>Factorial Points</em> - a full or fractional factorial design (typically Resolution V) that estimates the main effects and two factor interactions</li>
          <li style='text-align: justify;'><em>Axial Points</em> - a set of points places along each axis (between or beyond the extremes of the factors) to quadratic effects</li>
          <li style='text-align: justify;'><em>Center Points</em> - located at the midpoint of all factors; helps improve precision and estimate experimental error, while also preserving orthogonality and rotatability (which determines the number of center points)</li>
          </ol>
          <p style='text-align: justify;'>The positioning of the axial points depends on whether the design is Circumscribed, Inscribed, or Face-Centered.</p>
          <ul>
          <li style='text-align: justify;'><strong>Circumscribed </strong>- the axial points lie outside the range of the factorial design.</li>
          <li style='text-align: justify;'><strong>Inscribed </strong>- the axial points lie within the provided upper and lower bounds and the factorial design is 'shrunk' to fit within the range to preserve rotatability.&nbsp;</li>
          <li style='text-align: justify;'><strong>Face Centered </strong>- the axial points lie on the center of each face of the factorial design. These axial points do not exceed the bounds of the factorial design. Note: this option is not rotatable.</li>
          </ul>"
          ),
          style = "margin-top: 20px;"
        )
      })
    }
    if (input$type_of_design_select_id == "Latin Hypercube Sampling") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            '<p style="text-align: justify;"><strong>Latin Hypercube Sampling </strong>is a statistical method for creating designs that efficiently explore the design space using a set number of runs. It ensures that each factor is sampled uniformly across its range, allowing for a comprehensive exploration of the factor space with fewer experimental runs. It is particularly useful for exploring complex, high-dimensional spaces when other designs would be too experimentally onerous.</p>
           <p style="text-align: justify;">Note: There is some randomness involved in the LHS method, so a somewhat different design may be generated from the same set of inputs, but each design remains optimized for balanced, space-filling coverage.</p>'
          ),
          style = "margin-top: 20px;"
        )
      })
    }
    if (input$type_of_design_select_id == "D-Optimal") {
      output$selected_design_text_id <- renderUI({
        p(
          HTML(
            '<p style="text-align: justify;"><strong>D-Optimal Designs </strong>are an efficient alternative to classical designs, especially when traditional approaches require too many experimental runs or face constraints within the design space. Based on a specified model that includes the necessary effects and interactions, the D-Optimal algorithm selects a subset of experiments that provides the most precise parameter estimates. This software supports D-Optimal designs for any number of factors, levels, and factor types, including continuous, discrete, categorical, and blocking factors</p>'
          ),
          style = "margin-top: 20px;"
        )
      })
    }
  })

  # Begins the workflow for making the selected design
  # Sets the limits/requirements for the selected design
  # Gathers the factors' details
  observeEvent(input$make_design_id, {
    selected_design <- input$type_of_design_select_id

    if (selected_design == "") {
      showNotification(strong("A type of design must be selected."), duration = 15, closeButton = TRUE, type = "error")
    } else {
      if (selected_design == "Plackett-Burman") {
        # Settings for Plackett-Burman Designs
        factors_string <- paste0(strong("Plackett-Burman Designs"), " are supported for
                                4-23 factors.  These factors can be continuous or categorial.
                                For each factor, only high and low levels are considered.
                                Please provide the levels as comma separated values, e.g., '-1, 1'.")
        min_factors <- 4
        max_factors <- 23
        set_num_levels <- 2
        factor_df <- data.frame(
          Factor = rep("", min_factors),
          Type = factor(rep("Continuous", min_factors), levels = c("Continuous", "Discrete", "Categorical")),
          Level_Values = rep("-1, 1", min_factors)
        )
        colnames(factor_df) <- c("Factor Names", "Factor Type", "Level Values")
        column_widths <- c(200, 200, 700)
      } else if (selected_design == "Box-Behnken") {
        # Settings for Box-Behnken Designs
        factors_string <- paste0(strong("Box-Behnken Designs"), " are supported for
                                3-12 factors.  These factors must be continuous.
                                For each factor, three levels are considered: a high, mid, and low level.
                                Please provide the levels as comma separated values, e.g., '-1, 0, 1'.
                                If only two values are provided for the levels, the mid level will be
                                automatically calculated as the midpoint of the provided values.")
        min_factors <- 3
        max_factors <- 12
        set_num_levels <- 3
        factor_df <- data.frame(
          Factor = rep("", min_factors),
          Level_Values = rep("-1, 0, 1", min_factors)
        )
        colnames(factor_df) <- c("Factor Names", "Level Values")
        column_widths <- c(200, 800)
      } else if (selected_design == "Full Factorial") {
        # Settings for Full Factorial Designs
        factors_string <- paste0(strong("Full Factorial Designs"), " are supported for
                                1-25 factors, which can be continuous, discrete, categorical or blocking.
                                For each factor, 2-10 levels can be considered, select the appropriate
                                number of levels for the factor from the corresponding dropdown.
                                Please provide the values for the levels as comma separated values, e.g., '5, 10, 15'.")
        min_factors <- 1
        max_factors <- 25
        num_levels <- 10
        set_num_levels <- NULL
        factor_df <- data.frame(
          Factor = rep("", 2),
          Type = factor(rep("Continuous", 2), levels = c("Continuous", "Discrete", "Categorical")),
          Levels = factor(rep(2, 2), levels = 2:num_levels),
          Level_Values = rep("-1, 1", 2)
        )
        colnames(factor_df) <- c("Factor Names", "Factor Type", "Num Levels", "Level Values")
        column_widths <- c(150, 150, 100, 600)
      } else if (selected_design == "Fractional Factorial") {
        # Settings for Fractional Factorial Designs
        factors_string <- paste0(strong("Fractional Factorial"), "
                                are supported for 3 to 25 factors, which can be continuous,
                                categorical, or discrete.  Each factor is limited
                                to two levels. Please provide the values for the
                                levels as comma separated values, e.g., '-1, 1'.
                                When setting up the design, select
                                the required resolution from the dropdown menu.
                                Note: Resolution III requires at least 3 factors,
                                Resolution IV requires at least 4 factors, and
                                Resolution V requires at least 5 factors.
                                are supported for
                                1-25 factors.  These factors can be continuous or categorical.
                                For each factor, 2-10 levels can be considered, select the appropriate
                                number of levels for the factor from the corresponding dropdown.
                                ")
        min_factors <- 3
        max_factors <- 25
        set_num_levels <- 2
        factor_df <- data.frame(
          Factor = rep("", min_factors),
          Type = factor(rep("Continuous", min_factors), levels = c("Continuous", "Discrete", "Categorical")),
          Level_Values = rep("-1, 1", min_factors)
        )
        colnames(factor_df) <- c("Factor Names", "Factor Type", "Level Values")
        column_widths <- c(200, 200, 700)
      } else if (selected_design == "Central Composite") {
        # Settings for Central Composite Design
        factors_string <- paste0(strong("Central Composite"), " are supported for
                                2 more factors.  These factors must be continuous.
                                For each factor, three levels are considered.
                                Please provide the upper and lower limits for these levels as comma separated values, e.g., '-1, 1'.
                                The mid level will be
                                automatically calculated as the midpoint of the provided values.")
        min_factors <- 2
        max_factors <- Inf
        set_num_levels <- 2
        factor_df <- data.frame(
          Factor = rep("", min_factors),
          Level_Values = rep("-1, 1", min_factors)
        )
        colnames(factor_df) <- c("Factor Names", "Level Values")
        column_widths <- c(200, 800)
      } else if (selected_design == "Latin Hypercube Sampling") {
        # Settings for Latin Hypercube Designs
        factors_string <- paste0(strong("Latin Hypercube Designs"), " are supported for
                                any number of factors or runs.  These factors can be continuous
                                discrete, categorical, or blocks.  For continuous factors, provide
                                the upper and lower limits as comma separated values, e.g., -1, 1.
                                For discrete, categorical, or blocking factors, provide the allowable values
                                as comma separated values.  Note: if blocking factors are included,
                                it is recommended to choose a number of runs that is divisible by the number of blocks.")
        min_factors <- 1
        max_factors <- Inf
        set_num_levels <- NULL
        factor_df <- data.frame(
          Factor = rep("", min_factors),
          Type = factor(rep("Continuous", min_factors), levels = c("Continuous", "Discrete", "Categorical")),
          Level_Values = rep("-1, 1", min_factors)
        )
        colnames(factor_df) <- c("Factor Names", "Factor Type", "Level Values")
        column_widths <- c(200, 200, 700)
      } else if (selected_design == "D-Optimal") {
        # Settings for D-Optimal Designs
        factors_string <- paste0(strong("D-Optimal Designs"), " are supported for
                                any number of factors, which can be continuous
                                discrete, categorical, or blocks.  For continuous factors, provide
                                the upper and lower limits as comma separated values, e.g., -1, 1.
                                For discrete, categorical, or blocking factors, provide the allowable values
                                as comma separated values.  Note: if blocking factors are included,
                                it is recommended to choose a number of runs that is divisible by the number of blocks.")
        min_factors <- 1
        max_factors <- Inf
        set_num_levels <- NULL
        factor_df <- data.frame(
          Factor = rep("", min_factors),
          Type = factor(rep("Continuous", min_factors), levels = c("Continuous", "Discrete", "Categorical")),
          Level_Values = rep("-1, 1", min_factors)
        )
        colnames(factor_df) <- c("Factor Names", "Factor Type", "Level Values")
        column_widths <- c(200, 200, 700)
      }

      column_widths_reactive(column_widths)
      max_factors_reactive(max_factors)
      min_factors_reactive(min_factors)
      set_num_levels_reactive(set_num_levels)
      factors_string_reactive(factors_string)

      showModal(modalDialog(
        title = h3(strong("Provide Factor Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(
            12,
            p(HTML(factors_string))
          ),
          column(
            12,
            p(strong("All entiries in the table below must be completed.  Factor Names can only contain letters, numbers, and underscores.  Do NOT use spaces or special characters in the factor names."), style = "margin-bottom: 25px;")
          )
        ),
        p("Note:  Add or remove rows by right clicking on a row."),
        div(
          style = "height: 220px;",
          rHandsontableOutput("Factor_details_table_id", height = "100%")
        ),
        uiOutput("warning_message_id"),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2,
              offset = 8,
              actionButton(inputId = "open_blocking_modal_id", "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")
            ),
            column(
              2,
              actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
            )
          )
        )
      ))
      output$Factor_details_table_id <- renderRHandsontable({
        rhandsontable(factor_df, selectCallback = TRUE, width = 1150, height = 250) %>%
          hot_cols(colWidths = column_widths) %>%
          hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
      })
    }
  })

  # Ensure the user remains within the allowed number of factors for a given design
  observe({
    dim(input$Factor_details_table_id)[1]
    current_factor_data <- hot_to_r(input$Factor_details_table_id)
    if (!is.null(dim(current_factor_data)[1])) {
      if (dim(current_factor_data)[1] > max_factors_reactive()) {
        showNotification(strong(paste0("The maximum number of factors for this design is ", max_factors_reactive(), ".")), duration = 20, closeButton = TRUE, type = "error")
      }
      if (dim(current_factor_data)[1] < min_factors_reactive()) {
        showNotification(strong(paste0("The minimum number of factors for this design is ", min_factors_reactive(), ".")), duration = 20, closeButton = TRUE, type = "error")
      }
    }
  })

  # User clicks next after providing factor details, opens Design Details modal
  # This confirms the design including blocking options, if available
  observeEvent(input$open_blocking_modal_id, {
    # Get the factor data entered by the user
    current_factor_data <- hot_to_r(input$Factor_details_table_id)
    factor_data_reactive(current_factor_data)

    # Reviews the users entries to ensure all requirements are met
    warning_msg <- c()

    # Loops through each factor
    for (i in seq_len(dim(current_factor_data)[1])) {
      warning_msg_temp <- c()

      # Checks for Missing Factor Name
      factor_name <- current_factor_data$`Factor Names`[i]
      if (is.na(factor_name) || grepl("^\\s*$", factor_name) || factor_name == "") {
        warning_msg_temp <- c(warning_msg_temp, paste0("factor name is missing"))
      }
      
      # Check factor name is safe
      if (is.na(factor_name) ||
          grepl("^\\s*$", factor_name) ||
          factor_name == "" ||
          factor_name != make.names(factor_name)) {
        
        warning_msg_temp <- c(warning_msg_temp, "factor name cannot contain spaces or special characters")
      }

      # If Factor Type is a required field, ensure it has been specified otherwise, assume continuous
      if ("Factor Type" %in% colnames(current_factor_data)) {
        factor_type <- current_factor_data$`Factor Type`[i]
        if (is.na(as.character(factor_type))) {
          warning_msg_temp <- c(warning_msg_temp, paste0("factor type is missing"))
        }
      } else {
        factor_type <- "Continuous"
      }

      # If the number of levels must be specified, ensure the user has provided it
      if ("Num Levels" %in% colnames(current_factor_data)) {
        num_levels <- current_factor_data$`Num Levels`[i]
        if (is.na(as.character(num_levels))) {
          warning_msg_temp <- c(warning_msg_temp, paste0("number of levels is missing"))
        }
      } else {
        num_levels <- set_num_levels_reactive()
      }

      # Checks the user has provided level values, that they are in the appropriate format (numeric for continuous), and that the appropriate number has been provided
      level_values <- current_factor_data$`Level Values`[i]
      if (is.na(level_values) || grepl("^\\s*$", level_values) || level_values == "") {
        warning_msg_temp <- c(warning_msg_temp, paste0("level values are missing"))
      } else {
        level_values_vector <- gsub("^\\s+|\\s+$", "", unlist(strsplit(level_values, ",")))
        if (!is.null(num_levels)) {
          if (length(level_values_vector) != num_levels) {
            warning_msg_temp <- c(warning_msg_temp, paste0("the number of level values provided does not match the number of levels selected or required by the design"))
          }
        }

        if (factor_type == "Continuous") {
          if (!all(!is.na(as.numeric(level_values_vector)))) {
            warning_msg_temp <- c(warning_msg_temp, paste0("the level values for a Continuous factor must be numeric"))
          }
        }
      }

      # Combines all of the warning messages per row
      if (length(warning_msg_temp) != 0) {
        warning_msg_temp <- paste0("Row ", i, ": ", paste0(warning_msg_temp, collapse = "; "), ".")
        warning_msg <- c(warning_msg, warning_msg_temp)
      }
    }

    # If there are warning messages (i.e. key info is missing), print the warning messages and stop the user from proceeding
    if (length(warning_msg) != 0) {
      output$warning_message_id <- renderUI({
        lapply(seq_len(length(warning_msg)), function(j) {
          p(strong(warning_msg[j]), class = "text-red")
        })
      })
    }
    req(length(warning_msg) == 0, cancelOutput = TRUE)

    # If all the info is inputted, close the model and start the blocking process
    removeModal()

    selected_design <- input$type_of_design_select_id
    n_factors <- dim(current_factor_data)[1]

    # Open a modal with the options for blocking and randomization
    showModal(modalDialog(
      title = h3(strong("Design Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      fluidRow(
        column(
          6,
          p(HTML(paste0(strong("Selected Design: "), selected_design))),
          DT::dataTableOutput("summary_factors_table_id")
        ),
        column(
          6,
          uiOutput("design_specific_selections_id")
        )
      ),
      easyClose = FALSE,
      footer = tagList(
        fluidRow(
          column(2,
            offset = 6,
            actionButton(inputId = "reopen_factors_modal_id", "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")
          ),
          column(
            2,
            actionButton(inputId = "create_design_id", "Create", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")
          ),
          column(
            2,
            actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
          )
        )
      )
    ))

    # Show the summary of the factors
    output$summary_factors_table_id <- DT::renderDataTable(
      current_factor_data,
      escape = FALSE, selection = "none", server = FALSE,
      extensions = c("FixedColumns", "FixedHeader", "Scroller"),
      options = list(
        dom = "t",
        pageLength = 10, paging = TRUE,
        headerCallback = DT::JS("function(thead) {", "  $(thead).css('font-size', '1em');", "}")
      ),
      rownames = FALSE, filter = "none"
    )

    # Based on the selected design, checks if blocking is available and presents options to user
    # Shows the total number of shows and provides and option to randomize runs
    if (selected_design == "Plackett-Burman") {
      design <- Plackett_Burman_Designs(n_factors)
      num_runs <- dim(design)[1]

      output$design_specific_selections_id <- renderUI({
        fluidRow(
          p(HTML(paste0(strong("Number of Runs:  "), num_runs))),
          p("Blocking is not available for this design"),
          materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
        )
      })
    } else if (selected_design == "Full Factorial") {
      factor_levels <- as.numeric(as.character(current_factor_data$`Num Levels`))
      design <- Full_Factorial_Designs(factor_levels)
      num_runs <- dim(design)[1]

      output$design_specific_selections_id <- renderUI({
        fluidRow(
          p(HTML(paste0(strong("Number of Runs:  "), num_runs))),
          p("Blocking is not available for this design"),
          materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
        )
      })
    } else if (selected_design == "Fractional Factorial") {
      n_factors <- dim(current_factor_data)[1]
      min_power <- ceiling(log2(n_factors + 1))
      max_power <- n_factors - 1
      possible_powers <- min_power:max_power
      possible_runs <- sort(2^possible_powers)
      num_runs <- possible_runs[1]

      min_runs_per_block_power <- min_power
      max_runs_per_block_power <- log2(num_runs / 2)

      if (max_runs_per_block_power < min_runs_per_block_power) {
        possible_blocks <- c(1)
      } else if (2^max_runs_per_block_power == num_runs) {
        possible_blocks <- c(1)
      } else {
        possible_blocks <- sort(num_runs / 2^(min_runs_per_block_power:max_runs_per_block_power))
      }

      output$design_specific_selections_id <- renderUI({
        fluidRow(
          column(
            5,
            selectInput("fractional_factorial_runs_select_id", "Select number of runs:", choices = possible_runs, selected = possible_runs[1], width = "100%")
          ),
          column(
            1,
            div(
              style = "margin-left: -25px; margin-top: 15px;",
              dropdownButton(
                inputId = "dropdownButton_id",
                p(paste0(
                  "These are the possible numbers of runs for a fractional factorial containing ",
                  n_factors, " factors.  Selecting the desired number of runs will update the possible
                                         blocking options and resolution of the design"
                )),
                circle = TRUE, icon = icon("info"), tooltip = tooltipOptions(title = "Click for details.")
              )
            ),
          ),
          column(
            5,
            selectInput("fractional_factorial_blocks_select_id", "Select number of blocks", choices = possible_blocks, selected = possible_blocks[1], width = "100%")
          ),
          column(
            1,
            div(
              style = "margin-left: -25px; margin-top: 15px;",
              dropdownButton(
                inputId = "dropdownButton_id",
                p(paste0(
                  "These are the possible numbers of blocks for a fractional factorial containing ",
                  n_factors, " factors and the chosen number of runs.  A value of 1 indicates a single block, i.e.
                                         no blocking is implemented."
                )),
                circle = TRUE, icon = icon("info"), tooltip = tooltipOptions(title = "Click for details.")
              )
            ),
          ),
          column(
            12,
            uiOutput("fractional_factorial_resolution_ui_id")
          ),
          column(
            12,
            materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        )
      })
    } else if (selected_design == "Box-Behnken") {
      n_factors <- dim(current_factor_data)[1]

      design <- Box_Behnken_Designs(n_factors)
      blocks <- design$blocks

      block_options <- unique(c(1, length(unique(blocks))))

      output$design_specific_selections_id <- renderUI({
        fluidRow(
          column(
            5,
            selectInput("fractional_factorial_blocks_select_id", "Select number of blocks", choices = block_options, selected = block_options[1], width = "100%")
          ),
          column(
            1,
            div(
              style = "margin-left: -25px; margin-top: 15px;",
              dropdownButton(
                inputId = "dropdownButton_id",
                p(paste0(
                  "These are the possible numbers of blocks for a Box-Behnken Design containing ",
                  n_factors, " factors and the chosen number of runs.  A value of 1 indicates a single block, i.e.
                                         no blocking is implemented."
                )),
                circle = TRUE, icon = icon("info"), tooltip = tooltipOptions(title = "Click for details.")
              )
            ),
          ),
          column(
            12,
            materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        )
      })
    } else if (selected_design == "Central Composite") {
      n_factors <- dim(current_factor_data)[1]
      design <- Central_Composite_Designs(n_factors, "Circumscribed")
      num_runs <- dim(design)[1]

      output$design_specific_selections_id <- renderUI({
        fluidRow(
          column(
            6,
            selectInput(inputId = "central_composite_type_id", label = "Select Central Composite Design Type:", choices = c("Circumscribed", "Inscribed", "Face Centered"), selected = "Circumscribed", width = "100%")
          ),
          column(
            6,
            div(
              style = "margin-left: -25px; margin-top: 15px;",
              dropdownButton(
                inputId = "dropdownButton_id",
                p(HTML(
                  "<ul>
                                      <li style='text-align: justify;'><strong>Circumscribed </strong>- the axial points lie outside the range of the factorial design.</li>
                                      <li style='text-align: justify;'><strong>Inscribed </strong>- the axial points lie within the provided upper and lower bounds and the factorial design is 'shrunk' to fit within the range to preserve rotatability.</li>
                                      <li style='text-align: justify;'><strong>Face Centered </strong>- the axial points lie on the center of each face of the factorial design. These axial points do not exceed the bounds of the factorial design. Note: this option is not rotatable.</li>
                                      </ul>"
                )),
                circle = TRUE, icon = icon("info"), width = "100%", tooltip = tooltipOptions(title = "Click for details.")
              ),
            )
          ),
          column(
            12,
            p(HTML(paste0(strong("Number of Runs:  "), num_runs))),
            p("Blocking is not available for this design"),
            materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        )
      })
    } else if (selected_design == "Latin Hypercube Sampling") {
      n_factors <- dim(current_factor_data)[1]

      output$design_specific_selections_id <- renderUI({
        fluidRow(
          column(
            6,
            numericInput(inputId = "LHS_runs_id", label = "Number of Experimental Runs:", min = 1, max = Inf, value = n_factors + 1, step = 1, width = "100%")
          ),
          column(
            6,
            div(
              style = "margin-left: -25px; margin-top: 15px;",
              dropdownButton(
                inputId = "dropdownButton_id",
                p(HTML(
                  "<p style='text-align: justify;'>A good rule of
                                      thumb for the minimum number of runs in a Latin
                                      Hypercube design is the number of factors plus 1
                                      (or the number of model parameters plus 1 if
                                      higher-order effects are modeled). While more
                                      runs will generally improve parameter estimation,
                                      the number is ultimately constrained by your experimental capacity.</p>"
                )),
                circle = TRUE, icon = icon("info"), width = "100%", tooltip = tooltipOptions(title = "Click for details.")
              ),
            )
          ),
          column(
            12,
            p("Blocking is not available for this design"),
            materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        )
      })
    } else if (selected_design == "D-Optimal") {
      n_factors <- dim(current_factor_data)[1]
      factor_names <- as.character(current_factor_data$`Factor Names`)
      factor_types <- as.character(current_factor_data$`Factor Type`)
      level_values <- as.character(current_factor_data$`Level Values`)
      level_values <- sapply(strsplit(level_values, ","), function(x) length(x))


      output$design_specific_selections_id <- renderUI({
        fluidPage(
          fluidRow(
            column(
              12,
              p("For a D-Optimal design, it's essential to identify the effects you
                wish to estimate in advance. Please use the interface below to select
                all relevant effects, including main effects, interactions, quadratic
                terms, and any other higher-order effects.")
            )
          ),
          fluidRow(
            column(
              4,
              actionButton(inputId = "DOD_main_effects_add_btn_id", label = "Main Effects", width = "100%", class = "Blue-Button")
            ),
            column(
              4,
              actionButton(inputId = "DOD_2way_effects_add_btn_id", label = "2-Way Interactions", width = "100%", class = "Blue-Button")
            ),
            column(
              4,
              actionButton(inputId = "DOD_quadratics_add_btn_id", label = "Quadratics", width = "100%", class = "Blue-Button")
            )
          ),
          div(style = "height: 10px;"),
          fluidRow(
            column(
              4,
              actionButton(inputId = "DOD_add_btn_id", label = "Add", width = "100%", class = "Blue-Button")
            ),
            column(
              4,
              actionButton(inputId = "DOD_cross_btn_id", label = "Create Interaction", width = "100%", class = "Blue-Button")
            ),
            column(
              4,
              actionButton(inputId = "DOD_remove_btn_id", label = "Remove", width = "100%", class = "Blue-Button")
            )
          ),
          div(style = "height: 10px;"),
          fluidRow(
            column(
              4,
              selectInput(inputId = "available_factors_DOD_id", label = "Factor:", choices = factor_names, size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
            ),
            column(
              8,
              selectInput(inputId = "selected_factors_DOD_id", label = "Selected Effects:", choices = NULL, size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
            )
          ),
          fluidRow(
            column(
              4,
              numericInput(inputId = "DOD_num_runs_id", label = "Number of Runs:", min = 1, max = Inf, value = NULL, step = 1, width = "100%")
            ),
            column(
              3,
              uiOutput("DOD_min_runs_ui_id")
            ),
            column(
              3,
              uiOutput("DOD_max_runs_ui_id")
            )
          ),
          fluidRow(
            column(
              12,
              materialSwitch(inputId = "randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
            )
          )
        )
      })
    }
  })

  # Re-opens the modal for specifying the factor information if the user hits back on the design details modal
  observeEvent(input$reopen_factors_modal_id, {
    removeModal()

    showModal(modalDialog(
      title = h3(strong("Provide Factor Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      fluidRow(
        column(
          12,
          h4(HTML(factors_string_reactive()), style = "margin-bottom: 25px;")
        )
      ),
      p("Note:  Add or Remove columns by right clicking on the table."),
      div(
        rHandsontableOutput("Factor_details_table_id")
      ),
      uiOutput("warning_message_id"),
      easyClose = FALSE,
      footer = tagList(
        fluidRow(
          column(2,
            offset = 8,
            actionButton(inputId = "open_blocking_modal_id", "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")
          ),
          column(
            2,
            actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
          )
        )
      )
    ))
    output$Factor_details_table_id <- renderRHandsontable({
      rhandsontable(factor_data_reactive(), selectCallback = TRUE, width = "1150px", height = "250px", stretchH = "all") %>%
        hot_cols(colWidths = column_widths_reactive()) %>%
        hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })
  })

  # Calculates the possible blocks for a fractional factorial design based on the number of runs selected by the user.
  observeEvent(input$fractional_factorial_runs_select_id, {
    current_factor_data <- factor_data_reactive()
    n_factors <- dim(current_factor_data)[1]

    min_power <- ceiling(log2(n_factors + 1))

    max_power <- n_factors - 1
    num_runs <- as.numeric(input$fractional_factorial_runs_select_id)

    min_runs_per_block_power <- min_power
    max_runs_per_block_power <- log2(num_runs / 2)

    if (max_runs_per_block_power < min_runs_per_block_power) {
      possible_blocks <- c(1)
    } else if (2^max_runs_per_block_power == num_runs) {
      possible_blocks <- c(1)
    } else {
      possible_blocks <- sort(num_runs / 2^(min_runs_per_block_power:max_runs_per_block_power))
    }
    if (!1 %in% possible_blocks) {
      possible_blocks <- c(1, possible_blocks)
    }

    updateSelectInput(session, inputId = "fractional_factorial_blocks_select_id", choices = possible_blocks, selected = possible_blocks[1])


    design <- Fractional_Factorial_Designs(num_runs, n_factors, 1)

    output$fractional_factorial_resolution_ui_id <- renderUI({
      p(HTML(paste0("This is a ", strong(design$Resolution), " design.")))
    })
  })

  # ===== For the D-Optimal Design
  # Helper function
  classify_terms <- function(term) {
    length(unlist(strsplit(term, ":")))
  }

  # Adds all main effects to the D-Opt Design
  observeEvent(input$DOD_main_effects_add_btn_id, {
    current_factor_data <- factor_data_reactive()
    factor_names <- as.character(current_factor_data$`Factor Names`)

    currently_selected_factors <- d_optimal_model_components()

    updated_selected_factors <- unique(c(currently_selected_factors, factor_names))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    d_optimal_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = sorted_terms)
  })

  # Adds all 2-factor interactions to the D-Opt Design
  observeEvent(input$DOD_2way_effects_add_btn_id, {
    current_factor_data <- factor_data_reactive()
    factor_names <- as.character(current_factor_data$`Factor Names`)
    req(length(factor_names) > 1, cancelOutput = TRUE)

    currently_selected_factors <- d_optimal_model_components()

    interactions <- combn(factor_names, 2, FUN = function(x) paste(x, collapse = ":"))

    updated_selected_factors <- unique(c(currently_selected_factors, interactions))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    d_optimal_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = sorted_terms)
  })

  # Adds all quadratics to the D-Opt Design
  observeEvent(input$DOD_quadratics_add_btn_id, {
    current_factor_data <- factor_data_reactive()
    factor_names <- as.character(current_factor_data$`Factor Names`)

    currently_selected_factors <- d_optimal_model_components()

    quadratics <- paste(factor_names, factor_names, sep = ":")

    updated_selected_factors <- unique(c(currently_selected_factors, quadratics))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    d_optimal_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = sorted_terms)
  })

  # Add selected factor to the D-Opt Design
  observeEvent(input$DOD_add_btn_id, {
    selected_available_factors <- input$available_factors_DOD_id

    if (is.null(selected_available_factors)) {
      showNotification("Select a factor to add", duration = 20, closeButton = TRUE, type = "warning")
    }
    req(!(is.null(selected_available_factors)), cancelOutput = TRUE)

    currently_selected_factors <- d_optimal_model_components()

    updated_selected_factors <- unique(c(currently_selected_factors, selected_available_factors))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    d_optimal_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = sorted_terms)
  })

  # Cross selected factors and add to the D-Opt Design
  observeEvent(input$DOD_cross_btn_id, {
    currently_selected_factors <- d_optimal_model_components()

    selected_available_factors <- input$available_factors_DOD_id
    selected_model_effects <- input$selected_factors_DOD_id

    if (is.null(selected_available_factors) || is.null(selected_model_effects)) {
      showNotification("Select at least one factor and one effect to create an interaction.", duration = 20, closeButton = TRUE, type = "warning")
    }
    req(!is.null(selected_available_factors), !is.null(selected_model_effects), cancelOutput = TRUE)

    combinations <- expand.grid(selected_available_factors, selected_model_effects)

    interactions <- paste(combinations[, 1], combinations[, 2], sep = ":")

    updated_selected_factors <- unique(c(currently_selected_factors, interactions))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    d_optimal_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = sorted_terms)
  })

  # Remove selected effect from the D-Opt Design
  observeEvent(input$DOD_remove_btn_id, {
    currently_selected_factors <- d_optimal_model_components()

    selected_model_effects <- input$selected_factors_DOD_id
    if (is.null(selected_model_effects)) {
      showNotification("Select an effect to remove", duration = 20, closeButton = TRUE, type = "warning")
    }
    req(!(is.null(selected_model_effects)), cancelOutput = TRUE)

    updated_selected_factors <- setdiff(currently_selected_factors, selected_model_effects)
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    d_optimal_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = sorted_terms)
  })

  # As the selected effects for the D-Opt Design change,
  # this calculates and updates the min and max number of runs
  DOD_min_runs_reactive <- reactiveVal(NULL)
  DOD_max_runs_reactive <- reactiveVal(NULL)
  
  observeEvent(d_optimal_model_components(), {
    req(!is.null(d_optimal_model_components()), cancelOutput = TRUE)
    
    components <- d_optimal_model_components()
    if (length(components) == 0) {
      output$Augment_min_runs_ui_id <- renderUI(p(strong("Min Runs: NA")))
      output$Augment_max_runs_ui_id <- renderUI(p(strong("Max Runs: NA")))
      return(NULL)
    }
    
    factor_data <- factor_data_reactive()
    
    factor_info <- lapply(seq_len(nrow(factor_data)), function(i) {
      fname <- factor_data$`Factor Names`[i]
      ftype <- factor_data$`Factor Type`[i]
      levels <- strsplit(factor_data$`Level Values`[i], ",")[[1]]
      list(name = fname, type = ftype, levels = levels)
    })
    
    names(factor_info) <- factor_data$`Factor Names`
    
    term_params <- function(term) {
      factors <- unlist(strsplit(term, ":"))
      tab <- table(factors)
      
      prod(sapply(names(tab), function(f) {
        info <- factor_info[[f]]
        if (info$type == "Continuous") {
          1
        } else {
          length(info$levels) - 1
        }
      }))
    }
    
    num_parameters <- sum(sapply(components, term_params)) + 1
    
    factor_degree <- function(f) {
      deg <- 1
      for (term in components) {
        parts <- unlist(strsplit(term, ":"))
        count <- sum(parts == f)
        if (count > deg) deg <- count
      }
      deg
    }
    
    max_runs <- prod(sapply(names(factor_info), function(f) {
      info <- factor_info[[f]]
      if (info$type == "Continuous") {
        deg <- factor_degree(f)
        deg + 1
      } else {
        length(info$levels)
      }
    }))
    
    if (num_parameters > max_runs) {
      num_parameters <- max_runs
    }
    
    DOD_min_runs_reactive(num_parameters)
    DOD_max_runs_reactive(max_runs)
    
    output$DOD_min_runs_ui_id <- renderUI(p(strong(paste0("Min Runs: ", num_parameters))))
    output$DOD_max_runs_ui_id <- renderUI(p(strong(paste0("Max Runs: ", max_runs))))
  })
  

  # ===== Creates the Design
  observeEvent(input$create_design_id, {
    # Get the factor data, the selected design, and if the runs should be randomized
    factor_data <- factor_data_reactive()
    selected_design <- input$type_of_design_select_id
    n_factors <- dim(factor_data)[1]
    factor_names <- factor_data$`Factor Names`

    randomize <- input$randomize_checkbox_id

    # Helper function to check that the number of runs selected by the
    # user is in the right format for LHS and D-Opt
    validate_num_runs <- function(num_runs) {
      if (is.character(num_runs) || is.factor(num_runs)) {
        suppressWarnings(num_runs <- as.numeric(as.character(num_runs)))
        return(FALSE)
      }
      
      # Check it's numeric
      if (!is.numeric(num_runs) || length(num_runs) != 1 || is.na(num_runs)) {
        showNotification("The 'Number of Experimental Runs' must be a single numeric integer.", duration = 20, closeButton = TRUE)
        return(FALSE)
      }

      # Check it's an integer (within floating-point tolerance)
      if (abs(num_runs - round(num_runs)) > .Machine$double.eps^0.5) {
        showNotification("The 'Number of Experimental Runs' must be an integer (no decimal values).", duration = 20, closeButton = TRUE)
        return(FALSE)
      }

      # Check it's positive
      if (num_runs <= 0) {
        showNotification("The 'Number of Experimental Runs' must be greater than zero.")
        return(FALSE)
      }
      
      # check it is within the min and max number of runs for a DOD
      if (!(num_runs >= DOD_min_runs_reactive() && num_runs <= DOD_max_runs_reactive())) {
        showNotification("The 'Number of Experimental Runs' must be within the minumum and maximum number of runs provided.")
        return(FALSE) 
      }
      TRUE
    }


    # Creates the design based on the users choices
    if (selected_design == "Plackett-Burman") {
      design <- Plackett_Burman_Designs(n_factors)
      colnames(design) <- factor_names
    } else if (selected_design == "Full Factorial") {
      factor_levels <- as.numeric(as.character(factor_data$`Num Levels`))
      design <- Full_Factorial_Designs(factor_levels)
      colnames(design) <- factor_names
    } else if (selected_design == "Fractional Factorial") {
      num_runs <- as.numeric(input$fractional_factorial_runs_select_id)
      num_blocks <- as.numeric(input$fractional_factorial_blocks_select_id)

      design <- Fractional_Factorial_Designs(num_runs, n_factors, num_blocks)$design

      if (num_blocks > 1) {
        blocks <- design$blocks
        design <- subset(design, select = -blocks)
        design <- cbind(design, blocks)
        colnames(design) <- c(factor_names, "Block")
      } else {
        colnames(design) <- factor_names
      }
    } else if (selected_design == "Box-Behnken") {
      num_blocks <- input$fractional_factorial_blocks_select_id
      design <- Box_Behnken_Designs(n_factors)

      if (num_blocks == 1) {
        if ("Blocks" %in% colnames(design)) {
          design <- subset(design, select = -blocks)
        }
        colnames(design) <- factor_names
      } else {
        blocks <- design$blocks
        design <- subset(design, select = -blocks)
        design <- cbind(design, blocks)
        colnames(design) <- c(factor_names, "Block")
      }
    } else if (selected_design == "Central Composite") {
      ccd_type <- input$central_composite_type_id
      design <- Central_Composite_Designs(n_factors, ccd_type)
      colnames(design) <- factor_names
    } else if (selected_design == "Latin Hypercube Sampling") {
      req(validate_num_runs(input$LHS_runs_id), cancelOutput = TRUE)
      num_runs <- as.integer(round(input$LHS_runs_id))

      design <- Latin_Hypercube_Designs(n_factors, num_runs)
      colnames(design) <- factor_names
    } else if (selected_design == "D-Optimal") {
      model_components <- d_optimal_model_components()


      req(validate_num_runs(input$DOD_num_runs_id), cancelOutput = TRUE)
      num_runs <- as.integer(round(input$DOD_num_runs_id))

      level_values <- as.character(factor_data$`Level Values`)
      level_values <- sapply(strsplit(level_values, ","), function(x) length(x))
      factor_types <- as.character(factor_data$`Factor Type`)
      design <- D_Optimal_Designs(level_values, model_components, factor_types, num_runs)
    }

    # Rescales the Design coded values to the levels provided by the user
    if (selected_design == "Central Composite") {
      bounds <- factor_data$`Level Values`
      for (i in seq_len(n_factors)) {
        lower_upper <- as.numeric(unlist(strsplit(bounds[i], ",")))
        lower <- min(lower_upper)
        upper <- max(lower_upper)
        design[, i] <- lower + (design[, i] + 1) * (upper - lower) / 2
      }
    } else if (selected_design == "Latin Hypercube Sampling") {
      bounds <- factor_data$`Level Values`
      factor_types <- factor_data$`Factor Type`
      for (i in seq_len(ncol(design))) {
        factor_type <- factor_types[i]
        bound_str <- bounds[i]

        if (factor_type == "Continuous") {
          lower_upper <- as.numeric(unlist(strsplit(bound_str, ",")))
          lower <- min(lower_upper)
          upper <- max(lower_upper)
          design[, i] <- lower + design[, i] * (upper - lower)
        } else if (factor_type == "Discrete") {
          allowable_values <- sort(as.numeric(unlist(strsplit(bound_str, ","))))
          min_val <- min(allowable_values)
          max_val <- max(allowable_values)

          scaled_values <- min_val + design[, i] * (max_val - min_val)

          design[, i] <- sapply(scaled_values, function(val) {
            allowable_values[which.min(abs(allowable_values - val))]
          })
        } else if (factor_type == "Categorical") {
          categories <- unlist(strsplit(bound_str, ","))
          num_categories <- length(categories)

          breaks <- seq(0, 1, length.out = num_categories + 1)
          design[, i] <- cut(design[, i], breaks = breaks, labels = categories, include.lowest = TRUE)
        }
      }
    } else {
      for (col_name in colnames(design)) {
        if (col_name %in% factor_data$`Factor Names`) {
          level_string <- factor_data$`Level Values`[factor_data$`Factor Names` == col_name]
          level_values <- unlist(strsplit(level_string, ","))
          level_values <- trimws(level_values)
          if (all(sapply(level_values, is.numeric))) {
            level_values <- sort(as.numeric(level_values))
          }

          unique_design_values <- sort(unique(design[[col_name]]))
          value_mapping <- setNames(level_values, unique_design_values)
          design[[col_name]] <- value_mapping[as.character(design[[col_name]])]
        }
      }
    }

    # Adds a Response Column
    new <- rep("", nrow(design))
    design[, ncol(design) + 1] <- new
    colnames(design)[ncol(design)] <- "Response"
    
    design[, ncol(design) + 1] <- new
    colnames(design)[ncol(design)] <- "Metadata"

    # Scans through each factor in the design and assigns a role (type) to it
    factor_types <- c()
    for (i in seq_along(colnames(design))) {
      factor_name <- colnames(design)[i]
      if (factor_name == "Block") {
        factor_types <- c(factor_types, "Block")
      } else if (grepl("Response", factor_name)) {
        factor_types <- c(factor_types, "Response")
      } else if (grepl("Metadata", factor_name)) {
        factor_types <- c(factor_types, "Metadata")
      } else {
        index <- which(factor_data$`Factor Names` == factor_name)
        factor_type <- as.character(factor_data$`Factor Type`[index])
        if (is.null(factor_type) || factor_type == "" || length(factor_type) == 0) {
          factor_types <- c(factor_types, "Continuous")
        } else {
          factor_types <- c(factor_types, factor_type)
        }
      }
    }

    factor_types_df <- as.data.frame(t(factor_types))
    colnames(factor_types_df) <- colnames(design)

    # Enfores continuous, discrete, and response columns are numeric, block and categorical will be character
    for (col_name in colnames(design)) {
      col_type <- factor_types_df[[col_name]]
      if (col_type %in% c("Continuous", "Discrete", "Response")) {
        design[[col_name]] <- as.numeric(design[[col_name]])
      } else if (col_type %in% c("Block", "Categorical", "Metadata")) {
        design[[col_name]] <- as.character(design[[col_name]])
      }
    }

    if (randomize) {
      design <- design[sample(nrow(design)), , drop = FALSE]
    }

    design_reactive(design)

    factor_types_reactive(factor_types_df)

    # Sets up Augment Design button in Enter/Edit Results tab
    output$augment_design_btn_ui_id <- renderUI({
      div(
        align = "center",
        actionButton(inputId = "augment_design_btn_id", label = "Augment Design", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
      )
    })

    # Sets up add, remove, update column buttons in Enter/Edit Results tab
    output$design_table_UI_id <- renderUI({
      fluidPage(
        fluidRow(
          column(
            4,
            div(
              align = "center",
              actionButton(inputId = "remove_column_id", label = "Remove Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
            )
          ),
          column(
            4,
            div(
              align = "center",
              actionButton(inputId = "add_column_id", label = "Add Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
            )
          ),
          column(
            4,
            div(
              align = "center",
              actionButton(inputId = "update_column_id", label = "Update Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
            )
          )
        ),
        rHandsontableOutput(outputId = "design_table_id")
      )
    })

    output$design_table_id <- renderRHandsontable({
      rhandsontable(design_reactive(), selectCallback = TRUE, width = "1050px", height = "350px", stretchH = "all") %>%
        hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })

    removeModal()
    updateTabItems(session, "tabs", "Enter_edit_results_id")
  })


  # Close the Modal and clear any warning message
  observeEvent(input$close_modal_id, {
    output$warning_message_id <- renderUI(p(""))
    removeModal()
  })


  # ===== Edit the Design ======================================================
  # Binds the factor types and design from the table and saves it to a .csv
  output$save_design_id <- downloadHandler(
    filename = function() {
      "DoEIY_Design.csv"
    },
    content = function(file) {
      write.csv(rbind(factor_types_reactive(), hot_to_r(input$design_table_id)), file, row.names = FALSE)
    }
  )

  # Select a file to load
  observeEvent(input$load_design_id, {
    showModal(modalDialog(
      title = h3(strong("Select Design File to Load"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      column(
        10,
        fileInput(inputId = "select_load_file_id", label = NULL, width = "100%", accept = c(".csv", "text/csv"))
      ),
      column(
        2,
        actionButton(inputId = "load_file_id", label = "Upload", class = "Blue-Button", width = "100%", style = "height: 34px; margin-bottom: 0px;")
      ),
      easyClose = FALSE,
      footer = tagList(
        column(2,
          offset = 10,
          actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      )
    ))
  })
  # Load the selected file
  observeEvent(input$load_file_id, {
    file_path <- input$select_load_file_id[[1, "datapath"]]
    design <- read.csv(file_path, header = TRUE, sep = ",")

    colnames(design) <- gsub("\\.", " ", colnames(design))

    design[is.na(design)] <- ""

    factor_types <- design[1, ]
    design <- design[-1, ]
    rownames(design) <- NULL

    factor_types_reactive(factor_types)
    design_reactive(design)

    removeModal()

    output$augment_design_btn_ui_id <- renderUI({
      div(
        align = "center",
        actionButton(inputId = "augment_design_btn_id", label = "Augment Design", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
      )
    })
    
    output$design_table_UI_id <- renderUI({
      fluidPage(
        fluidRow(
          column(
            4,
            div(
              align = "center",
              actionButton(inputId = "remove_column_id", label = "Remove Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
            )
          ),
          column(
            4,
            div(
              align = "center",
              actionButton(inputId = "add_column_id", label = "Add Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
            )
          ),
          column(
            4,
            div(
              align = "center",
              actionButton(inputId = "update_column_id", label = "Update Column", width = "200px", class = "Blue-Button", style = "margin-bottom: 10px; margin-top: 0px;")
            )
          )
        ),
        rHandsontableOutput(outputId = "design_table_id")
      )
    })

    output$design_table_id <- renderRHandsontable({
      rhandsontable(design_reactive(), selectCallback = TRUE, width = "1050px", height = "350px", stretchH = "all") %>%
        hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })
  })

  # Select column to remove
  observeEvent(input$remove_column_id, {
    design <- hot_to_r(input$design_table_id)

    columns <- colnames(design)

    showModal(modalDialog(
      title = h3(strong("Select columns to delete"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      fluidRow(
        column(
          6,
          pickerInput(inputId = "remove_selected_columns_id", label = NULL, choices = columns, multiple = TRUE, width = "100%")
        ),
        column(
          3,
          actionButton(inputId = "delete_columns_id", label = "Delete", class = "Blue-Button", width = "100%", style = "margin-bottom: 0px;")
        )
      ),
      easyClose = FALSE,
      footer = tagList(
        column(2,
          offset = 10,
          actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      )
    ))
  })
  # Remove column
  observeEvent(input$delete_columns_id, {
    selected_columns <- input$remove_selected_columns_id
    design <- hot_to_r(input$design_table_id)

    design <- design[, !(names(design) %in% selected_columns)]
    design_reactive(design)

    factor_types <- factor_types_reactive()
    factor_types <- factor_types[, !(names(factor_types) %in% selected_columns)]
    factor_types_reactive(factor_types)

    removeModal()

    output$design_table_id <- renderRHandsontable({
      rhandsontable(design_reactive(), selectCallback = TRUE, width = "1050px", height = "350px", stretchH = "all") %>%
        hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })
  })

  # Select a column to add
  observeEvent(input$add_column_id, {
    showModal(modalDialog(
      title = h3(strong("Add a column"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      fluidRow(
        column(
          3,
          textInput(inputId = "added_column_name_id", label = "Column Name:", width = "100%")
        ),
        column(
          3,
          selectInput(inputId = "added_column_type_id", label = "Column Type:", width = "100%", choices = c("Continuous", "Discrete", "Categorical", "Block", "Response", "Metadata"))
        ),
        column(
          3,
          selectInput(inputId = "added_column_position_id", label = "Add Before Column:", width = "100%", choices = c(colnames(design_reactive()), "Add to End"))
        ),
        column(
          3,
          actionButton(inputId = "add_column_modal_id", label = "Add Column", class = "Blue-Button", width = "100%", style = "margin-top: 17px;")
        )
      ),
      easyClose = FALSE,
      footer = tagList(
        column(2,
          offset = 10,
          actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      )
    ))
  })
  # Add the column
  observeEvent(input$add_column_modal_id, {
    column_name <- input$added_column_name_id
    if (column_name == "") {
      showNotification(strong("A column name must be provided."), duration = 15, closeButton = TRUE, type = "error")
    }
    req(column_name != "", cancelOutput = TRUE)

    design <- design_reactive()
    colnames_initial <- colnames(design)
    factor_types <- factor_types_reactive()

    if (any(column_name %in% colnames_initial)) {
      showNotification(strong("The column name cannot already be in use."), duration = 15, closeButton = TRUE, type = "error")
    }
    req(!(any(column_name %in% colnames_initial)), cancelOutput = TRUE)

    column_type <- input$added_column_type_id
    column_position <- input$added_column_position_id

    if (column_position == "Add to End") {
      insert_position <- length(colnames(design)) + 1
    } else {
      insert_position <- which(colnames(design) == column_position)
    }

    new <- rep("", nrow(design))
    new_factor <- column_type

    if (insert_position == 1) {
      design_updated <- cbind(new, design)
      factor_types_updated <- cbind(new_factor, factor_types)
      colnames_updated <- c(column_name, colnames_initial)
    } else if (insert_position == (ncol(design) + 1)) {
      design_updated <- cbind(design, new)
      factor_types_updated <- cbind(factor_types, new_factor)
      colnames_updated <- c(colnames_initial, column_name)
    } else {
      design_updated <- data.frame(
        design[, 1:(insert_position - 1)],
        New_Column = new,
        design[, insert_position:ncol(design)]
      )
      factor_types_updated <- data.frame(
        factor_types[, 1:(insert_position - 1)],
        New_Column = new_factor,
        factor_types[, insert_position:ncol(factor_types)]
      )
      colnames_updated <- c(colnames_initial[1:(insert_position - 1)], column_name, colnames_initial[insert_position:length(colnames_initial)])
    }

    colnames(design_updated) <- colnames_updated
    colnames(factor_types_updated) <- colnames_updated
    design_reactive(design_updated)
    factor_types_reactive(factor_types_updated)

    removeModal()
  })

  # Selects updates to column and specs to change
  observeEvent(input$update_column_id, {
    showModal(modalDialog(
      title = h3(strong("Update a column"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      fluidRow(
        selectInput(inputId = "column_to_update_id", label = "Select a Column to Update:", selected = NULL, choices = colnames(design_reactive()))
      ),
      fluidRow(
        column(
          3,
          textInput(inputId = "update_column_name_id", label = "Column Name:", width = "100%")
        ),
        column(
          3,
          selectInput(inputId = "update_column_type_id", label = "Column Type:", width = "100%", choices = c("Continuous", "Discrete", "Categorical", "Block", "Response", "Metadata"))
        ),
        column(
          3,
          selectInput(inputId = "update_column_position_id", label = "Move Before Column:", width = "100%", choices = c("", colnames(design_reactive()), "Move to End"))
        ),
        column(
          3,
          actionButton(inputId = "update_column_modal_id", label = "Update Column", class = "Blue-Button", width = "100%", style = "margin-top: 17px;")
        )
      ),
      easyClose = FALSE,
      footer = tagList(
        column(2,
          offset = 10,
          actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
        )
      )
    ))
  })
  # Provides the current name and factor type for the selected column for updating
  observeEvent(input$column_to_update_id, {
    updateTextInput(session, inputId = "update_column_name_id", value = input$column_to_update_id)
    updateSelectInput(session, inputId = "update_column_type_id", selected = as.character(factor_types_reactive()[[input$column_to_update_id]]))
  })
  # Assesses and implements changes to column
  observeEvent(input$update_column_modal_id, {

    column_name <- input$update_column_name_id
    if (column_name == "") {
      showNotification(strong("A column name must be provided."), duration = 15, closeButton = TRUE, type = "error")
    }
    req(column_name != "", cancelOutput = TRUE)

    design <- design_reactive()
    colnames_initial <- colnames(design)
    filtered_colnames <- colnames_initial[colnames_initial != column_name]
    factor_types <- factor_types_reactive()

    if (any(column_name %in% filtered_colnames)) {
      showNotification(strong("The column name cannot already be in use."), duration = 15, closeButton = TRUE, type = "error")
    }
    req(!(any(column_name %in% filtered_colnames)), cancelOutput = TRUE)

    column_type <- input$update_column_type_id
    column_position <- input$update_column_position_id

    if (column_position != "") {
      names(design)[names(design) == input$column_to_update_id] <- "Temp_column"
      names(factor_types)[names(factor_types) == input$column_to_update_id] <- "Temp_column"
      col_names_new <- colnames(design)

      if (column_position == "Move to End") {
        insert_position <- length(colnames(design)) + 1
      } else {
        insert_position <- which(colnames(design) == column_position)
      }

      column_values <- design$Temp_column

      if (insert_position == 1) {
        design_updated <- cbind(column_values, design)
        factor_types_updated <- cbind(column_type, factor_types)
        colnames_updated <- c(column_name, col_names_new)
      } else if (insert_position == (ncol(design) + 1)) {
        design_updated <- cbind(design, column_values)
        factor_types_updated <- cbind(factor_types, column_type)
        colnames_updated <- c(col_names_new, column_name)
      } else {
        design_updated <- data.frame(
          design[, 1:(insert_position - 1)],
          New_Column = column_values,
          design[, insert_position:ncol(design)]
        )
        factor_types_updated <- data.frame(
          factor_types[, 1:(insert_position - 1)],
          New_Column = column_type,
          factor_types[, insert_position:ncol(factor_types)]
        )
        colnames_updated <- c(col_names_new[1:(insert_position - 1)], column_name, col_names_new[insert_position:length(col_names_new)])
      }
      colnames(design_updated) <- colnames_updated
      colnames(factor_types_updated) <- colnames_updated

      design_updated$Temp_column <- NULL
      factor_types_updated$Temp_column <- NULL

      design_reactive(design_updated)
      factor_types_reactive(factor_types_updated)
    } else {
      factor_types[[input$column_to_update_id]] <- column_type
      names(design)[names(design) == input$column_to_update_id] <- column_name
      names(factor_types)[names(factor_types) == input$column_to_update_id] <- column_name

      design_reactive(design)
      factor_types_reactive(factor_types)
    }

    removeModal()
  })


  # Watches the user editable table; any time changes are made, the Design_Reactive is updated
  observe({
    req(input$design_table_id)
    design_reactive(hot_to_r(input$design_table_id))
  })
  
  
  augment_factor_data_reactive <- reactiveVal(NULL)
  augment_design_model_components <- reactiveVal(value = NULL)
  augment_original_factor_data_reactive <- reactiveVal(NULL)
  
  # Augment the Design
  observeEvent(input$augment_design_btn_id, {
    
    current_factor_types = factor_types_reactive()
    current_design = design_reactive()
    
    
    types_vec <- as.character(current_factor_types[1, ])
    names(types_vec) <- colnames(current_factor_types)
    
    keep_types <- c("Continuous", "Discrete", "Categorical")
    kept_names <- names(types_vec)[types_vec %in% keep_types]
    
    design_subset <- current_design[, kept_names, drop = FALSE]
    
    factor_names <- kept_names
    factor_types <- types_vec[kept_names]
    
    level_values <- sapply(factor_names, function(fn) {
      col <- design_subset[[fn]]
      
      if (factor_types[[fn]] == "Continuous") {
        # Range for continuous
        col <- as.numeric(col)
        rng <- range(col, na.rm = TRUE)
        paste0(rng[1], ",", rng[2])
        
      } else if (factor_types[[fn]] == "Discrete") {
        # Sorted unique numeric values
        col <- as.numeric(col)
        vals <- sort(unique(col))
        paste(vals, collapse = ",")
        
      } else if (factor_types[[fn]] == "Categorical") {
        # Sorted unique categories
        vals <- sort(unique(as.character(col)))
        paste(vals, collapse = ",")
      }
    })
    
    # Ensure incoming factor types match exactly
    valid_types <- c("Continuous", "Discrete", "Categorical")
    
    factor_types_clean <- trimws(factor_types)                  
    factor_types_clean <- gsub("[\r\n\t]", "", factor_types_clean)  
    factor_types_clean <- tools::toTitleCase(tolower(factor_types_clean))  
    
    factor_types_clean[!(factor_types_clean %in% valid_types)] <- "Continuous"
    
    
    current_factor_information <- data.frame(
      "Factor Names" = factor_names,
      "Factor Type" = factor_types_clean,
      "Level Values" = level_values,
      check.names = FALSE, stringsAsFactors = FALSE
    )
    
    augment_factor_data_reactive(current_factor_information)
    augment_original_factor_data_reactive(current_factor_information)
    
    showModal(
      modalDialog(
        title = h3(strong("Augment Existing Design"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        p("This feature is currently in beta, and results should be reviewed carefully before use.", style = "color:red;"), 
        fluidRow(
          column(
            12,
            h4(strong("Current Factors")),
            p("Edit types or level values.  Add or delete factors as needed before augmenting by right clicking on the table below and inserting a new row."),
            rHandsontableOutput("augment_factor_table_id")
          )
        ),
        fluidRow(
          column(4, offset = 8, 
            actionButton("augment_continue_btn_id", "Continue", class = "Blue-Button", width = "100%", style = "margin-top: 25px;")
          )
        ),

        
        easyClose = FALSE,
        footer = tagList(
          column(2, offset = 10,
            actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
          )
        )
      )
    )
  })
  
  output$augment_factor_table_id <- renderRHandsontable({
    req(augment_factor_data_reactive())
    
    
    rhandsontable(augment_factor_data_reactive(), rowHeaders = NULL, width = "1150", height = "250px") %>% 
      hot_cols(colWidths = c(200, 200, 700)) %>%
      hot_col(col = "Factor Type", type = "dropdown", source = c("Continuous", "Discrete", "Categorical"), allowInvalid = FALSE) %>%
      hot_context_menu(allowRowEdit = TRUE, allowColEdit = TRUE)
  })
  
  observeEvent(input$augment_continue_btn_id, {
    # Get the factor data entered by the user
    current_factor_data <- hot_to_r(input$augment_factor_table_id)
    factor_data_reactive(current_factor_data)
    # Reviews the users entries to ensure all requirements are met
    warning_msg <- c()
    
    # Loops through each factor
    for (i in seq_len(dim(current_factor_data)[1])) {
      warning_msg_temp <- c()
      
      # Checks for Missing Factor Name
      factor_name <- current_factor_data$`Factor Names`[i]
      if (is.na(factor_name) || grepl("^\\s*$", factor_name) || factor_name == "") {
        warning_msg_temp <- c(warning_msg_temp, paste0("factor name is missing"))
      }
      
      # Check factor name is safe
      if (is.na(factor_name) ||
          grepl("^\\s*$", factor_name) ||
          factor_name == "" ||
          factor_name != make.names(factor_name)) {
        
        warning_msg_temp <- c(warning_msg_temp, "factor name cannot contain spaces or special characters")
      }
      
      # If Factor Type is a required field, ensure it has been specified otherwise, assume continuous
      if ("Factor Type" %in% colnames(current_factor_data)) {
        factor_type <- current_factor_data$`Factor Type`[i]
        if (is.na(as.character(factor_type))) {
          warning_msg_temp <- c(warning_msg_temp, paste0("factor type is missing"))
        }
      } else {
        factor_type <- "Continuous"
      }
      
      # If the number of levels must be specified, ensure the user has provided it
      if ("Num Levels" %in% colnames(current_factor_data)) {
        num_levels <- current_factor_data$`Num Levels`[i]
        if (is.na(as.character(num_levels))) {
          warning_msg_temp <- c(warning_msg_temp, paste0("number of levels is missing"))
        }
      } else {
        num_levels <- set_num_levels_reactive()
      }
      
      # Checks the user has provided level values, that they are in the appropriate format (numeric for continuous), and that the appropriate number has been provided
      level_values <- current_factor_data$`Level Values`[i]
      if (is.na(level_values) || grepl("^\\s*$", level_values) || level_values == "") {
        warning_msg_temp <- c(warning_msg_temp, paste0("level values are missing"))
      } else {
        level_values_vector <- gsub("^\\s+|\\s+$", "", unlist(strsplit(level_values, ",")))
        if (!is.null(num_levels)) {
          if (length(level_values_vector) != num_levels) {
            warning_msg_temp <- c(warning_msg_temp, paste0("the number of level values provided does not match the number of levels selected or required by the design"))
          }
        }
        
        if (factor_type == "Continuous") {
          if (!all(!is.na(as.numeric(level_values_vector)))) {
            warning_msg_temp <- c(warning_msg_temp, paste0("the level values for a Continuous factor must be numeric"))
          }
        }
      }
      
      # Combines all of the warning messages per row
      if (length(warning_msg_temp) != 0) {
        warning_msg_temp <- paste0("Row ", i, ": ", paste0(warning_msg_temp, collapse = "; "), ".")
        warning_msg <- c(warning_msg, warning_msg_temp)
      }
    }

    # If there are warning messages (i.e. key info is missing), print the warning messages and stop the user from proceeding
    if (length(warning_msg) != 0) {
      output$warning_message_id <- renderUI({
        lapply(seq_len(length(warning_msg)), function(j) {
          p(strong(warning_msg[j]), class = "text-red")
        })
      })
    }
    req(length(warning_msg) == 0, cancelOutput = TRUE)
    
    # If all the info is inputted, close the model and start the blocking process
    removeModal()
    
    n_factors <- dim(current_factor_data)[1]

    # Open a modal with the options for blocking and randomization
    showModal(modalDialog(
      title = h3(strong("Design Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
      fluidRow(
        column(
          6,
          p(strong("Factor Details: ")),
          DT::dataTableOutput("summary_factors_table_id")
        ),
        column(
          6,
          uiOutput("augment_design_selections_id")
        )
      ),
      easyClose = FALSE,
      footer = tagList(
        fluidRow(
          column(offset = 8,
                 2,
                 actionButton(inputId = "augment_design_id", "Update Design", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")
          ),
          column(
            2,
            actionButton(inputId = "close_modal_id", strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;")
          )
        )
      )
    ))
    
    # Show the summary of the factors
    output$summary_factors_table_id <- DT::renderDataTable(
      current_factor_data,
      escape = FALSE, selection = "none", server = FALSE,
      extensions = c("FixedColumns", "FixedHeader", "Scroller"),
      options = list(
        dom = "t",
        pageLength = 10, paging = TRUE,
        headerCallback = DT::JS("function(thead) {", "  $(thead).css('font-size', '1em');", "}")
      ),
      rownames = FALSE, filter = "none"
    )

    n_factors <- dim(current_factor_data)[1]
    factor_names <- as.character(current_factor_data$`Factor Names`)
    factor_types <- as.character(current_factor_data$`Factor Type`)
    level_values <- as.character(current_factor_data$`Level Values`)
    level_values <- sapply(strsplit(level_values, ","), function(x) length(x))
    
    output$augment_design_selections_id <- renderUI({
      fluidPage(
        fluidRow(
          column(
            12,
            p("To augment the design, select the effects you
                wish to estimate. Please use the interface below to select
                all relevant effects, including main effects, interactions, quadratic
                terms, and any other higher-order effects.")
          )
        ),
        fluidRow(
          column(
            4,
            actionButton(inputId = "Augment_main_effects_add_btn_id", label = "Main Effects", width = "100%", class = "Blue-Button")
          ),
          column(
            4,
            actionButton(inputId = "Augment_2way_effects_add_btn_id", label = "2-Way Interactions", width = "100%", class = "Blue-Button")
          ),
          column(
            4,
            actionButton(inputId = "Augment_quadratics_add_btn_id", label = "Quadratics", width = "100%", class = "Blue-Button")
          )
        ),
        div(style = "height: 10px;"),
        fluidRow(
          column(
            4,
            actionButton(inputId = "Augment_add_btn_id", label = "Add", width = "100%", class = "Blue-Button")
          ),
          column(
            4,
            actionButton(inputId = "Augment_cross_btn_id", label = "Create Interaction", width = "100%", class = "Blue-Button")
          ),
          column(
            4,
            actionButton(inputId = "Augment_remove_btn_id", label = "Remove", width = "100%", class = "Blue-Button")
          )
        ),
        div(style = "height: 10px;"),
        fluidRow(
          column(
            4,
            selectInput(inputId = "available_factors_augment_id", label = "Factor:", choices = factor_names, size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
          ),
          column(
            8,
            selectInput(inputId = "selected_factors_augment_id", label = "Selected Effects:", choices = NULL, size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
          )
        ),
        fluidRow(
          column(
            4,
            numericInput(inputId = "Augment_num_runs_id", label = "Number of Runs:", min = 1, max = Inf, value = NULL, step = 1, width = "100%")
          ),
          column(
            3,
            uiOutput("Augment_min_runs_ui_id")
          ),
          column(
            3,
            uiOutput("Augment_max_runs_ui_id")
          )
        ),
        fluidRow(
          column(
            12,
            materialSwitch(inputId = "Augment_randomize_checkbox_id", label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        )
      )
    })
  })
  
  observeEvent(input$Augment_main_effects_add_btn_id, {

    current_factor_data <- factor_data_reactive()
    factor_names <- as.character(current_factor_data$`Factor Names`)
    currently_selected_factors <- augment_design_model_components()
    
    updated_selected_factors <- unique(c(currently_selected_factors, factor_names))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    augment_design_model_components(sorted_terms)

    updateSelectInput(session, inputId = "selected_factors_augment_id", choices = sorted_terms)
  })
  
  observeEvent(input$Augment_2way_effects_add_btn_id, {
    current_factor_data <- factor_data_reactive()
    factor_names <- as.character(current_factor_data$`Factor Names`)
    req(length(factor_names) > 1, cancelOutput = TRUE)
    
    currently_selected_factors <- augment_design_model_components()
    
    interactions <- combn(factor_names, 2, FUN = function(x) paste(x, collapse = ":"))
    
    updated_selected_factors <- unique(c(currently_selected_factors, interactions))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    augment_design_model_components(sorted_terms)
    
    updateSelectInput(session, inputId = "selected_factors_augment_id", choices = sorted_terms)
  })
  
  observeEvent(input$Augment_quadratics_add_btn_id, {
    current_factor_data <- factor_data_reactive()
    factor_names <- as.character(current_factor_data$`Factor Names`)
    
    currently_selected_factors <- augment_design_model_components()
    
    quadratics <- paste(factor_names, factor_names, sep = ":")
    
    updated_selected_factors <- unique(c(currently_selected_factors, quadratics))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    augment_design_model_components(sorted_terms)
    
    updateSelectInput(session, inputId = "selected_factors_augment_id", choices = sorted_terms)
  })
  
  observeEvent(input$Augment_add_btn_id, {
    selected_available_factors <- input$available_factors_augment_id
    
    if (is.null(selected_available_factors)) {
      showNotification("Select a factor to add", duration = 20, closeButton = TRUE, type = "warning")
    }
    req(!(is.null(selected_available_factors)), cancelOutput = TRUE)
    
    currently_selected_factors <- augment_design_model_components()
    
    updated_selected_factors <- unique(c(currently_selected_factors, selected_available_factors))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    augment_design_model_components(sorted_terms)
    
    updateSelectInput(session, inputId = "selected_factors_augment_id", choices = sorted_terms)
  })
  
  observeEvent(input$Augment_cross_btn_id, {
    currently_selected_factors <- augment_design_model_components()
    
    selected_available_factors <- input$available_factors_augment_id
    selected_model_effects <- input$selected_factors_augment_id
    
    if (is.null(selected_available_factors) || is.null(selected_model_effects)) {
      showNotification("Select at least one factor and one effect to create an interaction.", duration = 20, closeButton = TRUE, type = "warning")
    }
    req(!is.null(selected_available_factors), !is.null(selected_model_effects), cancelOutput = TRUE)
    
    combinations <- expand.grid(selected_available_factors, selected_model_effects)
    
    interactions <- paste(combinations[, 1], combinations[, 2], sep = ":")
    
    updated_selected_factors <- unique(c(currently_selected_factors, interactions))
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    augment_design_model_components(sorted_terms)
    
    updateSelectInput(session, inputId = "selected_factors_augment_id", choices = sorted_terms)
  })
  
  observeEvent(input$Augment_remove_btn_id, {
    currently_selected_factors <- augment_design_model_components()
    
    selected_model_effects <- input$selected_factors_augment_id
    if (is.null(selected_model_effects)) {
      showNotification("Select an effect to remove", duration = 20, closeButton = TRUE, type = "warning")
    }
    req(!(is.null(selected_model_effects)), cancelOutput = TRUE)
    
    updated_selected_factors <- setdiff(currently_selected_factors, selected_model_effects)
    term_types <- sapply(updated_selected_factors, classify_terms)
    terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
    sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
    sorted_terms <- sorted_terms_df$term
    augment_design_model_components(sorted_terms)
    
    updateSelectInput(session, inputId = "selected_factors_augment_id", choices = sorted_terms)
  })
  
  Augment_min_runs_reactive <- reactiveVal(NULL)
  Augment_max_runs_reactive <- reactiveVal(NULL)
  observeEvent(augment_design_model_components(), {
    req(augment_design_model_components(), cancelOutput = TRUE)
    
    components <- augment_design_model_components()
    
    if (length(components) == 0) {
      output$Augment_min_runs_ui_id <- renderUI(p(strong("Min Runs: NA")))
      output$Augment_max_runs_ui_id <- renderUI(p(strong("Max Runs: NA")))
      Augment_min_runs_reactive(NULL)
      Augment_max_runs_reactive(NULL)
      return(NULL)
    }

    # ============================================================
    # Load factor information (current)
    # ============================================================
    factor_data <- factor_data_reactive()
    
    factor_info <- lapply(seq_len(nrow(factor_data)), function(i) {
      fname  <- factor_data$`Factor Names`[i]
      ftype  <- factor_data$`Factor Type`[i]
      levels <- trimws(strsplit(factor_data$`Level Values`[i], ",")[[1]])
      list(name = fname, type = ftype, levels = levels)
    })
    names(factor_info) <- factor_data$`Factor Names`
    
    factor_names <- names(factor_info)

    # ============================================================
    # Build formula for user-selected model
    # ============================================================
    fix_component <- function(term) {
      parts <- strsplit(term, ":")[[1]]
      uniq  <- unique(parts)
      
      # pure polynomial (X:X, X:X:X, ...)
      if (length(uniq) == 1 && length(parts) > 1) {
        fname <- uniq[1]
        deg   <- length(parts)
        return(paste0("I(", fname, "^", deg, ")"))
      }
      
      # leave interaction terms intact
      term
    }
    
    components_fixed <- sapply(components, fix_component)
    formula <- as.formula(paste("~", paste(components_fixed, collapse = " + ")))
    
    # ============================================================
    # STEP 1: Build a full candidate set → compute p correctly
    # ============================================================
    
    factor_degree <- function(fname) {
      deg <- 1
      for (term in components) {
        parts <- strsplit(term, ":")[[1]]
        count <- sum(parts == fname)
        if (count > deg) deg <- count
      }
      deg
    }

    candidate_list <- lapply(factor_names, function(f) {
      info <- factor_info[[f]]
      
      if (info$type == "Continuous") {
        deg <- factor_degree(f)
        seq(-1, 1, length.out = deg + 1)
        
      } else if (info$type == "Discrete") {
        vals <- as.numeric(info$levels)
        sort(unique(vals))
        
      } else if (info$type == "Categorical") {
        factor(info$levels, levels = info$levels)
        
      } else {
        stop("Unknown factor type for ", f)
      }
    })

    candidate_set <- do.call(expand.grid, candidate_list)
    colnames(candidate_set) <- factor_names
    
    X_full <- model.matrix(formula, data = candidate_set)
    p <- ncol(X_full)   # total model parameters
    
    # ============================================================
    # STEP 2: Build padded existing design → compute rank correctly
    # ============================================================

    existing_full <- design_reactive()
    
    # Which factors are already in the existing design?
    existing_cols <- intersect(factor_names, colnames(existing_full))

    new_cols      <- setdiff(factor_names, colnames(existing_full))

    if (length(existing_cols) == 0) {
      # No existing factor is present → rank = 0
      rank_exist <- 0
      
    } else {
      
      base_existing <- existing_full[, existing_cols, drop = FALSE]

      # Build padded design frame
      typed_existing <- data.frame(
        row_id = seq_len(nrow(base_existing)),
        stringsAsFactors = FALSE
      )
      
      for (f in factor_names) {
        info <- factor_info[[f]]
        
        if (f %in% existing_cols) {
          # Convert existing factor to correct type
          col <- base_existing[[f]]
          if (info$type %in% c("Continuous", "Discrete")) {
            typed_existing[[f]] <- as.numeric(col)
          } else {
            typed_existing[[f]] <- factor(as.character(col),
                                          levels = info$levels)
          }
          
        } else {
          # Newly added factor → pad with constants
          if (info$type %in% c("Continuous", "Discrete")) {
            typed_existing[[f]] <- rep(0, nrow(base_existing))
          } else {
            typed_existing[[f]] <- factor(rep(info$levels[1], nrow(base_existing)),
                                          levels = info$levels)
          }
        }
      }
      
      X_exist <- model.matrix(formula, data = typed_existing)
      rank_exist <- qr(X_exist)$rank
    }

    # ============================================================
    # STEP 3: Runs needed and maximum runs
    # ============================================================
    
    required_new_runs <- max(0, p - rank_exist)
    
    max_runs <- prod(sapply(factor_names, function(f) {
      info <- factor_info[[f]]
      if (info$type == "Continuous") {
        deg <- factor_degree(f)
        deg + 1
      } else {
        length(info$levels)
      }
    }))

    max_additional_runs = max_runs - rank_exist
    
    if (required_new_runs > max_additional_runs) {
      required_new_runs <- max_additional_runs
    }
    

    
    # ============================================================
    # Update UI outputs and reactives
    # ============================================================
    
    Augment_min_runs_reactive(required_new_runs)
    Augment_max_runs_reactive(max_additional_runs)
    
    output$Augment_min_runs_ui_id <- renderUI(
      p(strong(paste0("Min Additional Runs: ", required_new_runs)))
    )
    
    output$Augment_max_runs_ui_id <- renderUI(
      p(strong(paste0("Max Additional Runs: ", max_additional_runs)))
    )
  })
  
  observeEvent(input$augment_design_id, {
    # Get the factor data, the selected design, and if the runs should be randomized
    factor_data <- factor_data_reactive()
    original_factor_types <- factor_types_reactive()
    n_factors <- dim(factor_data)[1]
    factor_names <- factor_data$`Factor Names`
    
    randomize <- input$Augment_randomize_checkbox_id
    
    # Helper function to check that the number of runs selected by the
    # user is in the right format for LHS and D-Opt
    validate_num_runs <- function(num_runs) {
      if (is.character(num_runs) || is.factor(num_runs)) {
        suppressWarnings(num_runs <- as.numeric(as.character(num_runs)))
        return(FALSE)
      }
      
      # Check it's numeric
      if (!is.numeric(num_runs) || length(num_runs) != 1 || is.na(num_runs)) {
        showNotification("The 'Number of Experimental Runs' must be a single numeric integer.", duration = 20, closeButton = TRUE)
        return(FALSE)
      }
      
      # Check it's an integer (within floating-point tolerance)
      if (abs(num_runs - round(num_runs)) > .Machine$double.eps^0.5) {
        showNotification("The 'Number of Experimental Runs' must be an integer (no decimal values).", duration = 20, closeButton = TRUE)
        return(FALSE)
      }
      
      # Check it's positive
      if (num_runs <= 0) {
        showNotification("The 'Number of Experimental Runs' must be greater than zero.")
        return(FALSE)
      }
      
      # check it is within the min and max number of runs for a DOD
      if (!(num_runs >= Augment_min_runs_reactive() && num_runs <= Augment_max_runs_reactive())) {
        showNotification("The 'Number of Experimental Runs' must be within the minumum and maximum number of runs provided.")
        return(FALSE) 
      }
      TRUE
    }
    

    # Creates the design based on the users choices
    model_components <- augment_design_model_components()
    
    existing_design_full <- design_reactive()
    existing_design <- existing_design_full[, colnames(existing_design_full) %in% factor_names, drop = FALSE]
    
    other_columns <- existing_design_full[, !(colnames(existing_design_full) %in% factor_names), drop = FALSE]
    
    
    req(validate_num_runs(input$Augment_num_runs_id), cancelOutput = TRUE)
    num_runs <- as.integer(round(input$Augment_num_runs_id))
    
    factor_types <- setNames(as.character(factor_data$`Factor Type`), factor_names)
    
    level_values <- setNames(
      lapply(strsplit(as.character(factor_data$`Level Values`), ","), function(x) trimws(x)),
      factor_names
    )
    
    design <- D_Optimal_Augment(existing_design, level_values, model_components, factor_types, num_runs, randomize)
    
    req(!is.null(design))
    
    # Add in any existing non-factor roles:
    n_old <- nrow(other_columns)
    n_new <- nrow(design) - n_old
    
    filler <- as.data.frame(
      matrix("",
             nrow = n_new,
             ncol = ncol(other_columns)
      ),
      stringsAsFactors = FALSE
    )
    
    colnames(filler) <- colnames(other_columns)
    
    other_extended <- rbind(other_columns, filler)
    design <- cbind(design, other_extended)
    
    other_factor_types <- original_factor_types[colnames(other_columns)]
    
    factor_types <- c(factor_types, other_factor_types)
    factor_types_df <- as.data.frame(factor_types)
    colnames(factor_types_df) <- colnames(design)
    
    for (col_name in colnames(design)) {
      col_type <- factor_types_df[[col_name]]
      if (col_type %in% c("Continuous", "Discrete", "Response")) {
        design[[col_name]] <- as.numeric(design[[col_name]])
      } else if (col_type %in% c("Block", "Categorical")) {
        design[[col_name]] <- as.character(design[[col_name]])
      }
    }
    
    design_reactive(design)
    
    factor_types_reactive(factor_types_df)
    
    output$design_table_id <- renderRHandsontable({
      rhandsontable(design_reactive(), selectCallback = TRUE, width = "1050px", height = "350px", stretchH = "all") %>%
        hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
    })
    
    removeModal()

  })
  
  # ===== Analyse the Design and Model Explorer ================================
  output$model_visualizer_output_ui_id <- renderUI(p("Create a model in ", strong("Analyze the Design"), " to initialize the visualizer."))
  
  factor_names_reactive <- reactiveVal(value = NULL) # Stores names of Factors (excl blocks)
  factor_blocks_reactive <- reactiveVal(value = NULL) # Stores names of Blocks
  responses_reactive <- reactiveVal(value = NULL) # Stores names of Responses
  
  design_matrix_reactive <- reactiveVal(NULL) # Stores the design matrix
  model_formula_reactive <- reactiveVal(NULL) # stores the model formula
  model_aov_reactive <- reactiveVal(NULL) # stores the model AoV
  
  selected_factors_for_model_reactive <- reactiveVal(value = NULL) # Store the effects selected by the user for the model


  # Triggers when the user creates, loads, or edits the design
  observeEvent(design_reactive(), {
    if (!(is.null(design_reactive()))) {
      design <- design_reactive()
      factor_types <- factor_types_reactive()

      factor_names <- names(factor_types)[sapply(factor_types, function(x) !x %in% c("Response", "Block", "Metadata"))]
      factor_blocks <- names(factor_types)[sapply(factor_types, function(x) x == "Block")]
      responses <- names(factor_types)[sapply(factor_types, function(x) x == "Response")]

      factor_names_reactive(factor_names)
      factor_blocks_reactive(factor_blocks)
      responses_reactive(responses)

      # Once a design exists, it renders the make a design UI
      output$analysis_table_ui_id <- renderUI({
        div(
          p("Choose a default model from the dropdown.  Factors to include in the model can be manually adjusted as needed.  "),
          fluidRow(
            column(
              4,
              selectInput(inputId = "select_model_id", label = "Select Default Model (optional):", choices = c("Main Effects", "Main Effects + Interactions", "Response Surface"), selected = NULL, multiple = FALSE, width = "100%")
            ),
            column(
              2,
              actionButton(inputId = "apply_selected_model_id", label = "Apply", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 20px;")
            )
          ),
          fluidRow(
            column(
              3,
              selectInput(inputId = "available_factors_for_model_id", label = "Select Factor(s):", choices = c(factor_names, factor_blocks), size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
            ),
            column(
              1,
              div(style = "height: 25px;"),
              fluidRow(
                actionButton(inputId = "add_factor_to_model_id", label = "Add", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
              ),
              fluidRow(
                actionButton(inputId = "remove_factor_to_model_id", label = "Remove", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
              ),
              fluidRow(
                actionButton(inputId = "cross_factors_for_model_id", label = "Cross", width = "100%", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
              )
            ),
            column(
              3,
              selectInput(inputId = "selected_factors_for_model_id", label = "Factor(s) for Model:", choices = selected_factors_for_model_reactive(), size = 8, selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%")
            ),
            column(
              3,
              selectInput(inputId = "selected_responses_for_model_id", label = "Select Response:", choices = responses, selected = NULL, multiple = FALSE, width = "100%"),
            )
          ),
          fluidRow(
            column(4,
              offset = 8,
              div(
                align = "right",
                actionButton(inputId = "fit_model_button_id", label = "Analyze", width = "200px", class = "Blue-Button", style = "margin-bottom: 0px; margin-top: 5px;")
              )
            )
          )
        )
      })
    }
  })

  # Apply the model selected by the user (main factors, response surface, etc.)
  observeEvent(input$apply_selected_model_id, {
    selected_model <- input$select_model_id

    factor_names <- factor_names_reactive()
    factor_blocks <- factor_blocks_reactive()
    responses <- responses_reactive()

    current_selected_model <- input$select_model_id
    current_selected_response <- input$selected_responses_for_model_id

    if (selected_model == "Main Effects") {
      main_effects <- c(factor_names, factor_blocks)
      selected_factors_for_model_reactive(main_effects)
      updateSelectInput(session, "select_model_id", selected = current_selected_model)
      updateSelectInput(session, "selected_responses_for_model_id", selected = current_selected_response)
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    } else if (selected_model == "Main Effects + Interactions") {
      main_effects <- c(factor_names, factor_blocks)
      interactions <- combn(factor_names, 2, FUN = paste0, collapse = " * ")

      all_effects <- c(main_effects, interactions)
      selected_factors_for_model_reactive(all_effects)
      updateSelectInput(session, "select_model_id", selected = current_selected_model)
      updateSelectInput(session, "selected_responses_for_model_id", selected = current_selected_response)
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    } else if (selected_model == "Response Surface") {
      main_effects <- c(factor_names, factor_blocks)
      interactions <- combn(factor_names, 2, FUN = paste0, collapse = " * ")
      quadratics <- unname(sapply(c(factor_names), function(x) paste(x, x, sep = " * ")))

      all_effects <- c(main_effects, interactions, quadratics)
      selected_factors_for_model_reactive(all_effects)
      updateSelectInput(session, "select_model_id", selected = current_selected_model)
      updateSelectInput(session, "selected_responses_for_model_id", selected = current_selected_response)
      updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
    }

  })

  # The use can manually create a model by adding, crossing, and/or removing factors
  # Add a factor to the model
  observeEvent(input$add_factor_to_model_id, {
    factor_to_add <- input$available_factors_for_model_id
    selected_factors <- selected_factors_for_model_reactive()

    update_selected_factors <- unique(c(selected_factors, factor_to_add))
    selected_factors_for_model_reactive(update_selected_factors)

    current_selected_model <- input$select_model_id
    current_selected_response <- input$selected_responses_for_model_id
    updateSelectInput(session, "select_model_id", selected = current_selected_model)
    updateSelectInput(session, "selected_responses_for_model_id", selected = current_selected_response)
    updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
  })
  # Remove a factor or interactions
  observeEvent(input$remove_factor_to_model_id, {
    factor_to_remove <- input$selected_factors_for_model_id
    selected_factors <- selected_factors_for_model_reactive()

    remaining_factors <- setdiff(selected_factors, factor_to_remove)
    selected_factors_for_model_reactive(remaining_factors)

    current_selected_model <- input$select_model_id
    current_selected_response <- input$selected_responses_for_model_id
    updateSelectInput(session, "select_model_id", selected = current_selected_model)
    updateSelectInput(session, "selected_responses_for_model_id", selected = current_selected_response)
    updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
  })
  # Create an interaction by crossing factors (works for quadratics too)
  observeEvent(input$cross_factors_for_model_id, {
    selected_factors_left <- input$available_factors_for_model_id

    factor_names <- factor_names_reactive()
    factor_blocks <- factor_blocks_reactive()
    responses <- responses_reactive()

    if (any(selected_factors_left %in% factor_blocks)) {
      showNotification(strong("Interactions involving blocks are not currently supported."), duration = 15, closeButton = TRUE, type = "error")
    }
    req(!any(selected_factors_left %in% factor_blocks), cancelOutput = TRUE)

    selected_factors_right <- input$selected_factors_for_model_id

    factor_combinations <- expand.grid(selected_factors_left, selected_factors_right)
    crossed_factors <- apply(factor_combinations, 1, function(row) paste(row, collapse = " * "))

    existing_selected_factors <- selected_factors_for_model_reactive()
    all_selected_factors <- unique(c(existing_selected_factors, crossed_factors))
    selected_factors_for_model_reactive(all_selected_factors)

    current_selected_model <- input$select_model_id
    current_selected_response <- input$selected_responses_for_model_id
    updateSelectInput(session, "select_model_id", selected = current_selected_model)
    updateSelectInput(session, "selected_responses_for_model_id", selected = current_selected_response)
    updateSelectInput(session, "selected_factors_for_model_id", choices = selected_factors_for_model_reactive())
  })


  # Fit the model, present the results, prepare the ANOVA, Design Diagnostics,
  # the Visualizer and the Interaction Tabs
  observeEvent(input$fit_model_button_id, {
    # ===== Fit the Model ======================================================
    # Get the user inputs
    selected_factors_all <- selected_factors_for_model_reactive()
    selected_response <- input$selected_responses_for_model_id

    # Get the design and factor details
    design <- design_reactive()
    factor_types <- factor_types_reactive()

    # Prep the design matrix based on the factor types (e.g. categorical )
    design_matrix <- design
    for (col in colnames(design_matrix)) {
      role <- factor_types[[col]]
      if (role == "Categorical" || role == "Block") {
        design_matrix[[col]] <- as.factor(design_matrix[[col]])
      } else if (role == "Continuous") {
        design_matrix[[col]] <- as.numeric(design_matrix[[col]])
      }
    }
    
    design_matrix_reactive(design_matrix)

    # Function to change the formatting of:
    # quadratics, from A * A to I(A)^2
    transform_term <- function(term) {
      if (grepl("\\*", term)) {
        vars <- unlist(strsplit(term, "\\s*\\*\\s*"))
        if (length(vars) == 2 && vars[1] == vars[2]) {
          # Quadratic term like "A * A"
          paste0("I(", vars[1], "^2)")
        } else {
          # Interaction term like "A * B"
          paste0(vars[1], ":", vars[2])
        }
      } else {
        term
      }
    }

    transformed_factors <- sapply(selected_factors_all, transform_term)

    # Prepare the model formula based on the selected response and factors and do the ANOVA
    model_formula <- as.formula(paste0(selected_response, " ~ ", paste0(transformed_factors, collapse = " + ")))
    model_formula_reactive(model_formula)
    
    model_aov <- aov(model_formula, design_matrix)
    model_aov_reactive(model_aov)
    
    anova_results <- broom::tidy(model_aov)
    colnames(anova_results) <- c("Factor", "Degrees of Freedom", "Sum of Squares", "Mean Square", "F-Statistic", "P-Value")


    # Get the design correlations (except if there is only 1 factor in the model)
    design_correlations <- model.matrix(model_formula, data = design_matrix)
    colnames(design_correlations) <- update_factor_names(colnames(design_correlations), design_matrix)
    design_correlations <- design_correlations[, colnames(design_correlations) != "(Intercept)"]
    if (is.null(dim(design_correlations))) {
      "N/A"
    } else {
      model_correlations <- cor(design_correlations, method = c("pearson"))
      cor_data <- as.data.frame(as.table(model_correlations))
    }

    # Get the linear model summary from the ANOVA (Coefficients for each term, standard error, etc.)
    estimates <- summary.lm(model_aov)
    estimates_df <- broom::tidy(estimates)
    colnames(estimates_df) <- c("Factor", "Estimate", "Standard Error", "t-Statistic", "P-Value")
    estimates_df$Factor <- gsub("\\(Intercept\\)", "Intercept", estimates_df$Factor)
    estimates_df$Factor <- update_factor_names(estimates_df$Factor, design_matrix)

    # Get the alias, if Null, i.e. nothing is aliased, state "model terms are not aliased"
    alias_complete <- alias(model_aov)$Complete
    if (is.null(alias_complete)) {
      alias_structure <- "Model terms are not aliased."
    } else {
      alias_structure <- data.frame(alias_complete)
      colnames(alias_structure) <- update_factor_names(colnames(alias_structure), design_matrix)
      rownames(alias_structure) <- update_factor_names(rownames(alias_structure), design_matrix)
      names(alias_structure) <- gsub("X.Intercept.", "Intercept", names(alias_structure))
    }

    # extract/calculate and round the goodness of fit metrics, R-squared, adjusted_rsq,
    # model f-statistic, model p-value
    rsq <- signif(estimates$r.squared, 3)
    adjusted_rsq <- signif(estimates$adj.r.squared, 3)
    f_statistic <- signif(unname(estimates$fstatistic[1]), 3)
    model_p_value <- signif(pf(estimates$fstatistic[1], estimates$fstatistic[2], estimates$fstatistic[3], lower.tail = FALSE), 3)

    # Extract the actual values into a dataframe and add the predicted values as a column
    all_data <- data.frame("actual_values" = design_matrix[[selected_response]])
    all_data$predicted_values <- predict(model_aov)


    # ===== Present the Model Results ==========================================
    # Set up the ui elements for presenting the results, include the goodness of fit metrics
    output$ANOVA_output_ui_id <- renderUI({
      div(
        fluidRow(
          column(
            6,
            h4(strong("Predicted vs. Actual")),
            box(
              id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
              plotOutput("pred_vs_actual_output_id")
            )
          ),
          column(
            6,
            h4(strong("Model Fit Metrics")),
            box(
              id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
              p(HTML(paste0(strong("R-Squared: "), rsq))),
              p(HTML(paste0(strong("Adjusted R-Squared: "), adjusted_rsq))),
              p(HTML(paste0(strong("F-statistic: "), f_statistic, ", degrees of freedom: ", unname(estimates$fstatistic[2]), ", ", unname(estimates$fstatistic[3])))),
              p(HTML(paste0(strong("p-value: "), model_p_value)))
            )
          )
        ),
        fluidRow(
          column(
            12,
            h4(strong("Model Summary")),
            box(
              id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
              tableOutput("ANOVA_results_output_id")
            )
          ),
          column(
            12,
            h4(strong("Estimates")),
            box(
              id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
              tableOutput("Estimates_output_id")
            )
          )
        )
      )
    })

    # Print the parity plot
    min_value <- as.numeric(min(c(all_data$actual_values, all_data$predicted_values), na.rm = TRUE))
    max_value <- as.numeric(max(c(all_data$actual_values, all_data$predicted_values), na.rm = TRUE))
    output$pred_vs_actual_output_id <- renderPlot({
      p <- ggplot(all_data, aes(x = as.numeric(actual_values), y = as.numeric(predicted_values))) +
        geom_point(color = "#0097a9") + # Add points
        geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black") + # x = y line
        labs(x = "Actual", y = "Predicted", title = "Parity Plot") +
        xlim(min_value, max_value) + # Set x-axis limits
        ylim(min_value, max_value) + # Set y-axis limits
        coord_fixed(ratio = 1)

      p
    })

    # render the ANOVA results in a table
    output$ANOVA_results_output_id <- renderTable({
      col_names <- colnames(anova_results)
      data_rounded <- as.data.frame(lapply(anova_results, function(x) {
        if (is.numeric(x)) {
          as.character({
            signif(x, digits = 4)
          })
        } else {
          x
        } # Leave non-numeric columns unchanged
      }))

      data_rounded$significance <- sapply(data_rounded$`P.Value`, function(p) {
        if (!is.na(p) && as.numeric(p) < 0.001) {
          "***"
        } else if (!is.na(p) && as.numeric(p) < 0.01) {
          "**"
        } else if (!is.na(p) && as.numeric(p) < 0.05) {
          "*"
        } else {
          "" # Leave blank if none of the conditions are met or p is NA
        }
      })

      data_rounded$`P.Value` <- sapply(data_rounded$`P.Value`, function(p) {
        p_val <- as.numeric(p)
        if (!is.na(p_val)) {
          if (p_val < 0.001) {
            "<0.001"
          } else {
            as.character(signif(p_val, digits = 4))
          }
        } else {
          p
        }
      })

      colnames(data_rounded) <- c(col_names, "")
      data_rounded
    })

    # render the estimates in a table
    output$Estimates_output_id <- renderTable({
      col_names <- colnames(estimates_df)
      data_rounded <- as.data.frame(lapply(estimates_df, function(x) {
        if (is.numeric(x)) {
          as.character({
            signif(x, digits = 4)
          })
        } else {
          x
        } # Leave non-numeric columns unchanged
      }))

      data_rounded$significance <- sapply(data_rounded$`P.Value`, function(p) {
        if (!is.na(p) && as.numeric(p) < 0.001) {
          "***"
        } else if (!is.na(p) && as.numeric(p) < 0.01) {
          "**"
        } else if (!is.na(p) && as.numeric(p) < 0.05) {
          "*"
        } else {
          "" # Leave blank if none of the conditions are met or p is NA
        }
      })

      data_rounded$`P.Value` <- sapply(data_rounded$`P.Value`, function(p) {
        p_val <- as.numeric(p)
        if (!is.na(p_val)) {
          if (p_val < 0.001) {
            "<0.001"
          } else {
            as.character(signif(p_val, digits = 4))
          }
        } else {
          p
        }
      })

      colnames(data_rounded) <- c(col_names, "")
      data_rounded
    })

    # Go to the ANOVA tab to show the results
    updateTabsetPanel(session, "Analysis_results_tabs_id", selected = "Model Fitting Results")


    # ===== Design Diagnostics Tab =============================================
    # Prep the ui elements for the design diagnostics tab
    output$Design_diagnostics_output_ui_id <- renderUI({
      div(
        fluidRow(
          column(
            12,
            h3(strong("Alias Structure:")),
            uiOutput("alias_table_output_id")
          )
        ),
        fluidRow(
          column(
            12,
            h3(strong("Correlations:")),
            uiOutput("correlations_output_ui_id")
          )
        )
      )
    })

    # adjust the alias output if there are no aliased terms
    output$alias_table_output_id <- renderUI({
      if (is.data.frame(alias_structure)) {
        tableOutput("alias_table")
      } else {
        verbatimTextOutput("alias_text")
      }
    })

    # render the alias table or text as appropriate
    output$alias_table <- renderTable(
      {
        alias_structure
      },
      rownames = TRUE
    )
    output$alias_text <- renderText({
      as.character(alias_structure)
    })

    # adjust the correlations output if there are is only 1 factor in the model
    output$correlations_output_ui_id <- renderUI({
      if (is.null(dim(design_correlations))) {
        verbatimTextOutput("correlations_text_output_id")
      } else {
        plotlyOutput("correlations_plot_output_id")
      }
    })

    # render the correlations as text or a plot as appropriate
    output$correlations_text_output_id <- renderText({
      "Correlation matrix not available; model contains only one factor."
    })
    output$correlations_plot_output_id <- renderPlotly({
      plot_ly(
        data = cor_data,
        x = ~Var1, # Assuming Var1 and Var2 are the variable names in the matrix
        y = ~Var2,
        z = ~Freq, # Frequency here represents the correlation values
        type = "heatmap",
        colors = colorRamp(c("blue", "white", "red")),
        text = ~ paste("Correlation: ", round(Freq, 2)), # Text to show on hover
        hoverinfo = "text+x+y", # Custom hover info
        zmin = -1,
        zmax = 1,
        width = 600, # Width of the plot
        height = 600 # Height of the plot, same as width to ensure square
      ) %>% layout(
        xaxis = list(title = "", tickangle = -45, autorange = TRUE, side = "top"),
        yaxis = list(title = "", autorange = "reversed", scaleanchor = "x", scaleratio = 1), # Enforcing the square aspect ratio
        autosize = FALSE
      )
    })


    # ===== Model Visualizer ===================================================
    # Preps the UI for the visualizer
    output$model_visualizer_output_ui_id <- renderUI({
      # Ensure that model_formula, design_matrix, and model_aov exist before proceeding
      model_formula = model_formula_reactive()
      design_matrix = design_matrix_reactive()
      model_aov = model_aov_reactive()
      
      req(model_formula, design_matrix, model_aov)

      # Extract the names of all terms (main effects + interactions) from the
      # model formula and filter to identify the unique factors included in the model
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]

      n <- length(main_effects)
      if (n == 0) {
        return(h4("No main effect terms found in model formula."))
      }

      # Set the grid layout based on the number of effects
      ncol <- ceiling(sqrt(n))
      nrow <- ceiling(n / ncol)

      # Create a plot with corresponding slider bar for each factor in the model
      # The format of the plot and slider depend on the factor type
      plot_slider_grid <- lapply(seq_along(main_effects), function(i) {
        factor_name <- main_effects[i]
        column_data <- design_matrix[[factor_name]]
        slider_id <- paste0("slider_", factor_name)
        factor_type <- factor_types[[factor_name]]

        slider_ui <- switch(factor_type,
          "Continuous" = sliderInput(
            inputId = slider_id,
            label = NULL,
            min = min(column_data, na.rm = TRUE),
            max = max(column_data, na.rm = TRUE),
            value = mean(column_data, na.rm = TRUE),
            step = (max(column_data, na.rm = TRUE) - min(column_data, na.rm = TRUE)) / 100
          ),
          "Discrete" = sliderTextInput(
            inputId = slider_id,
            label = NULL,
            choices = sort(unique(column_data)),
            selected = unique(column_data)[1],
            grid = TRUE
          ),
          "Categorical" = sliderTextInput(
            inputId = slider_id,
            label = NULL,
            choices = levels(column_data),
            selected = levels(column_data)[1],
            grid = TRUE
          ),
          h5(paste("Unknown factor type for", factor_name))
        )

        column(
          width = 12 / ncol,
          plotOutput(paste0("factor_plot_", factor_name), height = "250px"),
          slider_ui
        )
      })

      # render the plots
      tagList(
        p("Use the sliders below to adjust the values of each factor.
           The plots will update automatically, and the Predicted Response
           will reflect your current selections."),
        h4("Predicted Response:"),
        hr(),
        uiOutput("selected_value_output"),
        fluidRow(plot_slider_grid)
      )
    })

    # Get the current values from the sliders as a dataframe for prediction
    new_data <- reactive({
      model_formula = model_formula_reactive()
      design_matrix = design_matrix_reactive()
      req(model_formula, design_matrix)

      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]

      data_list <- list()
      for (factor_name in main_effects) {
        column_data <- design_matrix[[factor_name]]
        factor_type <- factor_types[[factor_name]]
        slider_id <- paste0("slider_", factor_name)

        if (factor_type == "Continuous") {
          data_list[[factor_name]] <- input[[slider_id]]
        } else {
          data_list[[factor_name]] <- factor(input[[slider_id]], levels = levels(column_data))
        }
      }

      as.data.frame(data_list, stringsAsFactors = FALSE)
    })

    # Gets a common y-range for all factors
    shared_yrange <- reactive({
      model_formula = model_formula_reactive()
      design_matrix = design_matrix_reactive()
      model_aov = model_aov_reactive()

      req(model_formula, design_matrix, model_aov, new_data())

      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]

      all_preds <- c()

      for (factor_name in main_effects) {
        factor_type <- factor_types[[factor_name]]
        column_data <- design_matrix[[factor_name]]
        current_data <- new_data()

        if (factor_type == "Continuous") {
          x_seq <- seq(min(column_data), max(column_data), length.out = 100)
        } else if (factor_type == "Discrete") {
          x_seq <- sort(unique(column_data))
        } else if (factor_type == "Categorical") {
          x_seq <- levels(column_data)
        }

        tmp_df <- do.call(rbind, lapply(x_seq, function(val) {
          row <- current_data
          row[[factor_name]] <- if (factor_type == "Categorical") {
            factor(val, levels = levels(column_data))
          } else {
            val
          }
          row
        }))

        preds <- predict(model_aov, newdata = tmp_df)
        all_preds <- c(all_preds, preds)
      }

      range(all_preds, na.rm = TRUE)
    })

    # Renders each plot
    observe({
      req(model_formula, design_matrix, model_aov, new_data())

      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]

      for (factor_name in main_effects) {
        local({
          fname <- factor_name
          output_id <- paste0("factor_plot_", fname)

          output[[output_id]] <- renderPlot({
            req(new_data(), shared_yrange())

            column_data <- design_matrix[[fname]]
            factor_type <- factor_types[[fname]]
            current_data <- new_data()

            if (factor_type == "Continuous") {
              x_seq <- seq(min(column_data), max(column_data), length.out = 100)
              tmp_df <- do.call(rbind, lapply(x_seq, function(x) {
                row <- current_data
                row[[fname]] <- x
                row
              }))
              preds <- predict(model_aov, newdata = tmp_df)

              plot(x_seq, preds,
                type = "l", lwd = 2,
                xlab = fname, ylab = selected_response,
                main = paste(fname),
                ylim = shared_yrange()
              )
              grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)

              points(current_data[[fname]],
                predict(model_aov, newdata = current_data),
                col = "#428bca", pch = 19, cex = 1.5
              )
            } else if (factor_type %in% c("Discrete", "Categorical")) {
              x_vals <- if (factor_type == "Categorical") {
                levels(column_data)
              } else {
                sort(unique(column_data))
              }

              tmp_df <- do.call(rbind, lapply(x_vals, function(val) {
                row <- current_data
                row[[fname]] <- if (factor_type == "Categorical") {
                  factor(val, levels = levels(column_data))
                } else {
                  val
                }
                row
              }))
              preds <- predict(model_aov, newdata = tmp_df)

              plot(x_vals, preds,
                type = "b", pch = 19,
                xlab = fname, ylab = selected_response,
                main = paste("Profile for", fname),
                ylim = shared_yrange()
              )
              grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted", lwd = 0.5)

              current_val <- as.character(current_data[[fname]])
              idx <- which(x_vals == current_val)
              points(idx, preds[idx], col = "#428bca", pch = 19, cex = 1.5)
            }
          })
        })
      }
    })

    # Displays the predicted response for the current values
    output$selected_value_output <- renderUI({
      req(new_data())
      prediction <- signif(predict(model_aov, newdata = new_data()), 3)
      p("The predicted ", strong(selected_response), " for the selected values is:  ", strong(prediction))
    })


    # ===== Interact with the Model ============================================

    output$model_interaction_output_ui_id <- renderUI({
      design_matrix = design_matrix_reactive()
      model_formula = model_formula_reactive()
      req(model_formula, design_matrix)

      # Extract factors in the model (e.g. exclude interactions, just keep their components)
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]

      # Set a starting-value input for each factor (default to the midpoint)
      inputs <- lapply(main_effects, function(fname) {
        ftype <- factor_types[[fname]]
        values <- design_matrix[[fname]]

        if (ftype == "Continuous") {
          numericInput(
            inputId = paste0("start_", fname),
            label = paste("Start", fname),
            value = mean(range(values, na.rm = TRUE))
          )
        } else if (ftype == "Discrete") {
          selectInput(
            inputId = paste0("start_", fname),
            label = paste("Start", fname),
            choices = sort(unique(values)),
            selected = unique(values)[1]
          )
        } else {
          selectInput(
            inputId = paste0("start_", fname),
            label = paste("Start", fname),
            choices = levels(values),
            selected = levels(values)[1]
          )
        }
      })

      # Render the UI elements
      tagList(
        h4("Select Optimization Objective"),
        radioButtons("opt_type", "Goal:",
          choices = c("Maximize", "Minimize", "Find Target"),
          selected = "Maximize",
          inline = TRUE
        ),
        conditionalPanel(
          condition = "input.opt_type == 'Find Target'",
          numericInput("target_response", "Target Response", value = 0)
        ),
        hr(),
        h4("Starting Values:"),
        inputs,
        actionButton("run_optim_btn", "Run", class = "btn-primary"),
        hr(),
        h4("Resulting Conditions:"),
        uiOutput("target_solution"),
        h4(paste0("Predicted ", selected_response, ":")),
        uiOutput("target_predicted_response")
      )
    })

    # Run the optimisation as indicated by the user
    # Note:  this may take some time if there are many non-continuous variables
    observeEvent(input$run_optim_btn, {
      req(model_formula, design_matrix, factor_types, model_aov, input$opt_type)

      # Extract factors in the model (e.g. exclude interactions, just keep their components)
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]

      # Partition variables by type
      continuous_vars <- main_effects[factor_types[main_effects] == "Continuous"]
      discrete_vars <- main_effects[factor_types[main_effects] == "Discrete"]
      categorical_vars <- main_effects[factor_types[main_effects] == "Categorical"]
      noncont_vars <- c(discrete_vars, categorical_vars)

      # Build choice sets, i.e. the list of allowed values for each non-continuous factor
      choices_list <- list()
      for (v in noncont_vars) {
        colv <- design_matrix[[v]]
        if (v %in% discrete_vars) {
          choices_list[[v]] <- sort(unique(colv))
        } else {
          choices_list[[v]] <- levels(colv)
        }
      }

      # Create a grid of all combinations for non-continuous factors
      if (length(choices_list) > 0) {
        combo_grid <- do.call(expand.grid, c(choices_list, stringsAsFactors = FALSE))
      } else {
        combo_grid <- data.frame(dummy = 1, stringsAsFactors = FALSE) # single row stub
      }

      # Safety guard against explosion
      max_combos <- 5000L
      if (nrow(combo_grid) > max_combos) {
        showNotification(
          paste(
            "Evaluating", nrow(combo_grid),
            "combinations is large; sampling", max_combos, "rows for speed."
          ),
          type = "warning", duration = 8
        )
        set.seed(1)
        combo_grid <- combo_grid[sample(nrow(combo_grid), max_combos), , drop = FALSE]
      }

      # Set defaults for continuous vars
      cont_start <- numeric(0)
      cont_lower <- numeric(0)
      cont_upper <- numeric(0)

      # Sets up all starting values and bounds for continuous variables
      for (v in continuous_vars) {
        colv <- design_matrix[[v]]
        # Try user-provided starting value; fallback to mid-range
        start_val <- input[[paste0("start_", v)]]
        if (is.null(start_val)) {
          r <- range(colv, na.rm = TRUE)
          start_val <- mean(r)
        }
        cont_start <- c(cont_start, start_val)
        cont_lower <- c(cont_lower, min(colv, na.rm = TRUE))
        cont_upper <- c(cont_upper, max(colv, na.rm = TRUE))
      }
      names(cont_start) <- continuous_vars

      # adjusts the prediction so that the objective is to always minimize
      objective_value <- function(pred, opt_type, target = NULL) {
        if (opt_type == "Maximize") {
          return(-pred)
        }
        if (opt_type == "Minimize") {
          return(pred)
        }
        if (opt_type == "Find Target") {
          validate(need(!is.null(target), "Target response is required."))
          return((pred - target)^2)
        }
        return(pred)
      }

      # Initialize variables to track the global best solution across the grid of non-cont combos
      best_obj <- Inf
      best_df <- NULL
      best_pred <- NA_real_

      # Iterate over each non-continuous combination
      for (i in seq_len(nrow(combo_grid))) {
        # Fix the non-continuous values based on the current combo
        fixed_vars <- list()
        if (length(noncont_vars) > 0) {
          row_i <- combo_grid[i, , drop = FALSE]
          for (v in noncont_vars) {
            if (v %in% categorical_vars) {
              fixed_vars[[v]] <- factor(row_i[[v]], levels = levels(design_matrix[[v]]))
            } else { # Discrete
              fixed_vars[[v]] <- row_i[[v]]
            }
          }
        }

        # This function is called by optim() to assess the 'goodness' of a set of conditions
        fn <- function(par) {
          # par is a numeric vector of current guesses for the continuous variables
          # Assemble full row: continuous (from par) + fixed (from grid)
          test_row <- as.list(par)
          names(test_row) <- continuous_vars
          full_row <- c(test_row, fixed_vars)

          # If there are no continuous vars, ensure we still build a row
          if (length(continuous_vars) == 0) {
            full_row <- fixed_vars
          }

          df <- as.data.frame(full_row, stringsAsFactors = FALSE)

          # Ensures non-continuous factors have the expected leves as factors
          for (n in names(df)) {
            if (is.factor(design_matrix[[n]])) {
              df[[n]] <- factor(df[[n]], levels = levels(design_matrix[[n]]))
            }
          }

          # Calculate the prediction based on the model and values
          pred <- as.numeric(predict(model_aov, newdata = df))
          # send the prediction to the objective function to adjust for minimization (reframes all optimizations as minimization problems)
          objective_value(pred, input$opt_type, target = input$target_response)
        }

        # If we have continuous vars, run bounded optimization; otherwise just evaluate once
        if (length(continuous_vars) > 0) {
          res <- optim(
            par     = cont_start,
            fn      = fn,
            method  = "L-BFGS-B",
            lower   = cont_lower,
            upper   = cont_upper,
            control = list(maxit = 200)
          )

          # Take the best continuous variable values optim found and combine them
          # with the fixed discrete/categorical values for this run to create one complete set of inputs
          # Note:  this is the set of conditions that provided the best response
          # for the given combo of non-cont factors
          sol_cont <- as.list(res$par)
          names(sol_cont) <- continuous_vars
          full_row <- c(sol_cont, fixed_vars)
        } else {
          # No continuous vars: just use the fixed combination
          full_row <- fixed_vars
        }

        # Convert this 'full_row' into a dataframe, calculate the prediciton, and get the objective
        final_df <- as.data.frame(full_row, stringsAsFactors = FALSE)
        for (n in names(final_df)) {
          if (is.factor(design_matrix[[n]])) {
            final_df[[n]] <- factor(final_df[[n]], levels = levels(design_matrix[[n]]))
          }
        }
        pred_val <- as.numeric(predict(model_aov, newdata = final_df))
        obj_val <- objective_value(pred_val, input$opt_type, target = input$target_response)

        # If the new set of conditions is better than previous sets, replace the values
        if (obj_val < best_obj) {
          best_obj <- obj_val
          best_df <- final_df
          best_pred <- pred_val
        }
      }

      # ===== Output best solution =============================================
      output$target_solution <- renderUI({
        vals <- as.list(best_df[1, , drop = FALSE])

        tags$div(
          lapply(names(vals), function(nm) {
            tags$p(
              tags$strong(paste0(nm, ": ")),
              as.character(signif(vals[[nm]]), 3)
            )
          })
        )
      })
      output$target_predicted_response <- renderUI({
        req(best_pred)
        tags$p(
          tags$strong("Predicted Response: "),
          signif(best_pred, 3)
        )
      })
    })
  })

  # Cleans up column names after analysis
  update_factor_names <- function(new_columns, design_matrix) {
    original_columns <- colnames(design_matrix)
    # Initialize a named list to track renamed columns
    renamed_columns <- list()

    # Step 1: Rename columns if they match the pattern and have 2 unique values
    for (col in new_columns) {
      for (original_col in original_columns) {
        # Create a pattern to match the original column name followed by a number
        pattern <- paste0("^", original_col, "[0-9]+$")

        # Check if the current column name matches this pattern
        if (grepl(pattern, col)) {
          # Find out how many unique values are in the corresponding column in design_matrix
          unique_values <- unique(design_matrix[[original_col]])
          num_unique_values <- length(unique_values)

          # If there are exactly 2 unique values, rename the column in new_columns
          if (num_unique_values == 2) {
            new_name <- original_col
            new_columns[new_columns == col] <- new_name
            renamed_columns[[col]] <- new_name
          }
        }
      }
    }

    # Step 2: Update interaction terms in new_columns based on renamed columns
    for (col in new_columns) {
      for (old_name in names(renamed_columns)) {
        # Replace interaction terms in the format "old_name:other" with "new_name:other"
        # but ensure only whole-word matches are replaced
        new_name <- renamed_columns[[old_name]]
        updated_col <- gsub(paste0("\\b", old_name, "\\b"), new_name, col)

        # Update the column name if it was changed
        if (updated_col != col) {
          new_columns[new_columns == col] <- updated_col
        }
      }
    }

    # Step 3: Loop through each column name in new_columns and update interaction notation
    for (col in new_columns) {
      if (grepl(":", col)) {
        # Split the column name by ":"
        terms <- unlist(strsplit(col, ":"))

        # Check if all possible combinations of terms match any name in original_columns
        valid_combination <- TRUE
        for (i in seq_along(terms)) {
          combined_term <- terms[i]

          # Start combining terms one by one to check for multi-part names
          for (j in (i + 1):length(terms)) {
            combined_term <- paste0(combined_term, ":", terms[j])

            # Check if this combined term exists in original_columns
            if (combined_term %in% original_columns) {
              terms <- c(terms[1:(i - 1)], combined_term, terms[(j + 1):length(terms)])
              break
            }
          }

          # If any part of the name isn't valid, stop the process
          if (!(terms[i] %in% original_columns)) {
            valid_combination <- FALSE
            break
          }
        }

        # If the combination is valid, replace ":" with "*" in the column name
        if (valid_combination) {
          new_name <- paste(terms, collapse = "*")
          new_columns[new_columns == col] <- new_name
        }
      }
    }

    for (col in new_columns) {
      if (grepl("\\.", col)) {
        # Split the column name by ":"
        terms <- unlist(strsplit(col, "\\."))

        # Check if all possible combinations of terms match any name in original_columns
        valid_combination <- TRUE
        for (i in seq_along(terms)) {
          combined_term <- terms[i]

          # Start combining terms one by one to check for multi-part names
          for (j in (i + 1):length(terms)) {
            combined_term <- paste0(combined_term, ".", terms[j])

            # Check if this combined term exists in original_columns
            if (combined_term %in% original_columns) {
              terms <- c(terms[1:(i - 1)], combined_term, terms[(j + 1):length(terms)])
              break
            }
          }

          # If any part of the name isn't valid, stop the process
          if (!(terms[i] %in% original_columns)) {
            valid_combination <- FALSE
            break
          }
        }

        # If the combination is valid, replace ":" with "*" in the column name
        if (valid_combination) {
          new_name <- paste(terms, collapse = "*")
          new_columns[new_columns == col] <- new_name
        }
      }
    }
    # Return the updated column names
    new_columns
  }
}

# ===== Start the App ==========================================================
shinyApp(ui, server)
