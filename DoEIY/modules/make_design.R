# Make a Design Shiny Module

make_design_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
      fluidRow(
        column(
          4,
          selectizeInput(inputId = ns("type_of_design_select_id"), label = "Select a Design:", 
                         choices = list("Screening Designs" = list("Plackett-Burman"), 
                                        "Factorial" = list("Full Factorial", "Fractional Factorial"), 
                                        "Response Surface Designs" = list("Box-Behnken", "Central Composite"), 
                                        "Space Filling" = list("Latin Hypercube Sampling"), 
                                        "Custom" = list("D-Optimal")), 
                         selected = NULL, width = "100%", 
                         options = list(
                           placeholder = "Please select a suitable design",
                           onInitialize = I('function() { this.setValue(""); }')
                         )),
          br(),
          actionButton(inputId = ns("make_design_id"), label = "Make Design", width = "100%", class = "Blue-Button", style = "height: 34px;")
        ),
        column(
          8,
          uiOutput(outputId = ns("selected_design_text_id"))
        )
      )
    )
  )
}

make_design_server <- function(id, design_reactive, factor_data_reactive, factor_types_reactive, d_optimal_model_components) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Internal reactives for make_design module
    column_widths_reactive <- reactiveVal(value = NULL)
    max_factors_reactive <- reactiveVal(value = NULL)
    min_factors_reactive <- reactiveVal(value = NULL)
    set_num_levels_reactive <- reactiveVal(value = NULL)
    factors_string_reactive <- reactiveVal(value = NULL)
    
    DOD_min_runs_reactive <- reactiveVal(NULL)
    DOD_max_runs_reactive <- reactiveVal(NULL)
    
    # Description texts based on selected design
    observeEvent(input$type_of_design_select_id, {
      if (input$type_of_design_select_id == "") {
        output$selected_design_text_id <- renderUI({ p("") })
      }
      if (input$type_of_design_select_id == "Plackett-Burman") {
        output$selected_design_text_id <- renderUI({
          p(HTML('<p style="text-align: justify;"><strong>Plackett-Burman Designs</strong> are efficient screening tools for estimating main effects in an experiment. These designs allow for independent estimation of main effects, though they may be confounded with higher-order interactions. This software offers 2-level Plackett-Burman designs suitable for experiments involving 4 to 23 factors.</p>'), style = "margin-top: 20px;")
        })
      }
      if (input$type_of_design_select_id == "Full Factorial") {
        output$selected_design_text_id <- renderUI({
          p(HTML('<p style="text-align: justify;"><strong>Full Factorial Designs</strong> investigate every possible combination of factor levels, using two or more levels for each factor. This comprehensive approach enables the independent estimation of all main effects and higher-order interactions, providing a complete understanding of factor relationships.</p>'), style = "margin-top: 20px;")
        })
      }
      if (input$type_of_design_select_id == "Fractional Factorial") {
        output$selected_design_text_id <- renderUI({
          p(HTML("<p style='text-align: justify;'><strong>Fractional Factorial Designs</strong> are efficient alternatives to full factorial designs, utilizing only a fraction of the total possible runs. They allow for the estimation of main effects and, depending on the design's resolution, certain interactions. This software provides 2-level fractional factorial designs with the following resolutions: </p>
          <ul>
          <li style='text-align: justify;'><strong>Resolution III</strong> - main effects are not confounded with each other but may be confounded with two-factor interactions.</li>
          <li style='text-align: justify;'><strong>Resolution IV</strong> - main effects are not confounded with other main effects or two-factor interactions; two-factor interactions may still be confounded with each other.</li>
          <li style='text-align: justify;'><strong>Resolution V</strong> - two-factor interactions are not confounded with main effects or other two-factor interactions, but may be confounded with three-factor interactions.</li>
          </ul>"), style = "margin-top: 20px;")
        })
      }
      if (input$type_of_design_select_id == "Box-Behnken") {
        output$selected_design_text_id <- renderUI({
          p(HTML("<p><strong>Box-Behnken Designs</strong> are used to estimate quadratic models by employing three levels for each factor. This design enables the independent estimation of main effects, two-factor interactions, and quadratic effects. This software offers Box-Behnken designs for experiments with 3 to 7 factors.</p>"), style = "margin-top: 20px;")
        })
      }
      if (input$type_of_design_select_id == "Central Composite") {
        output$selected_design_text_id <- renderUI({
          p(HTML("<p style='text-align: justify;'><strong>Central Composite Designs (CCD)</strong> are response surface designs used to estimate quadratic models. They consist of a factorial or fractional factorial design, additional 'star' (axial) points, and center points. This structure allows for the estimation of curvature and the fitting of second-order response surfaces. This software offers Circumscribed, Inscribed, and Face Centered central composite designs for experiments with 2 or more factors.</p>"), style = "margin-top: 20px;")
        })
      }
      if (input$type_of_design_select_id == "Latin Hypercube Sampling") {
        output$selected_design_text_id <- renderUI({
          p(HTML("<p style='text-align: justify;'><strong>Latin Hypercube Sampling (LHS)</strong> is a space-filling design used to explore multi-dimensional parameter spaces. Unlike traditional factorial designs, LHS ensures that each factor is uniformly sampled across its range by dividing the range of each factor into equal intervals and selecting one point from each interval. This software provides improved LHS designs, optimizing the space-filling properties, suitable for computer experiments and simulations with any number of factors.</p>"), style = "margin-top: 20px;")
        })
      }
      if (input$type_of_design_select_id == "D-Optimal") {
        output$selected_design_text_id <- renderUI({
          p(HTML("<p style='text-align: justify;'><strong>D-Optimal Designs</strong> are computer-generated designs optimized based on a specified model. They choose a subset of runs from a candidate set of all possible combinations to minimize the variance of the estimated model coefficients (specifically maximizing the determinant of the information matrix). This approach is highly flexible, allowing for custom model specifications, constraints on run sizes, and the inclusion of continuous, discrete, and categorical factors with unequal levels.</p>"), style = "margin-top: 20px;")
        })
      }
    })
    
    # Opens modal to gather factor details
    observeEvent(input$make_design_id, {
      selected_design <- input$type_of_design_select_id
      
      if (selected_design == "") {
        showNotification(strong("A type of design must be selected."), duration = 15, closeButton = TRUE, type = "error")
      } else {
        if (selected_design == "Plackett-Burman") {
          factors_string <- paste0(strong("Plackett-Burman Designs"), " are supported for 4-23 factors. These factors can be continuous or categorical. For each factor, only high and low levels are considered. Please provide the levels as comma separated values, e.g., '-1, 1'.")
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
          factors_string <- paste0(strong("Box-Behnken Designs"), " are supported for 3-12 factors. These factors must be continuous. For each factor, three levels are considered: a high, mid, and low level. Please provide the levels as comma separated values, e.g., '-1, 0, 1'. If only two values are provided for the levels, the mid level will be automatically calculated as the midpoint of the provided values.")
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
          factors_string <- paste0(strong("Full Factorial Designs"), " are supported for 1-25 factors, which can be continuous, discrete, categorical or blocking. For each factor, 2-10 levels can be considered, select the appropriate number of levels for the factor from the corresponding dropdown. Please provide the values for the levels as comma separated values, e.g., '5, 10, 15'.")
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
          factors_string <- paste0(strong("Fractional Factorial"), " are supported for 3 to 25 factors, which can be continuous, categorical, or discrete. Each factor is limited to two levels. Please provide the values for the levels as comma separated values, e.g., '-1, 1'. When setting up the design, select the required resolution from the dropdown menu. Note: Resolution III requires at least 3 factors, Resolution IV requires at least 4 factors, and Resolution V requires at least 5 factors.")
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
          factors_string <- paste0(strong("Central Composite"), " are supported for 2 or more factors. These factors must be continuous. For each factor, three levels are considered. Please provide the upper and lower limits for these levels as comma separated values, e.g., '-1, 1'. The mid level will be automatically calculated as the midpoint of the provided values.")
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
          factors_string <- paste0(strong("Latin Hypercube Designs"), " are supported for any number of factors or runs. These factors can be continuous, discrete, categorical, or blocks. For continuous factors, provide the upper and lower limits as comma separated values, e.g., -1, 1. For discrete, categorical, or blocking factors, provide the allowable values as comma separated values. Note: if blocking factors are included, it is recommended to choose a number of runs that is divisible by the number of blocks.")
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
          factors_string <- paste0(strong("D-Optimal Designs"), " are supported for any number of factors, which can be continuous, discrete, categorical, or blocks. For continuous factors, provide the upper and lower limits as comma separated values, e.g., -1, 1. For discrete, categorical, or blocking factors, provide the allowable values as comma separated values. Note: if blocking factors are included, it is recommended to choose a number of runs that is divisible by the number of blocks.")
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
            column(12, p(HTML(factors_string))),
            column(12, p(strong("All entries in the table below must be completed. Factor Names can only contain letters, numbers, and underscores. Do NOT use spaces or special characters in the factor names."), style = "margin-bottom: 25px;"))
          ),
          p("Note: Add or remove rows by right clicking on a row."),
          div(style = "height: 220px;", rHandsontableOutput(ns("Factor_details_table_id"), height = "100%")),
          uiOutput(ns("warning_message_id")),
          easyClose = FALSE,
          footer = tagList(
            fluidRow(
              column(2, offset = 8, actionButton(inputId = ns("open_blocking_modal_id"), "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
              column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
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
    
    # Close modal triggers
    observeEvent(input$close_modal_id, { removeModal() })
    
    # Enforces min/max factors
    observe({
      req(input$Factor_details_table_id)
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
    
    # Confirm factor configuration details and navigate to blocking modal
    observeEvent(input$open_blocking_modal_id, {
      req(input$Factor_details_table_id)
      current_factor_data <- hot_to_r(input$Factor_details_table_id)
      factor_data_reactive(current_factor_data)
      
      warning_msg <- c()
      
      for (i in seq_len(dim(current_factor_data)[1])) {
        warning_msg_temp <- c()
        factor_name <- current_factor_data$`Factor Names`[i]
        
        if (is.na(factor_name) || grepl("^\\s*$", factor_name) || factor_name == "") {
          warning_msg_temp <- c(warning_msg_temp, paste0("factor name is missing"))
        }
        if (is.na(factor_name) || grepl("^\\s*$", factor_name) || factor_name == "" || factor_name != make.names(factor_name)) {
          warning_msg_temp <- c(warning_msg_temp, "factor name cannot contain spaces or special characters")
        }
        
        if ("Factor Type" %in% colnames(current_factor_data)) {
          factor_type <- current_factor_data$`Factor Type`[i]
          if (is.na(as.character(factor_type))) {
            warning_msg_temp <- c(warning_msg_temp, paste0("factor type is missing"))
          }
        } else {
          factor_type <- "Continuous"
        }
        
        if ("Num Levels" %in% colnames(current_factor_data)) {
          num_levels <- current_factor_data$`Num Levels`[i]
          if (is.na(as.character(num_levels))) {
            warning_msg_temp <- c(warning_msg_temp, paste0("number of levels is missing"))
          }
        } else {
          num_levels <- set_num_levels_reactive()
        }
        
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
        
        if (length(warning_msg_temp) != 0) {
          warning_msg_temp <- paste0("Row ", i, ": ", paste0(warning_msg_temp, collapse = "; "), ".")
          warning_msg <- c(warning_msg, warning_msg_temp)
        }
      }
      
      # If duplicate factor names exist, trigger a warning
      if (length(unique(current_factor_data$`Factor Names`)) != length(current_factor_data$`Factor Names`)) {
        warning_msg <- c(warning_msg, "Factor Names must be unique. Duplicate Factor Names are not allowed.")
      }
      
      if (length(warning_msg) != 0) {
        output$warning_message_id <- renderUI({
          lapply(seq_len(length(warning_msg)), function(j) {
            p(strong(warning_msg[j]), class = "text-red")
          })
        })
      }
      req(length(warning_msg) == 0, cancelOutput = TRUE)
      
      removeModal()
      
      selected_design <- input$type_of_design_select_id
      n_factors <- dim(current_factor_data)[1]
      
      # Determine design specific selections inside the Details modal
      showModal(modalDialog(
        title = h3(strong("Design Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(6,
                 p(HTML(paste0(strong("Selected Design: "), selected_design))),
                 DT::dataTableOutput(ns("summary_factors_table_id"))
          ),
          column(6, uiOutput(ns("design_specific_selections_id")))
        ),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2, offset = 6, actionButton(inputId = ns("reopen_factors_modal_id"), "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("create_design_id"), "Create", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
          )
        )
      ))
      
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
      
      if (selected_design == "Plackett-Burman") {
        design <- Plackett_Burman_Designs(n_factors)
        num_runs <- dim(design)[1]
        output$design_specific_selections_id <- renderUI({
          fluidRow(
            p(HTML(paste0(strong("Number of Runs: "), num_runs))),
            p("Blocking is not available for this design"),
            materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        })
      } else if (selected_design == "Full Factorial") {
        factor_levels <- as.numeric(as.character(current_factor_data$`Num Levels`))
        design <- Full_Factorial_Designs(factor_levels)
        num_runs <- dim(design)[1]
        output$design_specific_selections_id <- renderUI({
          fluidRow(
            p(HTML(paste0(strong("Number of Runs: "), num_runs))),
            p("Blocking is not available for this design"),
            materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
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
        
        output$design_specific_selections_id <- renderUI({
          fluidRow(
            selectInput(inputId = ns("fractional_factorial_runs_select_id"), label = "Number of Experimental Runs:", choices = possible_runs, selected = num_runs, width = "100%"),
            uiOutput(ns("fractional_factorial_blocks_UI_id")),
            uiOutput(ns("fractional_factorial_resolution_ui_id")),
            materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        })
      } else if (selected_design == "Box-Behnken") {
        design <- Box_Behnken_Designs(n_factors)
        num_runs <- dim(design)[1]
        
        output$design_specific_selections_id <- renderUI({
          fluidRow(
            p(HTML(paste0(strong("Number of Runs: "), num_runs))),
            selectInput(inputId = ns("fractional_factorial_blocks_select_id"), label = "Number of Blocks:", choices = c(1, unique(design$Block)), selected = NULL, width = "100%"),
            materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        })
      } else if (selected_design == "Central Composite") {
        output$design_specific_selections_id <- renderUI({
          fluidRow(
            selectInput(inputId = ns("central_composite_type_id"), label = "Select Central Composite Type:", choices = c("Circumscribed", "Inscribed", "Face Centered"), selected = "Face Centered", width = "100%"),
            p("Blocking is not available for this design"),
            materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        })
      } else if (selected_design == "Latin Hypercube Sampling") {
        output$design_specific_selections_id <- renderUI({
          fluidRow(
            numericInput(inputId = ns("LHS_runs_id"), label = "Number of Experimental Runs:", value = 10, min = 1, width = "100%"),
            materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
          )
        })
      } else if (selected_design == "D-Optimal") {
        factor_names <- current_factor_data$`Factor Names`
        factor_types <- as.character(current_factor_data$`Factor Type`)
        
        # Populate DOD model components selection
        d_optimal_model_components(factor_names)
        
        # Open model building dialog first
        removeModal()
        
        showModal(modalDialog(
          title = h3(strong("Build Model"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
          fluidRow(
            column(12, p("Provide the components of the model that you wish to fit to the experimental results.  Select additional interaction or quadratic terms from the dropdown menus below.")),
            column(6,
                   h4(strong("Allowed Model Terms:")),
                   selectInput(inputId = ns("selected_factors_DOD_id"), label = NULL, choices = d_optimal_model_components(), selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%", size = 10)
            ),
            column(6,
                   br(),
                   fluidRow(
                     column(10, offset = 1, actionButton(inputId = ns("DOD_add_interactions_id"), "Add Interactions", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))
                   ),
                   br(),
                   fluidRow(
                     column(10, offset = 1, actionButton(inputId = ns("DOD_add_quadratics_id"), "Add Quadratics", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))
                   ),
                   br(),
                   fluidRow(
                     column(10, offset = 1, actionButton(inputId = ns("DOD_remove_interactions_id"), "Remove Selected Terms", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))
                   )
            )
          ),
          easyClose = FALSE,
          footer = tagList(
            fluidRow(
              column(2, offset = 6, actionButton(inputId = ns("reopen_factors_modal_id"), "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
              column(2, actionButton(inputId = ns("open_dod_modal_id"), "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
              column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
            )
          )
        ))
      }
    })
    
    # Back action to reopen the factor spec table modal
    observeEvent(input$reopen_factors_modal_id, {
      removeModal()
      showModal(modalDialog(
        title = h3(strong("Provide Factor Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(12, p(HTML(factors_string_reactive()))),
          column(12, p(strong("All entries in the table below must be completed. Factor Names can only contain letters, numbers, and underscores. Do NOT use spaces or special characters in the factor names."), style = "margin-bottom: 25px;"))
        ),
        p("Note: Add or remove rows by right clicking on a row."),
        div(style = "height: 220px;", rHandsontableOutput(ns("Factor_details_table_id"), height = "100%")),
        uiOutput(ns("warning_message_id")),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2, offset = 8, actionButton(inputId = ns("open_blocking_modal_id"), "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
          )
        )
      ))
      output$Factor_details_table_id <- renderRHandsontable({
        rhandsontable(factor_data_reactive(), selectCallback = TRUE, width = 1150, height = 250) %>%
          hot_cols(colWidths = column_widths_reactive()) %>%
          hot_context_menu(allowRowEdit = TRUE, allowColEdit = FALSE)
      })
    })
    
    # Dynamic blocks for Fractional Factorial
    observeEvent(input$fractional_factorial_runs_select_id, {
      n_factors <- dim(factor_data_reactive())[1]
      num_runs <- as.numeric(input$fractional_factorial_runs_select_id)
      
      max_block_power <- log2(num_runs / 2)
      possible_block_runs <- 2^(0:max_block_power)
      
      output$fractional_factorial_blocks_UI_id <- renderUI({
        selectInput(inputId = ns("fractional_factorial_blocks_select_id"), label = "Number of Blocks:", choices = possible_block_runs, selected = possible_block_runs[1], width = "100%")
      })
    })
    
    observe({
      req(input$fractional_factorial_runs_select_id, input$fractional_factorial_blocks_select_id)
      n_factors <- dim(factor_data_reactive())[1]
      num_runs <- as.numeric(input$fractional_factorial_runs_select_id)
      num_blocks <- as.numeric(input$fractional_factorial_blocks_select_id)
      
      design <- Fractional_Factorial_Designs(num_runs, n_factors, num_blocks)
      output$fractional_factorial_resolution_ui_id <- renderUI({
        p(HTML(paste0("This is a ", strong(design$Resolution), " design.")))
      })
    })
    
    # D-Optimal Model specification observers
    observeEvent(input$DOD_add_interactions_id, {
      current_factor_data <- factor_data_reactive()
      factor_names <- as.character(current_factor_data$`Factor Names`)
      req(length(factor_names) > 1, cancelOutput = TRUE)
      
      showModal(modalDialog(
        title = h3(strong("Add Interaction Term")),
        selectInput(inputId = ns("DOD_interaction_factor1_id"), label = "Factor 1", choices = factor_names),
        selectInput(inputId = ns("DOD_interaction_factor2_id"), label = "Factor 2", choices = factor_names),
        easyClose = FALSE,
        footer = tagList(
          actionButton(inputId = ns("DOD_interactions_add_btn_id"), "Add", class = "Blue-Button", style = "height:34px; margin:5px;"),
          actionButton(inputId = ns("dod_model_select_back_btn_id"), "Back", class = "Grey-Button", style = "height:34px; margin:5px;")
        )
      ))
    })
    
    observeEvent(input$DOD_interactions_add_btn_id, {
      req(input$DOD_interaction_factor1_id, input$DOD_interaction_factor2_id)
      f1 <- input$DOD_interaction_factor1_id
      f2 <- input$DOD_interaction_factor2_id
      
      if (f1 == f2) {
        showNotification("Interaction terms must be between different factors.", type = "warning")
      } else {
        new_interaction <- paste(sort(c(f1, f2)), collapse = ":")
        currently_selected_factors <- d_optimal_model_components()
        updated_selected_factors <- unique(c(currently_selected_factors, new_interaction))
        term_types <- sapply(updated_selected_factors, classify_terms)
        terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
        sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
        sorted_terms <- sorted_terms_df$term
        d_optimal_model_components(sorted_terms)
      }
      
      # Reopen Build Model Modal
      removeModal()
      showBuildModelModal()
    })
    
    observeEvent(input$DOD_add_quadratics_id, {
      current_factor_data <- factor_data_reactive()
      continuous_discrete_factors <- as.character(current_factor_data$`Factor Names`[current_factor_data$`Factor Type` %in% c("Continuous", "Discrete")])
      
      if (length(continuous_discrete_factors) == 0) {
        showNotification("Quadratic terms are only supported for Continuous or Discrete factors.", type = "warning")
      } else {
        showModal(modalDialog(
          title = h3(strong("Add Quadratic Term")),
          selectInput(inputId = ns("DOD_quadratic_factor_id"), label = "Factor", choices = continuous_discrete_factors),
          easyClose = FALSE,
          footer = tagList(
            actionButton(inputId = ns("DOD_quadratics_add_btn_id"), "Add", class = "Blue-Button", style = "height:34px; margin:5px;"),
            actionButton(inputId = ns("dod_model_select_back_btn_id"), "Back", class = "Grey-Button", style = "height:34px; margin:5px;")
          )
        ))
      }
    })
    
    observeEvent(input$DOD_quadratics_add_btn_id, {
      req(input$DOD_quadratic_factor_id)
      f <- input$DOD_quadratic_factor_id
      new_quadratic <- paste0(f, ":", f)
      
      currently_selected_factors <- d_optimal_model_components()
      updated_selected_factors <- unique(c(currently_selected_factors, new_quadratic))
      term_types <- sapply(updated_selected_factors, classify_terms)
      terms_df <- data.frame(term = updated_selected_factors, type = term_types, stringsAsFactors = FALSE)
      sorted_terms_df <- terms_df[order(terms_df$type, terms_df$term), ]
      sorted_terms <- sorted_terms_df$term
      d_optimal_model_components(sorted_terms)
      
      removeModal()
      showBuildModelModal()
    })
    
    observeEvent(input$DOD_remove_interactions_id, {
      req(input$selected_factors_DOD_id)
      terms_to_remove <- input$selected_factors_DOD_id
      currently_selected_factors <- d_optimal_model_components()
      updated_selected_factors <- setdiff(currently_selected_factors, terms_to_remove)
      
      d_optimal_model_components(updated_selected_factors)
      updateSelectInput(session, inputId = "selected_factors_DOD_id", choices = updated_selected_factors)
    })
    
    observeEvent(input$dod_model_select_back_btn_id, {
      removeModal()
      showBuildModelModal()
    })
    
    showBuildModelModal <- function() {
      showModal(modalDialog(
        title = h3(strong("Build Model"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(12, p("Provide the components of the model that you wish to fit to the experimental results. Select additional interaction or quadratic terms from the dropdown menus below.")),
          column(6,
                 h4(strong("Allowed Model Terms:")),
                 selectInput(inputId = ns("selected_factors_DOD_id"), label = NULL, choices = d_optimal_model_components(), selected = NULL, multiple = TRUE, selectize = FALSE, width = "100%", size = 10)
          ),
          column(6,
                 br(),
                 fluidRow(column(10, offset = 1, actionButton(inputId = ns("DOD_add_interactions_id"), "Add Interactions", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))),
                 br(),
                 fluidRow(column(10, offset = 1, actionButton(inputId = ns("DOD_add_quadratics_id"), "Add Quadratics", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;"))),
                 br(),
                 fluidRow(column(10, offset = 1, actionButton(inputId = ns("DOD_remove_interactions_id"), "Remove Selected Terms", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")))
          )
        ),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2, offset = 6, actionButton(inputId = ns("reopen_factors_modal_id"), "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("open_dod_modal_id"), "Next", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
          )
        )
      ))
    }
    
    # DOD parameter calculations & dynamic Runs Spec Modal
    observeEvent(input$open_dod_modal_id, {
      factor_data <- factor_data_reactive()
      components <- d_optimal_model_components()
      
      factor_info <- lapply(seq_len(nrow(factor_data)), function(i) {
        fname <- factor_data$`Factor Names`[i]
        ftype <- factor_data$`Factor Type`[i]
        levels <- strsplit(factor_data$`Level Values`[i], ",")[[1]]
        list(name = fname, type = ftype, levels = levels)
      })
      names(factor_info) <- factor_data$`Factor Names`
      
      num_parameters <- sum(sapply(components, function(t) term_params(t, factor_info))) + 1
      
      max_runs <- prod(sapply(names(factor_info), function(f) {
        info <- factor_info[[f]]
        if (info$type == "Continuous") {
          deg <- factor_degree(f, components)
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
      
      removeModal()
      
      showModal(modalDialog(
        title = h3(strong("Design Details"), style = "margin: 0px; padding: 0px; color: #0097a9; margin-bottom: 5px;"),
        fluidRow(
          column(6,
                 p(HTML(paste0(strong("Selected Design: "), "D-Optimal"))),
                 DT::dataTableOutput(ns("summary_factors_table_id"))
          ),
          column(6,
                 fluidRow(
                   uiOutput(ns("DOD_min_runs_ui_id")),
                   uiOutput(ns("DOD_max_runs_ui_id")),
                   numericInput(inputId = ns("DOD_num_runs_id"), label = "Number of Experimental Runs:", value = num_parameters, min = num_parameters, max = max_runs, width = "100%"),
                   materialSwitch(inputId = ns("randomize_checkbox_id"), label = "Randomize Runs:", value = TRUE, width = "100%", status = "primary")
                 )
          )
        ),
        easyClose = FALSE,
        footer = tagList(
          fluidRow(
            column(2, offset = 6, actionButton(inputId = ns("open_dod_modal_back_btn_id"), "Back", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("create_design_id"), "Create", width = "100%", class = "Blue-Button", style = "height: 34px; margin: 5px;")),
            column(2, actionButton(inputId = ns("close_modal_id"), strong("Cancel"), width = "100%", class = "Grey-Button", style = "height: 34px; margin: 5px;"))
          )
        )
      ))
      
      output$DOD_min_runs_ui_id <- renderUI(p(strong(paste0("Min Runs: ", num_parameters))))
      output$DOD_max_runs_ui_id <- renderUI(p(strong(paste0("Max Runs: ", max_runs))))
      
      output$summary_factors_table_id <- DT::renderDataTable(
        factor_data,
        escape = FALSE, selection = "none", server = FALSE,
        extensions = c("FixedColumns", "FixedHeader", "Scroller"),
        options = list(
          dom = "t",
          pageLength = 10, paging = TRUE,
          headerCallback = DT::JS("function(thead) {", "  $(thead).css('font-size', '1em');", "}")
        ),
        rownames = FALSE, filter = "none"
      )
    })
    
    observeEvent(input$open_dod_modal_back_btn_id, {
      removeModal()
      showBuildModelModal()
    })
    
    # Executes the design creation upon clicking "Create"
    observeEvent(input$create_design_id, {
      factor_data <- factor_data_reactive()
      selected_design <- input$type_of_design_select_id
      n_factors <- dim(factor_data)[1]
      factor_names <- factor_data$`Factor Names`
      
      randomize <- input$randomize_checkbox_id
      
      # Generates design based on design type
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
        num_blocks <- as.numeric(input$fractional_factorial_blocks_select_id)
        design <- Box_Behnken_Designs(n_factors)
        
        if (num_blocks == 1) {
          design <- design[, colnames(design) != "Block", drop = FALSE]
          colnames(design) <- factor_names
        } else {
          block_col <- design$Block
          design <- design[, colnames(design) != "Block", drop = FALSE]
          design <- cbind(design, Block = block_col)
          colnames(design) <- c(factor_names, "Block")
        }
      } else if (selected_design == "Central Composite") {
        ccd_type <- input$central_composite_type_id
        design <- Central_Composite_Designs(n_factors, ccd_type)
        colnames(design) <- factor_names
      } else if (selected_design == "Latin Hypercube Sampling") {
        val_res <- validate_num_runs(input$LHS_runs_id, DOD_min_runs_reactive(), DOD_max_runs_reactive())
        if (!val_res$valid) {
          showNotification(val_res$message, duration = 20, closeButton = TRUE)
        }
        req(val_res$valid, cancelOutput = TRUE)
        num_runs <- as.integer(round(input$LHS_runs_id))
        
        design <- Latin_Hypercube_Designs(n_factors, num_runs)
        colnames(design) <- factor_names
      } else if (selected_design == "D-Optimal") {
        model_components <- d_optimal_model_components()
        
        val_res <- validate_num_runs(input$DOD_num_runs_id, DOD_min_runs_reactive(), DOD_max_runs_reactive())
        if (!val_res$valid) {
          showNotification(val_res$message, duration = 20, closeButton = TRUE)
        }
        req(val_res$valid, cancelOutput = TRUE)
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
            continuous_values <- min_val + design[, i] * (max_val - min_val)
            snapped_values <- sapply(continuous_values, function(val) {
              allowable_values[which.min(abs(allowable_values - val))]
            })
            design[, i] <- snapped_values
          } else if (factor_type == "Categorical") {
            allowable_values <- trimws(unlist(strsplit(bound_str, ",")))
            num_categories <- length(allowable_values)
            intervals <- seq(0, 1, length.out = num_categories + 1)
            coded_values <- cut(design[, i], breaks = intervals, labels = allowable_values, include.lowest = TRUE)
            design[, i] <- as.character(coded_values)
          }
        }
      } else if (selected_design %in% c("Plackett-Burman", "Full Factorial", "Fractional Factorial", "Box-Behnken", "D-Optimal")) {
        bounds <- factor_data$`Level Values`
        
        if ("Factor Type" %in% colnames(factor_data)) {
          factor_types <- as.character(factor_data$`Factor Type`)
        } else {
          factor_types <- rep("Continuous", n_factors)
        }
        
        cols_to_scale <- setdiff(colnames(design), "Block")
        
        for (col_name in cols_to_scale) {
          col_idx <- which(factor_names == col_name)
          factor_type <- factor_types[col_idx]
          bound_str <- bounds[col_idx]
          
          if (factor_type == "Continuous") {
            lower_upper <- as.numeric(unlist(strsplit(bound_str, ",")))
            lower <- min(lower_upper)
            upper <- max(lower_upper)
            design[[col_name]] <- lower + (design[[col_name]] + 1) * (upper - lower) / 2
          } else if (factor_type == "Discrete") {
            allowable_values <- sort(as.numeric(unlist(strsplit(bound_str, ","))))
            lower <- min(allowable_values)
            upper <- max(allowable_values)
            continuous_values <- lower + (design[[col_name]] + 1) * (upper - lower) / 2
            snapped_values <- sapply(continuous_values, function(val) {
              allowable_values[which.min(abs(allowable_values - val))]
            })
            design[[col_name]] <- snapped_values
          } else if (factor_type == "Categorical") {
            allowable_values <- trimws(unlist(strsplit(bound_str, ",")))
            design[[col_name]] <- allowable_values[design[[col_name]]]
          }
        }
      }
      
      # Adds blocks if present
      if ("Block" %in% colnames(design)) {
        block_col <- design$Block
        design <- subset(design, select = -c(Block))
        design$Block <- as.factor(block_col)
      }
      
      # Save the factor types vector
      factor_types_vector <- setNames(as.character(factor_data$`Factor Type`), factor_names)
      if ("Block" %in% colnames(design)) {
        factor_types_vector <- c(factor_types_vector, Block = "Blocking")
      }
      
      # Randomize order if requested
      if (randomize) {
        design <- design[sample(nrow(design)), , drop = FALSE]
      }
      
      rownames(design) <- NULL
      
      # Set shared reactive states
      design_reactive(design)
      factor_data_reactive(factor_data)
      factor_types_reactive(factor_types_vector)
      
      removeModal()
      
      # Switch tabs to the Enter Results tab
      updateTabItems(session, inputId = "tabs", selected = "Enter_edit_results_id")
      showNotification(strong("Design created successfully! Enter responses on the next tab."), duration = 10, type = "message")
    })
  })
}
