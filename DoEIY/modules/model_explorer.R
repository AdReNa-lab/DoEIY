# Model Explorer Shiny Module

model_explorer_ui <- function(id) {
  ns <- NS(id)
  tagList(
    box(
      id = NULL, title = NULL, status = "primary", solidHeader = FALSE, collapsible = FALSE, headerBorder = FALSE, width = 12,
      tabBox(
        id = ns("Explorer_tabs_id"), title = NULL, side = "left", width = "100%",
        tabPanel(
          "Visualizer",
          uiOutput(ns("model_visualizer_output_ui_id"))
        ),
        tabPanel(
          "Response Optimizer",
          uiOutput(ns("model_interaction_output_ui_id"))
        )
      )
    )
  )
}

model_explorer_server <- function(id, design_matrix_reactive, factor_types_reactive, 
                                  model_formula_reactive, model_aov_reactive) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Observe reactive state to render visualizer UI
    output$model_visualizer_output_ui_id <- renderUI({
      model_formula <- model_formula_reactive()
      design_matrix <- design_matrix_reactive()
      model_aov <- model_aov_reactive()
      factor_types <- factor_types_reactive()
      
      if (is.null(model_formula) || is.null(design_matrix) || is.null(model_aov)) {
        return(p("Create and fit a model in ", strong("Analyze the Design"), " to initialize the visualizer."))
      }
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
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
      
      ncol <- ceiling(sqrt(n))
      
      # Build plot-slider UI elements
      plot_slider_grid <- lapply(seq_along(main_effects), function(i) {
        factor_name <- main_effects[i]
        column_data <- design_matrix[[factor_name]]
        slider_id <- ns(paste0("slider_", factor_name))
        factor_type <- types_vec[[factor_name]]
        
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
          plotOutput(ns(paste0("factor_plot_", factor_name)), height = "250px"),
          slider_ui
        )
      })
      
      tagList(
        p("Use the sliders below to adjust the values of each factor. The plots will update automatically, and the Predicted Response will reflect your current selections."),
        h4("Predicted Response:"),
        hr(),
        uiOutput(ns("selected_value_output")),
        fluidRow(plot_slider_grid)
      )
    })
    
    # Reconstruct new_data frame from slider values
    new_data <- reactive({
      model_formula <- model_formula_reactive()
      design_matrix <- design_matrix_reactive()
      factor_types <- factor_types_reactive()
      req(model_formula, design_matrix)
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]
      
      data_list <- list()
      for (factor_name in main_effects) {
        column_data <- design_matrix[[factor_name]]
        factor_type <- types_vec[[factor_name]]
        slider_id <- paste0("slider_", factor_name)
        
        # Check input exists
        req(input[[slider_id]])
        
        if (factor_type == "Continuous") {
          data_list[[factor_name]] <- as.numeric(input[[slider_id]])
        } else if (factor_type == "Discrete") {
          # Keep discrete values as numeric/character accordingly
          val <- input[[slider_id]]
          data_list[[factor_name]] <- if (is.numeric(column_data)) as.numeric(val) else val
        } else {
          data_list[[factor_name]] <- factor(input[[slider_id]], levels = levels(column_data))
        }
      }
      
      as.data.frame(data_list, stringsAsFactors = FALSE)
    })
    
    # Calculate a shared y range across all prediction profiles
    shared_yrange <- reactive({
      model_formula <- model_formula_reactive()
      design_matrix <- design_matrix_reactive()
      model_aov <- model_aov_reactive()
      factor_types <- factor_types_reactive()
      
      req(model_formula, design_matrix, model_aov, new_data())
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]
      
      all_preds <- c()
      
      for (factor_name in main_effects) {
        factor_type <- types_vec[[factor_name]]
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
    
    # Generate visualizer plots dynamically
    observe({
      model_formula <- model_formula_reactive()
      design_matrix <- design_matrix_reactive()
      model_aov <- model_aov_reactive()
      factor_types <- factor_types_reactive()
      
      req(model_formula, design_matrix, model_aov, new_data())
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]
      
      selected_response <- as.character(model_formula)[2]
      
      for (factor_name in main_effects) {
        local({
          fname <- factor_name
          output_id <- paste0("factor_plot_", fname)
          
          output[[output_id]] <- renderPlot({
            req(new_data(), shared_yrange())
            
            column_data <- design_matrix[[fname]]
            factor_type <- types_vec[[fname]]
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
    
    output$selected_value_output <- renderUI({
      model_formula <- model_formula_reactive()
      model_aov <- model_aov_reactive()
      req(model_formula, model_aov, new_data())
      selected_response <- as.character(model_formula)[2]
      
      prediction <- signif(predict(model_aov, newdata = new_data()), 3)
      p("The predicted ", strong(selected_response), " for the selected values is: ", strong(prediction))
    })
    
    # Renders the optimizer parameters panel
    output$model_interaction_output_ui_id <- renderUI({
      design_matrix <- design_matrix_reactive()
      model_formula <- model_formula_reactive()
      factor_types <- factor_types_reactive()
      
      if (is.null(model_formula) || is.null(design_matrix)) {
        return(p("Create and fit a model in ", strong("Analyze the Design"), " to initialize the response optimizer."))
      }
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]
      
      selected_response <- as.character(model_formula)[2]
      
      # Build starting values inputs
      inputs <- lapply(main_effects, function(fname) {
        ftype <- types_vec[[fname]]
        values <- design_matrix[[fname]]
        start_id <- ns(paste0("start_", fname))
        
        if (ftype == "Continuous") {
          numericInput(
            inputId = start_id,
            label = paste("Start", fname),
            value = mean(range(values, na.rm = TRUE))
          )
        } else if (ftype == "Discrete") {
          selectInput(
            inputId = start_id,
            label = paste("Start", fname),
            choices = sort(unique(values)),
            selected = unique(values)[1]
          )
        } else {
          selectInput(
            inputId = start_id,
            label = paste("Start", fname),
            choices = levels(values),
            selected = levels(values)[1]
          )
        }
      })
      
      tagList(
        h4("Select Optimization Objective"),
        radioButtons(ns("opt_type"), "Goal:",
                     choices = c("Maximize", "Minimize", "Find Target"),
                     selected = "Maximize",
                     inline = TRUE
        ),
        conditionalPanel(
          condition = sprintf("input['%s'] == 'Find Target'", ns("opt_type")),
          numericInput(ns("target_response"), "Target Response", value = 0)
        ),
        hr(),
        h4("Starting Values:"),
        inputs,
        actionButton(ns("run_optim_btn"), "Run", class = "btn-primary"),
        hr(),
        h4("Resulting Conditions:"),
        uiOutput(ns("target_solution")),
        h4(paste0("Predicted ", selected_response, ":")),
        uiOutput(ns("target_predicted_response"))
      )
    })
    
    # Executes optim numerical optimization
    observeEvent(input$run_optim_btn, {
      model_formula <- model_formula_reactive()
      design_matrix <- design_matrix_reactive()
      factor_types <- factor_types_reactive()
      model_aov <- model_aov_reactive()
      
      req(model_formula, design_matrix, factor_types, model_aov, input$opt_type)
      
      if (input$opt_type == "Find Target") {
        validate(need(!is.null(input$target_response) && !is.na(input$target_response), "Target response is required."))
      }
      
      types_vec <- as.character(factor_types[1, ])
      names(types_vec) <- colnames(factor_types)
      
      term_labels <- attr(terms(model_formula), "term.labels")
      main_effects <- term_labels[
        !grepl(":", term_labels) &
          grepl("^[[:alnum:]_.]+$", term_labels) &
          term_labels %in% colnames(design_matrix)
      ]
      
      selected_response <- as.character(model_formula)[2]
      
      # Partition variables by type
      continuous_vars <- main_effects[types_vec[main_effects] == "Continuous"]
      discrete_vars <- main_effects[types_vec[main_effects] == "Discrete"]
      categorical_vars <- main_effects[types_vec[main_effects] == "Categorical"]
      noncont_vars <- c(discrete_vars, categorical_vars)
      
      choices_list <- list()
      for (v in noncont_vars) {
        colv <- design_matrix[[v]]
        if (v %in% discrete_vars) {
          choices_list[[v]] <- sort(unique(colv))
        } else {
          choices_list[[v]] <- levels(colv)
        }
      }
      
      if (length(choices_list) > 0) {
        combo_grid <- do.call(expand.grid, c(choices_list, stringsAsFactors = FALSE))
      } else {
        combo_grid <- data.frame(dummy = 1, stringsAsFactors = FALSE)
      }
      
      max_combos <- 5000L
      if (nrow(combo_grid) > max_combos) {
        showNotification(
          paste("Evaluating", nrow(combo_grid), "combinations is large; sampling", max_combos, "rows for speed."),
          type = "warning", duration = 8
        )
        set.seed(1)
        combo_grid <- combo_grid[sample(nrow(combo_grid), max_combos), , drop = FALSE]
      }
      
      cont_start <- numeric(0)
      cont_lower <- numeric(0)
      cont_upper <- numeric(0)
      
      for (v in continuous_vars) {
        colv <- design_matrix[[v]]
        start_val <- input[[paste0("start_", v)]]
        if (is.null(start_val)) {
          start_val <- mean(range(colv, na.rm = TRUE))
        }
        cont_start <- c(cont_start, start_val)
        cont_lower <- c(cont_lower, min(colv, na.rm = TRUE))
        cont_upper <- c(cont_upper, max(colv, na.rm = TRUE))
      }
      names(cont_start) <- continuous_vars
      
      best_obj <- Inf
      best_df <- NULL
      best_pred <- NA_real_
      
      # Evaluate combinations
      for (i in seq_len(nrow(combo_grid))) {
        fixed_vars <- list()
        if (length(noncont_vars) > 0) {
          row_i <- combo_grid[i, , drop = FALSE]
          for (v in noncont_vars) {
            if (v %in% categorical_vars) {
              fixed_vars[[v]] <- factor(row_i[[v]], levels = levels(design_matrix[[v]]))
            } else {
              fixed_vars[[v]] <- row_i[[v]]
            }
          }
        }
        
        # Local objective minimization wrapper
        fn <- function(par) {
          test_row <- as.list(par)
          names(test_row) <- continuous_vars
          full_row <- c(test_row, fixed_vars)
          
          if (length(continuous_vars) == 0) {
            full_row <- fixed_vars
          }
          
          df <- as.data.frame(full_row, stringsAsFactors = FALSE)
          for (n in names(df)) {
            if (is.factor(design_matrix[[n]])) {
              df[[n]] <- factor(df[[n]], levels = levels(design_matrix[[n]]))
            }
          }
          
          pred <- as.numeric(predict(model_aov, newdata = df))
          objective_value(pred, input$opt_type, target = input$target_response)
        }
        
        if (length(continuous_vars) > 0) {
          res <- optim(
            par     = cont_start,
            fn      = fn,
            method  = "L-BFGS-B",
            lower   = cont_lower,
            upper   = cont_upper,
            control = list(maxit = 200)
          )
          sol_cont <- as.list(res$par)
          names(sol_cont) <- continuous_vars
          full_row <- c(sol_cont, fixed_vars)
        } else {
          full_row <- fixed_vars
        }
        
        final_df <- as.data.frame(full_row, stringsAsFactors = FALSE)
        for (n in names(final_df)) {
          if (is.factor(design_matrix[[n]])) {
            final_df[[n]] <- factor(final_df[[n]], levels = levels(design_matrix[[n]]))
          }
        }
        pred_val <- as.numeric(predict(model_aov, newdata = final_df))
        obj_val <- objective_value(pred_val, input$opt_type, target = input$target_response)
        
        if (obj_val < best_obj) {
          best_obj <- obj_val
          best_df <- final_df
          best_pred <- pred_val
        }
      }
      
      output$target_solution <- renderUI({
        vals <- as.list(best_df[1, , drop = FALSE])
        tags$div(
          lapply(names(vals), function(nm) {
            tags$p(
              tags$strong(paste0(nm, ": ")),
              as.character(signif(vals[[nm]], 3))
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
}
