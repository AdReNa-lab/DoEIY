FROM rocker/shiny:4.4.0

# Install system dependencies required by R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libxt-dev \
    zlib1g-dev \
    libglpk-dev \
    make \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
# Explicitly use the Posit Package Manager (P3PM) binary URL to download precompiled Linux packages
RUN R -e "install.packages(c( \
    'shiny', \
    'shinydashboard', \
    'shinydashboardPlus', \
    'shinyWidgets', \
    'gfonts', \
    'waiter', \
    'ggplot2', \
    'rhandsontable', \
    'data.table', \
    'mltools', \
    'broom', \
    'plotly', \
    'rlang', \
    'FrF2', \
    'rsm', \
    'lhs', \
    'DoE.base', \
    'AlgDesign', \
    'dplyr' \
    ), repos='https://packagemanager.posit.co/cran/__linux__/jammy/latest')"

# Set working directory
WORKDIR /app

# Copy the app files into the container
COPY DoEIY/ /app/

# Expose port 3838 for the Shiny application
EXPOSE 3838

# OpenShift compatibility: Allow the root group (GID 0) to read/write in /app
RUN chgrp -R 0 /app && chmod -R g=u /app

# Run the Shiny application directly
CMD ["R", "-e", "shiny::runApp('/app', host='0.0.0.0', port=3838)"]
