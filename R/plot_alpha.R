#' Diagnostic Plot
#'
#' Number of Components, Component ranking and trace plot of alpha values.
#'
#' @param x object class 'dpmm_fit'.
#' @param nburn Iterations to be discarded in the plot of component rankings
#' @param thinning Iterations 
#' @param ... other parameters used by ggplot2.
#'
#' @return A 'ggplot' panel plot.
#'
#' @import cowplot
#' @importFrom purrr discard
#'
#' @examples
#' \dontrun{
#' ## load dataset
#' data(dataset_1)
#' 
#' ## fit model
#' posteriors <- runModel(dataset_1, 
#'                        mcmc_iterations = 100,
#'                        L = 6, 
#'                        mcmc_chains = 2, 
#'                        standardise = TRUE)
#'                
#'                               
#' ## diagnostic plot
#' plot_alpha(posteriors, nburn = 50, thinning = 1)
#' }
#'
#'
#' @export
plot_alpha <- function(x, nburn, thinning, ...) {
  
  #:---------------------------------------------------------------
  ## Check before running function
  if(missing(x)) stop("'x' needs to be supplied")
  if (!inherits(x, "dpmm_fit")) stop("'x' needs to be of class 'dpmm_fit'")
    
  if(missing(nburn)) stop("'nburn' must be provided")
  if(!missing(nburn)) {
    if (is.numeric(nburn) == FALSE) {
      stop("'nburn' must be numeric")
    }
  }
  
  if(missing(thinning)) stop("'thinning' needs to be supplied")
  if(!missing(thinning)) {
    if (is.numeric(thinning) == FALSE) {
      stop("'thinning' must be numeric")
    }
  }
  
  #:---------------------------------------------------------------
  ## Generate objects that will contain the information for plots
  summary.clusters_complete <- NULL
  postComp_complete <- NULL
  postAlpha_complete <- NULL
  
  #:---------------------------------------------------------------
  ## Iterate through all chains
  for (chain in 1:x$mcmc_chains) {
    
    title_chain <- paste0("Chain ",chain)
    
    if (x$mcmc_chains == 1) {
      samples <- x$samples %>%
        as_tibble()
    } else {
      samples <- x$samples[[chain]] %>%
        as_tibble()
    }
      
    iterations <- seq(nburn,nrow(samples), thinning)
    
    
    #:---------------------------------------------------------------
    ## Plot A - Number of Components
    postComp <- samples %>% 
      select(starts_with("z")) %>%
      apply(1, function(x)  length(unique(x))) %>%
      as.data.frame() %>%
      cbind(Iteration = rep(1:nrow(samples), 1)) %>%
      cbind(Chain = rep(title_chain, nrow(samples)))
    
    postComp_complete <- rbind(postComp_complete, postComp)
    
    #:---------------------------------------------------------------
    ## Plot B - Average number of individuals for ranked components
    samples_summary <- samples %>%
      select(starts_with("z"))
    
    summary.clusters <- NULL
    for (i in iterations) {
      row_summary <- samples_summary[i,] %>% t() %>% factor() %>% summary()
      summary.clusters <- dplyr::bind_rows(summary.clusters, row_summary)
    }
    
    summary.clusters <- summary.clusters %>%
      apply(1, function(x) x[order(x, decreasing = TRUE)]) %>%
      t() %>%
      as.data.frame() %>%
      discard(~all(is.na(.x))) %>%
      mutate_all(~replace(., is.na(.), 0)) %>%
      colMeans() %>%
      as.data.frame() %>%
      setNames(c("value"))
    
    summary.clusters <- summary.clusters %>%
      mutate(key = rep(as.character(1:nrow(summary.clusters)),1)) %>%
      mutate(key = factor(key, levels = as.character(1:nrow(summary.clusters)))) %>%
      cbind(Chain = rep(title_chain, nrow(summary.clusters)))
    
    summary.clusters_complete <- rbind(summary.clusters_complete, summary.clusters) 
    
    #:---------------------------------------------------------------
    ## Plot C - Alpha values
    postAlpha <- samples %>%
      select("alpha") %>%
      mutate(Iteration = rep(1:nrow(samples), 1)) %>%
      cbind(Chain = rep(title_chain, nrow(samples)))
    
    postAlpha_complete <- rbind(postAlpha_complete, postAlpha)
    
      
  }
  
  #:---------------------------------------------------------------
  ## Formatting datasets
  summary.clusters_complete <- summary.clusters_complete %>%
    mutate(Chain = factor(Chain))
  
  postComp_complete <- postComp_complete %>%
    mutate(Chain = factor(Chain))
  
  postAlpha_complete <- postAlpha_complete %>%
    mutate(Chain = factor(Chain))
  
  
  #:---------------------------------------------------------------
  ## Plot
  plot <- plot_grid(
    
    plot_grid(
      
      postComp_complete %>%
        ggplot() +
        geom_vline(xintercept = nburn, colour = "black", linetype = "dashed") +
        geom_path(aes(x = Iteration, y = `.`, colour = Chain), linewidth = 0.4, alpha = 0.7) +
        theme_bw() +
        scale_y_continuous(breaks = seq(1,200, by =1)) +
        labs(title = "Number of Components",
             x = "Iterations",
             y = "Components") +
        theme(legend.position = "none")
      
      ,
      
      summary.clusters_complete %>%
        ggplot() +
        geom_col(aes(x = key, y = value, colour = Chain, fill = Chain), linewidth = 0.2, position = "dodge2") +
        theme_bw() +
        labs(title = "Average Number of Individuals for Ranked Components",
             x = "Component Ranking",
             y = "Individuals") +
        theme(legend.position = "none")
      ,
      
      ncol = 2, rel_widths = c(15,28), labels = c("A","B")
      
      
    )
    
    ,
    
    postAlpha_complete %>%
      ggplot() +
      geom_vline(xintercept = nburn, colour = "black", linetype = "dashed") +
      geom_line(aes(x = Iteration, y = alpha, colour = Chain), alpha = 0.7) +
      theme_bw() +
      labs(title = "Alpha values",
           x = "Iterations") +
      theme(axis.title.y = element_blank()) +
      theme(legend.position = "none")
    ,
    
    ncol = 1, nrow = 2, rel_heights = c(23,16), labels = c("","C")
    
    
  )
  
  return(plot)
}



