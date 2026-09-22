#' Universal Network Analyzer Class
#'
#' A flexible and extensible R6 class for network construction, analysis, and comparison
#' that works with various input formats without dependency on specific packages.
#'
#' @import R6
#' @import igraph
#' @export

library(R6)
library(igraph)

NetworkAnalyzer <- R6::R6Class(
  "NetworkAnalyzer",

  # ============================================================================
  # PUBLIC FIELDS
  # ============================================================================
  public = list(
    # Network object
    network = NULL,

    # Original data
    data = NULL,

    # Metadata
    metadata = NULL,

    # Analysis results
    node_attributes = NULL,
    edge_attributes = NULL,
    modules = NULL,
    network_properties = NULL,
    node_roles = NULL,

    # ========================================================================
    # INITIALIZATION
    # ========================================================================

    #' @description
    #' Initialize the NetworkAnalyzer object
    #'
    #' @param data Input data (matrix, data.frame, or igraph object)
    #' @param metadata Optional metadata for nodes/samples
    #' @param data_type Type of input: "abundance", "correlation", "adjacency", "edgelist", "igraph"
    #' @param filter_threshold Abundance threshold for filtering (if applicable)
    #'
    #' @return A new NetworkAnalyzer object
    initialize = function(
      data = NULL,
      metadata = NULL,
      data_type = c(
        "abundance",
        "correlation",
        "adjacency",
        "edgelist",
        "igraph"
      ),
      filter_threshold = 0
    ) {
      data_type <- match.arg(data_type)

      # Store original data
      self$data <- data
      self$metadata <- metadata

      # Process data based on type
      if (!is.null(data)) {
        if (data_type == "igraph") {
          self$network <- data
        } else if (data_type == "abundance") {
          # Filter low abundance features if threshold provided
          if (filter_threshold > 0) {
            data <- private$filter_abundance(data, filter_threshold)
          }
          self$data <- data
        } else if (data_type == "adjacency") {
          self$network <- private$adjacency_to_network(data)
        } else if (data_type == "edgelist") {
          self$network <- private$edgelist_to_network(data)
        }
      }

      message("NetworkAnalyzer initialized successfully")
    },

    # ========================================================================
    # NETWORK CONSTRUCTION
    # ========================================================================

    #' @description
    #' Build network from data
    #'
    #' @param method Network construction method
    #' @param cor_method Correlation method (if applicable)
    #' @param cor_threshold Correlation coefficient threshold
    #' @param p_threshold P-value threshold
    #' @param use_abs Use absolute values for correlation (TRUE keeps both positive and negative)
    #' @param keep_negative Keep negative correlations (only if use_abs=FALSE)
    #' @param ... Additional parameters for specific methods
    #'
    #' @return Self (invisibly)
    build_network = function(
      method = c("correlation", "sparcc", "spieceasi", "custom"),
      cor_method = "pearson",
      cor_threshold = 0.6,
      p_threshold = 0.05,
      use_abs = TRUE,
      keep_negative = TRUE,
      ...
    ) {
      method <- match.arg(method)

      if (method == "correlation") {
        self$network <- private$build_correlation_network(
          cor_method = cor_method,
          cor_threshold = cor_threshold,
          p_threshold = p_threshold,
          use_abs = use_abs,
          keep_negative = keep_negative
        )
      } else if (method == "sparcc") {
        self$network <- private$build_sparcc_network(
          cor_threshold = cor_threshold,
          ...
        )
      } else if (method == "spieceasi") {
        self$network <- private$build_spieceasi_network(...)
      } else if (method == "custom") {
        # Allow custom network building function
        custom_func <- list(...)$custom_func
        if (!is.null(custom_func)) {
          self$network <- custom_func(self$data, ...)
        }
      }

      message(sprintf(
        "Network built with %d nodes and %d edges",
        vcount(self$network),
        ecount(self$network)
      ))

      invisible(self)
    },

    # ========================================================================
    # MODULE DETECTION
    # ========================================================================

    #' @description
    #' Detect modules/communities in the network
    #'
    #' @param method Clustering algorithm
    #' @param use_abs_weight Use absolute values of edge weights (for algorithms that don't support negative weights)
    #' @param ... Additional parameters for the clustering method
    #'
    #' @return Self (invisibly)
    detect_modules = function(
      method = c("louvain", "fast_greedy", "walktrap", "label_prop", "leiden"),
      use_abs_weight = TRUE,
      ...
    ) {
      if (is.null(self$network)) {
        stop("Network not built yet. Run build_network() first.")
      }

      method <- match.arg(method)

      # Handle negative weights
      temp_network <- self$network
      has_negative_weights <- any(E(temp_network)$weight < 0, na.rm = TRUE)

      if (has_negative_weights) {
        if (method %in% c("louvain", "fast_greedy", "leiden")) {
          if (use_abs_weight) {
            message(
              "Detected negative weights. Using absolute values for ",
              method,
              " algorithm."
            )
            E(temp_network)$weight <- abs(E(temp_network)$weight)
          } else {
            warning(
              method,
              " does not support negative weights. Consider using 'walktrap' or set use_abs_weight=TRUE"
            )
          }
        }
      }

      # Select clustering method
      cluster_func <- switch(
        method,
        "louvain" = igraph::cluster_louvain,
        "fast_greedy" = igraph::cluster_fast_greedy,
        "walktrap" = igraph::cluster_walktrap,
        "label_prop" = igraph::cluster_label_prop,
        "leiden" = igraph::cluster_leiden
      )

      # Detect communities
      communities <- cluster_func(temp_network, ...)

      # Store module information in original network
      V(self$network)$module <- paste0("M", membership(communities))
      self$modules <- communities

      # Store sign information if using absolute weights
      if (has_negative_weights && use_abs_weight) {
        message(
          "Note: Module detection used absolute weights. Original edge signs are preserved in edge attributes."
        )
      }

      message(sprintf(
        "Detected %d modules",
        length(unique(membership(communities)))
      ))

      invisible(self)
    },

    # ========================================================================
    # NETWORK PROPERTIES
    # ========================================================================

    #' @description
    #' Calculate network properties and node attributes
    #'
    #' @param calculate_roles Calculate node roles (requires modules)
    #'
    #' @return Self (invisibly)
    calculate_properties = function(calculate_roles = TRUE) {
      if (is.null(self$network)) {
        stop("Network not built yet. Run build_network() first.")
      }

      # Check if modules are detected when calculating roles
      if (calculate_roles && is.null(self$modules)) {
        warning(
          "Modules not detected. Run detect_modules() first to calculate node roles. Skipping role calculation."
        )
        calculate_roles <- FALSE
      }

      # Check for negative weights (affects centrality calculations)
      has_negative_weights <- any(E(self$network)$weight < 0, na.rm = TRUE)
      temp_network <- self$network

      if (has_negative_weights) {
        message(
          "Note: Network has negative weights. Using absolute values for centrality calculations."
        )
        E(temp_network)$weight <- abs(E(temp_network)$weight)
      }

      # Node-level metrics
      node_df <- data.frame(
        node = V(self$network)$name,
        degree = degree(self$network), # Degree doesn't use weights
        betweenness = betweenness(temp_network, normalized = TRUE), # Use temp network
        closeness = closeness(temp_network, normalized = TRUE), # Use temp network
        eigenvector = eigen_centrality(temp_network)$vector, # Use temp network
        clustering_coef = transitivity(self$network, type = "local"), # Doesn't use weights
        stringsAsFactors = FALSE
      )

      # Replace any NA values with 0 in centrality measures
      node_df$closeness[is.na(node_df$closeness)] <- 0
      node_df$clustering_coef[is.na(node_df$clustering_coef)] <- 0

      # Add module information if available
      if (!is.null(self$modules)) {
        node_df$module <- V(self$network)$module
      }

      # Calculate node roles if requested
      if (calculate_roles && !is.null(self$modules)) {
        roles <- private$calculate_node_roles()
        node_df <- cbind(node_df, roles)
      }

      self$node_attributes <- node_df

      # Network-level properties
      # Use temp_network for weighted metrics
      avg_path_length <- tryCatch(
        {
          mean_distance(temp_network, directed = FALSE)
        },
        error = function(e) {
          NA
        }
      )

      self$network_properties <- data.frame(
        nodes = vcount(self$network),
        edges = ecount(self$network),
        density = edge_density(self$network),
        transitivity = transitivity(self$network),
        diameter = diameter(temp_network, directed = FALSE),
        avg_path_length = avg_path_length,
        avg_degree = mean(degree(self$network)),
        modularity = if (!is.null(self$modules)) {
          modularity(self$modules)
        } else {
          NA
        }
      )

      message("Network properties calculated")

      invisible(self)
    },

    # ========================================================================
    # NETWORK COMPARISON
    # ========================================================================

    #' @description
    #' Compare this network with another NetworkAnalyzer object
    #'
    #' @param other Another NetworkAnalyzer object
    #' @param permutations Number of permutations for statistical testing
    #'
    #' @return Comparison results
    compare_with = function(other, permutations = 1000) {
      if (!inherits(other, "NetworkAnalyzer")) {
        stop("other must be a NetworkAnalyzer object")
      }

      if (is.null(self$network) || is.null(other$network)) {
        stop("Both networks must be built before comparison")
      }

      # Compare network properties
      prop_comparison <- private$compare_properties(other)

      # Compare node roles (if available)
      if (!is.null(self$node_roles) && !is.null(other$node_roles)) {
        role_comparison <- private$compare_node_roles(other)
      } else {
        role_comparison <- NULL
      }

      # Compare degree distributions
      degree_comparison <- private$compare_degree_distributions(
        other,
        permutations
      )

      # Compare module structures
      if (!is.null(self$modules) && !is.null(other$modules)) {
        module_comparison <- private$compare_modules(other)
      } else {
        module_comparison <- NULL
      }

      # Identify common and unique nodes
      nodes_self <- V(self$network)$name
      nodes_other <- V(other$network)$name

      node_overlap <- list(
        common = intersect(nodes_self, nodes_other),
        unique_to_self = setdiff(nodes_self, nodes_other),
        unique_to_other = setdiff(nodes_other, nodes_self),
        jaccard_index = length(intersect(nodes_self, nodes_other)) /
          length(union(nodes_self, nodes_other))
      )

      results <- list(
        properties = prop_comparison,
        node_roles = role_comparison,
        degree_distribution = degree_comparison,
        modules = module_comparison,
        node_overlap = node_overlap
      )

      class(results) <- "NetworkComparison"
      return(results)
    },

    # ========================================================================
    # SUBSET NETWORK
    # ========================================================================

    #' @description
    #' Extract a subnetwork based on criteria
    #'
    #' @param nodes Vector of node names to include
    #' @param module Module name to extract
    #' @param min_degree Minimum degree threshold
    #' @param remove_isolates Remove isolated nodes
    #'
    #' @return A new NetworkAnalyzer object with the subnetwork
    subset_network = function(
      nodes = NULL,
      module = NULL,
      min_degree = NULL,
      remove_isolates = TRUE
    ) {
      if (is.null(self$network)) {
        stop("Network not built yet")
      }

      # Start with all nodes
      selected_nodes <- V(self$network)$name

      # Filter by specified nodes
      if (!is.null(nodes)) {
        selected_nodes <- intersect(selected_nodes, nodes)
      }

      # Filter by module
      if (!is.null(module) && !is.null(V(self$network)$module)) {
        module_nodes <- V(self$network)$name[V(self$network)$module == module]
        selected_nodes <- intersect(selected_nodes, module_nodes)
      }

      # Filter by degree
      if (!is.null(min_degree)) {
        high_degree_nodes <- V(self$network)$name[
          degree(self$network) >= min_degree
        ]
        selected_nodes <- intersect(selected_nodes, high_degree_nodes)
      }

      # Create subgraph
      subgraph <- induced_subgraph(self$network, selected_nodes)

      # Remove isolates if requested
      if (remove_isolates) {
        isolated <- which(degree(subgraph) == 0)
        if (length(isolated) > 0) {
          subgraph <- delete_vertices(subgraph, isolated)
        }
      }

      # Create new NetworkAnalyzer object
      subnet <- NetworkAnalyzer$new(
        data = subgraph,
        data_type = "igraph"
      )

      message(sprintf(
        "Subnetwork created with %d nodes and %d edges",
        vcount(subgraph),
        ecount(subgraph)
      ))

      return(subnet)
    },

    # ========================================================================
    # EXPORT FUNCTIONS
    # ========================================================================

    #' @description
    #' Get node attributes table
    #'
    #' @return Data frame of node attributes
    get_node_table = function() {
      if (is.null(self$node_attributes)) {
        message("Calculating properties first...")
        self$calculate_properties()
      }
      return(self$node_attributes)
    },

    #' @description
    #' Get edge attributes table
    #'
    #' @return Data frame of edge attributes
    get_edge_table = function() {
      if (is.null(self$network)) {
        stop("Network not built yet")
      }

      edge_list <- as_edgelist(self$network)
      edge_df <- data.frame(
        from = edge_list[, 1],
        to = edge_list[, 2],
        weight = E(self$network)$weight,
        stringsAsFactors = FALSE
      )

      return(edge_df)
    },

    #' @description
    #' Get adjacency matrix
    #'
    #' @param sparse Return sparse matrix
    #'
    #' @return Adjacency matrix
    get_adjacency_matrix = function(sparse = FALSE) {
      if (is.null(self$network)) {
        stop("Network not built yet")
      }

      return(as_adjacency_matrix(self$network, sparse = sparse))
    },

    #' @description
    #' Export network to file
    #'
    #' @param filename Output filename
    #' @param format File format ("graphml", "gexf", "edgelist", "pajek")
    #'
    #' @return NULL (invisibly)
    save_network = function(
      filename,
      format = c("graphml", "gexf", "edgelist", "pajek")
    ) {
      if (is.null(self$network)) {
        stop("Network not built yet")
      }

      format <- match.arg(format)

      if (format == "graphml") {
        write_graph(self$network, filename, format = "graphml")
      } else if (format == "gexf") {
        if (requireNamespace("rgexf", quietly = TRUE)) {
          rgexf::write.gexf(igraph::get.edgelist(self$network), filename)
        } else {
          stop("Package 'rgexf' required for GEXF format")
        }
      } else if (format == "edgelist") {
        write_graph(self$network, filename, format = "edgelist")
      } else if (format == "pajek") {
        write_graph(self$network, filename, format = "pajek")
      }

      message(sprintf("Network saved to %s", filename))
      invisible(NULL)
    },

    # ========================================================================
    # VISUALIZATION
    # ========================================================================

    #' @description
    #' Plot the network
    #'
    #' @param layout Layout algorithm or coordinates
    #' @param color_by Node attribute for coloring
    #' @param size_by Node attribute for sizing
    #' @param method Plotting method ("igraph", "ggraph", "visNetwork")
    #' @param palette Palette name (from palette_system.R) or color vector.
    #'   NULL uses default palette. Only used when color_by is set.
    #' @param theme_use ggplot2 theme function (ggraph method only).
    #'   NULL uses theme_pub_stat if available, otherwise theme_minimal.
    #' @param ... Additional parameters passed to plotting function
    #'
    #' @return Plot object
    plot_network = function(
      layout = "fr",
      color_by = "module",
      size_by = "degree",
      method = c("igraph", "ggraph", "visNetwork"),
      palette = NULL,
      theme_use = NULL,
      ...
    ) {
      if (is.null(self$network)) {
        stop("Network not built yet")
      }

      method <- match.arg(method)

      if (method == "igraph") {
        p <- private$plot_igraph(layout, color_by, size_by, palette, ...)
      } else if (method == "ggraph") {
        p <- private$plot_ggraph(layout, color_by, size_by, palette, theme_use, ...)
      } else if (method == "visNetwork") {
        p <- private$plot_visnetwork(color_by, size_by, palette, ...)
      }

      return(p)
    },

    # ========================================================================
    # UTILITY FUNCTIONS
    # ========================================================================

    #' @description
    #' Print summary of the network
    print = function() {
      cat("NetworkAnalyzer Object\n")
      cat("======================\n\n")

      if (!is.null(self$network)) {
        cat(sprintf(
          "Network: %d nodes, %d edges\n",
          vcount(self$network),
          ecount(self$network)
        ))

        if (!is.null(self$modules)) {
          cat(sprintf(
            "Modules: %d detected\n",
            length(unique(membership(self$modules)))
          ))
        }

        if (!is.null(self$network_properties)) {
          cat("\nNetwork Properties:\n")
          print(self$network_properties)
        }
      } else {
        cat("Network: Not built yet\n")
      }

      invisible(self)
    }
  ),

  # ============================================================================
  # PRIVATE METHODS
  # ============================================================================
  private = list(
    # Filter low abundance features
    filter_abundance = function(data, threshold) {
      if (is.matrix(data) || is.data.frame(data)) {
        # Calculate relative abundance
        rel_abund <- sweep(data, 2, colSums(data), "/")
        # Keep features above threshold in at least one sample
        keep <- apply(rel_abund, 1, max) >= threshold
        data <- data[keep, , drop = FALSE]
        message(sprintf("Filtered to %d features above threshold", sum(keep)))
      }
      return(data)
    },

    # Convert adjacency matrix to network
    adjacency_to_network = function(adj_matrix) {
      graph_from_adjacency_matrix(
        adj_matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )
    },

    # Convert edge list to network
    edgelist_to_network = function(edgelist) {
      if (ncol(edgelist) >= 3) {
        graph_from_data_frame(edgelist, directed = FALSE)
      } else {
        graph_from_edgelist(as.matrix(edgelist[, 1:2]), directed = FALSE)
      }
    },

    # Build correlation-based network
    build_correlation_network = function(
      cor_method,
      cor_threshold,
      p_threshold,
      use_abs,
      keep_negative = TRUE
    ) {
      # Calculate correlations
      cor_matrix <- cor(
        t(self$data),
        method = cor_method,
        use = "pairwise.complete.obs"
      )

      # Calculate p-values
      n <- ncol(self$data)
      t_stat <- cor_matrix * sqrt(n - 2) / sqrt(1 - cor_matrix^2)
      p_matrix <- 2 * pt(abs(t_stat), n - 2, lower.tail = FALSE)

      # Apply thresholds
      if (use_abs) {
        # Use absolute correlation values
        sig_cor <- abs(cor_matrix) >= cor_threshold & p_matrix < p_threshold

        if (!keep_negative) {
          # Only keep positive correlations
          sig_cor <- sig_cor & cor_matrix > 0
          message(
            "Note: Only positive correlations are kept (set keep_negative=TRUE to include negative correlations)"
          )
        }
      } else {
        # Only positive correlations
        sig_cor <- cor_matrix >= cor_threshold & p_matrix < p_threshold
      }

      # Set diagonal to 0
      diag(sig_cor) <- FALSE

      # Create adjacency matrix
      adj_matrix <- cor_matrix * sig_cor

      # Build network
      network <- graph_from_adjacency_matrix(
        adj_matrix,
        mode = "undirected",
        weighted = TRUE,
        diag = FALSE
      )

      # Add edge sign attribute
      E(network)$sign <- ifelse(E(network)$weight > 0, "positive", "negative")

      # Count positive and negative edges
      n_pos <- sum(E(network)$sign == "positive")
      n_neg <- sum(E(network)$sign == "negative")

      if (n_neg > 0) {
        message(sprintf(
          "Network has %d positive and %d negative edges",
          n_pos,
          n_neg
        ))
      }

      return(network)
    },

    # Build SparCC network (placeholder - requires SpiecEasi or similar)
    build_sparcc_network = function(cor_threshold, ...) {
      stop(
        "SparCC method requires additional package. Please implement or use correlation method."
      )
    },

    # Build SPIEC-EASI network (placeholder)
    build_spieceasi_network = function(...) {
      stop(
        "SPIEC-EASI method requires SpiecEasi package. Please implement or use correlation method."
      )
    },

    # Calculate node roles (Zi-Pi)
    calculate_node_roles = function() {
      module_membership <- membership(self$modules)

      # Calculate within-module degree (Zi)
      Zi <- sapply(1:vcount(self$network), function(i) {
        module <- module_membership[i]
        module_nodes <- which(module_membership == module)

        # Skip if module has only one node
        if (length(module_nodes) <= 1) {
          return(0)
        }

        module_subgraph <- induced_subgraph(self$network, module_nodes)

        # Get degree of current node in the module
        node_name <- as.character(V(self$network)$name[i])
        ki <- degree(module_subgraph)[node_name]

        # Handle case where node not found
        if (is.na(ki)) {
          return(0)
        }

        k_mean <- mean(degree(module_subgraph))
        k_sd <- sd(degree(module_subgraph))

        # Handle zero or NA standard deviation
        if (is.na(k_sd) || k_sd == 0 || is.infinite(k_sd)) {
          return(0)
        }

        zi_value <- (ki - k_mean) / k_sd

        # Handle infinite or NA results
        if (is.na(zi_value) || is.infinite(zi_value)) {
          return(0)
        }

        return(zi_value)
      })

      # Calculate participation coefficient (Pi)
      Pi <- sapply(1:vcount(self$network), function(i) {
        neighbors_idx <- neighbors(self$network, i)

        # No neighbors means Pi = 0
        if (length(neighbors_idx) == 0) {
          return(0)
        }

        module_dist <- table(module_membership[neighbors_idx])
        ki <- length(neighbors_idx)

        # Handle edge case
        if (ki == 0) {
          return(0)
        }

        sum_term <- sum((module_dist / ki)^2)
        pi_value <- 1 - sum_term

        # Handle NA or infinite values
        if (is.na(pi_value) || is.infinite(pi_value)) {
          return(0)
        }

        return(pi_value)
      })

      # Classify roles
      roles <- rep("Peripherals", vcount(self$network))

      # Replace NA values with 0
      Zi[is.na(Zi)] <- 0
      Pi[is.na(Pi)] <- 0

      # Classify based on thresholds
      roles[Zi >= 2.5 & Pi < 0.62] <- "Module hubs"
      roles[Zi < 2.5 & Pi >= 0.62] <- "Connectors"
      roles[Zi >= 2.5 & Pi >= 0.62] <- "Network hubs"

      return(data.frame(
        Zi = Zi,
        Pi = Pi,
        role = roles,
        stringsAsFactors = FALSE
      ))
    },

    # Compare network properties
    compare_properties = function(other) {
      self_prop <- self$network_properties
      other_prop <- other$network_properties

      comparison <- data.frame(
        property = names(self_prop),
        network1 = as.numeric(self_prop),
        network2 = as.numeric(other_prop),
        difference = as.numeric(self_prop) - as.numeric(other_prop),
        pct_change = ((as.numeric(self_prop) - as.numeric(other_prop)) /
          as.numeric(other_prop)) *
          100
      )

      return(comparison)
    },

    # Compare node roles
    compare_node_roles = function(other) {
      self_roles <- table(self$node_attributes$role)
      other_roles <- table(other$node_attributes$role)

      all_roles <- union(names(self_roles), names(other_roles))

      comparison <- data.frame(
        role = all_roles,
        network1 = as.numeric(self_roles[all_roles]),
        network2 = as.numeric(other_roles[all_roles])
      )
      comparison[is.na(comparison)] <- 0

      return(comparison)
    },

    # Compare degree distributions
    compare_degree_distributions = function(other, permutations) {
      deg1 <- degree(self$network)
      deg2 <- degree(other$network)

      # Kolmogorov-Smirnov test
      ks_test <- ks.test(deg1, deg2)

      # Permutation test for mean degree
      observed_diff <- mean(deg1) - mean(deg2)

      perm_diffs <- replicate(permutations, {
        combined <- c(deg1, deg2)
        shuffled <- sample(combined)
        mean(shuffled[1:length(deg1)]) -
          mean(shuffled[(length(deg1) + 1):length(combined)])
      })

      perm_p <- mean(abs(perm_diffs) >= abs(observed_diff))

      return(list(
        ks_test = ks_test,
        mean_diff = observed_diff,
        permutation_p = perm_p,
        distribution1 = deg1,
        distribution2 = deg2
      ))
    },

    # Compare module structures
    compare_modules = function(other) {
      list(
        n_modules_1 = length(unique(membership(self$modules))),
        n_modules_2 = length(unique(membership(other$modules))),
        modularity_1 = modularity(self$modules),
        modularity_2 = modularity(other$modules)
      )
    },

    # Plot with igraph
    plot_igraph = function(layout, color_by, size_by, palette = NULL, ...) {
      # Handle negative weights for layout calculation
      temp_network <- self$network
      has_negative_weights <- any(E(temp_network)$weight < 0, na.rm = TRUE)

      if (has_negative_weights) {
        message("Note: Using absolute edge weights for layout calculation.")
        E(temp_network)$weight <- abs(E(temp_network)$weight)
      }

      # Get layout
      if (is.character(layout)) {
        # Handle different layout naming conventions
        layout_func <- tryCatch(
          {
            # Try new naming convention (layout_with_*)
            get(
              paste0("layout_with_", layout),
              mode = "function",
              envir = asNamespace("igraph")
            )
          },
          error = function(e) {
            # Try old naming convention (layout.*)
            tryCatch(
              {
                get(
                  paste0("layout.", layout),
                  mode = "function",
                  envir = asNamespace("igraph")
                )
              },
              error = function(e2) {
                # If both fail, use a default layout
                warning(
                  "Layout '",
                  layout,
                  "' not found. Using Fruchterman-Reingold instead."
                )
                igraph::layout_with_fr
              }
            )
          }
        )
        coords <- layout_func(temp_network)
      } else {
        coords <- layout
      }

      # Set colors using project palette system when available
      has_palette_sys <- exists("get_colors", mode = "function")
      if (
        !is.null(color_by) && color_by %in% names(vertex_attr(self$network))
      ) {
        colors <- vertex_attr(self$network, color_by)
        if (is.numeric(colors)) {
          if (has_palette_sys) {
            pal <- get_colors(palette, n = 100, type = "continuous")
            vertex_colors <- pal[cut(colors, 100, labels = FALSE)]
          } else {
            vertex_colors <- heat.colors(100)[cut(colors, 100)]
          }
        } else {
          unique_vals <- unique(colors)
          n_vals <- length(unique_vals)
          if (has_palette_sys) {
            color_palette <- get_colors(palette, n = n_vals)
          } else {
            color_palette <- rainbow(n_vals)
          }
          vertex_colors <- color_palette[match(colors, unique_vals)]
        }
      } else {
        vertex_colors <- "lightblue"
      }

      # Set sizes
      if (!is.null(size_by) && size_by %in% names(vertex_attr(self$network))) {
        sizes <- vertex_attr(self$network, size_by)
        vertex_sizes <- scales::rescale(sizes, to = c(3, 15))
      } else {
        vertex_sizes <- 5
      }

      plot(
        self$network,
        layout = coords,
        vertex.color = vertex_colors,
        vertex.size = vertex_sizes,
        vertex.label = V(self$network)$name,
        vertex.label.cex = 0.7,
        edge.arrow.size = 0.5,
        ...
      )
    },

    # Plot with ggraph
    plot_ggraph = function(layout, color_by, size_by, palette = NULL, theme_use = NULL, ...) {
      if (!requireNamespace("ggraph", quietly = TRUE)) {
        stop("Package 'ggraph' required for this plotting method")
      }
      if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' required for this plotting method")
      }

      # Handle negative weights for layout calculation
      temp_network <- self$network
      has_negative_weights <- any(E(temp_network)$weight < 0, na.rm = TRUE)

      if (has_negative_weights) {
        message("Note: Using absolute edge weights for layout calculation.")
        E(temp_network)$weight <- abs(E(temp_network)$weight)
      }

      # Pre-calculate layout coordinates to avoid ggraph internal layout call
      if (is.character(layout)) {
        # Get layout function
        layout_func <- tryCatch(
          {
            get(
              paste0("layout_with_", layout),
              mode = "function",
              envir = asNamespace("igraph")
            )
          },
          error = function(e) {
            tryCatch(
              {
                get(
                  paste0("layout.", layout),
                  mode = "function",
                  envir = asNamespace("igraph")
                )
              },
              error = function(e2) {
                warning(
                  "Layout '",
                  layout,
                  "' not found. Using Fruchterman-Reingold instead."
                )
                igraph::layout_with_fr
              }
            )
          }
        )

        # Calculate coordinates with temp_network (absolute weights)
        coords <- layout_func(temp_network)

        # Use manual layout with original network
        p <- ggraph::ggraph(
          self$network,
          layout = "manual",
          x = coords[, 1],
          y = coords[, 2]
        ) +
          ggraph::geom_edge_link(aes(alpha = 0.5), show.legend = FALSE) +
          ggraph::geom_node_point(aes(color = .data[[color_by]], size = .data[[size_by]])) +
          ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
          ggraph::theme_graph() +
          ggplot2::theme(legend.position = "right")
      } else {
        # Layout is already coordinates
        p <- ggraph::ggraph(self$network, layout = layout) +
          ggraph::geom_edge_link(aes(alpha = 0.5), show.legend = FALSE) +
          ggraph::geom_node_point(aes(color = .data[[color_by]], size = .data[[size_by]])) +
          ggraph::geom_node_text(aes(label = name), repel = TRUE, size = 3) +
          ggraph::theme_graph() +
          ggplot2::theme(legend.position = "right")
      }

      # Apply project color scale if palette system is available
      has_palette_sys <- exists("scale_fill_pub_d", mode = "function")
      if (has_palette_sys) {
        p <- p + scale_color_pub_d(palette)
      }

      # Apply user-specified theme, project theme, or keep ggraph default
      if (is.function(theme_use)) {
        p <- p + theme_use()
      } else if (exists("theme_pub_stat", mode = "function")) {
        p <- p + theme_pub_stat()
      }

      return(p)
    },

    # Plot with visNetwork
    plot_visnetwork = function(color_by, size_by, palette = NULL, ...) {
      if (!requireNamespace("visNetwork", quietly = TRUE)) {
        stop("Package 'visNetwork' required for this plotting method")
      }


      # Prepare nodes data
      nodes_df <- data.frame(
        id = 1:vcount(self$network),
        label = V(self$network)$name
      )

      if (
        !is.null(color_by) && color_by %in% names(vertex_attr(self$network))
      ) {
        nodes_df$group <- vertex_attr(self$network, color_by)
        # Apply project palette for group colors when available
        if (exists("get_colors", mode = "function")) {
          n_groups <- length(unique(nodes_df$group))
          group_colors <- get_colors(palette, n = n_groups)
          nodes_df$color <- group_colors[match(
            nodes_df$group,
            unique(nodes_df$group)
          )]
        }
      }

      if (!is.null(size_by) && size_by %in% names(vertex_attr(self$network))) {
        nodes_df$value <- vertex_attr(self$network, size_by)
      }

      # Prepare edges data
      edges_df <- as_data_frame(self$network, what = "edges")
      edges_df$from <- match(edges_df$from, V(self$network)$name)
      edges_df$to <- match(edges_df$to, V(self$network)$name)

      visNetwork::visNetwork(nodes_df, edges_df, ...) %>%
        visNetwork::visOptions(
          highlightNearest = TRUE,
          nodesIdSelection = TRUE
        )
    }
  )
)


# ============================================================================
# S3 METHODS FOR NetworkComparison
# ============================================================================

#' Print method for NetworkComparison
#' @export
print.NetworkComparison <- function(x, ...) {
  cat("Network Comparison Results\n")
  cat("===========================\n\n")

  cat("Network Properties:\n")
  print(x$properties)

  if (!is.null(x$node_roles)) {
    cat("\nNode Roles Distribution:\n")
    print(x$node_roles)
  }

  cat("\nNode Overlap:\n")
  cat(sprintf("  Common nodes: %d\n", length(x$node_overlap$common)))
  cat(sprintf("  Jaccard index: %.3f\n", x$node_overlap$jaccard_index))

  if (!is.null(x$degree_distribution)) {
    cat("\nDegree Distribution Comparison:\n")
    cat(sprintf("  Mean difference: %.3f\n", x$degree_distribution$mean_diff))
    cat(sprintf(
      "  Permutation p-value: %.3f\n",
      x$degree_distribution$permutation_p
    ))
    cat(sprintf(
      "  KS test p-value: %.3f\n",
      x$degree_distribution$ks_test$p.value
    ))
  }

  invisible(x)
}

#' Plot method for NetworkComparison
#'
#' @param palette Palette name or color vector for fills
#' @param theme_use ggplot2 theme function. NULL uses theme_pub_stat if available.
#' @export
plot.NetworkComparison <- function(
  x,
  type = c("properties", "roles", "degree"),
  palette = NULL,
  theme_use = NULL,
  ...
) {
  type <- match.arg(type)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }

  # Determine theme: user-specified > project > minimal
  has_proj_theme <- exists("theme_pub_stat", mode = "function")
  use_theme <- if (is.function(theme_use)) {
    theme_use
  } else if (has_proj_theme) {
    theme_pub_stat
  } else {
    ggplot2::theme_minimal
  }

  # Determine color scale when palette system is available
  has_palette_sys <- exists("scale_fill_pub_d", mode = "function")

  p <- NULL
  if (type == "properties") {
    p <- ggplot2::ggplot(x$properties, ggplot2::aes(x = property)) +
      ggplot2::geom_col(
        ggplot2::aes(y = network1, fill = "Network 1"),
        position = ggplot2::position_dodge(),
        alpha = 0.7
      ) +
      ggplot2::geom_col(
        ggplot2::aes(y = network2, fill = "Network 2"),
        position = ggplot2::position_dodge(),
        alpha = 0.7
      ) +
      ggplot2::labs(
        title = "Network Properties Comparison",
        y = "Value",
        fill = "Network"
      ) +
      use_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (has_palette_sys) {
      p <- p + scale_fill_pub_d(palette)
    }
  } else if (type == "roles" && !is.null(x$node_roles)) {
    roles_long <- tidyr::pivot_longer(
      x$node_roles,
      cols = c(network1, network2),
      names_to = "network",
      values_to = "count"
    )

    p <- ggplot2::ggplot(
      roles_long,
      ggplot2::aes(x = role, y = count, fill = network)
    ) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::labs(
        title = "Node Roles Distribution Comparison",
        x = "Role",
        y = "Count"
      ) +
      use_theme() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

    if (has_palette_sys) {
      p <- p + scale_fill_pub_d(palette)
    }
  } else if (type == "degree" && !is.null(x$degree_distribution)) {
    deg_data <- data.frame(
      degree = c(
        x$degree_distribution$distribution1,
        x$degree_distribution$distribution2
      ),
      network = rep(
        c("Network 1", "Network 2"),
        c(
          length(x$degree_distribution$distribution1),
          length(x$degree_distribution$distribution2)
        )
      )
    )

    p <- ggplot2::ggplot(
      deg_data,
      ggplot2::aes(x = degree, fill = network)
    ) +
      ggplot2::geom_density(alpha = 0.5) +
      ggplot2::labs(
        title = "Degree Distribution Comparison",
        x = "Degree",
        y = "Density"
      ) +
      use_theme()

    if (has_palette_sys) {
      p <- p + scale_fill_pub_d(palette)
    }
  }

  if (is.null(p)) {
    stop(
      "No data available for type = '", type, "' in this NetworkComparison ",
      "object (e.g. 'roles' needs node_roles from calculate_properties(calculate_roles = TRUE))"
    )
  }

  return(p)
}
