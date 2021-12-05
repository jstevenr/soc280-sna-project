create_sociomatrix <- function(df, sender_id, receiver_id) {
  
  require(dplyr)
  # counting the nodes (N), creating empty sociomatrix with N x N dimensions
  nodes <- c(df$sender_id, df$receiver_id) %>% unique() %>% sort()
  N <- length(nodes)
  mat <- matrix(0, nrow = N, ncol = N)
  
  # adding rownames and column names to the matrix that match the IDs
  # (since they may not actually be indexed by number)
  rownames(mat) <- nodes
  colnames(mat) <- nodes
  
  
  for (i in 1:N) {
    # temporary data.frame that stores all the receivers for the ith sender
    df_ties <- df %>% 
      filter(sender_id == nodes[i] & tie == 1) %>%
      select(sender_id, receiver_id)
    
    # if the sender has no ties, move on
    if(nrow(df_ties) == 0) {
      next
    }
    
    # iterating across all rows of _this_ data.frame
    # we then add a "1" to the sociomatrix 
    # for each receiver_id for which the sender has a tie
    for (j in 1:nrow(df_ties)) {
      # we enter the sender/receiver ids as a string, rather than by index
      # this is so that we can add the ties by the specific node label
      mat[as.character(nodes[i]), as.character(df_ties$receiver_id[j])] <- 1
    }
  }
  
  return(mat)
}

reorder_blockmodel <- function(bm) {
  bmR <- bm
  bmR$blocked.data[order(bm$block.membership),order(bm$block.membership)]
  bmR$blocked.data <- bm$blocked.data[order(bm$block.membership),order(bm$block.membership)]
  bmR$block.membership <- sort(bmR$block.membership)
  bmR$plabels <- colnames(bmR$blocked.data)
  
  return(bmR)  
}



plot_max_cc <- function(matrix, displaylabels = F, title = "", legend = T, coord = NULL) {
  
  require(sna)
  require(dplyr)
  # symmetrize the matrix
  mat_sym <- symmetrize(matrix, rule = "weak")
  colors <- colors()
  # get clique_census
  clique_census <- clique.census(mat_sym)
  clique_count <- clique_census$clique.count
  
  # get maximum clique size
  max_clique <- apply((clique_count[, -1] > 0) *
                        as.numeric(rownames(clique_count)), 2, max) %>% as.data.frame() %>%
    tibble::rownames_to_column() %>%
    rename(node = "rowname", size = ".") %>%
    mutate(size = as.factor(size))
  sizes <- max_clique$size %>% unique() %>% as.numeric() %>% sort()
  
  # creating a table of nodes per maximum clique size
  table(max_clique$size) %>% as.data.frame() %>% 
    mutate(Proportion = (Freq / sum(Freq)) %>% round(2) ) %>%
    rename("Max Clique Size" = "Var1") %>%
    knitr::kable(caption = "Nodes per max clique size")
  
  # plotting by maximum clique
  
  # if coordinates are provided, use them, otherwise dont
  if (is.matrix(coord)) {
    plot <- gplot(mat_sym,
                  displaylabels = T,
                  usearrows = F, 
                  vertex.col = max_clique$size,
                  vertex.cex = 1.5, label.cex = 0.7, pad = 1,
                  coord = coord)
  }
  else {
    plot <- gplot(mat_sym,
                  displaylabels = T,
                  usearrows = F, 
                  vertex.col = max_clique$size,
                  vertex.cex = 1.5, label.cex = 0.7, pad = 1)
    
  }
  title(title)
  # adding legend if specified (defaults to true)
  if (legend == TRUE) {
    legend("topleft", 
           legend = sizes,
           col = sizes,
           fill = F,
           xjust= 0,
           pch = 19, border = "white",
           title = "Max clique size")
  }
  
  return(plot)
}


plot_kcores <- function(matrix, title = "", displaylabels = F, 
                        coord = NA, sym = FALSE, mode = "graph",
                        cmode = "freeman") {
  
  # symmetrize the network
  require(sna)
  require(dplyr)
  
  if (sym == TRUE) {
    matrix <- symmetrize(matrix, rule = "weak") 
    kcor <- kcores(matrix, mode = "graph", 
                   cmode = cmode)
  }
  else {
    kcor <- kcores(matrix, mode = "digraph", cmode = cmode)
  }
    #as.data.frame()  %>%
    #tibble::rownames_to_column() %>%
    #rename(node = "rowname", kcore = ".") %>%
    #mutate(kcore = as.factor(kcore))
    # kcores <- df_kcor$kcore %>% unique() %>% as.numeric() %>% sort()
  gplot(matrix, displaylabels, 
        usearrows=F, 
        vertex.cex=1.5)
  title(title)
  legend("topleft", 
         legend = unique(kcor),
         col = kcor,
         fill = F,
         pch = 19, border = "white",
         title = "k")
  
}


plot_fastgreedy_cd <- function(matrix, title = "", 
                               layout = NULL, legend = TRUE,
                               vertex.size = 15) 
  {
  
  require(igraph)
  require(dplyr)
  # Create igraph from a sociomatrix
  message("Note: fast-greedy community algorithm is for undirected graphs.")
  inet <- graph.adjacency(matrix, mode = "undirected", diag = F)
  
  # fast.greedy communtiy detection algorithm
  fg <- fastgreedy.community(inet)
  colbar <- rainbow(max(membership(fg))+1)
  V(inet)$color <- colbar[membership(fg) + 1] # setting colors
  
  # getting a membership data.frame
  membership <- fg %>% membership() %>% 
    as.matrix() %>% as.data.frame() %>% 
    tibble::rownames_to_column() %>%
    rename(node = rowname, membership = V1)
  # plotting the sociogram
  # if a layout matrix is provided, use it
  # can be a matrix of vertex coordinates as given by sna::gplot()
  if (is.matrix(layout)) {
    plot.igraph(inet, layout = layout, vertex.label = NA, 
                vertex.size = vertex.size)
  }
  # otherwise, use the fructerman-reingold layout
  else {
    fr <- layout.fruchterman.reingold(inet)
    plot.igraph(inet, layout = fr, 
                vertex.label = NA, vertex.size = vertex.size)  
    
  }
  title(title)
  if (legend == TRUE) {
    legend('topright',
           legend = membership$membership %>% unique() %>% sort(),
           col = colbar,
           fill = F, pch = 19, border = "white",
           title = "Community")
  }
  
}

plot_eigenvector_cd <- function(matrix, title = "", fg, displaylabels = F) {
  
  require(igraph)
  require(network)
  
  inet <- graph.adjacency(matrix, mode = "undirected", diag = F)
  
  #  here's another community detection algorithm: leading eigenvector algorithm
  #  automatically extracts community structure with largest modularity
  lec <- leading.eigenvector.community(inet)
  #  plot the network using colors for memberships
  colbar <- rainbow(max(lec$membership)+1) # identify one color for each community
  V(inet)$color <- colbar[lec$membership+1] # set the colors
  
  # find the density within and between lec communities
  net <- as.network(matrix)
  network::set.vertex.attribute(net,"lc",lec$membership)
  
  gplot(net, usearrows=F, 
        vertex.col=colbar[lec$membership+1], 
        vertex.cex= 1.5,
        displaylabels)
  title(title)
  
  return(lec)
}

# function to calculate indegree, outdegree, betweennness, and eigenvector
get_centralities <- function(matrix, cmode = "directed") {
  
  require(sna)
  
  outdegree <- sna::degree(matrix, cmode = "outdegree")
  indegree <- sna::degree(matrix, cmode = "indegree") 
  betweenness <- sna::betweenness(matrix, cmode) # directed
  eigenvectors <- eigen(matrix)
  # The vector corresponding to the biggest eigenvalue is the 1st column
  eig <- as.numeric(eigenvectors$vectors[,1])
  
  df <- data.frame(outdegree, indegree, betweenness, eigenvector = eig)
  df$id <- row.names(matrix) 
  
  return(df)
  
}


# function to find the correlation of centrality scores within one network
centrality_correlations <- function(matrix, centralities) {
  
  names <- c("outdegree", "indegree", "betweenness", "eigenvector")
  # calculate correlation of ith column vs all other columns
  cormat <- matrix(data = NA, nrow = 4, ncol = 4,
                   dimnames = list(names, names))
  for (i in 1:4) {
    for(j in 1:4) {
      cormat[i,j] <- cor(centralities[,i], centralities[,j])
    }
  }
  
  return(cormat)
}


mc_sim <- function(sociomatrix, n, alpha = 0.05, test) {
  
  require(sna)
  require(ggplot2)
  # monte carlo simulation
  # assumes directed graph
  # creating n graphs of the same dimensions of the input sociomatrix
  # and the same number of edges
  # message("Generating ", n, " networks of size ", sum(sociomatrix), "...")
  
  sim <- rgnm(n, nv = dim(sociomatrix)[1], m = sum(sociomatrix))
  
  if (test == "mutuality") {
    # message("Testing observed level of mutuality...")
    
    # census of observed network
    dyad_census <- sna::dyad.census(sociomatrix)
    # observed mutual ties in the network
    mut_obs <- dyad_census[1] 
    cat("The observed number of mutual ties is ", mut_obs, ".\n", sep = "")
    # census of simulated networks
    sim_dyad_census <- sna::dyad.census(sim)
    # vector of simulated mutual ties, length n
    mut_sim <- sim_dyad_census[,1] 
    
    # count of mutual ties is our test statistic
    # message("Plotting the distribution of the simulations vs the observed value...")
    df <- data.frame(network = 1:n, mutual_ties = mut_sim)
    plot <- ggplot(df, aes(x = mutual_ties)) + 
      geom_histogram(aes(y = ..density..),
                     binwidth = 1,
                     fill = "black",
                     colour = "black",
                     alpha = 0.25) + 
      geom_density(fill = "black", colour = "black", alpha = 0.25) + 
      geom_vline(xintercept = mut_obs, lty = 2, col = "red",
                 show.legend = T) + 
      ggtitle("Distribution of n simulations for # of mutual ties") + 
      ylab("Density of the curve") + 
      xlab("Mutual ties")
    
    # p-value
    # 1 - proportion of networks with more observed mutual times than simulated
    p <- 1 - (sum(mut_obs > mut_sim) / n)
    
  } # end of if statement
  
  else if (test == "transitivity") {
    # message("Testing observed level of transitivity...")
    
    trans_obs <- gtrans(sociomatrix)
    cat("The observed clustering coefficient of the sociomatrix is ", trans_obs,
        ".\n",
        sep = "")
    trans_sim <- gtrans(sim)
    
    # message("Plotting the distribution of the simulations vs the observed value...")
    df <- data.frame(network = 1:n, transitivity = trans_sim)
    plot <- ggplot(df, aes(x = transitivity)) +
      geom_density(fill = "black", colour = "black", alpha = 0.25) + 
      geom_vline(xintercept = trans_obs, lty = 2, col = "red",
                 show.legend = T) + 
      ggtitle("Distribution of n simulations for transitivity (clustering coefficient)") +
      ylab("Density of the curve") + 
      xlab("Transivity")
    
    # p-value
    # 1 - proportion of networks with higher transitivity index than simulated
    p <- 1 - (sum(trans_obs > trans_sim) / n)
    
    
    
  } # end of else if 
  
  else if (test == "three-cycle") {
    # message("Testing observed number of three-cycles...")
    
    # observed cycle census
    cycle_census <- kcycle.census(sociomatrix, maxlen = 3)
    # number of 3-cycles
    three_cycle_obs <- cycle_census$cycle.count["3","Agg"]
    cat("The observed number of three-cycles is ", three_cycle_obs, 
        ".\n",
        sep = "")
    # cycle census of simulated networks
    sim_cycle_census <- kcycle.census(sim, maxlen = 3)
    # initializing empty vector of counts of 3 cycles
    three_cycle_sim <- c()
    # collecting all 3-cycle counts
    for (i in 1:n) {
      three_cycle_sim[i] <- sim_cycle_census[[i]]$cycle.count["3", "Agg"]
    }
    
    # message("Plotting the distribution of the simulations vs the observed value...")
    # data frame of count of three-cycles by network
    df <- data.frame(network = 1:n, three_cycles = three_cycle_sim)
    # plotting histogram and density curve of simulated three-cycle amounts
    plot <- ggplot(df, aes(x = three_cycles)) + 
      # adding histogram
      geom_histogram(aes(y = ..density..),
                     binwidth = 1, 
                     fill = "black",
                     colour = "black",
                     alpha = 0.25) + 
      # adding density curve
      geom_density() + 
      # adding line for the observed value on the x-axis
      geom_vline(xintercept = three_cycle_obs, lty= 2, col = "red",
                 show.legend = T) +
      ggtitle("Distribution of n simulations for # of three-cycle triads") + 
      ylab("Density of the curve") + 
      xlab("Three-cycles")
    
    # p-value
    # 1 - proportion of networks with more observed three cycles than simulated
    p <- 1 - (sum(three_cycle_obs > three_cycle_sim) / n)
    
    
  } # end of else if 
  
  # comparing p-value to significance level alpha, defaults to alpha = 0.05
  # if p-value is less than alpha, reject H0
  # otherwise, fail to reject H0
  
  if (p < alpha) {
    cat("p-value ", p, " < ", alpha, ", we reject the null hypothesis", sep = "")
  }
  else {
    cat("p-value ", p, " > ", alpha, ", we fail to reject the null hypothesis", sep = "")
  }
  
  return(plot)
  
}

# get aic function
get_from_summary <- function(ergm_model, object) {
    
  sum <- summary(ergm_model)  
  
  if (object == "aic") {
    y <- sum$aic[1]
  }
  else if (object == "formula") {
    y <- sum$formula
  }
  
  return(y)
}
