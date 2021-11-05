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
                        coord = NA, sym = FALSE, mode = "graph") {
  
  # symmetrize the network
  require(sna)
  require(dplyr)
  
  if (sym == TRUE) {
    matrix <- symmetrize(matrix, rule = "weak")  
  }
  
  kcor <- kcores(matrix, mode) %>% 
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
         legend = kcor,
         col = kcor,
         fill = F,
         pch = 19, border = "white",
         title = "k")
  
}


plot_fastgreedy_cd <- function(matrix, title = "", 
                               layout = NULL, legend = TRUE) 
  {
  
  require(igraph)
  require(dplyr)
  # Create igraph from a sociomatrix
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
    plot.igraph(inet, layout = layout, vertex.label = NA)
  }
  # otherwise, use the fructerman-reingold layout
  else {
    fr <- layout.fruchterman.reingold(inet)
    plot.igraph(inet, layout = fr, vertex.label = NA)  
    
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
  
  outdegree <- degree(matrix, cmode = "outdegree")
  indegree <- degree(matrix, cmode = "indegree") 
  betweenness <- betweenness(matrix, cmode) # directed
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
