graphToPhylo<-function (graph) 
{
    df <- igraph::as_data_frame(graph)
    node_counts <- table(c(df$to, df$from))
    tips <- names(node_counts)[node_counts == 1]
    nodes <- names(node_counts)[node_counts > 1]
    attr <- igraph::vertex_attr(graph)
    seqs <- attr$sequence
    names(seqs) <- attr$name
    germline <- tips[tips %in% df$from]
    if (length(germline) > 0) {
        ucanode <- paste0(germline, "_UCA")
        nodes <- c(ucanode, nodes)
        df[df$from == germline, ]$from <- ucanode
        row <- c(ucanode, germline, 0)
        names(row) <- c("from", "to", "weight")
        df <- rbind(df, row)
        seqs <- c(seqs, seqs["Germline"])
        names(seqs)[length(seqs)] = paste0(germline, "_UCA")
    }
    tipn <- 1:length(tips)
    names(tipn) <- tips
    noden <- (length(tips) + 1):(length(tips) + length(nodes))
    names(noden) <- nodes
    renumber <- c(tipn, noden)
    df$from <- as.numeric(renumber[df$from])
    df$to <- as.numeric(renumber[df$to])
    phylo <- list()
    phylo$edge <- matrix(cbind(df$from, df$to), ncol = 2)
    phylo$edge.length <- as.numeric(df$weight)
    phylo$tip.label <- tips
    phylo$node.label <- nodes
    phylo$Nnode <- length(nodes)
    phylo$node.label <- nodes
    class(phylo) <- "phylo"
    nnodes <- length(renumber)
    phylo$nodes <- lapply(1:nnodes, function(x) {
        n <- list()
        n$id <- names(renumber[renumber == x])
        n$sequence <- seqs[n$id]
        n
    })
    phylo = ape::ladderize(phylo, right = FALSE)
    return(phylo)
}