# Function for kmer enrichment with wilcoxon
kmer.wilcox <- function(freq, data) {
  freq <- freq[nchar(rownames(freq)) %in% c(2,3,4,5,6,7), ]
  freq[freq > 0] <- 1
  freq <- as.data.frame(t(freq))

  freq$Score <- data[match(rownames(freq), data$Motif), ]$log2FC.t

# Perform test
  result <- lapply(colnames(freq), function(x) {
    kmer_pos <- freq[freq[, x] == 1, "Score"]
    kmer_neg <- freq[freq[, x] == 0, "Score"]

    if (length(kmer_pos) >= 5 && length(kmer_neg) >= 5) {
      formula_str <- paste("Score ~ ", as.character(x))
      f <- as.formula(formula_str)
    
      w <- wilcox.test(formula = f, data = freq)
    
      data.frame(kmer = x,
               p.val = w$p.value,
               stat = w$statistic,
               estimate = mean(kmer_pos) - mean(kmer_neg) / sd(kmer_pos),
               effect.i = sqrt(w$statistic^2 / (w$statistic^2 + length(c(kmer_neg, kmer_neg)))))
    }
  })
  result <- do.call(rbind, result)
  result$p.val.adj <- p.adjust(result$p.val, method = 'fdr')

# Plot
  result$label <- ifelse(result$p.val.adj < 0.1, result$kmer, "")
  result$significance <- ifelse(result$p.val.adj < 0.1, "FDR < 0.1", "FDR > 0.1")

  g <- ggplot(result,
       aes(y = -log10(p.val.adj), x = estimate,
       label = label)) +
    geom_label_repel(size = 2, max.overlaps = ) +
    geom_point(aes(color = significance), alpha = 0.8, size = 2) +
    scale_color_manual(values = c("FDR < 0.1" = "#6699FF", "FDR > 0.1" = "black")) +
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    labs(title = "K-mer enrichment by Mann-Whitney test")

  print(g)
  return(result)
}

# Function for kmer enrichment with binomial
kmer.binomial <- function(data, kmer_summary,FC_threshold, p.val_threshold) {
  # HITS PPM
  kmers <- kmer_summary[kmer_summary$Score > threshold, ]

  list <- lapply(unique(kmers$Kmer), function(x) {
    kmer.sub <- kmers[kmers$Length == nchar(x), ]
    total <- length(unique(kmer.sub$Motif))
    count <- length(unique(kmer.sub[kmer.sub$Kmer == x, ]$Motif))
    
    data.frame(Motif = x,
               Length = nchar(x),
               Kmer = x,
               Count = count,
               Prob = count / total)
  })
  PPM.hits <- do.call(rbind, list)

# BACKGROUND PPM
  kmers <- kmer_summary
  list <- lapply(unique(kmers$Kmer), function(x) {
    kmer.sub <- kmers[kmers$Length == nchar(x), ]
    total <- length(unique(kmer.sub$Motif))
    count <- length(unique(kmer.sub[kmer.sub$Kmer == x, ]$Motif))
    
    data.frame(Motif = x,
               Length = nchar(x),
               Kmer = x,
               Count = count,
               Prob = count / total)
  })
  PPM.background <- do.call(rbind, list)

# Perform Binomial to determine enrichment of kmers
  common_kmers <- intersect(PPM.hits$Kmer, PPM.background$Kmer)
  
  list <- lapply(common_kmers, function(x) {
    hits_count <- sum(PPM.hits[PPM.hits$Kmer == x, ]$Count)
    background_count <- sum(PPM.background[PPM.background$Kmer == x, ]$Count)
    n <- sum(PPM.hits$Count)
    prob <- background_count / sum(PPM.background$Count)
    
    binom <- binom.test(hits_count, n, prob)
    
    data.frame(Position = p, Kmer = x, p.val = binom$p.value, estimate = binom$estimate)
    })
  binomial_test <- do.call(rbind, list)

  binomial_test$p.val.adj <- p.adjust(binomial_test$p.val, method = 'fdr')

# Visualize
  binomial_test$label <- ifelse(binomial_test$p.val.adj < p.val_threshold,
                              binomial_test$Kmer, "")
  binomial_test$significance <- ifelse(binomial_test$p.val.adj < p.val_threshold, "FDR < 0.1", "FDR > 0.1")

  g <- ggplot(binomial_test[nchar(binomial_test$Kmer)  %in% c(2,3,4,5,6), ],
       aes(y = -log10(p.val), x = log10(estimate + 1)*10000,
       label = label)) +
    geom_label_repel(size = 2) +
    geom_vline(xintercept = mean(log10(binomial_test$estimate+1)*10000), linetype = "dotted") +
    geom_vline(xintercept = mean(log10(binomial_test$estimate+1)*10000) + 1*sd(log10(binomial_test$estimate+1)*10000), linetype = "dotted") +
    geom_point(aes(color = significance), alpha = 0.8, size = 2) +
    scale_color_manual(values = c("FDR < 0.1" = "#6699FF", "FDR > 0.1" = "black")) +
    #annotate("text", x = 1, y = 31, label = paste0("Mean: ", round(log10(binomial_test$estimate+1)*10000), 3)) +
    theme_bw() + 
    theme(panel.grid = element_blank()) + 
    labs(title = "K-mer enrichment by binomial test", subtitle = "Threshold: 1 SD from mean")

  print(g)
  return(binomial_test)
}


# Function for applying weights to kmer matrix
weight <- function(freq, test, size, p.val, threshold) {
  # convert to binary
  freq <- freq[nchar(rownames(freq)) %in% size, ]
  freq[freq > 0] <- 1
  freq <- t(freq)

  # Apply score associated enrichment weights
  weights <- freq
  for (m in rownames(weights)) {
    kmers <- colnames(weights)
    kmers <- intersect(kmers, test[test$p.val.adj < p.val, ]$Kmer)
    for (kmer in kmers) {
      if (weights[m, kmer] == 1) {
        enrich_score <- test[test$Kmer == kmer, ]$estimate
        weights[m, kmer] <- as.numeric(weights[m, kmer]) * enrich_score
      }
      }
  }
  #Removed kmers not owned by any kmer
  weights <- weights[, colSums(weights) != 0]

  return(weights)
}


# UMAP and plotting function
umap.motifs <- function(Matrix, data, k, threshold, metric, label) {
  umap <- umap(kmer_weights, n_neighbors = k, metric = metric, approx_pow = TRUE, ret_nn = TRUE)

  umap <- as.data.frame(umap$embedding)
  umap$Score <- data[data$Motif %in% rownames(umap), ]$log2FC.t
  umap$Score <- umap$Score + abs(min(data[data$Motif %in% rownames(umap), ]$log2FC.t))
  
  umap$Score.i <- data[data$Motif %in% rownames(umap), ]$log2FC
  umap$Score.i <- umap$Score + abs(min(data[data$Motif %in% rownames(umap), ]$log2FC))
  
  umap$label <- ifelse(umap$Score > threshold, as.character(rownames(umap)), "")

  g <- ggplot(umap, aes(x = umap$V1, y = umap$V2, color = Score, label = label)) +
    geom_point(alpha = 0.6, size = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("Dim 1") + ylab("Dim 2") +
    scale_color_gradient2(low = "#CCCCCC", mid = "#CCCCCC", high = "#100CCF",
                          midpoint = mean(umap$Score)) + 
    labs(title = paste0("k: ", k))
  if (label == TRUE) {
    g <- g + geom_label_repel(max.overlaps = Inf, size = 2, color = "black", force = 2)
  }
  
  print(g)
  return(umap)
}

# K-nn clustering with options for k, resolution, and output matrix
knn.clustering <- function(matrix, umap, k, r) {
  knn <- get.knn(matrix, k = k)
  m <- Matrix(0, nrow = nrow(matrix), ncol = nrow(matrix), sparse = TRUE)

  for (i in 1:nrow(knn$nn.index)) {
    for (j in 1:ncol(knn$nn.index)) {
      if (knn$nn.index[i, j] != 0) {
        m[i, knn$nn.index[i, j]] = 1
    }
  }
}
  m %>% 
    graph_from_adjacency_matrix(mode = "undirected") %>%
    cluster_leiden(resolution_parameter = r, n_iterations = 1000) -> leiden
  umap$leiden <- as.factor(leiden$membership)

  g <- ggplot(umap, aes(x = umap$V1, y = umap$V2, color = leiden)) +
    geom_point(alpha = 0.6, size = 2) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    xlab("Dim 1") + ylab("Dim 2") + 
    labs(title = paste0("r: ", r))
  
  print(g)
  return(umap)
}

# Function for determining the number of nt differences of a string to reference
find_diff <- function(input, reference) {
  num_pos <- nchar(reference)
  num_changes <- 0
  for (p in 1:num_pos) {
    if (substr(input, p, p) != substr(reference, p, p)) {
      num_changes <- num_changes + 1
    }
  }
  return(as.numeric(num_changes))
}
