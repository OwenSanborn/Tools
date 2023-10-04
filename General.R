### Inputs a list of tables and merges them by binding columns with row matching rownames.
multimerge <- function (mylist) {
 
  unames = unique(unlist(lapply(mylist, rownames)))
 
  n = length(unames)
 
  out = lapply(mylist, function(df) {
 
    tmp = matrix(nr = n, nc = ncol(df), dimnames = list(unames,colnames(df)))
    tmp[rownames(df), ] = as.matrix(df)
    rm(df); gc()
 
    return(tmp)
  })
 
  stopifnot( all( sapply(out, function(x) identical(rownames(x), unames)) ) )
 
  bigout = do.call(cbind, out)
  colnames(bigout) = paste(rep(names(mylist), sapply(mylist, ncol)), unlist(sapply(mylist, colnames)), sep = "_")
  return(data.frame(bigout))
}

### Inputs a list of files and a filepath and outputs a table. Designed for EditR html reports. File must be named Group#.Rep# (ei. 1.1)
html_to_table <- function (mylst, filepath) {
  for (file in mylist) {
      html <- read_html(paste0(filepath, file, ".html"))
      xpath_value <- html_text(html_element(html, xpath = '//*[@id="for-use-in-r"]/pre[2]/code'))
      xpath_value <- gsub("\n", "", xpath_value)
      xpath_value <- gsub("\\\\", "", xpath_value)
     
      data_frame <- eval(parse(text = xpath_value))

      group <- gsub("\\..*", "", file)
      replicate <- gsub("^.*\\.","", file )
  
      data_frame$Group <- group
      data_frame$Rep <- replicate
      data_frames[[file]] <- data_frame
   }
   combined_df <- do.call(rbind, data_frames)
   return(combined_df)
}

 ### Generates a random nucleotide (DNA) for a nucleotide given. Returns nt at frequency freq and a random nt at (1 - freq) / 3                      
random_base <- function(nt, freq) {
  nucleotides <- c("A", "T", "C", "G")
  nucleotides <- nucleotides[nucleotides != nt]
  freq.i <- (1 - freq) / 3
  r <- runif(1)
  if (r < freq) {
    return(nt)
  } else if (r < (freq + freq.i) ) {
    return(nucleotides[[1]])
  } else if (r < (freq + freq.i*2) ) {
    return(nucleotides[[2]])
  } else {
    return(nucleotides[[3]])
  }
}

 ### Find the number of nt changes in a string.                       
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

# Output position of nucleotide changes in a string
find_diff_pos <- function(input, reference) {
  num_pos <- nchar(reference)
  loc_changes <- list()
  for (p in 1:num_pos) {
    if (substr(input, p, p) != substr(reference, p, p)) {
      loc_changes[[p]] <- p
    }
  }
  loc_changes <- Filter(function(x) !is.null(x), loc_changes)
  return(as.numeric(unlist(loc_changes)))
}
                         
### Takes Genes x Cells counts matrix and calculates log(CV^2)                         
raw_variance <- function(counts) {
  gene_names <- rownames(counts)
  counts.raw <- lapply(gene_names, function(gene_name) {
    cts <- counts[gene_name, ]
    if(sum(cts != 0) > 2) {
      CV <- log10((sd(cts[cts != 0])/mean(cts[cts != 0])^2))
      Mean <- mean(cts[cts != 0])
      data.frame(gene = gene_name, CV = CV, Mean = Mean)} })

  Raw.counts <- do.call("rbind", counts.raw)
  return(Raw.counts) 
}

### Takes Genes x Cells counts matrix and calculates relative counts per million and log(CV^2)                         
adj_variance <- function(raw) {
  cell_names <- colnames(raw)
  
  raw.list <- lapply(cell_names, function(cell) {
    cts <- raw[ , cell]
    total.transcripts <- sum(cts)
    cts.rc.adj <- (cts / total.transcripts) * 1e6
    data.frame(row.names = rownames(raw), cell = cts.rc.adj) })
  
  raw.RC <- do.call("cbind", raw.list)
  
  colnames(raw.RC) <- cell_names
  
  counts.adj <- lapply(rownames(raw.RC), function(gene_name) {
    cts <- raw.RC[gene_name, ]
    if(length(cts != 0) > 2) {
      x <- as.numeric(cts[cts != 0])
      CV <- log10((sd(x)/mean(x))^2)
      Mean <- mean(x)
      
      data.frame(gene = gene_name, CV = CV, Mean = Mean)}})
  
  adj.counts <- do.call("rbind", counts.adj)
  print("CV Adjustment Complete")
  return(adj.counts)
}
### Applies spline model to relative cpm matrix and normalizes varaiance with residuals
std_variance <- function(adj.counts) {
  adj.counts <- adj.counts[!(adj.counts$CV %in% c(Inf, -Inf, NaN, NA)), ]
  adj.counts <- adj.counts[!(adj.counts$Mean %in% c(Inf, -Inf, NaN, NA)), ]
  
  model <- lm(data = adj.counts, formula = CV ~ ns(Mean, df = 20))
  
  #model <- smooth.spline(adj.counts$Mean, adj.counts$CV, cv = TRUE)

  residuals <- residuals(model)
  mse <- sum(residuals^2) / (length(residuals) - length(model$coefficients))
  standardized_residuals <- residuals / sqrt(mse)

  adj.counts$std_res <- standardized_residuals
  print("Residual Calculations Complete")
  return(adj.counts)
}
