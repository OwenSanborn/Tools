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
