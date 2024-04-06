#@ test input
x <- c("CAACAACTCCAGTACA-1_2_2","TATCCTAGTAATACCC-1_2_1","GACAGCCCACTAAACC-1_2_1")
x <- c("wt_CAACAACTCCAGTACA-1", "mut_TATCCTAGTAATACCC-1")

#* The goal of this function is remove redundant prefix or suffix added to  
#* single cell barcodes in a given list of barocdes. These barcodes were likely
#* reformatted by Seurat merge function that either add project_id as 
#* prefix("projectID_") or add suffix as "_{1-9}+" to barcode. 

#* This function requires a list of barocdes with 16 bp nucleotide. Special
#* characters that separate prefix and suffix from the 16 bp. It will return 
#* a list of barcode with no prefix and suffix that is non-redundant. 

reformatBarcode <- function(x, prefix="_", suffix="-"){
  
  # split barcode into 16bp nt, prefix and suffix 
  # barcode_df <- data.frame(x)%>%
  #   separate(x,c("prefix","tmp"),sep=prefix)%>%
  #   separate(tmp,c("BP16","suffix"),sep=suffix)
  
  #* split barcode into 16bp nt, prefix and suffix and 
  #* return a table of three columns in this order 
  barcode_df <- do.call(rbind, lapply(x, function(x){
    # split a barcode in to 3 parts
    BP16 <- str_extract(x, "[ATCG]{16}")
    
    return(c(BP16, str_split(x,BP16)[[1]]))

  })) %>% data.frame() 
  names(barcode_df) <- c("BP16", "prefix", "suffix")
  ##################################################################################################################################
  
  # get order of samples by prefix
  prefix_list <- unique(barcode_df$prefix)
  uniq_pre <- unique(barcode_df$prefix) # projectID with prefix the special character

  ###* Set new suffix by order of prefix, or 
  if ( length(uniq_pre)!=1 && 
       all(grepl(prefix, uniq_pre)) ){
    
    # get new suffix according to order of samples in prefix_list
    barcode_df$new_suffix <- 0
    
    for (i in 1:length(prefix_list)){
      barcode_df[barcode_df$prefix==prefix_list[i],"new_suffix"] <- i
    }
  } else {
    
  #################################################################
  ###* Set new suffix as part of suffix
    
    ##### remove common string in start of suffix
    sfx_all <- barcode_df$suffix
    
    # create pattern by input suffix 
    # extension_pattern_start <- "^[-|_][0-9]+" 
    extension_pattern_start <- paste0("^", suffix ,"[0-9]+")
    
    # trim start by substring 
    common_pattern <- ""
    while(!is.na(common_pattern)){
      # extract the real common pattern in given format
      common_pattern <- unique(str_extract(sfx_all, extension_pattern_start))
      
      if (length(common_pattern)==1){
        sfx_all <- substring(sfx_all,
                             nchar(common_pattern)+1 ,
                             max(nchar(sfx_all)))
      }
      common_pattern <- unique(str_extract(sfx_all,extension_pattern_start))
      
      if (length(common_pattern)!=1){
        common_pattern <-NA
      }
    }
    
    ##### remove common string in end of suffix
    # extension_pattern_end <- "[-|_][0-9]+$"
    extension_pattern_end <- paste0(suffix,"[0-9]+$")
    n <- max(nchar(sfx_all))
    common_pattern <- substring(sfx_all[1], n-1, n)
    
    # trim end by str_sub 
    # while(!is.na(common_pattern)){
    #   common_pattern <- unique(str_extract(sfx_all, extension_pattern_end))
      
    #   if (length(common_pattern)==1){
    #     sfx_all <- str_sub(sfx_all, end=-nchar(common_pattern))
    #   }
    #   common_pattern <- unique(str_extract(sfx_all,extension_pattern_end))
      
    #   if (length(common_pattern)!=1){
    #     common_pattern <-NA
    #   }
    # }

    while(length(common_pattern) == 1 && !is.na(common_pattern[[1]])) {

      sfx_all <- str_sub(sfx_all, end = -nchar(common_pattern[[1]]))
      common_pattern <- unique(str_extract(sfx_all, extension_pattern_end))
      if (length(common_pattern) != 1) { break }
    }

    # sfx_all is the trimmed/updated new suffix
    uniq_suffice <- unique(sfx_all)
    if (length(uniq_suffice)==1){
      # all suffix is the same 
      barcode_df$new_suffix <- "-1"
    } else {
      #* check if sfx_all at certain index is part of barcode_df$suffix at exact certain index
      match_logic <- mapply(function(X,Y) {
        sapply(length(X), function(index) grepl(X[index], Y[index]))
      }, 
      X=sfx_all, Y=barcode_df$suffix)
      if (all(match_logic)){
        barcode_df$new_suffix <- sfx_all
      }
    }
  }    
  ##################################################################################################################################
  
  # export 
  if (any(barcode_df$new_suffix ==0)){
    print("Error in barcode format")
    return(NULL)
    
  } else {
    #* Check if the new suffix contains special character in the 
    #* beginning. If not, add "-". 
    if ( all(grepl('^[_|@|-]+',barcode_df$new_suffix)) ){
      barcode_df$out_suffix <- gsub('^[_|@]', '-',barcode_df$new_suffix)
    } else {
      barcode_df$out_suffix <- gsub('^', '-',barcode_df$new_suffix)
    }
    #* Check if the new suffix contains special character in the 
    #* end. If so, remove it.
    if (any(grepl('[-|_|@]$',barcode_df$out_suffix))){
      barcode_df$out_suffix <- gsub('[-|_|@]$','',barcode_df$out_suffix, perl=T)
    }
    
    # create new_barocde 
    output_barcode <- barcode_df %>%
      unite(new_barocde, c('BP16','out_suffix'),sep="")%>%
      select(new_barocde)%>%
      as.list()%>%unlist()
    
    if (length(output_barcode)!=length(x)){
      print("Error in barcode format")
      return(NULL)
    } else{
      return(output_barcode)
    }
  }
}
