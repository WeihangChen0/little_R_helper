# require(BiocManager)
# require(stringr)

installPkgs <- function(list_packages, update_R_pkgs=F, 
                        install_to_libPath=.libPaths()[1]){
  
  # get bioconductor version installed for current session 
  tmp <- BiocManager::repositories()
  bioconductor_version <- stringr::str_extract(tmp["BioCsoft"], "\\/[0-9]+\\.[0-9]+\\/")
  bioconductor_version <- gsub("/$","", gsub("^/","", bioconductor_version))
  
  
  # get list of packages that need to be installed  
  if (!update_R_pkgs){
    new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
  } else {
    new.packages <- list_packages
  }
  
  # 
  p_err<-c()
  for(pkg in new.packages){
    try(
      # install pkg 
      BiocManager::install(pkg,version=bioconductor_version,
                            force = T, update = update_R_pkgs,
                            lib=install_to_libPath)
      )
    
    if (! pkg%in%installed.packages()[,"Package"]){
      message(paste0("Error in installing : ", pkg))
      p_err<-c(p_err, pkg)
    }
  }
  return(p_err)
}

# test run 
list.of.packages <- c("testHere","KEGGgraph","siggenes","MSnbase","multtest","RBGL")
check_pkg <- installPkgs(list.of.packages)
print(check_pkg)
