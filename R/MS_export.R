###################################################
#    Functions used for export of mass spectra    #
#    from alignment list to matrix of vectors     #
#    with relative mass intensities as variables. #
#    Then used for similarity calculations.       #
##################################################

MS_rewrite <- function(MS){
  tryCatch(out <- {
    #function for transforming a selected spectrum defiend as "(m/z1|intensity1)(m/z2|intensity2)"
    #in to a vector of intensities
    Zk_new<-vector()
    #spliting stringinto a vector of variables where variables on even position
    #are intensities and variables with odd possition are m/z
    Zk <- as.numeric(unlist(strsplit(MS, split = "\\)\\(|\\)|\\(|\\|")))[-1]
    for (x in 1:length(Zk)){
      #selecting variables on even position
      if ((x%%2)==0){
        Zk_new <- append(Zk_new, Zk[x])
        #selecting variables on odd position
      }else{
        #number of zeros that should be appended between detected intensities for undetected m/z
        #for example: if first m/z = 29 we will first append 28 zeros and then first intensity value
        if (x==1){
          Zk_new <- append(Zk_new, numeric(Zk[x]-1))
        }else{
          # for any other place then first position we will append as many zeros
          #as is the difference between previews and next m/z minus one: previos29 next 31 -> 31-29-1=1
          #we will append one zero
          Zk_new <- append(Zk_new, numeric((Zk[x]-Zk[x-2])-1))
        }
      }
    }
    return(Zk_new)
  },
  error = function(e){
    print(
      sprintf("An error occurred in MS_rewrite at %s : %s",
              Sys.time(),
              e)
    )
  })

}


MS_export_tile <- function(Alingment_list, F.names = Alingment_list[,1]){
  tryCatch(out <- {
    #export from alignment list in created by Chromatof Tile
    #with mass spectra export option selected
    #applying MS_rewrite function on every single row of Alingment_list
    # vectors are saved as list
    Tab_export <- lapply(Alingment_list[,ncol(Alingment_list)] , MS_rewrite)
    #making every vector the same length as vector with maximum length
    Tab_export <- Map(function(x, y) append(y, numeric(x)),
                      max(lengths(Tab_export))-lengths(Tab_export),
                      Tab_export)
    #this list is exported to dataframe
    Tab_export <- as.data.frame(do.call(rbind, Tab_export))
    #nameing rows of Tab_export
    rownames(Tab_export) <- F.names
    #nameming columns of Tab_export
    MZ <- 1:ncol(Tab_export)
    colnames(Tab_export) <- paste("m/z", MZ, sep = "_")
    return(Tab_export)
  },
  error = function(e){
    print(
      sprintf("An error occurred in MS_export_tile at %s : %s",
              Sys.time(),
              e)
    )
  })
}


Tab_zk <-read.csv("Smelodi_covid4 - Copy.csv")

MS_export_tile(Tab_zk)
