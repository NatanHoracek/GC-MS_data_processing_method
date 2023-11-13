SMF <- function(vec1, vec2){

  tryCatch(out <- {
    # function for calculating simple match factor between two vectors
    smf <- ((sum((vec1^(1/2))*(vec2^(1/2))))^2)/(sum(vec1)*sum(vec2))
    return(smf)
  },
  error = function(e){
    print(
      sprintf("An error occurred in SMF at %s : %s",
              Sys.time(),
              e)
    )
  })
}

SMF_mat <- function(x){

  tryCatch(out <- {
    # function for calculating simple match factor form a matrix
    co = array(0, c(ncol(x), ncol(x)))
    f = colnames(x)
    dimnames(co) = list(f, f)
    for (i in 2:ncol(x)) {
      for (j in 1:(i - 1)) {
        co[i, j] = SMF(x[, i], x[, j])
      }
    }
    co = co + t(co)
    diag(co) = 1
    return(as.matrix(co))
  },
  error = function(e){
    print(
      sprintf("An error occurred in SMF_mat at %s : %s",
              Sys.time(),
              e)
    )
  })
}

#SMF based clusts and f. selection
SMF_based_extr <- function(align_MS,
                           MS_int_cutoff = 50, # cutoff for m/z intensity values
                           Smf_cutoff = 0.65, # cutoff for spectral similarity
                           Transitivity_cutoff = 0.8,# cutoff for Transitivity of groups
                           Group_size = 3){# cutoff for minimum number of members in group
  tryCatch(out <- {

    #SMF distance calculation
    LECO_export <- align_MS
    LECO_export[LECO_export < 50] <- 0# drop values of intensities lower then MS_int_cutoff
    LECO_export <- LECO_export[, colSums(LECO_export!=0)>0]# drop all zero columns
    row.names(LECO_export) <- rownames(align_MS)# adding row names
    Smf <- SMF_mat(t(LECO_export))# calculating SMF
    Smf[Smf < Smf_cutoff] <- 0# dropping SMF values lower than Smf_cutoff
    diag(Smf) <- 0# deleting diagonal of 1
    #molecular betwork plotng
    g <- igraph::graph_from_adjacency_matrix(Smf, mode = "lower", weighted = "weight")#creating network
    eb <- igraph::cluster_edge_betweenness(g)# calculating edge betweenness
    #combining grouping with mass spectra
    LECO_export_group <- cbind(eb$membership, LECO_export)# adding groups to spectral list
    #Transitivity calculation
    Trans_loc <- igraph::transitivity(g, type = "local") # calculating local transitivity of all nodes
    Trans_loc <- data.frame(Transitivity = Trans_loc, Group = eb$membership)# making adding transitivity to every node with group information
    Trans_loc <- Trans_loc[!is.na(Trans_loc$Transitivity), ]# deleting all nodes with transitivity = NA
    Trans_loc <- aggregate(.~Group, data = Trans_loc, mean) # calculating a mean transitivity of groups
    eb_memb_tab <- as.data.frame(table(eb$membership))# calculating number of nodes in each group
    eb_memb_tab <- eb_memb_tab[eb_memb_tab$Var1 %in% Trans_loc$Group, ]# selecting only those with transitivity
    Trans_loc <- cbind(Trans_loc, eb_memb_tab$Freq)# combinign mean transitivity with group name and number of nodes in group
    #selecting only groups with transitivyty > Transitivity_cutoff and number of nodes >= Group_size
    number_of_ms <- Trans_loc[Trans_loc$Transitivity> Transitivity_cutoff & eb_memb_tab$Freq >= Group_size, 1]
    # creating output list
    output <- cbind(Group_membership = eb$membership, F_ID = rownames(align_MS))
    SMF_score_extracted <- ifelse(output[,1]%in%number_of_ms, "Suspicios_feature", "regular_one")
    output <- cbind(output, Test_results = SMF_score_extracted)
    output <- list("Table" = output, "Num_of_edges" = igraph::gsize(g))
    return(output)

  },
  error = function(e){
    print(
      sprintf("An error occurred in MS_rewrite at %s : %s",
              Sys.time(),
              e)
    )
  })
}


#combining SMF and cos group search

Cos_SMF_comb <- function(cos_extract, SMF_extract){
  tryCatch(out <- {
    ID_in_cos <- cos_extract[cos_extract[,3]=="Suspicios_feature",2]
    SMF_g_in_cos <- SMF_extract[SMF_extract[,2] %in% ID_in_cos, 1]
    SMF_group_names <- as.data.frame(table(SMF_g_in_cos))[,1]
    ID_of_suspitios_features <- unique(c(ID_in_cos, SMF_extract[SMF_extract[,1] %in% SMF_group_names,2]))
    return(ID_of_suspitios_features)
  },
  error = function(e){
    print(
      sprintf("An error occurred in Cos_SMF_comb at %s : %s",
              Sys.time(),
              e)
    )
  })
}



tryCatch(out <- {

},
error = function(e){
  print(
    sprintf("An error occurred in MS_rewrite at %s : %s",
            Sys.time(),
            e)
  )
})

