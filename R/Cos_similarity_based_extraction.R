

#' Cos_sim_based_extr
#' @title  Function detects features that are possibly artefacts produced by
#' column bleed, solvent, machine and method imperfection,
#' sorbent decomposition etc.
#'
#' @description Mass spectra of secondary and primary metabolites and other
#' natural products are usually les similar compound to compound and more
#' similar between groups of compounds then mass spectra of impurities witch
#' are highly similar in group and less similar between groups. This is used
#' in this function which searches for groups of mass spectra in the sample that
#' are highly similar inside group and with low similarity outside group.
#'
#'
#'
#' @param align_MS Matrix in which rows are detected peaks and columns are
#' individual masses. Value is mass intensity.
#' @param MS_int_cutoff Mass intensity threshold every intensity below this
#' value will be set to zero.
#' @param Cos_cutoff Threshold for cosine similarity between two mass spectra
#' everything below is set to zero.
#' @param Transitivity_cutoff Transitivity threshold everything below is set to zero
#' @param Group_size Minimum number of detected peaks in group detected in
#' network any group with less members then this number will be drop.
#'
#' @return return table where rows are individual detected peaks second column
#' is number of its group and third column marks if keep the peak or not
#' @export
#'
#' @examples data("Tab_align_ms")
#' @examples Cos_res <- Cos_sim_based_extr(Tab_align_ms, Cos_cutoff = 0.8)
#' @examples summary(Cos_res)
#'
#'
#'
Cos_sim_based_extr <- function(align_MS,
                               MS_int_cutoff = 50, # cutoff for m/z intensity values
                               Cos_cutoff = 0.8,  # cutoff for spectral similarity
                               Transitivity_cutoff = 0.9, # cutoff for Transitivity of groups
                               Group_size = 3){# cutoff for minimum number of members in group
  # Function clustering network based on cosine similarity with edge betweenness algorithm
  #cosine distance calculation
  LECO_export <- align_MS
  LECO_export[is.na(LECO_export)] <- 0
  LECO_export[LECO_export < MS_int_cutoff] <- 0 # drop values of intensities lower then MS_int_cutoff
  LECO_export <- LECO_export[, colSums(LECO_export!=0)>0]# drop all zero columns
  row.names(LECO_export) <- rownames(align_MS)# adding row names
  Cos <- lsa::cosine(t(LECO_export)) # calculating cosine similarity
  Cos[Cos < Cos_cutoff] <- 0 # dropping cosine values lower than Cos_cutoff
  Cos[Cos > 0] <- 1 #setting remainig values > 0 to 1
  diag(Cos) <- 0 # deleting diagonal of 1
  #molecular network
  g <- igraph::graph_from_adjacency_matrix(Cos, mode = "lower", weighted = NULL) #creating network
  eb <- igraph::cluster_edge_betweenness(g) # calculating edge betweenness
  #combining grouping with mass spectra
  LECO_export_group <- cbind(eb$membership, LECO_export) # adding groups to spectral list
  #Transitivity calculation
  Trans_loc <- igraph::transitivity(g, type = "local") # calculating local transitivity of all nodes
  Trans_loc <- data.frame(Transitivity = Trans_loc, Group = eb$membership) # making adding transitivity to every node with group information
  Trans_loc <- Trans_loc[!is.na(Trans_loc$Transitivity), ] # deleting all nodes with transitivity = NA
  Trans_loc <- aggregate(.~Group, data = Trans_loc, mean) # calculating a mean transitivity of groups
  eb_memb_tab <- as.data.frame(table(eb$membership)) # calculating number of nodes in each group
  eb_memb_tab <- eb_memb_tab[eb_memb_tab$Var1 %in% Trans_loc$Group, ]# selecting only those with transitivity
  Trans_loc <- cbind(Trans_loc, eb_memb_tab$Freq)# combinign mean transitivity with group name and number of nodes in group
  #selecting only groups with transitivyty > Transitivity_cutoff and number of nodes >= Group_size
  number_of_ms <- Trans_loc[Trans_loc$Transitivity> Transitivity_cutoff & eb_memb_tab$Freq >= Group_size, 1]
  # creating output list
  output <- cbind(Group_membership = eb$membership, F_ID = rownames(align_MS))
  Cos_score_extracted <- ifelse(output[,1]%in%number_of_ms, "Suspicious_feature", "regular_one")
  output <- cbind(output, Test_results = Cos_score_extracted)
  return(output)
}
