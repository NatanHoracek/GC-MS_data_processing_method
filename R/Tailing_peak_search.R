###################################################
#   Tailing function detecting and combining
#   tailing peaks split by peak detection
#
###################################################


Tailing_peak_selection <- function(Clique_list,  #list of vectors of peak names belonging to one clique
                                   sample_list,  #list of smaples with 1d and 2d retention times and peak names as rownames
                                   fig = "G23_52_WB_results_0.9cuttoff.pdf", #name of figure with ploted cliques and wrongones highlighted
                                   fiq_w = 60, #diametrs of clique
                                   fig_h = 40){ #diametrs of clique
  #setting up PDF file size and layout
  pdf(file = fig, width = fiq_w, height = fig_h)
  layout(matrix(c(1:55), nrow = 5, ncol = 11)) ####################################
  Cliques_selected <- list() # feature vector of selected cliques
  list_of_models <- list() #list of fitted models information
  #loop looking for tailing peaks
  for (x in 1:length(Clique_list)){
    RT <- as.vector(names(Clique_list[[x]])) #selecting vector of peak names in one clique
    RT <- RT[!is.na(RT)]#deleting NA if present
    Vec <- sample_list[rownames(sample_list) %in% RT,1:2]#getting retention times of selected peaks
    x_axis <- seq(min(Vec[,1]), max(Vec[,1]), length=10)#selecting x axis bins
    y_axis_diff <- max(Vec[,2]) - min(Vec[,2])#getting y axis length
    # estimating coeficients for model
    theta.0 <- min(Vec$RT_2st_dim) * 0.5
    model.0 <- lm(log(RT_2st_dim - theta.0) ~ RT_1st_dim, data=Vec)
    alpha.0 <- exp(coef(model.0)[1])
    beta.0 <- coef(model.0)[2]
    Start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)
    # testing if nonlinear model is possible
    ## if it is not possible
    if (class(try(nls(RT_2st_dim ~ alpha * exp(beta * RT_1st_dim) + theta, data = Vec, start = Start), silent = TRUE))== "try-error"){
      # fitting linear model
      fit_1 <- lm(RT_2st_dim ~ RT_1st_dim, data=Vec)
      Pred <-  predict(fit_1, data.frame(RT_1st_dim=x_axis))
      Summary_1 <- summary(fit_1)
      Summary <- paste("RSE=", round((Summary_1$sigma/y_axis_diff), 6), sep = "")
      # preparing plot description of linear model
      A <- coef(fit_1)[1]
      A <- paste("a=", signif(A, digits = 3), sep = "")
      B <- coef(fit_1)[2]
      B <- paste("b=", signif(B, digits = 3), sep = "")
      lbl <- paste(x, Summary, sep = " | ")
      lbl <- paste(lbl, A, sep = " | ")
      lbl <- paste(lbl, B, sep = " | ")
      # ploting linear model
      plot(Vec, main=lbl)
      lines(x_axis, Pred, col='red')
      text(Vec, labels=row.names(Vec), pos = 4, cex = 0.5)
      # if linear model fits parametrs green box is added around the plot and clique is added to list fo bad cliques
      if ((Summary_1$sigma/y_axis_diff) < 0.1 & beta.0 < 0 & Summary_1$coefficients[2, 4] < 0.001){
        box(col = "green", lwd = 6)
        Cliques_selected <- append(Cliques_selected,  list(RT))
        list_of_models <- append(list_of_models, list(coef(fit_1)))
      }
      # use non linear model if possible
    }else{
      # fitting nonlinear model
      fit_1 <- nls(RT_2st_dim ~ alpha * exp(beta * RT_1st_dim) + theta, data = Vec, start = Start)
      Pred <-  predict(fit_1, data.frame(RT_1st_dim=x_axis))
      Summary_1 <- summary(fit_1)
      Summary <- paste("RSE=", round(Summary_1$sigma, 6), sep = "")
      lbl <- paste(x, Summary, sep = "-")
      # ploting nonlinear mondel
      A <- coef(fit_1)[1]
      A <- paste("a=", signif(A, digits = 3), sep = "")
      B <- coef(fit_1)[2]
      B <- paste("b=", signif(B, digits = 3), sep = "")
      Th <- coef(fit_1)[3]
      Th <- paste("Th=", signif(Th, digits = 3), sep = "")
      lbl <- paste(x, Summary, sep = " | ")
      lbl <- paste(lbl, A, sep = " | ")
      lbl <- paste(lbl, B, sep = " | ")
      lbl <- paste(lbl, Th, sep = " | ")
      plot(Vec, main=lbl)
      lines(x_axis, Pred, col='red')
      text(Vec, labels=row.names(Vec), pos = 4, cex = 0.5)
      # adding green box if model fits parameters and clique is added to list fo bad cliques
      if (Summary_1$sigma < 0.035 & beta.0 < 0){
        box(col = "green", lwd = 6)
        Cliques_selected <- append(Cliques_selected,  list(RT))
        list_of_models <- append(list_of_models, list(coef(fit_1)))
      }
    }
  }
  dev.off()
  output <- list(Cliques_selected, list_of_models)
  return(output)
}
