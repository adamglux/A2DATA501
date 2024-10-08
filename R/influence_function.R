#' Calculate and Plot Measures of Influence for Linear Models
#'
#' @description
#' This function calculates and evaluates key influence measures, including Cook's Distance, DFFITS, and Hadi's Influence, for linear or generalized linear models by computing diagnostic metrics to identify influential observations. The function has an optional setting that will generate plots to visualize these influence measures, aiding in the detection of data points that may impact the model's fit. This tool is helpful for model assessment and ensuring the reliability of statistical inferences.
#'
#'
#' @param model a glm or lm object
#' @param data a dataframe
#' @param plot optional parameter to include plots in the output, default is "none", options include "none", "all", "Cooks", "DFFITS", "Hadi"
#'
#' @return Returns a dataframe with all influential measures. If `plot` parameter is set to anything other than "none", the function will return a list with dataframe and `ggplot2` plot objects.
#'
#' @import ggplot2
#' @importFrom stats model.matrix
#' @export
#'
#' @examples
#' #import dataset mtcars
#' data1 <- mtcars
#' model1 <- lm(mpg ~ wt + hp, data = data1)
#'
#' ## run function to create influence measures + plots
#' sample1 <- influenceR(model1, data1, "all")
#'
#' ## show dataframe of measures
#' sample1$influence_measures
#'
#' ## show plots
#' sample1$plots
#'
#' ## show individual plots
#' sample1$plots$Cooks
#' sample1$plots$DFFITS
#' sample1$plots$Hadi
#'
#' ## run function to create dataframe only (no plots)
#' sample2 <- influenceR(model1, data1)
#' sample2
#'
#' ## view / manipulate individual influence measure
#' sample2$DFFITs


influenceR <- function(model, data, plot = "none") {

  ##### input validation
  if (!is.data.frame(data)) {
    stop("Input data must be a data.frame.")
  }

  if (!inherits(model, c("lm", "glm"))) {
    stop("Input model must be an lm or glm object.")
  }

  if (nrow(data) < 5) {
    stop("Dataset is too small")
  }

  if (ncol(data) < 2) {
    stop("Input data must have at least 2 columns.")
  }

  if (nrow(data) < ncol(data)) {
    stop("Dataset is the wrong dimension")
  }

  if (sum(is.na(data)) / prod(dim(data)) * 100 > 5) {
    stop("Input data contains more than 5% NA values.")
  }

  if (any(is.infinite(as.matrix(data)))) {
    stop("Input data contains infinite values.")
  }

  if (!plot %in% c("none", "all", "cooks", "DFFITS", "Hadi") ) {
    stop("Input for plot must be either: none, all, cooks, DFFITS, or Hadi")
  }

  ##### obtain model variables
  # get y
  y <- model$model[1]
  y <- data.matrix(y)
  # get X matrix
  X <- model.matrix(model, data = data)

  # compute Beta hat
  XtX <- t(X) %*% X
  XtX_inv <- solve(XtX)
  beta_hat <- XtX_inv %*% t(X) %*% y

  #get y-hat (fitted values)
  y_hat <- X %*% beta_hat
  #get residuals
  residuals <- y - y_hat

  #get them Hat matrix (and diagonals)
  H <- X %*% XtX_inv %*% t(X)
  leverage <- diag(H)

  #identity matrix


  #get studentised residuals
  n <- nrow(X)
  p <- ncol(X)
  Id <- diag(rep(1, n))
  SSE <- as.numeric(t(y) %*% (Id - H) %*% y)
  RSE <- sqrt(sum(residuals^2) / (n - p))
  MSE <- sum(residuals^2) / (n - p)
  studentized_residuals <- residuals / (RSE * sqrt(1 - leverage))
  normalised_residuals <- residuals / sqrt(SSE)
  high_residuals <- studentized_residuals[abs(studentized_residuals)>= 2]

  ##### distance measures

  # Cook's D
  cooks_distance <- (residuals^2 / (p * MSE)) * (leverage / (1 - leverage)^2)

  # DFFITs
  dffits <- studentized_residuals * sqrt(leverage) / (1 - leverage)

  # Hadi's
  H1 <- leverage / (1 - leverage)
  H2 <- p+1 / (1 - leverage) #(residuals^2) / (MSE * (1 - leverage))
  H3 <- normalised_residuals^2 / (1 - normalised_residuals^2)
  Hadi_influence <- H1 + (H2 * H3)

  ##### put it all into a dataframe
  influencers <- data.frame(cooks_distance,dffits,Hadi_influence)
  names(influencers) <- c("Cook's D", "DFFITs", "Hadi's Influence")

  #### Plotting
  # options: "none", "all", "Cooks", "DFFITS", "Hadi"

  #cooks threshold line:
  cooks_cutoff <- 4/nrow(data)
  dffits_cutoff <- 2*sqrt((p)/(n-p-1))
  hadi_cutoff <- 2

  # labels
  cooks_label <- rep(NA,length(cooks_distance))
  dffits_label <- rep(NA,length(dffits))
  hadi_label <- rep(NA,length(Hadi_influence))
  for (i in 1:nrow(data)){
    cooks_label[i] <- ifelse(abs(cooks_distance[i]) >= cooks_cutoff, i, "")
    dffits_label[i] <- ifelse(abs(dffits[i]) > dffits_cutoff, i, "")
    hadi_label[i] <- ifelse(abs(Hadi_influence[i]) > hadi_cutoff, i, "")
  }



  # empty list of functions
  plots <- list()

  cook_guides <- data.frame(1:length(cooks_distance),
                            cooks_distance, cooks_label)
  dffits_guides <- data.frame(1:length(dffits),
                              dffits, dffits_label)
  hadi_guides <- data.frame(1:length(Hadi_influence),
                            Hadi_influence, hadi_label)


  if (plot == "cooks" | plot == "all") {

    plots$Cooks <- ggplot(cook_guides, aes(x = cook_guides[,1], y = cooks_distance)) +
      geom_segment(aes(xend = cook_guides[,1], yend = 0)) +
      geom_hline(yintercept = cooks_cutoff, linetype = "dashed", color = "darkred") +
      geom_text(aes(label = cook_guides[,3]), vjust = 0.75, check_overlap = T,
                family = "Courier", nudge_x = 0.75) +
      labs(title = "Cook's Distance", y = "Cook's Distance", x = "Observation")
  }

  if (plot == "DFFITS" | plot == "all") {
    plots$DFFITS <- ggplot(dffits_guides, aes(x = dffits_guides[,1], y = dffits)) +
      geom_segment(aes(xend = dffits_guides[,1], yend = 0)) +
      geom_hline(yintercept = dffits_cutoff, linetype = "dashed", color = "darkred") +
      geom_hline(yintercept = -dffits_cutoff, linetype = "dashed", color = "darkred") +
      geom_text(aes(label = dffits_guides[,3]),  vjust = 0.75, check_overlap = T,
                family = "Courier", nudge_x = 0.75) +
      labs(title = "DFFITS", y = "DFFITS", x = "Observation")
  }

  if (plot == "Hadi" | plot == "all") {
    plots$Hadi <- ggplot(hadi_guides, aes(x = hadi_guides[,1], y = Hadi_influence)) +
      geom_segment(aes(xend = hadi_guides[,1], yend = 0)) +
      geom_hline(yintercept = hadi_cutoff, linetype = "dashed", color = "darkred") +
      geom_text(aes(label = hadi_guides[,3]),  vjust = 0.75, check_overlap = T,
                family = "Courier", nudge_x = 0.75) +
      labs(title = "Hadi's Influence", y = "Hadi's Influence", x = "Observation")
  }

  #### return
  if (plot == "none") {
    return(influencers)
  }

  else {
    return(list(influence_measures = influencers, plots = plots))
    }

}
