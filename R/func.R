 #' @title pMCMC Function
 #' @param x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
 #' @param null A numeric value decsribing what the null hypothesis should be
 #' @param twotail Whether to conduct a one-tailed hypothesis or a two-tailed hypotheses. Default = true indicating a two-tailed test will be done.
pmcmc <- function(x, null = 0, twotail = TRUE){
  if(twotail){
    2*(1 - max(table(x<=null) / length(x)))
  } else{
    (1 - max(table(x<=null) / length(x)))
  }
}


#' @title overall_repeatability
#' @description Will take the posterior distributions from a multi-response model and calculate the posterior distribution of repeatability
#' @param id_sd These are the SD estimates for ID
#' @param clutch_sd These are the SD estimates for clutch ID
#' @param sigma_sd These are the SD estimates for within individual / residual 
#' @param trait A character string with the name of the trait in the data set
#' 
#' 
id_sd = id_deli
clutch_sd = clutch_deli
sigma_sd = sigma_deli
repeatability <- function(id_sd, clutch_sd, sigma_sd, trait) {
			    id_v <- id_sd[,grep(trait, colnames(id_sd))]^2
			clutch_v <- clutch_sd[,grep(trait, colnames(clutch_sd))]^2
			 sigma_v <- sigma_sd[,grep(trait, colnames(sigma_sd))]^2
    
	# Calculate R
	    R  <- id_v / (id_v + clutch_v + sigma_v)

    return(data.frame(R = mean(R),
					  `L 95% CI` = quantile(R, 0.025)[1],
					  `U 95% CI` = quantile(R, 0.975)[1], check.names = FALSE, row.names = NULL))
}
