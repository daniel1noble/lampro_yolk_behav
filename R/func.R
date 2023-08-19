#' @title brms_model_check
#' @description Checks a 'brms' model by plotting a histogram of residuals and a scatterplot of the observed and predcited log mass
#' @param model The 'brms' model object that models log (mass (grams))
#' @param main Title of the plot
#' @param xlab Label of the x-axis which defaults to 'Residuals'
#' @return Returns z-transformed age (in SD units and cenetred on 0)
brms_model_check <- function(model, main = NULL, xlab = "Residuals"){
  
  # Histogram of residuals - assumed normal - pretty good to me
      resid <-  model$data$lnMass - predict(model, summary = TRUE)[,"Estimate"]
 
  # Look at the residuals - should be normally distributed
     hist(resid, main = main, xlab = xlab)
  
  # We already know roughly from R2 that model does good job predciting observed response, but lets have a look. Little bit of over/underpredciting but nothing serious
      plot(model$data$lnMass ~ predict(model, summary = TRUE)[,"Estimate"], ylab = "Observed lnMass", xlab = "Predicted lnMass", main = main)
      abline(0,1, col = "red")
}

#' @title extract_post
#' @param model The 'brms' model object
#' @param trait A character string with the name of the trait in the data set
#' @return A data frame with the posterior distribution of the estimated means for each level of the fixed effects
 extract_post <- function(model, trait){
		post <- posterior_samples(model, pars = paste0("^b_", trait))

		# Calculate means for each level from fixed effects
				       intercept <- paste0("b_", trait, "_Intercept")
				egg_treatcontrol <- paste0("b_", trait, "_egg_treatcontrol")
				         temphot <- paste0("b_", trait, "_temphot")
	    egg_treatcontrol_temphot <- paste0("b_", trait, "_temphot:egg_treatcontrol")

           A_23 <- post[,intercept]
           A_28 <- post[,intercept] + post[,temphot]
           C_23 <- post[,intercept] + post[,egg_treatcontrol]
           C_28 <- post[,intercept] + post[,temphot] + post[,egg_treatcontrol] + post[,egg_treatcontrol_temphot]

		   return(data.frame(A23 = A_23,
							 A28 = A_28,
							 C23 = C_23,
							 C28 = C_28))
}

#' @title contrast_post_main
#' @param posterior_data The dataframe that contains the posterior distribution of the estimated means for each level of the fixed effects
#' @return A data frame with the contrast along with 95% CI, pmcmc and contrast name for the main effects (pooled posteriors). 
contrast_post_main <- function(posterior_data){

	      temp <- c(posterior_data[,"C23"], posterior_data[,"A23"]) - (c(posterior_data[,"C28"], posterior_data[,"A28"]))
		invest <- c(posterior_data[,"C23"], posterior_data[,"C28"]) - (c(posterior_data[,"A23"], posterior_data[,"A28"]))
         
	   table <- data.frame(rbind(temp  %>% posterior_summary(), invest %>% posterior_summary()))
	   table <- table %>% mutate(pmcmc = c(pmcmc(temp), pmcmc(invest)),
	   							contrast  = c("Temp (23-28)", "Invest (C-A)"))  %>% select(contrast, everything())
	return(table)
}

#' @title contrast_post
#' @param posterior_data The dataframe that contains the posterior distribution of the estimated means for each level of the fixed effects
#' @return A data frame with the contrast along with 95% CI, pmcmc and contrast name. 
contrast_post <- function(posterior_data){

	    cold <- (posterior_data[,"C23"] - posterior_data[,"A23"])  
         hot <- (posterior_data[,"C28"] - posterior_data[,"A28"])
	   inter <- hot-cold

	   table <- data.frame(rbind(cold  %>% posterior_summary(), hot %>% posterior_summary(), inter %>% posterior_summary()))
	   table <- table %>% mutate(pmcmc = c(pmcmc(cold), pmcmc(hot), pmcmc(inter)),
	   							contrast  = c("C23 - A23", "C28 - A28", "(C23 - A23) - (C28 - A28) (interaction)"))  %>% select(contrast, everything())
	return(table)
}

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
