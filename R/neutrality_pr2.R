#' Creating a PR2 plot and neutrality plot based on the table with the get_RSCU function
#'
#' @param get_RSCU_out the table was created using the get_RSCU function
#' @param select the "PR2_plot" or "neutrality_plot" output
#'
#' @return A plot object.
#' @export
#'

neutrality_pr2 <- function(get_RSCU_out=get_RSCU_out,select=""){
  
  if (base::missing(get_RSCU_out)) {
    stop("The get_RSCU_out predictions are required. Please provide a valid argument.",
         call. = FALSE)
  }
  
  get_RSCU_out$A1 <- NA
  get_RSCU_out$A2 <- NA
  get_RSCU_out$A3 <- NA
  get_RSCU_out$C1 <- NA
  get_RSCU_out$C2 <- NA
  get_RSCU_out$C3 <- NA
  get_RSCU_out$T1 <- NA
  get_RSCU_out$T2 <- NA
  get_RSCU_out$T3 <- NA
  get_RSCU_out$G1 <- NA
  get_RSCU_out$G2 <- NA
  get_RSCU_out$G3 <- NA
  
  
  get_RSCU_out <- get_RSCU_out %>%
    mutate(
      A1 = ifelse(str_sub(codon, 1, 1) == "a", eff, 0),
      A2 = ifelse(str_sub(codon, 2, 2) == "a", eff, 0),
      A3 = ifelse(str_sub(codon, 3, 3) == "a", eff, 0),
      T1 = ifelse(str_sub(codon, 1, 1) == "t", eff, 0),
      T2 = ifelse(str_sub(codon, 2, 2) == "t", eff, 0),
      T3 = ifelse(str_sub(codon, 3, 3) == "t", eff, 0),
      C1 = ifelse(str_sub(codon, 1, 1) == "c", eff, 0),
      C2 = ifelse(str_sub(codon, 2, 2) == "c", eff, 0),
      C3 = ifelse(str_sub(codon, 3, 3) == "c", eff, 0),
      G1 = ifelse(str_sub(codon, 1, 1) == "g", eff, 0),
      G2 = ifelse(str_sub(codon, 2, 2) == "g", eff, 0),
      G3 = ifelse(str_sub(codon, 3, 3) == "g", eff, 0),
    ) %>% 
    group_by(Species) %>% summarise(sum_A1=sum(A1),
                                    sum_A2=sum(A2),
                                    sum_A3=sum(A3),
                                    sum_C1=sum(C1),
                                    sum_C2=sum(C2),
                                    sum_C3=sum(C3),
                                    sum_T1=sum(T1),
                                    sum_T2=sum(T2),
                                    sum_T3=sum(T3),
                                    sum_G1=sum(G1),
                                    sum_G2=sum(G2),
                                    sum_G3=sum(G3)) %>%
    mutate(
      sum3ATGC = sum_A3 + sum_T3 + sum_G3 + sum_C3,
      GC3 = ifelse(sum3ATGC == 0, 0, (sum_G3 + sum_C3) / sum3ATGC),
      sum1ATGC = sum_A1 + sum_T1 + sum_G1 + sum_C1,
      GC1 = ifelse(sum1ATGC == 0, 0, (sum_G1 + sum_C1) / sum1ATGC),
      sum2ATGC = sum_A2 + sum_T2 + sum_G2 + sum_C2,
      GC2 = ifelse(sum2ATGC == 0, 0, (sum_G2 + sum_C2) / sum2ATGC),
      GC12 = (GC1 + GC2) / 2,
      sum3AT = sum_A3 + sum_T3,
      A3_AT3 = ifelse(sum3AT == 0, 0, sum_A3 / sum3AT),
      sum3GC = sum_G3 + sum_C3,
      G3_GC3 = ifelse(sum3GC == 0, 0, sum_G3 / sum3GC)
      ) %>% select(-c(sum3ATGC,sum1ATGC,sum2ATGC,sum3AT,sum3GC))
  
  if (select=="neutrality_plot"){
  
  model <- lm(GC12 ~ GC3, data = get_RSCU_out)
  r_squared <- summary(model)$r.squared
  equation <- paste("y =", round(coef(model)[2], 2), "x +", round(coef(model)[1], 2))
  r_squared_label <- paste("R² =", round(r_squared, 3))
  
  neutrality_plot <- get_RSCU_out %>%
    ggplot(aes(x = GC3, y = GC12)) +
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "GC3", y = "GC12", title = paste0(r_squared_label,", ",equation)) +
    theme_bw() 

  return(list(table = get_RSCU_out, plot = neutrality_plot))
  }
  
  if (select=="PR2_plot"){
  model <- lm(A3_AT3 ~ G3_GC3, data = get_RSCU_out) 
  r_squared <- summary(model)$r.squared
  equation <- paste("y =", round(coef(model)[2], 2), "x +", round(coef(model)[1], 2))
  r_squared_label <- paste("R² =", round(r_squared, 3))
  
  PR2_plot <- get_RSCU_out %>%
    ggplot(aes(x = G3_GC3, y = A3_AT3)) +
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) + 
    labs(x = "G3/(G3+C3)", y = "A3/(A3+T3)", title = paste0(r_squared_label,", ",equation)) +
    theme_bw()
    
  return(list(table = get_RSCU_out, plot = PR2_plot))
  }
  base::message(base::paste0("Success"))
}