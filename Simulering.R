
# SIMULERINGSSTUDIEN


set.seed(1337)


# Inställningar

R <- 1000
n_train <- 500
n_cal   <- 500
n_test  <- 1000

gamma_values <- c(0, 0.5, 1, 2)
delta_values <- c(0, 0.5, 1, 2)


# Weighted conformal-kvantil


Viktad_q <- function(scores, Kalk_vikt, test_vikt, alpha) {
  
  ord <- order(scores)
  scores <- scores[ord]
  Kalk_vikt <- Kalk_vikt[ord]
  
  cum_w <- cumsum(Kalk_vikt) /
    (sum(Kalk_vikt) + test_vikt)
  
  i <- which(cum_w >= 1 - alpha)[1]
  
  if (is.na(i)) Inf else scores[i]
}



# Generell Monte Carlo-iteration


en_sim <- function(gamma = 0,
                          delta = 0,
                          alpha = 0.10,
                          weighted = TRUE) {
  
  # Data
  train <- data.frame(x = rnorm(n_train))
  cal   <- data.frame(x = rnorm(n_cal))
  test  <- data.frame(x = rnorm(n_test, mean = delta))
  
  train$y <- 2*train$x + gamma*train$x^2 + rnorm(n_train)
  cal$y   <- 2*cal$x   + gamma*cal$x^2   + rnorm(n_cal)
  test$y  <- 2*test$x  + gamma*test$x^2  + rnorm(n_test)
  
  
  # LM
  lin <- lm(y ~ x, data = train)
  
  cal_pred  <- predict(lin, cal)
  test_pred <- predict(lin, test)
  
  scores <- abs(cal$y - cal_pred)
  
  
  
  # Split CP
  
  k <- ceiling((n_cal + 1) * (1 - alpha))
  q_split <- sort(scores)[k]
  
  split_lower <- test_pred - q_split
  split_upper <- test_pred + q_split
  
  
  
  # LM
 
  
  lm_int <- predict(
    lin,
    test,
    interval = "prediction",
    level = 1 - alpha
  )
  
  
 
  # Resultat Grund

  
  Resultat <- data.frame(
    gamma = gamma,
    delta = delta,
    alpha = alpha,
    
    Split_covered = mean(
      test$y >= split_lower &
        test$y <= split_upper
    ),
    
    Split_bredd = median(
      split_upper - split_lower
    ),
    
    LM_covered = mean(
      test$y >= lm_int[, "lwr"] &
        test$y <= lm_int[, "upr"]
    ),
    
    LM_bredd = median(
      lm_int[, "upr"] -
        lm_int[, "lwr"]
    )
  )
  
  

  # viktad CP
 
  
  if (weighted) {
    
    Kalk_vikt <- exp(
      delta * cal$x - delta^2 / 2
    )
    
    test_vikts <- exp(
      delta * test$x - delta^2 / 2
    )
    
    q_viktad <- sapply(
      test_vikts,
      \(w) Viktad_q(
        scores,
        Kalk_vikt,
        w,
        alpha
      )
    )
    
    Viktad_bredd <- 2 * q_viktad
    
    Resultat$Weighted_covered <- mean(
      test$y >= test_pred - q_viktad &
        test$y <= test_pred + q_viktad
    )
    
    Resultat$Viktad_bredd <- median(
      Viktad_bredd[is.finite(Viktad_bredd)]
    )
    
    Resultat$Weighted_infinite_rate <- mean(
      is.infinite(Viktad_bredd)
    )
  }
  
  Resultat
}



# Funktion hela experiment


Simulering_n <- function(gamma,
                           delta,
                           alpha = 0.10,
                           R = 1000,
                           weighted = TRUE,
                           paired = FALSE) {
  
  # Experimentdesign
  if (paired) {
    
    # Endast parade gamma/delta
    design <- expand.grid(
      scenario = seq_along(gamma),
      alpha = alpha
    )
    
    design$gamma <- gamma[design$scenario]
    design$delta <- delta[design$scenario]
    design$scenario <- NULL
    
  } else {
    
    # Alla komb
    design <- expand.grid(
      gamma = gamma,
      delta = delta,
      alpha = alpha
    )
  }
  
  
  # Kör sim
  resultaten <- do.call(
    rbind,
    lapply(1:nrow(design), function(i) {
      
      d <- design[i, ]
      
      temp <- do.call(
        rbind,
        lapply(1:R, function(r) {
          
          ut <- en_sim(
            gamma = d$gamma,
            delta = d$delta,
            alpha = d$alpha,
            weighted = weighted
          )
          
          ut$rep <- r
          ut
        })
      )
      
      temp
    })
  )
  
  resultaten
}



# HUVUDSIM


huvud_resultaten <- Simulering_n(
  gamma = gamma_values,
  delta = delta_values,
  alpha = 0.10,
  R = R
)


# MONTE CARLO STANDARDFEL

mcse_mal <- 0.005

# MCSE för scenario och metod (alla)
mcse_resultaten <- aggregate(
  cbind(
    Split_covered,
    Weighted_covered,
    LM_covered
  ) ~ gamma + delta,
  data = huvud_resultaten,
  FUN = function(x) sd(x) / sqrt(length(x))
)

# Största 
max_mcse <- max(
  mcse_resultaten[
    c(
      "Split_covered",
      "Weighted_covered",
      "LM_covered"
    )
  ],
  na.rm = TRUE
)

# Antal faktiska replikationer
R_nu <- length(unique(huvud_resultaten$rep))

# Uppskattat R som krävs för MSCE 0.005
R_behov <- ceiling(
  R_nu * (max_mcse / mcse_mal)^2
)

# Resultat
cat(
  "Nuvarande R:", R_nu,
  "\nStörsta MCSE:", round(max_mcse, 4),
  "\nMål för MCSE:", mcse_mal,
  "\n"
)

if (max_mcse <= mcse_mal) {
  
  cat(
    "x R =", R_nu,
    "är tillräckligt för det valda MCSE-målet.\n"
  )
  
} else {
  
  cat(
    "✗ R =", R_nu,
    "är inte tillräckligt.\n",
    "Rekommenderat R är minst cirka",
    R_behov, "replikationer.\n"
  )
}


# SENSITIVITET


sensitivity_resultaten <- Simulering_n(
  gamma = c(0, 1, 2),
  delta = c(0, 1, 2),
  alpha = c(0.20, 0.10, 0.05),
  R = R,
  paired = TRUE
)

saveRDS(huvud_resultaten, "huvud_resultaten.rds")
saveRDS(sensitivity_resultaten, "sensitivity_resultaten.rds")
