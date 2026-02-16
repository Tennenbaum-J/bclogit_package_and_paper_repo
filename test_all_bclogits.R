pacman::p_load(bclogit, data.table)
prior_types = c("Naive", "G prior", "PMP", "Hybrid")
concordant_methods = c("GLM", "GEE", "GLMM")
data("fhs")
Nmax = 500
fhs = fhs[1 : Nmax, ]

res = data.table(expand.grid(prior_types, concordant_methods))
colnames(res) = c("prior", "method")
res[, lower := NA_real_]
res[, upper := NA_real_]

level = 0.95
for (prior_type in prior_types){
  for (concordant_method in concordant_methods){
    cat("\n\n ===== bclogit prior_type: ", prior_type, " with method: ", concordant_method, "\n")
    fit <- bclogit(
      formula = PREVHYP ~ TOTCHOL + CIGPDAY + BMI + HEARTRTE, 
      data = fhs, 
      treatment = PERIOD, 
      strata = RANDID,
      prior_type = prior_type,
      concordant_method = concordant_method
    )
    for (type in c("HPD_one", "HPD_many", "CR")){
      print(summary(fit, inference_method = type, level = level))
      print(confint(fit, type = type, level = level))
      if (type == "HPD_one"){
        conf_tab = confint(fit, type = type, level = level)
        res[prior == prior_type & method == concordant_method, lower := conf_tab[1, 1]]
        res[prior == prior_type & method == concordant_method, upper := conf_tab[1, 2]]
      }
    }
  }  
}

res
