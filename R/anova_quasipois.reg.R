#[export]
anova_quasipois.reg <- function(mod0, mod1, n) {
   phi <- mod1$phi   
   df1 <- length(mod1$be) - length(mod0$be) 
   stat <- (mod0$devi - mod1$devi) / df1 / phi
   df2 <- n - length(mod1$be)
   pvalue <- pf(stat, df1, df2, lower.tail = FALSE)
   res <- c(stat, pvalue, df1, df2)
   names(res) <- c("statistic", "p-value", "Numerator dof", "Denominator dof")
   res
}

