#[export]
anova_qpois.reg <- function(mod, poia = NULL) {
   b <- mod$info[, 1]
   vb <- mod$varb
   if (is.null(poia)) {
     poia <- 2:length(b)
     stat <- b[poia] %*% solve(vb[poia, poia], b[poia])
   } else if (length(poia) == 1) {
     stat <- mod$info[poia, 3]
   } else stat <- b[poia] %*% solve(vb[poia, poia], b[poia])
   dof <- length(poia)
   pvalue <- pchisq(stat, dof, lower.tail = FALSE)
   res <- c(stat, pvalue, dof)
   names(res) <- c("statistic", "p-value", "degrees of freedom")
   res
}