library(Rfast)

group <- c(2, 1, 3, 1, 2, 3, 2, 3, 2, 3)
x <- cbind(c(5, 2, 3, 4, 8, 4, 10, 15, 30, 38),
           c(1, 2, 3, 4, 5, 13, 6, 7, 8, 22))

reference <- function(x, group, center) {
  group <- factor(group)
  t(apply(x, 2L, function(column) {
    deviations <- abs(column - ave(column, group, FUN = center))
    fit <- anova(lm(deviations ~ group))
    c(fit[1L, "F value"], fit[1L, "Pr(>F)"])
  }))
}

for (type in c("levene", "bf")) {
  center <- if (type == "levene") mean else median
  expected <- reference(x, group, center)
  actual <- vartests(x, group, type = type)
  stopifnot(isTRUE(all.equal(unname(actual), expected)))

  single <- vartests(x[, 1L, drop = FALSE], factor(group), type = type,
                     logged = TRUE)
  stopifnot(isTRUE(all.equal(unname(single[, 1L]), expected[1L, 1L])),
            isTRUE(all.equal(unname(single[, 2L]), log(expected[1L, 2L]))))
}
