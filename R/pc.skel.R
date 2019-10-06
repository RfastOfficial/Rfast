#[export]
pc.skel <- function(dataset, method = "pearson", alpha = 0.01, R = 1, stat = NULL, ini.pvalue = NULL) {
  is.init.vals <- c(!is.null(stat), !is.null(ini.pvalue))
  if (is.null(stat)) stat <- matrix();
  if (is.null(ini.pvalue)) ini.pvalue <- matrix();
  .Call(Rfast_pc_skel,dataset, method, alpha, R, stat, ini.pvalue, is.init.vals)
}
