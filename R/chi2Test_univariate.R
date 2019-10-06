#[export]
chi2Test_univariate <- function(data,dc) {
    .Call(Rfast_chi2Test_univariate,data,dc)
}