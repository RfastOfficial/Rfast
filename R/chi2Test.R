#[export]
chi2Test <- function(data,x,y,cs,dc) {
    .Call(Rfast_chi2Test,data,x,y,cs,dc)
}


#[export]
chi2Test_univariate <- function(data,dc) {
    .Call(Rfast_chi2Test_univariate,data,dc)
}


#[export]
chi2tests <- function(data,x,y,dc) {
    .Call(Rfast_chi2tests,data,x,y,dc)
}