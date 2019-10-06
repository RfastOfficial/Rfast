#[export]
chi2tests <- function(data,x,y,dc) {
    .Call(Rfast_chi2tests,data,x,y,dc)
}