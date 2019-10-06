#[export]
chi2Test <- function(data,x,y,cs,dc) {
    .Call(Rfast_chi2Test,data,x,y,cs,dc)
}