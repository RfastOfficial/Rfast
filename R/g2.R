#[export]
g2tests_perm <- function(data,x,y,dc,nperm) {
    .Call(Rfast_g2tests_perm,data,x,y,dc,nperm)
}

#[export]
g2Test_univariate <- function(data,dc) {
    .Call(Rfast_g2Test_univariate,data,dc)
}

#[export]
g2Test_perm <- function(data,x,y,cs,dc,nperm) {
    .Call(Rfast_g2Test_perm,data,x,y,cs,dc,nperm)
}
#[export]
g2tests <- function(data,x,y,dc) {
    .Call(Rfast_g2tests,data,x,y,dc)
}
#[export]
g2Test_univariate_perm <- function(data,dc,nperm) {
    .Call(Rfast_g2Test_univariate_perm,data,dc,nperm)
}
#[export]
g2Test <- function(data,x,y,cs,dc,parallel = FALSE) {
    .Call(Rfast_g2Test,data,x,y,cs,dc,parallel)
}