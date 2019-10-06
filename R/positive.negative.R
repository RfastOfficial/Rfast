
#[export]
positive.negative <- function(x,method = "min"){
    .Call(Rfast_positive_negative,x,method)
}