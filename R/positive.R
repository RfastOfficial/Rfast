
#[export]
positive <- function(x,method = "min"){
    .Call(Rfast_positive,x,method)
}