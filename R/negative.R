#[export]
negative <- function(x,method = "min"){
    .Call(Rfast_negative,x,method)
}