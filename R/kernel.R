
#[export]
kernel <- function(x, h) {
    .Call(Rfast_kernel,t(x),h)
}