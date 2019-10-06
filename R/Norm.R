#[export]
Norm<- function(x, type = "F") {
  .Call(Rfast_Norm,x,type)
}