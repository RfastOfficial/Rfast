#[export]
colprods<- function(x,method = "direct"){
  .Call(Rfast_col_prods,x,method)
}

#[export]
rowprods<- function(x){
  .Call(Rfast_row_prods,x)
}