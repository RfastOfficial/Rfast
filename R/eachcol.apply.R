#[export]
eachcol.apply<-function(x,y,indices = NULL,oper = "*",apply = "sum", parallel = FALSE, cores = 0){
	.Call(Rfast_eachcol_apply,x,y,indices,oper,apply, parallel)
}
