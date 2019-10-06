#[export]
eachcol.apply<-function(x,y,indices = NULL,oper = "*",apply = "sum"){
	.Call(Rfast_eachcol_apply,x,y,indices,oper,apply)
}