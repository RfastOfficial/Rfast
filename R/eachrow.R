#[export]
eachrow <- function(x,y,oper = "*",method = NULL){
	.Call(Rfast_eachrow,x,y,if(oper=="==") "=" else oper,method)
}