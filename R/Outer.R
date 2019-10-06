#[export]
Outer<-function(x,y,oper="*"){
	if(identical(oper,"%%")) oper<-"%"
	.Call(Rfast_Outer,x,y,oper)
}