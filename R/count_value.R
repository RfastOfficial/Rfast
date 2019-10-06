#[export]
count_value <- function(x,value) {
	.Call(Rfast_count_value,x,value)
}

#[export]
colCountValues<-function(x,values,parallel = FALSE){
	if(parallel){
		.Call(Rfast_col_count_values_p,x,values)
	}else{
		.Call(Rfast_col_count_values,x,values)
	}
}

#[export]
rowCountValues<-function(x,values,parallel = FALSE){
	if(parallel){
		.Call(Rfast_row_count_values_p,x,values)
	}else{
		.Call(Rfast_row_count_values,x,values)
	}
}