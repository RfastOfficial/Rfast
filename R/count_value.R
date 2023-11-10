#[export]
count_value <- function(x,value) {
	.Call(Rfast_count_value,x,value)
}

#[export]
colCountValues<-function(x,values,parallel = FALSE,cores = 0){
	if(parallel){
		.Call(Rfast_col_count_values_p,x,values,cores)
	}else{
		.Call(Rfast_col_count_values,x,values)
	}
}

#[export]
rowCountValues<-function(x,values,parallel = FALSE,cores = 0){
	if(parallel){
		.Call(Rfast_row_count_values_p,x,values,cores)
	}else{
		.Call(Rfast_row_count_values,x,values)
	}
}