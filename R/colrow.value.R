#[export]
colrow.value <- function(x,value=0) {
	.Call(Rfast_col_row_value,x,value)
}