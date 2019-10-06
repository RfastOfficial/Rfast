#[export]
apply.condition <- function(x,method = "+",oper = ">",cond.val = 0) {
  .Call(Rfast_apply_condition,x,method,oper,cond.val)
}