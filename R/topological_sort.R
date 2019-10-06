#[export]
topological_sort <- function (dag) {
  a <- .Call(Rfast_topological_sort, dag) + 1
  if ( sort_unique.length(a) != dim(dag)[2] ) a <- NA
  a
}





# topological_sort <- function(dag) {
#    ## TOPOLOGICAL_SORT Return the nodes in topological order (parents before children).
#    ## order = topological_sort(adj_mat)
#    ## Base Code is from BNT
#   n = dim(dag)[1];
#   indeg = numeric(n)
#   zero_indeg = NULL    ## a stack of nodes with no parents
#     ##  for ( i in 1:n ) {
#     ##    indeg[i] = length( which(dag[, i] > 0) )  ## parents(A, i)
#     ##    if ( indeg[i] == 0 )   zero_indeg = c(i, zero_indeg)
#     ##  }
#   indeg = Rfast::colsums(dag)  ## parents(A, i)
#   for ( i in 1:n )   if ( indeg[i] == 0 )   zero_indeg = c(i, zero_indeg)
#   i = 1
#   ord = numeric(n)
#   while ( !is.null(zero_indeg) & i <= n ) {
#     v = zero_indeg[1]   ##  pop v
#     zero_indeg = zero_indeg[-1]
#     ord[i] = v
#     i = i + 1;
#     cs = which(dag[v, ] > 0)  ##  children(A, v)  
#     if ( length(cs) > 0 ) {
#       for ( j in 1:length(cs) ) {
#         m = cs[j]
#         indeg[m] = indeg[m] - 1
#         if ( indeg[m] == 0 ) {
#           zero_indeg = c(m, zero_indeg)   ## push m 
#         }
#       }
# 
#     }  ## end  if ( length(cs) > 0 )
# 
#   }  ## end while ( !is.null(zero_indeg) & i <= n )
#   
#   ord
# }
