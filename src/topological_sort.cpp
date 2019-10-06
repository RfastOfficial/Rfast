
#include <vector>
#include <stack>
#include <RcppArmadillo.h>
//#include "mn.h"

using namespace Rcpp;

using std::count;
using std::stack;
using std::vector;

static vector<int> mywhich(IntegerVector& x,int c){
  vector<int> f(c);
  for(int i=0,k=0;i<x.size();++i)
    if(x[i]>0)
      f[k++]=i;
    return f;
}

//[[Rcpp::export]]
IntegerVector topological_sort(IntegerMatrix dag) {
  const int n = dag.nrow();
  IntegerVector indeg = colSums(dag(Range(0,n-1),Range(0,n-1)));
  stack<int> zero_indeg; // a stack of nodes with no parents
  vector<int> cs;
  int i = 0;
  IntegerVector ord(n),tmp(dag.ncol());
  int v,c;
  double m;
  unsigned int j;
  for(int i=0;i<indeg.size();++i)
    if(indeg[i]==0)
      zero_indeg.push(i);
  for(i=0;zero_indeg.size()>0 && i < n; ++i) {
    v = zero_indeg.top();  // get top element
    zero_indeg.pop(); //  pop top element
    ord[i] = v;
    tmp=dag.row(v);
    c=count(tmp.begin(),tmp.end(),1);
    if(c>0){
      cs = mywhich(tmp,c); //  children(A, v) 
      for (m=j=0;j<cs.size();++j) {
        m = cs[j];
        indeg[m]--;
        if ( indeg[m] == 0 ) {
          zero_indeg.push(m);   // push m 
        }
      }
    }
  }  
  return ord;
}


static IntegerVector topological_sort2(IntegerMatrix dag) {
    const int n = dag.nrow();
    IntegerVector indeg = colSums(dag(Range(0,n-1),Range(0,n-1)));
    vector<int> cs,zero_indeg; // a stack of nodes with no parents
    int i = 0;
    IntegerVector ord(n),tmp(dag.ncol());
    int v,c;
    double m;
    unsigned int j;
    for(int i=0;i<indeg.size();++i){
        if(indeg[i]==0){
            zero_indeg.emplace_back(i);
        }
    }
    for(i=0;zero_indeg.size()>0 && i < n; ++i) {
        v = zero_indeg.back();  // get top element
        zero_indeg.pop_back();
        ord[i] = v;
        tmp=dag.row(v);
        c=count(tmp.begin(),tmp.end(),1);
        if(c>0){
            cs = mywhich(tmp,c); //  children(A, v) 
            for (m=j=0;j<cs.size();++j) {
                m = cs[j];
                indeg[m]--;
                if ( indeg[m] == 0 ) {
                    zero_indeg.emplace_back(m);   // push m 
                }
            }
        }
    }  
    return ord;
}






RcppExport SEXP Rfast_topological_sort(SEXP dagSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< IntegerMatrix >::type dag(dagSEXP);
    __result = topological_sort(dag);
    return __result;
END_RCPP
}
