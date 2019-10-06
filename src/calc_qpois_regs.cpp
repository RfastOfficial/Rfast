// Author: Manos Papadakis

#include "calc_qpois_regs.h"

arma::vec calc_qpois_regs(arma::mat& x,arma::vec& y,const double tol,const double ylogy,const double my) {
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=2;
  unsigned int i;
  colvec b_old(d),b_new(d),L1(d),yhat(n);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones);
  vec m(n),z_col_1(n);
  arma::colvec f(pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(my);
  for(i=0;i<pcols;++i){
    b_old(0)=lgmeany;
    b_old(1)=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    for(dif=1.0;dif > 0.000000001 && dif > tol;){
      sm=szm=sz2m=0.0;
      yhat=z*b_old;
      m=exp(yhat);
      L1=z_tr*(y-m);
      sm=sum(m);
      szm=sum(m%z_col_1);
      sz2m=sum(m%square(z_col_1));
      t=1.0/(sm*sz2m-szm*szm);
      inv_L2.at(0,0)=sz2m*t;
      inv_L2.at(0,1)=inv_L2.at(1,0)=-szm*t;
      inv_L2.at(1,1)=sm*t;
      b_new=b_old+inv_L2*L1;
      dif=sum(abs(b_new-b_old));
      b_old=b_new;
    }
  	const double phi=arma::sum(arma::pow(y-m,2)/m)/(n-d); 
	f[i]=(b_new[1]*b_new[1])/(phi*inv_L2.at(1,1));
  }
  return f;
}
