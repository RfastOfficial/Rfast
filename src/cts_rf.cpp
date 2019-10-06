#include "cts_rf.h"

#define DEBUG 0
#define db_print(...) \
	do { if (DEBUG) Rprintf(__VA_ARGS__); } while (0)

// Author: Manos Papadakis

double calc_med_rf(std::vector<double>& x){
  	double F;
  	const int sz=x.size(),middle=sz/2-1;
  	if(sz%2==0){
   		std::nth_element(x.begin(),x.begin()+middle,x.end());
    	F=(x[middle]+*(std::min_element(x.begin()+middle+1,x.end())))/2.0;
  	}else{
    	std::nth_element(x.begin(),x.begin()+middle+1,x.end());
    	F=x[middle+1];
  	}
  	return F;
}

static arma::mat sqrt_mat(arma::mat x);

arma::mat calc_dist_rf(arma::mat xnew, arma::mat x, const bool is_euclidean) {
  	const int n=x.n_cols,nu=xnew.n_cols;
  	arma::mat disa(n,nu,arma::fill::zeros);
  	if(is_euclidean) 
    	for (int i=0;i<nu;++i)
      		disa.col(i)=sqrt_mat(sum(square(x.each_col() - xnew.col(i)),0));
  	else 
    	for(int i=0;i<nu;++i)
      		disa.col(i)=(sum(abs(x.each_col() - xnew.col(i)),0)).t();
  	return disa;
}

static arma::mat sqrt_mat(arma::mat x){
	arma::colvec f(x.n_elem);
  	for(double *start=&x[0],*startf=&f[0],*end=&(*x.end());start!=end;++start,++startf){
    	*startf=std::sqrt(*start);    
  	}
  	return f;
}

static arma::uvec get_k_indices(arma::rowvec x, const int k);

arma::umat calc_dist_mem_eff_rf(arma::mat xnew, arma::mat x, const unsigned int k, const bool is_euclidean) {
  const unsigned int nu=xnew.n_cols;
  arma::umat disa(k,nu);
  if(is_euclidean){
    for(unsigned int i=0;i<nu;++i)
      disa.col(i)=get_k_indices(arma::sum(arma::square(x.each_col() - xnew.col(i)),0),k);
  }else{
    for(unsigned int i=0;i<nu;++i)
      disa.col(i)=get_k_indices(arma::sum(arma::abs(x.each_col() - xnew.col(i)),0),k);
  }
  return disa;
}

static arma::uvec get_k_indices(arma::rowvec x, const int k){
	arma::uvec ind=arma::linspace<arma::uvec>(1,x.size(),x.size());
	std::sort(ind.begin(),ind.end(),[&](int i,int j){return x[i-1]<x[j-1];});
  return ind(arma::span(0,k-1));
}

// Author: Giorgos Borboudakis

class TestResult {
public: double pvalue;
    double logpvalue;
    double stat;
    int df;

    TestResult(double _pvalue, double _stat, double _logpvalue, int _df) {
        pvalue=_pvalue;
        stat=_stat;
        logpvalue=_logpvalue;
        df=_df;
    }
};

TestResult g2_test(arma::mat& data, const unsigned int x, const unsigned int y, 
		arma::uvec& cs, const unsigned int ncs, arma::uvec& dc);

static double g2_statistic(arma::uvec& counts, 
		const unsigned int xdim, const unsigned int ydim);

TestResult g2_test(arma::mat& data, const unsigned int x, const unsigned int y, arma::uvec& dc);

Rcpp::List g2_test_univ(arma::mat& data, arma::uvec& dc) {
	const unsigned int nvars = data.n_cols;
    unsigned int nout = nvars * (nvars - 1) / 2;
	arma::uvec xout(nout);
	arma::uvec yout(nout);
	arma::vec statistics(nout);
	arma::vec df(nout);

    unsigned int idx = 0;
    for(unsigned int i = 0; i < nvars; ++i) {
        for(unsigned int j = i + 1; j < nvars; ++j) {
            TestResult result = g2_test(data, i, j, dc);
            xout.at(idx) = i;
            yout.at(idx) = j;
            statistics.at(idx) = result.stat;
            df.at(idx) = (dc.at(i) - 1) * (dc.at(j) - 1);
            ++idx;
        }
    }

    Rcpp::List out;
    out["statistic"] = statistics;
    out["x"] = xout;
    out["y"] = yout;
    out["df"] = df;
    return out;
}

TestResult g2_test(arma::mat& data, const unsigned int x, const unsigned int y, arma::uvec& dc) {
	const unsigned int xdim = dc.at(x);
	const unsigned int ydim = dc.at(y);
	arma::uvec counts(xdim * ydim, arma::fill::zeros);

	for (unsigned int i = 0; i < data.n_rows; ++i) {
		const unsigned int curx = (unsigned int) data.at(i, x);
		const unsigned int cury = (unsigned int) data.at(i, y);
		counts(cury * xdim + curx)++;
	}
	const int df = (xdim - 1) * (ydim - 1);
	const double statistic = g2_statistic(counts, xdim, ydim);

	return TestResult(0, statistic, 0, df);
}

static double g2_statistic(arma::uvec& counts, 
		const unsigned int xdim, const unsigned int ydim) {
	if (arma::all(counts == 0)) {
		return 0;
	}
	double statistic = 0;
	int countsXY = 0;
	arma::uvec countsX(xdim, arma::fill::zeros);
	arma::uvec countsY(ydim, arma::fill::zeros);

	for (unsigned int x = 0; x < xdim; ++x) {
		for (unsigned int y = 0; y < ydim; ++y) {
			const unsigned int curcounts = counts.at(y * xdim + x);
			countsXY += curcounts;
			countsX.at(x) += curcounts;
			countsY.at(y) += curcounts;
		}
	}

	for (unsigned int x = 0; x < xdim; ++x) {
		if (countsX.at(x) != 0) {
			for (unsigned int y = 0; y < ydim; ++y) {
				const unsigned int curcounts = counts(y * xdim + x);
				if (countsY.at(y) && curcounts) {
					statistic += curcounts * (std::log(((double) curcounts * countsXY) / ((double) countsX.at(x) * countsY.at(y))));
				}
			}
		}
	}
	return 2 * statistic;
}

Rcpp::List g2_test(arma::mat& data, const unsigned int x, const unsigned int y, 
		arma::uvec& cs, arma::uvec& dc) {
    TestResult result = g2_test(data, x, y, cs, cs.size(), dc);
    Rcpp::List out;
    out["statistic"] = result.stat;
    out["df"] = result.df;
    return out;
}

TestResult g2_test(arma::mat& data, const unsigned int x, const unsigned int y, 
		arma::uvec& cs, const unsigned int ncs, arma::uvec& dc) {
	if (!ncs) {
		return g2_test(data, x, y, dc);
	}
	const unsigned int xdim = dc.at(x);
	const unsigned int ydim = dc.at(y);
	const unsigned int nsamples = data.n_rows;
	arma::uvec prod(ncs + 1);
	prod.at(0) = 1;
	for (unsigned int i = 1; i <= ncs; ++i) {
		prod.at(i) = prod.at(i - 1) * dc.at(cs.at(i - 1));
	}

	const unsigned int size = prod.at(ncs);
	arma::umat counts(xdim * ydim, size, arma::fill::zeros);
	for (unsigned int i = 0; i < nsamples; ++i) {
		unsigned int key = 0;
		for (unsigned int j = 0; j < ncs; ++j) {
			key += (unsigned int) data.at(i, cs.at(j)) * prod.at(j);
		}
		const unsigned int curx = (unsigned int) data.at(i, x);
		const unsigned int cury = (unsigned int) data.at(i, y);
		counts(cury * xdim + curx, key)++;
	}

	double statistic = 0;
	for (unsigned int i = 0; i < size; ++i) {
		arma::uvec tmp = counts.col(i);
		statistic += g2_statistic(tmp, xdim, ydim);
	}
	const unsigned int df = (xdim - 1) * (ydim - 1) * prod.at(ncs);

	return TestResult(0, statistic, 0, df);
}

void random_contigency_table(int* matrix, const int* nrowt, const int* ncolt, 
		const unsigned int nrow, const unsigned int ncol, const double* logfact, int* jwork, 
		const int ntotal, std::mt19937&	rng);

static int total_counts(arma::uvec& counts, const unsigned int xdim, 
		const unsigned int ydim);

static void col_counts(arma::uvec& counts, const unsigned int xdim,
		const unsigned int ydim, int* counts_y);

static void row_counts(arma::uvec& counts, const unsigned int xdim, 
		const unsigned int ydim, int* counts_x);

TestResult perm_g2_test(arma::mat& data, const unsigned int x, const unsigned int y,
		arma::uvec& cs, const unsigned int ncs, arma::uvec& dc, const unsigned int nperm);

Rcpp::List g2_test_perm(arma::mat& data, const unsigned int x, const unsigned int y,
		arma::uvec& cs, arma::uvec& dc, const unsigned int nperm) {
	TestResult result = perm_g2_test(data, x, y, cs, cs.size(), dc, nperm);
	Rcpp::List out;
	out["statistic"] = result.stat;
	out["pvalue"] = result.pvalue;
	out["x"] = x;
	out["y"] = y;
	out["df"] = result.df;
	return out;
}

TestResult perm_g2_test(arma::mat& data, const unsigned int x, const unsigned int y,
		arma::uvec& cs, const unsigned int ncs, arma::uvec& dc, const unsigned int nperm) {
	const unsigned int xdim = dc.at(x);
	const unsigned int ydim = dc.at(y);

	const unsigned int nsamples = data.n_rows;
	arma::uvec prod(ncs + 1);
	prod.at(0) = 1;
	for (unsigned int i = 1; i <= ncs; i++) {
		prod.at(i) = prod.at(i - 1) * dc.at(cs.at(i - 1));
	}

	const unsigned int size = prod.at(ncs);
	arma::umat counts(xdim * ydim, size, arma::fill::zeros);
	for (unsigned int i = 0; i < nsamples; i++) {
		unsigned int key = 0;
		for (unsigned int j = 0; j < ncs; j++) {
			key += data.at(i, cs.at(j)) * prod.at(j);
		}
		const unsigned int curx = data.at(i, x);
		const unsigned int cury = data.at(i, y);
		counts.at(cury * xdim + curx, key)++;
	}

	double statistic = 0;
	for (unsigned int i = 0; i < size; i++) {
		arma::uvec tmp = counts.col(i);	
		statistic += g2_statistic(tmp, dc.at(x), dc.at(y));
	}
	const int df = (dc.at(x) - 1) * (dc.at(y) - 1) * prod.at(ncs);

	if (!nperm) {
		return TestResult(0, statistic, 0, df);
	}

	arma::vec permstats(nperm, arma::fill::zeros);

	std::random_device rd;
	std::mt19937 rng(rd());

	int* jwork = new int[ydim - 1];
	int* ct = new int[xdim * ydim];
	int* rowcounts = new int[xdim];
	int* colcounts = new int[ydim];
	int* totals = new int[size];
	double* nrc = new double[xdim * ydim];

	int maxtotal = 0;
	for (unsigned int i = 0; i < size; i++) {
		arma::uvec tmp = counts.col(i);
		totals[i] = total_counts(tmp, xdim, ydim);
		maxtotal = (totals[i] > maxtotal) ? totals[i] : maxtotal;
	}

	double* plog = new double[1 + maxtotal];
	double* logfact = new double[1 + maxtotal];
	plog[0] = 0;
	logfact[0] = 0;
	for (int j = 1; j <= maxtotal; j++) {
		plog[j] = std::log((double) j);
		logfact[j] = logfact[j - 1] + plog[j];
	}

	for (unsigned int i = 0; i < size; i++) {
		const int ntotal = totals[i];
		if (ntotal > 0) {
			arma::uvec tmp = counts.col(i);
			row_counts(tmp, xdim, ydim, rowcounts);
			col_counts(tmp, xdim, ydim, colcounts);
			int ctr = 0;
			for (unsigned int x = 0; x < xdim; x++) {
				for (unsigned int y = 0; y < ydim; y++) {
					nrc[ctr++] = plog[ntotal] - plog[rowcounts[x]] - plog[colcounts[y]];
				}
			}

			for (unsigned int p = 0; p < nperm; p++) {
				memcpy(jwork, colcounts, (ydim - 1) * sizeof(int));
				random_contigency_table(ct, rowcounts, colcounts, xdim, ydim, 
						logfact, jwork, ntotal, rng);

				double curstat = 0;
				int cti = 0;
				for (unsigned int x = 0; x < xdim; x++) {
					if (rowcounts[x]) {
						for (unsigned int y = 0; y < ydim; y++) {
							curstat += ct[cti] * (plog[ct[cti]] + nrc[cti]);
							cti++;
						}
					}
					else {
						cti += ydim;
					}
				}
				permstats.at(p) += (2 * curstat);
			}
		}
	}
		
	double pvalue = 1;
	for (unsigned int p = 0; p < nperm; p++) {
		pvalue += (permstats.at(p) >= statistic);
	}
	pvalue /= (nperm + 1);

	return TestResult(pvalue, statistic, std::log(pvalue), df);
}

static void row_counts(arma::uvec& counts, const unsigned int xdim, 
		const unsigned int ydim, int* counts_x) {
	if (arma::all(counts == 0)) {
		return;
	}
  	memset(counts_x, 0, xdim * sizeof(int));
	for (unsigned int x = 0; x < xdim; x++) {
		for (unsigned int y = 0; y < ydim; y++) {
			counts_x[x] += counts.at(y * xdim + x);
		}
	}
}

static void col_counts(arma::uvec& counts, const unsigned int xdim,
		const unsigned int ydim, int* counts_y) {
	if (arma::all(counts == 0)) {
		return;
	}
	memset(counts_y, 0, ydim * sizeof(int));
	for (unsigned int x = 0; x < xdim; x++) {
		for (unsigned int y = 0; y < ydim; y++) {
			counts_y[y] += counts.at(y * xdim + x);
		}
	}
}

static int total_counts(arma::uvec& counts, const unsigned int xdim, 
		const unsigned int ydim) {
	if (arma::all(counts == 0)) {
		return 0;
	}
	int counts_xy = 0;
	for (unsigned int x = 0; x < xdim; x++) {
		for (unsigned int y = 0; y < ydim; y++) {
			counts_xy += counts.at(y * xdim + x);
		}
	}
	return counts_xy;
}

void random_contigency_table(int* matrix, const int* nrowt, const int* ncolt, const unsigned int nrow, const unsigned int ncol, const double* logfact, int* jwork, 
		const int ntotal, std::mt19937& rng) {
	std::uniform_real_distribution<> dist(0, 1);
	int jc = ntotal;
	int ib = 0;
  
	for (unsigned int l = 0; l < nrow - 1; ++l) {
		int ia = nrowt[l];
		int ic = jc;
		jc -= ia;
		for (unsigned int m = 0; m < ncol - 1; ++m) {
			int id = jwork[m];
			int ie = ic;
			ib = ie - ia;
			int ii = ib - id;
			ic -= id;
		  
			if (ie == 0) {
				ia = 0;
				memset(&matrix[l * ncol + m], 0, (ncol - m) * sizeof(int));
				break;
			}
		  
			//  Compute the conditional expected value of MATRIX(L,M)
		  	bool done = false;
		  	int curnlm, nlm;
		  	curnlm = nlm = (int)(((double)ia * id) / ie + 0.5);
		  	double x = exp(logfact[ia] + logfact[ib] + logfact[ic] + logfact[id] - 
					logfact[ie] - logfact[nlm] - logfact[id - nlm] - logfact[ia - nlm] - logfact[ii + nlm]);
			for (double r = dist(rng), sumprb = x; !done && r > x; r = sumprb * dist(rng)) {
				bool lsp = false;
				int nll;
				double curx, y;
				curnlm = nll = nlm;
				curx = y = sumprb = x;
				
				// Increment entry in row L, column M.
				do {
					int j = (id - curnlm) * (ia - curnlm);
					if (j == 0) {
						lsp = true;
						for (j = nll * (ii + nll); j != 0; j = nll * (ii + nll)) {
							--nll;
							y = (y * j) / ((id - nll) * (ia - nll));
							sumprb = sumprb + y;

							if (r <= sumprb) {
							curnlm = nll;
							done = true;
							break;
							}
						}
					}
					else {
						++curnlm;
						curx = (curx * j) / (curnlm * (ii + curnlm));
						sumprb = sumprb + curx;

						if (r <= sumprb) {
							done = true;
							break;
						}
						for (j = nll * (ii + nll); j != 0; j = nll * (ii + nll)) {
							--nll;
							y = (y * j) / ((id - nll) * (ia - nll));
							sumprb = sumprb + y;
					  
							if (r <= sumprb) {
								curnlm = nll;
								done = true;
								break;
							}
							if (!lsp) {
								break;
							}
						}
					}
				} while (!done && !lsp);
			}
			matrix[l * ncol + m] = curnlm;
			ia = ia - curnlm;
			jwork[m] = jwork[m] - curnlm;
		}
    	matrix[l * ncol + ncol - 1] = ia;
  	}
  	memcpy(&matrix[(nrow - 1) * ncol], jwork, (ncol - 1) * sizeof(int));
  	matrix[(nrow - 1) * ncol + ncol - 1] = ib - matrix[(nrow - 1) * ncol + ncol - 2];
}

// Author: Manos Papadakis

NumericVector logistic_only(NumericMatrix& X,NumericVector& Y,const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=2;
  unsigned int j,i;
  colvec b_old(d),b_new(d),L1(d),yhat(n),expyhat,y(Y.begin(),n,false),W(n,fill::zeros),p(n),z_tr_i(n),x2_col;
  mat z(n,2,fill::ones),inv_L2(d,d),z_tr(2,n,fill::ones),x(X.begin(),n,pcols,false);
  NumericVector F(pcols);
  colvec::iterator expyhatiter=expyhat.begin();
  double dif,s,t,sw=0.0,szw=0.0,sz2w=0.0,lgmy=log(my/(1-my));
  for(i=0;i<pcols;++i){
    b_old(0)=lgmy;
    b_old(1)=0;
    z_tr_i=x.col(i);
    z.col(1)=z_tr_i;
    z_tr.row(1)=mat(z_tr_i.begin(),1,n,false);
    x2_col=square(z_tr_i);
    for(dif=1.0,s=0.0;dif>0.000000001;){
      sw=szw=sz2w=0.0;
      yhat = z*b_old;
      expyhat=exp(yhat);
      p = expyhat / ( 1 + expyhat );
      for(j=0;j<n;j++){
        t=p.at(j);
        W.at(j)=t*(1-t);
        sw+=W.at(j);
      }
      szw=sum(W%z_tr_i);
      sz2w=sum(W%x2_col);
      L1=z_tr*(y-p);
      t=1.0/(sw*sz2w-szw*szw);
      inv_L2.at(0,0)=sz2w*t;
      inv_L2.at(0,1)=inv_L2.at(1,0)=-szw*t;
      inv_L2.at(1,1)=sw*t;
      b_new=b_old+inv_L2*L1;
      dif=sum(abs(b_new-b_old));
      b_old=b_new;
    }  
    for(expyhatiter=expyhat.begin();expyhatiter!=expyhat.end();++expyhatiter)
      s+=log1p(*expyhatiter);
    F[i]=2.0*(s-sum(y%yhat));
  }
  return F;
}

NumericVector poisson_only(NumericMatrix& X,NumericVector& Y,const double ylogy,const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=2;
  unsigned int i;
  colvec b_old(d),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones),x(X.begin(),n,pcols,false);
  vec m(n),z_col_1(n);
  NumericVector F(pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(my);
  for(i=0;i<pcols;++i){
    b_old(0)=lgmeany;
    b_old(1)=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    for(dif=1.0;dif>0.000000001;){
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
    F[i]=2.0*(ylogy-sum(y%yhat));
  }
  return F;
}

double glm_logistic(NumericMatrix& X,NumericVector& Y,const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols+1;
  unsigned int j;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),expyhat,y(Y.begin(),n,false),W(n,fill::zeros);
  mat L2,x(X.begin(),n,pcols,false),x_tr(n,pcols+1);
  vec p(n);
  colvec::iterator expyhatiter=expyhat.begin();
  double dif,s=0.0,t;
  x.insert_cols(0,ones(n));
  x_tr=x.t();
  b_old(0)=log(my)-log(1-my);
  for(dif=1.0;dif>0.000000001;){
    yhat = x*b_old;
    expyhat=exp(yhat);
    p = expyhat / ( 1 + expyhat );
    for(j=0;j<n;++j){
      t=p.at(j);
      W.at(j)=t-t*t;
    }
    L1=x_tr*(y-p);
    L2=x.each_col()%W;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  t=sum(y%yhat);
  for(expyhatiter=expyhat.begin();expyhatiter!=expyhat.end();++expyhatiter)
    s+=log1p(*expyhatiter);
  return 2.0*(s-t);
}

double arma_glm_logistic(mat x,vec y,const double my) {
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=pcols+1;
  unsigned int j;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),expyhat,W(n,fill::zeros);
  mat L2,x_tr(n,pcols+1);
  vec p(n);
  colvec::iterator expyhatiter=expyhat.begin();
  double dif,s=0.0,t;
  x.insert_cols(0,ones(n));
  x_tr=x.t();
  b_old(0)=log(my)-log(1-my);
  for(dif=1.0;dif>0.000000001;){
    yhat = x*b_old;
    expyhat=exp(yhat);
    p = expyhat / ( 1 + expyhat );
    for(j=0;j<n;++j){
      t=p.at(j);
      W.at(j)=t-t*t;
    }
    L1=x_tr*(y-p);
    L2=x.each_col()%W;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  t=sum(y%yhat);
  for(expyhatiter=expyhat.begin();expyhatiter!=expyhat.end();++expyhatiter)
    s+=log1p(*expyhatiter);
  return 2.0*(s-t);
}

double glm_poisson(NumericMatrix& X,NumericVector& Y,const double ylogy,const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols+1;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false),m(n);
  mat L2,x(X.begin(),n,pcols,false),x_tr(n,pcols+1);
  double dif;
  x.insert_cols(0,ones(n));
  b_old(0)=log(my);
  x_tr=x.t();
  for(dif=1.0;dif>0.000000001;){
    yhat=x*b_old;
    m=exp(yhat);
    L1=x_tr*(y-m);
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  return 2.0*(ylogy-sum(y%yhat));
}

double arma_glm_poisson(mat x,vec y,const double ylogy,const double my){
  const unsigned int n=x.n_rows,pcols=x.n_cols,d=pcols+1;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),m(n);
  mat L2,x_tr(n,pcols+1);
  double dif;
  x.insert_cols(0,ones(n));
  b_old(0)=log(my);
  x_tr=x.t();
  for(dif=1.0;dif>0.000000001;){
    yhat=x*b_old;
    m=exp(yhat);
    L1=x_tr*(y-m);
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  return 2.0*(ylogy-sum(y%yhat));
}

static double calc_stat(arma::mat& l2, arma::colvec& u, arma::colvec& p, arma::colvec& b, 
		const double dof, const unsigned int nrows, const unsigned int ncols) {
	arma::mat l2_inv = arma::inv(l2);
	double sum = 0;
	for (unsigned int i = 0; i < nrows; ++i) {
		sum += std::pow(u[i], 2) / p[i] / (1 - p[i]);
	}
	const unsigned int last_ind = ncols - 1;
	return std::pow(b[last_ind], 2) / (sum / dof * l2_inv.at(last_ind, last_ind));
}

static Rcpp::NumericVector finalize(arma::colvec& f, arma::colvec& b, 
		const unsigned int n, const unsigned int size) {
	double fmax = 0;
	double findex = -1;
	for (unsigned int i = 0; i < size; ++i) {
	  if (f[i] > fmax) {
		  fmax = f[i];
		  findex = i;
	  }
	}
	const double bic = b[findex] + 2 * std::log(n);
	return Rcpp::NumericVector::create(bic, fmax, findex);
}

static double calc_bic(arma::colvec& y, arma::colvec& p, const unsigned int n) {
	double sum_a = 0;
	double sum_b = 0;
	for (unsigned int j = 0; j < n; ++j) {
		if (y[j] && p[j]) {
			sum_a += y[j] * std::log(y[j] / p[j]);
		}
		if (y[j] != 1 && p[j] != 1) {
			sum_b += (1 - y[j]) * std::log((1 - y[j]) / (1 - p[j]));
		}
	}
	return 2 * sum_a + 2 *sum_b;
}

NumericVector qs_binom_only(NumericMatrix& X, NumericVector& Y, const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=2;
  unsigned int j,i;
  colvec b_old(d),b_new(d),L1(d),yhat(n),expyhat,y(Y.begin(),n,false),W(n,fill::zeros),p(n),z_tr_i(n),x2_col,u(n,fill::zeros),u2(n);
  mat z(n,2,fill::ones),inv_L2(d,d),z_tr(2,n,fill::ones),x(X.begin(),n,pcols,false);
  arma::colvec f(pcols);
  arma::colvec b(pcols);
  double dif,s,t,sw=0.0,szw=0.0,sz2w=0.0,lgmy=log(my/(1-my)),vb;
  for(i=0;i<pcols;++i){
    b_old(0)=lgmy;
    b_old(1)=0;
    z_tr_i=x.col(i);
    z.col(1)=z_tr_i;
    z_tr.row(1)=mat(z_tr_i.begin(),1,n,false);
    x2_col=square(z_tr_i);
    for(dif=1.0,s=0.0;dif>0.000000001;){
      sw=szw=sz2w=0.0;
      yhat = z*b_old;
      expyhat=exp(yhat);
      p = expyhat / ( 1 + expyhat );
      for(j=0;j<n;j++){
        t=p.at(j);
        W.at(j)=t*(1-t);
        sw+=W.at(j);
      }
      szw=sum(W%z_tr_i);
      sz2w=sum(W%x2_col);
      u=y-p;
      L1=z_tr*u;
      t=1.0/(sw*sz2w-szw*szw);
      inv_L2.at(0,0)=sz2w*t;
      inv_L2.at(0,1)=inv_L2.at(1,0)=-szw*t;
      inv_L2.at(1,1)=sw*t;
      b_new=b_old+inv_L2*L1;
      dif=sum(abs(b_new-b_old));
      b_old=b_new;
    }  
    u2=square(u);
    s=sum(u2/W)/(n-2);
    vb=s*sw/(sw*sz2w-szw*szw);
    f[i]=b_new[1]*b_new[1]/vb;
	b[i] = calc_bic(y,p,n);
  }
  return finalize(f,b,n,pcols);
}

NumericVector qs_poisson_only(NumericMatrix& X, NumericVector& Y, const double ylogy, const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=2;
  unsigned int i;
  colvec b_old(d),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false);
  mat z(n,2,fill::ones),inv_L2(d,d),ytr=y.t(),z_tr(2,n,fill::ones),x(X.begin(),n,pcols,false);
  vec m(n),z_col_1(n);
  arma::colvec f(pcols);
  arma::colvec b(pcols);
  double dif,sm=0.0,szm=0.0,sz2m=0.0,t,lgmeany=log(my);
  for(i=0;i<pcols;++i){
    b_old(0)=lgmeany;
    b_old(1)=0;
    z_col_1=x.col(i);
    z.col(1)=z_col_1;
    z_tr.row(1)=mat(z_col_1.begin(),1,n,false);
    for(dif=1.0;dif>0.000000001;){
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
	f[i]=b_new[1]*b_new[1]/phi/inv_L2.at(1,1);
    b[i]=2.0*(ylogy-sum(y%yhat));
  }
  return finalize(f,b,n,pcols);
}

NumericVector glm_qs_binom(NumericMatrix& X,NumericVector& Y,const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols+1;
  unsigned int j,itters=0;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),expyhat,y(Y.begin(),n,false),W(n,fill::zeros);
  mat L2,x(X.begin(),n,pcols,false),x_tr(n,pcols+1);
  vec p(n);
  double dif,t;
  x.insert_cols(0,ones(n));
  x_tr=x.t();
  b_old(0)=log(my)-log(1-my);
  for(dif=1.0;dif>0.000000001;++itters){
    yhat = x*b_old;
    expyhat=exp(yhat);
    p = expyhat / ( 1 + expyhat );
    for(j=0;j<n;++j){
      t=p.at(j);
      W.at(j)=t-t*t;
    }
    L1=x_tr*(y-p);
    L2=x.each_col()%W;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  arma::colvec u=y-p;
  const double dof=n-pcols-1;
  const double stat=calc_stat(L2,u,p,b_new,dof,n,pcols+1);
  const double bic=calc_bic(y,p,n);
  return Rcpp::NumericVector::create(bic,stat);
}

NumericVector glm_qs_poisson(NumericMatrix& X,NumericVector& Y,const double ylogy,const double my){
  const unsigned int n=X.nrow(),pcols=X.ncol(),d=pcols+1;
  colvec b_old(d,fill::zeros),b_new(d),L1(d),yhat(n),y(Y.begin(),n,false),m(n);
  mat L2,x(X.begin(),n,pcols,false),x_tr(n,pcols+1);
  double dif;
  x.insert_cols(0,ones(n));
  b_old(0)=log(my);
  x_tr=x.t();
  for(dif=1.0;dif>0.000000001;){
    yhat=x*b_old;
    m=exp(yhat);
    L1=x_tr*(y-m);
    L2=x.each_col()%m;
    L2=x_tr*L2;
    b_new=b_old+solve(L2,L1);
    dif=sum(abs(b_new-b_old));
    b_old=b_new;
  }
  arma::mat l2_inv=arma::inv(L2);
  const double phi=arma::sum(arma::pow(y-m,2)/m)/(n-d); 
  const double stat=b_new[b_new.size()-1]*b_new[b_new.size()-1]/phi/l2_inv.at(l2_inv.n_rows-1,l2_inv.n_cols-1);
  const double bic=2.0*(ylogy-sum(y%yhat))+d*std::log(n);
  return Rcpp::NumericVector::create(bic, stat);
}
