// Author:  Marios Dimitriadis
// Contact: kmdimitriadis@gmail.com

#include "sw_regs.h"

#define DEBUG 0
#define db_print(...) \
	do { if (DEBUG) Rprintf(__VA_ARGS__); } while (0)

#define DEBUG_L 0
#define dbl_print(...) \
	do { if (DEBUG_L) Rprintf(__VA_ARGS__); } while (0)

#define RES_ROWS 1
#define RES_COLS 3
#define IDX_POS 0
#define VAL_POS 1

struct ms_t {
	double sum;
	double gt_one_lsum;
	double gt_one_gsum;
	double mean;
	double add;
};

// Alters ms
static double calc_base_dist(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const bool is_logistic, ms_t& ms) {
	double base_dist = 0;
	ms.sum = 0;
	ms.gt_one_lsum = 0;
	ms.gt_one_gsum = 0;
	ms.mean = 0;
	ms.add = 0;
	if (is_logistic) {
		for (int i = 0; i < y.size(); ++i) {
			ms.sum += y[i];
		}
		ms.mean = ms.sum / ds.nrow();
		base_dist = -2  * (ms.sum * std::log(ms.mean) + (ds.nrow() - ms.sum) * std::log(1 - ms.mean));
	}
	else {
		for (int i = 0; i < y.size(); ++i) {
			ms.sum += y[i];
			if (y[i] > 0) {
				ms.gt_one_lsum += (y[i] * std::log(y[i]));
				ms.gt_one_gsum += std::lgamma(y[i] + 1); }
		}
		ms.mean = ms.sum / y.size();
		double lsum = 0;
		for (int i = 0; i < y.size(); ++i) {
			if (y[i] > 0) {
				lsum += y[i] * std::log(y[i] / ms.mean);
			}
		}
		ms.add = -2 * (ms.gt_one_lsum - ms.sum - ms.gt_one_gsum);
		base_dist = 2 * lsum; 
	}
	return base_dist;
}

static Rcpp::NumericVector calc_type_only(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double gt_one_lsum, const double mean, const bool is_logistic) {
	if (is_logistic) {
		return logistic_only(ds, y, mean);
	}
	return poisson_only(ds, y, gt_one_lsum, mean);
}

static Rcpp::NumericVector calc_min(Rcpp::NumericVector& data) {
	if (!data.size()) {
		Rcpp::stop("Data size invalid.\n");
	}
	Rcpp::NumericVector min_min_col(2);
	min_min_col[0] = -1;
	for (int i = 0; i < data.size(); ++i) {
		if (min_min_col[0] == -1 || data[i] < min_min_col[1]) {
			min_min_col[0] = i;
			min_min_col[1] = data[i];
		}
	}
	return min_min_col;
}

static Rcpp::IntegerVector app_val(Rcpp::IntegerVector& data, const int elem) {
	Rcpp::IntegerVector app_data(data.size() + 1);
	int i;
	for (i = 0; i < data.size(); ++i) {
		app_data[i] = data[i];
	}
	app_data[i] = elem;
	return app_data;
}

static Rcpp::NumericMatrix form_cmat(Rcpp::NumericMatrix& ds,
		Rcpp::IntegerVector& rows, Rcpp::IntegerVector& cols) {
	Rcpp::NumericMatrix formed_ds(rows.size(), cols.size());
	for (int i = 0; i < rows.size(); ++i) {
		for (int j = 0; j < cols.size(); ++j) {
			formed_ds(i, j) = ds(rows[i], cols[j]);
		}
	}
	return formed_ds;
}

static double calc_glm_type(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds,  
		const double gt_one_lsum, const double mean, const bool is_logistic) {
	if (is_logistic) {
		return glm_logistic(ds, y, mean);
	}
	return glm_poisson(ds, y, gt_one_lsum, mean);
}

/*** _calc_bic_fs_reg_cpp_ ***/

static Rcpp::NumericMatrix finalize_bfs(Rcpp::IntegerVector& idxs, 
		Rcpp::NumericVector& bics, Rcpp::NumericVector& y, const double add);

// Alters prev_dist, idxs, bics
static bool update_vals_end_bfs(Rcpp::NumericVector& min_min_col, std::vector<bool>& used_cols,  
		double& prev_dist, const double tol,  const double log_n, const int step,
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& bics);	

// Alters ms
static bool calc_base_bfs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, const bool is_logistic, 
		std::vector<bool>& used_cols, double& prev_dist, 
		const double tol, const double log_n, const int step, ms_t& ms, 
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& bics);

Rcpp::NumericMatrix calc_bic_fs_reg(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double tol, const std::string type) {
	const bool is_logistic = !type.compare("logistic");							   
	const int adj_nrows = ds.nrow() - 15;
	std::vector<bool> used_cols(ds.ncol());
	Rcpp::IntegerVector idxs;
	Rcpp::NumericVector bics;
	double log_n = std::log(ds.nrow());
	int step = 1;
	bool end = false;
	double prev_dist = 0;
	ms_t ms;
	Rcpp::IntegerVector rows = Rcpp::Range(0, ds.nrow() - 1);

	db_print("Calculating base.\n");
	end = calc_base_bfs(y, ds, is_logistic, used_cols, prev_dist, tol, log_n, step, ms, idxs, bics);
	++step;
	if (end) {
		db_print("end\n");
		return finalize_bfs(idxs, bics, y, ms.add);
	}
	db_print("Entering main loop.\n");
	while (step < adj_nrows) {
		db_print("step < adj_nrows\n");
		db_print("step == %d\n", step);
		Rcpp::NumericVector min_min_col(2);
		min_min_col[0] = -1;
		for (int i = 0; i < ds.ncol(); ++i) {
			if (used_cols[i]) {
				continue;
			}
			Rcpp::IntegerVector cols = app_val(idxs, i); 
			Rcpp::NumericMatrix rh = form_cmat(ds, rows, cols);
			const double dist_i = calc_glm_type(y, rh, ms.gt_one_lsum, ms.mean, is_logistic);
			if (min_min_col[0] == -1 || dist_i < min_min_col[1]) {
				min_min_col[0] = i;
				min_min_col[1] = dist_i;
			}
		}
		db_print("Updating values.\nChecking end conditions.\n");
		end = update_vals_end_bfs(min_min_col, used_cols, 
				prev_dist, tol, log_n, step, idxs, bics);
		if (end) {
			db_print("end\n");
			return finalize_bfs(idxs, bics, y, ms.add);
		}
		++step;
	}
	db_print("Exiting main loop.\n");
	return finalize_bfs(idxs, bics, y, ms.add);
}

static bool calc_base_bfs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, const bool is_logistic, 
		std::vector<bool>& used_cols, double& prev_dist, const double tol, const double log_n, 
		const int step, ms_t& ms, 
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& bics) {
	db_print("Calculating base dist.\n");
	prev_dist = calc_base_dist(y, ds, is_logistic, ms);
	db_print("Calculating dist.\n");
	Rcpp::NumericVector dist = calc_type_only(y, ds, ms.gt_one_lsum, ms.mean, is_logistic);
	db_print("Finding min of dist.\n");
	Rcpp::NumericVector min_min_col = calc_min(dist);
	db_print("Updating values.\nChecking end conditions.\n");
	return update_vals_end_bfs(min_min_col, used_cols, prev_dist, tol, log_n, step, 
			idxs, bics);
}

static bool update_vals_end_bfs(Rcpp::NumericVector& min_min_col, std::vector<bool>& used_cols,  
		double& prev_dist, const double tol, const double log_n, const int step,
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& bics) {	
	if (min_min_col[0] == -1) {
		db_print("min_min_col[0] == -1\n");
		return true;
	}
	db_print("Calculating BIC.\n");
	const double bic = min_min_col[1] + (step + 1) * log_n;
	if (bics.size() && (bics[bics.size() - 1] - bic) <= tol) {
		db_print("bics.size() && (bics[bics.size() - 1] - bic) <= tol\n");
		return true;
	}
	db_print("Updating.\n");
	used_cols[min_min_col[0]] = true;
	idxs.push_back(min_min_col[0]);
	bics.push_back(bic);
	prev_dist = min_min_col[1];
	return false;
}					 

static Rcpp::NumericMatrix finalize_bfs(Rcpp::IntegerVector& idxs, Rcpp::NumericVector& bics, 
		Rcpp::NumericVector& y, const double add) {
	Rcpp::NumericMatrix ret(idxs.size(), 2);
	for (int i = 0; i < idxs.size(); ++i) {
		ret.at(i, 0) = idxs[i] + 1;
		ret.at(i, 1) = bics[i] + add;
	}
	return ret;
}

/*** _calc_fs_reg_cpp_ ***/

static Rcpp::NumericMatrix finalize_fs(Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics,
		Rcpp::NumericVector& y, const double add);

// Alters prev_dist, idxs, stats, pvalues, bics
static bool update_vals_end_fs(Rcpp::NumericVector& min_min_col, std::vector<bool>& used_cols,  
		double& prev_dist, const double log_sig, 
		const double tol,  const double log_n, const int step,
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics);	

// Alters ms
static bool calc_base_fs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, const bool is_logistic, 
		std::vector<bool>& used_cols, double& prev_dist, 
		const double log_sig, const double tol, const double log_n, 
		const int step, ms_t& ms, 
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats,
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics);

Rcpp::NumericMatrix calc_fs_reg_st(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double sig, const double tol, const std::string type) {
	const bool is_logistic = !type.compare("logistic");							   
	const double log_sig = std::log(sig);
	const int adj_nrows = ds.nrow() - 15;
	std::vector<bool> used_cols(ds.ncol());
	Rcpp::IntegerVector idxs;
	Rcpp::NumericVector stats;
 	Rcpp::NumericVector pvalues;
	Rcpp::NumericVector bics;
	double log_n = std::log(ds.nrow());
	int step = 1;
	bool end = false;
	double prev_dist = 0;
	ms_t ms;
	Rcpp::IntegerVector rows = Rcpp::Range(0, ds.nrow() - 1);

	db_print("Calculating base.\n");
	end = calc_base_fs(y, ds, is_logistic, used_cols, prev_dist, log_sig, tol, log_n, step, ms, idxs, stats, pvalues, bics);
	++step;
	if (end) {
		db_print("end\n");
		return finalize_fs(idxs, stats, pvalues, bics, y, ms.add);
	}
	db_print("Entering main loop.\n");
	while (step < adj_nrows) {
		db_print("step < adj_nrows\n");
		db_print("step == %d\n", step);
		Rcpp::NumericVector min_min_col(2);
		min_min_col[0] = -1;
		for (int i = 0; i < ds.ncol(); ++i) {
			if (used_cols[i]) {
				continue;
			}
			Rcpp::IntegerVector cols = app_val(idxs, i); 
			Rcpp::NumericMatrix rh = form_cmat(ds, rows, cols);
			const double dist_i = calc_glm_type(y, rh, ms.gt_one_lsum, ms.mean, is_logistic);
			if (min_min_col[0] == -1 || dist_i < min_min_col[1]) {
				min_min_col[0] = i;
				min_min_col[1] = dist_i;
			}
		}
		db_print("Updating values.\nChecking end conditions.\n");
		end = update_vals_end_fs(min_min_col, used_cols, 
				prev_dist, log_sig, tol, log_n, step, 
				idxs, stats, pvalues, bics);
		if (end) {
			db_print("end\n");
			return finalize_fs(idxs, stats, pvalues, bics, y, ms.add);
		}
		++step;
	}
	db_print("Exiting main loop.\n");
	return finalize_fs(idxs, stats, pvalues, bics, y, ms.add);
}

static bool calc_base_fs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, const bool is_logistic, 
		std::vector<bool>& used_cols, double& prev_dist, 
		const double log_sig, const double tol, const double log_n, 
		const int step, ms_t& ms, 
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats,
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics) {
	db_print("Calculating base dist.\n");
	prev_dist = calc_base_dist(y, ds, is_logistic, ms);
	db_print("Calculating dist.\n");
	Rcpp::NumericVector dist = calc_type_only(y, ds, ms.gt_one_lsum, ms.mean, is_logistic);
	db_print("Finding min of dist.\n");
	Rcpp::NumericVector min_min_col = calc_min(dist);
	db_print("Updating values.\nChecking end conditions.\n");
	return update_vals_end_fs(min_min_col, used_cols, prev_dist, log_sig, tol, log_n, step, 
			idxs, stats, pvalues, bics);
}

static bool update_vals_end_fs(Rcpp::NumericVector& min_min_col, std::vector<bool>& used_cols,  
		double& prev_dist, const double log_sig, 
		const double tol,  const double log_n, const int step,
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics) {	
	if (min_min_col[0] == -1) {
		db_print("min_min_col[0] == -1\n");
		return true;
	}
	db_print("Calculating statistic.\n");
	const double stat = prev_dist - min_min_col[1];
	db_print("Calculating pvalue.\n");
	const double pvalue = R::pchisq(stat, 1, false, true);
	if (pvalue >= log_sig) {
		db_print("pvalue >= log_sig\n");
		return true;
	}
	db_print("Calculating BIC.\n");
	const double bic = min_min_col[1] + (step + 1) * log_n;
	if (bics.size() && (bics[bics.size() - 1] - bic) <= tol) {
		db_print("bics.size() && (bics[bics.size() - 1] - bic) <= tol\n");
		return true;
	}
	db_print("Updating.\n");
	used_cols[min_min_col[0]] = true;
	idxs.push_back(min_min_col[0]);
	bics.push_back(bic);
	stats.push_back(stat);
	pvalues.push_back(pvalue);
	prev_dist = min_min_col[1];
	return false;
}					 

static Rcpp::NumericMatrix finalize_fs(Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics, 
		Rcpp::NumericVector& y, const double add) {
	Rcpp::NumericMatrix ret(idxs.size(), 4);
	for (int i = 0; i < idxs.size(); ++i) {
		ret(i, 0) = idxs[i] + 1;
		ret(i, 1) = pvalues[i];
		ret(i, 2) = stats[i];
		ret(i, 3) = bics[i] + add;
	}
	return ret;
}

static const int g_BIC_POS = 0;
static const int g_STAT_POS = 1;
static const int g_INDEX_POS = 2;

static Rcpp::NumericVector calc_glm_type_efs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds,  
		const double gt_one_lsum, const double mean, const bool is_quasi_logistic);

// Alters used_cols, idxs, stats, pvalues, bics
static bool update_vals_end_efs(Rcpp::NumericVector& vals, std::vector<bool>& used_cols,  
		const double log_sig, const double tol,  const int step,
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics);

static Rcpp::NumericVector calc_type_only_efs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double gt_one_lsum, const double mean, const bool is_quasi_logistic);

// Alters ms
static void setup_ms_efs(Rcpp::NumericVector& y, const double log_n, const bool is_quasi_logistic, ms_t& ms);

static bool calc_base_efs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, const bool is_quasi_logistic, 
		std::vector<bool>& used_cols, 
		const double log_sig, const double tol, const double log_n, 
		const int step, ms_t& ms, 
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats,
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics);

Rcpp::NumericMatrix calc_fs_reg_ext(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double sig, const double tol, const std::string type) {
	const bool is_quasi_logistic = !type.compare("quasilogistic");							   
	const double log_sig = std::log(sig);
	const int adj_nrows = ds.nrow() - 15;
	std::vector<bool> used_cols(ds.ncol());
	Rcpp::IntegerVector idxs;
	Rcpp::NumericVector stats;
 	Rcpp::NumericVector pvalues;
	Rcpp::NumericVector bics;
	double log_n = std::log(ds.nrow());
	int step = 1;
	bool end = false;
	ms_t ms;
	Rcpp::IntegerVector rows = Rcpp::Range(0, ds.nrow() - 1);

	db_print("Calculating base.\n");
	end = calc_base_efs(y, ds, is_quasi_logistic, used_cols, 
			log_sig, tol, log_n, step, ms, idxs, stats, pvalues, bics);
	++step;
	if (end) {
		db_print("end\n");
		return finalize_fs(idxs, stats, pvalues, bics, y, ms.add);
	}
	db_print("Entering main loop.\n");
	while (step < adj_nrows) {
		db_print("step < adj_nrows\n");
		db_print("step == %d\n", step);
		Rcpp::NumericVector vals = Rcpp::NumericVector::create(0, 0, -1);
		for (int i = 0; i < ds.ncol(); ++i) {
			if (used_cols[i]) {
				continue;
			}
			Rcpp::IntegerVector cols = app_val(idxs, i); 
			Rcpp::NumericMatrix rh = form_cmat(ds, rows, cols);
			Rcpp::NumericVector tmp_vals = calc_glm_type_efs(y, rh, ms.gt_one_lsum, ms.mean, is_quasi_logistic);
			if (vals[g_INDEX_POS] == -1 || tmp_vals[g_STAT_POS] > vals[g_STAT_POS]) {
				vals[g_BIC_POS] = tmp_vals[g_BIC_POS];
				vals[g_STAT_POS] = tmp_vals[g_STAT_POS];
				vals[g_INDEX_POS] = i;
			}
		}
		db_print("Updating values.\nChecking end conditions.\n");
		end = update_vals_end_efs(vals, used_cols, log_sig, tol, step, 
				idxs, stats, pvalues, bics);
		if (end) {
			db_print("end\n");
			return finalize_fs(idxs, stats, pvalues, bics, y, ms.add);
		}
		++step;
	}
	db_print("Exiting main loop.\n");
	return finalize_fs(idxs, stats, pvalues, bics, y, ms.add);
}

static bool calc_base_efs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, const bool is_quasi_logistic, 
		std::vector<bool>& used_cols, 
		const double log_sig, const double tol, const double log_n, const int step, ms_t& ms, 
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics) {
	db_print("Setting up ms.\n");
	setup_ms_efs(y, log_n, is_quasi_logistic, ms);
	db_print("Calculating type only.\n");
	Rcpp::NumericVector vals = calc_type_only_efs(y, ds, ms.gt_one_lsum, ms.mean, is_quasi_logistic);
	db_print("Updating values.\nChecking end conditions.\n");
	return update_vals_end_efs(vals, used_cols, log_sig, tol, step, 
			idxs, stats, pvalues, bics);
}

static void setup_ms_efs(Rcpp::NumericVector& y, const double log_n, const bool is_quasi_logistic, ms_t& ms) {
	ms.sum = 0;
	ms.gt_one_lsum = 0;
	ms.gt_one_gsum = 0;
	ms.mean = 0;
	ms.add = log_n;
	if (is_quasi_logistic) {
		for (int i = 0; i < y.size(); ++i) {
			ms.sum += y[i];
			ms.mean = ms.sum / y.size();
		}
	}
	else {
		for (int i = 0; i < y.size(); ++i) {
			ms.sum += y[i];
			if (y[i] > 0) {
				ms.gt_one_lsum += (y[i] * std::log(y[i]));
				ms.gt_one_gsum += std::lgamma(y[i] + 1); }
		}
		ms.mean = ms.sum / y.size();
		ms.add = -2 * (ms.gt_one_lsum - ms.sum - ms.gt_one_gsum);
	}
}

static Rcpp::NumericVector calc_type_only_efs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds, 
		const double gt_one_lsum, const double mean, const bool is_quasi_logistic) {
	if (is_quasi_logistic) {
		return qs_binom_only(ds, y, mean);
	}
	return qs_poisson_only(ds, y, gt_one_lsum, mean);
}

static bool update_vals_end_efs(Rcpp::NumericVector& vals, std::vector<bool>& used_cols,  
		const double log_sig, const double tol, const int step,
		Rcpp::IntegerVector& idxs, Rcpp::NumericVector& stats, 
		Rcpp::NumericVector& pvalues, Rcpp::NumericVector& bics) {	
	const double pvalue = R::pchisq(vals[g_STAT_POS], 1, false, true);
	if (pvalue >= log_sig) {
		db_print("pvalue < log_sig\n");
		return true;
	}
	if (bics.size() && (bics[bics.size() - 1] - vals[g_BIC_POS]) <= tol) {
		db_print("bics.size() && (bics[bics.size() - 1] - bic) <= tol\n");
		return true;
	}
	db_print("Updating.\n");
	used_cols[vals[g_INDEX_POS]] = true;
	idxs.push_back(vals[g_INDEX_POS]);
	bics.push_back(vals[g_BIC_POS]);
	stats.push_back(vals[g_STAT_POS]);
	pvalues.push_back(pvalue);
	return false;
}					 

static Rcpp::NumericVector calc_glm_type_efs(Rcpp::NumericVector& y, Rcpp::NumericMatrix& ds,  
		const double gt_one_lsum, const double mean, const bool is_quasi_logistic) {
	if (is_quasi_logistic) {
		return glm_qs_binom(ds, y, mean);
	}
	return glm_qs_poisson(ds, y, gt_one_lsum, mean);
}

/*** _bs_reg_ ***/

static double calc_stat_bs(arma::vec& y, arma::mat& ds, const double dist, ms_t ms, const bool is_logistic);

static arma::uvec adj_idxs_bs(arma::uvec& idxs, std::vector<bool>& idxs_used, const unsigned int idxs_used_cntr);

static void adj_vals_bs(arma::vec& types_gen, std::vector<unsigned int>& prev_min_idxs, const double val);

static arma::vec get_idx_min_bs(arma::vec& types_gen);

static double gen_type_bs(arma::vec& y, arma::mat& ds, const ms_t ms, const bool is_logistic);

static arma::vec gen_types_bs(arma::vec& y, arma::mat& ds, arma::uvec& idxs, 
		std::vector<bool> idxs_used, arma::vec& types_gen, 
		std::vector<unsigned int> rm_idxs, const ms_t ms, const bool is_logistic);

static void upd_ms_bs(arma::vec& y, ms_t& ms, const bool is_logistic);

static bool is_type_bs(const std::string type);

Rcpp::List calc_bs_reg(arma::vec& y, arma::mat& ds, const double sig, const std::string type) {
	const bool is_logistic = is_type_bs(type);
	ms_t ms = {0};
	upd_ms_bs(y, ms, is_logistic);
	const double sig_log = std::log(sig);
	arma::uvec idxs(ds.n_cols);
	std::iota(idxs.begin(), idxs.end(), 0);
	std::vector<bool> idxs_used(ds.n_cols);
	std::fill(idxs_used.begin(), idxs_used.end(), false);
	unsigned int idxs_used_cntr = 0;
	arma::vec types_gen(ds.n_cols);
	arma::mat res(RES_ROWS, RES_COLS, arma::fill::zeros);

	db_print("Generating types.\n");
	double type_gen = gen_type_bs(y, ds, ms, is_logistic);
	gen_types_bs(y, ds, idxs, idxs_used, types_gen, std::vector<unsigned int>(), ms, is_logistic);
	arma::vec idx_min = get_idx_min_bs(types_gen);
	double stat = idx_min[VAL_POS] - type_gen;
	db_print("The pre-loop stat is: %f\n", stat);
	double pvalue = R::pchisq(stat, 1, 0, 1);
	db_print("The pre-loop pvalue is: %f\n", pvalue);
	if (pvalue >= sig_log) {
		db_print("pvalue >= sig_log\n");
		unsigned int j = 0;
		res(j, 0) = idx_min[IDX_POS] + 1; res(j, 1) = stat; res(j, 2) = pvalue;
		idxs_used[idx_min[IDX_POS]] = true;
		++idxs_used_cntr;
		std::vector<unsigned int> prev_min_idxs = {(unsigned int) idx_min[IDX_POS]};
		db_print("Entering main loop.\n");
		while (pvalue > sig_log && j < (ds.n_cols - 2) && idxs_used.size() >= idxs_used_cntr) {
			++j;
			type_gen = idx_min[VAL_POS];
			adj_vals_bs(types_gen, prev_min_idxs, type_gen + 500);
			gen_types_bs(y, ds, idxs, idxs_used, types_gen, prev_min_idxs, ms, is_logistic);
			idx_min = get_idx_min_bs(types_gen);
			stat = idx_min[VAL_POS] - type_gen;
			db_print("An in-loop stat is: %f\n", stat);
			pvalue = R::pchisq(stat, 1, 0, 1);
			db_print("An in-loop pvalue is: %f\n", pvalue);
			if (pvalue > sig_log) {
				db_print("pvalue > sig_log\n");
				const unsigned int extra_row = res.n_rows;
				res.resize(extra_row + 1, res.n_cols);
				res(extra_row, 0) = idx_min(IDX_POS) + 1; res(extra_row, 1) = stat; res(extra_row, 2) = pvalue;
				idxs_used[idx_min[IDX_POS]] = true;
				++idxs_used_cntr;
				prev_min_idxs.push_back(idx_min[IDX_POS]);
			}
			if (idxs_used_cntr == (idxs.size() - 1)) {
				db_print("idxs_used_cntr == (idxs.size() - 1)\n");
				arma::uvec left_idxs = adj_idxs_bs(idxs, idxs_used, idxs_used_cntr);
				arma::mat tmp_ds = ds.col(left_idxs[0] - 1);
				const double tmp_dist = gen_type_bs(y, tmp_ds, ms, is_logistic);
				stat = calc_stat_bs(y, ds, tmp_dist, ms, is_logistic);
				pvalue = R::pchisq(stat, 1, 0, 1);
				if (pvalue > sig_log) {
					const unsigned int extra_row = res.n_rows;
					res.resize(extra_row + 1, res.n_cols);
					res(extra_row, 0) = left_idxs[0]; res(extra_row, 1) = stat; res(extra_row, 2) = pvalue;
					++idxs_used_cntr;
				}
			}
		}
	}
	Rcpp::List ret;
	arma::uvec left_idxs = adj_idxs_bs(idxs, idxs_used, idxs_used_cntr);
	ret["info"] = res; ret["vars"] = left_idxs;
	return ret;
}

static bool is_type_bs(const std::string type) {
	if (!type.compare("logistic")) {
		return true;
	}
	else if (!type.compare("poisson")) {
		return false;
	}
	stop("Type input invalid. Exiting...\n|");
}

static void upd_ms_bs(arma::vec& y, ms_t& ms, const bool is_logistic) {
	for (unsigned int i = 0; i < y.size(); i++) {
		ms.sum += y[i];
		if (!is_logistic && y[i] > 0) {
			ms.gt_one_lsum += (y[i] * std::log(y[i]));
		}
	}
	ms.mean = ms.sum / y.size();
}

static double gen_type_bs(arma::vec& y, arma::mat& ds, const ms_t ms, const bool is_logistic) {
	if (is_logistic) {
		return arma_glm_logistic(ds, y, ms.mean);
	}
	return arma_glm_poisson(ds, y, ms.gt_one_lsum, ms.mean);
}

static arma::vec gen_types_bs(arma::vec& y, arma::mat& ds, arma::uvec& idxs, 
		std::vector<bool> idxs_used, arma::vec& types_gen, 
		std::vector<unsigned int> rm_idxs, const ms_t ms, const bool is_logistic) {
	for (unsigned int i = 0; i < idxs.size(); ++i) {
		if (idxs_used[i]) {
			continue;
		}
		arma::uvec minus_idxs;
		if (!rm_idxs.size()) {
			minus_idxs = {idxs[i]};
		}
		else {
			std::vector<unsigned int> tmp_idxs = rm_idxs;
			tmp_idxs.push_back(idxs[i]);
			minus_idxs = arma::conv_to<arma::uvec>::from(tmp_idxs);
		}
		arma::mat adj_ds = rm_cols(ds, minus_idxs);
		if (is_logistic) {
			types_gen[idxs[i]] = arma_glm_logistic(adj_ds, y, ms.mean);
		}
		else {
			types_gen[idxs[i]] = arma_glm_poisson(adj_ds, y, ms.gt_one_lsum, ms.mean);
		}
	}
	return types_gen;
}

static arma::vec get_idx_min_bs(arma::vec& types_gen) {
	arma::vec idx_min(2);
	idx_min[IDX_POS] = -1;
	for (unsigned int i = 0; i < types_gen.size(); ++i) {
		if (idx_min[IDX_POS] == -1) {
			idx_min[IDX_POS] = i;
			idx_min[VAL_POS] = types_gen[i];
		}
		else if (idx_min[VAL_POS] > types_gen[i]) {
			idx_min[IDX_POS] = i;
			idx_min[VAL_POS] = types_gen[i];
		}
	}
	return idx_min;
}

static void adj_vals_bs(arma::vec& types_gen, std::vector<unsigned int>& prev_min_idxs, const double val) {
	for (unsigned int i = 0; i < prev_min_idxs.size(); i++) {
		types_gen[prev_min_idxs[i]] = val;
	}
}

static arma::uvec adj_idxs_bs(arma::uvec& idxs, std::vector<bool>& idxs_used, const unsigned int idxs_used_cntr) {
	int idxs_act_cntr = idxs.size() - idxs_used_cntr;
	idxs_act_cntr = std::abs(idxs_act_cntr);
	arma::uvec left_idxs(idxs_act_cntr);
	unsigned int i = 0;
	unsigned int j = 0;
	while (idxs_act_cntr) {
		if (!idxs_used[i]) {
			left_idxs[j++] = idxs[i] + 1;
			--idxs_act_cntr;
		}
		++i;
	}
	return left_idxs;
}

static double calc_stat_bs(arma::vec& y, arma::mat& ds, const double dist, ms_t ms, const bool is_logistic) {
	double tmp = 0;
	double init = 0;
	if (is_logistic) {
		tmp = ms.sum / ds.n_rows;
		init = -2 * (ds.n_rows * tmp * std::log(tmp) + (ds.n_rows - ds.n_rows * tmp) * std::log(1 - tmp));
	}
	else {
		tmp = ms.sum / ds.n_rows;
		init = 2 * ms.gt_one_lsum - 2 * ds.n_rows * tmp * std::log(tmp);
	}
	return init - dist;
}
