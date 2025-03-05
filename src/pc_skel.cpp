// Author:  Marios Dimitriadis
// Contact: kmdimitriadis@gmail.com

#include "pc_skel.h"

#define DEBUG 0
#define db_print(...) \
	do { if (DEBUG) Rprintf(__VA_ARGS__); } while (0)

#define DEBUG_L 0
#define dbl_print(...) \
	do { if (DEBUG_L) Rprintf(__VA_ARGS__); } while (0)

#define INIT_STAT_POS 0
#define INIT_PVALUE_POS 1

struct mth_t {
	bool is_pearson_method;
	bool is_spearman_method;
	bool is_cat_method;
};

static std::vector<unsigned int> loc_inactv_pos(arma::mat& sig_pairs,
		arma::mat& st, const unsigned int row);

// Alters sig_pairs
static unsigned int upd_state(arma::mat& st, arma::mat& sig_pairs, arma::mat& ch, 
		const unsigned int k, const unsigned int el, Rcpp::List& sep);

// Alters st, ch
static void update_st_ch(arma::mat& st, arma::mat& ch, arma::mat& sam, 
		arma::vec& cat_condi, const unsigned int spo_i0, const unsigned int spo_i1, 
		const unsigned int curr_row, const unsigned int m, const unsigned int k);

static arma::vec calc_cat_condi(arma::mat& ds, arma::mat& cor_ds, arma::uvec& max_min, 
		arma::mat& sam, const unsigned int spo_i0, const unsigned int spo_i1, 
		const unsigned int m, const unsigned int k,
		const bool is_cat_method, const std::string method, const unsigned int r);

// Alters sam
static void calc_sam(arma::mat& pvalues, std::vector<unsigned int>& adj,
		const unsigned int spo_row, const unsigned int k, arma::mat& sam);

// Alters xadj, yadj
static void calc_adj(arma::mat& st, const unsigned int xrow, const unsigned int yrow,
		std::vector<unsigned int>& xadj, std::vector<unsigned int>& yadj);

static unsigned int link_vars(arma::mat& ds, arma::mat& cor_ds, arma::uvec& max_min, arma::mat& pvalues, 
		arma::mat& st, arma::mat& sig_pairs, arma::mat& ch, 
		const mth_t mth, const std::string method, const double sig_log, 
		const unsigned int k, const unsigned int r);

// Alters pvalues, st
static arma::mat form_st_sig_pairs(const unsigned int ds_nrows, const unsigned int ds_ncols, const double sig_log, 
		arma::mat& test_stats, arma::mat& pvalues, arma::mat& st);

// Alters test_stats_abs, test_stats, pvalues
static void init_cat_data(arma::mat& ds, arma::uvec& max_min, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues);

// Alters test_stats_abs, test_stats, pvalues
static void init_pearson_spearman_yp_data(arma::mat& ds, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues, const unsigned int r);

// Alters test_stats_abs
static void calc_test_stats(arma::mat& cor_ds, const unsigned int ds_nrows, const double div, arma::mat& test_stats_abs);

// Alters cor_ds, test_stats_abs, test_stats, pvalues
static void init_pearson_spearman_np_data(arma::mat& ds, arma::mat& cor_ds, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues, 
		const bool is_pearson_method);

static void init_data(arma::mat& ds, arma::uvec& max_min, arma::mat& cor_ds, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues, 
		const mth_t mth, const unsigned int r, 
		arma::mat& stats_init, arma::mat& pvalues_init, arma::ivec& is_init_vals);

// Alters ds
static void check_NAs(arma::mat& ds, const unsigned int nrows, const unsigned int ncols, const bool is_cat_method);

Rcpp::List calc_pc_skel(arma::mat& ds, const std::string method, const double sig, const unsigned int r, 
		arma::mat& stats_init, arma::mat& pvalues_init, arma::ivec& is_init_vals) {
	mth_t mth;
	mth.is_pearson_method = !method.compare("pearson");
	mth.is_spearman_method = !method.compare("spearman");
	mth.is_cat_method = !method.compare("cat");
	if (!mth.is_pearson_method && !mth.is_spearman_method && !mth.is_cat_method) {
		Rcpp::stop("Method not recognised.\n");
	}
	const double sig_log = std::log(sig);
	const unsigned int ds_nrows = ds.n_rows;
	const unsigned int ds_ncols = ds.n_cols;

	db_print("Checking for NA values.\n");
	check_NAs(ds, ds_nrows, ds_ncols, mth.is_cat_method);

	if (mth.is_spearman_method) {
		db_print("is_spearman_method\n");
		db_print("Calculating rank.\n");
		ds = calc_rank(ds);
	}
	arma::uvec max_min;
	if (mth.is_cat_method) {
		db_print("is_cat_method\n");
		db_print("Calculating max_min.\n");
		max_min = sub_col_max_min(ds, false);
	}
	
	db_print("Initializing clock.\n");
	clock_t begin = clock();

	db_print("Initializing data based on arguments.\n");
	arma::mat cor_ds(ds_ncols, ds_ncols);
	arma::mat test_stats_abs(ds_ncols, ds_ncols, arma::fill::zeros);
	arma::mat test_stats(ds_ncols, ds_ncols, arma::fill::zeros);
	arma::mat pvalues(ds_ncols, ds_ncols, arma::fill::zeros);
	init_data(ds, max_min, cor_ds, test_stats_abs, test_stats, pvalues, mth, r, 
			stats_init, pvalues_init, is_init_vals);
	pvalues_init = pvalues;

	db_print("Forming matrix sig_pairs/updating pvalues.\n");
	arma::mat st(ds_ncols, ds_ncols, arma::fill::zeros);
	arma::mat sig_pairs = form_st_sig_pairs(ds_nrows, ds_ncols, sig_log, 
			test_stats, pvalues, st);

	Rcpp::List sep;
	std::vector<unsigned int> tests;
	tests.push_back(ds_ncols * (ds_ncols - 1) / 2);
	unsigned int k = 0;

	if (!sig_pairs.n_rows) {
		db_print("!sig_pairs.n_rows\n");
		db_print("Adjusting diagonal of st.\nExiting...\n");
		adj_diag(st, 0);
		Rcpp::List ret; ret["kappa"] = k; ret["G"] = st;
		return ret;
	}
	else {
		db_print("sig_pairs.n_rows\n");
		db_print("Entering main loop.\n");
		unsigned int el = 2;
		while (k < el && sig_pairs.n_rows > 0) {
			k++;
			db_print("k = %u\n", k);
			arma::mat ch(sig_pairs.n_rows, k + 2, arma::fill::zeros);
			db_print("Linking vars.\n");
			unsigned int ntests = link_vars(ds, cor_ds, max_min, pvalues, st, sig_pairs, ch,
					mth, method, sig_log, k, r);
			db_print("Updating state.\n");
			el = upd_state(st, sig_pairs, ch, k, el, sep);
			if (ntests) {
				tests.push_back(ntests);
			}
		}
		db_print("Exited main loop.\nFinalizing...\n");
		st = st / 2;
		adj_diag(st, 0);
		clock_t end = clock();
		double time_spent = (double) (end - begin) / CLOCKS_PER_SEC;
		Rcpp::List ret; 
		ret["stat"] = test_stats_abs; ret["ini.pvalue"] = pvalues_init; ret["pvalue"] = pvalues; 
		ret["runtime"] = time_spent; ret["kappa"] = tests.size() - 1; ret["n.tests"] = tests; 
		ret["G"] = st; ret["sepset"] = sep;
		return ret;
	}
	return NILSXP;
}

static void check_NAs(arma::mat& ds, const unsigned int nrows, const unsigned int ncols, const bool is_cat_method) {
	bool found_NA = false;
	if (!is_cat_method) {
		found_NA = adj_med_NAs(ds);
	}
	else {
		found_NA = adj_freq_NAs(ds);
	}
	if (found_NA) {
		Rcpp::Rcout << "NA values were found in the dataset. "
			<< "They were replaced by the median of the (numeric) column or by the most frequent mode."
			<< std::endl;
	}
}

static void init_data(arma::mat& ds, arma::uvec& max_min, arma::mat& cor_ds, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues, 
		const mth_t mth, const unsigned int r, 
		arma::mat& stats_init, arma::mat& pvalues_init, arma::ivec& is_init_vals) {
	if (!mth.is_cat_method) {
		db_print("!mth.is_cat_method\n");
		db_print("Calculating cor_ds.\n");
		cor_ds = arma::cor(ds);
	}
	if (is_init_vals(INIT_STAT_POS) && is_init_vals(INIT_PVALUE_POS)) {
		db_print("is_init_vals(INIT_STAT_POS) && is_init_vals(INIT_PVALUE_POS)\n");
		pvalues = pvalues_init;
		test_stats_abs = stats_init;
		test_stats = stats_init / (ds.n_rows - 3);
		return;
	}
	if (!mth.is_cat_method) {
		db_print("!is_cat_method\n");
		if (r == 1) {
			db_print("r == 1\n");
			init_pearson_spearman_np_data(ds, cor_ds, test_stats_abs, test_stats, pvalues, mth.is_pearson_method);
		}
		else {
			db_print("r != 1\n"); 
			init_pearson_spearman_yp_data(ds, test_stats_abs, test_stats, pvalues, r);
		}
	}
	else {
		db_print("is_cat_method\n");
		init_cat_data(ds, max_min, test_stats_abs, test_stats, pvalues);
	}
}
		
static void init_pearson_spearman_np_data(arma::mat& ds, arma::mat& cor_ds, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues, 
		const bool is_pearson_method) { 
	db_print("Calculating test_stats_abs.\n");
	if (is_pearson_method) {
		calc_test_stats(cor_ds, ds.n_rows, 1, test_stats_abs);
	}
	else {
		calc_test_stats(cor_ds, ds.n_rows, 1.029563, test_stats_abs);
	}
	db_print("Calculating test_stats and pvalues.\n");
	test_stats = test_stats_abs / (ds.n_rows - 3);
	pvalues = calc_pt(test_stats_abs, ds.n_rows - 3, false, true, std::log(2));
}

static void calc_test_stats(arma::mat& cor_ds, const unsigned int ds_nrows, const double div, arma::mat& test_stats_abs) {
	for (unsigned int i = 0; i < cor_ds.n_rows; i++) {
		for (unsigned int j = 0; j < cor_ds.n_cols; j++) {
			if (i == j) {
				continue;
			}
			const double curr = cor_ds(i, j);
			test_stats_abs(i, j) = std::abs(0.5 * std::log((1 + curr) / (1 - curr)) * std::sqrt(ds_nrows - 3)) / div;
		}
	}
}

static void init_pearson_spearman_yp_data(arma::mat& ds, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues, const unsigned int r) {
	db_print("Calculating test_stats_abs, test_stats, and pvalues.\n");
	arma::mat tmp_test_stats(ds.n_cols, ds.n_cols, arma::fill::zeros);
	arma::mat tmp_pvalues(ds.n_cols, ds.n_cols, arma::fill::zeros);
	for (unsigned int i = 0; i < ds.n_cols - 1; i++) {
		for (unsigned int j = i + 1; j < ds.n_cols; j++) {
			arma::vec ds_c0 = ds.col(i);
			arma::vec ds_c1 = ds.col(j);
			arma::vec pc = calc_perm_cor(ds_c0, ds_c1, r);
			tmp_test_stats(i, j) = pc(0);
			tmp_pvalues(i, j) = std::log(pc(1));
		}
	}
	test_stats_abs = tmp_test_stats + arma::trans(tmp_test_stats);
	test_stats = test_stats_abs / (ds.n_rows - 3);
	pvalues = tmp_pvalues + arma::trans(tmp_pvalues);
}

static void init_cat_data(arma::mat& ds, arma::uvec& max_min, 
		arma::mat& test_stats_abs, arma::mat& test_stats, arma::mat& pvalues) {
	arma::mat tmp_test_stats(ds.n_cols, ds.n_cols, arma::fill::zeros);
	arma::mat tmp_pvalues(ds.n_cols, ds.n_cols, arma::fill::zeros);
	arma::mat tmp_dfs(ds.n_cols, ds.n_cols, arma::fill::zeros);

	db_print("Calculating g2_test.\n");
	Rcpp::List g2_test = g2_test_univ(ds, max_min);
	arma::vec g2_test_stats = g2_test["statistic"];
	arma::vec g2_test_dfs = g2_test["df"];
	arma::uvec rows = g2_test["x"];
	arma::uvec cols = (g2_test["y"]);
	const unsigned int g2_test_stats_size = g2_test_stats.size();
	for (unsigned int i = 0; i < g2_test_stats_size; i++) {
		const unsigned int row = rows(i);
		const unsigned int col = cols(i);
		const double stat = g2_test_stats(i);
		const double df = g2_test_dfs(i);
		tmp_test_stats(row, col) = stat;
		tmp_pvalues(row, col) = R::pchisq(stat, df, false, true);
		tmp_dfs(row, col) = df;
	}
	db_print("Calculating test_stat_abs, test_stats, and pvalues.\n");
	test_stats_abs = tmp_test_stats + arma::trans(tmp_test_stats);
	adj_diag(test_stats_abs, 0);
	test_stats = test_stats_abs / tmp_dfs;
	pvalues = tmp_pvalues + arma::trans(tmp_pvalues);
	
}

static arma::mat form_st_sig_pairs(const unsigned int ds_nrows, const unsigned int ds_ncols, const double sig_log, 
		arma::mat& test_stats, arma::mat& pvalues, arma::mat& st) {
	arma::mat tmp_pvalues(ds_ncols, ds_ncols);
	cp_lt(pvalues, tmp_pvalues, 2);
	arma::mat tmp_sig_pairs(ds_ncols * ds_ncols, 4, arma::fill::zeros);
	unsigned int sp_row = 0;
	bool pvalue_lt_sig = false;
	for (unsigned int i = 0; i < ds_ncols; i++) {
		for (unsigned int j = 0; j < ds_ncols; j++) {
			pvalues(i, j) <= sig_log ? pvalue_lt_sig = true : pvalue_lt_sig = false;
			if (i == j) {
				st(i, j) = -100;
				pvalues(i, j) = 0;
				tmp_pvalues(i, j) = 0;
			}
			else if (pvalue_lt_sig) {
				st(i, j) = 2;
			}
			if (tmp_pvalues(i, j) <= sig_log) {
				tmp_sig_pairs(sp_row, 0) = i;
				tmp_sig_pairs(sp_row, 1) = j;
				tmp_sig_pairs(sp_row, 2) = test_stats(i, j);
				tmp_sig_pairs(sp_row, 3) = tmp_pvalues(i, j);
				sp_row++;
			}
		}
	}
	tmp_sig_pairs.resize(sp_row, 4);
	arma::vec tmp_sig_pairs_c3 = tmp_sig_pairs.col(3);
	arma::uvec ind_order = arma::sort_index(-tmp_sig_pairs_c3);
	arma::mat sig_pairs(sp_row, 4);
	for (unsigned int i = 0; i < sp_row; i++) {
		const unsigned int curr_pos = ind_order(i);
		sig_pairs(i, 0) = tmp_sig_pairs(curr_pos, 0);
		sig_pairs(i, 1) = tmp_sig_pairs(curr_pos, 1);
		sig_pairs(i, 2) = tmp_sig_pairs(curr_pos, 2);
		sig_pairs(i, 3) = tmp_sig_pairs(curr_pos, 3);
	}
	return sig_pairs;
}

static unsigned int link_vars(arma::mat& ds, arma::mat& cor_ds, arma::uvec& max_min, arma::mat& pvalues, 
		arma::mat& st, arma::mat& sig_pairs, arma::mat& ch, 
		const mth_t mth, const std::string method, const double sig_log, 
		const unsigned int k, const unsigned int r) {
	unsigned int ntests = 0;
	for (unsigned int i = 0; i < sig_pairs.n_rows; i++) {
		dbl_print("i = %u\n", i);

		const unsigned int spo_i0 = sig_pairs(i, 0);
		const unsigned int spo_i1 = sig_pairs(i, 1);
		dbl_print("spo_i0 = %u\n", spo_i0);
		dbl_print("spo_i1 = %u\n", spo_i1);

		dbl_print("Calculating xadj, yadj.\n");
		std::vector<unsigned int> xadj;
		std::vector<unsigned int> yadj;
		calc_adj(st, spo_i0, spo_i1,
				xadj, yadj);

		arma::mat xsam;
		arma::mat ysam;
		bool xsam_null = true;
		bool ysam_null = true;
		if (xadj.size() >= k) {
			dbl_print("x_size >= k\n");
			dbl_print("Calculating xsam.\n");
			calc_sam(pvalues, xadj, spo_i0, k, xsam);
			xsam_null = false;
		}
		if (yadj.size() >= k) {
			dbl_print("y_size >= k\n");
			dbl_print("Calculating ysam.\n");
			calc_sam(pvalues, yadj, spo_i1, k, ysam);
			ysam_null = false;
		}

		dbl_print("Calculating sam.\n");
		arma::mat sam = rbind_uniq(xsam, ysam, !xsam_null, !ysam_null);

		dbl_print("Calculating inter_res.\n");
		arma::vec lh(2); lh[0] = (double) spo_i0; lh[1] = (double) spo_i1;
		arma::vec rh = to_vec(sam);
		std::vector<double> inter_res = inter(lh, rh);

		if (inter_res.size()) {
			dbl_print("inter_res.size()\n");
			dbl_print("Calculating rem_rows.\n");
			arma::uvec rem_rows = arma::conv_to<arma::uvec>::from(index_row_eq(sam, inter_res));

			dbl_print("Removing rem_rows from sam.\n");
			sam = rm_rows(sam, rem_rows);
		}

		if (sam.n_rows > 1) {
			dbl_print("sam.n_rows > 1\n");
			sam = order_col(sam, 1);
		}
		if (sam.n_rows) {
			dbl_print("sam.n_rows\n");
			std::vector<double> pvs;
			arma::vec cat_condi = calc_cat_condi(ds, cor_ds, max_min, 
					sam, spo_i0, spo_i1, 0, k, mth.is_cat_method, method, r);
			pvs.push_back(cat_condi(1));
			if (cat_condi(1) > sig_log) {
				dbl_print("cat_condi[1] > sig_log\n");
				update_st_ch(st, ch, sam, cat_condi, spo_i0, spo_i1, i, 0, k);
				ntests++;
			}
			else {
				dbl_print("cat_condi[1] <= sig_log\n");
				unsigned int m = 0;
				while (cat_condi(1) < sig_log && m < (sam.n_rows - 1)) {
					m++;
					cat_condi = calc_cat_condi(ds, cor_ds, max_min, 
							sam, spo_i0, spo_i1, m, k, mth.is_cat_method, method, r);
					pvs.push_back(cat_condi(1));
					ntests++;
				}
				if (cat_condi(1) > sig_log) {
					dbl_print("cat_condi[1] > sig_log\n");
					update_st_ch(st, ch, sam, cat_condi, spo_i0, spo_i1, i, m, k);
				}
			}
			const double curr_pv = pvalues(spo_i0, spo_i1);
			pvs.push_back(curr_pv);
			const double pvs_max = *std::max_element(std::begin(pvs), std::end(pvs));
			pvalues(spo_i0, spo_i1) = pvs_max;
			pvalues(spo_i1, spo_i0) = pvs_max;
		}
	}
	return ntests;
}
	
static void calc_adj(arma::mat& st, const unsigned int xrow, const unsigned int yrow,
		std::vector<unsigned int>& xadj, std::vector<unsigned int>& yadj) {
	for (unsigned int j = 0; j < st.n_cols; j++) {
		if (st(xrow, j) == 2) {
			xadj.push_back(j);
		}
		if (st(yrow, j) == 2) {
			yadj.push_back(j);
		}
	}
}

static void calc_sam(arma::mat& pvalues, std::vector<unsigned int>& adj,
		const unsigned int spo_row, const unsigned int k, arma::mat& sam) {
	arma::mat info(adj.size(), 2);
	for (unsigned int j = 0; j < adj.size(); j++) {
		const unsigned int curr = adj.at(j);
		info(j, 0) = curr;
		info(j, 1) = pvalues(spo_row, curr);
	}
	arma::vec tmp_info = info.col(1);
	arma::uvec ind_order = arma::sort_index(-tmp_info);
	arma::mat info_ordered(adj.size(), 2);
	for (unsigned int j = 0; j < adj.size(); j++) {
		const unsigned int ord_row = ind_order(j);
		info_ordered(j, 0) = info(ord_row, 0);
		info_ordered(j, 1) = info(ord_row, 1);
	}
	if (info_ordered.n_rows == 1) {
		sam = info_ordered;
	}
	else {
		arma::vec tmp_info_ordered_0 = info_ordered.col(0);
		arma::vec tmp_info_ordered_1 = info_ordered.col(1);
		arma::mat combn_adj = find_combn(tmp_info_ordered_0, k);
		arma::mat combn_pvl = find_combn(tmp_info_ordered_1, k);
		sam = cbind_tran_mat(combn_adj, combn_pvl);
	}
}

static arma::vec calc_cat_condi(arma::mat& ds, arma::mat& cor_ds, arma::uvec& max_min, 
		arma::mat& sam, const unsigned int spo_i0, const unsigned int spo_i1, 
		const unsigned int m, const unsigned int k,
		const bool is_cat_method, const std::string method, const unsigned int r) {
	arma::uvec cols(k);
	std::iota(cols.begin(), cols.end(), 0);
	arma::uvec cs = form_vec(sam, m, cols);
	if (is_cat_method) {
		return cat_ci(spo_i0, spo_i1, cs, ds, max_min, r);
	}
	return calc_condi(spo_i0, spo_i1, cs, ds, cor_ds, method, r);
}

static void update_st_ch(arma::mat& st, arma::mat& ch, arma::mat& sam, 
		arma::vec& cat_condi, const unsigned int spo_i0, const unsigned int spo_i1, 
		const unsigned int curr_row, const unsigned int m, const unsigned int k) {
	st(spo_i0, spo_i1) = 0;
	st(spo_i1, spo_i0) = 0;
	arma::uvec cols(k);
	std::iota(cols.begin(), cols.end(), 0);
	arma::vec vals(2); vals[0] = cat_condi[0]; vals[1] = cat_condi[1];
	arma::vec tmp = form_vec_wvals(sam, m, cols, vals);
	append_row(ch, curr_row, tmp);
}

static unsigned int upd_state(arma::mat& st, arma::mat& sig_pairs, arma::mat& ch, 
		const unsigned int k, const unsigned int el, Rcpp::List& sep) {
	dbl_print("Updating el.\n");
	unsigned int max_len = el;
	for (unsigned int i = 0; i < sig_pairs.n_rows; i++) {
		std::vector<unsigned int> lens = loc_inactv_pos(sig_pairs, st, i);
		unsigned int curr_max_len = std::max(lens.at(0), lens.at(1));
		if (curr_max_len > max_len) {
			max_len = curr_max_len;
		}
	}
	dbl_print("el = %u\n", el);

	dbl_print("Forming sep.\n");
	std::vector<unsigned int> idxs = rsum_gt_zero_idxs(ch);
	std::ostringstream oss;
	oss << k;
	arma::mat res;
	if (idxs.size() == 1) {
		dbl_print("idxs.size() == 1\n");
		std::vector<unsigned int> lh_cols = {0, 1};
		arma::mat lh = form_rmat_std(sig_pairs, idxs, lh_cols);
		arma::rowvec rh = ch.row(idxs.at(0));
		res = form_cmat_vec(lh, rh);
	}
	else {
		dbl_print("idxs.size() != 1\n");
		std::vector<unsigned int> lh_cols = {0, 1};
		std::vector<unsigned int> rh_cols(ch.n_cols);
		std::iota(rh_cols.begin(), rh_cols.end(), 0);
		arma::mat lh = form_rmat_std(sig_pairs, idxs, lh_cols);
		arma::mat rh = form_rmat_std(ch, idxs, rh_cols);
		res = cbind_mat(lh, rh);
	}
	if (!res.is_empty()) {
		if (res.n_cols > 1) {
			res.cols(0, k + 1) += 1;
		}
		else {
			res.rows(0, k + 1) += 1;
		}
		if (res.n_cols == 1) {
			res.reshape(res.n_cols, res.n_rows);
		}
		sep[oss.str()] = res;
	}
	if (idxs.size()) {
		dbl_print("idxs.size()\n");
		dbl_print("Adjusting sig_pairs\n");
		sig_pairs = rm_rows_std(sig_pairs, idxs);
		if (sig_pairs.n_rows == 1) {
			dbl_print("sig_pairs.n_rows == 1\n");
			sig_pairs = adj_cols(sig_pairs, 4);
		}
	}
	else {
		dbl_print("!idxs.size()\n");
		sig_pairs.reset();
	}
	return max_len;
}

static std::vector<unsigned int> loc_inactv_pos(arma::mat& sig_pairs,
		arma::mat& st, const unsigned int row) {
	std::vector<unsigned int> xa;
	std::vector<unsigned int> ya;
	const unsigned int xrow = sig_pairs(row, 0);
	const unsigned int yrow = sig_pairs(row, 1);
	for (unsigned int j = 0; j < st.n_cols; j++) {
		if (st(xrow, j) == 2) {
			xa.push_back(j);
		}
		if (st(yrow, j) == 2) {
			ya.push_back(j);
		}
	}
	std::vector<unsigned int> lens;
	lens.push_back(xa.size());
	lens.push_back(ya.size());
	return lens;
}

