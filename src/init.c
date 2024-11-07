#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP Rfast_add_to_namespace(SEXP, SEXP, SEXP);
SEXP Rfast_as_integer(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_apply_condition(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_bcdcor(SEXP, SEXP);
SEXP Rfast_binarysearch(SEXP, SEXP);
SEXP Rfast_bic_fs_reg(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_bs_reg(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_bincomb(SEXP);
SEXP Rfast_coeff(SEXP, SEXP);
SEXP Rfast_coeff_vec(SEXP, SEXP);
SEXP Rfast_col_all(SEXP);
SEXP Rfast_col_count_values(SEXP, SEXP);
SEXP Rfast_col_cum_maxs(SEXP);
SEXP Rfast_col_cum_sums(SEXP);
SEXP Rfast_col_cum_mins(SEXP);
SEXP Rfast_col_cum_prods(SEXP);
SEXP Rfast_col_meds(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_min_indices(SEXP);
SEXP Rfast_col_min(SEXP, SEXP, SEXP);
SEXP Rfast_col_sums(SEXP, SEXP, SEXP);
SEXP Rfast_col_min_max(SEXP, SEXP, SEXP);
SEXP Rfast_col_max_indices(SEXP);
SEXP Rfast_col_max(SEXP, SEXP, SEXP);
SEXP Rfast_col_means(SEXP, SEXP, SEXP);
SEXP Rfast_col_nth(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_len_sort_un_int(SEXP);
SEXP Rfast_col_ranks(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_tabulate(SEXP, SEXP);
SEXP Rfast_col_shuffle(SEXP);
SEXP Rfast_col_true(SEXP);
SEXP Rfast_col_diffs(SEXP);
SEXP Rfast_col_false(SEXP);
SEXP Rfast_col_vars(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_order(SEXP, SEXP, SEXP);
SEXP Rfast_col_pmax(SEXP, SEXP);
SEXP Rfast_col_pmin(SEXP, SEXP);
SEXP Rfast_col_prods(SEXP, SEXP);
SEXP Rfast_columns(SEXP, SEXP);
SEXP Rfast_chi2Test(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_chi2Test_univariate(SEXP, SEXP);
SEXP Rfast_chi2tests(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_check_namespace(SEXP, SEXP, SEXP);
SEXP Rfast_check_aliases(SEXP, SEXP, SEXP);
SEXP Rfast_check_true_false(SEXP, SEXP);
SEXP Rfast_check_usage(SEXP, SEXP, SEXP);
SEXP Rfast_col_any(SEXP);
SEXP Rfast_col_anovas(SEXP, SEXP);
SEXP Rfast_cholesky(SEXP);
SEXP Rfast_cholesky_par(SEXP);
SEXP Rfast_col_mads(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_true_false(SEXP);
SEXP Rfast_count_value(SEXP, SEXP);
SEXP Rfast_Choose(SEXP, SEXP);
SEXP Rfast_dvar(SEXP);
SEXP Rfast_dcor(SEXP, SEXP);
SEXP Rfast_dcov(SEXP, SEXP);
SEXP Rfast_diag_matrix_fill_scalar(SEXP, SEXP);
SEXP Rfast_diag_matrix_fill_vec(SEXP, SEXP);
SEXP Rfast_diag_fill_scalar(SEXP, SEXP);
SEXP Rfast_diag_fill_vec(SEXP, SEXP);
SEXP Rfast_design_matrix(SEXP, SEXP);
SEXP Rfast_dist(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_dist_vec(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_Digamma(SEXP);
SEXP Rfast_design_matrix_big(SEXP);
SEXP Rfast_dista(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_diri_nr_type2(SEXP a1SEXP, SEXP a2SEXP, SEXP maSEXP, SEXP, SEXP);
SEXP Rfast_edist(SEXP, SEXP);
SEXP Rfast_eachcol_apply(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_eachrow(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_floyd_john(SEXP, SEXP);
SEXP Rfast_frame_to_matrix(SEXP);
SEXP Rfast_fs_reg(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_gaussian_nb(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_g2Test(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_g2Test_univariate_perm(SEXP, SEXP, SEXP);
SEXP Rfast_g2Test_univariate(SEXP, SEXP);
SEXP Rfast_g2Test_perm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_g2tests_perm(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_g2tests(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_group(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_group_sum(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_glm_logistic(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_glm_poisson(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_hash2list(SEXP, SEXP);
SEXP Rfast_hash_find(SEXP, SEXP);
SEXP Rfast_Hash_list(SEXP, SEXP);
SEXP Rfast_Hash_key_multi(SEXP, SEXP, SEXP);
SEXP Rfast_is_element(SEXP, SEXP);
SEXP Rfast_is_element_string(SEXP, SEXP);
SEXP Rfast_is_integer(SEXP);
SEXP Rfast_k_nn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_k_nn_cv(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_lowerbound(SEXP, SEXP);
SEXP Rfast_Lgamma(SEXP);
SEXP Rfast_len_sort_unique_int(SEXP);
SEXP Rfast_Log(SEXP);
SEXP Rfast_logistic_only(SEXP, SEXP, SEXP);
SEXP Rfast_logistic_only_b(SEXP, SEXP, SEXP);
SEXP Rfast_Lbeta(SEXP, SEXP);
SEXP Rfast_lower_tri(SEXP, SEXP);
SEXP Rfast_lower_tri_assign(SEXP, SEXP, SEXP);
SEXP Rfast_lower_tri_b(SEXP, SEXP, SEXP);
SEXP Rfast_Lchoose(SEXP, SEXP);
SEXP Rfast_mahaCpp(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_min_freq_d(SEXP, SEXP);
SEXP Rfast_min_freq_i(SEXP, SEXP);
SEXP Rfast_max_freq_d(SEXP, SEXP);
SEXP Rfast_max_freq_i(SEXP, SEXP);
SEXP Rfast_Match(SEXP, SEXP);
SEXP Rfast_mad2(SEXP, SEXP, SEXP);
SEXP Rfast_min_max(SEXP, SEXP);
SEXP Rfast_mat_mat(SEXP, SEXP);
SEXP Rfast_med(SEXP, SEXP);
SEXP Rfast_min_max_perc(SEXP);
SEXP Rfast_negative(SEXP, SEXP);
SEXP Rfast_nth(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_nth_int(SEXP, SEXP);
SEXP Rfast_Norm(SEXP, SEXP);
SEXP Rfast_Order(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_odds_helper(SEXP);
SEXP Rfast_Outer(SEXP, SEXP, SEXP);
SEXP Rfast_positive(SEXP, SEXP);
SEXP Rfast_positive_negative(SEXP, SEXP);
SEXP Rfast_poisson_only(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_partial_sort(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_partial_sort_index(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_prop_reg(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_prop_regs(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_pmin(SEXP, SEXP, SEXP);
SEXP Rfast_pmax(SEXP, SEXP, SEXP);
SEXP Rfast_pmin_pmax(SEXP, SEXP, SEXP);
SEXP Rfast_pc_skel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_poisson_only_b(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_permutation(SEXP, SEXP);
SEXP Rfast_permutation_next(SEXP, SEXP);
SEXP Rfast_permutation_prev(SEXP, SEXP);
SEXP Rfast_perm_cor(SEXP, SEXP, SEXP);
SEXP Rfast_qpois_reg(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_qpois_regs(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_rank(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_remove_from_namespace(SEXP, SEXP);
SEXP Rfast_rbing(SEXP, SEXP);
SEXP Rfast_rows(SEXP, SEXP);
SEXP Rfast_row_any(SEXP);
SEXP Rfast_row_means(SEXP);
SEXP Rfast_row_max(SEXP);
SEXP Rfast_row_meds(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_min(SEXP);
SEXP Rfast_row_len_sort_un_int(SEXP);
SEXP Rfast_row_ranks(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_sums(SEXP, SEXP, SEXP);
SEXP Rfast_rmdp(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_tabulate(SEXP, SEXP);
SEXP Rfast_rep_col(SEXP, SEXP);
SEXP Rfast_rep_row(SEXP, SEXP);
SEXP Rfast_row_nth(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_min_max(SEXP);
SEXP Rfast_row_shuffle(SEXP);
SEXP Rfast_Round(SEXP, SEXP, SEXP);
SEXP Rfast_rvmf(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_rvonmises(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_all(SEXP);
SEXP Rfast_row_true(SEXP);
SEXP Rfast_row_prods(SEXP);
SEXP Rfast_row_false(SEXP);
SEXP Rfast_row_order(SEXP, SEXP, SEXP);
SEXP Rfast_row_true_false(SEXP);
SEXP Rfast_read_examples(SEXP, SEXP);
SEXP Rfast_row_count_values(SEXP, SEXP);
SEXP Rfast_row_mads(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_vars(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_max_indices(SEXP);
SEXP Rfast_row_min_indices(SEXP);
SEXP Rfast_sort_unique_double(SEXP);
SEXP Rfast_sort_mat(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_submatrix(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_sum_XopY(SEXP, SEXP, SEXP);
SEXP Rfast_sum_XopX(SEXP, SEXP);
SEXP Rfast_sum_lower_tri(SEXP, SEXP);
SEXP Rfast_sum_upper_tri(SEXP, SEXP);
SEXP Rfast_sort_unique_int(SEXP);
SEXP Rfast_symmetric(SEXP);
SEXP Rfast_Sort(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_Sort_na_first(SEXP, SEXP, SEXP);
SEXP Rfast_sort_string(SEXP, SEXP, SEXP);
SEXP Rfast_sort_int(SEXP);
SEXP Rfast_stable_sort(SEXP, SEXP, SEXP);
SEXP Rfast_spat_med(SEXP, SEXP);
SEXP Rfast_squareform_c(SEXP);
SEXP Rfast_total_dists(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_total_dista(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_topological_sort(SEXP);
SEXP Rfast_Trigamma(SEXP);
SEXP Rfast_table_c(SEXP, SEXP);
SEXP Rfast_table_with_names(SEXP);
SEXP Rfast_table_sign(SEXP, SEXP, SEXP);
SEXP Rfast_table2_c(SEXP, SEXP, SEXP);
SEXP Rfast_table2_with_names(SEXP, SEXP, SEXP);
SEXP Rfast_transpose(SEXP);
SEXP Rfast_Unique(SEXP, SEXP);
SEXP Rfast_upper_tri(SEXP, SEXP);
SEXP Rfast_upper_tri_assign(SEXP, SEXP, SEXP);
SEXP Rfast_upper_tri_b(SEXP, SEXP, SEXP);
SEXP Rfast_var(SEXP, SEXP, SEXP);
SEXP Rfast_comb_n(SEXP, SEXP, SEXP);
SEXP Rfast_varcomps_mle(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_vecdist(SEXP);
SEXP Rfast_which_is(SEXP, SEXP);
SEXP Rfast_col_row_value(SEXP, SEXP);

SEXP Rfast_colrint_mle(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_eigs_sym_c(SEXP, SEXP, SEXP);
SEXP Rfast_geom_regs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_dir_knn(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_multinom_regs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_normlog_regs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_normlog_reg(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_rint_reg(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_rint_regs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_rint_mle(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_weibull_mle(SEXP, SEXP, SEXP);
SEXP Rfast_weib_reg(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_spml_mle(SEXP, SEXP, SEXP);
SEXP Rfast_spml_regs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_spml_reg(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_colweibull_mle(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_quasi_poisson_only(SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP Rfast_col_all_p(SEXP, SEXP);
SEXP Rfast_col_count_values_p(SEXP, SEXP, SEXP);
SEXP Rfast_col_nth_p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_order_p(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_col_sums_p(SEXP, SEXP);
SEXP Rfast_mat_mult_p(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_all_p(SEXP, SEXP);
SEXP Rfast_row_count_values_p(SEXP, SEXP, SEXP);
SEXP Rfast_row_nth_p(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_order_p(SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_ranks_p(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP Rfast_row_sums_p(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"Rfast_add_to_namespace", (DL_FUNC)&Rfast_add_to_namespace, 3},
    {"Rfast_apply_condition", (DL_FUNC)&Rfast_apply_condition, 4},
    {"Rfast_as_integer", (DL_FUNC)&Rfast_as_integer, 4},
    {"Rfast_bcdcor", (DL_FUNC)&Rfast_bcdcor, 2},
    {"Rfast_binarysearch", (DL_FUNC)&Rfast_binarysearch, 2},
    {"Rfast_bic_fs_reg", (DL_FUNC)&Rfast_bic_fs_reg, 4},
    {"Rfast_bs_reg", (DL_FUNC)&Rfast_bs_reg, 4},
    {"Rfast_bincomb", (DL_FUNC)&Rfast_bincomb, 1},
    {"Rfast_coeff", (DL_FUNC)&Rfast_coeff, 2},
    {"Rfast_coeff_vec", (DL_FUNC)&Rfast_coeff_vec, 2},
    {"Rfast_col_all", (DL_FUNC)&Rfast_col_all, 1},
    {"Rfast_col_count_values", (DL_FUNC)&Rfast_col_count_values, 2},
    {"Rfast_col_cum_maxs", (DL_FUNC)&Rfast_col_cum_maxs, 1},
    {"Rfast_col_cum_sums", (DL_FUNC)&Rfast_col_cum_sums, 1},
    {"Rfast_col_cum_mins", (DL_FUNC)&Rfast_col_cum_mins, 1},
    {"Rfast_col_cum_prods", (DL_FUNC)&Rfast_col_cum_prods, 1},
    {"Rfast_col_meds", (DL_FUNC)&Rfast_col_meds, 4},
    {"Rfast_col_min_indices", (DL_FUNC)&Rfast_col_min_indices, 1},
    {"Rfast_col_min", (DL_FUNC)&Rfast_col_min, 3},
    {"Rfast_col_sums", (DL_FUNC)&Rfast_col_sums, 3},
    {"Rfast_col_min_max", (DL_FUNC)&Rfast_col_min_max, 3},
    {"Rfast_col_max_indices", (DL_FUNC)&Rfast_col_max_indices, 1},
    {"Rfast_col_max", (DL_FUNC)&Rfast_col_max, 3},
    {"Rfast_col_means", (DL_FUNC)&Rfast_col_means, 3},
    {"Rfast_col_nth", (DL_FUNC)&Rfast_col_nth, 6},
    {"Rfast_col_len_sort_un_int", (DL_FUNC)&Rfast_col_len_sort_un_int, 1},
    {"Rfast_col_ranks", (DL_FUNC)&Rfast_col_ranks, 6},
    {"Rfast_col_tabulate", (DL_FUNC)&Rfast_col_tabulate, 2},
    {"Rfast_col_shuffle", (DL_FUNC)&Rfast_col_shuffle, 1},
    {"Rfast_col_true", (DL_FUNC)&Rfast_col_true, 1},
    {"Rfast_col_diffs", (DL_FUNC)&Rfast_col_diffs, 1},
    {"Rfast_col_false", (DL_FUNC)&Rfast_col_false, 1},
    {"Rfast_col_order", (DL_FUNC)&Rfast_col_order, 3},
    {"Rfast_col_pmax", (DL_FUNC)&Rfast_col_pmax, 2},
    {"Rfast_col_pmin", (DL_FUNC)&Rfast_col_pmin, 2},
    {"Rfast_col_prods", (DL_FUNC)&Rfast_col_prods, 2},
    {"Rfast_col_vars", (DL_FUNC)&Rfast_col_vars, 5},
    {"Rfast_columns", (DL_FUNC)&Rfast_columns, 2},
    {"Rfast_chi2Test_univariate", (DL_FUNC)&Rfast_chi2Test_univariate, 2},
    {"Rfast_chi2Test", (DL_FUNC)&Rfast_chi2Test, 5},
    {"Rfast_chi2tests", (DL_FUNC)&Rfast_chi2tests, 4},
    {"Rfast_check_namespace", (DL_FUNC)&Rfast_check_namespace, 3},
    {"Rfast_check_aliases", (DL_FUNC)&Rfast_check_aliases, 3},
    {"Rfast_check_true_false", (DL_FUNC)&Rfast_check_true_false, 2},
    {"Rfast_check_usage", (DL_FUNC)&Rfast_check_usage, 3},
    {"Rfast_col_any", (DL_FUNC)&Rfast_col_any, 1},
    {"Rfast_col_anovas", (DL_FUNC)&Rfast_col_anovas, 2},
    {"Rfast_cholesky", (DL_FUNC)&Rfast_cholesky, 1},
    {"Rfast_cholesky_par", (DL_FUNC)&Rfast_cholesky_par, 1},
    {"Rfast_col_mads", (DL_FUNC)&Rfast_col_mads, 5},
    {"Rfast_col_true_false", (DL_FUNC)&Rfast_col_true_false, 1},
    {"Rfast_count_value", (DL_FUNC)&Rfast_count_value, 2},
    {"Rfast_Choose", (DL_FUNC)&Rfast_Choose, 2},
    {"Rfast_dvar", (DL_FUNC)&Rfast_dvar, 1},
    {"Rfast_dcor", (DL_FUNC)&Rfast_dcor, 2},
    {"Rfast_dcov", (DL_FUNC)&Rfast_dcov, 2},
    {"Rfast_diag_matrix_fill_scalar", (DL_FUNC)&Rfast_diag_matrix_fill_scalar, 2},
    {"Rfast_diag_matrix_fill_vec", (DL_FUNC)&Rfast_diag_matrix_fill_vec, 2},
    {"Rfast_diag_fill_scalar", (DL_FUNC)&Rfast_diag_fill_scalar, 2},
    {"Rfast_diag_fill_vec", (DL_FUNC)&Rfast_diag_fill_vec, 2},
    {"Rfast_design_matrix", (DL_FUNC)&Rfast_design_matrix, 2},
    {"Rfast_dist", (DL_FUNC)&Rfast_dist, 5},
    {"Rfast_dist_vec", (DL_FUNC)&Rfast_dist_vec, 5},
    {"Rfast_Digamma", (DL_FUNC)&Rfast_Digamma, 1},
    {"Rfast_design_matrix_big", (DL_FUNC)&Rfast_design_matrix_big, 1},
    {"Rfast_dista", (DL_FUNC)&Rfast_dista, 8},
    {"Rfast_diri_nr_type2", (DL_FUNC)&Rfast_diri_nr_type2, 5},
    {"Rfast_edist", (DL_FUNC)&Rfast_edist, 2},
    {"Rfast_eachrow", (DL_FUNC)&Rfast_eachrow, 4},
    {"Rfast_eachcol_apply", (DL_FUNC)&Rfast_eachcol_apply, 6},
    {"Rfast_floyd_john", (DL_FUNC)&Rfast_floyd_john, 2},
    {"Rfast_frame_to_matrix", (DL_FUNC)&Rfast_frame_to_matrix, 1},
    {"Rfast_fs_reg", (DL_FUNC)&Rfast_fs_reg, 5},
    {"Rfast_gaussian_nb", (DL_FUNC)&Rfast_gaussian_nb, 7},
    {"Rfast_g2Test_univariate_perm", (DL_FUNC)&Rfast_g2Test_univariate_perm, 3},
    {"Rfast_g2Test_univariate", (DL_FUNC)&Rfast_g2Test_univariate, 2},
    {"Rfast_g2Test_perm", (DL_FUNC)&Rfast_g2Test_perm, 6},
    {"Rfast_g2Test", (DL_FUNC)&Rfast_g2Test, 5},
    {"Rfast_g2tests_perm", (DL_FUNC)&Rfast_g2tests_perm, 5},
    {"Rfast_g2tests", (DL_FUNC)&Rfast_g2tests, 4},
    {"Rfast_group", (DL_FUNC)&Rfast_group, 6},
    {"Rfast_group_sum", (DL_FUNC)&Rfast_group_sum, 4},
    {"Rfast_glm_logistic", (DL_FUNC)&Rfast_glm_logistic, 4},
    {"Rfast_glm_poisson", (DL_FUNC)&Rfast_glm_poisson, 4},
    {"Rfast_hash2list", (DL_FUNC)&Rfast_hash2list, 2},
    {"Rfast_hash_find", (DL_FUNC)&Rfast_hash_find, 2},
    {"Rfast_Hash_list", (DL_FUNC)&Rfast_Hash_list, 2},
    {"Rfast_Hash_key_multi", (DL_FUNC)&Rfast_Hash_key_multi, 3},
    {"Rfast_is_element", (DL_FUNC)&Rfast_is_element, 2},
    {"Rfast_is_element_string", (DL_FUNC)&Rfast_is_element_string, 2},
    {"Rfast_is_integer", (DL_FUNC)&Rfast_is_integer, 1},
    {"Rfast_k_nn", (DL_FUNC)&Rfast_k_nn, 9},
    {"Rfast_k_nn_cv", (DL_FUNC)&Rfast_k_nn_cv, 10},
    {"Rfast_lowerbound", (DL_FUNC)&Rfast_lowerbound, 2},
    {"Rfast_Lgamma", (DL_FUNC)&Rfast_Lgamma, 1},
    {"Rfast_len_sort_unique_int", (DL_FUNC)&Rfast_len_sort_unique_int, 1},
    {"Rfast_Log", (DL_FUNC)&Rfast_Log, 1},
    {"Rfast_logistic_only_b", (DL_FUNC)&Rfast_logistic_only_b, 3},
    {"Rfast_logistic_only", (DL_FUNC)&Rfast_logistic_only, 3},
    {"Rfast_Lbeta", (DL_FUNC)&Rfast_Lbeta, 2},
    {"Rfast_lower_tri", (DL_FUNC)&Rfast_lower_tri, 2},
    {"Rfast_lower_tri_assign", (DL_FUNC)&Rfast_lower_tri_assign, 3},
    {"Rfast_lower_tri_b", (DL_FUNC)&Rfast_lower_tri_b, 3},
    {"Rfast_Lchoose", (DL_FUNC)&Rfast_Lchoose, 2},
    {"Rfast_mahaCpp", (DL_FUNC)&Rfast_mahaCpp, 4},
    {"Rfast_min_freq_d", (DL_FUNC)&Rfast_min_freq_d, 2},
    {"Rfast_min_freq_i", (DL_FUNC)&Rfast_min_freq_i, 2},
    {"Rfast_max_freq_d", (DL_FUNC)&Rfast_max_freq_d, 2},
    {"Rfast_max_freq_i", (DL_FUNC)&Rfast_max_freq_i, 2},
    {"Rfast_Match", (DL_FUNC)&Rfast_Match, 2},
    {"Rfast_mad2", (DL_FUNC)&Rfast_mad2, 3},
    {"Rfast_min_max", (DL_FUNC)&Rfast_min_max, 2},
    {"Rfast_mat_mat", (DL_FUNC)&Rfast_mat_mat, 2},
    {"Rfast_med", (DL_FUNC)&Rfast_med, 2},
    {"Rfast_min_max_perc", (DL_FUNC)&Rfast_min_max_perc, 1},
    {"Rfast_negative", (DL_FUNC)&Rfast_negative, 2},
    {"Rfast_nth", (DL_FUNC)&Rfast_nth, 6},
    {"Rfast_nth_int", (DL_FUNC)&Rfast_nth_int, 2},
    {"Rfast_Norm", (DL_FUNC)&Rfast_Norm, 2},
    {"Rfast_Order", (DL_FUNC)&Rfast_Order, 4},
    {"Rfast_Outer", (DL_FUNC)&Rfast_Outer, 3},
    {"Rfast_odds_helper", (DL_FUNC)&Rfast_odds_helper, 1},
    {"Rfast_positive", (DL_FUNC)&Rfast_positive, 2},
    {"Rfast_positive_negative", (DL_FUNC)&Rfast_positive_negative, 2},
    {"Rfast_poisson_only", (DL_FUNC)&Rfast_poisson_only, 4},
    {"Rfast_partial_sort", (DL_FUNC)&Rfast_partial_sort, 4},
    {"Rfast_partial_sort_index", (DL_FUNC)&Rfast_partial_sort_index, 4},
    {"Rfast_prop_reg", (DL_FUNC)&Rfast_prop_reg, 4},
    {"Rfast_prop_regs", (DL_FUNC)&Rfast_prop_regs, 5},
    {"Rfast_pmin", (DL_FUNC)&Rfast_pmin, 3},
    {"Rfast_pmax", (DL_FUNC)&Rfast_pmax, 3},
    {"Rfast_pmin_pmax", (DL_FUNC)&Rfast_pmin_pmax, 3},
    {"Rfast_pc_skel", (DL_FUNC)&Rfast_pc_skel, 7},
    {"Rfast_poisson_only_b", (DL_FUNC)&Rfast_poisson_only_b, 4},
    {"Rfast_permutation", (DL_FUNC)&Rfast_permutation, 2},
    {"Rfast_permutation_next", (DL_FUNC)&Rfast_permutation_next, 2},
    {"Rfast_permutation_prev", (DL_FUNC)&Rfast_permutation_prev, 2},
    {"Rfast_perm_cor", (DL_FUNC)&Rfast_perm_cor, 3},
    {"Rfast_qpois_reg", (DL_FUNC)&Rfast_qpois_reg, 5},
    {"Rfast_qpois_regs", (DL_FUNC)&Rfast_qpois_regs, 5},
    {"Rfast_rank", (DL_FUNC)&Rfast_rank, 5},
    {"Rfast_remove_from_namespace", (DL_FUNC)&Rfast_remove_from_namespace, 2},
    {"Rfast_rbing", (DL_FUNC)&Rfast_rbing, 2},
    {"Rfast_rows", (DL_FUNC)&Rfast_rows, 2},
    {"Rfast_row_count_values", (DL_FUNC)&Rfast_row_count_values, 2},
    {"Rfast_row_any", (DL_FUNC)&Rfast_row_any, 1},
    {"Rfast_row_means", (DL_FUNC)&Rfast_row_means, 1},
    {"Rfast_row_max", (DL_FUNC)&Rfast_row_max, 1},
    {"Rfast_row_meds", (DL_FUNC)&Rfast_row_meds, 4},
    {"Rfast_row_min", (DL_FUNC)&Rfast_row_min, 1},
    {"Rfast_row_len_sort_un_int", (DL_FUNC)&Rfast_row_len_sort_un_int, 1},
    {"Rfast_row_ranks", (DL_FUNC)&Rfast_row_ranks, 4},
    {"Rfast_row_sums", (DL_FUNC)&Rfast_row_sums, 3},
    {"Rfast_rmdp", (DL_FUNC)&Rfast_rmdp, 5},
    {"Rfast_row_tabulate", (DL_FUNC)&Rfast_row_tabulate, 2},
    {"Rfast_rep_col", (DL_FUNC)&Rfast_rep_col, 2},
    {"Rfast_rep_row", (DL_FUNC)&Rfast_rep_row, 2},
    {"Rfast_row_nth", (DL_FUNC)&Rfast_row_nth, 6},
    {"Rfast_row_min_max", (DL_FUNC)&Rfast_row_min_max, 1},
    {"Rfast_row_shuffle", (DL_FUNC)&Rfast_row_shuffle, 1},
    {"Rfast_Round", (DL_FUNC)&Rfast_Round, 3},
    {"Rfast_rvmf", (DL_FUNC)&Rfast_rvmf, 4},
    {"Rfast_rvonmises", (DL_FUNC)&Rfast_rvonmises, 4},
    {"Rfast_row_all", (DL_FUNC)&Rfast_row_all, 1},
    {"Rfast_row_true", (DL_FUNC)&Rfast_row_true, 1},
    {"Rfast_row_prods", (DL_FUNC)&Rfast_row_prods, 1},
    {"Rfast_row_false", (DL_FUNC)&Rfast_row_false, 1},
    {"Rfast_row_order", (DL_FUNC)&Rfast_row_order, 3},
    {"Rfast_row_true_false", (DL_FUNC)&Rfast_row_true_false, 1},
    {"Rfast_read_examples", (DL_FUNC)&Rfast_read_examples, 2},
    {"Rfast_row_mads", (DL_FUNC)&Rfast_row_mads, 5},
    {"Rfast_row_vars", (DL_FUNC)&Rfast_row_vars, 5},
    {"Rfast_row_max_indices", (DL_FUNC)&Rfast_row_max_indices, 1},
    {"Rfast_row_min_indices", (DL_FUNC)&Rfast_row_min_indices, 1},
    {"Rfast_sort_unique_double", (DL_FUNC)&Rfast_sort_unique_double, 1},
    {"Rfast_sort_mat", (DL_FUNC)&Rfast_sort_mat, 6},
    {"Rfast_submatrix", (DL_FUNC)&Rfast_submatrix, 5},
    {"Rfast_sum_XopY", (DL_FUNC)&Rfast_sum_XopY, 3},
    {"Rfast_sum_XopX", (DL_FUNC)&Rfast_sum_XopX, 2},
    {"Rfast_sum_lower_tri", (DL_FUNC)&Rfast_sum_lower_tri, 2},
    {"Rfast_sum_upper_tri", (DL_FUNC)&Rfast_sum_upper_tri, 2},
    {"Rfast_sort_unique_int", (DL_FUNC)&Rfast_sort_unique_int, 1},
    {"Rfast_symmetric", (DL_FUNC)&Rfast_symmetric, 1},
    {"Rfast_Sort", (DL_FUNC)&Rfast_Sort, 4},
    {"Rfast_Sort_na_first", (DL_FUNC)&Rfast_Sort_na_first, 3},
    {"Rfast_sort_string", (DL_FUNC)&Rfast_sort_string, 3},
    {"Rfast_stable_sort", (DL_FUNC)&Rfast_stable_sort, 3},
    {"Rfast_sort_int", (DL_FUNC)&Rfast_sort_int, 1},
    {"Rfast_spat_med", (DL_FUNC)&Rfast_spat_med, 2},
    {"Rfast_squareform_c", (DL_FUNC)&Rfast_squareform_c, 1},
    {"Rfast_total_dists", (DL_FUNC)&Rfast_total_dists, 5},
    {"Rfast_total_dista", (DL_FUNC)&Rfast_total_dista, 7},
    {"Rfast_topological_sort", (DL_FUNC)&Rfast_topological_sort, 1},
    {"Rfast_Trigamma", (DL_FUNC)&Rfast_Trigamma, 1},
    {"Rfast_table_c", (DL_FUNC)&Rfast_table_c, 2},
    {"Rfast_table_sign", (DL_FUNC)&Rfast_table_sign, 3},
    {"Rfast_table_with_names", (DL_FUNC)&Rfast_table_with_names, 1},
    {"Rfast_table2_c", (DL_FUNC)&Rfast_table2_c, 3},
    {"Rfast_table2_with_names", (DL_FUNC)&Rfast_table2_with_names, 3},
    {"Rfast_transpose", (DL_FUNC)&Rfast_transpose, 1},
    {"Rfast_Unique", (DL_FUNC)&Rfast_Unique, 2},
    {"Rfast_upper_tri", (DL_FUNC)&Rfast_upper_tri, 2},
    {"Rfast_upper_tri_assign", (DL_FUNC)&Rfast_upper_tri_assign, 3},
    {"Rfast_upper_tri_b", (DL_FUNC)&Rfast_upper_tri_b, 3},
    {"Rfast_var", (DL_FUNC)&Rfast_var, 3},
    {"Rfast_varcomps_mle", (DL_FUNC)&Rfast_varcomps_mle, 4},
    {"Rfast_vecdist", (DL_FUNC)&Rfast_vecdist, 1},
    {"Rfast_comb_n", (DL_FUNC)&Rfast_comb_n, 3},
    {"Rfast_which_is", (DL_FUNC)&Rfast_which_is, 2},
    {"Rfast_col_row_value", (DL_FUNC)&Rfast_col_row_value, 2},

    {"Rfast_col_all_p", (DL_FUNC)&Rfast_col_all_p, 2},
    {"Rfast_col_count_values_p", (DL_FUNC)&Rfast_col_count_values_p, 3},
    {"Rfast_col_nth_p", (DL_FUNC)&Rfast_col_nth_p, 6},
    {"Rfast_col_sums_p", (DL_FUNC)&Rfast_col_sums_p, 2},
    {"Rfast_col_order_p", (DL_FUNC)&Rfast_col_order_p, 4},
    {"Rfast_mat_mult_p", (DL_FUNC)&Rfast_mat_mult_p, 4},
    {"Rfast_row_all_p", (DL_FUNC)&Rfast_row_all_p, 2},
    {"Rfast_row_count_values_p", (DL_FUNC)&Rfast_row_count_values_p, 3},
    {"Rfast_row_nth_p", (DL_FUNC)&Rfast_row_nth_p, 6},
    {"Rfast_row_order_p", (DL_FUNC)&Rfast_row_order_p, 4},
    {"Rfast_row_ranks_p", (DL_FUNC)&Rfast_row_ranks_p, 5},
    {"Rfast_row_sums_p", (DL_FUNC)&Rfast_row_sums_p, 2},

    {"Rfast_colrint_mle", (DL_FUNC)&Rfast_colrint_mle, 6},
    {"Rfast_eigs_sym_c", (DL_FUNC)&Rfast_eigs_sym_c, 3},
    {"Rfast_geom_regs", (DL_FUNC)&Rfast_geom_regs, 7},
    {"Rfast_normlog_regs", (DL_FUNC)&Rfast_normlog_regs, 8},
    {"Rfast_dir_knn", (DL_FUNC)&Rfast_dir_knn, 6},
    {"Rfast_multinom_regs", (DL_FUNC)&Rfast_multinom_regs, 6},
    {"Rfast_normlog_reg", (DL_FUNC)&Rfast_normlog_reg, 4},
    {"Rfast_rint_reg", (DL_FUNC)&Rfast_rint_reg, 6},
    {"Rfast_rint_regs", (DL_FUNC)&Rfast_rint_regs, 7},
    {"Rfast_rint_mle", (DL_FUNC)&Rfast_rint_mle, 5},
    {"Rfast_weibull_mle", (DL_FUNC)&Rfast_weibull_mle, 3},
    {"Rfast_weib_reg", (DL_FUNC)&Rfast_weib_reg, 4},
    {"Rfast_spml_mle", (DL_FUNC)&Rfast_spml_mle, 3},
    {"Rfast_spml_regs", (DL_FUNC)&Rfast_spml_regs, 6},
    {"Rfast_spml_reg", (DL_FUNC)&Rfast_spml_reg, 5},
    {"Rfast_colweibull_mle", (DL_FUNC)&Rfast_colweibull_mle, 4},
    {"Rfast_quasi_poisson_only", (DL_FUNC)&Rfast_quasi_poisson_only, 5},
    {NULL, NULL, 0}};

void R_init_Rfast(DllInfo *info)
{
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}
