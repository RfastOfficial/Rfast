#include <Rinternals.h>
#include <R_ext/Rdynload.h>



SEXP Rfast_add_to_namespace(SEXP dir_to_exportSEXP,SEXP dir_to_fileSEXP);
SEXP Rfast_as_integer(SEXP x,SEXP sortedSEXP,SEXP initSEXP);
SEXP Rfast_apply_condition(SEXP x,SEXP methodSEXP,SEXP operSEXP,SEXP valSEXP);
SEXP Rfast_bcdcor(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast_binarysearch(SEXP xSEXP,SEXP vSEXP);
SEXP Rfast_bic_fs_reg(SEXP ySEXP,SEXP dsSEXP,SEXP tolSEXP,SEXP typeSEXP);
SEXP Rfast_bs_reg(SEXP ySEXP,SEXP dsSEXP,SEXP sigSEXP,SEXP typeSEXP);
SEXP Rfast_bincomb(SEXP xSEXP);
SEXP Rfast_col_all(SEXP x);
SEXP Rfast_col_count_values(SEXP xSEXP,SEXP valuesSEXP);
SEXP Rfast_col_cum_maxs(SEXP x);
SEXP Rfast_col_cum_sums(SEXP x);
SEXP Rfast_col_cum_mins(SEXP x);
SEXP Rfast_col_cum_prods(SEXP x);
SEXP Rfast_col_meds(SEXP xSEXP,SEXP na_rmSEXP,SEXP parallelSEXP);
SEXP Rfast_col_min_indices(SEXP xSEXP);
SEXP Rfast_col_min(SEXP x);
SEXP Rfast_col_sums(SEXP xSEXP,SEXP SEXPindices);
SEXP Rfast_col_min_max(SEXP x);
SEXP Rfast_col_max_indices(SEXP xSEXP);
SEXP Rfast_col_max(SEXP x);
SEXP Rfast_col_means(SEXP xSEXP);
SEXP Rfast_col_nth(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP);
SEXP Rfast_col_len_sort_un_int(SEXP xSEXP);
SEXP Rfast_col_ranks(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP);
SEXP Rfast_col_tabulate(SEXP xSEXP,SEXP nrowwSEXP);
SEXP Rfast_col_shuffle(SEXP xSEXP);
SEXP Rfast_col_true(SEXP xSEXP);
SEXP Rfast_col_diffs(SEXP x);
SEXP Rfast_col_false(SEXP xSEXP);
SEXP Rfast_col_vars(SEXP xSEXP,SEXP parallelSEXP);
SEXP Rfast_col_order(SEXP xSEXP,SEXP stableSEXP,SEXP descendingSEXP);
SEXP Rfast_col_pmax(SEXP x,SEXP y);
SEXP Rfast_col_pmin(SEXP x,SEXP y);
SEXP Rfast_col_prods(SEXP xSEXP,SEXP methodSEXP);
SEXP Rfast_columns(SEXP x,SEXP ind);
SEXP Rfast_chi2Test(SEXP dataSEXP,SEXP xSEXP,SEXP ySEXP,SEXP csSEXP,SEXP dcSEXP);
SEXP Rfast_chi2Test_univariate(SEXP dataSEXP,SEXP dcSEXP);
SEXP Rfast_chi2tests(SEXP dataSEXP,SEXP xSEXP,SEXP ySEXP,SEXP dcSEXP);
SEXP Rfast_check_namespace(SEXP dir_to_exportSEXP,SEXP dir_to_fileSEXP);
SEXP Rfast_check_aliases(SEXP dir_to_manSEXP,SEXP dir_to_fileSEXP);
SEXP Rfast_check_true_false(SEXP path_manSEXP);
SEXP Rfast_check_usage(SEXP path_manSEXP,SEXP path_rfSEXP);
SEXP Rfast_col_any(SEXP x);
SEXP Rfast_col_anovas(SEXP ySEXP,SEXP xSEXP);
SEXP Rfast_cholesky(SEXP x);
SEXP Rfast_cholesky_par(SEXP x);
SEXP Rfast_col_mads(SEXP xSEXP,SEXP methodSEXP,SEXP na_rmSEXP,SEXP parallelSEXP);
SEXP Rfast_col_true_false(SEXP xSEXP);
SEXP Rfast_count_value(SEXP xSEXP,SEXP valueSEXP);
SEXP Rfast_Choose(SEXP x,SEXP kSEXP);
SEXP Rfast_dvar(SEXP xSEXP);
SEXP Rfast_dcor(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast_dcov(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast_diag_matrix_fill_scalar(SEXP xSEXP,SEXP vSEXP); 
SEXP Rfast_diag_matrix_fill_vec(SEXP lenSEXP,SEXP vSEXP);
SEXP Rfast_diag_fill_scalar(SEXP lenSEXP,SEXP vSEXP);
SEXP Rfast_diag_fill_vec(SEXP lenSEXP,SEXP vSEXP);
SEXP Rfast_design_matrix(SEXP xSEXP,SEXP onesSEXP);
SEXP Rfast_dist(SEXP xSEXP,SEXP methodSEXP,SEXP sqrSEXP,SEXP pSEXP);
SEXP Rfast_dist_vec(SEXP xSEXP,SEXP methodSEXP,SEXP sqrSEXP,SEXP pSEXP);
SEXP Rfast_Digamma(SEXP x);
SEXP Rfast_design_matrix_big(SEXP xSEXP);
SEXP Rfast_dista(SEXP XnewSEXP,SEXP XSEXP,SEXP sqrSEXP,SEXP typeSEXP);
SEXP Rfast_dista_index(SEXP XnewSEXP,SEXP XSEXP,SEXP kSEXP,SEXP typeSEXP);
SEXP Rfast_dista_values(SEXP XnewSEXP,SEXP XSEXP,SEXP kSEXP,SEXP sqrSEXP,SEXP typeSEXP);
SEXP Rfast_diri_nr_type2(SEXP a1SEXP,SEXP a2SEXP,SEXP maSEXP,SEXP pSEXP,SEXP tolSEXP);
SEXP Rfast_edist(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast_eachcol_apply(SEXP xSEXP,SEXP ySEXP,SEXP ind,SEXP operSEXP,SEXP methodSEXP);
SEXP Rfast_eachrow(SEXP x,SEXP y,SEXP operSEXP,SEXP method);
SEXP Rfast_floyd_john(SEXP nSEXP,SEXP xSEXP);
SEXP Rfast_frame_to_matrix(SEXP xSEXP);
SEXP Rfast_fs_reg(SEXP ySEXP,SEXP dsSEXP,SEXP sigSEXP,SEXP tolSEXP,SEXP typeSEXP);
SEXP Rfast_g2Test(SEXP dataSEXP,SEXP xSEXP,SEXP ySEXP,SEXP csSEXP,SEXP dcSEXP);
SEXP Rfast_g2Test_univariate_perm(SEXP dataSEXP,SEXP dcSEXP,SEXP npermSEXP);
SEXP Rfast_g2Test_univariate(SEXP dataSEXP,SEXP dcSEXP);
SEXP Rfast_g2Test_perm(SEXP dataSEXP,SEXP xSEXP,SEXP ySEXP,SEXP csSEXP,SEXP dcSEXP,SEXP npermSEXP);
SEXP Rfast_g2tests_perm(SEXP dataSEXP,SEXP xSEXP,SEXP ySEXP,SEXP dcSEXP,SEXP npermSEXP);
SEXP Rfast_g2tests(SEXP dataSEXP,SEXP xSEXP,SEXP ySEXP,SEXP dcSEXP);
SEXP Rfast_group(SEXP xSEXP,SEXP methodSEXP,SEXP groupSEXP,SEXP method_varSEXP,SEXP minSEXP,SEXP maxSEXP);
SEXP Rfast_group_sum(SEXP xSEXP,SEXP groupSEXP,SEXP minn,SEXP maxx);
SEXP Rfast_group_all(SEXP xSEXP,SEXP groupSEXP,SEXP nSEXP);
SEXP Rfast_group_any(SEXP xSEXP,SEXP groupSEXP,SEXP nSEXP);
SEXP Rfast_group_mad(SEXP xSEXP,SEXP groupSEXP,SEXP methodSEXP);
SEXP Rfast_group_max(SEXP xSEXP,SEXP groupSEXP,SEXP minn,SEXP maxx);
SEXP Rfast_group_mean(SEXP xSEXP,SEXP groupSEXP,SEXP max_nSEXP);
SEXP Rfast_group_med(SEXP xSEXP,SEXP groupSEXP,SEXP length_unique);
SEXP Rfast_group_min(SEXP xSEXP,SEXP groupSEXP,SEXP max_nSEXP);
SEXP Rfast_group_min_max(SEXP xSEXP,SEXP groupSEXP,SEXP max_nSEXP);
SEXP Rfast_group_var(SEXP xSEXP,SEXP groupSEXP,SEXP nSEXP);
SEXP Rfast_glm_logistic(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast_glm_poisson(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP);
SEXP Rfast_hash2list(SEXP xSEXP,SEXP sortingSEXP);
SEXP Rfast_hash_find(SEXP xSEXP,SEXP valueSEXP);
SEXP Rfast_Hash_list(SEXP keySEXP,SEXP xSEXP);
SEXP Rfast_Hash_key_multi(SEXP xSEXP,SEXP valueSEXP,SEXP sepSEXP);
SEXP Rfast_is_element(SEXP xSEXP,SEXP elSEXP);
SEXP Rfast_is_element_string(SEXP xSEXP,SEXP elSEXP);
SEXP Rfast_is_integer(SEXP xSEXP);
SEXP Rfast_k_nn(SEXP ds_extraSEXP,SEXP ySEXP,SEXP dsSEXP,SEXP idxsSEXP,SEXP dist_typeSEXP,SEXP typeSEXP,SEXP methodSEXP,SEXP freq_optionSEXP,SEXP mem_eff);
SEXP Rfast_k_nn_cv(SEXP foldsSEXP,SEXP ySEXP,SEXP dsSEXP,SEXP idxsSEXP,SEXP dist_typeSEXP,SEXP typeSEXP,SEXP methodSEXP,SEXP freq_optionSEXP,SEXP pred_retSEXP,SEXP mem_eff);
SEXP Rfast_lowerbound(SEXP xSEXP,SEXP vSEXP);
SEXP Rfast_Lgamma(SEXP x);
SEXP Rfast_len_sort_unique_int(SEXP xSEXP);
SEXP Rfast_Log(SEXP x);
SEXP Rfast_logistic_only(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP);
SEXP Rfast_logistic_only_b(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP);
SEXP Rfast_Lbeta(SEXP x,SEXP y);
SEXP Rfast_lower_tri(SEXP xSEXP,SEXP dgSEXP);
SEXP Rfast_lower_tri_assign(SEXP xSEXP,SEXP vSEXP,SEXP dgSEXP);
SEXP Rfast_lower_tri_b(SEXP nclSEXP, SEXP nrwSEXP,SEXP dgSEXP);
SEXP Rfast_Lchoose(SEXP x,SEXP kSEXP);
SEXP Rfast_mahaCpp(SEXP X, SEXP mu, SEXP sigma, SEXP isChol);
SEXP Rfast_min_freq_d(SEXP xSEXP,SEXP na_rmSEXP);
SEXP Rfast_min_freq_i(SEXP xSEXP,SEXP na_rmSEXP);
SEXP Rfast_max_freq_d(SEXP xSEXP,SEXP na_rmSEXP); 
SEXP Rfast_max_freq_i(SEXP xSEXP,SEXP na_rmSEXP);
SEXP Rfast_Match(SEXP xSEXP,SEXP keySEXP);
SEXP Rfast_mad2(SEXP xSEXP,SEXP methodSEXP,SEXP na_rmSEXP);
SEXP Rfast_min_max(SEXP x,SEXP indexSEXP);
SEXP Rfast_mat_mat(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast_med(SEXP x,SEXP na_rmSEXP);
SEXP Rfast_min_max_perc(SEXP x);
SEXP Rfast_negative(SEXP xSEXP,SEXP methodSEXP);
SEXP Rfast_nth(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP);
SEXP Rfast_nth_int(SEXP xSEXP,SEXP elemSEXP);
SEXP Rfast_Norm(SEXP xSEXP,SEXP typeSEXP);
SEXP Rfast_Order(SEXP xSEXP,SEXP stableSEXP,SEXP descendSEXP);
SEXP Rfast_odds_helper(SEXP x);
SEXP Rfast_Outer(SEXP x,SEXP y,SEXP operSEXP);
SEXP Rfast_positive(SEXP xSEXP,SEXP methodSEXP);
SEXP Rfast_positive_negative(SEXP xSEXP,SEXP methodSEXP);
SEXP Rfast_poisson_only(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP);
SEXP Rfast_partial_sort(SEXP x,SEXP nSEXP,SEXP descendSEXP);
SEXP Rfast_partial_sort_index(SEXP x,SEXP nSEXP,SEXP descendSEXP);
SEXP Rfast_prop_reg(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast_prop_regs(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP varbSEXP,SEXP maxitersSEXP);
SEXP Rfast_pmin(SEXP x,SEXP y,SEXP na_rmSEXP);
SEXP Rfast_pmax(SEXP x,SEXP y,SEXP na_rmSEXP);
SEXP Rfast_pmin_pmax(SEXP x,SEXP y,SEXP na_rmSEXP);
SEXP Rfast_pc_skel(SEXP dsSEXP,SEXP methodSEXP,SEXP sigSEXP,SEXP rSEXP,SEXP stats_initSEXP,SEXP pvalues_initSEXP,SEXP is_init_valsSEXP);
SEXP Rfast_poisson_only_b(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP);
SEXP Rfast_permutation(SEXP xSEXP,SEXP npermSEXP);
SEXP Rfast_permutation_next(SEXP xSEXP,SEXP npermSEXP);
SEXP Rfast_permutation_prev(SEXP xSEXP,SEXP npermSEXP);
SEXP Rfast_perm_cor(SEXP xSEXP,SEXP ySEXP,SEXP rSEXP);
SEXP Rfast_qpois_reg(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast_qpois_regs(SEXP xSEXP,SEXP ySEXP,SEXP tolSEXP,SEXP ylogySEXP,SEXP mySEXP);
SEXP Rfast_rank(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP);
SEXP Rfast_remove_from_namespace(SEXP dir_to_exportSEXP,SEXP files_to_removeSEXP);
SEXP Rfast_rbing(SEXP nSEXP,SEXP lamSEXP);
SEXP Rfast_rows(SEXP x,SEXP ind);
SEXP Rfast_row_any(SEXP xSEXP);
SEXP Rfast_row_means(SEXP xSEXP);
SEXP Rfast_row_max(SEXP xSEXP);
SEXP Rfast_row_meds(SEXP xSEXP,SEXP na_rmSEXP,SEXP parallelSEXP);
SEXP Rfast_row_min(SEXP xSEXP);
SEXP Rfast_row_len_sort_un_int(SEXP xSEXP);
SEXP Rfast_row_ranks(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP);
SEXP Rfast_row_sums(SEXP xSEXP,SEXP SEXPindices);
SEXP Rfast_rmdp(SEXP ySEXP,SEXP hSEXP,SEXP rndSEXP,SEXP itertimeSEXP);
SEXP Rfast_row_tabulate(SEXP xSEXP,SEXP ncollSEXP);
SEXP Rfast_rep_col(SEXP xSEXP,SEXP nSEXP);
SEXP Rfast_rep_row(SEXP xSEXP,SEXP nSEXP);
SEXP Rfast_row_nth(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP);
SEXP Rfast_row_min_max(SEXP x);
SEXP Rfast_row_shuffle(SEXP xSEXP);
SEXP Rfast_Round(SEXP x,SEXP dgSEXP,SEXP na_rmSEXP);
SEXP Rfast_read_directory(SEXP pathSEXP);
SEXP Rfast_rvmf_h(SEXP xSEXP,SEXP caSEXP,SEXP d1SEXP,SEXP x0SEXP,SEXP mSEXP,SEXP kSEXP,SEXP bSEXP);
SEXP Rfast_row_all(SEXP xSEXP);
SEXP Rfast_row_true(SEXP xSEXP);
SEXP Rfast_row_prods(SEXP xSEXP);
SEXP Rfast_row_false(SEXP xSEXP);
SEXP Rfast_row_order(SEXP xSEXP,SEXP stableSEXP,SEXP descendingSEXP);
SEXP Rfast_row_true_false(SEXP xSEXP);
SEXP Rfast_read_examples(SEXP path_manSEXP);
SEXP Rfast_row_count_values(SEXP xSEXP,SEXP valuesSEXP);
SEXP Rfast_row_mads(SEXP xSEXP,SEXP methodSEXP,SEXP na_rmSEXP,SEXP parallelSEXP);
SEXP Rfast_row_vars(SEXP xSEXP,SEXP parallelSEXP);
SEXP Rfast_row_max_indices(SEXP xSEXP);
SEXP Rfast_row_min_indices(SEXP xSEXP);
SEXP Rfast_sort_unique_double(SEXP xSEXP);
SEXP Rfast_sort_mat(SEXP xSEXP,SEXP descendSEXP,SEXP by_rowSEXP,SEXP stableSEXP,SEXP parallelSEXP);
SEXP Rfast_submatrix(SEXP xSEXP,SEXP rowstartSEXP,SEXP rowendSEXP,SEXP colstartSEXP,SEXP colendSEXP);
SEXP Rfast_sum_XopY(SEXP x,SEXP y,SEXP operSEXP);
SEXP Rfast_sum_XopX(SEXP x,SEXP operSEXP);
SEXP Rfast_sum_lower_tri(SEXP xSEXP,SEXP dgSEXP);
SEXP Rfast_sum_upper_tri(SEXP xSEXP,SEXP dgSEXP);
SEXP Rfast_sort_unique_int(SEXP xSEXP);
SEXP Rfast_symmetric(SEXP xSEXP);
SEXP Rfast_Sort(SEXP x,SEXP descendSEXP,SEXP na);
SEXP Rfast_Sort_na_first(SEXP xSEXP,SEXP descendSEXP);
SEXP Rfast_sort_string(SEXP xSEXP,SEXP descendSEXP);
SEXP Rfast_sort_int(SEXP xSEXP);
SEXP Rfast_stable_sort(SEXP x,SEXP descendSEXP);
SEXP Rfast_spat_med(SEXP xSEXP,SEXP tolSEXP);
SEXP Rfast_squareform_c(SEXP xSEXP);
SEXP Rfast_total_dists(SEXP xSEXP,SEXP methodSEXP,SEXP sqrSEXP,SEXP pSEXP);
SEXP Rfast_total_dista(SEXP xSEXP,SEXP ySEXP,SEXP sqrSEXP);
SEXP Rfast_topological_sort(SEXP dagSEXP);
SEXP Rfast_Trigamma(SEXP x);
SEXP Rfast_table_c(SEXP x,SEXP use_naSEXP);
SEXP Rfast_table_with_names(SEXP x);
SEXP Rfast_table_sign(SEXP xSEXP,SEXP naSEXP,SEXP namesSEXP);
SEXP Rfast_table2_c(SEXP x,SEXP y,SEXP rm_zerosSEXP);
SEXP Rfast_table2_with_names(SEXP x,SEXP y,SEXP rm_zerosSEXP);
SEXP Rfast_transpose(SEXP xSEXP);
SEXP Rfast_upper_tri(SEXP xSEXP,SEXP dgSEXP);
SEXP Rfast_upper_tri_assign(SEXP xSEXP,SEXP vSEXP,SEXP dgSEXP);
SEXP Rfast_upper_tri_b(SEXP nclSEXP, SEXP nrwSEXP,SEXP dgSEXP);
SEXP Rfast_var(SEXP xSEXP,SEXP stdSEXP,SEXP na_rmSEXP);
SEXP Rfast_comb_n(SEXP dataSEXP,SEXP nSEXP,SEXP simplifySEXP);
SEXP Rfast_varcomps_mle(SEXP xSEXP,SEXP inaSEXP,SEXP nSEXP,SEXP tolSEXP);
SEXP Rfast_vecdist(SEXP x);
SEXP Rfast_which_is(SEXP xSEXP,SEXP methodSEXP);
SEXP Rfast_col_row_value(SEXP xSEXP,SEXP vSEXP);

SEXP Rfast_colrint_mle(SEXP XSEXP,SEXP idSEXP,SEXP ranefSEXP,SEXP tolSEXP,SEXP maxitersSEXP,SEXP parallelSEXP);
SEXP Rfast_eigs_sym_c(SEXP XSEXP,SEXP kSEXP,SEXP vectorsSEXP);
SEXP Rfast_geom_regs(SEXP YSEXP,SEXP XSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP typeSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast_dir_knn(SEXP tXnewSEXP,SEXP tXSEXP,SEXP YSEXP,SEXP KSEXP,SEXP typeSEXP,SEXP parallelSEXP);
SEXP Rfast_multinom_regs(SEXP Y0SEXP,SEXP X0SEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP maxitersSEXP,SEXP parallelSEXP);
SEXP Rfast_normlog_regs(SEXP YSEXP,SEXP XSEXP,SEXP BESEXP,SEXP conSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast_normlog_reg(SEXP YSEXP,SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast_rint_reg(SEXP YSEXP,SEXP XSEXP,SEXP BESEXP,SEXP tolSEXP,SEXP ranefSEXP,SEXP loggedSEXP);
SEXP Rfast_rint_regs(SEXP XSEXP,SEXP YSEXP,SEXP idSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast_rint_mle(SEXP XSEXP,SEXP idSEXP,SEXP ranefSEXP,SEXP tolSEX,SEXP maxitersSEXP);
SEXP Rfast_weibull_mle(SEXP XSEXP,SEXP tolSEX,SEXP maxitersSEXP);
SEXP Rfast_weib_reg(SEXP YSEXP,SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast_spml_mle(SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast_spml_regs(SEXP YSEXP,SEXP X0SEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP maxitersSEXP,SEXP parallelSEXP);
SEXP Rfast_spml_reg(SEXP YSEXP, SEXP XSEXP,SEXP tolSEXP,SEXP sebSEXP,SEXP maxitersSEXP);
SEXP Rfast_colweibull_mle(SEXP XSEXP,SEXP tolSEXP,SEXP maxitersSEXP,SEXP parallelSEXP);
SEXP Rfast_quasi_poisson_only(SEXP xSEXP,SEXP ySEXP,SEXP ylogySEXP,SEXP tolSEXP,SEXP maxitersSEXP);

SEXP Rfast_col_all_p(SEXP xSEXP);
SEXP Rfast_col_count_values_p(SEXP xSEXP,SEXP valuesSEXP);
SEXP Rfast_col_max_p(SEXP nSEXP);
SEXP Rfast_col_mean_p(SEXP xSEXP);
SEXP Rfast_col_min_p(SEXP nSEXP);
SEXP Rfast_col_nth_p(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP);
SEXP Rfast_col_order_p(SEXP xSEXP,SEXP stableSEXP,SEXP descendingSEXP);
SEXP Rfast_col_ranks_p(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP);
SEXP Rfast_col_sums_p(SEXP xSEXP);
SEXP Rfast_mat_mult_p(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast_row_all_p(SEXP xSEXP);
SEXP Rfast_row_count_values_p(SEXP xSEXP,SEXP valuesSEXP);
SEXP Rfast_row_nth_p(SEXP xSEXP,SEXP ySEXP,SEXP descendSEXP,SEXP na_rmSEXP,SEXP indexSEXP);
SEXP Rfast_row_order_p(SEXP xSEXP,SEXP stableSEXP,SEXP descendingSEXP);
SEXP Rfast_row_ranks_p(SEXP xSEXP,SEXP methodSEXP,SEXP descendSEXP,SEXP stableSEXP);
SEXP Rfast_row_sums_p(SEXP xSEXP);

static const R_CallMethodDef CallEntries[] = {
  {"Rfast_add_to_namespace", (DL_FUNC) &Rfast_add_to_namespace, 2},
  {"Rfast_apply_condition", (DL_FUNC) &Rfast_apply_condition, 4},
  {"Rfast_as_integer", (DL_FUNC) &Rfast_as_integer, 3},
  {"Rfast_bcdcor", (DL_FUNC) &Rfast_bcdcor, 2},
  {"Rfast_binarysearch", (DL_FUNC) &Rfast_binarysearch, 2},
  {"Rfast_bic_fs_reg", (DL_FUNC) &Rfast_bic_fs_reg, 4},
  {"Rfast_bs_reg", (DL_FUNC) &Rfast_bs_reg, 4},
  {"Rfast_bincomb", (DL_FUNC) &Rfast_bincomb, 1},
  {"Rfast_col_all", (DL_FUNC) &Rfast_col_all, 1},
  {"Rfast_col_count_values", (DL_FUNC) &Rfast_col_count_values, 2},
  {"Rfast_col_cum_maxs", (DL_FUNC) &Rfast_col_cum_maxs, 1},
  {"Rfast_col_cum_sums", (DL_FUNC) &Rfast_col_cum_sums, 1},
  {"Rfast_col_cum_mins", (DL_FUNC) &Rfast_col_cum_mins, 1},
  {"Rfast_col_cum_prods", (DL_FUNC) &Rfast_col_cum_prods, 1},
  {"Rfast_col_meds", (DL_FUNC) &Rfast_col_meds, 3},
  {"Rfast_col_min_indices", (DL_FUNC) &Rfast_col_min_indices, 1},
  {"Rfast_col_min", (DL_FUNC) &Rfast_col_min, 1},
  {"Rfast_col_sums", (DL_FUNC) &Rfast_col_sums, 2},
  {"Rfast_col_min_max", (DL_FUNC) &Rfast_col_min_max, 1},
  {"Rfast_col_max_indices", (DL_FUNC) &Rfast_col_max_indices, 1},
  {"Rfast_col_max", (DL_FUNC) &Rfast_col_max, 1},
  {"Rfast_col_means", (DL_FUNC) &Rfast_col_means, 1},
  {"Rfast_col_nth", (DL_FUNC) &Rfast_col_nth, 5},
  {"Rfast_col_len_sort_un_int", (DL_FUNC) &Rfast_col_len_sort_un_int, 1},
  {"Rfast_col_ranks", (DL_FUNC) &Rfast_col_ranks, 4},
  {"Rfast_col_tabulate", (DL_FUNC) &Rfast_col_tabulate, 2},
  {"Rfast_col_shuffle", (DL_FUNC) &Rfast_col_shuffle, 1},
  {"Rfast_col_true", (DL_FUNC) &Rfast_col_true, 1},
  {"Rfast_col_diffs", (DL_FUNC) &Rfast_col_diffs, 1},
  {"Rfast_col_false", (DL_FUNC) &Rfast_col_false, 1},
  {"Rfast_col_order", (DL_FUNC) &Rfast_col_order, 3},
  {"Rfast_col_pmax", (DL_FUNC) &Rfast_col_pmax, 2},
  {"Rfast_col_pmin", (DL_FUNC) &Rfast_col_pmin, 2},
  {"Rfast_col_prods", (DL_FUNC) &Rfast_col_prods, 2},
  {"Rfast_col_vars", (DL_FUNC) &Rfast_col_vars, 4},
  {"Rfast_columns", (DL_FUNC) &Rfast_columns, 2},
  {"Rfast_chi2Test_univariate", (DL_FUNC) &Rfast_chi2Test_univariate, 2},
  {"Rfast_chi2Test", (DL_FUNC) &Rfast_chi2Test, 5},
  {"Rfast_chi2tests", (DL_FUNC) &Rfast_chi2tests, 4},
  {"Rfast_check_namespace", (DL_FUNC) &Rfast_check_namespace, 2},
  {"Rfast_check_aliases", (DL_FUNC) &Rfast_check_aliases, 2},
  {"Rfast_check_true_false", (DL_FUNC) &Rfast_check_true_false, 1},
  {"Rfast_check_usage", (DL_FUNC) &Rfast_check_usage, 2},
  {"Rfast_col_any", (DL_FUNC) &Rfast_col_any, 1},
  {"Rfast_col_anovas", (DL_FUNC) &Rfast_col_anovas, 2},
  {"Rfast_cholesky", (DL_FUNC) &Rfast_cholesky, 1},
  {"Rfast_cholesky_par", (DL_FUNC) &Rfast_cholesky_par, 1},
  {"Rfast_col_mads", (DL_FUNC) &Rfast_col_mads, 4},
  {"Rfast_col_true_false", (DL_FUNC) &Rfast_col_true_false, 1},
  {"Rfast_count_value", (DL_FUNC) &Rfast_count_value, 2},
  {"Rfast_Choose", (DL_FUNC) &Rfast_Choose, 2},
  {"Rfast_dvar", (DL_FUNC) &Rfast_dvar, 1},
  {"Rfast_dcor", (DL_FUNC) &Rfast_dcor, 2},
  {"Rfast_dcov", (DL_FUNC) &Rfast_dcov, 2},
  {"Rfast_diag_matrix_fill_scalar", (DL_FUNC) &Rfast_diag_matrix_fill_scalar, 2},
  {"Rfast_diag_matrix_fill_vec", (DL_FUNC) &Rfast_diag_matrix_fill_vec, 2},
  {"Rfast_diag_fill_scalar", (DL_FUNC) &Rfast_diag_fill_scalar, 2},
  {"Rfast_diag_fill_vec", (DL_FUNC) &Rfast_diag_fill_vec, 2},
  {"Rfast_design_matrix", (DL_FUNC) &Rfast_design_matrix, 2},
  {"Rfast_dist", (DL_FUNC) &Rfast_dist, 4},
  {"Rfast_dist_vec", (DL_FUNC) &Rfast_dist_vec, 4},
  {"Rfast_Digamma", (DL_FUNC) &Rfast_Digamma, 1},
  {"Rfast_design_matrix_big", (DL_FUNC) &Rfast_design_matrix_big, 1},
  {"Rfast_dista", (DL_FUNC) &Rfast_dista, 4},
  {"Rfast_dista_index", (DL_FUNC) &Rfast_dista_index, 4},
  {"Rfast_dista_values", (DL_FUNC) &Rfast_dista_values, 5},
  {"Rfast_diri_nr_type2", (DL_FUNC) &Rfast_diri_nr_type2, 5},
  {"Rfast_edist", (DL_FUNC) &Rfast_edist, 2},
  {"Rfast_eachrow", (DL_FUNC) &Rfast_eachrow, 4},
  {"Rfast_eachcol_apply", (DL_FUNC) &Rfast_eachcol_apply, 5},
  {"Rfast_floyd_john", (DL_FUNC) &Rfast_floyd_john, 2},
  {"Rfast_frame_to_matrix", (DL_FUNC) &Rfast_frame_to_matrix, 1},
  {"Rfast_fs_reg", (DL_FUNC) &Rfast_fs_reg, 5},
  {"Rfast_g2Test_univariate_perm", (DL_FUNC) &Rfast_g2Test_univariate_perm, 3},
  {"Rfast_g2Test_univariate", (DL_FUNC) &Rfast_g2Test_univariate, 2},
  {"Rfast_g2Test_perm", (DL_FUNC) &Rfast_g2Test_perm, 6},
  {"Rfast_g2Test", (DL_FUNC) &Rfast_g2Test, 5},
  {"Rfast_g2tests_perm", (DL_FUNC) &Rfast_g2tests_perm, 5},
  {"Rfast_g2tests", (DL_FUNC) &Rfast_g2tests, 4},
  {"Rfast_group", (DL_FUNC) &Rfast_group, 6},
  {"Rfast_group_sum", (DL_FUNC) &Rfast_group_sum, 4},
  {"Rfast_group_all", (DL_FUNC) &Rfast_group_all, 3},
  {"Rfast_group_any", (DL_FUNC) &Rfast_group_any, 3},
  {"Rfast_group_mad", (DL_FUNC) &Rfast_group_mad, 3},
  {"Rfast_group_max", (DL_FUNC) &Rfast_group_max, 4},
  {"Rfast_group_mean", (DL_FUNC) &Rfast_group_mean, 3},
  {"Rfast_group_med", (DL_FUNC) &Rfast_group_med, 3},
  {"Rfast_group_min", (DL_FUNC) &Rfast_group_min, 3},
  {"Rfast_group_min_max", (DL_FUNC) &Rfast_group_min_max, 3},
  {"Rfast_group_var", (DL_FUNC) &Rfast_group_var, 3},
  {"Rfast_glm_logistic", (DL_FUNC) &Rfast_glm_logistic, 4},
  {"Rfast_glm_poisson", (DL_FUNC) &Rfast_glm_poisson, 4},
  {"Rfast_hash2list", (DL_FUNC) &Rfast_hash2list, 2},
  {"Rfast_hash_find", (DL_FUNC) &Rfast_hash_find, 2},
  {"Rfast_Hash_list", (DL_FUNC) &Rfast_Hash_list, 2},
  {"Rfast_Hash_key_multi", (DL_FUNC) &Rfast_Hash_key_multi, 3},
  {"Rfast_is_element", (DL_FUNC) &Rfast_is_element, 2},
  {"Rfast_is_element_string", (DL_FUNC) &Rfast_is_element_string, 2},
  {"Rfast_is_integer", (DL_FUNC) &Rfast_is_integer, 1},
  {"Rfast_k_nn", (DL_FUNC) &Rfast_k_nn, 9},
  {"Rfast_k_nn_cv", (DL_FUNC) &Rfast_k_nn_cv, 10},
  {"Rfast_lowerbound", (DL_FUNC) &Rfast_lowerbound, 2},
  {"Rfast_Lgamma", (DL_FUNC) &Rfast_Lgamma, 1},
  {"Rfast_len_sort_unique_int", (DL_FUNC) &Rfast_len_sort_unique_int, 1},
  {"Rfast_Log", (DL_FUNC) &Rfast_Log, 1},
  {"Rfast_logistic_only_b", (DL_FUNC) &Rfast_logistic_only_b, 3},
  {"Rfast_logistic_only", (DL_FUNC) &Rfast_logistic_only, 3},
  {"Rfast_Lbeta", (DL_FUNC) &Rfast_Lbeta, 2},
  {"Rfast_lower_tri", (DL_FUNC) &Rfast_lower_tri, 2},
  {"Rfast_lower_tri_assign", (DL_FUNC) &Rfast_lower_tri_assign, 3},
  {"Rfast_lower_tri_b", (DL_FUNC) &Rfast_lower_tri_b, 3},
  {"Rfast_Lchoose", (DL_FUNC) &Rfast_Lchoose, 2},
  {"Rfast_mahaCpp", (DL_FUNC) &Rfast_mahaCpp, 4},
  {"Rfast_min_freq_d", (DL_FUNC) &Rfast_min_freq_d, 2},
  {"Rfast_min_freq_i", (DL_FUNC) &Rfast_min_freq_i, 2},
  {"Rfast_max_freq_d", (DL_FUNC) &Rfast_max_freq_d, 2},
  {"Rfast_max_freq_i", (DL_FUNC) &Rfast_max_freq_i, 2},
  {"Rfast_Match", (DL_FUNC) &Rfast_Match, 2},
  {"Rfast_mad2", (DL_FUNC) &Rfast_mad2, 3},
  {"Rfast_min_max", (DL_FUNC) &Rfast_min_max, 2},
  {"Rfast_mat_mat", (DL_FUNC) &Rfast_mat_mat, 2},
  {"Rfast_med", (DL_FUNC) &Rfast_med, 2},
  {"Rfast_min_max_perc", (DL_FUNC) &Rfast_min_max_perc, 1},
  {"Rfast_negative", (DL_FUNC) &Rfast_negative, 2},
  {"Rfast_nth", (DL_FUNC) &Rfast_nth, 5},
  {"Rfast_nth_int", (DL_FUNC) &Rfast_nth_int, 2},
  {"Rfast_Norm", (DL_FUNC) &Rfast_Norm, 2},
  {"Rfast_Order", (DL_FUNC) &Rfast_Order, 3},
  {"Rfast_Outer", (DL_FUNC) &Rfast_Outer, 3},
  {"Rfast_odds_helper", (DL_FUNC) &Rfast_odds_helper, 1},
  {"Rfast_positive", (DL_FUNC) &Rfast_positive, 2},
  {"Rfast_positive_negative", (DL_FUNC) &Rfast_positive_negative, 2},
  {"Rfast_poisson_only", (DL_FUNC) &Rfast_poisson_only, 4},
  {"Rfast_partial_sort", (DL_FUNC) &Rfast_partial_sort, 3},
  {"Rfast_partial_sort_index", (DL_FUNC) &Rfast_partial_sort_index, 3},
  {"Rfast_prop_reg", (DL_FUNC) &Rfast_prop_reg, 4},
  {"Rfast_prop_regs", (DL_FUNC) &Rfast_prop_regs, 5},
  {"Rfast_pmin", (DL_FUNC) &Rfast_pmin, 3},
  {"Rfast_pmax", (DL_FUNC) &Rfast_pmax, 3},
  {"Rfast_pmin_pmax", (DL_FUNC) &Rfast_pmin_pmax, 3},
  {"Rfast_pc_skel", (DL_FUNC) &Rfast_pc_skel, 7},
  {"Rfast_poisson_only_b", (DL_FUNC) &Rfast_poisson_only_b, 4},
  {"Rfast_permutation", (DL_FUNC) &Rfast_permutation, 2},
  {"Rfast_permutation_next", (DL_FUNC) &Rfast_permutation_next, 2},
  {"Rfast_permutation_prev", (DL_FUNC) &Rfast_permutation_prev, 2},
  {"Rfast_perm_cor", (DL_FUNC) &Rfast_perm_cor, 3},
  {"Rfast_qpois_reg", (DL_FUNC) &Rfast_qpois_reg, 5},
  {"Rfast_qpois_regs", (DL_FUNC) &Rfast_qpois_regs, 5},
  {"Rfast_rank", (DL_FUNC) &Rfast_rank, 3},
  {"Rfast_remove_from_namespace", (DL_FUNC) &Rfast_remove_from_namespace, 2},
  {"Rfast_rbing", (DL_FUNC) &Rfast_rbing, 2},
  {"Rfast_rows", (DL_FUNC) &Rfast_rows, 2},
  {"Rfast_row_count_values", (DL_FUNC) &Rfast_row_count_values, 2},
  {"Rfast_row_any", (DL_FUNC) &Rfast_row_any, 1},
  {"Rfast_row_means", (DL_FUNC) &Rfast_row_means, 1},
  {"Rfast_row_max", (DL_FUNC) &Rfast_row_max, 1},
  {"Rfast_row_meds", (DL_FUNC) &Rfast_row_meds, 3},
  {"Rfast_row_min", (DL_FUNC) &Rfast_row_min, 1},
  {"Rfast_row_len_sort_un_int", (DL_FUNC) &Rfast_row_len_sort_un_int, 1},
  {"Rfast_row_ranks", (DL_FUNC) &Rfast_row_ranks, 4},
  {"Rfast_row_sums", (DL_FUNC) &Rfast_row_sums, 2},
  {"Rfast_rmdp", (DL_FUNC) &Rfast_rmdp, 4},
  {"Rfast_row_tabulate", (DL_FUNC) &Rfast_row_tabulate, 2},
  {"Rfast_rep_col", (DL_FUNC) &Rfast_rep_col, 2},
  {"Rfast_rep_row", (DL_FUNC) &Rfast_rep_row, 2},
  {"Rfast_row_nth", (DL_FUNC) &Rfast_row_nth, 5},
  {"Rfast_row_min_max", (DL_FUNC) &Rfast_row_min_max, 1},
  {"Rfast_row_shuffle", (DL_FUNC) &Rfast_row_shuffle, 1},
  {"Rfast_Round", (DL_FUNC) &Rfast_Round, 3},
  {"Rfast_read_directory", (DL_FUNC) &Rfast_read_directory, 1},
  {"Rfast_rvmf_h", (DL_FUNC) &Rfast_rvmf_h, 7},
  {"Rfast_row_all", (DL_FUNC) &Rfast_row_all, 1},
  {"Rfast_row_true", (DL_FUNC) &Rfast_row_true, 1},
  {"Rfast_row_prods", (DL_FUNC) &Rfast_row_prods, 1},
  {"Rfast_row_false", (DL_FUNC) &Rfast_row_false, 1},
  {"Rfast_row_order", (DL_FUNC) &Rfast_row_order, 3},
  {"Rfast_row_true_false", (DL_FUNC) &Rfast_row_true_false, 1},
  {"Rfast_read_examples", (DL_FUNC) &Rfast_read_examples, 1},
  {"Rfast_row_mads", (DL_FUNC) &Rfast_row_mads, 4},
  {"Rfast_row_vars", (DL_FUNC) &Rfast_row_vars, 4},
  {"Rfast_row_max_indices", (DL_FUNC) &Rfast_row_max_indices, 1},
  {"Rfast_row_min_indices", (DL_FUNC) &Rfast_row_min_indices, 1},
  {"Rfast_sort_unique_double", (DL_FUNC) &Rfast_sort_unique_double, 1},
  {"Rfast_sort_mat", (DL_FUNC) &Rfast_sort_mat, 5},
  {"Rfast_submatrix", (DL_FUNC) &Rfast_submatrix, 5},
  {"Rfast_sum_XopY", (DL_FUNC) &Rfast_sum_XopY, 3},
  {"Rfast_sum_XopX", (DL_FUNC) &Rfast_sum_XopX, 2},
  {"Rfast_sum_lower_tri", (DL_FUNC) &Rfast_sum_lower_tri, 2},
  {"Rfast_sum_upper_tri", (DL_FUNC) &Rfast_sum_upper_tri, 2},
  {"Rfast_sort_unique_int", (DL_FUNC) &Rfast_sort_unique_int, 1},
  {"Rfast_symmetric", (DL_FUNC) &Rfast_symmetric, 1},
  {"Rfast_Sort", (DL_FUNC) &Rfast_Sort, 3},
  {"Rfast_Sort_na_first", (DL_FUNC) &Rfast_Sort_na_first, 2},
  {"Rfast_sort_string", (DL_FUNC) &Rfast_sort_string, 2},
  {"Rfast_stable_sort", (DL_FUNC) &Rfast_stable_sort, 2},
  {"Rfast_sort_int", (DL_FUNC) &Rfast_sort_int, 1},
  {"Rfast_spat_med", (DL_FUNC) &Rfast_spat_med, 2},
  {"Rfast_squareform_c", (DL_FUNC) &Rfast_squareform_c, 1},
  {"Rfast_total_dists", (DL_FUNC) &Rfast_total_dists, 4},
  {"Rfast_total_dista", (DL_FUNC) &Rfast_total_dista, 3},
  {"Rfast_topological_sort", (DL_FUNC) &Rfast_topological_sort, 1},
  {"Rfast_Trigamma", (DL_FUNC) &Rfast_Trigamma, 1},
  {"Rfast_table_c", (DL_FUNC) &Rfast_table_c, 2},
  {"Rfast_table_sign", (DL_FUNC) &Rfast_table_sign, 3},
  {"Rfast_table_with_names", (DL_FUNC) &Rfast_table_with_names, 1},
  {"Rfast_table2_c", (DL_FUNC) &Rfast_table2_c, 3},
  {"Rfast_table2_with_names", (DL_FUNC) &Rfast_table2_with_names, 3},
  {"Rfast_transpose", (DL_FUNC) &Rfast_transpose, 1},
  {"Rfast_upper_tri", (DL_FUNC) &Rfast_upper_tri, 2},
  {"Rfast_upper_tri_assign", (DL_FUNC) &Rfast_upper_tri_assign, 3},
  {"Rfast_upper_tri_b", (DL_FUNC) &Rfast_upper_tri_b, 3},
  {"Rfast_var", (DL_FUNC) &Rfast_var, 3},
  {"Rfast_varcomps_mle", (DL_FUNC) &Rfast_varcomps_mle, 4},
  {"Rfast_vecdist", (DL_FUNC) &Rfast_vecdist, 1},
  {"Rfast_comb_n", (DL_FUNC) &Rfast_comb_n, 3},
  {"Rfast_which_is", (DL_FUNC) &Rfast_which_is, 2},
  {"Rfast_col_row_value", (DL_FUNC) &Rfast_col_row_value, 2},


  {"Rfast_col_all_p", (DL_FUNC) &Rfast_col_all_p, 1},
  {"Rfast_col_count_values_p", (DL_FUNC) &Rfast_col_count_values_p, 2},
  {"Rfast_col_max_p", (DL_FUNC) &Rfast_col_max_p, 1},
  {"Rfast_col_mean_p", (DL_FUNC) &Rfast_col_mean_p, 1},
  {"Rfast_col_min_p", (DL_FUNC) &Rfast_col_min_p, 1},
  {"Rfast_col_nth_p", (DL_FUNC) &Rfast_col_nth_p, 5},
  {"Rfast_col_ranks_p", (DL_FUNC) &Rfast_col_ranks_p, 4},
  {"Rfast_col_sums_p", (DL_FUNC) &Rfast_col_sums_p, 1},
  {"Rfast_col_order_p", (DL_FUNC) &Rfast_col_order_p, 3},
  {"Rfast_mat_mult_p", (DL_FUNC) &Rfast_mat_mult_p, 2},
  {"Rfast_row_all_p", (DL_FUNC) &Rfast_row_all_p, 1},
  {"Rfast_row_count_values_p", (DL_FUNC) &Rfast_row_count_values_p, 2},
  {"Rfast_row_nth_p", (DL_FUNC) &Rfast_row_nth_p, 5},
  {"Rfast_row_order_p", (DL_FUNC) &Rfast_row_order_p, 3},
  {"Rfast_row_ranks_p", (DL_FUNC) &Rfast_row_ranks_p, 4},
  {"Rfast_row_sums_p", (DL_FUNC) &Rfast_row_sums_p, 1},

  {"Rfast_colrint_mle", (DL_FUNC) &Rfast_colrint_mle, 6},
  {"Rfast_eigs_sym_c", (DL_FUNC) &Rfast_eigs_sym_c, 3},
  {"Rfast_geom_regs", (DL_FUNC) &Rfast_geom_regs, 7},
  {"Rfast_normlog_regs", (DL_FUNC) &Rfast_normlog_regs, 8},
  {"Rfast_dir_knn", (DL_FUNC) &Rfast_dir_knn, 6},
  {"Rfast_multinom_regs", (DL_FUNC) &Rfast_multinom_regs, 6},
  {"Rfast_normlog_reg", (DL_FUNC) &Rfast_normlog_reg, 4},
  {"Rfast_rint_reg", (DL_FUNC) &Rfast_rint_reg, 6},
  {"Rfast_rint_regs", (DL_FUNC) &Rfast_rint_regs, 7},
  {"Rfast_rint_mle", (DL_FUNC) &Rfast_rint_mle, 5},
  {"Rfast_weibull_mle", (DL_FUNC) &Rfast_weibull_mle, 3},
  {"Rfast_weib_reg", (DL_FUNC) &Rfast_weib_reg, 4},
  {"Rfast_spml_mle", (DL_FUNC) &Rfast_spml_mle, 3},
  {"Rfast_spml_regs", (DL_FUNC) &Rfast_spml_regs, 6},
  {"Rfast_spml_reg", (DL_FUNC) &Rfast_spml_reg, 5},
  {"Rfast_colweibull_mle", (DL_FUNC) &Rfast_colweibull_mle, 4},
  {"Rfast_quasi_poisson_only", (DL_FUNC) &Rfast_quasi_poisson_only, 5},
  {NULL, NULL, 0}
};


void R_init_Rfast(DllInfo *info)
{
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

