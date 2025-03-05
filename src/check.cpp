// Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <algorithm>
#include "system_files.h"

using namespace Rcpp;
using namespace std;

using std::binary_search;
using std::getline;
using std::ifstream;
using std::remove;
using std::sort;
using std::string;
using std::vector;

vector<string> check_namespace(const string dir_to_export, const string dir_to_file, const bool full_paths = false)
{
	List all_functions = read_functions_and_signatures(dir_to_file, full_paths)["export"];
	vector<string> all_r_functions = all_functions["functions"];
	int which_string_has_export = 0, len_which_not_exp = 1;
	vector<string> which_undefined_function, all_exported_files;
	if (all_r_functions.empty())
	{
		stop("Warning: empty folder.\n");
	}
	vector<string> data_export = readNamespaceFile(dir_to_export, which_string_has_export);
	if (which_string_has_export == -1)
	{
		stop("Error. can't find \"export\" function in NAMESPACE file.\n");
	}
	string exported_files = data_export[which_string_has_export];
	exported_files.erase(exported_files.end() - 1);
	exported_files.erase(exported_files.begin(), exported_files.begin() + 7);
	all_exported_files = split_words(exported_files, ",");
	sort(all_r_functions.begin(), all_r_functions.end());
	for (unsigned int i = 0; i < all_exported_files.size(); ++i)
	{
		if (binary_search(all_r_functions.begin(), all_r_functions.end(), all_exported_files[i]) == false)
		{
			which_undefined_function.resize(len_which_not_exp);
			which_undefined_function[len_which_not_exp - 1] = all_exported_files[i];
			++len_which_not_exp;
		}
	}
	return which_undefined_function;
}

RcppExport SEXP Rfast_check_namespace(SEXP dir_to_exportSEXP, SEXP dir_to_fileSEXP, SEXP full_pathsSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const string>::type dir_to_export(dir_to_exportSEXP);
	traits::input_parameter<const string>::type dir_to_file(dir_to_fileSEXP);
	traits::input_parameter<const bool>::type full_paths(full_pathsSEXP);
	__result = check_namespace(dir_to_export, dir_to_file, full_paths);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////////////////

List check_true_false(string path_to_man, const bool full_paths = false)
{
	List ex = read_examples(path_to_man, full_paths), L;
	CharacterVector names = ex["files"];
	vector<string> exams = as<vector<string>>(ex["examples"]);
	string s;
	CharacterVector trues, falses;
	for (unsigned int i = 0; i < exams.size(); ++i)
	{
		s = exams[i];
		[[maybe_unused]] auto temp = remove(s.begin(), s.end(), ' ');
		if (find_string(s, "=T)") || find_string(s, "=T,"))
		{
			trues.push_back(names[i]);
		}
		else if (find_string(s, "=F)") || find_string(s, "=F,"))
		{
			falses.push_back(names[i]);
		}
	}
	trues = sort_unique(trues);
	falses = sort_unique(falses);
	if (trues.size())
		L["TRUE"] = trues;
	if (falses.size())
		L["FALSE"] = falses;
	if (ex.containsElementNamed("dont read"))
		L["dont read"] = ex["dont read"];
	return L;
}

RcppExport SEXP Rfast_check_true_false(SEXP path_manSEXP, SEXP full_pathsSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type path_man(path_manSEXP);
	traits::input_parameter<const bool>::type full_paths(full_pathsSEXP);
	__result = check_true_false(path_man, full_paths);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////////////
using std::remove;

//[[Rcpp::export]]
List check_aliases(const string path_man, const string path_rf, const bool full_paths = false)
{
	ifstream file;
	List data = read_functions_and_signatures(path_rf, full_paths);
	List all_functions = data["export"];
	vector<string> aliases, all_r_functions = all_functions["functions"],
							all_s3method = all_functions["s3"], tmp, dontread_rd;
	Files all_rd_files = readDirectory(path_man);
	all_r_functions.reserve(all_r_functions.size() + all_s3method.size());
	all_r_functions.insert(all_r_functions.end(), all_s3method.begin(), all_s3method.end());
	for (auto &rd_file : all_rd_files)
	{
		file.open(rd_file.filename(true));
		if (!file.is_open())
		{
			Rcout << "Can't open file " << rd_file.filename(full_paths) << ".\n";
			continue;
		}
		if (check_read_file(file, '%'))
		{
			tmp = read_aliases(file);
			aliases.reserve(aliases.size() + tmp.size());
			aliases.insert(aliases.end(), tmp.begin(), tmp.end());
		}
		else
		{
			DEBUG("Find attribute dont read file with name: " + rd_file.filename(full_paths));
			dontread_rd.push_back(rd_file.filename(full_paths));
		}
		file.close();
	}

	sort(aliases.begin(), aliases.end());
	sort(all_r_functions.begin(), all_r_functions.end());
	List ls;
	auto r_to_al = find_which(all_r_functions, aliases);
	auto al_to_r = find_which(aliases, all_r_functions);
	auto dupAl = find_duplis(aliases);

	if (!r_to_al.empty())
		ls["Missing Man aliases"] = r_to_al;
	if (!al_to_r.empty())
		ls["Missing R functions"] = al_to_r;
	if (!dupAl.empty())
		ls["Duplicate alias"] = dupAl;
	if (!dontread_rd.empty())
		ls["dont read"] = List::create(_["Rd"] = dontread_rd);
	return ls;
}

RcppExport SEXP Rfast_check_aliases(SEXP dir_to_manSEXP, SEXP dir_to_fileSEXP, SEXP full_pathsSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const string>::type dir_to_man(dir_to_manSEXP);
	traits::input_parameter<const string>::type dir_to_file(dir_to_fileSEXP);
	traits::input_parameter<const bool>::type full_paths(full_pathsSEXP);
	__result = check_aliases(dir_to_man, dir_to_file, full_paths);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

using std::remove;

//[[Rcpp::export]]
List check_usage(string path_man, string path_rf, const bool full_paths = false)
{
	DEBUG("START");
	ifstream file_rd, file_r;
	vector<string> dontread_rd, dontread_r;
	Files all_rd_files = readDirectory(path_man);
	std::vector<string> exported_functions_names, not_exported_functions_names;
	vector<string> missing_functions, aliases, functions_usage, missmatch_functions, name_of_functions_in_usage, aliases_with_lines_more_than_90;
	string r_file, function_signature;
	List aliases_per_rd_with_lines_more_than_90;

	List functions = read_functions_and_signatures(path_rf, full_paths), functions_signatures = functions["signatures"];

	for (auto &rd_file : all_rd_files)
	{
		file_rd.open(rd_file.filename(true));
		if (!file_rd.is_open())
		{
			// Rcout << "Can't open file " << rd_file.filename(full_paths) << ".\n";
			continue;
		}
		if (!check_read_file(file_rd, '%'))
		{
			dontread_rd.push_back(rd_file.filename(full_paths));
			DEBUG("Find attribute dont read file with name: " + rd_file.filename(full_paths));
		}
		else
		{
			DEBUG("file: ",rd_file.filename(true));
			aliases = read_aliases(file_rd);
			functions_usage = read_usage(file_rd, aliases_with_lines_more_than_90);

			if (!aliases_with_lines_more_than_90.empty())
			{
				aliases_per_rd_with_lines_more_than_90[rd_file.filename()] = aliases_with_lines_more_than_90;
				aliases_with_lines_more_than_90.clear();
			}

			string curr_func, func_from_r_file;
			for (auto &al : aliases)
			{
				DEBUG("\t",al);
				for (auto &tmp : functions_usage)
				{ // sigourevo oti gia to alias iparxei h antistoixi sinartisi sto usage
					if (tmp.compare(0, al.size(), al) == 0 and tmp[al.size()] == '(')
					{ // otan to onoma einai idio akrivos kai amesos meta ksekinaei "("
						curr_func = tmp;
						break;
					}
				}
				if (curr_func.empty())
				{
					missing_functions.push_back("<" + al + "> not in <" + rd_file.filename(full_paths) + ">"); // aliase not in usage
				}
				else
				{
					if (functions_signatures.containsElementNamed(al.c_str()))
					{ //  an iparxei to alias
						List functions_details = functions_signatures[al];
						// r_file = as<string>(functions_details["filename"]); // to onoma toy arxeiou pou iparxei to aliase
						function_signature = as<string>(functions_details["signature"]); // ipografi tis sinartiseis

						//DEBUG("current: "+curr_func+" , fromRfile: "+function_signature);
						if (curr_func != function_signature)
						{
							DEBUG(curr_func + " : " + function_signature + " [" + al + "]");
							missmatch_functions.push_back("signature of <" + al + "> missmatch with usage in <" + rd_file.filename(full_paths) + ">");
						}
					}
					else
					{ // an to alias den iparxei stis sinartiseis simainei oti leipei h R sinartisi
						missing_functions.push_back("<" + al + "> not in any R file");
					}
					curr_func.clear();
				}
			}
		}
		file_rd.close();
	}
	List L, r_rd;
	if (!missing_functions.empty())
		L["missing functions"] = missing_functions;
	if (!missmatch_functions.empty())
		L["missmatch_functions"] = missmatch_functions;
	if (functions.containsElementNamed("dont read"))
		r_rd["R"] = functions["dont read"];
	if (!dontread_rd.empty())
		r_rd["Rd"] = dontread_rd;
	if (r_rd.size() != 0)
		L["dont read"] = r_rd;
	if (aliases_per_rd_with_lines_more_than_90.size() != 0)
		L["usage lines wider than 90 characters"] = aliases_per_rd_with_lines_more_than_90;
	if (functions.containsElementNamed("hidden functions"))
		L["hidden functions"] = functions["hidden functions"];
	DEBUG("END");
	return L;
}

RcppExport SEXP Rfast_check_usage(SEXP path_manSEXP, SEXP path_rfSEXP, SEXP full_pathsSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type path_man(path_manSEXP);
	traits::input_parameter<string>::type path_rf(path_rfSEXP);
	traits::input_parameter<const bool>::type full_paths(full_pathsSEXP);
	__result = check_usage(path_man, path_rf, full_paths);
	return __result;
	END_RCPP
}
