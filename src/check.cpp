// Author: Manos Papadakis

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

vector<string> check_namespace(const string dir_to_export, const string dir_to_file)
{
	int which_string_has_export = 0, len_which_not_exp = 1;
	vector<string> allfiles = readDirectory(dir_to_file, 2), which_undefined_function, all_exported_files;
	if (allfiles.empty())
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
	sort(allfiles.begin(), allfiles.end());
	for (unsigned int i = 0; i < all_exported_files.size(); ++i)
	{
		if (binary_search(allfiles.begin(), allfiles.end(), all_exported_files[i]) == false)
		{
			which_undefined_function.resize(len_which_not_exp);
			which_undefined_function[len_which_not_exp - 1] = all_exported_files[i];
			++len_which_not_exp;
		}
	}
	return which_undefined_function;
}

RcppExport SEXP Rfast_check_namespace(SEXP dir_to_exportSEXP, SEXP dir_to_fileSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const string>::type dir_to_export(dir_to_exportSEXP);
	traits::input_parameter<const string>::type dir_to_file(dir_to_fileSEXP);
	__result = check_namespace(dir_to_export, dir_to_file);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////////////////

List check_true_false(string path_to_man)
{
	List ex = read_examples(path_to_man), L;
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

RcppExport SEXP Rfast_check_true_false(SEXP path_manSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type path_man(path_manSEXP);
	__result = check_true_false(path_man);
	return __result;
	END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////////////////
using std::remove;
//[[Rcpp::export]]
List read_export(const string path_rf)
{
	return read_functions_and_signatures(path_rf);
}

//[[Rcpp::export]]
List check_aliases(const string path_man, const string path_rf)
{
	ifstream file;
	List data = read_functions_and_signatures(path_rf);
	List all_functions = data["export"];
	vector<string> aliases, all_r_functions = all_functions["functions"],
							all_s3method = all_functions["s3"], all_rd_files = readDirectory(path_man, 3), tmp, dontread_rd;
	all_r_functions.reserve(all_r_functions.size() + all_s3method.size());
	all_r_functions.insert(all_r_functions.end(), all_s3method.begin(), all_s3method.end());
	for (auto &rd_file : all_rd_files)
	{
		file.open(path_man + rd_file + ".Rd");
		if (!file.is_open())
		{
			Rcout << "Can't open file " << rd_file << ".\n";
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
			DEBUG("Find attribute dont read file with name: " + rd_file);
			dontread_rd.push_back(rd_file);
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

RcppExport SEXP Rfast_check_aliases(SEXP dir_to_manSEXP, SEXP dir_to_fileSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const string>::type dir_to_man(dir_to_manSEXP);
	traits::input_parameter<const string>::type dir_to_file(dir_to_fileSEXP);
	__result = check_aliases(dir_to_man, dir_to_file);
	return __result;
	END_RCPP
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct File : public ifstream
{
	string name;
	void fopen(string path, string name)
	{
		this->name = name;
		this->open(path + name);
	}
	void fclose()
	{
		this->close();
		name = "";
	}
};
using std::remove;

//[[Rcpp::export]]
List check_usage(string path_man, string path_rf)
{
	DEBUG("START");
	File file_rd, file_r;
	vector<string> all_rd_files = read_directory(path_man), dontread_rd, dontread_r;
	std::vector<string> exported_functions_names, not_exported_functions_names;
	vector<string> missing_functions, aliases, functions_usage, missmatch_functions, name_of_functions_in_usage;
	string r_file, function_signature;

	List functions = read_functions_and_signatures(path_rf), functions_signatures = functions["signatures"];

	for (unsigned int i = 0; i < all_rd_files.size(); ++i)
	{
		file_rd.fopen(path_man, all_rd_files[i]);
		if (!file_rd.is_open())
		{
			Rcout << "Can't open file " << all_rd_files[i] << ".\n";
			continue;
		}
		if (!check_read_file(file_rd, '%'))
		{
			dontread_rd.push_back(file_rd.name);
			DEBUG("Find attribute dont read file with name: " + file_rd.name);
		}
		else
		{

			aliases = read_aliases(file_rd);
			functions_usage = read_usage(file_rd);

			string curr_func, func_from_r_file;
			for (auto &al : aliases)
			{
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
					missing_functions.push_back("<" + al + "> not in <" + file_rd.name + ">"); // aliase not in usage
				}
				else
				{
					if (functions_signatures.containsElementNamed(al.c_str()))
					{ //  an iparxei to alias
						List functions_details = functions_signatures[al];
						// r_file = as<string>(functions_details["filename"]); // to onoma toy arxeiou pou iparxei to aliase
						function_signature = as<string>(functions_details["signature"]); // ipografi tis sinartiseis

						// DEBUG("current: "+curr_func+" , fromRfile: "+function_signature);
						if (curr_func != function_signature)
						{
							DEBUG(curr_func + " : " + function_signature + " [" + al + "]");
							missmatch_functions.push_back("signature of <" + al + "> missmatch with usage in <" + file_rd.name + ">");
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
		file_rd.fclose();
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
	if (functions.containsElementNamed("hidden functions"))
		L["hidden functions"] = functions["hidden functions"];
	DEBUG("END");
	return L;
}

RcppExport SEXP Rfast_check_usage(SEXP path_manSEXP, SEXP path_rfSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<string>::type path_man(path_manSEXP);
	traits::input_parameter<string>::type path_rf(path_rfSEXP);
	__result = check_usage(path_man, path_rf);
	return __result;
	END_RCPP
}
