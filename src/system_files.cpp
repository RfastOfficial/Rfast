// Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "system_files.h"
#include <algorithm>

using namespace std;
using namespace Rcpp;
using std::binary_search;
using std::count;
using std::endl;
using std::find_if;
using std::getline;
using std::ifstream;
using std::lower_bound;
using std::ofstream;
using std::remove;
using std::sort;
using std::string;
using std::strtok;
using std::vector;

void print_error() {}

void reset_file(ifstream &file)
{
    file.clear();
    file.seekg(0, ios::beg);
}

vector<string> split_words(string x, const char *sep = ",")
{
    x.erase(remove(x.begin(), x.end(), ' '), x.end());
    int n = count(x.begin(), x.end(), sep[0]) + 1;
    vector<string> y(n);
    x += sep;
    int i = 0;
    const char *split = sep;
    char *token = strtok(&x[0], split);
    while (token != NULL)
    {
        y[i++] = token;
        token = strtok(NULL, split);
    }
    return y;
}

array<string, 2> split_words_in_half(string x, const char sep)
{
    x.erase(remove(x.begin(), x.end(), ' '), x.end());
    int index_to_sep = find(x.begin(), x.end(), sep) - x.begin();
    return {x.substr(0, index_to_sep), x.substr(index_to_sep + 1, x.size() - 1)};
}

void writeFile(vector<string> f, string path)
{
    ofstream oput(path.c_str());
    if (!oput.is_open())
    {
        Rcpp::stop("can't open file\n");
    }
    for (unsigned int i = 0; i < f.size(); ++i)
    {
        oput << f[i] << endl;
    }
}

vector<string> readNamespaceFile(string path, int &which_string_has_export)
{
    ifstream input(path.c_str());
    string s;
    vector<string> f;
    which_string_has_export = -1;
    bool found_export = false;
    while (getline(input, s))
    {
        if (is_namespace_export(s) and !found_export)
        { // oso briskei to export kai den to exei ksanavrei
            which_string_has_export = f.size();
            found_export = true;
        }
        f.push_back(s);
    }
    return f;
}

bool is_namespace_export(string x)
{
    return x.size() > sizeof("export") and x[0] == 'e' and x[1] == 'x' and x[2] == 'p' and x[3] == 'o' and x[4] == 'r' and x[5] == 't';
}

Files readDirectory(const fs::path path, const string extension)
{
    bool checkExtension = !extension.empty();
    Files files;
    if (fs::exists(path))
    {
        for (const auto& entry : fs::recursive_directory_iterator(path)) { // automatically skips special directories . and ..
            if (fs::is_regular_file(entry) && (!checkExtension || (checkExtension && entry.path().extension() == extension))) {
                files.push_back(entry.path());
            }
        }
    }
    else
    {
        Rcpp::stop("Error: Could not open directory with path \"" + path.generic_string() + "\"");
    }
    return files;
}

vector<string> find_which(vector<string> big, vector<string> small)
{
    vector<string> f;
    for (unsigned int i = 0; i < big.size(); ++i)
    {
        if (binary_search(small.begin(), small.end(), big[i]) == false)
        {
            f.push_back(big[i]);
        }
    }
    return f;
}

vector<string> find_duplis(vector<string> x)
{
    x.push_back("@");
    vector<string>::iterator a = x.begin(), b = a + 1;
    vector<string> f;
    int s = 0;
    for (; b != x.end(); ++b)
    {
        if (*a != *b)
        {
            if (s)
            {
                f.push_back(*a);
            }
            a = b;
            s = 0;
        }
        else
        {
            ++s;
        }
    }
    return f;
}

bool binary_help(vector<string>::iterator first, vector<string>::iterator last, string &val, vector<string>::iterator &res)
{
    vector<string>::iterator t = lower_bound(first, last, val);
    int tt = t - first + 1;
    bool found = false;
    if (tt != last - first and val >= *first)
    {
        res = t;
        found = true;
    }
    return found;
}

void dont_read_man(vector<string> &all_rd_files, vector<string> &no_read)
{
    if (no_read[0] != "")
    {
        sort(all_rd_files.begin(), all_rd_files.end());
        vector<string>::iterator fv;
        for (unsigned int i = 0; i < no_read.size(); ++i)
        {
            if (binary_help(all_rd_files.begin(), all_rd_files.end(), no_read[i], fv))
            {
                all_rd_files.erase(fv);
            }
        }
    }
}

bool is_alias(string &s)
{
    return (s.size() > 5 and s[0] == '\\' and s[1] == 'a' and s[2] == 'l' and s[3] == 'i' and s[4] == 'a' and s[5] == 's');
}

bool is_title(string &s)
{
    return (s.size() > 5 and s[0] == '\\' and s[1] == 't' and s[2] == 'i' and s[3] == 't' and s[4] == 'l' and s[5] == 'e');
}

void remove_alias_and_spaces(string &s)
{
    DEBUG("Start remove_alias_and_spaces");
    s.erase(s.end() - 1);
    s.erase(s.begin(), s.begin() + 7);
    remove_spaces_from_begin_end(s);
    if (is_R_operator(s.substr(0, 2)) or find_string(s, "<-"))
    {
        s = "\"" + s + "\"";
    }
    DEBUG("End remove_alias_and_spaces");
}

vector<string> read_aliases(ifstream &file)
{
    DEBUG("Start read_aliases");
    reset_file(file);
    vector<string> als;
    string s;
    do
    {
        getline(file, s);
        if (is_alias(s))
        {
            remove_alias_and_spaces(s);
            DEBUG(s);
            als.push_back(s);
        }
    } while (!is_title(s));
    DEBUG("End read_aliases");
    return als;
}

bool is_dont_runtest(string &s)
{
    return ((s.size() >= (sizeof("\\dontrun") - 1) or s.size() >= (sizeof("\\donttest")) - 1) and
            s[0] == '\\' and s[1] == 'd' and s[2] == 'o' and s[3] == 'n' and s[4] == 't' and
            ((s[5] == 'r' and s[6] == 'u' and s[7] == 'n') or (s[5] == 't' and s[6] == 'e' and s[7] == 's' and s[8] == 't')));
}

void pass_dont_run(ifstream &file)
{
    string tmp;
    while (getline(file, tmp))
    {
        if (tmp == "}")
        {
            break; // pass all the lines from dont run
        }
    }
}

bool is_example(const char *s, size_t len)
{
    return (len >= (sizeof("\\examples") - 1) and s[0] == '\\' and s[1] == 'e' and s[2] == 'x' and s[3] == 'a' and s[4] == 'm' and s[5] == 'p' and s[6] == 'l' and s[7] == 'e' and s[8] == 's');
}

int get_example(ifstream &file, string &res)
{
    string s;
    int is_e = 0; // not example
    if (getline(file, s))
    {
        is_e = is_example(s.c_str(), s.size());
        res = is_e ? s : "";
    }
    else
    {
        is_e = -1; // failed to read. Maybe EOF found
    }
    return is_e;
}

string read_example(ifstream &file, int &long_lines)
{
    string als;
    string s;
    int found_example;
    unsigned int count_curly_bracket = 1; /*at least there will be an empty example section*/
    do
    {
        found_example = get_example(file, s);
    } while (found_example == 0); // while not found example
    if (found_example > 0)
    {
        do
        {
            if (!getline(file, s))
            {
                break;
            }
            if (is_dont_runtest(s))
            {
                pass_dont_run(file);
                s.clear();
            }
            if (s.size() > 99)
            { // 100 max lines
                ++long_lines;
            }
            for (auto &symbol : s)
            {
                if (symbol == '{')
                    ++count_curly_bracket;
                else if (symbol == '}')
                    --count_curly_bracket;
            }
            s += "\n";
            als += s;
        } while (count_curly_bracket > 0); /* check for {} and extract the example correct*/
        als[als.size() - 2] = '\n';        // replace } with new line
        als.erase(als.end() - 1);          // remove extra new line
    }
    return als;
}

bool is_usage(string &s)
{
    return (s.size() > 5 and s[0] == '\\' and s[1] == 'u' and s[2] == 's' and s[3] == 'a' and s[4] == 'g' and s[5] == 'e');
}

bool get_usage(ifstream &file, string &res)
{
    // DEBUG("Start get_usage");
    string s;
    getline(file, s);
    bool is_e = is_usage(s);
    res = is_e ? s : "";
    // DEBUG("End get_usage");
    return is_e;
}

bool is_method(string &s)
{
    return s.size() >= sizeof("\\method") and s[0] == '\\' and s[1] == 'm' and
           s[2] == 'e' and s[3] == 't' and s[4] == 'h' and s[5] == 'o' and s[6] == 'd';
}

string convert_method_to_function(string s)
{
    int position_of_assign;

    for (position_of_assign = s.size() - 1; position_of_assign >= 0; --position_of_assign)
    {
        if (s[position_of_assign] == '-' and s[position_of_assign - 1] == '<')
        {
            --position_of_assign;
            break;
        }
        else if (s[position_of_assign] == ')')
        {
            position_of_assign = -1;
        }
    }
    string function_name, class_name;
    if (position_of_assign >= 0)
    {
        string argument_name = s.substr(position_of_assign + 2, s.size() - position_of_assign); //+2 epeidi eimai sto <-
        s.erase(s.begin() + position_of_assign - 1, s.end());
        s += "," + argument_name + ")";
        function_name = "<-"; // stin periptosi pou exo to <- tote prepei na to valo kai sto onoma tis sinartisis
        // function(x)<-value
    }
    int position_start_name = 0, position_end_name;

    if (is_method(s))
    {
        position_start_name = s.find('{') + 1;
        position_end_name = s.find('}', position_start_name);
        function_name = s.substr(position_start_name, position_end_name - position_start_name) + function_name;
        position_start_name = s.find('{', position_end_name) + 1;
        position_end_name = s.find('}', position_start_name);
        class_name = s.substr(position_start_name, position_end_name - position_start_name);
        s.erase(s.begin(), s.begin() + position_end_name + 1);
    }
    else
    { // periptosi pou einai function(x)<-value
        position_end_name = s.find('(');
        function_name = s.substr(position_start_name, position_end_name - position_start_name) + function_name;
        s.erase(s.begin(), s.begin() + position_end_name);
        return "\"" + function_name + "\"" + s;
    }

    return is_R_operator(function_name) or position_of_assign >= 0 ? "\"" + function_name + "." + class_name + "\"" + s : function_name + "." + class_name + s;
}

vector<string> read_usage(ifstream &file, vector<string> &aliases_with_lines_more_than_90)
{
    DEBUG("START read_usage");
    vector<string> usg;
    string s;
    bool sinexeia_apo_kato_grammi = false;
    reset_file(file);
    while (!get_usage(file, s))
        ;
    do
    {
        getline(file, s);

        if (s.size() > 90)
        {
            aliases_with_lines_more_than_90.push_back(s);
        }

        remove_spaces(s);
        if (s != "" and sinexeia_apo_kato_grammi)
        {
            DEBUG("ektelesi tin sinexeia stin apo kato grammi");
            sinexeia_apo_kato_grammi = false;
            usg[usg.size() - 1] += s;
        }
        else if (s != "}" and s[s.size() - 1] != '}' and s != "")
        { //  keni grammi, sketo "}", "sinartisi}"
            usg.push_back(s);
        }
        if (s != "" and !find_string(s, "<-") and s[s.size() - 1] != ')')
        { //  periptosi pou i sinartisi paei kai stin kato grammi kai den prepei na einai assign function
            DEBUG("BBrike sinexeia stin apo kato grammi");
            sinexeia_apo_kato_grammi = true;
        }
    } while (s[s.size() - 1] != '}');
    if (s.size() > 1 and s[s.size() - 1] == '}')
    { //  periptosi "sinartisi}"
        s.erase(s.end() - 1);
        usg.push_back(s);
    }
    for (auto &v : usg)
    {
        if (is_method(v) or find_string(v, "<-"))
        {
            v = convert_method_to_function(v);
        }
    }
    // for(auto& v : usg)
    //   name_of_functions_in_usage.push_back(v.substr(0,v.find(')')));
    DEBUG("END read_usage");
    return usg;
}

void remove_spaces(string &s)
{
    s.erase(remove_if(s.begin(), s.end(), [&](char &x)
                      { return isspace(x); }),
            s.end());
}

void remove_spaces_from_begin_end(string &s)
{
    auto it = find_if(s.begin(), s.end(), [&](char &c)
                      { return !isspace(c); });
    s.erase(s.begin(), it);
    int it2 = find_if(s.rbegin(), s.rend(), [&](char &c)
                      { return !isspace(c); }) -
              s.rbegin();
    s.erase(s.end() - it2, s.end());
}

string read_function_from_r_file(ifstream &file)
{
    string func;
    string s;
    size_t bg;
    reset_file(file);
    DEBUG("START read_function_from_r_file");
    do
    {
        getline(file, s);
    } while (s[0] == '#'); // oso briskei sxolia

    DEBUG(s);
    remove_spaces(s);
    func = s;
    if (!find_string(s, "){"))
    { // periptosi pou paei kai se alli grammi i sinartisi
        do
        {
            getline(file, s);
            remove_spaces(s);
            func += s;
        } while (!find_string(s, "{"));
    }
    DEBUG("function: ", func);
    string keyword_function1 = "<-function";
    string keyword_function2 = "=function";
    int keyword_function_size = keyword_function1.size();
    bg = func.find(keyword_function1);
    DEBUG(bg);
    if (bg == string::npos)
    {
        DEBUG("trying operator =.");
        bg = func.find(keyword_function2);
        keyword_function_size = keyword_function2.size();
    }
    func.erase(func.begin() + bg, func.begin() + bg + keyword_function_size);
    DEBUG(func);
    func.erase(func.end() - 1);
    DEBUG(func);
    DEBUG("END read_function_from_r_file");
    return func;
}

bool check_read_file(ifstream &file, char attr)
{
    DEBUG("Start checking file");
    string s;
    bool ret = true;
    while (getline(file, s))
    {
        if (is_dont_read(s, attr))
        {
            ret = false;
            break;
        }
        else if (!isspace(s[0]))
        {
            break;
        }
    }
    DEBUG("End checking file");
    return ret;
}

bool is_dont_read(string &s, char attr)
{
    return s[0] == attr and s.size() >= sizeof("[dontread]") and // sizeof("[dontread]") 11 + \0
           s[1] == '[' and s[2] == 'd' and s[3] == 'o' and s[4] == 'n' and s[5] == 't' and s[6] == ' ' and s[7] == 'r' and
           s[8] == 'e' and s[9] == 'a' and s[10] == 'd' and s[11] == ']';
}

bool is_export(string &s)
{
    return s[0] == '#' and s.size() >= sizeof("[export]") and // sizeof("[export]") 8 + \0
           s[1] == '[' and s[2] == 'e' and s[3] == 'x' and s[4] == 'p' and
           s[5] == 'o' and s[6] == 'r' and s[7] == 't' and s[8] == ']';
}

bool is_export_s3(string &s)
{
    return s[0] == '#' and s.size() >= sizeof("[exports3]") and // sizeof("[exports3]") 8 + \0
           s[1] == '[' and s[2] == 'e' and s[3] == 'x' and s[4] == 'p' and
           s[5] == 'o' and s[6] == 'r' and s[7] == 't' and s[8] == 's' and
           s[9] == '3' and s[10] == ']';
}

bool is_export_special(string &s)
{
    return s[0] == '#' and s.size() >= sizeof("[exportspecial]") and // sizeof("[exportspecial]") 8 + \0
           s[1] == '[' and s[2] == 'e' and s[3] == 'x' and s[4] == 'p' and
           s[5] == 'o' and s[6] == 'r' and s[7] == 't' and s[8] == 's' and
           s[9] == 'p' and s[10] == 'e' and s[11] == 'c' and s[12] == 'i' and
           s[13] == 'a' and s[14] == 'l';
}

bool is_s3method(string &s)
{
    return s.size() >= sizeof("S3method") and s[0] == 'S' and s[1] == '3' and s[2] == 'm' and s[3] == 'e' and s[4] == 't' and s[5] == 'h' and s[6] == 'o' and s[7] == 'd';
}

bool is_R_operator(string s)
{
    return s[0] == '[' or s[0] == ']' or s[0] == '+' or s[0] == '-' or
           s[0] == '&' or s[0] == '|' or s[0] == '/' or s[0] == '!' or
           s == "!=" or s == "==" or s == "*" or s == "and" or
           s == "||";
}

bool is_hidden_function(string &s)
{
    return s.size() > 1 and s[0] == '.';
}

string read_current_signature_function_from_r_file(string &line, string keyword_function, ifstream &file, const int position_of_function_key)
{
    string func = line;
    DEBUG("START read_function_from_r_file");
    remove_spaces(line);
    if (!find_string(line, "){"))
    { // periptosi pou paei kai se alli grammi i sinartisi
        do
        {
            getline(file, line);
            remove_spaces(line);
            func += line;
        } while (!find_string(line, "{"));
        //++depth_scope;
        line = func;
        // ayto edo einai akros xrisimo dioti an to anoigma tou scope den einai mazi me
        // tin ipografi tis sinartisi tote den tha metrithei sosta.
    }
    DEBUG("function: ", func);
    func.erase(func.begin() + position_of_function_key, func.begin() + position_of_function_key + keyword_function.size());
    DEBUG(func);
    func.erase(func.end() - 1);
    DEBUG(func);
    DEBUG("END read_function_from_r_file");
    return func;
}

// proipothesi na einai dilomeno h sinartisi tou stil a<-function h a=function
void read_functions_from_r_file(
    Path& filename,
    const bool full_paths,
    vector<string> &exported_functions_names,
    vector<string> &exported_functions_s3,
    vector<string> &not_exported_functions_names,
    vector<string> &hidden_functions_names,
    vector<string> &exported_special_functions,
    List &signatures,
    bool &found_dont_read)
{
    ifstream file(filename.filename(true));
    size_t position_of_function_key1 = 0, position_of_function_key2 = 0;
    string line;
    int depth_scope = 0, number_of_line = 0, number_of_export_line = 0;
    bool found_export = false, found_export_s3 = false, found_export_special = false;

    while (getline(file, line) and isspace(line[0]))
    {
        ++number_of_line;
    }

    if (is_dont_read(line, '#'))
    {
        found_dont_read = true;
        DEBUG("found dont read: " + number_of_line);
    }
    else
    {
        auto read_name_from_function = [&](string &line)
        {
            string function_name;
            bool read_constant_string = false;
            DEBUG(line);

            for (unsigned int i = 0; i < line.size(); ++i)
            {
                char ch = line[i];

                if (isspace(ch))
                    continue;
                else if (ch == '\"')
                { // kathe fora pou vrisko aytakia tote energopoio-apenergopoio ton mixanismo tou string
                    read_constant_string = !read_constant_string;
                }
                else if (ch == '{' and !read_constant_string)
                { // an eimai mesa se string tote den metrao to {
                    ++depth_scope;
                    DEBUG(to_string(__LINE__) + ": depth-> " + to_string(depth_scope));
                    continue;
                }
                else if (ch == '}' and !read_constant_string)
                { // an eimai mesa se string tote den metrao to }
                    --depth_scope;
                    DEBUG(to_string(__LINE__) + ": depth-> " + to_string(depth_scope));
                    continue;
                }
                if (ch == '<' and i + 9 < line.size())
                { // an bro < kai exo akoma 9 theseis na psakso
                    position_of_function_key1 = line.find("<-function");
                    // an brika to function kai eimai sto global scope
                    if (position_of_function_key1 == i and depth_scope == 0)
                    {
                        if (found_export)
                        {
                            DEBUG("<-function export: " + function_name);
                            exported_functions_names.push_back(function_name);
                            signatures[function_name] = List::create(_["signature"] = read_current_signature_function_from_r_file(line, "<-function", file, position_of_function_key1), _["filename"] = filename.filename(full_paths), _["export type"] = "export");
                            found_export = false;
                        }
                        else if (found_export_s3)
                        {
                            DEBUG("<-function export s3: " + function_name);
                            exported_functions_s3.push_back(function_name);
                            signatures[function_name] = List::create(_["signature"] = read_current_signature_function_from_r_file(line, "<-function", file, position_of_function_key1), _["filename"] = filename.filename(full_paths), _["export type"] = "export s3");
                            found_export_s3 = false;
                        }
                        else if (found_export_special)
                        {
                            DEBUG("<-function export special: " + function_name);
                            exported_special_functions.push_back(function_name);
                            signatures[function_name] = List::create(_["signature"] = read_current_signature_function_from_r_file(line, "<-function", file, position_of_function_key1), _["filename"] = filename.filename(full_paths), _["export type"] = "export special");
                            found_export_special = false;
                        }
                        else if (is_hidden_function(function_name))
                        {
                            hidden_functions_names.push_back(function_name);
                        }
                        else
                        {
                            DEBUG("<-function: " + function_name);
                            not_exported_functions_names.push_back(function_name);
                        }
                    }
                }
                else if (ch == '=' and i + 8 < line.size())
                { // an bro = kai exo akoma 8 theseis na psakso
                    position_of_function_key2 = line.find("=function");
                    // an brika to function kai eimai sto global scope
                    if (position_of_function_key2 == i and depth_scope == 0)
                    {
                        DEBUG("=function: " + function_name);
                        if (found_export)
                        {
                            DEBUG("<-function export: " + function_name);
                            exported_functions_names.push_back(function_name);
                            signatures[function_name] = List::create(_["signature"] = read_current_signature_function_from_r_file(line, "=function", file, position_of_function_key2), _["filename"] = filename.filename(full_paths), _["export type"] = "export");
                            found_export = false;
                        }
                        else if (found_export_s3)
                        {
                            DEBUG("<-function export s3: " + function_name);
                            exported_functions_s3.push_back(function_name);
                            signatures[function_name] = List::create(_["signature"] = read_current_signature_function_from_r_file(line, "=function", file, position_of_function_key2), _["filename"] = filename.filename(full_paths), _["export type"] = "export s3");
                            found_export_s3 = false;
                        }
                        else if (found_export_special)
                        {
                            DEBUG("<-function export s3: " + function_name);
                            exported_special_functions.push_back(function_name);
                            signatures[function_name] = List::create(_["signature"] = read_current_signature_function_from_r_file(line, "=function", file, position_of_function_key2), _["filename"] = filename.filename(full_paths), _["export type"] = "export s3");
                            found_export_special = false;
                        }
                        else if (is_hidden_function(function_name))
                        {
                            hidden_functions_names.push_back(function_name);
                        }
                        else
                        {
                            DEBUG("<-function: " + function_name);
                            not_exported_functions_names.push_back(function_name);
                        }
                    }
                }
                function_name += ch;
            }
            if ((found_export or found_export_s3 or found_export_special) and number_of_export_line == number_of_line - 1)
            { // an exo vrei export  stin porigoumeni grammi kai h epomeni einai space
                Rcout << "Warning: In file '" << filename.filename(false) << "' unused [export] attribute in line " << number_of_export_line << ".\n";
            }
            found_export = false;
            found_export_s3 = false;
            found_export_special = false;
        };

        do
        {
            ++number_of_line;
            remove_spaces(line);
            // mono oi sinartiseis poy einai sto global scope epitrepontai na einai export
            if (is_export(line) and depth_scope == 0)
            {
                DEBUG("found export: in line " + line + "in depth " + to_string(depth_scope));
                number_of_export_line = number_of_line;
                found_export = true;
                continue;
            }
            else if (is_export_s3(line) and depth_scope == 0)
            {
                DEBUG("found export s3: in line " + line + "in depth " + to_string(depth_scope));
                number_of_export_line = number_of_line;
                found_export_s3 = true;
                continue;
            }
            else if (is_export_special(line) and depth_scope == 0)
            {
                DEBUG("found export special: in line " + line + "in depth " + to_string(depth_scope));
                number_of_export_line = number_of_line;
                found_export_special = true;
                continue;
            }
            else if (line[0] == '#')
            { // pass comments
                DEBUG("found comments: in line " + line + "in depth " + to_string(depth_scope));
                continue;
            }

            if (find_string(line, "#"))
            { // an iparxei miksei kodika kai sxolion (} ## aaaa)
                line.erase(line.find("#"));
            }
            read_name_from_function(line);
        } while (getline(file, line));
    }
}

List read_functions_and_signatures(string path,const bool full_paths)
{
    vector<string> exported_functions_names, exported_functions_s3, hidden_functions_names,
        not_exported_functions_names, dont_read, exported_special_functions;
    Files files = readDirectory(path, ".R");
    exported_functions_names.reserve(500);
    exported_functions_s3.reserve(50);
    not_exported_functions_names.reserve(500);
    bool found_dont_read = false;
    List signatures;
    for (auto &file : files)
    {
        read_functions_from_r_file(
            file,
            full_paths,
            exported_functions_names,
            exported_functions_s3,
            not_exported_functions_names,
            hidden_functions_names,
            exported_special_functions,
            signatures,
            found_dont_read);
        if (found_dont_read)
        {
            found_dont_read = false;
            dont_read.push_back(file.filename(full_paths));
        }
    }
    List l, exp;
    if (!dont_read.empty())
        l["dont read"] = dont_read;
    if (!not_exported_functions_names.empty())
        l["without export"] = not_exported_functions_names;
    if (signatures.size() != 0)
        l["signatures"] = signatures;
    if (!hidden_functions_names.empty())
        l["hidden functions"] = hidden_functions_names;
    if (!exported_special_functions.empty())
        exp["special"] = exported_special_functions;
    if (!exported_functions_names.empty())
        exp["functions"] = exported_functions_names;
    if (!exported_functions_s3.empty())
        exp["s3"] = exported_functions_s3;
    else
        exp["s3"] = std::vector<string>();
    if (exp.size() != 0)
        l["export"] = exp;

    return l;
}
