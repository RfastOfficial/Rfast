//Author: Manos Papadakis

#ifndef SYSTEM_FILES
#define SYSTEM_FILES

#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <dirent.h>
#include <chrono>
#include <algorithm>
#include <iterator>
#include <Rcpp.h>

using Rcpp::List;
using std::vector;
using std::array;
using std::string;
using std::ifstream;

class Timer {
    std::chrono::steady_clock sc;
    bool new_measure=false;
    std::chrono::time_point<std::chrono::steady_clock> start,end;
public:
    void Start(){
        if(this->new_measure){
          Rcpp::stop("Error: you haven't stop this measure and the results will be wrong.\n");
        }
        this->new_measure=true;
        this->start=sc.now();
    }
    void Stop(){
        this->end=sc.now();
        if(!this->new_measure){
          Rcpp::stop("Error: you must Start measure first.\n");
        }
        this->new_measure=false;
    }
    double getTime(){
        /*
        string unit;
        long unit_number;
        if(result>=1){
            unit="seconds";
            unit_number=1; //nothing
        }else if(result>=1e-2 and result<1){
            unit="milliseconds";
            unit_number=1000; // 10^3
        }else if(result>=1e-3 and result<10e-2){
            unit="milliseconds";
            unit_number=10000; // 10^4
        }else if(result>=1e-6 and result<10e-3){
            unit="microseconds";
            unit_number=1000000; //10^6
        }else{
            stop("Error: times is weird...\n");
        }
        
        return to_string(result*unit_number)+" "+unit;*/
        return static_cast<std::chrono::duration<double>>(this->end-this->start).count();
    }
};


void print_error();

template<class T,class... Args>
void print_error(T values,Args... args){
	Rcpp::Rcout<<values<<"\n";
	print_error(args...);
}


//#define PRINT_ERRORS
#ifdef PRINT_ERRORS
#define DEBUG print_error
#else
#define DEBUG(...);
#endif


vector<string> split_words(string,const char*);
array<string,2> split_words_in_half(string,const char);
void writeFile(vector<string>,string);
vector<string> readNamespaceFile(string,int&);
bool is_namespace_export(string);
vector<string> readDirectory(const string,const int);
bool is_alias(const char *s,int);
bool next_alias(ifstream &,string &);
vector<string> read_aliases(ifstream &);
vector<string> find_which(vector<string>,vector<string>);
vector<string> find_duplis(vector<string>);
bool is_example(const char *,int);
int get_example(ifstream&,string&);
vector<string> read_directory(string);
string read_example(ifstream &,int&);
bool binary_help(vector<string>::iterator,vector<string>::iterator,string&,vector<string>::iterator&);
vector<string> read_usage(ifstream &, vector<string> &);
string read_function_from_r_file(ifstream &);
void remove_spaces(string&);
void remove_spaces_from_begin_end(string&);
List read_examples(string);
bool check_read_file(ifstream&,char);
void dont_read_man(vector<string>&,vector<string>&);
void reset_file(ifstream& file);

bool is_dont_read(string& s,char attr);
bool is_export(string& s);
string read_current_signature_function_from_r_file(string& ,string ,ifstream &,const int );

void read_functions_from_r_file(
    const string,
    vector<string> &,
    vector<string> &,
    vector<string> &,
    vector<string> &,
    vector<string> &,
    List& ,
    bool& );        

List read_functions_and_signatures(string path);
bool is_export_s3(string&);
bool is_export_special(string&);
bool is_s3method(string&);
bool is_R_operator(string);
bool is_hidden_function(string&);

template<class T>
bool find_string(string& s,T f){
  return s.find(f)!=string::npos;
}

#endif
