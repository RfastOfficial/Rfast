
//Author: Manos Papadakis
#include <RcppArmadillo.h>
#include "system_files.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace Rcpp;

using std::ifstream;
using std::vector;
using std::string;

List read_examples(string path_man, const bool full_paths = false){
  ifstream file;
  vector<string> examples,files_long_lines,dontread_rd;
  Files all_rd_files=readDirectory(path_man);
  string tmp;
  int longlines=0,i=0;
  for(auto &rd_file : all_rd_files){
	    file.open(rd_file.filename(true));
	    if(!file.is_open()){
	      stop("Can't open file \"%s\".",rd_file.filename(full_paths));
	    }
	    if(check_read_file(file,'%')){
	      longlines=0;
	      tmp=read_example(file,longlines);
	      if(longlines){
	      	files_long_lines.push_back(rd_file.filename(full_paths));
	      }
	      if(!tmp.empty())
	        examples.push_back(tmp);
	    }else{
	      DEBUG("Find attribute dont read file with name: "+rd_file.filename(full_paths));
	      dontread_rd.push_back(rd_file.filename(full_paths));
	      all_rd_files.erase(all_rd_files.begin()+i);
	      --i;
	    }
      ++i;
	    file.close();
  }
  List l;
  if(!examples.empty())
    l["examples"]=examples;
  if(!all_rd_files.empty())
    l["files"]=all_rd_files;
  if(!files_long_lines.empty())
    l["long_lines"]=files_long_lines;
  if(!dontread_rd.empty())
    l["dont read"]=List::create(_["Rd"]=dontread_rd);
  return l;
}

RcppExport SEXP Rfast_read_examples(SEXP path_manSEXP, SEXP full_pathsSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< string >::type path_man(path_manSEXP);
	  traits::input_parameter<const bool>::type full_paths(full_pathsSEXP);
    __result = read_examples(path_man,full_paths);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////
