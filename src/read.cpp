
//Author: Manos Papadakis
#include <RcppArmadillo.h>
#include <dirent.h>
#include "system_files.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace Rcpp;

using std::vector;
using std::string;


vector<string> read_directory(string path){
  DIR *dir=NULL;
  struct dirent *ent;
  vector<string> files;
  string textf;
  if((dir = opendir(path.c_str())) != NULL) {
    readdir(dir);
    readdir(dir);
    while((ent = readdir(dir)) != NULL) {
      textf=ent->d_name;
      files.push_back(textf);
    }
    closedir(dir);
  }else{
    stop("Error: Could not open directory with path \""+path+"\"");
  }
  return files;
}


RcppExport SEXP Rfast_read_directory(SEXP pathSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const string >::type path(pathSEXP);
    __result = read_directory(path);
    return __result;
END_RCPP
}


/////////////////////////////////////////////////////////////////////////

int is_regular_file(const char *path)
{
    struct stat path_stat;
    stat(path, &path_stat);
    return S_ISREG(path_stat.st_mode);
}

using std::ifstream;
using std::vector;
using std::string;

List read_examples(string path_man){
  ifstream file;
  vector<string> examples,all_rd_files=read_directory(path_man),files_long_lines,dontread_rd;
  string tmp;
  int longlines=0;
  for(unsigned int i=0;i<all_rd_files.size();++i){
  	string filename = path_man+all_rd_files[i];
  	if(is_regular_file(filename.c_str())){
	    file.open(filename);
	    if(!file.is_open()){
	      stop("Can't open file \"%s\".",all_rd_files[i]);
	    }
	    if(check_read_file(file,'%')){
	      longlines=0;
	      tmp=read_example(file,longlines);
	      if(longlines){
	      	files_long_lines.push_back(all_rd_files[i]);
	      }
	      if(!tmp.empty())
	        examples.push_back(tmp);
	    }else{
	      DEBUG("Find attribute dont read file with name: "+all_rd_files[i]);
	      dontread_rd.push_back(all_rd_files[i]);
	      all_rd_files.erase(all_rd_files.begin()+i);
	      --i;
	    }
	    file.close();
	}
  }
  List l;
  l["examples"]=examples;
  l["files"]=all_rd_files;
  l["long_lines"]=files_long_lines;
  l["dont read"]=List::create(_["Rd"]=dontread_rd);
  return l;
}

RcppExport SEXP Rfast_read_examples(SEXP path_manSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< string >::type path_man(path_manSEXP);
    __result = read_examples(path_man);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////
