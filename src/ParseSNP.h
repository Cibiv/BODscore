/*
 * ParseSNP.h
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#ifndef PARSESNP_H_
#define PARSESNP_H_

#include <regex>

#include <string>
#include <vector>
#include <map>
#include <assert.h>     /* assert */
#include "Global.h"
#include <sqlite3pp.h>
#include "SqliteDb.h"
#include <memory>

size_t const range =200;
// int const step=5;
using namespace std;

template <class T>
void deleteInVector(vector<T*>* deleteme) {
     while(!deleteme->empty()) {
          delete deleteme->back();
          deleteme->pop_back();
     }

     delete deleteme;
}
// ================================================================================
class ParseSNP: public BasisClass{
private:
    bool IMMEDIATE = false;
    int verbose = 0;
    string snpfile;
    string read_filename;
    string reffile;
    string output;
    string plot_file;
    string db_file;
    string sample_label;
    bool db_flag = false;
    // unique_ptr<sqlite3pp::database> 
    SqliteDb * db;
    int read_length;
    int range;
    int num_test = 0; // for testing purposes: maximal number of snps per contig to process.
    int n_snp = 0;

//    vector<vector<Coverage*> > genome;
//    map<string,Coverage*> covs;
    map<string,int>chromosome_vector;

    vector<string> exclude_contigs;//    vector<vector<int> >  snps;

    char * buffer;
    size_t buffer_size;
    bool process_snp(Coverage* cov, string & ref, Parser * mapped_file, const size_t &cc, const size_t & chr_bam);
    void print();
    void parseVCF();
    void init();

public:
	ParseSNP(){
		read_length=0;
		buffer_size = 3000;
		buffer = new char[buffer_size];
	}
	~ParseSNP(){
		delete buffer;
	}
	void parse_cmd_line(int argc, char *argv[]);
    void compute();
};



#endif /* PARSESNP_H_ */
