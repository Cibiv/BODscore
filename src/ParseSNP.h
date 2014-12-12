/*
 * ParseSNP.h
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#ifndef PARSESNP_H_
#define PARSESNP_H_

#include <string>
#include <vector>
#include <map>
#include <assert.h>     /* assert */
#include "Global.h"

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
    string snpfile;
    string read_filename;
    string reffile;
    string output;
    string plot_file;
    int read_length;
    int range;

    vector<vector<Coverage*> > genome;
    map<string,Coverage*> covs;
    map<string,int>chrs;

    vector<vector<int> >  snps;

    char * buffer;
    size_t buffer_size;
    void process_snp(Coverage* cov, string & ref, Parser * mapped_file, const size_t &cc);
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
