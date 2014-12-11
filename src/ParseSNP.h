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

#include "Global.h"

//int const range =200;
//int const step=5;
using namespace std;

struct coverage{
	int start_pos;
	int*cov_100;
	int*cov_90;
	int chr;
	float score;
	string seq;
};

class ParseSNP: public BasisClass{
private:
	string snpfile;
	string read_filename;
	string reffile;
	string output;
	string plot_file;
	int read_length;
	int range;

	map<string,coverage*> covs;
	map<string,int>chrs;

	vector<vector<int> >  snps;

	char * buffer;
	size_t buffer_size;
	void process_alignemnts(vector< vector<Alignment *> > &alignments,FastaParser * fasta);
	void print();
	void compute_cov(Alignment * aln, coverage* cov);
	void parseVCF();
	void init();
	float conf_val(int start_cov,int stop_cov,int max_cov);
	float estimate(int * cov_100,int * cov_90, int pos);
public:
	ParseSNP(){
		read_length=0;
		buffer_size = 3000;
		buffer = new char[buffer_size];
		range=0;
	}
	~ParseSNP(){
		delete [] buffer;
		covs.clear();
		snps.clear();
	}
	void parse_cmd_line(int argc, char *argv[]);
	void compute();

};



#endif /* PARSESNP_H_ */
