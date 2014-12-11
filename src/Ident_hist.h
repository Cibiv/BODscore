/*
 * Indent_hist.h
 *
 *  Created on: Sep 27, 2012
 *      Author: fritz
 */

#ifndef INDENT_HIST_H_
#define INDENT_HIST_H_

#include <string>
#include <vector>

#include "math.h"
#include "Global.h"

using namespace std;



class Ident_hist: public BasisClass{
private:
	string read_filename;
	string reffile;
	string output;
	vector<int >  ident;


	void process_alignemnts(vector< vector<Alignment *> > &alignments,FastaParser * fasta);
	void print();
	void compute_cov(Alignment * aln);
public:
	Ident_hist(){
		ident.resize(101,0);
	}
	~Ident_hist(){

	}
	void parse_cmd_line(int argc, char *argv[]);
	void compute();

};




#endif /* INDENT_HIST_H_ */
