/*
 * FataParser.h
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#ifndef FATAPARSER_H_
#define FATAPARSER_H_


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>


using namespace std;


class FastaParser{
public:

	FastaParser(string file);

	vector< string> contig_name;

	~FastaParser(){
		filename.clear();
		filepos.clear();
		delete [] buffer;
	}
	string getChr(int id);

private:
	string filename;
	vector< streampos > filepos;
	size_t buffer_size;
	char*buffer;
	ifstream myfile;


};
#endif /* FATAPARSER_H_ */
