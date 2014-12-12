/*
 * SamParser.h
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#ifndef SAMPARSER_H_
#define SAMPARSER_H_


#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include "api/BamAux.h"
#include "Alignment.h"
#include "Parser.h"


using namespace BamTools;
using namespace std;


class SamParser : public Parser{
private:
	size_t buffer_size;
	char*buffer;
	ifstream myfile;
	vector<string> chrs;
	vector<CigarOp> convertCigar(string cigar);
//	string convertRead(vector<CigarOp>cigar, string read);
	int getRefID(string chr);
public:
	SamParser(string file);
	~SamParser(){
		delete [] buffer;
		myfile.close();
		chrs.clear();
	}
	Alignment * parseRead(uint16_t mappingQv);
	void parseReadFast(uint16_t mappingQv,Alignment *& align);
    bool SetRegion(const int & RefId, const int & leftPos, const int & rightPos);

};
#endif /* SAMPARSER_H_ */
