/*
 * Alignments.h
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#ifndef ALIGNMENTS_H_
#define ALIGNMENTS_H_

#include <string.h>
#include <vector>
#include "api/BamAux.h"
#include "api/BamMultiReader.h"
#include "api/BamWriter.h"
using namespace BamTools;
using namespace std;

class Alignment {

private:
	BamAlignment * al;
	pair<string,string> alignment;
	bool is_computed;
	int32_t orig_length;
//	static int num;
public:
	Alignment(){
//		num++;
		orig_length=0;
		al=NULL;
		is_computed=false;
	}
	~Alignment(){
//		num--;
//		cout<<"DEL Align "<<num<<endl;
		alignment.first.clear();
		alignment.second.clear();
		delete al;
	}
	void setAlignment(BamAlignment * al);
	void setRef(string sequence);
	void computeAlignment();

	pair<string,string> getSequence();
	int32_t getPosition();
	int32_t getRefID();
	bool getStrand();
	uint16_t getMappingQual();
	string getName();
	vector<CigarOp> getCigar();
	string getQualitValues();
	size_t getRefLength();
	size_t getOrigLen();
	BamAlignment * getAlignment();
	float getIdentity();
	void initAlignment();
	int getAlignmentFlag();
	string getQueryBases();
	string getQualities();
	string getTagData();
};


#endif /* ALIGNMENTS_H_ */
