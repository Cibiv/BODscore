/*
 * parser.cpp
 *
 *  Created on: Apr 17, 2012
 *      Author: fritz
 */

#include "BamParser.h"

BamParser::BamParser(string file){
	vector<string > tmps;
	tmps.push_back(file);

	if(!reader.Open(tmps)){
		cerr<<"BAM Parser: could not open file: "<<file<<endl;
		exit(0);
	}

   tmps[0] = file + ".bai";
   if(!reader.OpenIndexes(tmps)){
        cerr<<"BAM Parser: could not open index for : "<<file<<endl;
    //    exit(0);
    }

}

bool BamParser::SetRegion(const int & RefId, const int & leftPos, const int & rightPos){
    return reader.SetRegion(RefId, leftPos, RefId, rightPos );
}

Alignment* BamParser::parseRead(uint16_t mappingQv){

	Alignment *align = new Alignment();
	BamAlignment* al = new BamAlignment();
	while(reader.GetNextAlignmentCore(al[0])){
		if( al->IsMapped() && al->MapQuality >= mappingQv){
			al->BuildCharData();
			align->setAlignment(al);
			return align;
		}
	}
	return align;

}
void BamParser::parseReadFast(uint16_t mappingQv,Alignment*& align){

//	Alignment *align = new Alignment();
	BamAlignment* al = align->getAlignment();

	while(reader.GetNextAlignmentCore(al[0])){
		if( al->IsMapped() && al->MapQuality > mappingQv){
			al->BuildCharData();
			align->setAlignment(al);
			return;
		}
	}
}
RefVector BamParser::get_refInfo(){
	return reader.GetReferenceData();
}

string BamParser::get_header(){
	return reader.GetHeaderText();
}

int BamParser::GetReferenceID(const string& refName){
    return reader.GetReferenceID( refName );
}
