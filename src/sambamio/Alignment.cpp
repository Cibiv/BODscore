/*
 * Alignments.cpp
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#include "Alignment.h"
//int Alignment::num = 0;
void Alignment::setRef(string sequence) {
	alignment.second = sequence;
}

void Alignment::initAlignment() {
	//al = new BamAlignment();
}

void Alignment::setAlignment(BamAlignment * align) {
	al = align;
	alignment.first.clear();
	alignment.second.clear();
	is_computed = false;

	orig_length = al->QueryBases.size();
	for (size_t i = 0; i < al->QueryBases.size(); i++) {
		alignment.first += toupper(al->QueryBases[i]);
	}
}

void Alignment::computeAlignment() {
	int pos = 0;
	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'I') {
			for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
				alignment.second.insert(pos, "-");
				alignment.second.erase(alignment.second.size() - 1, 1);
				pos++;
			}
		} else if (al->CigarData[i].Type == 'D') {
			for (uint32_t t = 0; t < al->CigarData[i].Length; t++) {
				alignment.first.insert(pos, "-");
				pos++;
			}
		} else if (al->CigarData[i].Type == 'S') {

			if (pos == 0) { //front side
				alignment.second.erase(((int) alignment.second.size()) - al->CigarData[i].Length,
						al->CigarData[i].Length);
			} else { //backside
				alignment.second.erase(pos, al->CigarData[i].Length);
			}
			alignment.first.erase(pos, al->CigarData[i].Length);

		} else if (al->CigarData[i].Type == 'M') {
			pos += al->CigarData[i].Length;
		} else if (al->CigarData[i].Type == 'H') {

		}else if (al->CigarData[i].Type == 'N') {
			alignment.second.erase(pos, al->CigarData[i].Length);
		}

	}
	for (size_t i = 0; i < alignment.first.size(); i++) {
		if (alignment.first[i] == '=') {
			alignment.first[i] = alignment.second[i];
		}
	}

	is_computed = true;

	if (alignment.first.size() != alignment.second.size()) {
		cerr << "Error alignment has different length" << endl;
		cerr << " ignoring alignment " << al->Name << endl;
		cerr << al->Position << endl;

		cerr << endl;
		cerr << "read: " << alignment.first << endl;
		cerr << endl;
		cerr << " ref: " << alignment.second << endl;
		cerr << endl;
		cerr << orig_length << endl;
		vector<CigarOp> cig = getCigar();

		for (size_t i = 0; i < cig.size(); i++) {
			cerr << cig[i].Length << cig[i].Type << " ";
		}
		cerr << endl;
		exit(0);
		return;
	}

}

void Alignment::processAlignment(string &ref){ // by dima

int32_t al_pos = this->getPosition();
string ref_al;

if (al_pos + getRefLength() >= ref.size()) {
           ref_al = ref.substr((size_t) al_pos,  (int32_t) ref.size() - al_pos);

            while (ref_al.size() < this->getRefLength()) { ref_al += 'N'; }
            } else {
                ref_al = ref.substr((size_t) al_pos, this->getRefLength());
            }
// set stuff
this->setRef(ref_al);
this->computeAlignment();
}


int32_t Alignment::getPosition() {
	return al->Position;
}
int32_t Alignment::getRefID() {
	return al->RefID;
}
bool Alignment::getStrand() {
	return !al->IsReverseStrand();
}
vector<CigarOp> Alignment::getCigar() {
	return al->CigarData;
}
string Alignment::getQualitValues() {
	return al->Qualities;
}

size_t Alignment::getRefLength() {
	size_t len = orig_length;

	for (size_t i = 0; i < al->CigarData.size(); i++) {
		if (al->CigarData[i].Type == 'D') {
			len += al->CigarData[i].Length;
		}else if(al->CigarData[i].Type == 'N'){
			len += al->CigarData[i].Length;
		}
	}
	return len;
}

size_t Alignment::getOrigLen() {
	return orig_length;
}
pair<string, string> Alignment::getSequence() {
	return alignment;
}
BamAlignment * Alignment::getAlignment() {
	return al;
}
string Alignment::getName() {
	return al->Name;
}
uint16_t Alignment::getMappingQual() {
	return al->MapQuality;
}
float Alignment::getIdentity() {
	if (is_computed) {
		float match = 0;
		for (size_t i = 0; i < alignment.first.size(); i++) {
			if (alignment.first[i] == alignment.second[i]) {
				match++;
			}
		}
		return match / (float) alignment.first.size();
	}
	return -1;
}

int Alignment::getAlignmentFlag() {
	return al->AlignmentFlag;
}
string Alignment::getQueryBases() {
	return al->QueryBases;
}

string Alignment::getQualities() {
	return al->Qualities;
}

string convertInt(int number){
	stringstream ss; //create a stringstream
	ss << number; //add number to the stream
	return ss.str(); //return a string with the contents of the stream
}

string Alignment::getTagData() {
	vector<string> tags;

	uint32_t i = 0;
	if (al->GetTag("AS", i)) {
		string tmp = "AS:i:";
		tmp += convertInt(i);
		tags.push_back(tmp);

	}
	i = 0;
	if (al->GetTag("NM", i)) {
		string tmp = "NM:i:";
		tmp += convertInt(i);
		tags.push_back(tmp);
	}

	string md;
	if (al->GetTag("MD", md)) {
		string tmp = "MD:Z:";
		tmp += md;
		tags.push_back(tmp);
	}

	i = 0;
	if (al->GetTag("UQ", i)) {
		string tmp = "UQ:i:";
		tmp += convertInt(i);
		tags.push_back(tmp);
	}

	string res;
	for(size_t i=0;i<tags.size();i++){
		res+=tags[i];
		if(i+1<tags.size()){
			res+='\t';
		}
	}
	return res;
}
