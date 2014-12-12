/*
 * SamParser.cpp
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */
#include "SamParser.h"

SamParser::SamParser(string file) {
	myfile.open(file.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Sam Parser: could not open file: " << file.c_str() << endl;
		exit(0);
	}

	buffer_size = 3000;
	buffer = new char[buffer_size];
	myfile.getline(buffer, buffer_size);

	string match = "SN:";
	while (!myfile.eof() && buffer[0] == '@') {
		for (size_t i = 0; i < buffer_size; i++) {

			if (strncmp(&buffer[i], match.c_str(), 3) == 0) {
				i += 3;	//SN:
				string chr;
				while (buffer[i] != '\t') {
					chr += buffer[i];
					i++;
				}
				chrs.push_back(chr);
			}
		}
		myfile.getline(buffer, buffer_size);
	}
	if (chrs.size() == 0) {
		cerr << "Error no header detected in SAM file which is required!" << endl;
		exit(0);
	}
}

Alignment * SamParser::parseRead(uint16_t mappingQv) {
	//	myfile.getline(buffer, buffer_size);
	Alignment * aln = new Alignment();
	parseReadFast(mappingQv,aln);

	return aln;
}

void SamParser::parseReadFast(uint16_t mappingQv,Alignment *& aln){
	while (aln->getSequence().first.empty() && !myfile.eof()) {
		BamAlignment * bam = new BamAlignment();
		int count = 0;
		string chr;
		string cigar;

		size_t i = 0;

		while (i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n') {
			switch (count) {
			case 0:
				if (buffer[i] != '\t') {
					bam->Name += buffer[i];
				}
				break;
			case 1:
				if (buffer[i - 1] == '\t') {
					int flag = atoi(&buffer[i]);
					bam->SetIsDuplicate(!(flag & 0x400));
					//					Sets value of "PCR duplicate" flag to ok.
					bam->SetIsFailedQC(!(flag & 0x200));
					//					Sets "failed quality control" flag to ok.
					bam->SetIsFirstMate(!(flag & 0x40));
					//					Sets "alignment is first mate" flag to ok.
					bam->SetIsMapped(!(flag & 0x4));
					//					Sets "alignment is mapped" flag to ok.
					bam->SetIsMateMapped(!(flag & 0x8));
					//					Sets "alignment's mate is mapped" flag to ok.
					bam->SetIsMateReverseStrand(!(flag & 0x20));
					//					Sets "alignment's mate mapped to reverse strand" flag to ok.
					bam->SetIsPaired(!(flag & 0x1));
					//					Sets "alignment part of paired-end read" flag to ok.
					bam->SetIsPrimaryAlignment(!(flag & 0x100));
					//					Sets "position is primary alignment" flag to ok.
					bam->SetIsProperPair(!(flag & 0x2));
					//					Sets "alignment is part of read that satisfied paired-end resolution" flag to ok.
					bam->SetIsReverseStrand(!(flag & 0x10));
					//					Sets "alignment mapped to reverse strand" flag to ok.
					bam->SetIsSecondMate((flag & 0x40));
					//					Sets "alignment is second mate on read" flag to ok.
					//					bam->SetIsMateUnmapped ( !(flag &0x8));
					////					Complement of using SetIsMateMapped().
					//					bam->SetIsSecondaryAlignment ( (flag &0x100));
					////					Complement of using SetIsPrimaryAlignment().
					//					bam->SetIsUnmapped ( !(flag &0x4));
					//					Complement of using SetIsMapped().
				}
				break;
			case 2:
				if (buffer[i] != '\t') {
					chr += buffer[i];
				}
				break;
			case 3:
				if (buffer[i - 1] == '\t') {
					bam->Position = atoi(&buffer[i]);
				}
				break;
			case 4:
				if (buffer[i - 1] == '\t') {
					bam->MapQuality = atoi(&buffer[i]);
				}
				break;
			case 5:
				if (buffer[i] != '\t') {
					cigar += buffer[i];
				}
				break;
			case 9:
				if (buffer[i] != '\t') {
					bam->QueryBases += toupper(buffer[i]);
				}
				break;
			case 10:
				if (buffer[i] != '\t') {
					bam->Qualities += buffer[i];
				}
				break;
				//			default:
				//				cout<<"else"<<i<<" "<<buffer_size<<endl;;
				//				bam->TagData+=buffer[i];
			}

			if (buffer[i] == '\t') {
				count++;
			}
			i++;
		}
		if (bam->IsMapped()){ //) && bam->MapQuality >= mappingQv) {
			bam->Length = bam->QueryBases.size();
			bam->RefID = getRefID(chr);
			if(bam->RefID == -1){
				cerr<<"No RefID assigned"<<endl;
			}
			bam->CigarData = convertCigar(cigar);
			aln->setAlignment(bam);
		} else {
			delete bam;
		}
		chr.clear();
		cigar.clear();
		myfile.getline(buffer, buffer_size);
	}
}

int SamParser::getRefID(string chr) {

	for (size_t i = 0; i < chrs.size(); i++) {
		if (strcmp(chrs[i].c_str(), chr.c_str()) == 0) {
			return (int) i;
		}
	}
	return -1;
}
vector<CigarOp> SamParser::convertCigar(string cigar) {

	int num = -1;
	vector<CigarOp> tmp;
	for (size_t i = 0; i < cigar.size(); i++) {
		if (num == -1) {
			num = atoi(&cigar[i]);
		} else if (num != -1 && atoi(&cigar[i]) == 0 && cigar[i] != '0') {
			CigarOp op;
			op.Length = num;
			op.Type = cigar[i];
			tmp.push_back(op);
			num = -1;
		}
	}
	return tmp;
}

bool SamParser::SetRegion(const int & RefId, const int & leftPos, const int & rightPos){
    return false;
}

