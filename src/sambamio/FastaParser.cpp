/*
 * FastaParser.cpp
 *
 *  Created on: May 25, 2012
 *      Author: fritz
 */

#include "FastaParser.h"


FastaParser::FastaParser(string file){
	filename=file;
	myfile.open(filename.c_str(), ifstream::in);
	if (!myfile.good()) {
		cout << "Fasta Parser: could not open file: " << filename.c_str() << endl;
		exit(0);
	}

	buffer_size = 2000;
	buffer = new char[buffer_size];

	myfile.getline(buffer, buffer_size);

	while(!myfile.eof()){
		if(buffer[0]=='>'){
			filepos.push_back(myfile.tellg());
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
}

string FastaParser::getChr(int id){
	if((size_t)id > filepos.size()){
		cerr<<"Wrong number of Chrs! Maybe you choose the wrong Ref file"<<endl;
		exit(0);
	}
	myfile.open(filename.c_str(), ifstream::in);
	myfile.seekg(filepos[id]);
	myfile.getline(buffer, buffer_size);
	string seq;
	while(!myfile.eof() && buffer[0]!='>'){
		for(size_t i=0;i<buffer_size && buffer[i]!='\0'&&buffer[i]!='\n'&& buffer[0]!='>';i++){
			seq+=toupper(buffer[i]);
		}
		myfile.getline(buffer, buffer_size);
	}
	myfile.close();
	return seq;
}
