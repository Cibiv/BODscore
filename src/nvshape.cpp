//============================================================================
// Name        : VShape.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : artistic license
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "BasisClass.h"
#include "ParseSNP.h"
// #include "Ident_hist.h"


#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {

//	if (argc < 1) {
//			cout << "choose" << endl;
//			cout<<"1: summarize Identity"<<endl;
//			cout<<"2: Filter SNPs"<<endl;
//			exit(0);
//		}

    // BasisClass * tool =NULL;
    ParseSNP * tool = new ParseSNP();
	tool->parse_cmd_line(argc, &argv[0]);
//	tool->compute_jump();

	delete tool;
	return 0;
}
