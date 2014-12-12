//============================================================================
// Name        : VShape.cpp
// Author      : Fritz Sedlazeck
// Version     :
// Copyright   : artistic license
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "BasisClass.h"
#include "ParseSNP.h"
#include "Ident_hist.h"


#include <iostream>
using namespace std;

int main(int argc, char *argv[]) {

	if (argc < 2) {
			cout << "choose" << endl;
			cout<<"1: summarize Identity"<<endl;
			cout<<"2: Filter SNPs"<<endl;
			exit(0);
		}

		BasisClass * tool =NULL;
		switch(atoi(argv[1])){
		case 1:
			tool=new Ident_hist();
			tool->parse_cmd_line(argc-1,&argv[1]);
			tool->compute();
			break;
		case 2:
			tool=new ParseSNP();
			tool->parse_cmd_line(argc-1,&argv[1]);
			tool->compute();
			break;
		default:
			cerr<<"No available option"<<endl;
			break;
		}



		delete tool;
		return 0;
}
