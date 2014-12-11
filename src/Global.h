/*
 * global.h
 *
 *  Created on: May 30, 2012
 *      Author: fritz
 */

#ifndef GLOBAL_H_
#define GLOBAL_H_

#include "string.h"
#include "vector"
#include <iostream>


#include "tclap/CmdLine.h"
#include "sambamio/BamParser.h"
#include "sambamio/FastaParser.h"
#include "sambamio/Parser.h"
#include "sambamio/SamParser.h"
#include "sambamio/Alignment.h"


#include "BasisClass.h"
using namespace TCLAP;
using namespace std;
int const step=5;
int const MAX_NUM_ALIGNMENTS = 5000000;

//string convertInt(int number) {
//	if (number == 0)
//		return "0";
//	string temp = "";
//	string returnvalue = "";
//	while (number > 0) {
//		temp += number % 10 + 48;
//		number /= 10;
//	}
//	for (int i = 0; i < (int) temp.length(); i++)
//		returnvalue += temp[temp.length() - i - 1];
//	return returnvalue;
//}

#endif /* GLOBAL_H_ */
