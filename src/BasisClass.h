/*
 * BasisClass.h
 *
 *  Created on: May 30, 2012
 *      Author: fritz
 */

#ifndef BASISCLASS_H_
#define BASISCLASS_H_

#include "Global.h"

class BasisClass {

protected:
	string version;

public:
	virtual ~BasisClass(){};
	virtual void parse_cmd_line(int argc, char *argv[])=0;
	virtual void compute()=0;
};
#endif /* BASISCLASS_H_ */
