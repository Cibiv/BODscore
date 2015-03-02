/*
 * Parser.h
 *
 *  Created on: May 29, 2012
 *      Author: fritz
 */

#ifndef PARSER_H_
#define PARSER_H_

#include "Alignment.h"

class Parser {

public:
	virtual ~Parser(){};
	virtual Alignment * parseRead(uint16_t mappingQv) = 0;
	virtual void parseReadFast(uint16_t mappingQv,Alignment *& align) = 0;
    virtual bool SetRegion(const int & RefId, const int & leftPos, const int & rightPos) = 0;
    virtual  string get_header() = 0 ;
    virtual int GetReferenceID(const string& refName) = 0;
//    virtual  RefVector get_refInfo();
};

#endif /* PARSER_H_ */
