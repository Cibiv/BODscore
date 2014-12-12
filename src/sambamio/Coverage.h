/*
 * Coverage.h
 *
 *  Created on: Dec 4, 2014
 *      Author: dima
 */

#ifndef COVERAGE_H_
#define COVERAGE_H_

#include <string.h>
#include <vector>
using namespace std;

class Coverage{

private:
//	static int num;
public:
     int start_pos;
     int*cov_100;
     int chr;
     float score;
     int pos;
     int*cov_90;

    Coverage(int &x, int &r): pos(x), range(r)
    {
        cov_100 = new int[range * 2];
        cov_90 = new int[range * 2];
        memset(cov_100, 0, range * 2 * sizeof(int));
        memset(cov_90, 0, range * 2 * sizeof(int));
	}
    ~Coverage(){
//		num--;
//		cout<<"DEL Align "<<num<<endl;
	    delete cov_90;
        delete cov_100;
	}

};


#endif /* ALIGNMENTS_H_ */
