#ifndef Coverage_H_
#define Coverage_H_

#include<string>
#include<Alignment.h>
#include <stdexcept>
#include <sqlite3pp.h>
#include "QdrArray.h"
// #include <stdexcept>

using namespace std;

int const step=5;
/*
size_t const LO = 0;
size_t const HI = 1;
size_t const FW = 0;
size_t const RV = 1;
*/

class Coverage{
private:
    float conf_val(const int & start_cov, const int & stop_cov, const int & max_cov);
    int buf_len;

    void allocate_arrays(){
        tot_cov     = new QdrArray( range * 2 );
        snp_cov     = new QdrArray( range * 2 );
        aln_centres = new QdrArray( range * 2 );
        // clog << "snp_cov[0][0][1] : " << (*snp_cov)[0][0][1] << endl;
    }
    static const int ZERO;

    const int & max(const int & a, const int & b){   return (a>b ? a : b);   }
    const int & min(const int & a, const int & b){   return (a<b ? a : b);   }
    float mean(const size_t & x_t, const size_t & t, const float & old_mean){
        float ft = (float) t;
        return ( ft-1 ) / ft * old_mean + 1/ft * (float) x_t;
    }

    size_t end_pos(){  return (size_t) start_pos + range * 2;   }
    
    void (Coverage::*pb)(FILE *file, const int * var) = NULL;

    void (Coverage::*pbdb_int)(int *outstr, const int * var) = NULL;

    void print_block(FILE *file, const int * var );
    void print_char_block(FILE *file, const int * var );

public:
    std::string contig;
    int pos ;
    int start_pos ; 
//    int end_pos;
    QdrArray * tot_cov;
    QdrArray * snp_cov;
    QdrArray * aln_centres;
 
    int chr;
	float score;
    int const range;

    int altCounts = 0;
    int totCounts = 0;
    float snp_ratio = 0;

    Coverage(): range(1), pos(0), start_pos(0){ 
        allocate_arrays();
    };
    
    Coverage(int &x, int &r): range(r), pos(x),  start_pos(  ( x-r > 0) ? (x-r) : 0 ) {
        allocate_arrays();
     }

    Coverage(const Coverage &obj): range(obj.range) {
        pos = obj.pos ;
        start_pos = obj.start_pos;
//        tot_cov[0][0] = obj.tot_cov[0][0];
//        snp_cov[0][0] = obj.snp_cov[0][0];
//        aln_centres[0][0] = obj.aln_centres[0][0];
        clear_arrays(); // memset all to zeros
    }

     Coverage(const int &r):  range(r) {
        allocate_arrays();
        // clog << "coverage arrays allocated for " << range*2<< " ints" << endl;
     }

    ~Coverage(){ 
      delete [] snp_cov;
      delete [] tot_cov;
      delete [] aln_centres;
    }

    void clear_arrays(){ // memset all to zeros
        snp_cov->memset_subarrays( );
        tot_cov->memset_subarrays( );
        aln_centres->memset_subarrays( );
    }

    float estimate( int & read_length );
    void compute_cov( Alignment * aln );
    void print_cov( const int & cc, FILE *file );
    bool within( const int & x );

    void reset(std::string current_chr, int x) //: contig(current_chr), pos(x)
    {
        contig = current_chr;
        pos = x;
        start_pos = pos > range ? pos - range : 0;
        // clog << "snp_cov[0][0][1] : " << (*snp_cov)[0][0][1] << endl;
        clear_arrays();
    }
    void reset(int cc, int x)//:  chr(cc), pos(x)
    {
        chr = cc;
        pos = x;
        start_pos = pos > range ? pos - range : 0;
        clear_arrays();
    }
    void reset(std::string current_chr, int cc, int x) //: contig(current_chr), pos(x)
    {
        contig = current_chr;
        chr = cc;
        pos = x;
        start_pos = pos > range ? pos - range : 0;
        clear_arrays();
    }

};

#endif // Coverage_H_
