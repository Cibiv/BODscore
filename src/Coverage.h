#ifndef Coverage_H_
#define Coverage_H_

#include<string>
#include<Alignment.h>
#include <stdexcept>

using namespace std;

int const step=5;
size_t const HI = 0;
size_t const LO = 1;
size_t const FW = 0;
size_t const RV = 1;

class Coverage{
private:
    float conf_val(const int & start_cov, const int & stop_cov, const int & max_cov);

    void allocate_arrays(){
        tot_cov[0][0] = new int[range * 2];
        tot_cov[0][1] = new int[range * 2];
        tot_cov[1][0] = new int[range * 2];
        tot_cov[1][1] = new int[range * 2];
        
        snp_cov[0][0] = new int[range * 2];
        snp_cov[0][1] = new int[range * 2];
        snp_cov[1][0] = new int[range * 2];
        snp_cov[1][1] = new int[range * 2];
 
        cov_lo_f = new int[range * 2];
        cov_lo_r = new int[range * 2];
        cov_hi_f = new int[range * 2];
        cov_hi_r = new int[range * 2];

        aln_centres[HI][FW] = new int[range * 2];
        aln_centres[HI][RV] = new int[range * 2];
        aln_centres[LO][FW] = new int[range * 2];
        aln_centres[LO][RV] = new int[range * 2];
        
        snp_cov_lo_f = new int[range * 2];
        snp_cov_lo_r = new int[range * 2];
        snp_cov_hi_f = new int[range * 2];
        snp_cov_hi_r = new int[range * 2];
    
        clear_arrays(); // memset all to zeros;
    }
    static const int ZERO;

    const int & max(const int & a, const int & b){   return (a>b ? a : b);   }
    const int & min(const int & a, const int & b){   return (a<b ? a : b);   }
    float mean(const size_t & x_t, const size_t & t, const float & old_mean){
        float ft = (float) t;
        return ( ft-1 ) / ft * old_mean + 1/ft * (float) x_t;
    }

    size_t end_pos(){  return (size_t) start_pos + range * 2;   }
    void print_block(FILE *file, const int * var );

public:
    int pos ;
    int start_pos ; 
//    int end_pos;
    int* cov_hi_f;
    int* cov_hi_r;
	int* cov_lo_f;
	int* cov_lo_r;
    
    int * tot_cov[2][2];
    int * snp_cov[2][2];
    int * aln_centres[2][2];
    
    int* snp_cov_lo_f;
    int* snp_cov_lo_r;
    int* snp_cov_hi_f;
    int* snp_cov_hi_r;
    int chr;
	float score;
    int const range;

    Coverage(): range(1), pos(0), start_pos(0){ 
        allocate_arrays();
        clog << "allocated for " << range*2<< " ints" << endl;
    };
    
    Coverage(int &x, int &r): range(r), pos(x),  start_pos(  ( x-r > 0) ? (x-r) : 0 ) {
        allocate_arrays();
     }

    Coverage(const Coverage &obj): range(obj.range) {
        pos = obj.pos ;
        start_pos = obj.start_pos;
        cov_lo_f = obj.cov_lo_f;
        cov_lo_r = obj.cov_lo_r;
        cov_hi_f = obj.cov_hi_f;
        cov_hi_r = obj.cov_hi_r;
        clear_arrays(); // memset all to zeros
    }

     Coverage(const int &r):  range(r) {
        allocate_arrays();
     }

    void delete_subarrays(int *p[2][2] ){
      delete[] p[0][0];
      delete[] p[0][1];
      delete[] p[1][0];
      delete[] p[1][1];
    }

    ~Coverage(){ 
      delete_subarrays( snp_cov );
      delete_subarrays( tot_cov );
      delete_subarrays( aln_centres );

      delete[] cov_hi_f;
      delete[] cov_hi_r;
      delete[] cov_lo_f;
      delete[] cov_lo_r;

      delete[] snp_cov_lo_f;
      delete[] snp_cov_lo_r;
      delete[] snp_cov_hi_f;
      delete[] snp_cov_hi_r;
    }

    void memset_subarrays(int *p[2][2]){
        memset( p[0][0], 0, range * 2 * sizeof(int));
        memset( p[0][1], 0, range * 2 * sizeof(int));
        memset( p[1][0], 0, range * 2 * sizeof(int));
        memset( p[1][1], 0, range * 2 * sizeof(int));
    }

    void clear_arrays(){ // memset all to zeros
        
        memset_subarrays( snp_cov );
        memset_subarrays( tot_cov );
        memset_subarrays( aln_centres );

        memset(cov_lo_f, 0, range * 2 * sizeof(int));
        memset(cov_lo_r, 0, range * 2 * sizeof(int));
        memset(cov_hi_f, 0, range * 2 * sizeof(int));
        memset(cov_hi_r, 0, range * 2 * sizeof(int));
       
        memset(snp_cov_lo_f, 0, range * 2 * sizeof(int));
        memset(snp_cov_lo_r, 0, range * 2 * sizeof(int));
        memset(snp_cov_hi_f, 0, range * 2 * sizeof(int));
        memset(snp_cov_hi_r, 0, range * 2 * sizeof(int));
    }

    float estimate( int & read_length );
    void compute_cov( Alignment * aln );
    void print_cov( const int & cc, FILE *file );
    bool within( const int & x );

};

#endif // Coverage_H_
