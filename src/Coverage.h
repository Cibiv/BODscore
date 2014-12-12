#ifndef Coverage_H_
#define Coverage_H_

#include<string>
#include<Alignment.h>
#include <stdexcept>

using namespace std;

int const step=5;

class Coverage{
private:
    float conf_val(const int & start_cov, const int & stop_cov, const int & max_cov);

    void allocate_arrays(){
        cov_90 = new int[range * 2];
        cov_100 = new int[range * 2];
        cov_diff =  new int[range * 2];
        alnm_centres = new int[range * 2];
        alnm_centres_mut = new int[range * 2];
        loc_cov_90 = new int[range * 2];
        loc_cov_100 = new int[range * 2];
    
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
    int* cov_100;
	int* cov_90;
    int* cov_diff;
    int* alnm_centres;
    int* alnm_centres_mut;
    int* loc_cov_90;
    int* loc_cov_100;
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
        cov_90 = obj.cov_90;
        cov_100 = obj.cov_100;
        cov_diff = obj.cov_diff;
        clear_arrays(); // memset all to zeros
    }

     Coverage(const int &r):  range(r) {
        allocate_arrays();
     }

    ~Coverage(){ 
      delete[] cov_100;
      delete[] cov_90;
      delete[] cov_diff;
      delete[] alnm_centres;
      delete[] alnm_centres_mut;
      delete[] loc_cov_90;
      delete[] loc_cov_100;
    }

    void clear_arrays(){ // memset all to zeros
        memset(cov_90, 0, range * 2 * sizeof(int));
        memset(cov_100, 0, range * 2 * sizeof(int));
        memset(cov_diff, 0, range * 2 * sizeof(int));
        memset(alnm_centres, 0, range * 2 * sizeof(int));
        memset(alnm_centres_mut, 0, range * 2 * sizeof(int));
        memset(loc_cov_90, 0, range * 2 * sizeof(int));
        memset(loc_cov_100, 0, range * 2 * sizeof(int));
    }

    float estimate(int & read_length);
    void compute_cov(Alignment * aln);
    void print_cov(const int & cc, FILE *file);
    bool within(const int & x);

};

#endif // Coverage_H_
