#ifndef QdrArray_H_
#define QdrArray_H_

// #include<string>
// #include <stdexcept>

using namespace std;

size_t const LO = 0;
size_t const HI = 1;
size_t const FW = 0;
size_t const RV = 1;

class QdrArray{
private:
    size_t subarr_len;
    size_t buf_len;
   
    int * qarr[2][2];

public:
    QdrArray( size_t subarr_len_ )
    : subarr_len(subarr_len_) , buf_len(subarr_len_ + 1)
    {
        qarr[0][0] = new int[ subarr_len ];
        qarr[0][1] = new int[ subarr_len ];
        qarr[1][0] = new int[ subarr_len ];
        qarr[1][1] = new int[ subarr_len ];
        
        memset_subarrays( );
        //  clog << "qarr[0][0][1] : " << (*this)[0][0][1] << endl;
    }

    ~QdrArray(){
      delete[] qarr[0][0];
      delete[] qarr[0][1];
      delete[] qarr[1][0];
      delete[] qarr[1][1];
    }

    void memset_subarrays(){
        memset( qarr[0][0], 0, subarr_len * sizeof(int));
        memset( qarr[0][1], 0, subarr_len * sizeof(int));
        memset( qarr[1][0], 0, subarr_len * sizeof(int));
        memset( qarr[1][1], 0, subarr_len * sizeof(int));
    }

    int ** operator[] (const size_t nIndex)
    {
        // if (nIndex <=1)
        return qarr[nIndex];
    }

   template<typename TT>
   const void sprint_subarrays( TT * buf ) const {
        sprint_char_block(buf            , qarr[HI][FW] );
        sprint_char_block(buf + buf_len  , qarr[HI][RV] );
        sprint_char_block(buf + buf_len*2, qarr[LO][FW] );
        sprint_char_block(buf + buf_len*3, qarr[LO][RV] );
    }

    template<typename TT>
    const void sprint_char_block(  TT outstr[],  const int * arr ) const {
        size_t buf_len = subarr_len + 1;
        for (int j = 0; j < buf_len - 1 ; j++) {
            outstr[j] = (TT) arr[j] + 1; // sprintf(outstr + 1 + j, "%C", var[j]+33);
        }
        // the last element is set to zero as a control point
        outstr[ buf_len - 1 ] = 0;
        return;
    }

    void print_subarrays(FILE *file ){
       print_block(file, qarr[HI][FW]);
       print_block(file, qarr[HI][RV]);
       print_block(file, qarr[LO][FW]);
       print_block(file, qarr[LO][RV]);
    }

    void print_block(FILE *file, const int * var ){
        fprintf(file, "|\t");
        for (size_t j = 0; j < (size_t) subarr_len; j++) {
            fprintf(file, "%i\t", var[j]);
        }
    }


};

// template<> const void sprint_subarrays( unsigned short int * buf ) const;

#endif // QdrArray_H_
