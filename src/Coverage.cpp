/*
by Dmytro S Lituiev, 8 Dec 2014;
based on Fritz Sedlazeck's;
*/

#include "Coverage.h"
const int Coverage::ZERO = 0;

float Coverage::conf_val(const int & start_cov, const int & stop_cov, const int & max_cov) {
        //always: y2-y1=1
        //    cout<<"\ts1: "<<(float)(stop_cov-start_cov)<<endl;
        //    cout<<"\tsmax: "<<(float)(max_cov-start_cov)<<endl;
        if (max_cov - start_cov == 0 || stop_cov == 0 || (stop_cov - start_cov) == 0) {
            return 2.0; //-1.95;
        }

        return (((float) (stop_cov - start_cov) / (float) (max_cov - start_cov)) / step );
}

float Coverage::estimate(int & read_length) {

        int cov_hi[2*range];
        int cov_90[2*range];
        for (size_t ii=0; ii< 2*range; ii++){
            cov_hi[ii] = tot_cov[HI][RV][ii] + tot_cov[HI][FW][ii];
            cov_90[ii] = cov_hi[ii] + tot_cov[LO][RV][ii] + tot_cov[LO][FW][ii];
            }

        int max_right = 0;
        int max_left = 0;
        int max_cov_90 = 0;
        // cov_xx[range] -- coverage in the very SNP position
        for (int i = 0; i < read_length; i++) {
            if (i < range && cov_hi[i + range] > max_right) {
                max_right = cov_hi[i + range];
            }
            if (i - range >= 0 && cov_hi[range - i] > max_left) {
                max_left = cov_hi[range - i];
            }
            if (i < range  && (cov_90[i + range] > max_cov_90) ) {
                max_cov_90 = (cov_90[i + range] );
            }
            if (i - range >= 0 && (cov_90[range - i]  > max_cov_90)) {
                max_cov_90 = (cov_90[range - i] );
            }
        }

        //check the slope of the left and right angle:
        float r1 = conf_val(cov_hi[range], cov_hi[range + step], max_right);
        float r2 = conf_val(cov_hi[range], cov_hi[range - step], max_left);
        float r_tot; //=1-(0.05-(r1+r2)/2);

        if (r1 < 0 || r2 < 0) {
            r_tot = 0;
        } else {
            r_tot = 1 - ((r1 + r2) / 2);
        }

        //check the cov at position of cov_90:
        float r_cov = ((float) (cov_90[range]) / (float) max_cov_90);
        //    cout<<"\tr: "<<(r_cov+r_tot)/2<<endl;
        if ((r_cov + r_tot) / 2 < 0) {
            score = 0;
        } else {
            score =  (r_cov + r_tot) / 2;
        }
        if (score > 1){ 
            char buffer [50];
            sprintf(buffer, "score exceeds one: %4.2f", score );
            throw range_error( buffer ); 
        }
        return score;
    }

void Coverage::compute_cov(Alignment * aln) {

//    clog<<aln->getSequence().first<<endl;
//    clog<<"P: "<<aln->getIdentity()<<" "<<aln->getRefID()<<endl;
//    clog<<aln->getSequence().second<<endl;
    int const al_pos = aln->getPosition() - start_pos;
    //   =========================================================================
    int pos_snp_in_aln = pos - aln->getPosition();

    int fw_rv = (int) (aln->IsReverseStrand());
    int hi_lo = (int) (aln->getIdentity() == 1);
    int snp_hi_lo = (int) ( (aln->getIdentity() == 1) && aln->getIdentity( pos_snp_in_aln ) );

    bool within = pos_snp_in_aln > 0 && pos_snp_in_aln <  aln->getSequence().first.size(); 
   
    size_t end_point = min(aln->getSequence().first.size(), 2*range - al_pos);
    size_t offset = 0;
    int jj = al_pos + max(ZERO, -al_pos);

    if (aln->getIdentity() == 1) {
        for (size_t i = max(ZERO, -al_pos); i < end_point ; i++) {
            tot_cov[hi_lo][fw_rv][jj]++;
            if (within){ 
                snp_cov[snp_hi_lo][fw_rv][jj]++;
                }
            jj++;
            if (i > 2*range){ break; }
        }
    } else {
        for (size_t i = max(ZERO, -al_pos); i < end_point ; i++) {
             tot_cov[hi_lo][fw_rv][jj]++;
             if ( within ){ 
                 snp_cov[snp_hi_lo][fw_rv][jj]++;
                 }
             if (aln->getSequence().second[i] != '-') {
                 jj++;
                 if (jj > 2*range){ break; }
             }
        }
    }

//    int pos_snp_in_aln = pos - aln->getPosition(); 
    // ==== IF snp is within the alignment
    if ( within ){
        int centre_al =  range - pos_snp_in_aln + aln->getSequence().first.size() / 2  ;
        aln_centres[snp_hi_lo][fw_rv][centre_al]++; 
    }

};

void Coverage::print_cov(const int & cc, FILE *file){
    // clog << "printing ..." << endl;

    fprintf(file, "%u\t", cc  + 1 );
    fprintf(file, "%u\t", pos + 1 );
    fprintf(file, "%f\t", score);
    
    print_subarrays( file, tot_cov );
    
    print_subarrays( file, aln_centres);
    
    print_subarrays( file, snp_cov );
    
    fprintf(file, "%c", '\n');
}
/*
char * uint16_to_char(unsigned short int * const val){
char bytes[2];
    bytes[0] = (val) & 0xFF;  // low byte
    bytes[1] = (val >> 8) & 0xFF;  // high byte
    return bytes;
}

void uint16_to_char(char * out_buffer, unsigned short int * const val){
    out_buffer[0] = (val) & 0xFF;  // low byte
    out_buffer[1] = (val >> 8) & 0xFF;  // high byte
    return;
}

void array_uint16_to_char(char * out_buffer, unsigned short int * arr, size_t len){ 
    for (size_t ii = 0; ii<len; ii++){
        uint16_to_char( out_buffer + 2*ii, arr[ii]);
    }
    return;
}
*/
void Coverage::print_cov_db(const char * table_name, sqlite3pp::database & db){
                        // clog << "printing ..." << endl;
//                    cout << "REACHED 1!" << endl;
                
    buf_len = range*2 + 1 ; // class int 
    int buf_tot_len = 4 * buf_len;

    char aln_ctr_str[buf_tot_len] ;
    sprint_subarrays(aln_ctr_str, aln_centres);

    char tot_cov_str[buf_tot_len] ;
    sprint_subarrays(tot_cov_str, tot_cov);

    char snp_cov_str[buf_tot_len] ;
    sprint_subarrays(snp_cov_str, snp_cov);

    int altCounts = tot_cov[LO][FW][range] + tot_cov[LO][RV][range];
    int totCounts = altCounts + tot_cov[HI][FW][range] + tot_cov[HI][RV][range];
    float snp_ratio = (float) altCounts / (float) totCounts;

/*  clog << "blob length: "<< tot_cov_str.length() << endl;
    if (tot_cov_str.length() < buf_tot_len){
        cerr << "short string!!" << endl << tot_cov_str.c_str() << endl;
    }
*/
    for (size_t i = 1; i<=4; i++){
        if (tot_cov_str[i*buf_len - 1] != 0){
            cerr << "corrupted byte ending # " << i ;
            cerr << " value " << (int) tot_cov_str[i*buf_len] << endl;
        }
    }


    char sql[1024];
    sprintf (sql, "INSERT OR REPLACE INTO %s "  \
         "(pos, totCounts, refCounts, snp_ratio, score, totCov, snpCov, alnCtr) " \
         "VALUES (:pos, :totCounts, :refCounts, :snp_ratio, :score, :tot_cov, :snp_cov, :aln_ctr) ", \
          table_name );

    sqlite3pp::command cmd( db, sql);
    cmd.bind(":pos", pos + 1 );
    cmd.bind(":totCounts", totCounts);
    cmd.bind(":refCounts", totCounts - altCounts);
    cmd.bind(":snp_ratio", snp_ratio);
    cmd.bind(":score", score);
    cmd.bind(":tot_cov", (void const*) tot_cov_str, buf_tot_len );
    cmd.bind(":snp_cov", (void const*) snp_cov_str, buf_tot_len);
    cmd.bind(":aln_ctr", (void const*) aln_ctr_str, buf_tot_len);
    
    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }
//                    cout << "REACHED 2!" << endl;
//      print_subarrays( file, tot_cov );
    
//      print_subarrays( file, aln_centres);
    
//      print_subarrays( file, snp_cov );
    
//      fprintf(file, "%c", '\n');
}
/*
std::string Coverage::col_names(const char * base){
    char buf [(strlen(base) + 8)*4 ] ;
    int n = 0;
    n = n+ sprintf(buf+n,  "%s_hi_fw, ", base);
    n = n+ sprintf(buf+n,  "%s_hi_rv, ", base);
    n = n+ sprintf(buf+n,  "%s_lo_fw, ", base);
    n = n+ sprintf(buf+n,  "%s_lo_rv ", base);
    std::string out = buf;
}
*/
void Coverage::sprint_subarrays( char * buf, int * const  p[2][2] ){

        pbdb = &Coverage::sprint_char_block;
//        pb = &Coverage::print_block;
        (this->*Coverage::pbdb)(buf            , p[HI][FW]);
        (this->*Coverage::pbdb)(buf + buf_len  , p[HI][RV]);
        (this->*Coverage::pbdb)(buf + buf_len*2, p[LO][FW]);
        (this->*Coverage::pbdb)(buf + buf_len*3, p[LO][RV]);
        
        // std::string my_string( buf, 4 * buf_len );
        // return my_string;

     }

void Coverage::print_subarrays(FILE *file, int * const  p[2][2] ){
        
        pb = &Coverage::print_char_block;
//        pb = &Coverage::print_block;
        (this->*Coverage::pb)(file, p[HI][FW]);
        (this->*Coverage::pb)(file, p[HI][RV]);
        (this->*Coverage::pb)(file, p[LO][FW]);
        (this->*Coverage::pb)(file, p[LO][RV]);
     }

void Coverage::sprint_char_block(  char outstr[], const int * var ){
    for (int j = 0; j < 2 * range; j++) {
        outstr[j] = var[j] + 1; // sprintf(outstr + 1 + j, "%C", var[j]+33);
    }
    outstr[ 2*range ] = 0;
    return;
}

void Coverage::print_char_block(FILE *file, const int * var ){
//    fprintf(file, "\t");
    char buf[range*2 + 2 ];
    sprint_char_block( buf, var);
//    printf ( "%.*s", range*2,  buf );
    fprintf(file, "%.*s", range*2,  buf );
}

void Coverage::print_block(FILE *file, const int * var ){
    fprintf(file, "|\t");
    for (size_t j = 0; j < (size_t) range * 2; j++) {
        fprintf(file, "%i\t", var[j]);
    }
}
bool Coverage::within(const int & x){
    return ( (pos - x) < range && (x - pos) < range ) ? true : false ;
}
