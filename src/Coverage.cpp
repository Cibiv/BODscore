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
            cov_hi[ii] = (*tot_cov)[HI][RV][ii] + (*tot_cov)[HI][FW][ii];
            cov_90[ii] = cov_hi[ii] + (*tot_cov)[LO][RV][ii] + (*tot_cov)[LO][FW][ii];
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
            (*tot_cov)[hi_lo][fw_rv][jj]++;
            if (within){ 
                (*snp_cov)[snp_hi_lo][fw_rv][jj]++;
                }
            jj++;
            if (i > 2*range){ break; }
        }
    } else {
        for (size_t i = max(ZERO, -al_pos); i < end_point ; i++) {
             (*tot_cov)[hi_lo][fw_rv][jj]++;
             if ( within ){ 
                 (*snp_cov)[snp_hi_lo][fw_rv][jj]++;
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
        (*aln_centres)[snp_hi_lo][fw_rv][centre_al]++; 
    }
    // summary 
    altCounts = (*tot_cov)[LO][FW][range] + (*tot_cov)[LO][RV][range];
    totCounts = altCounts + (*tot_cov)[HI][FW][range] + (*tot_cov)[HI][RV][range];
    snp_ratio = (float) altCounts / (float) totCounts;

};

void Coverage::print_cov(const int & cc, FILE *file){
    // clog << "printing ..." << endl;

    fprintf(file, "%u\t", cc  + 1 );
    fprintf(file, "%u\t", pos + 1 );
    fprintf(file, "%f\t", score);
    
    tot_cov->print_subarrays( file );
    aln_centres->print_subarrays( file );
    snp_cov->print_subarrays( file );
    
    fprintf(file, "%c", '\n');
}

bool Coverage::within(const int & x){
    return ( (pos - x) < range && (x - pos) < range ) ? true : false ;
}
