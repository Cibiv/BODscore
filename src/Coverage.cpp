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

        int max_right = 0;
        int max_left = 0;
        int max_cov_90 = 0;
        // cov_xx[range] -- coverage in the very SNP position
        for (int i = 0; i < read_length; i++) {
            if (i < range && cov_100[i + range] > max_right) {
                max_right = cov_100[i + range];
            }
            if (i - range >= 0 && cov_100[range - i] > max_left) {
                max_left = cov_100[range - i];
            }
            if (i < range  && cov_90[i + range] > max_cov_90) {
                max_cov_90 = cov_90[i + range];
            }
            if (i - range >= 0 && cov_90[range - i] > max_cov_90) {
                max_cov_90 = cov_90[range - i];
            }
        }

        //check the slope of the left and right angle:
        float r1 = conf_val(cov_100[range], cov_100[range + step], max_right);
        float r2 = conf_val(cov_100[range], cov_100[range - step], max_left);
        float r_tot; //=1-(0.05-(r1+r2)/2);

        if (r1 < 0 || r2 < 0) {
            r_tot = 0;
        } else {
            r_tot = 1 - ((r1 + r2) / 2);
        }

        //check the cov at position of cov_90:
        float r_cov = ((float) cov_90[range] / (float) max_cov_90);

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
//    clog << ":";
 
//    if (al_pos >= 0 && al_pos < 2 * range) {

//    size_t end_pos = (size_t) start_pos + range * 2; // MAKE GLOBAL CONSTANT
        
/*       int jj = al_pos + max(ZERO, -al_pos);
       for (size_t i = max(ZERO, -al_pos); i < end_point ; i++) {
            if (aln->getSequence().second[i] != '-') {
                if (aln->getIdentity() == 1) {
                    cov_100[jj]++;
                    } else {
                    cov_diff[jj]++; 
                    if ( cov_diff + al_pos == alnm_centres ) { cerr << " collision! al_pos: " << al_pos << endl;}
                    }
                cov_90[ j]++;
                jj++;
                if (jj > 2*range){ break; }
             } else {
                 cov_90[jj]++;
             }
        }
*/
//   =========================================================================
    int pos_snp_in_aln = pos - aln->getPosition();
    bool within = pos_snp_in_aln > 0 && pos_snp_in_aln <  aln->getSequence().first.size(); 

    size_t end_point = min(aln->getSequence().first.size(), 2*range - al_pos);
    size_t offset = 0;
    int jj = al_pos + max(ZERO, -al_pos);

    if (aln->getIdentity() == 1) {
        for (size_t i = max(ZERO, -al_pos); i < end_point ; i++) {
            cov_100[jj]++;
            cov_90[jj]++;
            if (within){ loc_cov_100[jj]++; loc_cov_90[jj]++; }
            jj++;
            if (i > 2*range){ break; }
        }
    } else {
        for (size_t i = max(ZERO, -al_pos); i < end_point ; i++) {
          if (aln->getSequence().second[i] != '-') {
             cov_diff[jj]++; 
             cov_90[jj]++;
             if (within){ loc_cov_90[jj]++; }
             jj++;
             if (jj > 2*range){ break; }
          } else {
              cov_90[jj]++; 
              cov_diff[jj]++;
              if (within){ loc_cov_90[jj]++; }
          }
          // if ( cov_diff + al_pos == alnm_centres ) { cerr << " collision! al_pos: " << al_pos << endl;}
        }
    }

//    int pos_snp_in_aln = pos - aln->getPosition(); 
    // ==== IF snp is within the alignment
    if ( within ){
        int centre_al =  range - pos_snp_in_aln + aln->getSequence().first.size() / 2  ;
        alnm_centres[centre_al]++;
        if (aln->getIdentity() < 1) { alnm_centres_mut[centre_al]++;}
    }

};

void Coverage::print_cov(const int & cc, FILE *file){

                        fprintf(file, "%u\t", cc  + 1 );
                        fprintf(file, "%u\t", pos + 1 );
                        fprintf(file, "%f\t", score);

                        print_block(file, cov_100);
                        print_block(file, cov_90);
                        
                        print_block(file, alnm_centres);
                        print_block(file, alnm_centres_mut);

                        print_block(file, loc_cov_100);
                        print_block(file, loc_cov_90);

                        fprintf(file, "%c", '\n');
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
