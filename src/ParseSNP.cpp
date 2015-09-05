/*
  ParseSNP.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#include "ParseSNP.h"
#include <sstream>
#include <stdexcept>
// #include <unistd.h>     //for the function usleep

#define NA 18446744073709551615

void ParseSNP::parse_cmd_line(int argc, char *argv[]) {

    CmdLine cmdss("computes local coverage and a derived BOD quality score", ' ', version.c_str());
    
    ValueArg<string> read_filename_arg("q", "mapped_reads",
            "Mapped read file (sam/bam)", true, "#", "string");
    cmdss.add(read_filename_arg);
    //
    ValueArg<string> ref_file_arg("r", "ref_file", "Reference file (fasta)", true, "#", "string");
    cmdss.add(ref_file_arg);
    //
    ValueArg<string> snp_file_arg("s", "snp_file", "VCF file", true, "#", "string");
    cmdss.add(snp_file_arg);
    //
    ValueArg<int> leng("l", "read_length", "maximum read length", true, -1, "int");
    cmdss.add(leng);
    //
    ValueArg<string> db_file_arg("d", "db_file", "Sqlite file to write the output", false, "#", "string");
    //
    ValueArg<string> plot_file_arg("p", "plot_file",
            "Name of the file for the data required for the Python script",
            false, "#", "string");
    //
    cmdss.xorAdd( db_file_arg, plot_file_arg ); // <- whether Sqlite DB or plain text
    //
    ValueArg<string> sample_label_arg("t", "sample_label",
            "sample label for the database",
            false, "#", "string");
    cmdss.add(sample_label_arg);
    //
    MultiArg<string> exclude_contigs_arg("x", "exclude", "exclude contigs from analysis", false, "string");
    cmdss.add(exclude_contigs_arg);
    //
    MultiSwitchArg verbose_arg("v", "verbose", "print snp-by-snp progress", 0);
    cmdss.add(verbose_arg);
    //
    ValueArg<int> num_test_arg("b", "break_after", "max number of position per chromosome [test mode]", false, 0, "int");
    cmdss.add(num_test_arg);
 
    try {
        cmdss.parse(argc, argv); //parse arguments
        read_filename = read_filename_arg.getValue();
/*        output = out.getValue();
        if (output[0] == '#') {
            output = read_filename;
            output += ".mot";
        }
*/
        reffile = ref_file_arg.getValue();
        snpfile = snp_file_arg.getValue();
        read_length = leng.getValue();
        range =  3 * read_length / 2 ;
     
        if ( plot_file_arg.isSet() )
            plot_file = plot_file_arg.getValue();
            //db_file = '#';
        else if ( db_file_arg.isSet() )
            // plot_file = '#';
            plot_file = db_file_arg.getValue();
        //

       if ( sample_label_arg.isSet() )
            sample_label = sample_label_arg.getValue();
        else 
            sample_label = "[" + db_file_arg.getValue() + "]";

        clog << "table base name: " << sample_label << endl;

        exclude_contigs = exclude_contigs_arg.getValue();

        verbose = verbose_arg.getValue();
        if (verbose){
            num_test = num_test_arg.getValue();
            }
        //
        db_flag = db_file_arg.isSet();

        init();
        parseVCF();

    } catch (ArgException &e) {
        StdOutput out;
        out.failure(cmdss, e);
    }
}

string convert(int number) {
    if (number == 0)
        return "0";
    string temp = "";
    string returnvalue = "";
    while (number > 0) {
        temp += number % 10 + 48;
        number /= 10;
    }
    for (int i = 0; i < (int) temp.length(); i++)
        returnvalue += temp[temp.length() - i - 1];
    return returnvalue;
}

string gen_key(int chr, int pos) {
    string chr_ref = convert(chr+1);
    chr_ref += ":";
    chr_ref += convert(pos);
    return chr_ref;
}

///////////////////////////////////////////////////////////////////////////////////
void ParseSNP::init() {
    // Process Fasta File  ================================================
    ifstream fastaFile;
    fastaFile.open(reffile.c_str(), ifstream::in);
    if (!fastaFile.good()) {
        cerr << "Fasta Parser: could not open file: `" << reffile.c_str() << "`" <<endl;
        exit(0);
    }

        clog << "reading the reference FASTA file : `" << reffile.c_str() << "`"<< endl;

    safeGetline(fastaFile ,  buffer, buffer_size);
    // fastaFile.getline(buffer, buffer_size);
    int chr_ref = 0;
        while (!fastaFile.eof()) {

            if (buffer[0] == '>') {
            string chr;
            chr.clear();
            for (size_t i = 1;
                    i < buffer_size && buffer[i] != '\0' && ( buffer[i] != '\r' )  && buffer[i] != '\n'
                            && buffer[i] != ' '; i++) {
                chr += buffer[i];
            }
            chromosome_vector[chr.c_str()] = chr_ref; // fill in a mapping : chromosome name [char] -> chr_ref [size_t]
            chr_ref++;
        }
        fastaFile.getline(buffer, buffer_size);
    }
    clog << "In reference, " << chromosome_vector.size() << " contigs have been detected." << endl;
    fastaFile.close();
//    std::sort( chromosome_vector.begin(), chromosome_vector.end());
    // read_register_table();
    
   for (std::vector<string>::iterator it = exclude_contigs.begin(); it != exclude_contigs.end(); ++it){
          if (chromosome_vector.erase( *it ) )    clog << endl << "Excluded contig: \t" << *it << endl;
    }
}

///////////////////////////////////////////////////////////////////////////////////
void ParseSNP::parseVCF() {
    size_t COMMIT_EACH = 200;
    bool IMMEDIATE = false;

    ifstream vcfFile;
    vcfFile.open(snpfile.c_str(), ifstream::in);
    if (!vcfFile.good()) {
        clog << "SNP Parser: could not open file: " << snpfile.c_str()
                << endl;
        exit(0);
    }
        
        clog << "reading the VCF file" << endl;
         
//    vcfFile.getline(buffer, buffer_size);

    vector<int> tmp;
    string broken_chromosome = "";
//    snps.resize(chromosome_vector.size(), tmp);
   
    Parser * mapped_file = 0;
    FastaParser * fasta = new FastaParser(reffile);
    string ref;
    
    FILE * plotFile = NULL;
    // sqlite3pp::transaction * xct;
    bool transaction_flag = false;
    if (db_flag){
        clog << "writing to database: " <<  plot_file.c_str() << endl;
        db = new SqliteDb( plot_file.c_str() , verbose );
        db->init_sql_table( sample_label );
        db->init_register_table();
        db->place_register_record_begin(  snpfile , read_filename);
        db->init_contig_table();
        if (!db) {  throw std::range_error("null pointer to the database");   
        } else { clog<< "table has been successfully initialized\t" << db << endl;} ;

    } else {
        remove(plot_file.c_str());
        plotFile = fopen(plot_file.c_str(), "a");
        if (plotFile == NULL) {
            cerr << "Error in printing: The file or path that you set "
//                << output.c_str()
                << " is not valid. It can be that there is no disc space available."
                << endl;
                ios_base::failure("cannot open output file for writing!");
                exit(0);
        }
   };

    if (read_filename.find(".bam") != string::npos) {
        mapped_file = new BamParser(read_filename);
    }
    clog << "======================================" << endl;
    clog << mapped_file->get_header() << endl;
    clog << "======================================" << endl;

//    clog << "num of chr " << genome.size() << endl;
 //   clog << "first chr size " << genome[0].size() << endl;
    clog << "`range` has been set to: " << range << endl;
    size_t chr_bam = 0;
    size_t temp_chr_ref = NA;
    size_t chr_ref = NA;
    int pos = 0;
    //    int n_snp = 0;  // <- global
    // ref = fasta->getChr(chr_ref);

    Coverage * cov;
    cov = new Coverage(range);

    while (!vcfFile.eof()) {

//        vcfFile.getline(buffer, buffer_size);
        safeGetline(vcfFile ,  buffer, buffer_size);
        if (buffer[0] == '#'){
            if (verbose){ clog << "skipping comment" << endl ;
                std::string s( buffer );
                if (s.find_last_of("\r") > 0){
                    clog << "this file contains Windows / MacOS specific end-of-line character \'\\r\'\n" <<  \
                        "@ position:\t" <<  s.find_last_of("\r") << "\n" << \
                        "consider running `dos2unix` or similar!" << endl;
                }
            }
            continue;
        }
        int field_count = 0;
        string current_chr;
        // read vcf file:
        // `i` -- runs for symbols
        // `field_count` -- runs for fields
        for (size_t i = 0; field_count< 2 && i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
            if (field_count == 0 && buffer[i] != '\t') {
                current_chr += buffer[i];
            }

            if (field_count == 1 && buffer[i - 1] == '\t') { //start: pos column
                pos = atoi(&buffer[i]) - 1;
                break;
            } // end: pos column
            if (buffer[i] == '\t') {
                field_count++;
            }
        }
        // process chromosome:
        if (verbose) { clog << "reading\t" << current_chr.c_str() << "\t" << pos << endl;}; 
        if ( chromosome_vector.count(current_chr.c_str()) > 0){ //found

            temp_chr_ref = chromosome_vector[current_chr.c_str()];
            if  (temp_chr_ref != chr_ref) { // new chromosome
                chr_ref = chromosome_vector[current_chr.c_str()];
                clog << endl; // "; getting next bam chromosome..." << endl;
                chr_bam = mapped_file->GetReferenceID( current_chr);
                if (chr_ref != chr_bam ){
                    cerr << "contig: " << current_chr << "; [fasta#:] " << chr_ref << "; [bam#:] " << chr_bam << ";  mismatch!" << endl;
                }
                if (verbose) {
                    cout << endl;
                    cout << "contig [fasta#:]\t " << chr_ref + 1 \
                    << "\t[bam#:]\t" << chr_bam + 1 \
                    << "\t[vcf:]\t" << current_chr.c_str() << "\t[fasta:]\t " << fasta->contig_name[chr_ref] << endl;
                } else {
                    clog << "contig # " << chr_ref+1 << endl;
                }
                ref = fasta->getChr(chr_ref);

                if (db_flag) {
                    if (transaction_flag){
                        clog << "commiting" << endl;
                        if (!db) {  throw std::range_error("null pointer to the database");    };
                        db->intermediate_commit();
                    } else {
                        if (!db) {  throw std::range_error("null pointer to the database");    };
                        db->new_transaction();
                        transaction_flag = true;
                    }
                    db->place_contig_table_record( chr_ref, current_chr);
                }
                n_snp = 0;
           }
             
        } else if (broken_chromosome.compare(current_chr)!=0){
            broken_chromosome = current_chr;
            cerr << endl << "Contig not found: \"\t" << broken_chromosome << "\t\"" << endl;
            continue;
        } else {
            if (verbose){ cerr << "\rskipping:\t" << current_chr << "\t" << pos << endl; }
            continue;
        }

        // process position
        if (verbose){
           if (num_test && (n_snp >= num_test )){
                continue;
           } else {
                n_snp ++;
           }
          //  clog << endl;
            // the info will be printed later in the `process_snp` routine
        } else {
            clog << "\r" << setfill(' ') << setw(8) << pos+1;
        }
   
        cov->reset(current_chr, chr_ref, pos); // cov->reset(chr_ref, pos);

        if (! process_snp(cov, ref, mapped_file, chr_ref, chr_bam) ){
            if (verbose>1) {
                throw std::logic_error(" error while processing a snp!!! ");
            } else {
                cerr << "  error while processing a snp! skipping... " << endl;
            }
        }
        if (db_flag){
            db->print_cov_db( *cov );
        } else { cov->print_cov(chr_ref, plotFile); }

        cov->estimate( read_length );

        if ( !(n_snp % COMMIT_EACH) && n_snp >0 ){
            clog << " | commiting" << endl;
            db->intermediate_commit();
       }
        // finally:       
        // vcfFile.getline(buffer, buffer_size);
    }
    clog << endl << "VCF file `" << snpfile.c_str() << "` has been successfully processed" << endl;
    vcfFile.close();

    if (db_flag){
        db->place_register_record( );
        db->composite_index();
         // xct->rollback(); // sqlite
        db->commit();
    } else  fclose(plotFile);
}
///////////////////////////////////////////////////////////////////////////////////

bool ParseSNP::process_snp(Coverage* cov, string & ref, Parser * mapped_file, const size_t &cc, const size_t & chr_bam){
    int leftPos = cov->start_pos;
    int rightPos =  cov->pos + range;
    int MIN_QUALITY = 20;
    int MIN_LENGTH = 20;
    
    std::ostringstream oss;
    if (verbose){       
        clog << "                                                                                                                            ";
        oss << "\r";
        oss << "chr # " << setfill(' ') << setw(2) << cc + 1 << " : " \
        << "(#" <<  setfill(' ') << setw(4) <<  n_snp  << ") " 
        << setfill(' ') << setw(8) << cov->pos + 1 << " >> " \
        << setfill(' ') << setw(8) << leftPos + 1 << " ... " \
        << setfill(' ') << setw(8) << rightPos + 1 ; } ;
    // set the region of interest
    if (!mapped_file->SetRegion( (int) chr_bam, leftPos, rightPos)){
        cerr << endl << "cannot jump to position " << cov->pos << " on chr " << cc << endl;
    }
    
    Alignment * tmp_aln = mapped_file->parseRead(MIN_QUALITY);
    vector<Alignment *> al_vect;
    size_t aln_count = 0;
    
    // these two loops CAN BE COMBINED !!!
    while(!tmp_aln->getSequence().first.empty()) {
                // collect alignments of the current SNP neighbourhood into `al_vect`                
                al_vect.push_back(tmp_aln);
                tmp_aln = mapped_file->parseRead(MIN_QUALITY);
                if (verbose >1) {
                    clog << oss.str() << " << aln. count:"<< setfill(' ') << setw(9) << aln_count ;
                    aln_count++;
                    }
    }

    if (verbose > 1){
         oss << " << aln. count:"<< setfill(' ') << setw(9) << aln_count ; // << "//" << al_vect.size();
   }

    std::ostringstream sc;
    for (size_t aa = 0; aa < al_vect.size(); aa++ )    {
    al_vect[aa]->processAlignment(ref); // a huge chunk has been factored out to the `processAlignment` method
        if ( (al_vect[aa]->getIdentity() >= 0.90) && (al_vect[aa]->getOrigLen() > MIN_LENGTH) ) {
             cov->compute_cov(al_vect[aa]); 
             if (verbose>1) { 
                 clog << "\r";
                 for (size_t jj = range; jj< range+50; jj=jj+3){
                    clog << setfill(' ') << setw(4) << (int) ( * cov->snp_cov)[0][0][jj];
//                    usleep(2000);
                 }
             }
             if (verbose == 1) { clog << oss.str() <<" <<  coverage: " << setfill(' ') << setw(9) << aa ;}
        }
    }
  if (verbose > 0){
    if (( al_vect.size() ) && (al_vect[0]->getRefID() != chr_bam) ){ 
       cerr << "| wrong aln ref id: " << al_vect[0]->getRefID() << ", chr_bam: " << chr_bam; 
    }
  }
    return true;
}

void ParseSNP::compute(){

}

//std::istream& ParseSNP::safeGetline(std::istream& is, std::string& t)
std::istream& ParseSNP::safeGetline(std::istream& is, char * t, const size_t & buf_len)
{
    //t.clear();

    // The characters in the stream are read one-by-one using a std::streambuf.
    // That is faster than reading them one-by-one using the std::istream.
    // Code that uses streambuf this way must be guarded by a sentry object.
    // The sentry object performs various tasks,
    // such as thread synchronization and updating the stream state.

    std::istream::sentry se(is, true);
    std::streambuf* sb = is.rdbuf();

    for(size_t i=0; i < buf_len ; i ++ ) {
        int c = sb->sbumpc();
        switch (c) {
        case '\n':
            t[i] = '\0';
            return is;
        case '\r':
            if(sb->sgetc() == '\n')
                sb->sbumpc();
            t[i] = '\0';
            return is;
        case EOF:
            // Also handle the case when the last line has no line ending
            t[i] = '\0';
            if(i==0) // t.empty())
                is.setstate(std::ios::eofbit);
            return is;
        default:
            t[i] = (char)c;
        }
    }
}
