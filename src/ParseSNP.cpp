/*
  ParseSNP.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#include "ParseSNP.h"
#include <sstream>
#include <stdexcept>

void ParseSNP::parse_cmd_line(int argc, char *argv[]) {

    CmdLine cmdss("Comp Motiv", ' ', version.c_str());
    
    ValueArg<string> read_file("q", "mapped_reads",
            "Mapped read file (sam/bam)", true, "#", "string");
    cmdss.add(read_file);
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
    cmdss.xorAdd( db_file_arg, plot_file_arg );
    //
   // ValueArg<string> out("o", "outputfile", "File to print the distripution to",
   //         true, "#", "string");
   // cmdss.add(out);
    //
    ValueArg<string> sample_label_arg("t", "sample_label",
            "sample label for the database",
            false, "#", "string");
    cmdss.add(sample_label_arg);
    //
    SwitchArg verbose_arg("v", "verbose", "print snp-by-snp progress", false);
    cmdss.add(verbose_arg);

    try {
        cmdss.parse(argc, argv); //parse arguments
        read_filename = read_file.getValue();
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

        verbose = verbose_arg.getValue();
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

void ParseSNP::init() {
    ifstream fastaFile;
    fastaFile.open(reffile.c_str(), ifstream::in);
    if (!fastaFile.good()) {
        cerr << "Fasta Parser: could not open file: " << reffile.c_str() << endl;
        exit(0);
    }

        clog << "reading the reference FASTA file : " << reffile.c_str() << endl;

    fastaFile.getline(buffer, buffer_size);
    int chr_ref = 0;
    while (!fastaFile.eof()) {
        if (buffer[0] == '>') {
            string chr;
            chr.clear();
            for (size_t i = 1;
                    i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'
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
}

void ParseSNP::parseVCF() {
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
    sqlite3pp::transaction * xct;
    bool transaction_flag = false;
    if (db_flag){
        clog << "writing to database: " <<  plot_file.c_str() << endl;
        db  = new sqlite3pp::database( plot_file.c_str() ) ;
        init_register_table();
        place_register_record_begin();

        init_tag_table();

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
    size_t chr_ref = 100000;
    int pos = 0;
    // ref = fasta->getChr(chr_ref);

    Coverage * cov;
    cov = new Coverage(range);

    while (!vcfFile.eof()) {

        vcfFile.getline(buffer, buffer_size);

        if (buffer[0] == '#'){
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
        // clog << current_chr.c_str() << "  "; 
//        int

        if ( chromosome_vector.count(current_chr.c_str()) > 0){ //found

            if  (chromosome_vector[current_chr.c_str()] != chr_ref) { // new chromosome
                chr_ref = chromosome_vector[current_chr.c_str()];
                clog << endl;
                if (chr_ref != mapped_file->GetReferenceID( current_chr) ){
                    cerr << "chromosome : " << current_chr << "[fasta #] " << chr_ref << "[bam #]" << mapped_file->GetReferenceID( current_chr)<< endl;
                    // throw std::logic_error("chromosome numbering in the fasta file and bam file do not match!");
                    cerr << "chromosome numbering in the fasta file and bam file do not match!" << endl;
                }
                if (verbose) {
                    cout << "chromosome # " << chr_ref+1 \
                    << "\t[bam:]\t" << mapped_file->GetReferenceID( current_chr ) + 1 \
                    << "\t" << current_chr.c_str() << " \t " << fasta->contig_name[chr_ref] << endl;
                } else {
                    clog << "chromosome # " << chr_ref+1 << endl;
                }
                ref = fasta->getChr(chr_ref);

                if (db_flag) {
                    if (transaction_flag){
                        clog << "commiting" << endl;
                        xct->commit();
                        delete [] xct;
                        //  xct->rollback(); // sqlite transaction;
                    }
                    xct = new sqlite3pp::transaction(*db);
                    init_sql_table(chr_ref, current_chr);
                }
            };
                      
        } else if (broken_chromosome.compare(current_chr)!=0){
            broken_chromosome = current_chr;
            cerr << endl << "Contig not found: \t" << broken_chromosome << endl;
            continue;
        } else {
            if (verbose){ cerr << "\r skipping: " << current_chr << " : " << pos; }
            continue;
        }

        // process position
        if (verbose){
            // clog << endl;
            // the info will be printed later in the `process_snp` routine
        } else {
            clog << "\r" << setfill(' ') << setw(8) << pos+1;
        }
     
        cov->pos = pos;
        cov->start_pos = pos > range ? pos - range : 0;
        cov->clear_arrays();
        process_snp(cov, ref, mapped_file, chr_ref);
        if (db_flag){
            cov->print_cov_db(table_name.c_str(), *db);
        } else { cov->print_cov(chr_ref, plotFile); }

        cov->estimate( read_length );

        // finally:       
        // vcfFile.getline(buffer, buffer_size);
    }
    clog << "VCF file has been successfully processed" << endl;
    place_register_record();

    xct->commit();
    // xct->rollback(); // sqlite
    vcfFile.close();
    if (db_flag){
        
    } else  fclose(plotFile);
}
///////////////////////////////////////////////////////////////////////////////////
void ParseSNP::read_register_table(){
// one needs to check whether the table `register` exists whatsoever here:

//    const char * qrycheck = "SELECT  count(*) FROM sqlite_master WHERE type='table' AND name='register';";
//    sqlite3pp::query qry(*db, qrycheck);
    try {
        sqlite3pp::query qry(*db, "SELECT name, notes  FROM register");

        for (int i = 0; i < qry.column_count(); ++i) {
          cout << qry.column_name(i) << "\t";
        }
    } catch (int n) {
        clog << "`register` table was not found" << endl;
    }
}
void ParseSNP::parseSQLite() {
    return; }
///////////////////////////////////////////////////////////////////////////////////
void ParseSNP::process_snp(Coverage* cov, string & ref, Parser * mapped_file, const size_t &cc){
    int leftPos = cov->start_pos;
    int rightPos =  cov->pos + range;
    int MIN_QUALITY = 0;
    
    std::ostringstream oss;
    if (verbose){ 
        oss << "chr # " << setfill(' ') << setw(2) << cc + 1 << " : " \
        << setfill(' ') << setw(8) << cov->pos + 1 << " >> " \
        << setfill(' ') << setw(8) << leftPos + 1<< " ... " \
        << setfill(' ') << setw(8) << rightPos + 1 << " << coverage: " ; } ;
    // set the region of interest
    if (!mapped_file->SetRegion( (int) cc, leftPos, rightPos)){
    cerr << "cannot jump to position " << cov->pos << " on chr " << cc << endl;
    }
    
    Alignment * tmp_aln = mapped_file->parseRead(MIN_QUALITY);
    vector<Alignment *> al_vect;
    size_t aln_count = 0;
    
    // these two loops CAN BE COMBINED !!!
    while(!tmp_aln->getSequence().first.empty()) {
                // collect alignments of the current SNP neighbourhood into `al_vect`                
                al_vect.push_back(tmp_aln);
                aln_count++;
                if (verbose) { clog << "\r" << oss.str() << setfill(' ') << setw(9) << aln_count ;}
                tmp_aln = mapped_file->parseRead(MIN_QUALITY);
    }
    for (size_t aa = 0; aa < al_vect.size(); aa++ )    {
        al_vect[aa]->processAlignment(ref); // a huge chunk has been factored out to the `processAlignment` method
        if (al_vect[aa]->getIdentity() >= 0.90) {
             cov->compute_cov(al_vect[aa]); 
        // } else { clog << " low identity! " ;
        }
    }
    // al_vect.clear();
}
    

void ParseSNP::compute(){

}

void ParseSNP::init_sql_table( size_t & chr_id, string & chr_name){
    table_name = sample_label + "__coverage_" + std::to_string( (int) chr_id + 1);

    char sql[256];
    sprintf(sql, "CREATE TABLE %s ("  \
                       "pos INT PRIMARY KEY     NOT NULL , " \
                       "totCounts INT, " \
                       "refCounts INT, " \
                       "snp_ratio FLOAT, " \
                       "score FLOAT , " \
                       "totCov BLOB, " \
                       "snpCov BLOB, " \
                       "alnCtr BLOB);", table_name.c_str() );
                      // "tot_cov CHAR(%u) );", (int) chr_ref+1, (range*2+1)*4 );
                      //    ", cov_hi    CHAR(50)" ", cov_lo    CHAR(50)"
                      //                            cout << sql[120] << endl;
    exec_sql_log(sql);
    place_tag_table_record(sample_label, table_name, chr_id, chr_name );
}

void ParseSNP::init_sql_table( size_t & chr_id){
    table_name = sample_label + "__coverage_" + std::to_string( (int) chr_id + 1);

    char sql[256];
    sprintf(sql, "CREATE TABLE %s ("  \
                       "pos INT PRIMARY KEY     NOT NULL , " \
                       "totCounts INT, " \
                       "refCounts INT, " \
                       "snp_ratio FLOAT, " \
                       "score FLOAT , " \
                       "totCov BLOB, " \
                       "snpCov BLOB, " \
                       "alnCtr BLOB);", table_name.c_str() );
                      // "tot_cov CHAR(%u) );", (int) chr_ref+1, (range*2+1)*4 );
                      //    ", cov_hi    CHAR(50)" ", cov_lo    CHAR(50)"
                      //                            cout << sql[120] << endl;
      exec_sql_log(sql);
}

void ParseSNP::init_register_table( ){
    const char * sql = "CREATE TABLE IF NOT EXISTS register ( "
                 "name TEXT PRIMARY KEY, "
            "vcf_file TEXT, "
            "mapper TEXT, "
            "nm INT, "
            "mq INT, "
            "nt_q INT, "
            "wt BOOL, "
            "notes TEXT ); ";
   exec_sql_log(sql);
}

void ParseSNP::place_register_record_begin(){
    std::string note =  "coverage analysis started : " + snpfile;
    sqlite3pp::command cmd( *db, "INSERT OR REPLACE INTO register (name, notes) VALUES (:name, :notes) ");
    cmd.bind(":name", sample_label.c_str());
    cmd.bind(":notes", note.c_str() );

// sqlite3pp::command cmd( *db, "INSERT OR REPLACE INTO register (name, vcf_file, notes) VALUES (:name, :vcf_file, :notes) ");
//    cmd.bind(":vcf_file", snpfile.c_str() );
    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }

}
void ParseSNP::place_register_record(){
    std::string note =  "coverage_done : " + snpfile;
    sqlite3pp::command cmd( *db, "INSERT OR REPLACE INTO register (name, notes) VALUES (:name, :notes) ");
    cmd.bind(":name", sample_label.c_str());
    cmd.bind(":notes", note.c_str() );

// sqlite3pp::command cmd( *db, "INSERT OR REPLACE INTO register (name, vcf_file, notes) VALUES (:name, :vcf_file, :notes) ");
//    cmd.bind(":vcf_file", snpfile.c_str() );
    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }

}

void ParseSNP::init_tag_table(){

    char sql[256];
    sprintf(sql, "CREATE TABLE %s ("  \
                       "chr_table TEXT PRIMARY KEY     NOT NULL , " \
                       "chr_name TEXT, " \
                       "data_type TEXT, " \
                       "chr_id INT, " \
                       "max_pos INT, " \
                       "records INT);", sample_label.c_str() );
    exec_sql_log(sql);
}

void ParseSNP::exec_sql_log(char const * sql){
    try {
        int exitcode = db->execute(sql);
        if (verbose){   clog << "[sqlite3:] " << sql << " [returned:]  " << exitcode << endl;}
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
}

void ParseSNP::place_tag_table_record( string & sample_label, string & table_name, size_t & chr_ref, string & chr_name ){

    char sql[256];
    sprintf(sql, "INSERT OR REPLACE INTO %s " \
        "( chr_table, chr_name, data_type, chr_id ) " \
        "VALUES " \
        "( :chr_table, :chr_name, :data_type, :chr_id);", sample_label.c_str() );
 
    sqlite3pp::command cmd( *db, sql);
    cmd.bind(":chr_table", table_name.c_str());
    cmd.bind(":data_type", "coverage");
    cmd.bind(":chr_name", chr_name.c_str() );
    cmd.bind(":chr_id", (int) chr_ref + 1 );

    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }
}

