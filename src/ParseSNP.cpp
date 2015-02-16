/*
  ParseSNP.cpp
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#include "ParseSNP.h"

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
    string id = convert(chr+1);
    id += ":";
    id += convert(pos);
    return id;
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
    int id = 0;
    while (!fastaFile.eof()) {
        if (buffer[0] == '>') {
            string chr;
            chr.clear();
            for (size_t i = 1;
                    i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'
                            && buffer[i] != ' '; i++) {
                chr += buffer[i];
            }
            chrs[chr.c_str()] = id;
            id++;
        }
        fastaFile.getline(buffer, buffer_size);
    }
    clog << "In reference, " << chrs.size() << " contigs have been detected." << endl;
    fastaFile.close();
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
//    snps.resize(chrs.size(), tmp);
   
    vector<Coverage> tmp_cov;
    Parser * mapped_file = 0;
    FastaParser * fasta = new FastaParser(reffile);
    string ref;
    
    FILE * plotFile = NULL;
    sqlite3pp::transaction * xct;
    bool transaction_flag = false;
    if (db_flag){
        clog << "writing to database: " <<  plot_file.c_str() << endl;
        db  = new sqlite3pp::database( plot_file.c_str() ) ;
    } else {
        remove(plot_file.c_str());
        plotFile = fopen(plot_file.c_str(), "a");
        if (plotFile == NULL) {
            cout << "Error in printing: The file or path that you set "
//                << output.c_str()
                << " is not valid. It can be that there is no disc space available."
                << endl;
                ios_base::failure("cannot open output file for writing!");
                exit(0);
        }
   };

    if (read_filename.find(".bam") != string::npos) {
    //      cout << "BAM File" << endl;
    mapped_file = new BamParser(read_filename);
    }
                            
//    clog << "num of chr " << genome.size() << endl;
 //   clog << "first chr size " << genome[0].size() << endl;
    clog << "`range` has been set to: " << range << endl;
    size_t id = 100000;
    int pos = 0;
    // ref = fasta->getChr(id);

    Coverage * cov;
    cov = new Coverage(range);

    while (!vcfFile.eof()) {

        vcfFile.getline(buffer, buffer_size);

        if (buffer[0] == '#'){
            continue;
        }
        int count = 0;
        string current_chr;
        // read vcf file:
        // `i` -- runs for symbols
        // `count` -- runs for fields
        for (size_t i = 0; count< 2 && i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
            if (count == 0 && buffer[i] != '\t') {
                current_chr += buffer[i];
            }

            if (count == 1 && buffer[i - 1] == '\t') { //start: pos column
                pos = atoi(&buffer[i]) - 1;
                break;
            } // end: pos column
            if (buffer[i] == '\t') {
                count++;
            }
        }
        // process chromosome:
        // clog << current_chr.c_str() << "  "; 

        if ( chrs.count(current_chr.c_str()) > 0){ //found

            if  (chrs[current_chr.c_str()]!= id) { // new chromosome
                id = chrs[current_chr.c_str()];
                clog << endl;
                clog << "chromosome # " << id+1 << endl;
                ref = fasta->getChr(id);

                if (db_flag) {
                    if (transaction_flag){
                        cout << "commiting" << endl;
                        xct->commit();
                        delete [] xct;
                        //  xct->rollback(); // sqlite transaction;
                    }
                    xct = new sqlite3pp::transaction(*db);
                    init_sql_table(id);
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
        clog << "\r" << pos+1;
     
        cov->pos = pos;
        cov->start_pos = pos - range;
        cov->clear_arrays();
        process_snp(cov, ref, mapped_file, id);
        if (db_flag){
            cov->print_cov_db(table_name.c_str(), *db);
        } else { cov->print_cov(id, plotFile); }

        cov->estimate( read_length );

        // finally:       
        // vcfFile.getline(buffer, buffer_size);
    }
    clog << "VCF file has been successfully processed" << endl;
    init_register_table();
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
    ifstream vcfFile;
    vcfFile.open(snpfile.c_str(), ifstream::in);
    if (!vcfFile.good()) {
        clog << "SNP Parser: could not open file: " << snpfile.c_str()
                << endl;
        exit(0);
    }
        
        clog << "reading the VCF file" << endl;
         
    vcfFile.getline(buffer, buffer_size);

    vector<int> tmp;
    string broken_chromosome = "";
//    snps.resize(chrs.size(), tmp);
   
    vector<Coverage> tmp_cov;
    Parser * mapped_file = 0;
    FastaParser * fasta = new FastaParser(reffile);
    string ref;
    
    FILE * plotFile = NULL;
    sqlite3pp::transaction * xct;
    bool transaction_flag = false;
    if (db_flag){
        cout << "writing to database: " <<  plot_file.c_str() << endl;
        db  = new sqlite3pp::database( plot_file.c_str() ) ;
    } else {
        remove(plot_file.c_str());
        plotFile = fopen(plot_file.c_str(), "a");
        if (plotFile == NULL) {
            cout << "Error in printing: The file or path that you set "
//                << output.c_str()
                << " is not valid. It can be that there is no disc space available."
                << endl;
                ios_base::failure("cannot open output file for writing!");
                exit(0);
        }
   };

    if (read_filename.find(".bam") != string::npos) {
    //      cout << "BAM File" << endl;
    mapped_file = new BamParser(read_filename);
    }
                            
//    clog << "num of chr " << genome.size() << endl;
 //   clog << "first chr size " << genome[0].size() << endl;
    clog << "`range` has been set to: " << range << endl;
  //  int old_pos  = -100000;
    size_t id = 100000;
    // ref = fasta->getChr(id);

    Coverage * cov;
    cov = new Coverage(range);

    while (!vcfFile.eof()) {
        if (buffer[0] == '#'){
            vcfFile.getline(buffer, buffer_size);
            continue;
        }
        int count = 0;
        string current_chr;
        int pos = 0;
        for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
            if (count == 0 && buffer[i] != '\t') {
                current_chr += buffer[i];
            }

            if (count == 1 && buffer[i - 1] == '\t') { //start: pos column
                pos = atoi(&buffer[i]) - 1;
                break;
            } // end: pos column
            if (buffer[i] == '\t') {
                count++;
            }
            // remove db unique_ptr
        }
        // process chromosome:
        // clog << current_chr.c_str() << "  "; 

        if ( chrs.count(current_chr.c_str()) > 0){ //found

            if  (chrs[current_chr.c_str()]!= id) { // new chromosome
                id = chrs[current_chr.c_str()];
                clog << "current_chr " << id  << endl;
                ref = fasta->getChr(id);
                // old_pos = -100000;

                if (db_flag) {
                    if (transaction_flag){
                        cout << "commiting" << endl;
                        xct->commit();
                        delete [] xct;
                        //  xct->rollback(); // sqlite transaction;
                    }
                    xct = new sqlite3pp::transaction(*db);
                    init_sql_table(id);
                }
            };
                      
        } else if (broken_chromosome.compare(current_chr)!=0){
            broken_chromosome = current_chr;
            cerr << endl << "Contig not found: \t" << broken_chromosome << endl;
            continue;
        } else {
            // cerr << "skipping" << endl;
            continue;
        }

        // process position
        // old_pos = pos;

        clog << "chr " <<  id + 1 << " pos " << pos + 1;
     
        cov->pos = pos;
        cov->start_pos = pos - range;
        cov->clear_arrays();
        process_snp(cov, ref, mapped_file, id);
        if (db_flag){
            cov->print_cov_db(table_name.c_str(), *db);
        } else { cov->print_cov(id, plotFile); }

        cov->estimate( read_length );

        // finally:       

        vcfFile.getline(buffer, buffer_size);
    }
    clog << "VCF file has been successfully processed" << endl;
    init_register_table();
    place_register_record();

    xct->commit();
    // xct->rollback(); // sqlite
    vcfFile.close();
    if (db_flag){
        
    } else  fclose(plotFile);
}
///////////////////////////////////////////////////////////////////////////////////
void ParseSNP::process_snp(Coverage* cov, string & ref, Parser * mapped_file, const size_t &cc){
    int leftPos = cov->start_pos;
    int rightPos =  cov->pos + range;
    int MIN_QUALITY = 0;
    // clog << ":" << cc << ":" << leftPos << "-" << rightPos ;
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

void ParseSNP::init_sql_table( size_t & id){
    table_name = sample_label + "__coverage_" + std::to_string( (int) id+1);

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
                      // "tot_cov CHAR(%u) );", (int) id+1, (range*2+1)*4 );
                      //    ", cov_hi    CHAR(50)" ", cov_lo    CHAR(50)"
                      //                            cout << sql[120] << endl;
    try {
        int exitcode = db->execute(sql);
        if (verbose){   clog << "[sqlite3:] " << sql << " [returned:]  " << exitcode ;}
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
    // clog << "table " << id+1 << " has been created" << endl;

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
   try {
        int exitcode = db->execute(sql);
        if (verbose){   clog << "[sqlite3:] " << sql << " [returned:]  " << exitcode ;}
   } catch (exception& ex) {
         cerr << ex.what() << endl;
   }
    // clog << "`register` table has been created" << endl;
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
