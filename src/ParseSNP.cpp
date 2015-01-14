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
    //cmdss.add(db_file_arg);
    //
    ValueArg<string> plot_file_arg("p", "plot_file",
            "Name of the file for the data required for the Python script",
            false, "#", "string");
    //cmdss.add(plot_file_arg);
    //
    ValueArg<string> out("o", "outputfile", "File to print the distripution to",
            true, "#", "string");
    cmdss.add(out);
    //
    cmdss.xorAdd( db_file_arg, plot_file_arg );

    try {
        cmdss.parse(argc, argv); //parse arguments
        read_filename = read_file.getValue();
        output = out.getValue();
        if (output[0] == '#') {
            output = read_filename;
            output += ".mot";
        }
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

        clog << "reading the reference FASTA file" << endl;

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
    clog << "In reference " << chrs.size() <<   " contigs detected." << endl;
    fastaFile.close();
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
         
    vcfFile.getline(buffer, buffer_size);


    vector<int> tmp;
    string broken_chromosome = "";
//    snps.resize(chrs.size(), tmp);
   
    vector<Coverage> tmp_cov;
    Parser * mapped_file = 0;
    FastaParser * fasta = new FastaParser(reffile);
    string ref;
    
    FILE * plotFile = NULL;
    if (db_flag){
        cout << "REACHED !" << endl;
        db  = new sqlite3pp::database( plot_file.c_str() ) ;
        //db.reset(new sqlite3pp::database( plot_file.c_str() ) );
    } else {
        remove(plot_file.c_str());
        plotFile = fopen(plot_file.c_str(), "a");
        if (plotFile == NULL) {
            cout << "Error in printing: The file or path that you set "
                << output.c_str()
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
    int old_pos  = -100000;
    size_t id = 0;
    size_t old_id = 100000;
    ref = fasta->getChr(id);

    clog << "allocating coverage";
    Coverage * cov;
    cov = new Coverage(range);
    clog << " done" << endl;

    while (!vcfFile.eof()) {
        if (buffer[0] != '#') {
            int count = 0;
            string current_chr;
            for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
                if (count == 0 && buffer[i] != '\t') {
                    current_chr += buffer[i];
                }
                if (count == 1 && buffer[i - 1] == '\t') {

                    int pos = atoi(&buffer[i]) - 1;
                    
                    if ( chrs.count(current_chr.c_str()) > 0){
                        if  (chrs[current_chr.c_str()]!= id) { //found
//                        clog << "current_chr " << id << " num of snps " << genome[id].size() << endl;
                            id = chrs[current_chr.c_str()];
                            ref = fasta->getChr(id);
                                                 // clog << "table " << id+1 << " has been created" << endl;
                           };
                      
                    } else if (broken_chromosome.compare(current_chr)!=0){
                        broken_chromosome = current_chr;
                        cerr << endl << "Contig not found: \t" << broken_chromosome << endl;
                        continue;
                    } else {
                        continue;
                    }

                    if (id != old_id) {
                        old_pos = -100000;
                        if (db_flag) {
                            init_sql_table(id);
                        }
                    }

                   old_id = id;
                   old_pos = pos;

                   clog << "chr " <<  id << " pos " << pos;
     
                   // *cov = Coverage(pos, range);
                   cov->pos = pos;
                   cov->start_pos = pos - range;
  //                      clog << " +++ ";
  //                      clog << cov->cov_90[0];
                   cov->clear_arrays();
                   process_snp(cov, ref, mapped_file, id);
                   if (db_flag){
                       cov->print_cov_db(id, *db);
                   } else { cov->print_cov(id, plotFile); }

                   cov->estimate( read_length );

                }
                if (buffer[i] == '\t') {
                    count++;
                }
            // remove db unique_ptr
            }
        }

        vcfFile.getline(buffer, buffer_size);
    }
    clog << "VCF file has been successfully processed" << endl;
    vcfFile.close();
    if (db_flag){
        
    } else  fclose(plotFile);
}

void ParseSNP::process_snp(Coverage* cov, string & ref, Parser * mapped_file, const size_t &cc){
    int leftPos= cov->start_pos;
    int rightPos =  cov->pos + range;
    int MIN_QUALITY = 0;
    clog << ":" << cc << ":" << leftPos << "-" << rightPos ;
    // set the region of interest
    if (!mapped_file->SetRegion( (int) cc, leftPos, rightPos)){
    cerr << "cannot jump to position " << cov->pos << " on chr " << cc << endl;
    }
    
    Alignment * tmp_aln = mapped_file->parseRead(MIN_QUALITY);
    clog << endl;
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
    

void ParseSNP::print() {
    cout << "printing:" << endl;
    FILE *file;
    file = fopen(output.c_str(), "w");
    if (file == NULL) {
        cout << "Error in printing: The file or path that you set "
                << output.c_str()
                << " is not valid. It can be that there is no disc space available."
                << endl;
        exit(0);
    }

    ifstream vcfFile;
    vcfFile.open(snpfile.c_str(), ifstream::in);
    if (!vcfFile.good()) {
        cout << "Could not open file: " << snpfile.c_str()
                << endl;
        exit(0);
    }

    vcfFile.getline(buffer, buffer_size);

    while (!vcfFile.eof()) {
        if (buffer[0] != '#') {
            int count = 0;
            string chr;
            int pos = 0;
            float score = -1;
            int id = 0;
            string line;
            string key;
            for (size_t i = 0; i < buffer_size && buffer[i] != '\0' && buffer[i] != '\n'; i++) {
                if (count == 0 && buffer[i] != '\t') {
                    chr += buffer[i];
                }

                if (count == 1 && buffer[i - 1] == '\t') {

                    pos = atoi(&buffer[i]);
                    if (chrs.find(chr.c_str()) != chrs.end()) { //found
                        id = chrs[chr.c_str()];
                    } else {
                        cerr << "Chr not found: " << chr << endl;
                        exit(0);
                    }
                    key = gen_key(id, pos);
                    if (covs.find(key.c_str()) != covs.end()) {
                        score = covs[key]->estimate(read_length);

                        if (score > 1 || score < 0) {
                            score = 0;
                        }
                        covs[key]->score = score;
                    } else {
                        score = -1; //default
                    }
                }

                if (count == 5 && buffer[i - 1] == '\t' && score != -1) { //SNP Quality
                    float snp_qual = atof(&buffer[i]);
                    score = score * 255;
                    score = (score + snp_qual) / 2;
                    covs[key]->score = score;
//                    cout << "HIT: " << snp_qual << " " << score << endl;
                    ostringstream ss;
                    ss << score;
                    string s(ss.str());
                    if (score > 99) {
                        line += s.substr(0, 5);
                    } else { //if(score >9){
                        line += s.substr(0, 4);
                    }
                    line += '\t';

                } else if (count != 5 || score == -1) {
                    line += buffer[i];
                }

                if (buffer[i] == '\t') {
                    count++;
                }
            }
            fprintf(file, "%s", line.c_str());
            fprintf(file, "%c", '\n');
            line.clear();
            key.clear();
        } else {
            fprintf(file, "%s", buffer);
            fprintf(file, "%c", '\n');
        }

        vcfFile.getline(buffer, buffer_size);
    }
    vcfFile.close();
    fclose(file);

    if (plot_file[0] != '#') {
        file = fopen(plot_file.c_str(), "w");
        if (file == NULL) {
            cout << "Error in printing: The file or path that you set "
                    << plot_file.c_str()
                    << " is not valid. It can be that there is no disc space available."
                    << endl;
            exit(0);
        }

        for (map<string, Coverage*>::iterator i = covs.begin(); i != covs.end(); i++) {
    //                 print_coverage(*i, file);
             }
                 fclose(file);
        }
     cout << "Printing finished" << endl;
    for (map<string, Coverage*>::iterator i = covs.begin(); i != covs.end(); i++) {
//        delete [] (*i).second->cov_hi;
//        delete [] (*i).second->cov_lo;
        delete (*i).second;
    }
}

void ParseSNP::compute(){

}

void ParseSNP::init_sql_table( size_t & id){
        char sql[128];
        sprintf(sql, "CREATE TABLE IF NOT EXISTS coverage_%u ("  \
                           "pos INT PRIMARY KEY     NOT NULL , " \
                           "score FLOAT );", (int) id+1);
                          //    ", cov_hi    CHAR(50)" ", cov_lo    CHAR(50)"
                          //                            cout << sql[120] << endl;
        try {
             clog << "[sqlite3 :] " << sql << " >>> " ;
             clog << db->execute(sql) << endl;
        } catch (exception& ex) {
             cerr << ex.what() << endl;
        }
        // clog << "table " << id+1 << " has been created" << endl;

}
