/*
 * ParseSNP.h
 *
 *  Created on: Sep 26, 2012
 *      Author: fritz
 */

#ifndef SQLITEDB_H_
#define SQLITEDB_H_

// #include <regex>
//#include <vector>
//#include <map>
#include <string>
#include <assert.h>     /* assert */
#include <sqlite3pp.h>
//#include <sqlite3ppext.h>
#include <memory>
#include "Coverage.h" 

using namespace std;

class SqliteDb {

private:
    bool IMMEDIATE = false;
    int verbose = 0;
    string db_file;
    string snpfile;
    string sample_label;
    bool db_flag = false;
    // unique_ptr<sqlite3pp::database> 
    sqlite3pp::database * db;
//    sqlite3pp::transaction * xct;
    string table_name;

    void exec_sql_log(char const * sql);

public:
	SqliteDb(const char * db_file){
        db = new sqlite3pp::database(db_file);
        clog << "writing to database: " << db_file << endl;
    }

	SqliteDb(const char * db_file, int verbose_)
    : verbose(verbose_) {
        db = new sqlite3pp::database(db_file);
//		buffer_size = 3000;
//		buffer = new char[buffer_size];
	}
	
    ~SqliteDb(){
        delete [] db; 
	}

    int intermediate_commit();
    int commit();
    void new_transaction();

    void print_cov_db( const Coverage & cov );
    void init_sql_table( std::string & sample_label_ );
    void init_contig_table();
    void init_register_table();
    void place_register_record( );
    void place_register_record_begin( const std::string & snpfile_ ,  const std::string & read_filename );
    void read_register_table();
    void place_contig_table_record( size_t & chr_ref, string & chr_name );
    void parseSQLite();
 
    void composite_index();
};



#endif /* PARSESNP_H_ */
