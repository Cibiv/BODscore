#ifndef SqliteDB_H_
#define SqliteDB_H_

#include<string>
#include <stdio.h>
#include <sqlite3.h>

using namespace std;

int const step=5;
size_t const RV = 1;

class SqliteDB{
  private:
    static int callback(void *NotUsed, int argc, char **argv, char **azColName){
      int i;
      for(i=0; i<argc; i++){
         printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
      }
      printf("\n");
      return 0;
    }

  public:
  int create_table(void) 
  int insert(char * HI, char * LO)

};

#endif 
