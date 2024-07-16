#include <iostream>
#include <string>
#include <stdio.h>
#include <sqlite3.h>

static int callback(void *NotUsed, int argc, char **argv, char **azColName){
  int i;
  for(i=0; i<argc; i++){
    printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
  }
  printf("\n");
   return 0;
}

int
read_stations(sqlite3 *db)
{
    int rc;
    char *zErrMsg = nullptr;
    std::string zSql = "SELECT * FROM stations";

    sqlite3_stmt *stmt;
    rc = sqlite3_prepare_v2(db, zSql.c_str(), zSql.length(), &stmt, nullptr);
    if (rc) {
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
        return 1;
    }

    while (sqlite3_step(stmt) != SQLITE_DONE) {
        int num_cols = sqlite3_column_count(stmt);
        std::string s_name = (char *) sqlite3_column_text(stmt, 0);
        int s_number = sqlite3_column_int(stmt, 1);
        std::cout << s_name << " " << s_number << std::endl;
    }

    sqlite3_finalize(stmt);

    return 0;
}

int
main(int argc, char **argv) {
    sqlite3 *db;
  
    char *zErrMsg = 0;
    int rc;

    if( argc != 2 ) {
        fprintf(stderr, "Usage: %s DATABASE\n", argv[0]);
        return 1;
    }
  
    rc = sqlite3_open(argv[1], &db);
    if(rc) {
        fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(db));
        sqlite3_close(db);
        return 1;
    }

    /*
    const char *get_stations = "SELECT * FROM stations";
    rc = sqlite3_exec(db, get_stations, callback, 0, &zErrMsg);
    if( rc!=SQLITE_OK ){
        fprintf(stderr, "SQL error: %s\n", zErrMsg);
        sqlite3_free(zErrMsg);
    }
    */

    read_stations(db);

    sqlite3_close(db);
    return 0;
}