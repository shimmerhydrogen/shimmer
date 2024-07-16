cmake_minimum_required(VERSION 3.23)
project(sqlitetest)

find_package(SQLite3 REQUIRED)

add_executable(stest main.cpp)
target_link_libraries(stest SQLite::SQLite3)

