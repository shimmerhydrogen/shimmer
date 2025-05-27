add_executable(ptd populate_test_db.cpp)
target_link_libraries(ptd shimmer ${LINK_LIBS})
target_include_directories(ptd PUBLIC src)

add_executable(ldb load_database.cpp)
target_link_libraries(ldb shimmer ${LINK_LIBS})
target_include_directories(ldb PUBLIC src)

add_executable(testprepro testprepro.cpp)
target_link_libraries(testprepro shimmer ${LINK_LIBS})
target_include_directories(testprepro PUBLIC src)
