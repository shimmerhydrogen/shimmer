add_executable(ptd populate_test_db.cpp)
target_link_libraries(ptd shimmer ${LINK_LIBS})
target_include_directories(ptd PUBLIC src)

#add_executable(ptd_qt1 populate_test_qt_db1.cpp)
#target_link_libraries(ptd_qt1 shimmer ${LINK_LIBS})
#target_include_directories(ptd_qt1 PUBLIC src)

add_executable(ldb load_database.cpp)
target_link_libraries(ldb shimmer ${LINK_LIBS})
target_include_directories(ldb PUBLIC src)

add_executable(testprepro testprepro.cpp)
target_link_libraries(testprepro shimmer ${LINK_LIBS})
target_include_directories(testprepro PUBLIC src)
