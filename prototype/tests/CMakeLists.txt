enable_testing()
find_package(GTest REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS} ../src)

add_executable(kernel kernel.cpp)

add_executable(hashgrid_tests hashgrid_tests.cpp)
target_link_libraries(hashgrid_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} -lpthread -ltbb)

add_executable(vector_tests vector_tests.cpp)
target_link_libraries(vector_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} -lpthread -ltbb)

add_executable(matrix_tests matrix_tests.cpp)
target_link_libraries(matrix_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} -lpthread)

#add_executable(unigrid_tests unigrid_tests.cpp)
#target_link_libraries(unigrid_tests ${GTEST_LIBRARY} ${GTEST_MAIN_LIBRARY} -lpthread -ltbb)


add_test(HashGridTests hashgrid_tests)
add_test(VectorTests vector_tests)
#add_test(UniGridTests unigrid_tests)
