add_executable(comparison_test comp_test.cpp)
target_include_directories(comparison_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(comparison_test ENGORGIO)

add_executable(sort_test sort_test.cpp)
target_include_directories(sort_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(sort_test ENGORGIO)

add_executable(hybrid_query_test hybrid_query_test.cpp)
target_include_directories(hybrid_query_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(hybrid_query_test ENGORGIO)

add_executable(relational_query_test relational_query_test.cpp)
target_include_directories(relational_query_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(relational_query_test ENGORGIO)

add_executable(vectorized_query_test vectorized_query_test.cpp)
target_include_directories(vectorized_query_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(vectorized_query_test ENGORGIO)

add_executable(topk_test topk_test.cpp)
target_include_directories(topk_test 
    PUBLIC ${CMAKE_SOURCE_DIR}/src
)
target_link_libraries(topk_test ENGORGIO)