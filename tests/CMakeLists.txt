set(TEST_COMMON_FILES test_common.c test_common.h)

set(TEST_INTERFACE_NAME "test_interface_${TEST_TYPE_NAME}")
add_library(${TEST_INTERFACE_NAME} INTERFACE )
target_compile_definitions(${TEST_INTERFACE_NAME} INTERFACE INCLUDE_HEADER=<${TEST_INCLUDE_HEADER}>)
target_include_directories(${TEST_INTERFACE_NAME} INTERFACE ../source)
target_link_libraries(${TEST_INTERFACE_NAME} INTERFACE jmtx)

set(BASE_BRM_NAME "base_brm_${TEST_TYPE_NAME}")
add_executable(${BASE_BRM_NAME} matrix_ops_tests/base_brm_test.c ${TEST_COMMON_FILES})
target_link_libraries(${BASE_BRM_NAME} PRIVATE ${TEST_INTERFACE_NAME})
add_test(NAME ${BASE_BRM_NAME} COMMAND ${BASE_BRM_NAME})


set(BASE_CCS_NAME "base_ccs_${TEST_TYPE_NAME}")
add_executable(${BASE_CCS_NAME} ${TEST_COMMON_FILES} matrix_ops_tests/base_ccs_test.c)
target_link_libraries(${BASE_CCS_NAME} PRIVATE ${TEST_INTERFACE_NAME})
add_test(NAME ${BASE_CCS_NAME} COMMAND ${BASE_CCS_NAME})
