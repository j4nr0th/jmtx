set(TEST_COMMON_FILES test_common.c test_common.h)

set(TEST_INTERFACE_NAME "test_interface_${TEST_TYPE_NAME}")
add_library(${TEST_INTERFACE_NAME} INTERFACE)
target_compile_definitions(${TEST_INTERFACE_NAME} INTERFACE INCLUDE_HEADER=<${TEST_INCLUDE_HEADER}>)
target_include_directories(${TEST_INTERFACE_NAME} INTERFACE ../source)
target_link_libraries(${TEST_INTERFACE_NAME} INTERFACE jmtx)

function(get_test_name OUTPUT_VAR FILENAME)
    get_filename_component(TEST_NAME ${FILENAME} NAME_WE)
    set(${OUTPUT_VAR} "${TEST_NAME}_${TEST_TYPE_NAME}" PARENT_SCOPE)
endfunction()

file(GLOB MATRIX_TEST_FILES "matrix_ops_tests/*.c")
file(GLOB SOLVER_TEST_FILES "solver_tests/*.c")
# file(GLOB TEST_FILES "conversion/*.c")
list(APPEND TEST_FILES
         ${MATRIX_TEST_FILES}
         ${SOLVER_TEST_FILES}

)

foreach (TEST_FILE ${TEST_FILES})
    get_test_name(TEST_NAME ${TEST_FILE})
    add_executable(${TEST_NAME} ${TEST_COMMON_FILES} ${TEST_FILE})
    target_link_libraries(${TEST_NAME} PRIVATE ${TEST_INTERFACE_NAME})
    add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach ()
