list(APPEND CONVERSION_SOURCE_FILES
        source/conversion/crs_conversion.c
        source/conversion/ccs_conversion.c
        source/conversion/cds_conversion.c
        source/conversion/brm_conversion.c)

list(APPEND CONVERSION_HEADER_FILES
        include/jmtx/conversion/crs_conversion.h
        include/jmtx/conversion/ccs_conversion.h
        include/jmtx/conversion/cds_conversion.h
        include/jmtx/conversion/brm_conversion.h)

if (CMAKE_PROJECT_NAME STREQUAL "jmtx")
add_executable(float_to_double_conversion_test
        tests/converison/float_to_double_conversion_test.c
        tests/float/test_common.c tests/double/test_common.c)
target_link_libraries(float_to_double_conversion_test PRIVATE jmtx)
add_test(NAME float_to_double_conversion_test COMMAND float_to_double_conversion_test)
endif()
