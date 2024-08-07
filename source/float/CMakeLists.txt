
list(APPEND MATRIX_TYPE_SOURCE_FILES
        source/float/matrices/sparse_row_compressed.c
        source/float/matrices/sparse_column_compressed.c
        source/float/matrices/sparse_conversion.c
        source/float/matrices/sparse_multiplication.c
        source/float/matrices/band_row_major.c
        source/float/matrices/sparse_diagonal_compressed.c
        source/float/matrices/dense_row_major.c
)

list(APPEND MATRIX_TYPE_HEADER_FILES
        include/jmtx/float/matrices/sparse_row_compressed.h
        include/jmtx/float/matrices/sparse_column_compressed.h
        include/jmtx/float/matrices/sparse_conversion.h
        include/jmtx/float/matrices/sparse_conversion_safe.h
        include/jmtx/float/matrices/sparse_multiplication.h
        include/jmtx/float/matrices/sparse_multiplication_safe.h
        include/jmtx/float/matrices/band_row_major.h
        include/jmtx/float/matrices/band_row_major_safe.h
        include/jmtx/float/matrices/sparse_diagonal_compressed.h
        include/jmtx/float/matrices/sparse_diagonal_compressed_safe.h
        include/jmtx/float/matrices/dense_row_major.h
)

list(APPEND INTERNAL_HEADER_FILES
        source/float/matrices/sparse_row_compressed_internal.h
        source/float/matrices/sparse_column_compressed_internal.h
        source/float/matrices/band_row_major_internal.h
        source/float/matrices/sparse_diagonal_compressed_internal.h
        source/float/solvers/gmres_internal.h
)

list(APPEND INTERNAL_SOURCE_FILES
        source/float/matrices/basic_io.c
        source/float/matrices/basic_io.h)

list(APPEND SOLVER_SOURCE_FILES
        source/float/solvers/jacobi_point_iteration.c
        source/float/solvers/gauss_seidel_iteration.c
        source/float/solvers/bicgstab_iteration.c
        source/float/solvers/conjugate_gradient_iteration.c
        source/float/solvers/lu_solving.c
        source/float/solvers/cholesky_solving.c
        source/float/solvers/generalized_minimum_residual_iteration.c
        source/float/solvers/recursive_generalized_minimum_residual_iteration.c
)

list(APPEND SOLVER_HEADER_FILES
        include/jmtx/float/solvers/jacobi_point_iteration.h
        include/jmtx/float/solvers/gauss_seidel_iteration.h
        include/jmtx/float/solvers/bicgstab_iteration.h
        include/jmtx/float/solvers/conjugate_gradient_iteration.h
        include/jmtx/float/solvers/lu_solving.h
        include/jmtx/float/solvers/cholesky_solving.h
        include/jmtx/float/solvers/generalized_minimum_residual_iteration.h
        include/jmtx/float/solvers/recursive_generalized_minimum_residual_iteration.h
)

list(APPEND DECOMPOSITIONS_SOURCE_FILES
        source/float/decompositions/incomplete_lu_decomposition.c
        source/float/decompositions/incomplete_cholesky_decomposition.c
        source/float/decompositions/band_lu_decomposition.c
        source/float/decompositions/dense_lu_decomposition.c
)

list(APPEND DECOMPOSITIONS_HEADER_FILES
        include/jmtx/float/decompositions/incomplete_lu_decomposition.h
        include/jmtx/float/decompositions/incomplete_cholesky_decomposition.h
        include/jmtx/float/decompositions/band_lu_decomposition.h
        include/jmtx/float/decompositions/dense_lu_decomposition.h
)

if (CMAKE_PROJECT_NAME STREQUAL "jmtx")

add_executable(base_crs_test_float tests/float/matrix_ops_tests/base_crs_test.c tests/float/test_common.c)
target_link_libraries(base_crs_test_float jmtx m)
add_test(NAME base_crs_float COMMAND base_crs_test_float)

add_executable(base_ccs_test_float tests/float/matrix_ops_tests/base_ccs_test.c tests/float/test_common.c)
target_link_libraries(base_ccs_test_float jmtx m)
add_test(NAME base_ccs_float COMMAND base_ccs_test_float)

add_executable(old_crs_demo_float tests/float/really_old_tests/run_crs_tests.c tests/float/test_common.c)
target_link_libraries(old_crs_demo_float PRIVATE jmtx m)
add_test(NAME old_crs_float COMMAND old_crs_demo_float)

add_executable(old_ccs_demo_float tests/float/really_old_tests/run_ccs_tests.c tests/float/test_common.c)
target_link_libraries(old_ccs_demo_float PRIVATE jmtx m)
add_test(NAME old_ccs_float COMMAND old_ccs_demo_float)

add_executable(add_crs_test_float tests/float/matrix_ops_tests/crs_addition_test.c tests/float/test_common.c)
target_link_libraries(add_crs_test_float PRIVATE jmtx m)
add_test(NAME addition_crs_float COMMAND add_crs_test_float)

add_executable(jacobi_point_test_float tests/float/solver_tests/jacobi_test.c tests/float/test_common.c source/float/solvers/jacobi_point_iteration.c)
target_link_libraries(jacobi_point_test_float PRIVATE jmtx m)
add_test(NAME jacobi_float COMMAND jacobi_point_test_float)

add_executable(omp_jacobi_test_float tests/float/solver_tests/omp_jacobi_test.c tests/float/test_common.c)
target_link_libraries(omp_jacobi_test_float PRIVATE jmtx)
add_test(NAME omp_jacobi_float COMMAND omp_jacobi_test_float)

add_executable(gauss_test_float tests/float/solver_tests/gauss_test.c tests/float/test_common.c)
target_link_libraries(gauss_test_float PRIVATE jmtx)
add_test(NAME gauss_float COMMAND gauss_test_float)

add_executable(omp_build_test_float tests/float/matrix_ops_tests/omp_build_test.c tests/float/test_common.c)
target_link_libraries(omp_build_test_float PRIVATE jmtx)
add_test(NAME omp_build_float COMMAND omp_build_test_float)

add_executable(omp_cg_test_float tests/float/solver_tests/omp_cg_test.c tests/float/test_common.c)
target_link_libraries(omp_cg_test_float PRIVATE jmtx)
add_test(NAME omp_cg_float COMMAND omp_cg_test_float)

add_executable(incomplete_lu_test_float tests/float/solver_tests/inclu_test.c tests/float/test_common.c)
target_link_libraries(incomplete_lu_test_float PRIVATE jmtx)
add_test(NAME solve_ilu_float COMMAND incomplete_lu_test_float)

add_executable(ilu_test_float tests/float/solver_tests/ilu_test.c tests/float/test_common.c)
target_link_libraries(ilu_test_float PRIVATE jmtx)
target_compile_options(ilu_test_float PRIVATE -fopenmp)

add_executable(icho_test_float tests/float/solver_tests/incho_test.c tests/float/test_common.c)
target_link_libraries(icho_test_float PRIVATE jmtx)
add_test(NAME solve_icho_float COMMAND icho_test_float)

add_executable(icpcg_test_float tests/float/solver_tests/icpcg_test.c tests/float/test_common.c)
target_link_libraries(icpcg_test_float PRIVATE jmtx)
add_test(NAME test_icpcg_float COMMAND icpcg_test_float)

add_executable(lu_solve_test_float tests/float/matrix_ops_tests/lu_solve_test.c tests/float/test_common.c)
target_link_libraries(lu_solve_test_float PRIVATE jmtx)
add_test(NAME lu_solve_test_float COMMAND lu_solve_test_float)

add_executable(base_io_test_float tests/float/matrix_ops_tests/io_test.c tests/float/test_common.c)
target_link_libraries(base_io_test_float PRIVATE jmtx)
add_test(NAME base_io_test_float COMMAND base_io_test_float)

add_executable(base_brm_test_float tests/float/matrix_ops_tests/base_brm_test.c tests/float/test_common.c)
target_link_libraries(base_brm_test_float PRIVATE jmtx)
add_test(NAME base_brm_test_float COMMAND base_brm_test_float)

add_executable(lu_solve_test_brm_float tests/float/matrix_ops_tests/lu_solve_test_brm.c tests/float/test_common.c)
target_link_libraries(lu_solve_test_brm_float PRIVATE jmtx)
add_test(NAME lu_solve_test_brm_float COMMAND lu_solve_test_brm_float)

add_executable(brm_decompose_test_float tests/float/solver_tests/brm_decompose_test.c tests/float/test_common.c)
target_link_libraries(brm_decompose_test_float PRIVATE jmtx)
add_test(NAME brm_decompose_test_float COMMAND brm_decompose_test_float)

add_executable(brm_lu_improvement_test_float tests/float/solver_tests/brm_lu_improvement_test.c tests/float/test_common.c)
target_link_libraries(brm_lu_improvement_test_float PRIVATE jmtx)
add_test(NAME brm_lu_improvement_test_float COMMAND brm_lu_improvement_test_float)

add_executable(base_cds_test_float tests/float/matrix_ops_tests/base_cds_test.c tests/float/test_common.c)
target_link_libraries(base_cds_test_float PRIVATE jmtx)
add_test(NAME base_cds_test_float COMMAND base_cds_test_float)

add_executable(cg_cds_test_float tests/float/solver_tests/cg_cds_test.c tests/float/test_common.c)
target_link_libraries(cg_cds_test_float PRIVATE jmtx)
add_test(NAME cg_cds_float COMMAND cg_cds_test_float)

add_executable(icho_cds_test_float tests/float/solver_tests/incho_cds_test.c tests/float/test_common.c)
target_link_libraries(icho_cds_test_float PRIVATE jmtx)
add_test(NAME icho_cds_test_float COMMAND icho_cds_test_float)

add_executable(bicgstab_cds_test_float tests/float/solver_tests/bicgstab_cds_test.c tests/float/test_common.c)
target_link_libraries(bicgstab_cds_test_float PRIVATE jmtx)
add_test(NAME bicgstab_cds_test_float COMMAND bicgstab_cds_test_float)

add_executable(gmres_test_float tests/float/solver_tests/gmres_test.c tests/float/test_common.c)
target_link_libraries(gmres_test_float PRIVATE jmtx)
add_test(NAME gmres_test_float COMMAND gmres_test_float)

add_executable(gmres_cds_test_float tests/float/solver_tests/gmres_cds_test.c tests/float/test_common.c)
target_link_libraries(gmres_cds_test_float PRIVATE jmtx)
add_test(NAME gmres_cds_test_float COMMAND gmres_cds_test_float)

add_executable(gmresr_cds_test_float tests/float/solver_tests/gmresr_cds_test.c tests/float/test_common.c)
target_link_libraries(gmresr_cds_test_float PRIVATE jmtx)
add_test(NAME gmresr_cds_test_float COMMAND gmresr_cds_test_float)

add_executable(base_drm_test_float tests/float/matrix_ops_tests/base_drm_test.c tests/float/test_common.c)
target_link_libraries(base_drm_test_float PRIVATE jmtx)
add_test(NAME base_drm_test_float COMMAND base_drm_test_float)

if (NOT CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    add_executable(drm_decompose_test_float tests/float/solver_tests/drm_decompose_test.c tests/float/test_common.c)
    target_link_libraries(drm_decompose_test_float PRIVATE jmtx)
    add_test(NAME drm_decompose_test_float COMMAND drm_decompose_test_float)
endif()

endif()