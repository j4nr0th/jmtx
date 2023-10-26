! Created by jan on 29.9.2023.
module f_crs_base
use, intrinsic :: ISO_C_Binding
implicit none

public :: binary_search, crs_multiply_with_vector, crs_get_element, crs_multiply_with_vector_row_exclude,&
        crs_multiply_with_vector_row

contains

    ! @brief
    !   Uses binary search to find a value inside a sorted array.
    !
    ! @param x value to find
    ! @param len length of array indices
    ! @param indices sorted array of indices of type int32_t
    !
    ! @return idx index of the element inside the array or zero if non-existant
    pure integer(C_INT32_T) function binary_search(x, len, indices) result(idx) bind(C, name="f_bin_search")
        integer(C_INT32_T), intent(in), value :: x
        integer(C_INT32_T), intent(in), value :: len
        integer(C_INT32_T), intent(in) :: indices(len)
        integer(C_INT32_T) :: pos, step, l, v
        l = len
        step = len / 2
        pos = 1

        do while (step /= 0)
            v = indices(pos + step)
            if (v > x) then
                l = step
                step = l / 2
            else if (v < x) then
                pos = pos + step
                l = l - step
                step = l / 2
            else
                idx = pos + step
                return
            end if
        end do

        if (indices(pos + step) == x) then
            idx = pos
        else
            idx = 0
        end if
    end function binary_search

    pure real(C_FLOAT) function crs_get_element(n_elements, dims, end_of_row_offsets, indices, values, i, j) result(v)
        integer(C_INT32_T), intent(in) :: dims, n_elements
        integer(C_INT32_T), intent(in) :: end_of_row_offsets(dims), indices(n_elements)
        real(C_FLOAT), intent(in) :: values(n_elements)
        integer(C_INT32_T), intent(in) :: i, j
        integer :: idx, i_begin, i_end
        if (i == 1) then
            i_begin = 1
            i_end = end_of_row_offsets(1)
        else
            i_begin = end_of_row_offsets(i-1) + 1
            i_end = end_of_row_offsets(i)
        end if
        v = 0

        if (i_begin - 1 == i_end) then
            return
        end if

        idx = binary_search(j - 1, i_end - i_begin + 1, indices(i_begin:i_end))
        if (idx /= 0) then
            v = values(i_begin + idx - 1)
        end if
    end function crs_get_element

    pure real(C_FLOAT) function crs_multiply_with_vector_row(n_elements, dim, indices, values, vector_x) result(v)! &
            !bind(C, name="f_crs_mul_vec_row")
        integer(C_INT32_T), intent(in), value :: n_elements, dim
        integer(C_INT32_T), intent(in) ::indices(n_elements)
        real(C_FLOAT), intent(in) :: values(n_elements), vector_x(dim)
        integer(C_INT32_T) :: n, idx

        v = 0

        do n = 1, n_elements
            idx = indices(n) + 1
            v = v + values(n) * vector_x(idx)
        end do
    end function crs_multiply_with_vector_row

    pure real(C_FLOAT) function crs_multiply_with_vector_row_exclude(n_elements, dim, indices, values, vector_x, exclude) result(v)! &
        !bind(C, name="f_crs_mul_vec_row")
        integer(C_INT32_T), intent(in), value :: n_elements, dim
        integer(C_INT32_T), intent(in) ::indices(n_elements), exclude
        real(C_FLOAT), intent(in) :: values(n_elements), vector_x(dim)
        integer(C_INT32_T) :: n, idx

        v = 0

        do n = 1, n_elements
            idx = indices(n) + 1
            if (idx /= exclude) v = v + values(n) * vector_x(idx)
        end do
    end function crs_multiply_with_vector_row_exclude

    impure integer(C_INT32_T) function crs_multiply_with_vector(n_elements, dims, end_of_row_offsets, indices, values, &
        vector_in, vector_out) result(v) bind(C, name="f_crs_mul_vec")
        integer(C_INT32_T), intent(in), value :: n_elements, dims
        integer(C_INT32_T), intent(in) ::indices(n_elements), end_of_row_offsets(dims)
        real(C_FLOAT), intent(in) :: values(n_elements), vector_in(dims)
        real(C_FLOAT), intent(out) :: vector_out(dims)
        integer(C_INT32_T) :: i_begin, i_end, n
        real(C_FLOAT) :: row_v

        !   Handeling the first row separately, because Fortran shits itselft if you manually bound check and write
        !   code that would be bad without that manual bound checking

        i_begin = 1
        i_end = end_of_row_offsets(1)
        if (i_end /= i_begin - 1) then
            row_v = crs_multiply_with_vector_row(i_end - i_begin + 1, dims, indices(i_begin:i_end), &
                    values(i_begin:i_end), vector_in)
            vector_out(1) = row_v
        end if

        do concurrent (n = 2:dims)
            i_begin = end_of_row_offsets(n - 1) + 1
            i_end = end_of_row_offsets(n)
            if (i_end /= i_begin - 1) then
                row_v = crs_multiply_with_vector_row(i_end - i_begin + 1, dims, indices(i_begin:i_end), &
                    values(i_begin:i_end), vector_in)
            else
                row_v = 0
            end if
            vector_out(n) = row_v
        end do
        v = 1

    end function crs_multiply_with_vector

end module f_crs_base