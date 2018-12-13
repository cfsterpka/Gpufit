#ifndef GPUFIT_QUADRATIC1D_CUH_INCLUDED
#define GPUFIT_QUADRATIC1D_CUH_INCLUDED

__device__ void calculate_quadratic1d(
    REAL const * parameters,
    int const n_fits,
    int const n_points,
    REAL * value,
    REAL * derivative,
    int const point_index,
    int const fit_index,
    int const chunk_index,
    char * user_info,
    std::size_t const user_info_size)
{
    // indices

    REAL * user_info_float = (REAL*) user_info;
    REAL x = 0;
    if (!user_info_float)
    {
        x = point_index;
    }
    else if (user_info_size / sizeof(REAL) == n_points)
    {
        x = user_info_float[point_index];
    }
    else if (user_info_size / sizeof(REAL) > n_points)
    {
        int const chunk_begin = chunk_index * n_fits * n_points;
        int const fit_begin = fit_index * n_points;
        x = user_info_float[chunk_begin + fit_begin + point_index];
    }

    // value
    value[point_index] = parameters[0]*(x - parameters[1])*(x - parameters[1]) + parameters[2]; //vertex form

    // derivatives

    REAL * current_derivatives = derivative + point_index;

    //vertex form derivatives
    current_derivatives[0 * n_points] = (x - parameters[1])*(x - parameters[1]); //why both of these negative?
    current_derivatives[1 * n_points] = 2.0*parameters[0]*(parameters[1] - x); //why both of these negative?
    current_derivatives[2 * n_points] = 1.0;
}

#endif
