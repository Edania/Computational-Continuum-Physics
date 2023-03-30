#include <iostream> 
#include <cmath>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#define PI 3.14159265359

namespace py = pybind11;

using namespace std;

//Initialize 2d matrix using basic C++ syntax
double ** create_2d_matrix(int rows, int columns){
    double** matrix = new double*[rows];

    for(int i = 0; i < rows; i++){
            matrix[i] = new double[columns];
    }
    //Set the matrix to zero everywhere just in case
    for (int n = 0; n < rows; n++){
        for(int m = 0; m < columns; m++){
            matrix[n][m] = 0;
        }
    }
    return matrix;    
}

void c_em_fdtd(double** E_y, double** B_z, int dim_t_, int dim_x, double tau, double h){
    //Initial conditions 
    for (int m = 0; m < dim_x; m++){
        E_y[0][m] = sin(20*PI*m);
        B_z[0][m] = 0; //sin(20*PI*m);
    }


    //Periodic and reflective boundary conditions?
    for (int n = 0; n < dim_t_; n++){
        E_y[n][0] = 0;
        E_y[n][dim_t_-1] = 0; 
    }

    //Explicit calculation
    for (int n = 0; n < dim_t_; n++){
        for (int m = 0; m < dim_x; m++){
            E_y[n][m] = B_z[n][m];
            B_z[n][m] = 1;
            //E_y[n+1][m] = tau*(E_y[n][m-1] + E_y[n][m+1] - 2*E_y[n][m])/(h*h) + E_y[n][m];           
        }
    }
}

//Function connecting to python
py::list em_fdtd(int dim_t_, int dim_x, double tau, double h){
    double** c_E_y = create_2d_matrix(dim_t_, dim_x);
    double** c_B_z = create_2d_matrix(dim_t_, dim_x);
    c_em_fdtd(c_E_y, c_B_z, dim_t_, dim_x, tau, h);
    py::array_t<double, py::array::c_style> E_y_arr({dim_t_, dim_x});
    py::array_t<double, py::array::c_style> B_z_arr({dim_t_, dim_x});
    
    //Allows direct writing access
    auto E_y = E_y_arr.mutable_unchecked();
    auto B_z = B_z_arr.mutable_unchecked();

    for(int n = 0; n < dim_t_; n++){
        for(int m = 0; m < dim_x; m++){
            E_y(n,m) = c_E_y[n][m];
            B_z(n,m) = c_B_z[n][m];
        } 
    }

    py::list result;
    result.append(E_y_arr);
    result.append(B_z_arr);
    return result;
}

PYBIND11_MODULE(fdtd, m){
    m.def("em_fdtd", &em_fdtd, py::return_value_policy::copy);
}
