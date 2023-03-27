#include <iostream> 
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;


double ** create_2d_matrix(int rows, int columns){
    double** result = new double*[rows];

    for(int i = 0; i < rows; i++){
            result[i] = new double[columns];
    }
    //Set the matrix to zero everywhere just in case
    for (int n = 0; n < rows; n++){
        for(int m = 0; m < columns; m++){
            result[n][m] = 0;
        }
    }
    return result;    
}


double * c_thomas_algorithm(double * a, double * b, double* c, double * f, int size){
    //u[0] is boundary condition
    
    double * x = new double[size]; double * y = new double[size];
    double * u = new double[size];
    
    x[0] = c[0]/b[0];
    y[0] = f[0]/b[0];

    for(int i = 1; i < size-1; i++){
        x[i] = c[i]/(b[i]-a[i]*x[i-1]);
        y[i] = (f[i] - a[i]*f[i-1])/(b[i] - a[i]*x[i-1]); 
    }
    y[size-1] = (f[size-2] - a[size-2]*f[size-3])/(b[size-2] - a[size-2]*x[size-3]); 

    u[size-1] = y[size-1];
    for(int i = size-2; i >= 0; i--){
        u[i] = y[i] - x[i]*u[i+1];
    }
    
    /*u[0] = f[0];
    x[size-2] = -a[size-1]/b[size-1];
    y[size-2] = -f[size-1]/b[size-1];
    for(int i = size-3; i > 0; i--){
        x[i] = -a[i+1]/(b[i+1] + c[i+1]*x[i+1]);
        y[i] = (f[i+1]-c[i+1]*y[i+1])/(b[i+1] + c[i+1]*x[i+1]);
        //u[i+1] = u[i]*x[i] + y[i];
    }
    x[0] = 0;
    y[0] = (f[0] - c[0]*y[1])/(b[0] + c[0]*x[1]);
    for(int i = 1; i < size-1; i++){
        u[i+1] = u[i]*x[i] + y[i];
    }*/
    return u;
}

double ** c_BCTS_heat(double** f, int rows, int columns, double tau, double h, double c){
    double** result = create_2d_matrix(rows, columns);

    double* diag = new double[columns];
    double* off_diag = new double[columns];
    double* g = new double[columns];

    off_diag[0] = 1; diag[0] = 1; off_diag[columns-1] = 1; diag[columns-1] = 1;

    for(int m = 1; m < columns-1; m++){
        off_diag[m] = -c*tau/(h*h);
        diag[m] = 2*c*tau/(h*h) + 1;
    }

    //Initial values (boundary values implicit in f)
    for(int m = 0; m < columns; m++){
        result[0][m] = h*m;
    }

    for(int n = 0; n < rows; n++){
        for(int m = 0; m < columns; m++){
            g[m] = f[n][m] - result[n][m];
        }
        result[n+1] = c_thomas_algorithm(off_diag, diag, off_diag, g, columns);
    }

    return result;
    //Create the mxm matrix A in the equation A y^(n+1) = f^n
    /*double ** A = create_2d_matrix(columns, columns);
    


    A[0][0] = 1; A[columns-1][columns-1] = 1;
    for(int m = 1; m < columns-1; m++){
        A[m][m-1] = -c*tau/(h*h);
        A[m][m] = 2*c*tau/(h*h) + 1;
        A[m][m+1] = -c*tau/(h*h);
    }*/

}

double ** c_FCTS_heat(int rows, int columns, double tau, double h, double c){
    double** result = new double*[rows];

    for(int i = 0; i < rows; i++){
            result[i] = new double[columns];
    }
    //Set the matrix to zero everywhere just in case
    for (int n = 0; n < rows; n++){
        for(int m = 0; m < columns; m++){
            result[n][m] = 0;
        }
    }
    //Initial and boundary conditions
    for (int m = 0; m < columns; m++){
        result[0][m] = m*h;
    }

    for (int n = 0; n < rows; n++){
        result[n][0] = 0;
        result[n][columns-1] = 0; 
    }

    //Explicit calculation
    for (int n = 0; n < rows; n++){
        for (int m = 0; m < columns; m++){
            result[n+1][m] = c*tau*(result[n][m-1] + result[n][m+1] - 2*result[n][m])/(h*h) + result[n][m];           
        }
    }

    return result;
}

double ** c_two_level_FDM_heat(int rows, int columns, double tau, double h, double xi){
    //Find some matrix equation :)
}

double ** c_advection(int rows, int columns, double tau, double h, double c) 
{
    double** result = new double*[rows];
    for(int i = 0; i < rows; i++){
        result[i] = new double[columns];
    }
    for (int n = 0; n < rows-1; n++){
        for(int m = 1; m < columns; m++){
            result[n][m] = 0;
        }
    }

    for (int m = 0; m < columns; m++){
        result[0][m] = m*h;
    }

    for (int n = 0; n < rows; n++){
        result[n][0] = -n*tau;
    }


    for (int n = 0; n < rows-1; n++){
        for(int m = 1; m < columns; m++){
            result[n+1][m] = -tau*c*(result[n][m] - result[n][m-1])/h + result[n][m];
        }
    }

    return result;
}



py::array_t<double> advection(int rows, int columns, double tau, double h, double c){
    double** c_array = c_advection(rows, columns, tau, h, c);
    py::array_t<double, py::array::c_style> arr({rows, columns});
    auto ra = arr.mutable_unchecked();
    for(int n = 0; n < rows; n++){
        for(int m = 0; m < columns; m++){
            ra(n,m) = c_array[n][m];
            //ra(n,m) = n+m;
        }
    } 
    return arr;
}

PYBIND11_MODULE(testing, m){
    m.def("advection", &advection, py::return_value_policy::copy);
}
