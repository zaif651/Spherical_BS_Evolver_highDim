#ifndef MATHUTILS_CPP_
#define MATHUTILS_CPP_

#include<math.h>
#include <iostream>
#include<vector>


using std::cout;
using std::cerr;
using std::endl;


//shorthand to force n into finite range (inclusive) of integers
int bound(int n, int lower, int upper)
{
    if (n < lower) n = lower;
    if (n > upper) n = upper;
    return n;
}

//given values of Y at X, return value of Lagrange polynomial passing through each at point x
double lagrange_interp(double x, const std::vector<double>& X, const std::vector<double>& Y)
{
    if (X.size() != Y.size())
    {
        cout << "ERROR: lagrange_interp called with point and value vectors of unqeual size" << endl;
        exit(1);
    }

    if (x == X[0])
        return Y[0];

    int n_points = X.size();

    for (int j = 1; j < n_points; j++)
    {
         if (X[j] <= X[j - 1])
        {
            cout << "ERROR: lagrange_interp called with invalid position vector; entries must be in strictly ascending order" << endl;
            exit(1);
        }

        if (x == X[j])
            return Y[j];
    }

    double y = 0.; //value we will return

    for (int j = 0; j < n_points; j++ ) //outer loop for sum
    {
        double term_j = Y[j];

        for (int k = 0; k < n_points; k++ ) //inner loop for product
        {
            if (k == j)
                continue;

            term_j *= ((x - X[k]) / (X[j] - X[k]) );
        }

        y += term_j;
    }

    return y;
}

//interpolates f at r using 4 gridpoints starting at j0*dr (r should be between (j0+1)*dr and (j0+2)*dr)
double cubic_interp(double r, double f0, double f1, double f2, double f3, int j0, double dr)
{
    double j1 = j0 + 1;
    double j2 = j1 + 1;
    double j3 = j2 + 1;

    std::vector<double> R = {j0*dr,j1*dr,j2*dr, j3*dr};
    std::vector<double> F = {f0,f1,f2,f3};

    return lagrange_interp(r, R, F);


   // return -(r - j1*dr) * (r - j2*dr ) * (r - j3*dr ) * f0 / (6. * dr * dr * dr) +  (r - j0*dr) * (r - j2*dr ) * (r - j3*dr ) * f1 / (2. * dr * dr * dr)
     //      -(r - j0*dr) * (r - j1*dr ) * (r - j3*dr ) * f2 / (2. * dr * dr * dr) + (r - j0*dr) * (r - j1*dr ) * (r - j2*dr ) * f3 / (6. * dr * dr * dr);

}

//returns the value of the derivative estimated with a five point centered stencil at f3 at given order (i.e. order = 2 -> 2nd derivative) up to 4
double fivePointDeriv(double step, int order, double f1, double f2, double f3, double f4, double f5)
{
    switch (order)
    {
        case 0:
            return f3;
        case 1:
            return ((1.0 / 12.0) * f1 - (2.0 / 3.0) * f2 + (2.0 / 3.0) * f4 - (1.0 / 12.0) * f5) / step; //first deriv, error at O(step^5)
        case 2:
            return ( (-1.0 / 12.0) * f1 + (4.0 / 3.0) * f2 - (5.0 / 2.0) * f3 + (4.0 / 3.0) * f4 - (1.0 / 12.0) * f5) / pow(step,2); //second deriv, error at O(step^4)
        case 3:
            return ((-1.0 / 2.0) * f1 + f2 -f4 + (1.0 / 2.0) * f5) / pow(step,3); //third deriv, error at O(step^2)
        case 4:
            return (f1 -4. * f2 + 6. * f3 - 4. * f4 + f5) / pow(step,4); //fourth deriv, error at O(step^2)
        default:
            printf("ERROR: invalid derivative order requested");
            abort();

    }
}

// Natural cubic spline solver
void computeSplineCoefficients(const std::vector<double>& y, std::vector<double>& M) {
    int n = y.size();
    M.assign(n, 0.0);
    if (n < 3) return;

    std::vector<double> a(n - 1), b(n), c(n - 1), d(n);

    // Setup the tridiagonal system
    for (int i = 1; i < n - 1; ++i) {
        a[i - 1] = 1.0;          // sub-diagonal
        b[i] = 4.0;              // main diagonal
        c[i] = 1.0;              // super-diagonal
        d[i] = 6.0 * (y[i + 1] - 2.0 * y[i] + y[i - 1]); // RHS
    }

    // Natural spline: second derivatives at ends = 0
    b[0] = b[n - 1] = 1.0;
    d[0] = d[n - 1] = 0.0;

    // Forward elimination
    for (int i = 1; i < n; ++i) {
        double m = a[i - 1] / b[i - 1];
        b[i] -= m * c[i - 1];
        d[i] -= m * d[i - 1];
    }

    // Back substitution
    M[n - 1] = d[n - 1] / b[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        M[i] = (d[i] - c[i] * M[i + 1]) / b[i];
    }
}


// Evaluates the spline at a given index
double evaluateSpline(const std::vector<double>& y, const std::vector<double>& secondDerivatives, int i, double x) {
    double h = 1.0; // uniform spacing
    double a = (i + 1 - x) / h;
    double b = (x - i) / h;
    return a * y[i] + b * y[i + 1] +
           ((a * a * a - a) * secondDerivatives[i] +
            (b * b * b - b) * secondDerivatives[i + 1]) * (h * h) / 6.0;
}

double evaluateSplineSegment(const std::vector<double>& x, const std::vector<double>& y,
                              const std::vector<double>& M, int seg, double xi) {
    double h = x[seg + 1] - x[seg];
    double a = (x[seg + 1] - xi) / h;
    double b = (xi - x[seg]) / h;

    return a * y[seg] + b * y[seg + 1] +
           ((a * a * a - a) * M[seg] + (b * b * b - b) * M[seg + 1]) * (h * h) / 6.0;
}

#endif /* MATHUTILS_CPP_ */

