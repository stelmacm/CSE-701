#include <initializer_list>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cassert>
#include <iomanip>

using namespace std;

//Will be used to read initial penalization parameter
class read_theta
{
    string filetheta;

public:
    read_theta(string filenametheta) : filetheta(filenametheta)
    {
    }
    double getData();
};

class read_matsize
{
    string fileName;

public:
    read_matsize(string filename) : fileName(filename)

    {
    }

    vector<double> getData();
};

//Matrix Class
class matrix
{
public:
    // Constructor to create a zero matrix
    // First argument: number of rows
    // Second argument: number of columns
    matrix(const size_t &, const size_t &);

    // Constructor to create a diagonal matrix
    // Argument: a vector containing the elements on the diagonal
    // Number of rows and columns is inferred automatically
    matrix(const vector<double> &);

    // Constructor to create a diagonal matrix
    // Argument: an initializer_list containing the elements on the diagonal
    // Number of rows and columns is inferred automatically
    matrix(const initializer_list<double> &);

    // Constructor to create a matrix and initialize it to the given elements
    // First argument: number of rows
    // Second argument: number of columns
    // Third argument: a vector containing the elements (flattened)
    matrix(const size_t &, const size_t &, const vector<double> &);

    // Constructor to create a matrix and initialize it to the given elements
    // First argument: number of rows
    // Second argument: number of columns
    // Third argument: an initializer_list containing the elements (flattened)
    matrix(const size_t &, const size_t &, const initializer_list<double> &);

    // Member function used to obtain (but not modify) the number of rows in the matrix
    size_t get_rows() const;

    // Member function used to obtain (but not modify) the number of columns in the matrix
    size_t get_cols() const;

    // Overloaded operator () used to access matrix elements WITHOUT range checking
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2
    // First version: returns a reference, thus allows modification of the element
    double &operator()(const size_t &, const size_t &);

    // Overloaded operator () used to access matrix elements WITHOUT range checking
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2
    // Second version: does not return a reference and declared as const, does not allow modification of the element
    double operator()(const size_t &, const size_t &) const;

    // Member function used to access matrix elements WITH range checking (throws out_of_range via vector::at)
    // The indices start from 0: m.at(0, 1) would be the element at row 1, column 2
    // First version: returns a reference, thus allows modification of the element
    double &at(const size_t &, const size_t &);

    // Member function used to access matrix elements WITH range checking (throws out_of_range via vector::at)
    // The indices start from 0: m.at(0, 1) would be the element at row 1, column 2
    // Second version: does not return a reference and declared as const, does not allow modification of the element
    double at(const size_t &, const size_t &) const;

    // Exception to be thrown if the number of rows or columns given to the constructor is zero
    class zero_size
    {
    };

    // Exception to be thrown if the vector of elements provided to the constructor is of the wrong size
    class initializer_wrong_size
    {
    };

    // Exception to be thrown if two matrices of different sizes are added or subtracted
    class incompatible_sizes_add
    {
    };

    // Exception to be thrown if two matrices of incompatible sizes are multiplied
    class incompatible_sizes_multiply
    {
    };

    class matrix_not_square
    {
    };

private:
    // The number of rows
    size_t rows{0};

    // The number of columns
    size_t cols{0};

    // A vector storing the elements of the matrix in flattened (1-dimensional) form
    vector<double> elements;
};

// Overloaded binary operator << used to easily print out a matrix to a stream
ostream &operator<<(ostream &, const matrix &);

// Overloaded binary operator + used to add two matrices
matrix operator+(const matrix &, const matrix &);

// Overloaded binary operator += used to add two matrices and assign the result to the first one
matrix operator+=(matrix &, const matrix &);

// Overloaded unary operator - used to take the negative of a matrix
matrix operator-(const matrix &);

// Overloaded binary operator - used to subtract two matrices
matrix operator-(const matrix &, const matrix &);

// Overloaded binary operator -= used to subtract two matrices and assign the result to the first one
matrix operator-=(matrix &, const matrix &);

// Overloaded binary operator * used to multiply two matrices
matrix operator*(const matrix &, const matrix &);

// Overloaded binary operator * used to multiply a scalar on the left and a matrix on the right
matrix operator*(const double &, const matrix &);

// Overloaded binary operator * used to multiply a matrix on the left and a scalar on the right
matrix operator*(const matrix &, const double &);

// Overloaded binary operator / used to divide a matrix by a scalar values
matrix operator/(const matrix &, const double &);

// Overloaded binary operator / used to divide matracies element wise
matrix operator/(const matrix &, const matrix &);
