#pragma once
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>

using namespace std;
using std::initializer_list;
using std::ostream, std::setw, std::unique_ptr, std::vector;

template <typename T>
class matrix
{
public:
    // First argument: number of rows
    // Second argument: number of columns
    matrix(const size_t &, const size_t &);

    // Constructor to create a diagonal matrix
    // Argument: a vector containing the elements on the diagonal
    // Number of rows and columns is inferred automatically
    matrix(const vector<T> &);

    // Constructor to create a diagonal matrix
    // Argument: an initializer_list containing the elements on the diagonal
    // Number of rows and columns is inferred automatically
    matrix(const initializer_list<T> &);

    // Constructor to create a matrix and initialize it to the given elements
    // First argument: number of rows
    // Second argument: number of columns
    // Third argument: a vector containing the elements (flattened)
    matrix(const size_t &, const size_t &, const vector<T> &);

    // Constructor to create a matrix and initialize it to the given elements
    // First argument: number of rows
    // Second argument: number of columns
    // Third argument: an initializer_list containing the elements (flattened)
    matrix(const size_t &, const size_t &, const initializer_list<T> &);

    // Copy constructor to create a new matrix with the same elements as an existing matrix
    matrix(const matrix &);

    // Move constructor to move the elements of an existing matrix to a new matrix
    matrix(matrix &&);

    //Memeber functions
    // Overloaded operator = to assign the elements of one matrix to another matrix
    matrix &operator=(const matrix &);

    // Overloaded operator = to move the elements of one matrix to another matrix
    matrix &operator=(matrix &&);

    // Member function used to obtain (but not modify) the number of rows in the matrix
    size_t get_rows() const;

    // Member function used to obtain (but not modify) the number of columns in the matrix
    size_t get_cols() const;

    // Overloaded operator () used to access matrix elements WITHOUT range checking
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2
    // First version: returns a reference, thus allows modification of the element
    T &operator()(const size_t &, const size_t &);

    // Overloaded operator () used to access matrix elements WITHOUT range checking
    // The indices start from 0: m(0, 1) would be the element at row 1, column 2
    // Second version: does not return a reference and declared as const, does not allow modification of the element
    T operator()(const size_t &, const size_t &) const;

    // Member function used to access matrix elements WITH range checking (throws index_out_of_range)
    // The indices start from 0: m.at(0, 1) would be the element at row 1, column 2
    // First version: returns a reference, thus allows modification of the element
    T &at(const size_t &, const size_t &);

    // Member function used to access matrix elements WITH range checking (throws index_out_of_range)
    // The indices start from 0: m.at(0, 1) would be the element at row 1, column 2
    // Second version: does not return a reference and declared as const, does not allow modification of the element
    T at(const size_t &, const size_t &) const;

    // Static member function used to set the character width of the elements when printing a matrix (will be used with std::setw)
    static void set_output_width(const int &);

    // === Friend functions ===

    // Overloaded binary operator << used to easily print out a matrix to a stream
    template <typename U>
    friend ostream &operator<<(ostream &, const matrix<U> &);

    // Overloaded binary operator << used to easily print out a vector to a stream
    template <typename U>
    friend ostream &operator<<(ostream &, const vector<U> &);

    // Overloaded binary operator + used to add two matrices
    template <typename U>
    friend matrix<U> operator+(const matrix<U> &, const matrix<U> &);

    // Overloaded binary operator += used to add two matrices and assign the result to the first one
    template <typename U>
    friend matrix<U> operator+=(matrix<U> &, const matrix<U> &);

    // Overloaded unary operator - used to take the negative of a matrix
    template <typename U>
    friend matrix<U> operator-(const matrix<U> &);

    // Overloaded binary operator - used to subtract two matrices
    template <typename U>
    friend matrix<U> operator-(const matrix<U> &, const matrix<U> &);

    // Overloaded binary operator -= used to subtract two matrices and assign the result to the first one
    template <typename U>
    friend matrix<U> operator-=(matrix<U> &, const matrix<U> &);

    // Overloaded binary operator * used to multiply two matrices
    template <typename U>
    friend matrix<U> operator*(const matrix<U> &, const matrix<U> &);

    // Overloaded binary operator * used to multiply a scalar on the left and a matrix on the right
    template <typename U>
    friend matrix<U> operator*(const U &, const matrix<U> &);

    // Overloaded binary operator * used to multiply a matrix on the left and a scalar on the right
    template <typename U>
    friend matrix<U> operator*(const matrix<U> &, const U &);

    // Because I want to in a niche situation
    // Overloaded binary operator * used to multiply a vector on the left and a matrix on the right
    template <typename U>
    friend matrix<U> operator*(const vector<U> &, const matrix<U> &);

    // Overloaded binary operator * used to multiply a matrix on the left and a vector on the right
    template <typename U>
    friend matrix<U> operator*(const matrix<U> &, const vector<U> &);

    // Overloaded binary operator / used to divide a matrix on the left and a scalar on the right (essentially inverse multiplication)
    template <typename U>
    friend matrix<U> operator/(const matrix<U> &, const U &);

    // Exception

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

    // Exception to be thrown if the requested index is out of range
    class index_out_of_range
    {
    };

    // Exception for case where we need to divide by 0
    class cannot_divide_by_zero
    {
    };

private:
    // The number of rows
    size_t rows{0};

    // The number of columns
    size_t cols{0};

    // A pointer to an array storing the elements of the matrix in flattened (1-dimensional) form
    T *elements{nullptr};

    // A smart pointer to manage the memory allocated for the elements
    unique_ptr<T[]> smart{nullptr};

    // The character width of the elements when printing a matrix (will be used with std::setw)
    static int output_width;
};

// Initialize output_width to have a default value of 5
template <typename T>
int matrix<T>::output_width{5};

// Constructors

template <typename T>
matrix<T>::matrix(const size_t &input_rows, const size_t &input_cols)
    : rows(input_rows), cols(input_cols)
{
    if (rows == 0 or cols == 0)
        throw zero_size{};
    smart.reset(new T[rows * cols]);
    elements = smart.get();
}

template <typename T>
matrix<T>::matrix(const vector<T> &input_diagonal)
    : rows(input_diagonal.size()), cols(input_diagonal.size())
{
    if (rows == 0)
        throw zero_size{};
    smart.reset(new T[rows * cols]);
    elements = smart.get();
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            elements[(cols * i) + j] = ((i == j) ? input_diagonal[i] : 0);
}

template <typename T>
matrix<T>::matrix(const initializer_list<T> &input_diagonal)
    : matrix(vector<T>{input_diagonal}) {}

template <typename T>
matrix<T>::matrix(const size_t &input_rows, const size_t &input_cols, const vector<T> &input_elements)
    : rows(input_rows), cols(input_cols)
{
    if (rows == 0 or cols == 0)
        throw zero_size{};
    if (input_elements.size() != rows * cols)
        throw initializer_wrong_size{};
    smart.reset(new T[rows * cols]);
    elements = smart.get();
    for (size_t i{0}; i < rows * cols; i++)
        elements[i] = input_elements[i];
}

template <typename T>
matrix<T>::matrix(const size_t &input_rows, const size_t &input_cols, const initializer_list<T> &input_elements)
    : matrix(input_rows, input_cols, vector<T>{input_elements}) {}

template <typename T>
matrix<T>::matrix(const matrix<T> &m)
    : rows(m.rows), cols(m.cols)
{
    smart.reset(new T[rows * cols]);
    elements = smart.get();
    for (size_t i{0}; i < rows * cols; i++)
        elements[i] = m.elements[i];
}

template <typename T>
matrix<T>::matrix(matrix<T> &&m)
    : rows(m.rows), cols(m.cols)
{
    smart = move(m.smart);
    elements = smart.get();
    m.rows = 0;
    m.cols = 0;
    m.elements = nullptr;
}

//Member functions

template <typename T>
matrix<T> &matrix<T>::operator=(const matrix<T> &m)
{
    rows = m.rows;
    cols = m.cols;
    smart.reset(new T[rows * cols]);
    elements = smart.get();
    for (size_t i{0}; i < rows * cols; i++)
        elements[i] = m.elements[i];
    return *this;
}

template <typename T>
matrix<T> &matrix<T>::operator=(matrix<T> &&m)
{
    rows = m.rows;
    cols = m.cols;
    smart = move(m.smart);
    elements = smart.get();
    m.rows = 0;
    m.cols = 0;
    m.elements = nullptr;
    return *this;
}

template <typename T>
inline size_t matrix<T>::get_rows() const
{
    return rows;
}

template <typename T>
inline size_t matrix<T>::get_cols() const
{
    return cols;
}

template <typename T>
inline T &matrix<T>::operator()(const size_t &row, const size_t &col)
{
    return elements[(cols * row) + col];
}

template <typename T>
inline T matrix<T>::operator()(const size_t &row, const size_t &col) const
{
    return elements[(cols * row) + col];
}

template <typename T>
inline T &matrix<T>::at(const size_t &row, const size_t &col)
{
    if (row >= rows or col >= cols)
        throw index_out_of_range{};
    return elements[(cols * row) + col];
}

template <typename T>
inline T matrix<T>::at(const size_t &row, const size_t &col) const
{
    if (row >= rows or col >= cols)
        throw index_out_of_range{};
    return elements[(cols * row) + col];
}

template <typename T>
inline void matrix<T>::set_output_width(const int &w)
{
    output_width = w;
}

//Friend functions

template <typename T>
ostream &operator<<(ostream &out, const matrix<T> &m)
{
    if (m.rows == 0 and m.cols == 0)
        out << "()\n";
    else
    {
        for (size_t i{0}; i < m.rows; i++)
        {
            out << "( ";
            for (size_t j{0}; j < m.cols; j++)
                out << setw(m.output_width) << m(i, j) << ' ';
            out << ")\n";
        }
        out << '\n';
    }
    return out;
}

template <typename T>
ostream &operator<<(ostream &out, const vector<T> &m)
{
    if (m.size() == 0)
    {
        out << "()\n";
    }
    else
    {
        for (size_t i{0}; i < m.size(); i++)
        {
            out << "( ";
            out << m[i] << " ";
            out << ")\n";
        }
    }
    return out;
}

template <typename T>
matrix<T> operator+(const matrix<T> &a, const matrix<T> &b)
{
    if ((a.rows != b.rows) or (a.cols != b.cols))
        throw typename matrix<T>::incompatible_sizes_add{};
    matrix<T> c(a.rows, a.cols);
    for (size_t i{0}; i < a.rows; i++)
        for (size_t j{0}; j < a.cols; j++)
            c(i, j) = a(i, j) + b(i, j);
    return c;
}

template <typename T>
inline matrix<T> operator+=(matrix<T> &a, const matrix<T> &b)
{
    a = a + b;
    return a;
}

template <typename T>
matrix<T> operator-(const matrix<T> &m)
{
    matrix<T> c(m.rows, m.cols);
    for (size_t i{0}; i < m.rows; i++)
        for (size_t j{0}; j < m.cols; j++)
            c(i, j) = -m(i, j);
    return c;
}

template <typename T>
matrix<T> operator-(const matrix<T> &a, const matrix<T> &b)
{
    if ((a.rows != b.rows) or (a.cols != b.cols))
        throw typename matrix<T>::incompatible_sizes_add{};
    matrix<T> c(a.rows, a.cols);
    for (size_t i{0}; i < a.rows; i++)
        for (size_t j{0}; j < a.cols; j++)
            c(i, j) = a(i, j) - b(i, j);
    return c;
}

template <typename T>
inline matrix<T> operator-=(matrix<T> &a, const matrix<T> &b)
{
    a = a - b;
    return a;
}

template <typename T>
matrix<T> operator*(const matrix<T> &a, const matrix<T> &b)
{
    if (a.cols != b.rows)
        throw typename matrix<T>::incompatible_sizes_multiply{};
    matrix<T> c(a.rows, b.cols);
    for (size_t i{0}; i < a.rows; i++)
        for (size_t j{0}; j < b.cols; j++)
        {
            c(i, j) = 0;
            for (size_t k{0}; k < a.cols; k++)
                c(i, j) += a(i, k) * b(k, j);
        }
    return c;
}

template <typename T>
matrix<T> operator*(const matrix<T> &a, const vector<T> &b)
{
    //Just turn it into matrix multiplication duh
    //Turn Vector into transposed vector
    matrix<T> d(b.size(), 1);
    for (size_t m{0}; m < b.size(); m++)
    {
        d(m, 0) = b[m];
    }
    //Cheeky
    //Check similar to normal matrix multiplication
    if (a.cols != d.rows)
        throw typename matrix<T>::incompatible_sizes_multiply{};
    matrix<T> c(a.rows, d.cols);
    for (size_t i{0}; i < a.rows; i++)
        for (size_t j{0}; j < d.cols; j++)
        {
            c(i, j) = 0;
            for (size_t k{0}; k < a.cols; k++)
                c(i, j) += a(i, k) * d(k, j);
        }
    return c;
}

template <typename T>
matrix<T> operator*(const vector<T> &a, const matrix<T> &b)
{
    //Just turn it into matrix multiplication since that is the intent
    matrix<T> d(1, a.size());
    for (size_t m{0}; m < a.size(); m++)
    {
        d(0, m) = a[m];
    }
    //Cheeky

    if (d.cols != b.rows)
        throw typename matrix<T>::incompatible_sizes_multiply{};
    matrix<T> c(d.rows, b.cols);
    for (size_t i{0}; i < d.rows; i++)
        for (size_t j{0}; j < b.cols; j++)
        {
            c(i, j) = 0;
            for (size_t k{0}; k < b.cols; k++)
                c(i, j) += d(i, k) * b(k, j);
        }
    return c;
}

template <typename T>
matrix<T> operator*(const T &s, const matrix<T> &m)
{
    matrix<T> c(m.rows, m.cols);
    for (size_t i{0}; i < m.rows; i++)
        for (size_t j{0}; j < m.cols; j++)
            c(i, j) = s * m(i, j);
    return c;
}

template <typename T>
inline matrix<T> operator*(const matrix<T> &m, const T &s)
{
    return s * m;
}

template <typename T>
matrix<T> operator/(const matrix<T> &m, const T &s)
{
    //Ensure we are not dividing by 0
    if (s == 0)
        throw typename matrix<T>::cannot_divide_by_zero{};
    matrix<T> c(m.rows, m.cols);
    for (size_t i{0}; i < m.rows; i++)
        for (size_t j{0}; j < m.cols; j++)
            c(i, j) = m(i, j) / s;
    return c;
}

// Class to read data from text file
class read_data
{
    string fileName;

public:
    read_data(string filename) : fileName(filename)

    {
    }

    vector<double> getData();

    //Stop program completely if rows > cols
    class matrix_long_not_eligible
    {
    };

    class insufficient_input
    {
    };
};

// Transpose matrix quickly and efficiently
class transposematrix
{
public:
    transposematrix(matrix<double_t> &input_X)
        : X{input_X}
    {
    }

    matrix<double_t> transpose();

    // Exception if matrix 1x1
    class nothing_to_transpose
    {
    };

private:
    matrix<double> X;
};

// Normalize matrix
class normalize
{
public:
    normalize(matrix<double_t> input_m)
        : M{input_m}
    {
    }
    matrix<double_t> normalize_mat();

    //Exception if matrix is actually just full of 0's
    class matrix_full_of_zeros
    {
    };

private:
    matrix<double_t> M;
};

//Randomize matrix
class randomize
{
public:
    randomize(const double &inputrows, const double &inputcols)
        : R{inputrows}, C{inputcols}
    {
    }
    matrix<double_t> generate_rand_mat();

    //Exceptions to matrix randomization
    class rows_not_eligible
    {
    };

    class cols_not_eligible
    {
    };

private:
    double R;
    double C;
};

//Find MSE of OLS
class solve_OLS
{
public:
    solve_OLS(matrix<double_t> input_X, matrix<double_t> input_Y)
        : X{input_X}, Y{input_Y}
    {
    }
    double_t getOLS();

    // Exeptions
    class size_must_match_for_matrix_multiplication
    {
    };

private:
    matrix<double_t> X;
    matrix<double_t> Y;
};

//Lasso class
class lasso
{
public:
    lasso(matrix<double_t> input_X, matrix<double_t> input_Y, double_t input_lambda)
        : X{input_X}, Y{input_Y}, lambda{input_lambda}
    {
    }
    vector<double_t> coordinate_descent();

private:
    matrix<double_t> X;
    matrix<double_t> Y;
    double_t lambda;
    size_t iterations = 100; //Hard coded values because 100 is sufficient
};

//Create log space values which coordinate descent will use
template <typename T>
class Logspace
{
public:
    Logspace(T first, T last, int num, T base = 10.0) : curValue(first), base(base)
    {
        step = (last - first) / (num - 1);
    }

    T operator()()
    {
        T retval = pow(base, curValue);
        curValue += step;
        return retval;
    }

private:
    T curValue, base, step;
};

//Validate data results
class kfoldCV
{
public:
    kfoldCV(matrix<double_t> &input_X, matrix<double_t> &input_Y, matrix<double_t> &input_theta)
        : X{input_X}, Y{input_Y}, theta{input_theta}
    {
    }

    vector<double_t> findCV();

private:
    matrix<double_t> X;
    matrix<double_t> Y;
    matrix<double_t> theta;
};