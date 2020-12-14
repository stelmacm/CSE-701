
#include <initializer_list>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cassert>
#include <iomanip>

#include "stelmacmclasses.hpp"

using namespace std;

vector<double> read_matsize::getData()
{
    ifstream file(fileName);

    vector<double> matsize;
    double value;

    while (file >> value)
    {
        matsize.push_back(value);
    }
    file.close();

    return matsize;
}
//Will be used to read initial penalization parameter
double read_theta::getData()
{
    ifstream fileT(filetheta);
    double thetaval;
    fileT >> thetaval;
    fileT.close();
    return thetaval;
}

//Throw matrix if the values of rows or columns is equal to zero
matrix::matrix(const size_t &input_rows, const size_t &input_cols)
    : rows(input_rows), cols(input_cols)
{
    if (rows == 0 or cols == 0)
        throw zero_size{};
    elements = vector<double>(rows * cols);
}

//Define the vector input to the undeclared be the diagonal
matrix::matrix(const vector<double> &input_diagonal)
    : rows(input_diagonal.size()), cols(input_diagonal.size())
{
    if (rows == 0)
        throw zero_size{};
    elements = vector<double>(rows * cols);
    for (size_t i{0}; i < rows; i++)
        elements[(cols * i) + i] = input_diagonal[i];
}

matrix::matrix(const initializer_list<double> &input_diagonal)
    : matrix(vector<double>{input_diagonal}) {}

matrix::matrix(const size_t &input_rows, const size_t &input_cols, const vector<double> &input_elements)
    : rows(input_rows), cols(input_cols), elements(input_elements)
{
    if (rows == 0 or cols == 0)
        throw zero_size{};
    if (input_elements.size() != rows * cols)
        throw initializer_wrong_size{};
}

matrix::matrix(const size_t &input_rows, const size_t &input_cols, const initializer_list<double> &input_elements)
    : matrix(input_rows, input_cols, vector<double>{input_elements}) {}

//Gets the rows of the matrix
size_t matrix::get_rows() const
{
    return rows;
}

//Gets the columns of the matrix
size_t matrix::get_cols() const
{
    return cols;
}

//Creating default for how matrices are stored
double &matrix::operator()(const size_t &row, const size_t &col)
{
    return elements[(cols * row) + col];
}

//Creating default for how matrices are stored
double matrix::operator()(const size_t &row, const size_t &col) const
{
    return elements[(cols * row) + col];
}

//Finds element point
double &matrix::at(const size_t &row, const size_t &col)
{
    return elements.at((cols * row) + col);
}

double matrix::at(const size_t &row, const size_t &col) const
{
    return elements.at((cols * row) + col);
}

//Defines a standard output for the matrix
ostream &operator<<(ostream &out, const matrix &m)
{
    size_t rows{m.get_rows()}, cols{m.get_cols()};
    for (size_t i{0}; i < rows; i++)
    {
        out << "( ";
        for (size_t j{0}; j < cols; j++)
            out << m(i, j) << '\t';
        out << ")\n";
    }
    out << '\n';
    return out;
}

//Defines a standard elementwise addition
matrix operator+(const matrix &a, const matrix &b)
{
    size_t rows{a.get_rows()}, cols{a.get_cols()};
    if ((rows != b.get_rows()) or (cols != b.get_cols()))
        throw matrix::incompatible_sizes_add{};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = a(i, j) + b(i, j);
    return c;
}

matrix operator+=(matrix &a, const matrix &b)
{
    a = a + b;
    return a;
}

//Creates a negative matrix
matrix operator-(const matrix &m)
{
    size_t rows{m.get_rows()}, cols{m.get_cols()};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = -m(i, j);
    return c;
}

//Defines matrix subtraction
matrix operator-(const matrix &a, const matrix &b)
{
    size_t rows{a.get_rows()}, cols{a.get_cols()};
    if ((rows != b.get_rows()) or (cols != b.get_cols()))
        throw matrix::incompatible_sizes_add{};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = a(i, j) - b(i, j);
    return c;
}

matrix operator-=(matrix &a, const matrix &b)
{
    a = a - b;
    return a;
}

//Defines element wise multiplication
matrix operator*(const matrix &a, const matrix &b)
{
    size_t a_rows{a.get_rows()}, a_cols{a.get_cols()};
    size_t b_rows{b.get_rows()}, b_cols{b.get_cols()};
    if (a_cols != b_rows)
        throw matrix::incompatible_sizes_multiply{};
    matrix c(a_rows, b_cols);
    for (size_t i{0}; i < a_rows; i++)
        for (size_t j{0}; j < b_cols; j++)
            for (size_t k{0}; k < a_cols; k++)
                c(i, j) += a(i, k) * b(k, j);
    return c;
}

//Defines a scalar times each element of the matrix
matrix operator*(const double &s, const matrix &m)
{
    size_t rows{m.get_rows()}, cols{m.get_cols()};
    matrix c(rows, cols);
    for (size_t i{0}; i < rows; i++)
        for (size_t j{0}; j < cols; j++)
            c(i, j) = s * m(i, j);
    return c;
}

//Same thing different order
matrix operator*(const matrix &m, const double &s)
{
    return s * m;
}

//Defines element wise matrix division by single scalar
matrix operator/(const matrix &m, const double &n)
{
    size_t rows{m.get_rows()};
    size_t cols{m.get_cols()};
    matrix w(rows, cols);
    //double r = (1 / n);
    for (size_t i{0}; i < rows; i++)
    {
        for (size_t j{0}; j < cols; j++)
        {
            w(i, j) = m(i, j) / n;
        }
    }

    return w;
}

//Defines matrix division as elementwise division of each element against one another
matrix operator/(const matrix &m, const matrix &n)
{
    size_t rows{m.get_rows()};
    size_t cols{m.get_cols()};
    if (rows != n.get_rows() or cols != n.get_cols())
    {
        throw matrix::incompatible_sizes_add{};
    }
    matrix w(rows, cols);
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {
            w(i, j) = m(i, j) / n(i, j);
        }
    }

    return w;
}

//Create Random values inside a matrix
matrix randomize(matrix p)
{
    random_device rd;
    mt19937 mt(rd());
    uniform_int_distribution<int> uid(0, 200); //Sampled from uniform distribution

    size_t rows{p.get_rows()}, col{p.get_cols()};
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < col; ++j)
        {
            p(i, j) = uid(mt); //Fills every value with a random one
        }
    }
    return p;
}

//Take transpose of matrix
matrix transpose(const matrix &p)
{
    size_t rows{p.get_rows()}, cols{p.get_cols()};
    matrix ptrans(cols, rows);
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {
            ptrans(j, i) = p(i, j); //Swaps all rows and columns as ya know...transpose would
        }
    }

    return ptrans;
}

//Augment matrix for Gauss Jordan method to find X^-1 (X inverse)
//Augment matrix. Takes matrix A[m][n] and puts the identity matrix right beside it
//A[m][2n] where the other matrix is I[m][n]
matrix augmentedmatrix(const matrix &N)
{
    size_t rows{N.get_rows()}, cols{N.get_cols()};
    //Creating elongated matrix of for A|I
    matrix aug(rows, 2 * cols);

    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {

            if (i == j)
            {
                aug(i, (j + cols)) = 1;
            }
            aug(i, j) = N(i, j);
        }
    }

    return aug;
}

matrix identitymatrix(const double &m)
{
    matrix N(m, m);
    size_t rows{N.get_rows()}, cols{N.get_cols()};

    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {
            //Fill diagonal with 1's
            if (i == j)
            {
                N(i, j) = 1;
            }
        }
    }

    return N;
}

//Swap rows so that we are able to complete reduced row echelon form
void swap_rows(matrix &m, size_t i, size_t j)
{
    size_t columns = m.get_cols();
    for (size_t column = 0; column < columns; ++column)
        std::swap(m(i, column), m(j, column));
}

//reduced for echelon form function
void rref(matrix &m)
{
    size_t rows = m.get_rows();
    size_t columns = m.get_cols();
    for (size_t row = 0, lead = 0; row < rows && lead < columns; ++row, ++lead)
    {
        size_t i = row;
        while (m(i, lead) == 0)
        {
            if (++i == rows)
            {
                i = row;
                if (++lead == columns)
                    return;
            }
        }
        swap_rows(m, i, row);
        if (m(row, lead) != 0)
        {
            float f = m(row, lead);
            for (size_t column = 0; column < columns; ++column)
                m(row, column) /= f;
        }
        for (size_t j = 0; j < rows; ++j)
        {
            if (j == row)
                continue;
            float f = m(j, lead);
            for (size_t column = 0; column < columns; ++column)
                m(j, column) -= f * m(row, column);
        }
    }
}

//Inverse the matrix once you get the matrix in rref form
matrix inverse(const matrix &m)
{
    assert(m.get_rows() == m.get_cols());
    if (m.get_rows() != m.get_cols())
    {
        throw matrix::matrix_not_square{};
    }

    size_t rows = m.get_rows();
    matrix tmp(rows, 2 * rows);             //temp matrix twice as big to add identity mat beside
    for (size_t row = 0; row < rows; ++row) //for every row in the set of rows
    {
        for (size_t column = 0; column < rows; ++column) //for every column in the set of rows
            tmp(row, column) = m(row, column);
        tmp(row, row + rows) = 1; //set diag to 0
    }
    rref(tmp); //reduced row echelon form
    matrix inv(rows, rows);
    for (size_t row = 0; row < rows; ++row) //for every row in the set of rows
    {
        for (size_t column = 0; column < rows; ++column) //for every column in the set of rows
            inv(row, column) = tmp(row, column + rows);  //make the right side of augmented matrix the final product
    }
    return inv;
}

matrix finderror(const matrix &X, const matrix &Y, const double &K)
{

    size_t rows = X.get_rows();
    size_t cols = X.get_cols();
    try
    {
        //Find transpose of X's
        matrix XT = transpose(X);
        //Multiply XT * X
        matrix XBXT = XT * X;
        //Create Identity matrix to penalize
        matrix I = identitymatrix(cols);
        //penalization
        matrix penalize = XBXT - (K * I);
        //cout << test;
        //Matrix inverse found using Gauss Jordan method
        matrix XBXTINV = inverse(penalize);
        //cout << XBXTINV;
        //Find inverse times transpose
        matrix invxtrans = XBXTINV * XT;
        //Find Beta hat
        matrix Beta = invxtrans * Y;
        //Find Y hat
        matrix yhat = X * Beta;
        //Define n
        double n = (1 / rows); //1/N ie 1 over number of rows
        //Define error
        matrix error = (Y - yhat);
        //MSE = 1/n e^2 = 1/n e^T * e
        matrix MSE = n * (transpose(error) * error); //Returns 0 as MSE which is odd.

        return error; //Returning error instead since it does not return as 0,
    }
    catch (const out_of_range &e)
    {
        cout << "Error: Matrix index out of range!\n";
    }
    catch (const matrix::incompatible_sizes_add &e)
    {
        cout << "Error: Two matrices can only be added or subtracted if they are of the same size!\n";
    }
    catch (const matrix::incompatible_sizes_multiply &e)
    {
        cout << "Error: Two matrices can only be multiplied if the number of columns in the first matrix is equal to the number of rows in the second matrix!\n";
    }
    catch (const matrix::matrix_not_square &e)
    {
        cout << "Error: The matrix is not square anf therefor the inverse cannot be taken!\n";
    }
}

//Find sumation of all values in a matrix
double sum(const matrix &m)
{
    size_t rows{m.get_rows()};
    size_t cols{m.get_cols()};
    double sum_of_matrix;
    double sum_init = 0;
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {
            sum_of_matrix = sum_init + m(i, j);
            sum_init = sum_of_matrix;
        }
    }
    return sum_of_matrix;
}

//Find the mean of the sum of the values in a matrix
double mean(const matrix &m)
{
    size_t rows{m.get_rows()};
    size_t cols{m.get_cols()};
    double n = (rows * cols);
    double meanofsum = sum(m) / n;

    return meanofsum;
}

//Find the sum of squares for a given
double sumofsquares(const matrix &m)
{
    size_t rows{m.get_rows()};
    size_t cols{m.get_cols()};
    double means = mean(m);
    //Expecting only a single row for Y
    double individ = 0;
    double individsqr = 0;
    double SST = 0;
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {
            individ = m(i, j) - means;
            individsqr = individ * individ;
            SST = SST + individsqr;
        }
    }

    return SST;
}

//Find the R^2 value of the data in question
double Rsqrd(const double &sst, const matrix &error)
{
    try
    {
        matrix ssr = transpose(error) * error;
        matrix FVU = ssr / sst;
        double rsqr = 1 - FVU(0, 0);
        return rsqr;
    }
    catch (const matrix::incompatible_sizes_multiply &e)
    {
        cout << "Error: Two matrices can only be multiplied if the number of columns in the first matrix is equal to the number of rows in the second matrix!\n";
    }
}

//Declare the variance inflation factor equation
double VIF(const double &rsqr)
{
    double Vif = 1 / (1 - rsqr);
    return Vif;
}

//Create a function that does all the
double findVIF(const matrix &X, const matrix &Y)
{
    double sumsqrs = sumofsquares(Y);
    matrix sserror = finderror(X, Y, 0);
    double rcubed = Rsqrd(sumsqrs, sserror);
    double varinf = VIF(rcubed);

    return varinf;
}