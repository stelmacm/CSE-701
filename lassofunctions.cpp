#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <memory>
#include <vector>
#include <random>
#include <fstream>
#include <algorithm>
#include <iterator>

#include "lassoclasses.hpp"

//Function that reads data
vector<double> read_data::getData()
{
    //Checking to see if file exists
    ifstream file{fileName};
    if (!file)
    {
        perror("Error opening input file");
    }

    vector<double> dataList;
    double list;

    //read data and input into vector list
    while (file >> list)
    {
        dataList.push_back(list);
    }
    file.close();

    if (dataList.size() < 2)
    {
        throw insufficient_input{};
    }
    if (dataList[0] > dataList[1])
    {
        throw matrix_long_not_eligible{};
    }

    return dataList;
}

//Transpose the matrix
matrix<double_t> transposematrix::transpose()
{
    size_t rows = X.get_rows();
    size_t cols = X.get_cols();
    //Check that there is a matrix to transpose
    if (rows == 1 && cols == 1)
    {
        throw nothing_to_transpose{};
    }

    matrix<double_t> Xtrans(cols, rows);
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; j < cols; ++j)
        {
            Xtrans(j, i) = X(i, j); //Swaps all rows and columns
        }
    }
    return Xtrans;
}

//Normalize matrix
matrix<double_t> normalize::normalize_mat()
{
    size_t rows = M.get_rows();
    size_t columns = M.get_cols();
    double_t suma2ij = 0;
    //Find the Norm of the matrix
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            double_t aij = M(i, j);
            double_t a2ij = aij * aij;
            suma2ij += a2ij;
        }
    }
    //Ensure we are not dividing by zero
    if (suma2ij == 0)
    {
        throw matrix_full_of_zeros{};
    }

    double_t anorm = sqrt(suma2ij);
    //Use division operator
    matrix<double_t> c = M / anorm;

    return c;
}

//Randomize data for matrix
matrix<double_t> randomize::generate_rand_mat()
{
    //Check there is infact a matrix
    if (R <= 0)
    {
        throw rows_not_eligible{};
    }

    if (C <= 0)
    {
        throw cols_not_eligible{};
    }

    //Generate random numbers
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_int_distribution<size_t> uid(0, 90);

    matrix<double_t> filled(R, C);

    for (size_t i{0}; i < R; ++i)
    {
        for (size_t j{0}; j < C; ++j)
        {
            filled(i, j) = uid(mt); //Fills every value with a random one
        }
    }

    return filled;
}

//Solve for original fit of data
double_t solve_OLS::getOLS()
{
    //Find matrix transpose
    transposematrix XTrans(X);
    matrix<double_t> XT = XTrans.transpose();
    //Finished finding matrix transpose
    matrix<double_t> XTX = XT * X;
    //Find matrix inverse
    if (XTX.get_rows() != XTX.get_cols())
    {
        throw solve_OLS::size_must_match_for_matrix_multiplication{};
    }
    //creates augmented matrix
    size_t rows = XTX.get_rows();
    matrix<double_t> tmp(rows, 2 * rows); //temp matrix twice as big to add identity mat beside
    //Set to 0 or else memory leaks can occur
    for (size_t row = 0; row < rows; ++row) //for every row in the set of rows
    {
        for (size_t column = 0; column < 2 * rows; ++column) //for every column in the set of rows
            tmp(row, column) = 0;                            //set
    }
    for (size_t row = 0; row < rows; ++row) //for every row in the set of rows
    {
        for (size_t column = 0; column < rows; ++column) //for every column in the set of rows
            tmp(row, column) = XTX(row, column);         //This should be XTX
        tmp(row, row + rows) = 1;                        //set diag to 1
    }

    //Reduce Row echelon the matrix
    size_t rows2 = tmp.get_rows();
    size_t columns = tmp.get_cols();
    for (size_t row = 0, lead = 0; row < rows2 && lead < columns; ++row, ++lead)
    {
        size_t i = row;
        while (tmp(i, lead) == 0)
        {
            if (++i == rows2)
            {
                i = row;
                if (++lead == columns)
                    break;
            }
        }

        //Swap rows around
        for (size_t column = 0; column < columns; ++column)
            std::swap(tmp(i, column), tmp(row, column));
        //Row Swaping ends

        if (tmp(row, lead) != 0)
        {
            double f = tmp(row, lead);
            for (size_t column = 0; column < columns; ++column)
                tmp(row, column) /= f;
        }
        for (size_t j = 0; j < rows2; ++j)
        {
            if (j == row)
                continue;
            double f = tmp(j, lead);
            for (size_t column = 0; column < columns; ++column)
                tmp(j, column) -= f * tmp(row, column);
        }
    }
    matrix<double_t> inv(rows, rows);
    for (size_t row = 0; row < rows; ++row) //for every row in the set of rows
    {
        for (size_t column = 0; column < rows; ++column) //for every column in the set of rows
            inv(row, column) = tmp(row, column + rows);  //make the right side of augmented matrix the final product
    }

    //Found matrix inverse
    //Now compute matrix estimators of B^hat
    matrix<double_t> XTXinvXT = inv * XT;
    matrix<double_t> Bhat = XTXinvXT * Y;

    //create test data
    randomize testdim(X.get_rows(), X.get_cols()); //confirm which one needs to change (needs to work with B)
    matrix<double_t> testmat{testdim.generate_rand_mat()};

    normalize testnorm(testmat);
    matrix<double_t> testdata{testnorm.normalize_mat()};

    //Find y^hat
    matrix<double_t> Yhat = testdata * Bhat;

    randomize testY(Y.get_rows(), 1);
    matrix<double_t> Ytest{testY.generate_rand_mat()};

    //Compute error
    matrix<double_t> error = Ytest - Yhat;
    matrix<double_t> errorT(error.get_cols(), error.get_rows());
    for (size_t i{0}; i < error.get_rows(); ++i)
    {
        for (size_t j{0}; j < error.get_cols(); ++j)
        {
            errorT(j, i) = error(i, j); //Swaps all rows and columns
        }
    }
    double_t N = X.get_rows();
    matrix<double_t> MSEmat = (errorT * error) / N;
    double_t MSE = MSEmat(0, 0); // convert the (1,1) matrix to a single point
    return MSE;
}

vector<double_t> lasso::coordinate_descent()
{
    size_t rows = X.get_rows();
    size_t columns = X.get_cols();

    //Have theta's inside rather than inputted since coord descent will always start at max theta
    vector<double_t> thetas(X.get_cols());
    for (size_t t = 0; t < thetas.size(); t++)
    {
        thetas[t] = 1;
    }
    //Begin coordinate descent for 100 iterations
    for (size_t i = 0; i < iterations; i++)
    {
        //For every column of the matrix (aka every predictor)
        for (size_t j = 0; j < columns; j++)
        {
            //Create column matrix
            matrix<double_t> X_j(rows, 1);

            //Need just the column j
            for (size_t m = 0; m < rows; m++)
            {
                X_j(m, 0) = X(m, j);
            }

            matrix<double_t> y_pred = X * thetas;

            matrix<double_t> X_jT(X_j.get_cols(), X_j.get_rows());
            for (size_t r{0}; r < X_j.get_rows(); ++r)
            {
                for (size_t c{0}; c < X_j.get_cols(); ++c)
                {
                    X_jT(c, r) = X_j(r, c); //Swaps all rows and columns
                }
            }
            //Calculate rho
            matrix<double_t> rho = X_jT * (Y - y_pred + thetas[j] * X_j); //Theta is a vector so vector multiplication
            double_t rhop = rho(0, 0);
            double_t rl;
            if (rhop + lambda < 0)
            {
                rl = rhop + lambda;
            }
            if (rhop - lambda > 0)
            {
                rl = rhop - lambda;
            }
            else
            {
                rl = 0;
            }
            //Store theta in a vector
            thetas[j] = rl;
        }
    }
    //Return vector of theta's
    return thetas;
}

// Writing my own sum function as the usual one doesn't add in the final element
double_t sum(vector<double> vec)
{
    size_t s{vec.size()};
    double_t sum_of_vector;
    double_t sum_init{vec[0]};
    for (size_t i = 0; i < s - 1; i++)
    {
        sum_of_vector = sum_init + vec[i + 1];
        sum_init = sum_of_vector;
    }
    return sum_of_vector;
}

//Function to find the closest value in a vector of a different value
double_t closest(vector<double_t> const &vec, double_t value)
{
    auto const it = lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end())
    {
        return -1;
    }

    return *it;
}

//CV function
vector<double_t> kfoldCV::findCV()
{
    size_t rows = X.get_rows();
    size_t cols = X.get_cols();
    double_t N = X.get_rows(); //Need a double type for division operator to work

    //Create test data to work on
    randomize testdim(rows, cols);
    matrix<double_t> testmat{testdim.generate_rand_mat()};

    normalize testnorm(testmat);
    matrix<double_t> testdata{testnorm.normalize_mat()};

    randomize testY(Y.get_rows(), 1);
    matrix<double_t> Ytest{testY.generate_rand_mat()};

    //Create empty vector
    vector<double_t> CVerrors;

    //Test data is now generated. Now we find y^hat = X*B^hat
    //So we need to get just B^hat first
    matrix<double_t> bhat(cols, 1);
    for (size_t i = 0; i < cols; i++) //for every column
    {
        for (size_t j = 0; j < rows; j++) //go through every row
        {
            bhat(j, 0) = theta(j, i); //create just a bhat column
        }

        matrix<double_t> yhat = testdata * bhat;

        matrix<double_t> error = Ytest - yhat;
        //Look into finding nice way to transpose matrices
        transposematrix errorT(error);
        matrix<double_t> errortrans{errorT.transpose()};
        matrix<double_t> errorsquared = (error * errortrans) / N; //works out conveniently
        //leave one out implies that all are being tested hence why we use 1/n

        CVerrors.push_back(errorsquared(0, 0)); //put the value of the (1,1) matrix into the vector
    }
    //create results vector
    vector<double_t> results;
    double_t minerror = *min_element(CVerrors.begin(), CVerrors.end());
    results.push_back(minerror); //results[0] = Smallest MSE

    for (size_t p = 0; p < CVerrors.size(); p++)
    {
        if (CVerrors[p] == minerror)
        {
            results.push_back(p); //results[1] = value of lambda @ smallest MSE
        }
    }
    //A little bit sloppy but works
    double_t SoS, DoS;
    double_t meanofCV = sum(CVerrors) / CVerrors.size();
    for (size_t i = 0; i < CVerrors.size(); i++)
    {
        DoS = (CVerrors[i] - meanofCV) * (CVerrors[i] - meanofCV);
        SoS += DoS;
    }
    double_t SD = sqrt(SoS / (CVerrors.size() - 1));
    double_t CVplusSD = minerror + SD;
    results.push_back(CVplusSD); // results[2] = Smallest MSE + SD

    double_t closestCV = closest(CVerrors, CVplusSD);
    results.push_back(closestCV); //results[3] = interpolation of smallest MSE + SD

    //Now need to find the index of the best lambda that we will use
    for (size_t h = 0; h < CVerrors.size(); h++)
    {
        if (CVerrors[h] == closestCV)
        {
            results.push_back(h); //result[4] = index of lambda @ smallest MSE + SD
        }
    }
    //Return the important results
    return results;
}