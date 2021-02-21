#include <initializer_list>
#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cassert>
#include <iomanip>
#include <fstream>
#include <ctime>

#include "lassoclasses.hpp"

using namespace std;

int main()
{
    //Start timing calculations
    clock_t tStart = clock();
    try
    {
        //Create data frame of predictor variables
        read_data X_dimensions("data.txt");
        vector<double> dim_X;
        dim_X = X_dimensions.getData();

        //Fill data with randomly generated data
        randomize train_random(dim_X[0], dim_X[1]);
        matrix<double_t> Xtrain{train_random.generate_rand_mat()};

        //Normalize data
        normalize Xnorm(Xtrain);
        matrix<double_t> Xnormalized{Xnorm.normalize_mat()};

        //Create Y Matrix
        randomize Ymat(dim_X[0], 1);
        matrix<double_t> Ytrain{Ymat.generate_rand_mat()};

        //Least squares accuracy estimate
        solve_OLS olstest(Xtrain, Ytrain);
        double_t ols = olstest.getOLS();
        cout << "The OLS Mean Squared Error is " << ols << "\n";

        //Create vector of lambda's
        vector<double_t> Blambda;
        double_t start = -3; //Base to exponent -3 (aka 10^-2)
        double_t stop = 6;   //Base to exponent 6 (aka 1 million)
        int num = 60;        //Hard coded value of 60 to ensure a sensible number of iterations
        vector<double_t> vals;
        //Generate list of log spaced values to test as lambda's
        generate_n(back_inserter(vals), num, Logspace<double_t>(start, stop, num));
        for (double_t num : vals)
        {
            lasso lassotest(Xnormalized, Ytrain, num);
            vector<double_t> lassovec = lassotest.coordinate_descent();
            for (size_t i = 0; i < lassovec.size(); i++)
            {
                Blambda.push_back(lassovec[i]);
            }
        }
        //Collect all the estimators from coordinate descent
        matrix<double_t> lambdaperiteration(vals.size(), Xnormalized.get_cols(), Blambda);
        //Need to transpose this matrix
        transposematrix trans(lambdaperiteration);
        matrix<double_t> lambdamatrix{trans.transpose()};
        //Now have it so each column is a B^hat for a respective lambda

        kfoldCV loocv(Xnormalized, Ytrain, lambdamatrix);
        vector<double_t> LOO = loocv.findCV();
        //Report findings
        cout << "The lowest Penatlized Mean Squared Error in the lasso regression is " << LOO[0] << "\n";

        cout << "This is the lambda that gave us that result = " << vals[LOO[1]] << "\n";

        cout << "We used the one-standard error rule to avoid overfitting which results in an MSE of " << LOO[2] << "\n";

        cout << "This is most similar to " << vals[LOO[4]] << " which gave an MSE of " << LOO[3] << "\n";
    }
    // Exceptions
    catch (const read_data::insufficient_input &e)
    {
        cout << "Please fill textfile with enough information.";
        return -1;
    }
    catch (const randomize::rows_not_eligible &e)
    {
        cout << "This was not a viable value";
        return -1;
    }
    catch (const randomize::cols_not_eligible &e)
    {
        cout << "This was not a viable value";
        return -1;
    }
    catch (const read_data::matrix_long_not_eligible &e)
    {
        cout << "The matrix needs to be square or have columns greater than rows for this program to work";
        return -1;
    }
    catch (const normalize::matrix_full_of_zeros &e)
    {
        cout << "The matrix is full of zeros and no normalization can occur!";
        return -1;
    }
    //Report time it took to run program
    cout << "Time taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s \n";
}