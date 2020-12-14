//Program will do ridge regression which will add a penalization parameter theta
//In order to optimize the following objective
//Objective = RSS + theta * Sum of squares coeff
//The theta gets determined by determining the minimized MSE
//This is determined through gradient descent
//Since the function of MSE is concave up (convex) we can solve for the minima
//by selecting any point theta and then increaing it or decreasing it therefore
//we are moving the MSE along the convex function until we reach the minima through iterative steps
//Once the minima is passed and the MSE increases again, we know that we have achiecve the desired MSE

#include <iostream>
#include <vector>
#include <initializer_list>
#include <random>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <iomanip>

#include "stelmacmfinal.hpp"

using namespace std;

int main()
{
    clock_t tStart = clock();

    //read_theta theta{"theta.txt"};
    //double T{theta.getData()};

    read_matsize matsize("matsize.txt");
    vector<double> matrixsize{matsize.getData()};

    try
    {

        //initialize matrix of size based on input
        matrix M(matrixsize[0], matrixsize[1]);
        //Fill matrix with random generated data
        matrix X = randomize(M);
        //Create set of Y
        matrix N(matrixsize[0], 1);
        matrix Y = randomize(N);

        //Create VIF verifier.
        double valueinflation = findVIF(X, Y);
        //Tell the user is there is collinearity or if there wont be any
        if (valueinflation < 10)
        {
            cout << "We have okay data. The program will return the appropriate MSE of the regression by nature of defaulting to no penalization.\n";
        }
        else if (valueinflation > 10)
        {
            cout << "We have lots of collinearity in our data. The program will find the optimal penalization.\n";
        }

        //Initialize penalization parameter
        double k = 1.1;
        matrix error = finderror(X, Y, k);

        double n = (1 / matrixsize[0]); //1/N ie 1 over number of rows
        //Calculate MSE
        matrix MSE = n * (transpose(error) * error);
        //Define tolerance that seems reasonable
        double tol = 1e-11;
        //MSE = 1/n (e^T * e)
        //MSE = 1/n (y^T * y - 2*\beta^T*x^T*y + \beta^T*x^T*x*\beta)

        //Take MSE and then do stepwise gradient descent to find minimum MSE
        if (MSE(0, 0) < tol)
        {
            cout << MSE << " is already the min";
        }
        else if (MSE(0, 0) > tol)
        {
            //Initialize which direction we will be going in
            //Ie where we are on the convex hull
            double stepsize = 0.1;
            double greaterk = k + stepsize; //1.6
            double lesserk = k - stepsize;  //1.4
            //Defined new penalization parameters to check
            //Now compare the MSE of the two
            matrix greatererror = finderror(X, Y, greaterk);
            matrix MSEgreat = n * (transpose(greatererror) * greatererror);
            matrix lessererror = finderror(X, Y, lesserk);
            matrix MSEless = n * (transpose(lessererror) * lessererror);
            matrix prevMSE = MSE;
            if (MSEgreat(0, 0) > MSEless(0, 0))
            {
                while (prevMSE(0, 0) > MSEless(0, 0)) // Looking for change of direction
                {
                    prevMSE = MSEless;
                    //Now we shift the penalization parameter by a step size and recalculate MSE
                    lesserk = lesserk - stepsize;
                    matrix lessererror = finderror(X, Y, lesserk);
                    matrix MSEless = n * (transpose(lessererror) * lessererror);
                }
                cout << "The optimal MSE has been reached by going left :" << MSEless;
                cout << "The optimal value of lambda is :" << lesserk << "\n";
            }
            if (MSEgreat(0, 0) < MSEless(0, 0))
            {

                while (prevMSE(0, 0) > MSEgreat(0, 0)) //Looking till change of direction
                {
                    prevMSE = MSEgreat;
                    //Now we shift the penalization parameter by a step size and recalculate MSE
                    greaterk = greaterk + stepsize;
                    matrix greatererror = finderror(X, Y, greaterk);
                    matrix MSEgreat = n * (transpose(greatererror) * greatererror);
                }
                cout << "The optimal MSE has been reached by going right :" << MSEgreat;
                cout << "The optimal value of lambda is :" << greaterk << "\n";
            }
        }
    }
    catch (const out_of_range &e)
    {
        cout << "Error: Matrix index out of range!\n";
    }
    catch (const matrix::zero_size &e)
    {
        cout << "Error: Cannot create a matrix with zero rows or columns!\n";
    }
    catch (const matrix::initializer_wrong_size &e)
    {
        cout << "Error: Initializer elements do not match the expected number of elements!\n";
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
    };

    cout << "Time taken: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << "s" << '\n';
}