# Finding optimal penalization parameter lambda for ridge regression - Martin Stelmach - December 13, 2020

## Summary of Project

The goal of the program created is to take user input regarding the size of the data, check if the data is multicollinear by checking the variance inflation factor, and if it is, using gradient descent to find the optimal value of lambda that minimizes the mean squared error and hence finds the optimal penalization parameter.

## Explanation of Algorithm

Regression is used for a multitude of things and is one of the most important things in todays world. For this reason there are so many different things one can do with regression to try and fit data better. One of those things is ridge regression. Ridge regression is used for data has lots of multicollinearity and there for is harder to predict. Typical dataset that ridge regression is used for is datasets that contain more factors than observations, making it difficult to interpret which are more important and should have higher coefficients. The use of ridge regression on such a multicollinear dataset is effective because ridge regression implements bias into the data. What happens is a shrinkage estimator is assigned to each coefficent to differentiate the more important coefficients. This is essentially a form of penalization to factors that should have less of a weight to the regression. The penalization parameter gives us a new equation for regression coefficients in the form of B = (X'X + aI)^-1 X' Y. The tricky part of ridge regression is selecting an appropriate alpha value. There is no definitive formula for solving for this parameter. Instead what is done is several values of alpha are computed and compared for goodness of fit by seeing which value minimizes the mean squared error. This essentially guessing and checking so the value is not necessarily precise nor the most accurate. We should recognize that the function of mean square error is convex for values of alpha. Knowing this, we are able to identify that there is an absolute minima which can be found. What happens is an alpha is selected originally and the respective mean squared error is computed. From there, the mean squared error of an alpha to the left and the right, by a defined stepsize, of the original alpha is computed. The one that is lower indicates which direction alpha should continue to head in that direction. This occurs while the mean sqaure error is decreasing; once it begins to increase again, the minima has occured. 

## Implementation

### Libraries 
For this program a multitude of libraries were used, including ctime, to time the length of the program, random to randomly generate data, algorithm to use functions such as 'swap', fstream to read in data and a few others for various purposes.
  
``` cpp

#include <iostream>
#include <vector>
#include <initializer_list>
#include <random>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <cassert>
#include <iomanip>

using namespace std;
```

### Reading input

Like most regression programs, this program requires data to run! To show that this method of regression is most effective we will be simply randomizing the data based on users request. The user must input two values that will dictate the size of the matrix. This is done by creating a class that reads the file and saves the values. In goal of being object oriented, a 'class read_matsize' was created to use a constructor that reads the inputted txt file and the appropriate member function 'getData()' reads it and stores it as a vector of type double. 

```cpp
class read_matsize
{
    string fileName;

public:
    read_matsize(string filename) : fileName(filename)
    {
    }

    vector<double> getData();
};
//Used to read in matrix size for user input
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
```

### Randomize the Data

So once the size of the matrix is inputted by the user, the values of the matrix are randomized using the uniform distribution. 

```cpp
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
```

### Matrix class and its operators

The matrix class is the most important class since it defines matrixies and their respective operators. Since it was created in class there is not neccessarily too much to comment on it other than the additions made to it. Two division operators were added to do division of a scalar into a matrix and to do elementwise division of a matrix respectively. Although it linear algebra matrices are obviously not divided, element wise division will be used in the program.

```cpp

matrix operator/(const matrix &m, const double &n)
{
    size_t rows{m.get_rows()};
    size_t cols{m.get_cols()};
    matrix w(rows, cols);
    for (size_t i{0}; i < rows; ++i)
    {
        for (size_t j{0}; i < cols; ++j)
        {
            w(i, j) = m(i, j) / n;
        }
    }

    return w;
}

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

```
