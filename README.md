# Comparing optimal penalization parameter lambda for Lasso Regression found using Coordinate Descent - Martin Stelmach - February 21, 2021

## Summary of Project

The goal of the program created is to take user find the optimal value of lambda for which we do penalized regression. This is done by taking a user defined randomly generated dataset and doing coordinate descent along log-spaced lambda values until the penalization has been performed. From there the optimal lambda is selected via cross validation and the penalized mean squared error is compared to the standard OLS proving that the lasso is a better fit but avoids overfitting, as we will further see. It is worth mentioning that lasso optimization is not a closed form solution so a forward stagewise algorithm is implemented. 

## Explanation of Algorithm

Regression is one of the most important tools in statistics. Often times, we fit data using a simple least squares approximation, however, this is not always sufficient. As data becomes larger and we have more predictors, we start to encounter more variance in our estimates for $$\beta$$. To avoid this, we apply a penalization parameter which regularizes the coefficients and reduces variance. By implementing L1 penalization, we are ensuring our coefficients do not blow up and get really big. Once we have done the coordinate descent with the lambda in question, we use the one standard rule for model selection. This means that the optimal lambda is selected one standard deviation from the optimal that is found via cross validation. This is because the optimal lambda in cross validation can often times be an overfit and result in certain estimators being very close to zero but not zero, and shifting by a standard deviation would "push" those over.

## Implementation

### Libraries 

The program created used many functions of the standard template libraries. Since so many were used I decided to use `using namespace std` instead of individually calling them based on the use. The algorithms standard library allows the program to use functions such as generate_n, which allowed me to create the log spaced sequence, while the vectors class allows the program to use vectors as containers. 

``` cpp

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
```

### Reading input

The program begins with reading the text file in the document to create the matrix size in order the generate the date. This is done in a class that reads each individual character in the file and stores it as a vector that can be easily accessed. 

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

### Randomly Generate data

Using the `random` library, the program is able to generate random values that will fill the matrix the user has designed. The values are selected from a uniform distribution and are then used to fill in the matrix that has been created. An X matrix and Y matrix are created respectively. The randomize class is also used further on in the cross validation class as we need to generate testing data to perform validation on. The randomize class also has exceptions if the user has somehow created a matrix with 0 rows or 0 columns, an error will occur. 

### Normalization

After the data is randomly generated, it is beneficial for Lasso that the data is normalized. In order to do this, a division operator had to be created. This allowed us to divide every individual component of the matrix by a given scalar. This normalization is what helps the program run more effectively. 

``` cpp
template <typename T>
matrix<T> operator/(const matrix<T> &m, const T &s)
{
    matrix<T> c(m.rows, m.cols);
    for (size_t i{0}; i < m.rows; i++)
        for (size_t j{0}; j < m.cols; j++)
            c(i, j) = m(i, j) / s;
    return c;
}
```

### Ordinary Least Squares Class

After the matrices are created and normalized, they are then used to find $\beta^{hat}$. This is how we would normally fit data. When performing matrix multiplication, we are always checking for exceptions to see that matrices can be multiplied to one another and transposed. In the class an augmented matrix is created and reduced row echelon is performed, after that the matrix is inverted and we have a matrix of our coefficients. As part of the class, the mean squared error of the prediction of the ordinary least squares was calculated to compare to the Lasso penalization. 

### Coordinate Descent

The coordinate descent begins by creating a list of lambda's. This list is created from a log spaced generated list. The list's start, end and size are hard coded since there is not particularly a need to specify these before hand, as creating to long a list can lead to aggresive overfitting. The following class defines the values in the list of these inputted values. This list is then created using the `generate_n` function from `<algorithms>` standard library. This ensures we are testing a good proportion of lambda values without testing too many that are very similar. 

```cpp

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
```

From there, the program tests every lambda in the list, we cycle through all the features one at a time, minimizing the cost function with respect to the coordinate. Lasso regression does not have a closed form solution as it becomes a single variable problem so the solution is defined in terms of a soft threshold, S(), which is as follows:

```cpp
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
```

In the scenario where rho + lambda > 0, then the estimator remains unchanged. The coordinate descent update rule is defined as follows. For every j = 0, 1, ... n , we compute rho where 
![Rho equation](https://github.com/stelmacm/CSE-701/blob/main/rho%20equation.png?raw=true)

After this we set theta_j = S(rho_j , lambda). This allows us perform step wise coordinate descent for every column or arguement of X. What is interesting about this is that theta is a vector, so the multiplication operator had to be adjusted and created to allow for matrix and vector compatible multiplication. The final result is a list of theta's that are the $\beta$ estimators of the matrix X. 

### Leave One Out Cross Validation

To determine which of the lambda's best fit the model, we perform cross validation on all the estimators in order to determine which will result in the lowest mean squared error or in this case penalized mean squared error. We begin by creating a test set which can be easily done since our data is randomly generated. Once a training and test sets have been assigned, we compute the mean squared error of each Beta estimator and compare them. We identify the one that has the lowest mean squared error. Tibeshirani, a world famous data scientist who invented Lasso, states that the lowest mean squared error is infact an overfit of lasso regression and that the correct lambda associated with the mean squared error that is one standard deviation away from the original. In order to do this, we need to compute the standard deviation of the mean squared errors from the cross validation and find the one that is closest to the CV(theta) + SD(theta). This is done using a function created in order to find the closest value to that of the one in question 
```cpp

double_t closest(vector<double_t> const &vec, double_t value)
{
    auto const it = lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end())
    {
        return -1;
    }

    return *it;
}
```
This is in form of double rather than in the `template <typename T>` because the result should always return a double format and the MSE would be incorrect if returned as anything else. These values are stored in a vector and incremented with `.pushback()`. Because the same values will always be in the same place, it will become very easy for us to reference those values and call upon them. 
  
## Acknowledgements

Special thanks to Professor Barak Shoshany of McMaster University for guidance and aid throughout the project. 
