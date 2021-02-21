# Comparing optimal penalization parameter lambda for Lasso Regression found using Coordinate Descent - Martin Stelmach - February 21, 2021

## Summary of Project

The goal of the program created is to take user find the optimal value of lambda for which we do penalized regression. This is done by taking a user defined randomly generated dataset and doing coordinate descent along log-spaced lambda values until the penalization has been performed. From there the optimal lambda is selected via cross validation and the penalized mean squared error is compared to the standard OLS proving that the lasso is a better fit but avoids overfitting, as we will further see. It is worth mentioning that lasso optimization is not a closed form solution so a forward stagewise algorithm is implemented. 

## Explanation of Algorithm

Regression is one of the most important tools in statistics. Often times, we fit data using a simple least squares approximation, however, this is not always sufficient. As data becomes larger and we have more predictors, we start to encounter more variance in our estimates for $\beta$. To avoid this, we apply a penalization parameter which regularizes the coefficients and reduces variance. By implementing L1 penalization, we are ensuring our coefficients do not blow up and get really big. Once we have done the coordinate descent with the lambda in question, we use the one standard rule for model selection. This means that the optimal lambda is selected one standard deviation from the optimal that is found via cross validation. This is because the optimal lambda in cross validation can often times be an overfit and result in certain estimators being very close to zero but not zero, and shifting by a standard deviation would "push" those over.

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

### Ordinary Least Squares

After the matrices are created 
  
## Acknowledgements

Special thanks to Professor Barak Shoshany of McMaster University for guidance and aid throughout the project. 
