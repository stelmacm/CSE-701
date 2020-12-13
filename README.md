# Finding optimal penalization parameter lambda for ridge regression - Martin Stelmach - December 13, 2020

## Summary of Project

The goal of the program created is to take user input regarding the size of the data, check if the data is multicollinear by checking the variance inflation factor, and if it is, using gradient descent to find the optimal value of lambda that minimizes the mean squared error and hence finds the optimal penalization parameter.

## Explanation of methods

Regression is used for a multitude of things and is one of the most important things in todays world. For this reason there are so many different things one can do with regression to try and fit data better. One of those things is ridge regression. Ridge regression is used for data has lots of multicollinearity and there for is harder to predict. Typical dataset that ridge regression is used for is datasets that contain more factors than observations, making it difficult to interpret which are more important and should have higher coefficients. The use of ridge regression on such a multicollinear dataset is effective because ridge regression implements bias into the data. What happens is a shrinkage estimator is assigned to each coefficent to differentiate the more important coefficients. This is essentially a form of penalization to factors that should have less of a weight to the regression. The penalization parameter gives us a new equation for regression coefficients in the form of $$\hat{\beta}$$
