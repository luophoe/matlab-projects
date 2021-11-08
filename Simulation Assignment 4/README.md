# Simulation Assignment 4
Develop your Derivative / Integration code that can do the following:  
a) Upload the data which will be in the form of x, f(x) (one file will be uploaded, test_1).  
b) Once you run the program, it will ask you what you want to perform Derivative or Integration (you select one of them)  
c) If you selected Derivative, then the program will do the following:  
i. Ask you to decide at what point, p you want to perform the derivative  
ii. If the point is from the data set and the spacing between the points is even, then the program will calculate the derivative using the CDD method.  
iii. If the point is not from the data set and/or the spacing between the points is not even, then the program first will use polynomial regression to estimate the function and then perform Derivative using CDD method with h = minimum of Δx from the data points.  
iv. show the solution for following cases of test1.txt  
1) p = 7  
2) p = 18.5  
d) If you selected Integration, then the program will do the following:  
i. Ask for the integration limit; p1, p2 where p2 > p1  
ii. Ask for the number of segments; n  
iii. Generate new data set with the limits p1, p2 and n evenly spaced segments. If all the values of the new data set belong to the original data set, then the program will calculate the Integration simply using the discrete Trapezoidal method.  
iv. If the limits (p1 and p2) are within the range of the original data set, but any of the values of the new data set do not belong to the original dataset, then the program first will use polynomial regression to estimate the function and then perform Integration using the Trapezoidal method. Select appropriate Dx.  
v. If at least one of the limits (p1 and p2) is out of the range of the original data set, display a message specifying the error (without calculating the integral).  
vi. show the solution for following cases of test1.txt  
1) p1 = 3, p2 = 7, n = 4  
2) p1 = 3, p2 = 7, n = 10  
