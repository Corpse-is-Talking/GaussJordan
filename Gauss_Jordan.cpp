#include<vector>
#include<iostream>
#include<iomanip>
#include"Matrix.h"
/*
Description:
Using basic algorithm that i learned in Linear algebra class, (Gaussian Jordan), I made a cpp that gets the 
inverse of specific matrix , even if the matrix is not invertible, it will get the rank of the matrix and tell
you if the matrix is invertible or not.
As this was the first thing that i've made with cpp, 
Im sorry that this program can have some errors and codes are not easy to see..
*/
using namespace std;
int main()
{
	std::vector < std::vector < double>> vec;
	vec = { {1,1,1,-1} ,{1,1,-1,1},{1,-1,1,1},{-1,1,1,1} };
	Matrix mat(vec, vec.capacity());
	mat.printMatrix(); 
	mat.printInverse();
	std::vector<std::vector<double>> vec2;
	vec2 = { {0,0,0},{0,1,0},{0,0,1} };
	Matrix mat2(vec2, vec2.capacity());
	mat2.printMatrix();
	mat2.printInverse();
	cout << "Rank is "<<mat2.getRank() << endl;
	cout << endl;
	cout << endl;
	
}