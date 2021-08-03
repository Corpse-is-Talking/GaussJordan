#include<vector>
#include<iostream>
#include<iomanip>
#include"Matrix.h"
/*
	Things to Consider more..
	1. swapping rows when leading 0's occur.. to avoid dividing by 0.. -? row swapping okay....
	2. Getting Rank of the Matrix-> Check if the matrix can be inverted.. okay..
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