/*
 * Matrix.cpp
 *
 *  Created on: Sep 21, 2016
 *      Author: linh
 */

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <time.h>

using namespace std;

#include "Point.h"

#include "Matrix.h"

template<typename T>
Matrix<T>::Matrix()
{
	rows = 0;
	cols = 0;
}
template<typename T>
Matrix<T>::Matrix(int nrows, int ncols)
{
	rows = nrows;
	cols = ncols;
	// Initialization the matrix
	Init();
}
template<typename T>
Matrix<T>::Matrix(int nrows, int ncols, T value)
{
	rows = nrows;
	cols = ncols;
	InitWithValue(value);
}
template<typename T>
Matrix<T>::Matrix(const Matrix<T> &tMatrix)
{
	this->rows = tMatrix.rows;
	this->cols = tMatrix.cols;
	Init();
	this->data = tMatrix.data;
}
template<typename T>
void Matrix<T>::Init()
{
	data.resize(getRows());
	for (int i = 0; i < getRows(); ++i)
		data[i].resize(getCols());
}

template<typename T>
void Matrix<T>::InitWithValue(T value)
{
	data.resize(getRows());
	for (int i = 0; i < getRows(); ++i)
		data[i].resize(getCols(), value);
}

template<typename T>
T Matrix<T>::getAtPosition(int r, int c)
{
	return data[r][c];
}

template<typename T>
void Matrix<T>::setAtPosition(int r, int c, T value)
{
	data[r][c] = value;
}

template<typename T>
void Matrix<T>::printMatrix()
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			//cout << data[i][j] << "\t";
		}
		cout << endl;
	}
}
template<typename T>
void Matrix<T>::operator=(const Matrix<T> &tMatrix)
{
	rows = tMatrix.rows;
	cols = tMatrix.cols;
	//Init();
	data = tMatrix.data;
}
// get rotation matrix with angle (in degree) and scale
ptr_DoubleMatrix getRotationMatrix2D(Point center, double angle, double scale)
{
	double alpha = cos(angle * M_PI / 180) * scale;
	double beta = sin(angle * M_PI / 180) * scale;

	ptr_DoubleMatrix rotateM = new Matrix<double>(2, 3, 0);

	rotateM->setAtPosition(0, 0, alpha);
	rotateM->setAtPosition(0, 1, beta);
	rotateM->setAtPosition(0, 2,
		(1 - alpha) * center.getX() - beta * center.getY());
	rotateM->setAtPosition(1, 0, -beta);
	rotateM->setAtPosition(1, 1, alpha);
	rotateM->setAtPosition(1, 2,
		beta * center.getX() + (1 - alpha) * center.getY());

	return rotateM;
}
// angle in degree
void rotateAPoint(int x, int y, Point center, double angle, double scale,
	int &xnew, int &ynew)
{
	xnew = 0;
	ynew = 0;
	ptr_DoubleMatrix rotationMatrix = getRotationMatrix2D(center, angle, scale);

	double a00 = rotationMatrix->getAtPosition(0, 0);
	double a01 = rotationMatrix->getAtPosition(0, 1);
	double b00 = rotationMatrix->getAtPosition(0, 2);
	double a10 = rotationMatrix->getAtPosition(1, 0);
	double a11 = rotationMatrix->getAtPosition(1, 1);
	double b10 = rotationMatrix->getAtPosition(1, 2);
	xnew = round(a00 * x + a01 * y + b00);
	ynew = round(a10 * x + a11 * y + b10);
	delete rotationMatrix;
}
template<typename T>
void Matrix<T>::rotation(Point center, double angle, double scale,
	T defaultValue)
{
	Matrix<T> result(rows, cols);
	result.InitWithValue(defaultValue);
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < cols; col++)
		{
			T value = data[row][col];

			int xnew = 0; //round(a00 * col + a01 * row + b00);
			int ynew = 0; //round(a10 * col + a11 * row + b10);
			rotateAPoint(col, row, center, angle, scale, xnew, ynew);
			if (xnew >= 0 && xnew < cols && ynew >= 0 && ynew < rows)
			{
				result.data[ynew][xnew] = value;
			}
		}
	}
	for (int row = 0; row < rows; row++)
	{
		for (int col = 1; col < cols - 1; col++)
		{
			T value = result.getAtPosition(row, col);
			T leftValue = result.getAtPosition(row, col - 1);
			T rightValue = result.getAtPosition(row, col + 1);
			if (value == 0 && leftValue > 0 && rightValue >= 0)
			{
				result.data[row][col] = (leftValue + rightValue) / 2;
			}
		}
	}
	data = result.data;
	//return result;
}
template<typename T>
Matrix<T> Matrix<T>::translate(int dx, int dy, T defaultValue)
{
	Matrix<T> result(rows, cols);
	result.InitWithValue(defaultValue);
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < cols; col++)
		{
			T value = data[row][col];

			int xnew = col + dx;
			int ynew = row + dy;
			if (xnew >= 0 && xnew < cols && ynew >= 0 && ynew < rows)
			{
				result.data[ynew][xnew] = value;
			}
		}
	}
	return result;
}
template<typename T>
Matrix<T> Matrix<T>::transposition(T defaultValue)
{
	Matrix<T> result(cols, rows);
	result.InitWithValue(defaultValue);
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < cols; col++)
		{
			T value = data[row][col];
			result.data[col][row] = value;
		}
	}
	return result;
}
template<typename T>
Matrix<T> Matrix<T>::multiply(Matrix<T> object, T defaultValue)
{
	int objRows = object.rows;
	int objCols = object.cols;
	Matrix<T> result(rows, objCols);
	result.InitWithValue(defaultValue);
	if (cols != objRows)
		return result;
	for (int row = 0; row < rows; row++)
	{
		for (int col = 0; col < objCols; col++)
		{
			T sum = defaultValue;
			for (int k = 0; k < cols; k++)
			{
				T value1 = data[row][k];
				T value2 = object.data[k][col];
				sum += value1 * value2;
			}
			result.data[row][col] = sum;
		}
	}
}

template class Matrix<int> ;
template class Matrix<double> ;
template class Matrix<float>;
template class Matrix<RGB> ;

