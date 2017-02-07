/*
 * Thresholds.cpp
 *
 *  Created on: Nov 26, 2016
 *      Author: linh
 */
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <vector>
#include <string.h>
#include <fstream>
#include <time.h>

using namespace std;
#include "../imageModel/Point.h"
#include "../imageModel/Line.h"
#include "../imageModel/Edge.h"
#include "../imageModel/Matrix.h"

#include "Thresholds.h"


ptr_IntMatrix binaryThreshold(ptr_IntMatrix inputMatrix, int tValue,
		int maxValue)
{
	int rows = inputMatrix->getRows();
	int cols = inputMatrix->getCols();

	ptr_IntMatrix result = new Matrix<int>(rows, cols, maxValue);
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			if (inputMatrix->getAtPosition(r, c) > tValue)
				result->setAtPosition(r, c, maxValue);
			else
				result->setAtPosition(r, c, 0);
		}
	}

	return result;
}
