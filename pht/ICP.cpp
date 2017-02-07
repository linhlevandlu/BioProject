/*
 * ICP.cpp
 *
 *  Created on: Jan 31, 2017
 *      Author: linh
 */
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <float.h>
#include <string.h>
using namespace std;

#include "../utils/Converter.h"

#include "../imageModel/Point.h"
#include "../imageModel/Line.h"
#include "../imageModel/Matrix.h"
#include "../imageModel/Image.h"
#include "../segmentation/Thresholds.h"
#include "../segmentation/Canny.h"
#include "../utils/Drawing.h"
#include "GHTInPoint.h"
//#include "Vector.h"
#include "SVD.h"

#include "ICP.h"

/*Point centerPoint(vector<Point> listPoints)
 {
 int totalX = 0, totalY = 0;
 int size = listPoints.size();
 for (int i = 0; i < size; i++)
 {
 totalX += listPoints.at(i).getX();
 totalY += listPoints.at(i).getY();
 }
 return Point(totalX / size, totalY / size);
 }
 void icpmethod(vector<Point> modelPoints, vector<Point> scenePoints)
 {
 Point mCenter = centerPoint(modelPoints);
 Point sCenter = centerPoint(scenePoints);
 int dx = sCenter.getX() - mCenter.getX();
 int dy = sCenter.getY() - mCenter.getY();
 // translate scene to model
 for(size_t i = 0; i < modelPoints.size();i++)
 {
 int x = modelPoints.at(i).getX();
 int y = modelPoints.at(i).getY();
 modelPoints.at(i).setX(x+dx);
 modelPoints.at(i).setY(x+dy);
 }
 }*/
double errorSquared(vector<Line> pairPoints)
{
	double result = 0;
	int size = pairPoints.size();
	for (int i = 0; i < size; i++)
	{
		Point begin = pairPoints.at(i).getBegin();
		Point end = pairPoints.at(i).getEnd();
		result += (pow((double) (end.getX() - begin.getX()), 2)
			+ pow((double) (end.getY() - begin.getY()), 2));
	}
	return sqrt(result / size);
}
Matrix<float> equivalent(Point center, vector<Point> listPoints)
{
	Matrix<float> result(2, (int) listPoints.size(), 0);
	Point pi;
	float dx = 0.0, dy = 0.0;
	for (size_t i = 0; i < listPoints.size(); i++)
	{
		pi = listPoints.at(i);
		dx = pi.getX() - center.getX();
		dy = pi.getY() - center.getY();
		result.setAtPosition(0, (int) i, dx);
		result.setAtPosition(1, (int) i, dy);
	}
	return result;
}

void icpMethod()
{
	Matrix<float> m(2,2,0.0);
	Matrix<float> v(2,2,0.0);
	Matrix<float> s(2,2,0.0);
	Matrix<float> u(2,2,0.0);
	svdMethod(m,u,s,v);
}

/*void icpmethod(Image modelGray, Image &sceneGray)
{
	int rows = modelGray.getGrayMatrix()->getRows();
	int cols = modelGray.getGrayMatrix()->getCols();
	vector<Point> modelPoints;
	ptr_IntMatrix modelGrandient = new Matrix<int>(rows, cols, -1);
	*modelGrandient = *(getGradientDMatrix(modelGray, modelPoints));
	vector<Point> modelPoints2(modelPoints);
	Point mPoint;
	Line mAxis = principalAxis(modelGrandient, mPoint);
	Point mend = mAxis.getEnd();

	vector<Point> scenePoints;
	ptr_IntMatrix sceneGrandient = new Matrix<int>(rows, cols, -1);
	*sceneGrandient = *(getGradientDMatrix(sceneGray, scenePoints));
	Point sPoint;
	Line sAxis = principalAxis(sceneGrandient, sPoint);

	double angle = mAxis.angleLines(sAxis); // in degree
	int dx = sPoint.getX() - mPoint.getX();
	int dy = sPoint.getY() - mPoint.getY();
	// translate scene to model
	for (size_t i = 0; i < modelPoints.size(); i++)
	{
		int x = modelPoints.at(i).getX();
		int y = modelPoints.at(i).getY();
		modelPoints.at(i).setX(x + dx);
		modelPoints.at(i).setY(y + dy);
	}
	// translate the endpoint of mAxis
	mend.setX(mend.getX() + dx);
	mend.setY(mend.getY() + dy);
	int xnew = 0, ynew = 0;
	rotateAPoint(mend.getX(), mend.getY(), sPoint, angle, 1, xnew, ynew);
	double distance1 = sAxis.perpendicularDistance(Point(xnew, ynew));
	xnew = ynew = 0;
	rotateAPoint(mend.getX(), mend.getY(), sPoint, -angle, 1, xnew, ynew);
	double distance2 = sAxis.perpendicularDistance(Point(xnew, ynew));
	if (distance2 < distance1)
		angle = -angle;
	cout << "\nAngle difference: " << angle << endl;
	for (size_t i = 0; i < modelPoints.size(); i++)
	{
		int x = modelPoints.at(i).getX();
		int y = modelPoints.at(i).getY();
		int xnew = 0, ynew = 0;
		rotateAPoint(x, y, sPoint, angle, 1, xnew, ynew);
		modelPoints.at(i).setX(xnew);
		modelPoints.at(i).setY(ynew);
	}
	RGB color;
	color.R = 255;
	color.G = color.B = 0;
	for (size_t i = 0; i < scenePoints.size(); i++)
	{
		Point pi = scenePoints.at(i);
		if (pi.getX() >= 0 && pi.getX() < cols && pi.getY() >= 0
			&& pi.getY() < rows)
		{
			sceneGray.getRGBMatrix()->setAtPosition(pi.getY(), pi.getX(), color);
		}
	}
	color.R = 0;
	color.B = 255;
	for (size_t i = 0; i < modelPoints.size(); i++)
	{
		Point pi = modelPoints.at(i);
		if (pi.getX() >= 0 && pi.getX() < cols && pi.getY() >= 0
			&& pi.getY() < rows)
		{
			sceneGray.getRGBMatrix()->setAtPosition(pi.getY(), pi.getX(), color);
		}
	}
	color.R = 255;
	fillCircle(*sceneGray.getRGBMatrix(), sPoint, 5, color);
	// finding closest points
	vector<Line> pairPoints;
	double minerror = DBL_MAX;
	double t;
	//do
	//{
	for (size_t i = 0; i < modelPoints.size(); i++)
	{
		Point pm = modelPoints.at(i);
		double minDist = DBL_MAX;
		Point pk;
		for (size_t j = 0; j < scenePoints.size(); j++)
		{
			Point ps = scenePoints.at(j);
			Line temp(pm, ps);
			if (temp.getLength() < minDist)
			{
				minDist = temp.getLength();
				pk = ps;
			}
		}
		pairPoints.push_back(Line(pm, pk));
		//drawingLine(*sceneGray.getRGBMatrix(), Line(pm, pk), color);
	}
	vector<Point> modelCorr, sceneCorr;
	Line l;
	for (size_t i = 0; i < pairPoints.size(); i++)
	{
		l = pairPoints.at(i);
		modelCorr.push_back(l.getBegin());
		sceneCorr.push_back(l.getEnd());
	}

	// apply svd
	Matrix<float> modelP = equivalent(sPoint, modelCorr);
	cout << "\nModel P: " << modelP.getRows() << "\t" << modelP.getCols();
	Matrix<float> sceneQ = equivalent(sPoint, sceneCorr);
	cout << "\nScene Q: " << sceneQ.getRows() << "\t" << sceneQ.getCols() << endl;
	Matrix<float> m = modelP.multiply(sceneQ.transposition(0), 0.0);
	Matrix<float> mT = m.transposition(0.0);
	Matrix<float> mTm = mT.multiply(m, 0);

	float mTm00 = mTm.getAtPosition(0, 0);
	float mTm01 = mTm.getAtPosition(0, 1);
	float mTm10 = mTm.getAtPosition(1, 0);
	float mTm11 = mTm.getAtPosition(1, 1);
	// pt bac 2
	float a = 1;
	float b = -(mTm00 + mTm11);
	float c = (mTm00 * mTm11) - (mTm01 * mTm10);
	vector<float> n0 = ptB2(a, b, c);
	Matrix<float> v(2, 2);
	Matrix<float> w(2, 2, 0.0);
	if (n0.size() != 0)
	{
		for (int i = 0; i < n0.size(); ++i)
		{
			float k = n0.at(i);
			w.setAtPosition(i, i, sqrt(k));
			float k1 = mTm00 - k;
			float k2 = -mTm01 / k1;
			float summ = (k2 * k2) + 1;
			v.setAtPosition(0, i, k2 / summ);
			v.setAtPosition(1, i, 1 / summ);
		}
	}
	cout << "\nV matrix: " << endl;
	for (int r = 0; r < 2; ++r)
	{
		for (int c = 0; c < 2; ++c)
		{
			cout << "\t" << v.getAtPosition(r, c);
		}
		cout << "\n";
	}

	Matrix<float> cv = mTm.multiply(v, 0.0);
	float sumc1 = 0, sumc2 = 0;

	sumc1 = cv.getAtPosition(0, 0) + cv.getAtPosition(1, 0);
	sumc2 = cv.getAtPosition(0, 1) + cv.getAtPosition(1, 1);
	cout << "\nU matrix: " << endl;
	Matrix<float> u(2, 2, 0.0);
	for (int r = 0; r < 2; ++r)
	{
		for (int c = 0; c < 2; ++c)
		{
			float rc = 0.0;
			if (c % 2 == 0)
				rc = cv.getAtPosition(r, c) / sumc1;
			else
				rc = cv.getAtPosition(r, c) / sumc2;
			cout << "\t" << rc;
			u.setAtPosition(r, c, rc);
		}
		cout << "\n";
	}
	Matrix<float> r = v.multiply(u.transposition(0.0), 0.0);
	for (size_t i = 0; i < modelPoints2.size(); i++)
	{
		int x = modelPoints2.at(i).getX();
		int y = modelPoints2.at(i).getY();
		int xnew = 0, ynew = 0;
		xnew = mPoint.getX()
			+ (r.getAtPosition(0, 0) * (x - sPoint.getX())
				+ r.getAtPosition(0, 1) * (y - sPoint.getY()));
		ynew = mPoint.getY()
			+ (r.getAtPosition(1, 0) * (x - sPoint.getX())
				+ r.getAtPosition(1, 1) * (y - sPoint.getY()));
		modelPoints2.at(i).setX(xnew);
		modelPoints2.at(i).setY(ynew);
	}
	for (size_t i = 0; i < modelPoints2.size(); i++)
	{
		Point pi = modelPoints2.at(i);
		if (pi.getX() >= 0 && pi.getX() < cols && pi.getY() >= 0
			&& pi.getY() < rows)
		{
			sceneGray.getRGBMatrix()->setAtPosition(pi.getY(), pi.getX(), color);
		}
	}

	// compute stop criteria
	double err = errorSquared(pairPoints);
	 t = minerror - err;
	 minerror = err;

//}
//while (t > 1);
}*/
