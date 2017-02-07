/*
 * GHTInPoint.cpp
 *
 *  Created on: Jan 17, 2017
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
#include "GHTInPoint.h"

ptr_IntMatrix getGradientDMatrix(Image grayImage, vector<Point> &edgePoints)
{
	int rows = grayImage.getGrayMatrix()->getRows();
	int cols = grayImage.getGrayMatrix()->getCols();
	ptr_IntMatrix mMatrix = new Matrix<int>(rows, cols, 0);
	*mMatrix = *(grayImage.getGrayMatrix());
	int mThresholdValue = (int) grayImage.getThresholdValue();
	ptr_IntMatrix mbinMatrix = new Matrix<int>(
		*binaryThreshold(mMatrix, mThresholdValue, 255));
	ptr_IntMatrix mgradirection = new Matrix<int>(rows, cols, -1);
	ptr_IntMatrix mcannyMatrix = new Matrix<int>(
		*cannyProcess2(mbinMatrix, mThresholdValue, 3 * mThresholdValue,
			mgradirection, edgePoints));
	delete mMatrix;
	delete mbinMatrix;
	delete mcannyMatrix;
	return mgradirection;
}
double angleVector(Point p1, Point p2, Point q1, Point q2)
{
	Point vp(p2.getX() - p1.getX(), p2.getY() - p1.getY());
	Point vq(q2.getX() - q1.getX(), q2.getY() - q1.getY());

	double lengthVP = sqrt(
		(double) (vp.getX() * vp.getX()) + (double) (vp.getY() * vp.getY()));
	double lengthVQ = sqrt(
		(double) (vq.getX() * vq.getX()) + (double) (vq.getY() * vq.getY()));

	double vpq = (double) (vp.getX() * vq.getX())
		+ (double) (vp.getY() * vq.getY());

	return acos(vpq / (lengthVP * lengthVQ));
}
RTable rTableConstruct(ptr_IntMatrix gradMatrix, Point center)
{
	RTable rtable;
	rtable.center.setX(center.getX());
	rtable.center.setY(center.getY());
	rtable.entriesTable.resize(361);

	// in default, rtable includes 360 rows for 360 degrees of gradient directions.
	int rows = gradMatrix->getRows();
	int cols = gradMatrix->getCols();

	int gd;
	REntry rentry;
	//ofstream off("rtable.txt");
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			gd = gradMatrix->getAtPosition(r, c);
			if (gd != -1)
			{
				PolarValue pv;
				pv.p.setX(c);
				pv.p.setY(r);
				Point hPoint(pv.p.getX() + 100, pv.p.getY());
				pv.angle = angleVector(pv.p, hPoint, pv.p, rtable.center);
				pv.distance = (Line(rtable.center, pv.p)).getLength();
				rentry = rtable.entriesTable.at((gd));

				(rtable.entriesTable.at((gd))).gradient = gd;
				(rtable.entriesTable.at((gd))).polarValues.push_back(pv);
				//off << "\n" << pv.angle << "\t" << pv.distance;
			}
		}
	}
	//off.close();
	return rtable;
}

Point houghSpace(ptr_IntMatrix gradMatrix, RTable rentries)
{
	int rows = gradMatrix->getRows();
	int cols = gradMatrix->getCols();
	ptr_IntMatrix acc = new Matrix<int>(rows, cols, 0);

	int gradient = 0;
	REntry entry;
	int maxValue = 0, maxXIndex = 0, maxYIndex = 0;
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			gradient = gradMatrix->getAtPosition(r, c);
			if (gradient != -1)
			{
				entry = rentries.entriesTable.at(gradient);
				if (entry.gradient == gradient)
				{
					vector<PolarValue> polarValues = entry.polarValues;
					if (polarValues.size() > 0)
					{
						//PolarValue pvalue;
						double avalue = 0, dvalue = 0;
						int xindex = 0, yindex = 0, accValue = 0;
						for (size_t k = 0; k < polarValues.size(); k++)
						{
							PolarValue pvalue = polarValues.at(k);
							avalue = pvalue.angle;
							dvalue = pvalue.distance;

							xindex = c + (dvalue * cos(avalue));
							yindex = r + (dvalue * sin(avalue));
							//cout << "\n" << xindex << "\t" << yindex;
							if (xindex >= 0 && yindex >= 0 && xindex < cols && yindex < rows)
							{
								// increase the value in accumulator
								accValue = acc->getAtPosition(yindex, xindex);
								accValue += 1;
								acc->setAtPosition(yindex, xindex, accValue);
							}
							if (accValue > maxValue && xindex > 0 && xindex < cols
								&& yindex > 0 && yindex < rows)
							{
								maxValue = accValue;
								maxXIndex = xindex;
								maxYIndex = yindex;
							}
						}
					}
				}
			}
		}
	}
	cout << "\n center: ";
	rentries.center.toString();
	cout << "\nMax Value:" << maxValue;
	cout << "\nMax X index:" << maxXIndex;
	cout << "\nMax Y index:" << maxYIndex << endl;
	return Point(maxXIndex, maxYIndex);
}

double avgDistance(vector<Point> listPoints, Line axis)
{
	double totalDist = 0;
	size_t nPoints = listPoints.size();
	for (size_t j = 0; j < nPoints; j++)
	{
		Point pj = listPoints.at(j);
		double distance = axis.perpendicularDistance(pj);
		totalDist += distance;
	}
	return totalDist / (int) nPoints;
}

// compute centroid point from the gradient matrix of image
// output: centroid point and list points of edge.
Point centroidPoint(ptr_IntMatrix gradMatrix, vector<Point> &listPoints)
{
	int rows = gradMatrix->getRows();
	int cols = gradMatrix->getCols();
	int totalX = 0, totalY = 0, count = 0;
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			if (gradMatrix->getAtPosition(r, c) != -1)
			{
				listPoints.push_back(Point(c, r));
				totalX += c;
				totalY += r;
				count++;
			}
		}
	}
	return Point(totalX / count, totalY / count);
}

Line principalAxis(ptr_IntMatrix gradMatrix, Point &cPoint)
{
	vector<Point> listOfPoints;
	cPoint = centroidPoint(gradMatrix, listOfPoints);
	Point sPoint;
	double minAvgDist = DBL_MAX;
	size_t nPoints = listOfPoints.size();
	for (size_t i = 0; i < nPoints; ++i)
	{
		Point pi = listOfPoints.at(i);
		Line l(cPoint, pi);
		double avgDist = avgDistance(listOfPoints, l);
		if (avgDist < minAvgDist)
		{
			minAvgDist = avgDist;
			sPoint.setX(pi.getX());
			sPoint.setY(pi.getY());
		}
	}

	return Line(cPoint, sPoint);
}
vector<Point> detectLandmarks(Point refPoint, Point ePoint,
	vector<Point> mlandmarks)
{
	vector<Point> esLandmarks;
	Point mlm;
	int deltaX = 0, deltaY = 0, xn = 0, yn = 0;
	for (size_t i = 0; i < mlandmarks.size(); i++)
	{
		mlm = mlandmarks.at(i);
		deltaX = ePoint.getX() - refPoint.getX();
		deltaY = ePoint.getY() - refPoint.getY();

		xn = mlm.getX() + deltaX;
		yn = mlm.getY() + deltaY;
		esLandmarks.push_back(Point(xn, yn));
	}
	return esLandmarks;
}

vector<Point> generalizingHoughTransform(ptr_IntMatrix mGradient,
	ptr_IntMatrix sGradient, vector<Point> mLandmarks, Point &ePoint,
	Point &mPoint, double &angle, Point &translation)
{
	int rows = mGradient->getRows();
	int cols = mGradient->getCols();
	Point center(cols / 2, rows / 2);
	RTable rentries = rTableConstruct(mGradient, center);
	Point sPoint = houghSpace(sGradient, rentries);
	vector<Point> eslm = detectLandmarks(center, sPoint, mLandmarks);

	Point mTemp;
	Line sLine = principalAxis(sGradient, ePoint);
	Line mLine = principalAxis(mGradient, mTemp);
	double anglek = sLine.angleLines(mLine); // compute angle difference between two images
	int dxold = sPoint.getX() - center.getX();
	int dyold = sPoint.getY() - center.getY();
	mTemp.setX(mTemp.getX() + dxold);
	mTemp.setY(mTemp.getY() + dyold);

	// move sPoint to mTemp
	int dx = mTemp.getX() - ePoint.getX();
	int dy = mTemp.getY() - ePoint.getY();
	cout << "\nTranslate: " << dx << "\t" << dy << endl;

	// checking the rotate direction
	sLine.getBegin().setX(mTemp.getX());
	sLine.getBegin().setY(mTemp.getY());
	sLine.getEnd().setX(sLine.getEnd().getX() + dx);
	sLine.getEnd().setY(sLine.getEnd().getY() + dy);
	int xpk = 0, ypk = 0;
	Point pk = sLine.getEnd();
	rotateAPoint(pk.getX(), pk.getY(), mTemp, anglek, 1, xpk, ypk);
	Point pk1(xpk, ypk);
	xpk = ypk = 0;
	rotateAPoint(pk.getX(), pk.getY(), mTemp, -anglek, 1, xpk, ypk);
	Point pk2(xpk, ypk);
	double angled1 = mLine.angleLines(Line(mTemp, pk1));
	double angled2 = mLine.angleLines(Line(mTemp, pk2));
	if (angled2 < angled1 && anglek <= 90)
		anglek = -anglek;
	if (angled2 < angled1 && anglek >= 90)
		anglek = -anglek;
	angle = anglek;
	if (isnan(angle))
		angle = 0;
	cout << "\nAngle: " << angle << endl;
	cout << "\nLandmarks by GHT: " << eslm.size() << endl;

	mPoint.setX(mTemp.getX());
	mPoint.setY(mTemp.getY());
	translation.setX(dx);
	translation.setY(dy);
	return eslm;
}
