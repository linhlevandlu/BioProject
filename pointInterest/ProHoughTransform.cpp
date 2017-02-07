/*
 * ProHoughTransform.cpp
 *
 *  Created on: Dec 2, 2016
 *      Author: linh
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#include "../imageModel/Point.h"
#include "../imageModel/Line.h"
#include "../imageModel/Edge.h"
#include "../imageModel/Matrix.h"
#include "../imageModel/Image.h"

#include "../segmentation/Canny.h"
#include "../segmentation/Thresholds.h"

#include "../pht/PHTEntry.h"
#include "../pht/PHoughTransform.h"
#include "../pht/GHTInPoint.h"

#include "Treatments.h"
#include "ProHoughTransform.h"

vector<Point> findLandmarks(Point refPoint, Point esPoint,
		vector<Point> refLandmarks, int width, int height, int &positive)
{
	vector<Point> esLandmarks;
	positive = 0;
	Point temp, lm;
	for (size_t t = 0; t < refLandmarks.size(); t++)
	{
		temp = refLandmarks.at(t);

		lm.setX(temp.getX());
		lm.setY(temp.getY());

		int px = refPoint.getX() - esPoint.getX();
		int py = refPoint.getY() - esPoint.getY();
		int x;
		int y;
		x = lm.getX() - px;
		if (py > 0)
			y = lm.getY() - py;
		else
			y = lm.getY() + py;
		esLandmarks.push_back(Point(x, y));
		if (x >= 0 && x <= width && y >= 0 && y <= height)
		{
			positive += 1;
		}
	}
	return esLandmarks;
}

Point refPointInScene(PHTEntry entry, vector<Line> matchLines,
		double &angleDiff, vector<Point> refLandmarks, int width, int height)
{

	Point inter(0, 0);
	Point refPoint(width / 2, height / 2);

	if (matchLines.size() <= 0)
		return inter;
	Line objl1 = matchLines.at(0);
	Line objl2 = matchLines.at(1);

	Line lineEntry1 = entry.getRefLine();
	Line lineEntry2 = entry.getObjLine();

	HoughSpace hs1, hs2;
	hs1 = entry.getListHoughSpace().at(0);
	hs2 = entry.getListHoughSpace().at(1);

	// display the angle
	double angle3 = lineEntry1.angleLines(objl1);
	double angle4 = lineEntry1.angleLines(objl2);

	vector<Point> intersects1 = objl1.interParallel(objl1, objl2, hs1.distance,
			hs2.distance, width, height);
	int max = 0;
	vector<Point> estLM;
	vector<double> angles;
	Point esPoint;
	vector<Point> lms;
	int positive = 0;
	for (size_t i = 0; i < intersects1.size(); i++)
	{
		esPoint = intersects1.at(i);
		positive = 0;
		lms = findLandmarks(refPoint, esPoint, refLandmarks, width, height,
				positive);
		if (positive > max)
		{
			estLM.clear();
			estLM.push_back(esPoint);
			angles.clear();
			angles.push_back(angle3);
			max = positive;
		}
		else
		{
			if (positive == max)
			{
				estLM.push_back(esPoint);
				angles.push_back(angle3);
			}
		}
	}

	vector<Point> intersects2 = objl1.interParallel(objl1, objl2, hs2.distance,
			hs1.distance, width, height);

	for (size_t i = 0; i < intersects2.size(); i++)
	{
		esPoint = intersects2.at(i);
		positive = 0;
		lms = findLandmarks(refPoint, esPoint, refLandmarks, width, height,
				positive);
		if (positive > max)
		{
			estLM.clear();
			estLM.push_back(esPoint);
			angles.clear();
			angles.push_back(angle4);
			max = positive;
		}
		else
		{
			if (positive == max)
			{
				estLM.push_back(esPoint);
				angles.push_back(angle4);
			}
		}
	}

	intersects1.clear();
	intersects2.clear();

	if (estLM.size() == 1)
	{
		inter = estLM.at(0);
		angleDiff = angles.at(0);
	}
	else
	{
		double angleds = 180;
		Point es;
		for (size_t i = 0; i < estLM.size(); i++)
		{
			es = estLM.at(i);
			double angle = angles.at(i);
			if (angle <= angleds)
			{
				angleds = angle;
				angleDiff = angleds;
				inter = es;
			}
		}
	}
	cout << "\n Reference point in scene: " << inter.getX() << ", "
			<< inter.getY();

	return inter;
}
bool similarPairLines(Line ref1, Line ref2, Line scene1, Line scene2)
{
	int cond1 = 1;
	int cond2 = 1;
	int cond3 = 2;

	double refAngle = ref1.angleLines(ref2);
	double rd1 = ref1.perpendicularDistance(ref2.getBegin());
	double rd2 = ref1.perpendicularDistance(ref2.getEnd());
	double rd = rd1 + rd2;

	double sceneAngle = scene1.angleLines(scene2);
	double sd1 = scene1.perpendicularDistance(scene2.getBegin());
	double sd2 = scene1.perpendicularDistance(scene2.getEnd());
	double sd = sd1 + sd2;

	if (abs(refAngle - sceneAngle) < cond1
			&& (abs(
					(ref1.getLength() / scene1.getLength())
							- (ref2.getLength() / scene2.getLength())) < cond2)
			&& (abs(rd - sd) < cond3))
	{
		return true;
	}
	return false;
}

PHTEntry findHoughSpace(vector<PHTEntry> entryTable, Line line1, Line line2)
{
	PHTEntry entry;
	Line ref1, ref2;
	for (size_t i = 0; i < entryTable.size(); i++)
	{
		ref1 = entryTable.at(i).getRefLine();
		ref2 = entryTable.at(i).getObjLine();
		if (similarPairLines(ref1, ref2, line1, line2))
		{
			entry = entryTable.at(i);
		}
	}
	return entry;
}

PHTEntry matchingInScene(vector<PHTEntry> entryTable, vector<Line> sceneLines,
		int width, int height, vector<Line> &maxVector)
{
	ptr_IntMatrix accumulator = new Matrix<int>(
			floor(sqrt(width * width + height * height)), 361);
	int maxValue = 0;
	PHTEntry maxEntry;
	Line objLine1;
	Line objLine2;
	PHTEntry entry;
	vector<HoughSpace> hspace;
	for (size_t i = 0; i < sceneLines.size(); i++)
	{
		objLine1 = sceneLines.at(i);
		for (size_t j = 0; j < sceneLines.size(); j++)
		{
			objLine2 = sceneLines.at(j);
			if (i != j && closetLine(objLine2, objLine1))
			{
				entry = findHoughSpace(entryTable, objLine1, objLine2);
				hspace = entry.getListHoughSpace();
				if (hspace.size() != 0)
				{
					HoughSpace hsp;
					for (size_t k = 0; k < hspace.size(); k++)
					{
						hsp = hspace.at(k);
						int angle = round(hsp.angle);
						int distance = round(hsp.distance);
						if (!isnan(angle) && !isnan(distance) && angle >= 0
								&& distance >= 0)
						{
							int value = accumulator->getAtPosition(distance,
									angle);
							accumulator->setAtPosition(distance, angle,
									value + 1);
							if (accumulator->getAtPosition(distance, angle)
									> maxValue)
							{
								maxVector.clear();
								maxVector.push_back(objLine1);
								maxVector.push_back(objLine2);
								maxValue = accumulator->getAtPosition(distance,
										angle);
								maxEntry.setRefLine(entry.getRefLine());
								maxEntry.setObjLine(entry.getObjLine());
								maxEntry.setListHoughSpace(
										entry.getListHoughSpace());
								//maxEntry= entry;
							}
							else
							{
								if (k == 0
										&& accumulator->getAtPosition(distance,
												angle) == maxValue)
								{
									maxVector.push_back(objLine1);
									maxVector.push_back(objLine2);
								}
							}
						}
					}
				}
			}
		}
	}

	delete accumulator;
	return maxEntry;
}
vector<Point> phtLandmarks(vector<PHTEntry> entriesTable, Point refPoint,
		vector<Line> sceneLines, int width, int height,
		vector<Point> mLandmarks, double &angleDiff, Point &ePoint)
{
	vector<Point> eLandmarks;
	vector<Line> maxVector;
	PHTEntry entry = matchingInScene(entriesTable, sceneLines, width, height,
			maxVector);
	if (maxVector.size() > 0)
	{

		ePoint = refPointInScene(entry, maxVector, angleDiff, mLandmarks, width,
				height);
		double angle1 = entry.getRefLine().angleLines(entry.getObjLine());
		double angle2 = maxVector.at(0).angleLines(maxVector.at(1));
		angleDiff += abs(angle1 - angle2);
		int positive = 0;
		eLandmarks = findLandmarks(refPoint, ePoint, mLandmarks, width, height,
				positive);
	}
	maxVector.clear();

	return eLandmarks;
}
ProHoughTransform::ProHoughTransform()
{
	// TODO Auto-generated constructor stub

}

ProHoughTransform::~ProHoughTransform()
{
	// TODO Auto-generated destructor stub
}
PHoughTransform ProHoughTransform::constructPHT()
{
	ptr_IntMatrix grayImage;
	grayImage = (Matrix<int> *) malloc(sizeof(Matrix<int> ));
	*grayImage = *(Treatments::refImage.getGrayMatrix());
	//grayImage = Treatments::refImage.getGrayMatrix();
	int width = grayImage->getCols();
	int height = grayImage->getRows();

	PHoughTransform pht;
	pht.setRefPoint(Point(width / 2, height / 2));
	pht.constructPHTTable(Treatments::refImage.getListOfLines());

	delete grayImage;
	return pht;
}

vector<Point> ProHoughTransform::estimateLandmarks(Image sImage,
		double &angleDiff, Point &ePoint)
{
	vector<Point> eLandmarks;
	ptr_IntMatrix mMatrix = Treatments::refImage.getGrayMatrix();
	int width = mMatrix->getCols();
	int height = mMatrix->getRows();

	vector<Point> mLandmarks = Treatments::refImage.getListOfManualLandmarks();
	vector<Line> mLines = Treatments::refImage.getListOfLines();
	vector<Line> sLines = sImage.getListOfLines();

	Point mPoint(width / 2, height / 2);

	PHoughTransform pht;
	pht.setRefPoint(mPoint);
	vector<PHTEntry> entryTable = pht.constructPHTTable(mLines);

	eLandmarks = phtLandmarks(entryTable, mPoint, sLines, width, height,
			mLandmarks, angleDiff, ePoint);

	return eLandmarks;

}

vector<Point> ProHoughTransform::generalTransform(Image &sImage, double &angle,
		Point &ePoint, Point &mPoint, ptr_IntMatrix &newScene)
{
	int rows = sImage.getGrayMatrix()->getRows();
	int cols = sImage.getGrayMatrix()->getCols();
	Image mImage = Treatments::refImage;
	ptr_IntMatrix mgradirection = new Matrix<int>(rows, cols, -1);
	vector<Point> modelPoints;
	*mgradirection = *(getGradientDMatrix(mImage,modelPoints));
	vector<Point> mLandmarks = mImage.getListOfManualLandmarks();

	ptr_IntMatrix gradirection = new Matrix<int>(rows, cols, -1);
	vector<Point> scenePoints;
	*gradirection = *(getGradientDMatrix(sImage,scenePoints));

	vector<Point> eslm;
	Point translation;
	eslm = generalizingHoughTransform(mgradirection, gradirection, mLandmarks,
			ePoint, mPoint, angle, translation);
	int dx = translation.getX();
	int dy = translation.getY();
	// move the scene to the same and rotate the model
	for (int r = 0; r < rows; r++)
	{
		for (int c = 0; c < cols; c++)
		{
			int value = sImage.getGrayMatrix()->getAtPosition(r, c);
			int xnew = c + dx;
			int ynew = r + dy;
			rotateAPoint(c + dx, r + dy, mPoint, angle, 1, xnew, ynew);
			if (xnew >= 0 && ynew >= 0 && ynew < rows && xnew < cols)
			{
				newScene->setAtPosition(ynew, xnew, value);
			}
		}
	}
	return eslm;
}
