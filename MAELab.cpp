/*
 * Analysis.cpp
 *
 *  Created on: Dec 7, 2016
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

#include "imageModel/Point.h"
#include "imageModel/Line.h"
#include "imageModel/Edge.h"
#include "imageModel/Matrix.h"
#include "imageModel/Image.h"

#include "histograms/ShapeHistogram.h"
#include "pht/PHTEntry.h"
#include "pht/PHoughTransform.h"
#include "correlation/CrossCorrelation.h"

#include "pointInterest/Treatments.h"
#include "pointInterest/Segmentation.h"
#include "pointInterest/GeometricHistgoram.h"
#include "pointInterest/ProHoughTransform.h"
#include "pointInterest/LandmarkDetection.h"

#include "MAELab.h"

vector<Line> segment(ptr_Treatments treatment, Image image, int minDistance)
{

	Segmentation sg;
	treatment = &sg;
	treatment->setRefImage(image);

	vector<Line> lines = sg.segment();
	return lines;
}

ShapeHistogram geomHistogram(ptr_Treatments treatment, Image image,
	AngleAccuracy angleAcc, int cols)
{
	GeometricHistgoram geoHistogram;
	treatment = &geoHistogram;
	treatment->setRefImage(image);
	return geoHistogram.geomHistogram(angleAcc, cols);
}
double bhattacharyyaDistance(ptr_Treatments treatment, Image sceneImage,
	AngleAccuracy angleAcc, int cols)
{
	GeometricHistgoram geoHistogram;
	treatment = &geoHistogram;
	return geoHistogram.bhattacharyyaDistance(sceneImage, angleAcc, cols);
}
PHoughTransform phtEntriesTable(ptr_Treatments treatment, Image image)
{
	ProHoughTransform pht;
	treatment = &pht;
	treatment->setRefImage(image);
	return pht.constructPHT();
}
vector<Point> estimatedLandmarks(ptr_Treatments treatment, Image sceneImage,
	AngleAccuracy acc, int cols, int templSize, int sceneSize, Point &ePoint,
	double angleDiff)
{
	LandmarkDetection lmd;
	lmd.setRefImage(treatment->getRefImage());
	//treatment = &lmd;
	return lmd.landmarksAutoDectect(sceneImage, acc, cols, templSize, sceneSize,
		ePoint, angleDiff);
}
double measureCentroidPoint(vector<Point> landmarks, Point &ebary)
{
	ebary.setX(0);
	ebary.setY(0);

	int totalX = 0;
	int totalY = 0;
	size_t lmSize = landmarks.size();

	Point pi;
	for (size_t i = 0; i < lmSize; i++)
	{
		pi = landmarks.at(i);
		totalX += pi.getX();
		totalY += pi.getY();
	}

	if (lmSize > 0)
	{
		ebary.setX(totalX / lmSize);
		ebary.setY(totalY / lmSize);
	}
	double totalDistance = 0;
	Point lm;
	for (size_t j = 0; j < lmSize; j++)
	{
		lm = landmarks.at(j);
		//ptr_Line line = new Line(ebary,lm);
		double distance = distancePoints(ebary,landmarks.at(j));
		//cout<<"\nDistance: \t"<<distance;
		totalDistance += (distance * distance);

	}
	//cout<<"\n"<<sqrt(totalDistance);
	return sqrt(totalDistance);
}

