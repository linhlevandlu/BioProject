/*
 * LandmarkDetection.cpp
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

#include "../imageModel/Point.h"
#include "../imageModel/Line.h"
#include "../imageModel/Edge.h"
#include "../imageModel/Matrix.h"
#include "../imageModel/Image.h"

#include "../segmentation/Canny.h"
#include "../segmentation/Thresholds.h"

#include "../histograms/ShapeHistogram.h"
#include "../pht/PHTEntry.h"
#include "../pht/PHoughTransform.h"
#include "../pht/GHTInPoint.h"
#include "../correlation/CrossCorrelation.h"

#include "Treatments.h"
#include "Segmentation.h"
#include "GeometricHistgoram.h"
#include "ProHoughTransform.h"
#include "LandmarkDetection.h"

LandmarkDetection::LandmarkDetection()
{
	// TODO Auto-generated constructor stub

}

LandmarkDetection::~LandmarkDetection()
{
	// TODO Auto-generated destructor stub
}

vector<Point> LandmarkDetection::landmarksAutoDectect(Image sceneImage,
	AngleAccuracy acc, int cols, int templSize, int sceneSize, Point &ePoint,
	double &angleDiff)
{
	vector<Point> result;
	Image modelImage = Treatments::refImage;
	vector<Point> manualLMs = modelImage.getListOfManualLandmarks();
	vector<Line> mLines = modelImage.getListOfLines();
	vector<Line> sLines = sceneImage.getListOfLines();
	int width = modelImage.getGrayMatrix()->getCols();
	int height = modelImage.getGrayMatrix()->getRows();

	ShapeHistogram mHistogram;
	vector<LocalHistogram> mLocalHist = mHistogram.constructPGH(mLines);
	mHistogram.constructPGHMatrix(mLocalHist, acc, cols);

	ShapeHistogram sHistogram;
	vector<LocalHistogram> sLocalHist = sHistogram.constructPGH(sLines);
	sHistogram.constructPGHMatrix(sLocalHist, acc, cols);

	double bhatt = bhattacharyyaMetric(mHistogram, sHistogram);
	cout << "\nBhattacharrya: " << bhatt << endl;
	//if (bhatt > 0.9)
	//{
	PHoughTransform mpht;
	mpht.setRefPoint(Point(width / 2, height / 2));

	vector<PHTEntry> entriesTable = mpht.constructPHTTable(mLines);
	vector<Point> phtEsLM = phtLandmarks(entriesTable, mpht.getRefPoint(), sLines,
		width, height, manualLMs, angleDiff, ePoint);

	cout << "\n Number of landmarks (pht): " << phtEsLM.size();
	cout << "\nAngle difference: " << angleDiff << endl;
	if (phtEsLM.size() > 0)
	{
		result = verifyLandmarks(modelImage, sceneImage, manualLMs, phtEsLM,
			templSize, sceneSize, angleDiff, ePoint);
		//result = phtEsLM;

	}
	entriesTable.clear();
	phtEsLM.clear();
	//}

	return result;
}

vector<Point> LandmarkDetection::landmarksAutoDectect2(Image &sceneImage,
	int templSize, int sceneSize)
{
	Image modelImage = Treatments::refImage;
	vector<Point> manualLMs = modelImage.getListOfManualLandmarks();
	ProHoughTransform proHT;
	proHT.setRefImage(modelImage);
	double angle = 0;
	Point ePoint, mPoint;
	int rows = sceneImage.getGrayMatrix()->getRows();
	int cols = sceneImage.getGrayMatrix()->getCols();
	ptr_IntMatrix newScene = new Matrix<int>(rows, cols, 0);

	vector<Point> phtEsLM = proHT.generalTransform(sceneImage, angle, ePoint,
		mPoint, newScene);
	vector<Point> result;

	if (phtEsLM.size() > 0)
	{
		result = verifyLandmarks2(modelImage.getGrayMatrix(), newScene, manualLMs,
			phtEsLM, templSize, sceneSize);
	}
	//sceneImage.setGrayMatrix(*newScene);
	phtEsLM.clear();
	delete newScene;
	//result = phtEsLM;
	// reverse the coordinate of estimated landmarks
	Point pi;
	int dx = ePoint.getX() - mPoint.getX();
	int dy = ePoint.getY() - mPoint.getY();
	for (size_t i = 0; i < result.size(); i++)
	{
		pi = result.at(i);
		int xnew = 0, ynew = 0;
		rotateAPoint(pi.getX(), pi.getY(), mPoint, -angle, 1, xnew, ynew);
		xnew += dx;
		ynew += dy;
		result.at(i).setX(xnew);
		result.at(i).setY(ynew);
	}

	return result;
}

void LandmarkDetection::landmarksOnDir(string modelName, string folderScene,
	vector<string> sceneImages, AngleAccuracy acc, int cols, int templSize,
	int sceneSize, Point &ePoint, double &angleDiff, string saveFolder)
{

	Image modelImage = Treatments::refImage;
	vector<Point> manualLMs = modelImage.getListOfManualLandmarks();
	vector<Line> mLines = modelImage.getListOfLines();
	int width = modelImage.getGrayMatrix()->getCols();
	int height = modelImage.getGrayMatrix()->getRows();

	ShapeHistogram mHistogram;
	vector<LocalHistogram> mLocalHist = mHistogram.constructPGH(mLines);
	mHistogram.constructPGHMatrix(mLocalHist, acc, cols);

	ShapeHistogram sHistogram;
	vector<LocalHistogram> sLocalHist;
	Image *sceneImage;
	for (size_t i = 0; i < sceneImages.size(); i++)
	{
		string sceneName = sceneImages.at(i);
		cout << "\n==============================================" << sceneName;
		vector<Point> result;
		angleDiff = 0;
		sceneImage = new Image(folderScene + "/" + sceneName);
		vector<Line> sLines = sceneImage->getListOfLines();

		sLocalHist = sHistogram.constructPGH(sLines);
		sHistogram.constructPGHMatrix(sLocalHist, acc, cols);

		double bhatt = bhattacharyyaMetric(mHistogram, sHistogram);
		cout << "\nBhattacharrya: " << bhatt << endl;
		PHoughTransform mpht;
		mpht.setRefPoint(Point(width / 2, height / 2));

		vector<PHTEntry> entriesTable = mpht.constructPHTTable(mLines);
		vector<Point> phtEsLM = phtLandmarks(entriesTable, mpht.getRefPoint(),
			sLines, width, height, manualLMs, angleDiff, ePoint);

		cout << "\n Number of landmarks (pht): " << phtEsLM.size();
		cout << "\nAngle difference: " << angleDiff << endl;

		if (phtEsLM.size() > 0)
		{
			result = verifyLandmarks(modelImage, *sceneImage, manualLMs, phtEsLM,
				templSize, sceneSize, angleDiff, ePoint);

		}
		entriesTable.clear();
		phtEsLM.clear();

		string saveFile = saveFolder + "/" + modelName + "_" + sceneName + ".TPS";
		ofstream inFile(saveFile.c_str());
		inFile << "LM=" << result.size() << "\n";
		Point pk;
		for (size_t k = 0; k < result.size(); k++)
		{
			pk = result.at(k);
			inFile << pk.getX() << "\t" << pk.getY() << "\n";
		}
		inFile << "IMAGE=" << saveFile << "\n";
		inFile.close();
		delete sceneImage;
	}

}
void LandmarkDetection::landmarksOnDir2(string modelName, string folderScene,
	vector<string> sceneImages, string saveFolder)
{

	Image modelImage = Treatments::refImage;

	int rows = modelImage.getGrayMatrix()->getRows();
	int cols = modelImage.getGrayMatrix()->getCols();
	ptr_IntMatrix mgradirection = new Matrix<int>(rows, cols, -1);
	vector<Point> modelPoints;
	*mgradirection = *(getGradientDMatrix(modelImage, modelPoints));
	vector<Point> mLandmarks = modelImage.getListOfManualLandmarks();

	for (size_t i = 0; i < sceneImages.size(); i++)
	{
		Image *sceneImage;
		string sceneName = sceneImages.at(i);
		cout << "\n==============================================" << sceneName;

		sceneImage = new Image(folderScene + "/" + sceneName);

		ptr_IntMatrix gradirection = new Matrix<int>(rows, cols, -1);
		vector<Point> scenePoints;
		*gradirection = *(getGradientDMatrix(*sceneImage, scenePoints));

		Point mPoint, ePoint, translation;
		double angle;
		vector<Point> eslm = generalizingHoughTransform(mgradirection, gradirection,
			mLandmarks, ePoint, mPoint, angle, translation);

		ptr_IntMatrix newScene = new Matrix<int>(rows, cols, 0);
		// move the model to the same and rotate the scene
		int dx = translation.getX();
		int dy = translation.getY();
		for (int r = 0; r < rows; r++)
		{
			for (int c = 0; c < cols; c++)
			{
				int value = sceneImage->getGrayMatrix()->getAtPosition(r, c);
				int xnew = c + dx;
				int ynew = r + dy;
				rotateAPoint(c + dx, r + dy, mPoint, angle, 1, xnew, ynew);
				if (xnew >= 0 && ynew >= 0 && ynew < rows && xnew < cols)
				{
					newScene->setAtPosition(ynew, xnew, value);
				}
			}
		}

		vector<Point> estLandmarks;

		if (eslm.size() > 0)
		{
			cout << "\nAngle difference: " << angle << endl;
			cout << "\n Number of landmarks (ght): " << eslm.size();
			estLandmarks = verifyLandmarks2(modelImage.getGrayMatrix(), newScene,
				mLandmarks, eslm, 100, 300);
		}

		cout << "\n Number of landmarks (matching): " << estLandmarks.size();

		// reverse the coordinate of estimated landmarks
		Point pi;
		int dx2 = ePoint.getX() - mPoint.getX();
		int dy2 = ePoint.getY() - mPoint.getY();
		for (size_t i = 0; i < estLandmarks.size(); i++)
		{
			pi = estLandmarks.at(i);
			int xnew = 0, ynew = 0;
			rotateAPoint(pi.getX(), pi.getY(), mPoint, -angle, 1, xnew, ynew);
			xnew += dx2;
			ynew += dy2;
			estLandmarks.at(i).setX(xnew);
			estLandmarks.at(i).setY(ynew);
		}

		string saveFile = saveFolder + "/" + modelName + "_" + sceneName + ".TPS";
		ofstream inFile(saveFile.c_str());
		inFile << "LM=" << estLandmarks.size() << "\n";
		Point epk;
		for (size_t k = 0; k < estLandmarks.size(); k++)
		{
			epk = estLandmarks.at(k);
			inFile << epk.getX() << " " << rows - epk.getY() << "\n";
		}
		inFile << "IMAGE=" << saveFile << "\n";
		inFile.close();

		estLandmarks.clear();
		eslm.clear();

		delete sceneImage;
		delete newScene;
	}

}
