/*
 * Analysis.h
 *
 *  Created on: Dec 7, 2016
 *      Author: linh
 */

#ifndef MAELAB_H_
#define MAELAB_H_
vector<Line> segment(ptr_Treatments treatment, Image sceneImage,
	int minDistance);
ShapeHistogram geomHistogram(ptr_Treatments treatment, Image image,
	AngleAccuracy angleAcc, int cols);
double bhattacharyyaDistance(ptr_Treatments treatment, Image sceneImage,
	AngleAccuracy angleAcc, int cols);
PHoughTransform phtEntriesTable(ptr_Treatments treatment, Image image);
vector<Point> estimatedLandmarks(ptr_Treatments treatment, Image sceneImage,
	AngleAccuracy acc, int cols, int templSize, int sceneSize,Point &ePoint, double angleDiff);
double measureCentroidPoint(vector<Point> landmarks, Point &ebary);
#endif /* ANALYSIS_H_ */
