/*
 * SVD.cpp
 *
 *  Created on: Feb 2, 2017
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
#include "../imageModel/Point.h"
#include "../imageModel/Matrix.h"

#include "SVD.h"

vector<float> ptB2(float a, float b, float c)
{
	float delta = (b * b) - (4 * a * c);
	vector<float> result;
	if (delta < 0)
		return result;
	if (delta == 0)
	{
		float n = -b / (2 * a);
		result.push_back(n);
		return result;
	}
	float n1 = (-b + sqrt(delta)) / (2 * a);
	float n2 = (-b - sqrt(delta)) / (2 * a);
	result.push_back(n1);
	result.push_back(n2);
	return result;

}
void svdMethod(Matrix<float> m, Matrix<float> &u, Matrix<float> &s,
	Matrix<float> &v)
{

	Matrix<float> mT = m.transposition(0);
	Matrix<float> mTm = mT.multiply(m, 0);
	int rows = mTm.getRows(); // normally, size of m is 2x2
	int cols = mTm.getCols();
	if (rows == 2 && cols == 2)
	{
		int m00 = mTm.getAtPosition(0, 0);
		int m01 = mTm.getAtPosition(0, 1);
		int m10 = mTm.getAtPosition(1, 0);
		int m11 = mTm.getAtPosition(1, 1);

		int detma = 1;
		int detmb = -(m00 + m11);
		int detmc = (m00 * m11) - (m01*m10);
		vector<float> ws = ptB2(detma,detmb,detmc);
		if(ws.size() > 0)
		{
			float k = 0;
			for (size_t i = 0; i < ws.size(); i++) {
				k = ws.at(i);
				float a1 = m00 - k;
				float b1 = m01;
				float a2 = m10;
				float b2 = m11 - k;
				float delta = (a1*b2) - (a2*b1);
				float x,y;
				if(delta ==0)
				{
					y = 1;
					x = (-b1 *y)/a1;
					float xy = (x * x) + (y*y);
					x = x/xy;
					y = y/xy;
					v.setAtPosition(0,i,x);
					v.setAtPosition(1,i,y);
				}
				s.setAtPosition(i,i,k);
			}
			Matrix<float> mv = m.multiply(v,0.0);
			float mva = mv.getAtPosition(0,0);
			float mvb = mv.getAtPosition(0,1);
			float mvc = mv.getAtPosition(1,0);
			float mvd = mv.getAtPosition(1,1);
			float wa = s.getAtPosition(0,0);
			float wb = s.getAtPosition(1,1);
			float ua = mva/wa;
			float ub = mvb/wb;
			float uc = mvc/wa;
			float ud = mvd/wb;
			float uac = (ua * ua) + (uc *uc);
			ua = ua/uac;
			uc = uc/uac;
			float ubd = (ub * ub) + (ud * ud);
			ub = ub/ubd;
			ud = ud/ubd;
			u.setAtPosition(0,0,ua);
			u.setAtPosition(0,1,ub);
			u.setAtPosition(1,0,uc);
			u.setAtPosition(1,1,ud);
		}
	}
}

