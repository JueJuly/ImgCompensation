
#pragma once

#include "common.h"

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include <opencv.hpp>

using namespace cv;
using namespace std;

#ifndef uchar
typedef unsigned char uchar;
#endif
#define WS_SWAP(a,b,t) ((t) = (a), (a) = (b), (b) = (t))
#define WS_MIN(a,b) ( ((a)>(b)) ? (b) : (a) )
#define WS_MAX(a,b) ( ((a)<(b)) ? (b) : (a) )
#define WS_ABS(x) ((x) >= 0 ? (x) : (-(x)))

#define UNDEFINEDCOLOR -1

void RGB2HSV(int *h, int *s, int *v, uchar r, uchar g, uchar b);

void HSV2RGB(uchar *r, uchar *g, uchar *b, int h, int s, int v);

//Calculate the average saturation of the pixels on the limit line
float CalSatMeanVal(uchar *r, uchar *g, uchar *b, int num);

//Image Saturation gain compensation process;
void ImgSatGainCompensate(uchar *rgb, int width, int height, float gainVal);


void testColorTransfer();


int ImgSatCompensationTest();

//对重合区域的中间线求取平均值，然后根据左右相机的作S通道的blending
int ImgSatCompenTest1();

//对整个重合区域求S的平均值；然后进行alpha 融合补偿
int ImgSatCompenTest2();

//对整个图像求取S的平均值，然后使用平均值进行替换补偿
int ImgSatCompenTest3();

//F_S(i, j) =  L_mean + (L_std / F_std) * (F_S(i, j)-  F_mean );
int ImgSatCompenTest4();

//使用OpenCV函数;
int ImgSatCompenTest5();

//使用OpenCV函数;
int ImgSatCompenTest6();

//按照shader语言的编程格式修改RGB2BGR
int Test7();

	

//shader language RGB <===> HSV
//vec3 rgb2hsv(vec3 c)
//{
//	vec3 hsv;
//	vec3 rgb = vec3(c.x / 255.0, c.y / 255.0, c.z / 255.0);
//	vec3 tempv = vec3(max(rgb.x, max(rgb.y, rgb.z)), min(rgb.x, min(rgb.y, rgb.z)), 0.0);
//	tempv.z = tempv.x - tempv.y;
//	if (tempv.x < -0.00001 || tempv.x > 0.00001)
//	{
//		if ((rgb.x > (tempv.x - 0.00001)) && (rgb.x < (tempv.x + 0.00001)))
//		{
//			hsv.x = (rgb.y - rgb.z) / tempv.z;
//		}
//		else if ((rgb.y > (tempv.x - 0.00001)) && (rgb.y < (tempv.x + 0.00001)))
//		{
//			hsv.x = 2.0 + (rgb.z - rgb.x) / tempv.z;
//		}
//		else
//		{
//			hsv.x = 4.0 + (rgb.x - rgb.y) / tempv.z;
//		}
//
//		hsv.x *= 60.0;
//		if (hsv.x < 0.0)
//		{
//			hsv.x += 360.0;
//		}
//
//		hsv = vec3(hsv.x, tempv.z / tempv.x*100.0, tempv.x*100.0);
//	}
//	else
//	{
//		hsv = vec3(-1.0, 0.0, tempv.x);
//	}
//	return hsv;
//}
//
//vec3 hsv2rgb(vec3 c)
//{
//	vec3 hsv = c;
//	vec3 rgb;
//	vec3 tempv = vec3(hsv.z*2.55, (hsv.z*2.55)*(100.0 - hsv.y) / 100.0, 0.0);
//	tempv.z = (tempv.x - tempv.y)*mod(hsv.x, 60.0) / 60.0;
//
//	if (hsv.x > -1.00001 && hsv.x < -0.99999)
//	{
//		rgb = vec3(0.0, 0.0, 0.0);
//	}
//
//	hsv.x = (hsv.x - mod(hsv.x, 60.0)) / 60.0;
//
//	if (hsv.x > -0.00001 && hsv.x < 0.00001)
//	{
//		rgb = vec3(tempv.x, tempv.y + tempv.z, tempv.y);
//	}
//	else if (hsv.x > 0.99999 && hsv.x < 1.00001)
//	{
//		rgb = vec3(tempv.x - tempv.z, tempv.x, tempv.y);
//	}
//	else if (hsv.x > 1.99999 && hsv.x < 2.00001)
//	{
//		rgb = vec3(tempv.y, tempv.x, tempv.y + tempv.z);
//	}
//	else if (hsv.x > 2.99999 && hsv.x < 3.00001)
//	{
//		rgb = vec3(tempv.y, tempv.x - tempv.z, tempv.x);
//	}
//	else if (hsv.x > 3.99999 && hsv.x < 4.00001)
//	{
//		rgb = vec3(tempv.y + tempv.z, tempv.y, tempv.x);
//	}
//	else if (hsv.x > 4.99999 && hsv.x < 5.00001)
//	{
//		rgb = vec3(tempv.y , tempv.y, tempv.x - tempv.z);
//	}
//	else
//	{
//		rgb = vec3(tempv.x, tempv.y, tempv.z);
//	}
//	rgb.x = rgb.x > 255.0 ? 255.0 : rgb.x;
//	rgb.y = rgb.y > 255.0 ? 255.0 : rgb.y;
//	rgb.z = rgb.z > 255.0 ? 255.0 : rgb.z;
//
//	rgb.x = rgb.x < 0.0 ? 0.0 : rgb.x;
//	rgb.y = rgb.y < 0.0 ? 0.0 : rgb.y;
//	rgb.z = rgb.z < 0.0 ? 0.0 : rgb.z;
//
//	return rgb;
//}

