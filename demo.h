
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

