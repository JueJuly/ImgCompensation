
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

//���غ�������м�����ȡƽ��ֵ��Ȼ����������������Sͨ����blending
int ImgSatCompenTest1();

//�������غ�������S��ƽ��ֵ��Ȼ�����alpha �ںϲ���
int ImgSatCompenTest2();

//������ͼ����ȡS��ƽ��ֵ��Ȼ��ʹ��ƽ��ֵ�����滻����
int ImgSatCompenTest3();

