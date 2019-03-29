
#include "demo.h"

int main(void)
{
	/*Gcout;*/ cout << string(32, '-') << endl;/* Wcout;*/
	/*Rcout;*/ cout << "对图像饱和度补偿测试!" << endl; /*Wcout;*/

	//ImgSatCompensationTest();
	//testColorTransfer();

	//ImgSatCompenTest1();
	ImgSatCompenTest2();
	//ImgSatCompenTest3();

	/*Gcout; */cout << string(32, '-') << endl; /*Wcout;*/
	system("pause");
	return 0;
}


int ImgSatCompensationTest()
{
	//string imgPath = "E:/matlab_project/mdlt/mdlt/images/case1/3.jpg";
	string strFrontImg = "E:/matlab_project/mdlt/mdlt/diffAngleImg/20DU_front.png";

	Mat srcHSV, sat, satAdj, dstMerge, dst;     //sat - saturation饱和度分量

#if 1
	Mat imageAwb, srcImg;
	srcImg = imread(strFrontImg, IMREAD_COLOR);

	//从原始图像中裁剪重合区域;
	Rect OverlapRect(0, 0, 320, 320);

	imageAwb = srcImg(OverlapRect).clone();
#else
	Mat imageAwb = imread(strFrontImg, IMREAD_COLOR);
#endif

	vector<Mat> channels, channels1;
	double p1, p2, p3;

	imshow("srcImg", imageAwb);


	cvtColor(imageAwb, srcHSV, CV_BGR2HSV);
	split(srcHSV, channels);
	split(srcHSV, channels1);
	sat = channels.at(1);
	Scalar m = mean(sat);

	if (m(0) <= 51.5)
	{
		p1 = -0.002714, p2 = 0.9498, p3 = -0.5073;  
		Ycout; cout << "High Color Enhancement!" << endl; Wcout;
	}//高
	else if (m(0) > 38.5 && m(0) <= 89.5)
	{
		p1 = -0.001578, p2 = 0.6238, p3 = -0.1514; 
		Ycout; cout << "Middle Color Enhancement!" << endl; Wcout;
	}//中
	else if (m(0) > 89.5 && m(0) <= 127.5)
	{
		p1 = -0.0006512, p2 = 0.2639, p3 = -0.9246; 
		Ycout; cout << "Low Color Enhancement!" << endl; Wcout;
	}//低
	else
	{
		p1 = 0, p2 = 0, p3 = 0; 
		Ycout; cout << "No Color Enhancement!" << endl; Wcout;
	}

	satAdj = sat;
	for (int i = 0; i < sat.rows; i++)
	{
		for (int j = 0; j < sat.cols; j++)
		{
			/*if (i > sat.rows / 2 && i < sat.rows * 3 / 4 &&
				j > sat.cols / 2 && j < sat.cols * 3 / 4)
			{
				uchar val = sat.at<uchar>(i, j);
				satAdj.at<uchar>(i, j) = (uchar)(val + p1 * val * val + p2 * val + p3);
			}*/
			uchar val = sat.at<uchar>(i, j);
			satAdj.at<uchar>(i, j) = (uchar)(val + p1 * val * val + p2 * val + p3);
		}
	}

	channels1.at(1) = satAdj;
	merge(channels1, dstMerge);
	cvtColor(dstMerge, dst, CV_HSV2BGR);

	imshow("dstImg", dst);
	waitKey(0);

	//imwrite("m_ImageCE.bmp", dst);
	destroyAllWindows();

	return 0;
}

/*
 @brief convert from RGB to HSV color
 @param inout	int   *h	hue 
 @param inout	int   *s	saturation 
 @param inout	int	  *v	value
 @param in		uchar r		red 
 @param in		uchar g		green 
 @param in		uchar b		blue 
 @note ref https://blog.csdn.net/hzdiy/article/details/8734229
	r,g,b values are from 0 to 255
	h = [0,360],s = [0,100], v = [0,100]
	if s == 0, then h = -1(undefined)
 */
void RGB2HSV(int *h, int *s, int *v, uchar r, uchar g, uchar b)
{

	float min, max, delta, tmp;
	float R = r / 255.0, G = g / 255.0, B = b / 255.0;
	float H, S, V;
	tmp = WS_MIN(R, G);
	min = WS_MIN(tmp, B);
	tmp = WS_MAX(R, G);
	max = WS_MAX(tmp, B);
	V = max; //v
	delta = max - min;

	if (max != 0 /*&& delta != 0*/)
	{ 
		S = delta / max; 

		if (R == max)
			H = (G - B) / delta;
		else if (G == max)
			H = 2 + (B - R) / delta;
		else
			H = 4 + (R - G) / delta;

		H *= 60;
		if (H < 0)
			H += 360;
	}
	else
	{
		//r = g = b = 0; //s = 0, v is undefined
		S = 0;
		H = UNDEFINEDCOLOR;
		return;
	}

	*h = (int)(H);
	*s = (int)(S * 100);
	*v = (int)(V * 100);
	
}

/*
 @brief convert from HSV to RGB  color
 @param inout	uchar   *r	red
 @param inout	uchar   *g	green
 @param inout	uchar	*b	blue
 @param in		int     h	hue
 @param in		int     s	saturation
 @param in		int     v	value
 @note ref https://blog.csdn.net/hzdiy/article/details/8734229
	R,G,B from 0-255, H from 0-360, S,V from 0-100
 */
void HSV2RGB(uchar *r, uchar *g, uchar *b, int h, int s, int v)
{

	int i;

	float RGB_min, RGB_max;
	RGB_max = v * 2.55f;
	RGB_min = RGB_max * (100 - s) / 100.0f;

	if (h == UNDEFINEDCOLOR) //black 
	{
		*r = 0;
		*g = 0;
		*b = 0;
		return;
	}

	i = h / 60;
	int difs = h % 60; // factorial part of h

	// RGB adjustment amount by hue 
	float RGB_Adj = (RGB_max - RGB_min) * difs / 60.0f;

	switch (i) {
	case 0:
		*r = (uchar)RGB_max;
		*g = (uchar)(RGB_min + RGB_Adj);
		*b = (uchar)RGB_min;
		break;
	case 1:
		*r = (uchar)(RGB_max - RGB_Adj);
		*g = (uchar)RGB_max;
		*b = (uchar)RGB_min;
		break;
	case 2:
		*r = (uchar)RGB_min;
		*g = (uchar)RGB_max;
		*b = (uchar)(RGB_min + RGB_Adj);
		break;
	case 3:
		*r = (uchar)RGB_min;
		*g = (uchar)(RGB_max - RGB_Adj);
		*b = (uchar)RGB_max;
		break;
	case 4:
		*r = (uchar)(RGB_min + RGB_Adj);
		*g = (uchar)RGB_min;
		*b = (uchar)RGB_max;
		break;
	default:		// case 5:
		*r = (uchar)RGB_max;
		*g = (uchar)RGB_min;
		*b = (uchar)(RGB_max - RGB_Adj);
		break;
	}
}


/*
 @brief Calculate the average saturation of the pixels on the limit line 
 @param in uchar *r		red pointer
 @param in uchar *g     green pointer
 @param in uchar *b     blue pointer
 @param in int   num    the number of the pixels on the limit line
 @return float type the average saturation value
*/
float CalSatMeanVal(uchar *r, uchar *g, uchar *b, int num)
{
	float meanVal = 0.0f;
	uchar R = 0, G = 0, B = 0;
	int H = 0, S = 0, V = 0;

	for (int i = 0; i < num; i++)
	{
		H = 0;
		S = 0;
		V = 0;
		R = *(r + i);
		G = *(g + i);
		B = *(b + i);
		RGB2HSV(&H, &S, &V, R, G, B);
		cout << i << " : " << S << endl;
		meanVal += (float)(S);

	}

	meanVal /= num;

	return meanVal;
}

/*
 @brief Image Saturation gain compensation process;
 @param inout uchar  *rgb     color image pixel value pointer
 @param in    int    width    image width
 @param in    int    height   image height
 @param in    float  gainVal  gain compensation value;
 @note
 */
void ImgSatGainCompensate(uchar *rgb, int width, int height, float gainVal)
{
	for (int y = 0; y < height; y++)
	{
		int indexY = y * width * 3;//get row index;

		for (int x = 0; x < width; x++)
		{
			int indexX = x * 3; //get col index;

			//get source image RGB pixel value
			uchar R = rgb[indexY + indexX + 0];
			uchar G = rgb[indexY + indexX + 1];
			uchar B = rgb[indexY + indexX + 2];

			int H = 0, S = 0, V = 0;

			//invalid region;
			if (R == 0 && G == 0 && B == 0)
				continue;

			//gain compensation process
			RGB2HSV(&H, &S, &V, R, G, B);
			S = (int)(gainVal);
			HSV2RGB(&R, &G, &B, H, S, V);

			//modify source image pixel value;
			rgb[indexY + indexX + 0] = R;
			rgb[indexY + indexX + 1] = G;
			rgb[indexY + indexX + 2] = B;
		}
	}
}

/*
 @brief test RGB2HSV and HSV2RGB function;
 */
void testColorTransfer()
{
	string imgPath = "E:/matlab_project/mdlt/mdlt/images/case1/3.jpg";

	//Mat srcHSV, sat, satAdj, dstMerge, dst;     //sat - saturation饱和度分量
	Mat srcImg = imread(imgPath, IMREAD_COLOR);
	Mat dstImg = Mat::zeros(srcImg.rows, srcImg.cols, CV_8UC3);

	for (int y = 0; y < srcImg.rows; y++)
	{
		for (int x = 0; x < srcImg.cols; x++)
		{
			//读取图像的像素值
			uchar b = srcImg.at<Vec3b>(y, x)[0]; // blue
			uchar g = srcImg.at<Vec3b>(y, x)[1]; // green
			uchar r = srcImg.at<Vec3b>(y, x)[2]; // red

			int H=0, S=0, V=0;

			RGB2HSV(&H, &S, &V, r, g, b);
			//S -= 50; //for testing;
			HSV2RGB(&r, &g, &b, H, S, V);

			dstImg.at<Vec3b>(y, x)[0] = b;
			dstImg.at<Vec3b>(y, x)[1] = g;
			dstImg.at<Vec3b>(y, x)[2] = r;
		}
	}

	imshow("srcImg", srcImg);
	imshow("dstImg", dstImg);
	waitKey(0);

	destroyAllWindows();

}

struct ValidRegion
{
	Point stPt;
	Point edPt;
};

//对重合区域的中间线求取平均值，然后根据左右相机的作S通道的blending
int ImgSatCompenTest1()
{
	string strLeftImg = "E:/matlab_project/mdlt/mdlt/diffAngleImg/20DU_left.png";
	string strFrontImg = "E:/matlab_project/mdlt/mdlt/diffAngleImg/20DU_front.png";

	//读原始图像;
	Mat leftImg = imread(strLeftImg, 1);
	Mat frontImg = imread(strFrontImg, 1);

	//从原始图像中裁剪重合区域;
	Rect OverlapRect(0, 0, 320, 320);

	vector<ValidRegion> leftImgValidReg;
	vector<ValidRegion> frontImgValidReg;

	Mat leftImgOverlap = leftImg(OverlapRect).clone();
	Mat frontImgOverlap = frontImg(OverlapRect).clone();

	Mat leftOverlapImgRGB, frontOverlapImgRGB;

	//cvtColor(leftImgOverlap, leftOverlapImgRGB, CV_BGR2RGB);
	//cvtColor(frontImgOverlap, frontOverlapImgRGB, CV_BGR2RGB);

	imshow("srcLeftImg", leftImgOverlap);
	imshow("srcFrontImg", frontImgOverlap);

	leftImgValidReg.clear();
	frontImgValidReg.clear();

	//得到有效区域的起始和结束点坐标;
	for (int y = 0; y < OverlapRect.height; y++)
	{
		bool stFlagL = false;
		bool edFlagL = false;
		ValidRegion validReg;

		for (int x = 3; x < OverlapRect.width-3; x++)
		{
			Vec3b PrePixel_L = leftImgOverlap.at<Vec3b>(y, x-1); //previous pixel for left image;
			Vec3b CurPixel_L = leftImgOverlap.at<Vec3b>(y, x);   //current pixel for left image;

			//得到左相机有效区域的起始、结束点;
			//----------------------------------------------
			if ((PrePixel_L[0] == 0 && PrePixel_L[1] == 0 && PrePixel_L[2] == 0) &&
				(CurPixel_L[0] != 0 || CurPixel_L[1] != 0 || CurPixel_L[2] != 0))
			{
				validReg.stPt.x = x;
				validReg.stPt.y = y;
				//cout << "start point : " << Point(x, y) << endl;
				stFlagL = true;
			}

			if ((PrePixel_L[0] != 0 || PrePixel_L[1] != 0 || PrePixel_L[2] != 0) &&
				(CurPixel_L[0] == 0 && CurPixel_L[1] == 0 && CurPixel_L[2] == 0))
			{
				validReg.edPt = Point(x, y);
				edFlagL = true;
			}

			if (stFlagL == true && edFlagL == true)
			{
				leftImgValidReg.push_back(validReg);
				stFlagL = false;
				edFlagL = false;
				break;
			}
				
		}
	}

	//-----------------------------------
	for (int y = 0; y < OverlapRect.height; y++)
	{
		bool stFlagF = false;
		bool edFlagF = false;
		ValidRegion validReg;

		for (int x = 3; x < OverlapRect.width - 3; x++)
		{
			Vec3b PrePixel_F = frontImgOverlap.at<Vec3b>(y, x - 1); //previous pixel for front image;
			Vec3b CurPixel_F = frontImgOverlap.at<Vec3b>(y, x);   //current pixel for front image;

			//得到前相机有效区域的起始、结束点;
			//-----------------------------------------------
			if ((PrePixel_F[0] == 0 && PrePixel_F[1] == 0 && PrePixel_F[2] == 0) &&
				(CurPixel_F[0] != 0 || CurPixel_F[1] != 0 || CurPixel_F[2] != 0))
			{
				validReg.stPt = Point(x, y);
				stFlagF = true;
			}

			if ((PrePixel_F[0] != 0 || PrePixel_F[1] != 0 || PrePixel_F[2] != 0) &&
				(CurPixel_F[0] == 0 && CurPixel_F[1] == 0 && CurPixel_F[2] == 0))
			{
				validReg.edPt = Point(x, y);
				edFlagF = true;
			}

			if (stFlagF == true && edFlagF == true)
			{
				frontImgValidReg.push_back(validReg);
				stFlagF = false;
				edFlagF = false;
				break;
			}

		}
	}

	vector<Point> limitLinePtL;
	vector<Point> limitLinePtF;

	limitLinePtF.clear();
	for (auto &i : frontImgValidReg)
	{
		Point midPt;
		midPt.x = (int)((i.stPt.x + i.edPt.x) / 2);
		midPt.y = (int)((i.stPt.y + i.edPt.y) / 2);

		limitLinePtF.push_back(midPt);
	}

	limitLinePtL.clear();
	for (auto &i : leftImgValidReg)
	{
		Point midPt;
		midPt.x = (int)((i.stPt.x + i.edPt.x) / 2);
		midPt.y = (int)((i.stPt.y + i.edPt.y) / 2);

		limitLinePtL.push_back(midPt);
	}

	//int num = WS_MAX(, limitLinePtL.size());
	int numF = limitLinePtF.size();
	int numL = limitLinePtL.size();

	//计算有效区域中间线的平均值;
	uchar *R_F = new uchar[numF];
	uchar *G_F = new uchar[numF];
	uchar *B_F = new uchar[numF];

	uchar *R_L = new uchar[numL];
	uchar *G_L = new uchar[numL];
	uchar *B_L = new uchar[numL];


	memset(R_F, 0, sizeof(uchar) * numF);
	memset(G_F, 0, sizeof(uchar) * numF);
	memset(B_F, 0, sizeof(uchar) * numF);

	memset(R_L, 0, sizeof(uchar) * numL);
	memset(G_L, 0, sizeof(uchar) * numL);
	memset(B_L, 0, sizeof(uchar) * numL);

	int k = 0;
	for (auto &i : limitLinePtF)
	{
		R_F[k] = frontImgOverlap.at<Vec3b>(i.y, i.x)[2];
		G_F[k] = frontImgOverlap.at<Vec3b>(i.y, i.x)[1];
		B_F[k] = frontImgOverlap.at<Vec3b>(i.y, i.x)[0];
		k++;
	}

	k = 0;
	for (auto &i : limitLinePtL)
	{
		R_L[k] = leftImgOverlap.at<Vec3b>(i.y, i.x)[2];
		G_L[k] = leftImgOverlap.at<Vec3b>(i.y, i.x)[1];
		B_L[k] = leftImgOverlap.at<Vec3b>(i.y, i.x)[0];
		k++;
	}

	float SMeanVal_F = CalSatMeanVal(R_F, G_F, B_F, numF);
	cout << "==================" << endl;
	float SMeanVal_L = CalSatMeanVal(R_L, G_L, B_L, numL);


	//ImgSatGainCompensate(leftOverlapImgRGB.data, leftOverlapImgRGB.cols, leftOverlapImgRGB.rows, SMeanVal_L);
	//ImgSatGainCompensate(frontOverlapImgRGB.data, frontOverlapImgRGB.cols, frontOverlapImgRGB.rows, SMeanVal_F);
	//根据limit line对重合区域做 饱和度过度补偿
	//公式: k*S1 + (1-k)*S2

	//---对左图像处理
	size_t validNum_L = leftImgValidReg.size();
	for (auto i = 0; i < validNum_L; i++)
	{
		int num = leftImgValidReg[i].edPt.x - leftImgValidReg[i].stPt.x + 1;
		float step = 1.0 / num;
		int y = leftImgValidReg[i].edPt.y;
		int stX = leftImgValidReg[i].stPt.x;
		int edX = leftImgValidReg[i].edPt.x;
		for (int x = stX; x <= edX; x++)
		{
			float k = (edX - x)*step;
			uchar r = leftImgOverlap.at<Vec3b>(y, x)[2];
			uchar g = leftImgOverlap.at<Vec3b>(y, x)[1];
			uchar b = leftImgOverlap.at<Vec3b>(y, x)[0];
			int h = 0, s = 0, v = 0;
			RGB2HSV(&h, &s, &v, r, g, b);
			s = (int)(k * SMeanVal_L + (1 - k)*SMeanVal_F);
			//s = (int)(SMeanVal_L + 0.5);
			HSV2RGB(&r, &g, &b, h, s, v);
			leftImgOverlap.at<Vec3b>(y, x)[2] = r;
			leftImgOverlap.at<Vec3b>(y, x)[1] = g;
			leftImgOverlap.at<Vec3b>(y, x)[0] = b;
		}
	}

	//对前相机图像进行处理;
	size_t validNum_F = frontImgValidReg.size();
	for (auto i = 0; i < validNum_F; i++)
	{
		int num = frontImgValidReg[i].edPt.x - frontImgValidReg[i].stPt.x + 1;
		float step = 1.0 / num;
		int y = frontImgValidReg[i].edPt.y;
		int stX = frontImgValidReg[i].stPt.x;
		int edX = frontImgValidReg[i].edPt.x;
		for (int x = stX; x <= edX; x++)
		{
			float k = (edX - x)*step;
			uchar r = frontImgOverlap.at<Vec3b>(y, x)[2];
			uchar g = frontImgOverlap.at<Vec3b>(y, x)[1];
			uchar b = frontImgOverlap.at<Vec3b>(y, x)[0];
			int h = 0, s = 0, v = 0;
			RGB2HSV(&h, &s, &v, r, g, b);
			s = (int)(k * SMeanVal_L + (1 - k)*SMeanVal_F );
			//s = (int)(SMeanVal_F + 0.5);
			HSV2RGB(&r, &g, &b, h, s, v);
			frontImgOverlap.at<Vec3b>(y, x)[2] = r;
			frontImgOverlap.at<Vec3b>(y, x)[1] = g;
			frontImgOverlap.at<Vec3b>(y, x)[0] = b;
		}
	}

	imshow("dstLeftImg", leftImgOverlap);
	imshow("dstFrontImg", frontImgOverlap);
	waitKey(0);

	destroyAllWindows();

	delete[]R_F; R_F = NULL;
	delete[]G_F; G_F = NULL;
	delete[]B_F; B_F = NULL;

	delete[]R_L; R_L = NULL;
	delete[]G_L; G_L = NULL;
	delete[]B_L; B_L = NULL;

	return 0;
}

/*
 @brief 对整个重合区域求S的平均值；然后进行alpha 融合补偿
 */
int ImgSatCompenTest2()
{
	string strLeftImg = "E:/matlab_project/mdlt/mdlt/diffAngleImg/20DU_left.png";
	string strFrontImg = "E:/matlab_project/mdlt/mdlt/diffAngleImg/20DU_front.png";

	//读原始图像;
	Mat leftImg = imread(strLeftImg, 1);
	Mat frontImg = imread(strFrontImg, 1);

	//从原始图像中裁剪重合区域;
	Rect OverlapRect(0, 0, 320, 320);

	Mat leftImgOverlap = leftImg(OverlapRect).clone();
	Mat frontImgOverlap = frontImg(OverlapRect).clone();

	imshow("srcLeftImg", leftImgOverlap);
	imshow("srcFrontImg", frontImgOverlap);

	vector<ValidRegion> leftImgValidReg;
	vector<ValidRegion> frontImgValidReg;

	//-------------------------------------------------
	//得到有效区域的起始和结束点坐标;
	for (int y = 0; y < OverlapRect.height; y++)
	{
		bool stFlagL = false;
		bool edFlagL = false;
		ValidRegion validReg;

		for (int x = 3; x < OverlapRect.width - 3; x++)
		{
			Vec3b PrePixel_L = leftImgOverlap.at<Vec3b>(y, x - 1); //previous pixel for left image;
			Vec3b CurPixel_L = leftImgOverlap.at<Vec3b>(y, x);   //current pixel for left image;

			//得到左相机有效区域的起始、结束点;
			//----------------------------------------------
			if ((PrePixel_L[0] == 0 && PrePixel_L[1] == 0 && PrePixel_L[2] == 0) &&
				(CurPixel_L[0] != 0 || CurPixel_L[1] != 0 || CurPixel_L[2] != 0))
			{
				validReg.stPt.x = x;
				validReg.stPt.y = y;
				//cout << "start point : " << Point(x, y) << endl;
				stFlagL = true;
			}

			if ((PrePixel_L[0] != 0 || PrePixel_L[1] != 0 || PrePixel_L[2] != 0) &&
				(CurPixel_L[0] == 0 && CurPixel_L[1] == 0 && CurPixel_L[2] == 0))
			{
				validReg.edPt = Point(x, y);
				edFlagL = true;
			}

			if (stFlagL == true && edFlagL == true)
			{
				leftImgValidReg.push_back(validReg);
				stFlagL = false;
				edFlagL = false;
				break;
			}

		}
	}

	//-----------------------------------
	for (int y = 0; y < OverlapRect.height; y++)
	{
		bool stFlagF = false;
		bool edFlagF = false;
		ValidRegion validReg;

		for (int x = 3; x < OverlapRect.width - 3; x++)
		{
			Vec3b PrePixel_F = frontImgOverlap.at<Vec3b>(y, x - 1); //previous pixel for front image;
			Vec3b CurPixel_F = frontImgOverlap.at<Vec3b>(y, x);   //current pixel for front image;

			//得到前相机有效区域的起始、结束点;
			//-----------------------------------------------
			if ((PrePixel_F[0] == 0 && PrePixel_F[1] == 0 && PrePixel_F[2] == 0) &&
				(CurPixel_F[0] != 0 || CurPixel_F[1] != 0 || CurPixel_F[2] != 0))
			{
				validReg.stPt = Point(x, y);
				stFlagF = true;
			}

			if ((PrePixel_F[0] != 0 || PrePixel_F[1] != 0 || PrePixel_F[2] != 0) &&
				(CurPixel_F[0] == 0 && CurPixel_F[1] == 0 && CurPixel_F[2] == 0))
			{
				validReg.edPt = Point(x, y);
				edFlagF = true;
			}

			if (stFlagF == true && edFlagF == true)
			{
				frontImgValidReg.push_back(validReg);
				stFlagF = false;
				edFlagF = false;
				break;
			}

		}
	}
	//-------------------------------------------------

	/*int numF = frontImgValidReg.size();
	int numL = leftImgValidReg.size();*/

	//计算有效区域中间线的平均值;
	//uchar *R_F = new uchar[numF];
	//uchar *G_F = new uchar[numF];
	//uchar *B_F = new uchar[numF];

	//uchar *R_L = new uchar[numL];
	//uchar *G_L = new uchar[numL];
	//uchar *B_L = new uchar[numL];


	//memset(R_F, 0, sizeof(uchar) * numF);
	//memset(G_F, 0, sizeof(uchar) * numF);
	//memset(B_F, 0, sizeof(uchar) * numF);

	//memset(R_L, 0, sizeof(uchar) * numL);
	//memset(G_L, 0, sizeof(uchar) * numL);
	//memset(B_L, 0, sizeof(uchar) * numL);

	float meanValF = 0;
	float meanValL = 0;
	int numL = 0;
	int numF = 0;

	for (int y = 0; y < OverlapRect.height; y++)
	{
		for (int x = 3; x < OverlapRect.width - 3; x++)
		{
			uchar r_F = frontImgOverlap.at<Vec3b>(y, x)[2];
			uchar g_F = frontImgOverlap.at<Vec3b>(y, x)[1];
			uchar b_F = frontImgOverlap.at<Vec3b>(y, x)[0];

			uchar r_L = leftImgOverlap.at<Vec3b>(y, x)[2];
			uchar g_L = leftImgOverlap.at<Vec3b>(y, x)[1];
			uchar b_L = leftImgOverlap.at<Vec3b>(y, x)[0];

			int h_F = 0, s_F = 0, v_F = 0;
			int h_L = 0, s_L = 0, v_L = 0;

			if (r_F != 0 || g_F != 0 || b_F != 0)
			{
				numF++;
				RGB2HSV(&h_F, &s_F, &v_F, r_F, g_F, b_F);
				meanValF += s_F;
			}

			if (r_L != 0 || g_L != 0 || b_L != 0)
			{
				numL++;
				RGB2HSV(&h_L, &s_L, &v_L, r_L, g_L, b_L);
				meanValL += s_L;
			}
		}
	}

	meanValF /= numF;
	meanValL /= numL;

	//for (int y = 0; y < OverlapRect.height; y++)
	//{
	//	for (int x = 3; x < OverlapRect.width - 3; x++)
	//	{
	//		uchar r_F = frontImgOverlap.at<Vec3b>(y, x)[2];
	//		uchar g_F = frontImgOverlap.at<Vec3b>(y, x)[1];
	//		uchar b_F = frontImgOverlap.at<Vec3b>(y, x)[0];

	//		uchar r_L = leftImgOverlap.at<Vec3b>(y, x)[2];
	//		uchar g_L = leftImgOverlap.at<Vec3b>(y, x)[1];
	//		uchar b_L = leftImgOverlap.at<Vec3b>(y, x)[0];

	//		int h_F = 0, s_F = 0, v_F = 0;
	//		int h_L = 0, s_L = 0, v_L = 0;

	//		if (r_F != 0 || g_F != 0 || b_F != 0)
	//		{
	//			//numF++;
	//			RGB2HSV(&h_F, &s_F, &v_F, r_F, g_F, b_F);
	//			s_F = meanValF;
	//			HSV2RGB(&r_F, &g_F, &b_F, h_F, s_F, v_F);
	//			frontImgOverlap.at<Vec3b>(y, x)[2] = r_F;
	//			frontImgOverlap.at<Vec3b>(y, x)[1] = g_F;
	//			frontImgOverlap.at<Vec3b>(y, x)[0] = b_F;
	//			//meanValF += s_F;
	//		}

	//		if (r_L != 0 || g_L != 0 || b_L != 0)
	//		{
	//			//numL++;
	//			RGB2HSV(&h_L, &s_L, &v_L, r_L, g_L, b_L);
	//			s_L = meanValL;
	//			HSV2RGB(&r_L, &g_L, &b_L, h_L, s_L, v_L);
	//			//meanValL += s_L;
	//			leftImgOverlap.at<Vec3b>(y, x)[2] = r_L;
	//			leftImgOverlap.at<Vec3b>(y, x)[1] = g_L;
	//			leftImgOverlap.at<Vec3b>(y, x)[0] = b_L;
	//		}
	//	}
	//}

	//---对左图像处理
	size_t validNum_L = leftImgValidReg.size();
	for (auto i = 0; i < validNum_L; i++)
	{
		int num = leftImgValidReg[i].edPt.x - leftImgValidReg[i].stPt.x + 1;
		float step = 1.0 / num;
		int y = leftImgValidReg[i].edPt.y;
		int stX = leftImgValidReg[i].stPt.x;
		int edX = leftImgValidReg[i].edPt.x;
		for (int x = stX; x <= edX; x++)
		{
			float k = (x-stX)*step;
			float diff = (0.5 - WS_ABS(k - 0.5))*(meanValF - meanValL);

			uchar r = leftImgOverlap.at<Vec3b>(y, x)[2];
			uchar g = leftImgOverlap.at<Vec3b>(y, x)[1];
			uchar b = leftImgOverlap.at<Vec3b>(y, x)[0];
			int h = 0, s = 0, v = 0;
			RGB2HSV(&h, &s, &v, r, g, b);
			s = (int)(s + diff);
			//s = (int)(k * meanValL + (1 - k)*meanValF);
			//s = (int)(SMeanVal_L + 0.5);
			HSV2RGB(&r, &g, &b, h, s, v);
			leftImgOverlap.at<Vec3b>(y, x)[2] = r;
			leftImgOverlap.at<Vec3b>(y, x)[1] = g;
			leftImgOverlap.at<Vec3b>(y, x)[0] = b;
		}
	}

	//对前相机图像进行处理;
	size_t validNum_F = frontImgValidReg.size();
	for (auto i = 0; i < validNum_F; i++)
	{
		int num = frontImgValidReg[i].edPt.x - frontImgValidReg[i].stPt.x + 1;
		float step = 1.0 / num;
		int y = frontImgValidReg[i].edPt.y;
		int stX = frontImgValidReg[i].stPt.x;
		int edX = frontImgValidReg[i].edPt.x;
		for (int x = stX; x <= edX; x++)
		{
			float k = (x-stX)*step;
			float diff = (0.5 - WS_ABS(k - 0.5))*(meanValF - meanValL);

			uchar r = frontImgOverlap.at<Vec3b>(y, x)[2];
			uchar g = frontImgOverlap.at<Vec3b>(y, x)[1];
			uchar b = frontImgOverlap.at<Vec3b>(y, x)[0];
			int h = 0, s = 0, v = 0;
			RGB2HSV(&h, &s, &v, r, g, b);
			s = (int)(s - diff);
			//s = (int)(k * SMeanVal_L + (1 - k)*SMeanVal_F);
			//s = (int)(SMeanVal_F + 0.5);
			HSV2RGB(&r, &g, &b, h, s, v);
			frontImgOverlap.at<Vec3b>(y, x)[2] = r;
			frontImgOverlap.at<Vec3b>(y, x)[1] = g;
			frontImgOverlap.at<Vec3b>(y, x)[0] = b;
		}
	}

	imshow("dstLeftImg", leftImgOverlap);
	imshow("dstFrontImg", frontImgOverlap);
	waitKey(0);

	destroyAllWindows();

	//delete[]R_F; R_F = NULL;
	//delete[]G_F; G_F = NULL;
	//delete[]B_F; B_F = NULL;

	//delete[]R_L; R_L = NULL;
	//delete[]G_L; G_L = NULL;
	//delete[]B_L; B_L = NULL;

	return 0;

}

/*
 @brief 对整个图像求取S的平均值，然后使用平均值进行替换补偿
 */
int ImgSatCompenTest3()
{
	string strLeftImg = "\\\\Sh-ws-share-01\\Development\\Data\\AVM_Joint_PC_Demo_Images\\fig-11-29\\Y0_0.bmp";
	string strFrontImg = "\\\\Sh-ws-share-01\\Development\\Data\\AVM_Joint_PC_Demo_Images\\fig-11-29\\Y0_1.bmp";

	//string strLeftImg = "/Sh-ws-share-01/Development/Data/AVM_Joint_PC_Demo_Images/fig-11-29/Y0_0.bmp";
	//string strFrontImg = "/Sh-ws-share-01/Development/Data/AVM_Joint_PC_Demo_Images/fig-11-29/Y0_1.bmp";

	//读原始图像;
	Mat leftImg = imread(strLeftImg, 1);
	Mat frontImg = imread(strFrontImg, 1);

	//从原始图像中裁剪重合区域;
	Rect OverlapRect(0, 0, 320, 320);

	Mat leftImgOverlap = leftImg.clone();
	Mat frontImgOverlap = frontImg.clone();

	namedWindow("srcLeftImg", 0);
	namedWindow("srcFrontImg", 0);

	namedWindow("dstLeftImg", 0);
	namedWindow("dstFrontImg", 0);

	imshow("srcLeftImg", leftImgOverlap);
	imshow("srcFrontImg", frontImgOverlap);

	float meanValF = 0;
	float meanValL = 0;
	int numL = 0;
	int numF = 0;

	for (int y = 0; y < leftImg.rows; y++)
	{
		for (int x = 0; x < leftImg.cols; x++)
		{
			uchar r_F = frontImgOverlap.at<Vec3b>(y, x)[2];
			uchar g_F = frontImgOverlap.at<Vec3b>(y, x)[1];
			uchar b_F = frontImgOverlap.at<Vec3b>(y, x)[0];

			uchar r_L = leftImgOverlap.at<Vec3b>(y, x)[2];
			uchar g_L = leftImgOverlap.at<Vec3b>(y, x)[1];
			uchar b_L = leftImgOverlap.at<Vec3b>(y, x)[0];

			int h_F = 0, s_F = 0, v_F = 0;
			int h_L = 0, s_L = 0, v_L = 0;

			if (r_F != 0 || g_F != 0 || b_F != 0)
			{
				numF++;
				RGB2HSV(&h_F, &s_F, &v_F, r_F, g_F, b_F);
				//meanValF += s_F; //更改S
				meanValF += v_F; //更改V
			}

			if (r_L != 0 || g_L != 0 || b_L != 0)
			{
				numL++;
				RGB2HSV(&h_L, &s_L, &v_L, r_L, g_L, b_L);
				//meanValL += s_L; //更改S
				meanValL += v_L; //更改V
			}
		}
	}

	meanValF /= numF;
	meanValL /= numL;

	for (int y = 0; y < frontImg.rows; y++)
	{

		for (int x = 0; x < frontImg.cols; x++)
		{
			uchar r_F = frontImgOverlap.at<Vec3b>(y, x)[2];
			uchar g_F = frontImgOverlap.at<Vec3b>(y, x)[1];
			uchar b_F = frontImgOverlap.at<Vec3b>(y, x)[0];

			uchar r_L = leftImgOverlap.at<Vec3b>(y, x)[2];
			uchar g_L = leftImgOverlap.at<Vec3b>(y, x)[1];
			uchar b_L = leftImgOverlap.at<Vec3b>(y, x)[0];

			int h_F = 0, s_F = 0, v_F = 0;
			int h_L = 0, s_L = 0, v_L = 0;

			if (r_F != 0 || g_F != 0 || b_F != 0)
			{
				//numF++;
				RGB2HSV(&h_F, &s_F, &v_F, r_F, g_F, b_F);
				//s_F = meanValF;
				v_F = meanValF;
				HSV2RGB(&r_F, &g_F, &b_F, h_F, s_F, v_F);
				frontImgOverlap.at<Vec3b>(y, x)[2] = r_F;
				frontImgOverlap.at<Vec3b>(y, x)[1] = g_F;
				frontImgOverlap.at<Vec3b>(y, x)[0] = b_F;
				//meanValF += s_F;
			}

			if (r_L != 0 || g_L != 0 || b_L != 0)
			{
				//numL++;
				RGB2HSV(&h_L, &s_L, &v_L, r_L, g_L, b_L);
				//s_L = meanValL;
				v_L = meanValL;
				HSV2RGB(&r_L, &g_L, &b_L, h_L, s_L, v_L);
				//meanValL += s_L;
				leftImgOverlap.at<Vec3b>(y, x)[2] = r_L;
				leftImgOverlap.at<Vec3b>(y, x)[1] = g_L;
				leftImgOverlap.at<Vec3b>(y, x)[0] = b_L;

			}
		}
	}

	imshow("dstLeftImg", leftImgOverlap);
	imshow("dstFrontImg", frontImgOverlap);
	waitKey(0);

	destroyAllWindows();

	return 0;
}