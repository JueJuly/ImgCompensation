#ifndef _COMMON_H
#define _COMMON_H

#if 1

#include <windows.h>
//#undef min

HANDLE handle = GetStdHandle(STD_OUTPUT_HANDLE);

#define Rcout SetConsoleTextAttribute( handle, 12 );
#define Gcout SetConsoleTextAttribute( handle, 10 );
#define Bcout SetConsoleTextAttribute( handle, 9 );
#define Ccout SetConsoleTextAttribute( handle, 11 );
#define Mcout SetConsoleTextAttribute( handle, 13 );
#define Ycout SetConsoleTextAttribute( handle, 14 );
#define Wcout SetConsoleTextAttribute( handle, 15 );
#define Kcout SetConsoleTextAttribute( handle, 0 );
#define GRcout SetConsoleTextAttribute( handle, 8 );

#endif

#endif