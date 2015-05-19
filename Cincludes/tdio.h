#ifndef TDHHDR
#define TDHHDR

#include <iostream>
#include <fstream>
#include <iomanip>
//#include <strstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

/*
   Standard output calls and functions for TopDraw plots.
*/

class output {
 protected:
  //  _IO_ostream_withassign file;
//    ostream_withassign file;
    ostream file;
 public:  
    output(void) : file(cout.rdbuf()) {}//{file=cout;}
    output(ofstream &os) : file(os.rdbuf()) {} // {file=os;}
    inline void pData(double x, double y, double dy) {
            file.setf(ios::scientific,ios::floatfield);
            file << setiosflags(ios::uppercase)
                 << setprecision(6)
                 <<"   "<<setw(13)<<x<<"   "<<setw(13)<<y<<"   "
                                    <<setw(13)<<dy<<"\n";}
    inline void pData(double x, double y) {
            file.setf(ios::scientific,ios::floatfield);
            file << setiosflags(ios::uppercase)
                 << setprecision(6)
                 <<"   "<<setw(13)<<x<<"   "<<setw(13)<<y<<"\n";}
    inline void pFunction(double (*fun)(double), const double &min,
                           const double &size, const int &num);
    inline void pFunctionLog(double (*fun)(double), const double &min,
                           const double &size, const int &num);
};


class topdraw : public output {
 private:
    double xl, xu, yl, yu;
    int xlog, ylog, sp;
    char *topTitle, *topCase, *leftTitle, *leftCase, *bottomTitle, *bottomCase;
 public:
    topdraw( ofstream &os, double wl=0.0, double wu=1.0, double zl=0.0, 
              double zu=1.0) : output(os) {
              xl=wl; xu=wu; yl=zl; yu=zu; topTitle=0; topCase=0;
              leftTitle=0; leftCase=0; bottomTitle=0; 
              bottomCase=0; xlog=0; ylog=0; sp=0;}
    topdraw( double wl=0.0, double wu=1.0, double zl=0.0, double zu=1.0) 
                               : output() {
              xl=wl; xu=wu; yl=zl; yu=zu; topTitle=0; topCase=0;
              leftTitle=0; leftCase=0; bottomTitle=0; 
              bottomCase=0; xlog=0; ylog=0; sp=0;}
    topdraw(topdraw &td) {xl=td.xl; xu=td.xu; yl=td.yl; yu=td.yu;
                          topTitle=td.topTitle; topCase=td.topCase;
                          leftTitle=td.leftTitle; leftCase=td.leftCase;
                          bottomTitle=td.bottomTitle; bottomCase=td.bottomCase;
                          xlog=0; ylog=0; sp=0;}    
    void left(char *Val) {leftTitle=Val;}
    void leftC(char *cVal) {leftCase=cVal;}
    void bottom(char *Val) {bottomTitle=Val;}
    void bottomC(char *cVal) {bottomCase=cVal;}
    void top(char *Val) {topTitle=Val;}
    void topC(char *cVal) {topCase=cVal;}
    inline void pString(char *Val) {file<<"\n  "<<Val<<"\n\n";}
    inline void xLog(void) {xlog=1;}
    inline void yLog(void) {ylog=1;}  
    inline void pHeader(void);
                      // Print out the TopDraw header.
    inline void pSolid(void) {file<<"\n  JOIN SOLID\n\n";}
    inline void pDots(void) {file<<"\n  JOIN DOTS\n\n";}
    inline void pDashes(void) {file<<"\n  JOIN DASHES\n\n";}
    inline void pDotdash(void) {file<<"\n  JOIN DOTDASH\n\n";}
    inline void pSymbol(int o) {file<<"  SET SYMBOL "<<o<<"O\n\n";}
    inline void pHistogram(void) {file<<"\n  HISTOGRAM\n\n";}
    inline void pHistogramDots(void) {file<<"\n  HISTOGRAM DOTS\n\n";}
    inline void pHistogramDashes(void) {file<<"\n  HISTOGRAM DASHES\n\n";}
    inline void pHistogramDotdash(void) {file<<"\n  HISTOGRAM DOTDASH\n\n";}
    inline void pHistogramSolid(void) {file<<"\n  HISTOGRAM SOLID\n\n";}
    inline void pPlot(void) {file<<"\n  PLOT\n\n";}
    inline void sidePost(void) {sp=1;}
    inline void pHorizontal(void) {
              pData(xl,0.0); pData(xu,0.0); pSolid();}
    inline void pVertical(void) {
              pData(0.0,yl); pData(0.0,yu); pSolid();}
};    


inline void topdraw::pHeader(void)
{
      file.setf(ios::floatfield);
  //      file.setf(0,ios::floatfield);
      file<<setiosflags(ios::uppercase);
      if (sp == 1 )
          file<<"  SET DEVICE ORIENTATION 3 POSTSCRIPT\n";
      file<<"  SET LIMITS X FROM "<<xl<<" TO "<<xu<<" Y FROM "<<yl<<" TO "
                                  <<yu<<"\n";
      file<<"  SET MODE VECTOR\n";      
      file<<"  SET FONT DUPLEX\n";
      file<<"  SET ORDER X Y DY\n";
      if (xlog == 1) 
	file<<"  SET SCALE X LOG\n";
      if (ylog == 1) 
	file<<"  SET SCALE Y LOG\n";
      if (topTitle != 0) 
          file<<"  TITLE TOP \'"<<topTitle<<"\'\n";      
      if (topCase != 0) 
          file<<"  CASE      \'"<<topCase<<"\'\n";      
      if (leftTitle != 0) 
          file<<"  TITLE LEFT \'"<<leftTitle<<"\'\n";
      if (leftCase != 0) 
          file<<"  CASE       \'"<<leftCase<<"\'\n";
      if (bottomTitle != 0) 
          file<<"  TITLE BOTTOM \'"<<bottomTitle<<"\'\n";
      if (bottomCase != 0) 
          file<<"  CASE         \'"<<bottomCase<<"\'\n";
      file<<"\n";
}


inline void output::pFunction(double (*fun)(double), const double &min,
                           const double &size, const int &num)
{
      double x = min - size;
      for(int i=0;i<=num; i++) {
          x += size;
          pData(x,fun(x));
	}
}
 
inline void output::pFunctionLog(double (*fun)(double), const double &min,
                           const double &size, const int &num)
{
      double x = min/size;
      for(int i=0;i<=num; i++) {
          x *= size;
          pData(x,fun(x));
	}
}
 
#endif
