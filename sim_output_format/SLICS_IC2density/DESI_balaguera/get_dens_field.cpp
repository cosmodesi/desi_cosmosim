/* Code to read IC's from the SLIC simulation. Andrés Balaguera-Antolínez IAC 2020 */


// This option is only available if FFTW is propery working

//#define USE_SUPERSAMPLING

#include <ctime>
#include <cmath>
#include <cctype>
#include <string>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <cassert>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <omp.h>
#include <unistd.h>
using namespace std;
#define USE_OMP
#define MAX_MAS_DEG static_cast<int>(3)
#define MAX_MAS_DEG_PSC static_cast<int>(4)
#define RED     "\033[31m"      /* Red */
#define RESET   "\033[0m"
#define ULONG unsigned long
#define real_prec double
#define num_1 static_cast<int>(1)
#define num_0_5 static_cast<real_prec>(0.5)
#define Number_headers static_cast<int>(1)  // Number of properties of IC (x,y,z,vx,vy,vz)
#define Nres_sim static_cast<int>(4)  //Number  of sub-cells (per dim) in which the simulation volume has ben divided
#define Lside_sim static_cast<real_prec>(3072.0)  //Lenght of simulation volume
#define Lside static_cast<real_prec>(505.0)  //Comoving lenght (in Mpc/h) of cosmological volume
#define column_x static_cast<int>(0)
#define column_y static_cast<int>(1)
#define column_z static_cast<int>(2)
#define Ncolumns static_cast<int>(6)  // Number of properties of IC (x,y,z,vx,vy,vz)
// **********************************************************************************************************************
inline ULONG index_3d(int i, int j, int k, int Nj, int Nk)
{
  return static_cast<ULONG>(k)+static_cast<ULONG>(Nk*j)+static_cast<ULONG>(Nk*Nj*i);
}
// **********************************************************************************************************************
void index2coords(int N, int index, int  &XG, int &YG, int &ZG )
{// Get the index in each direction : F-order cells (N,index, X,Y,Z). -Corder (N,index, Z,Y,X)
  XG=index % N;
  index = static_cast<int>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  YG=index % N;
  index = static_cast<int>(static_cast<real_prec>(index)/static_cast<real_prec>(N));
  ZG = index;
}


// **********************************************************************************************************************
// **********************************************************************************************************************
// **********************************************************************************************************************
inline real_prec MAS_TSC(real_prec x)
{
  x=fabs(x);
  real_prec ans;
  if(x<0.5)
      ans= (0.75-x*x);
  if(x>=0.5 && x<1.5)
  {
    real_prec r = 1.5 - x;
    r *= r;
    ans= (0.5*r);
  }
    else
      ans= 0;
    return ans;
}



// **********************************************************************************************************************
inline real_prec MAS_PCS(real_prec x)
{
   real_prec ans;
   real_prec y=fabs(x);
   if(y>=0 && y<1)
     ans=(1./6.)*(4.0- 6.0*x*x + 3.*pow(y,3));
   else if(y>=1. && y<2.)
     ans= (1./6.)*pow(2.0-y, 3);
   else
    ans=0.;
  return ans;
}
// **********************************************************************************************************************
void getDensity_PCS(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec deltax, real_prec deltay, real_prec deltaz, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, vector<real_prec>&delta)
{

#ifndef  _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif
  real_prec rdelta_x=1./deltax;
  real_prec rdelta_y=1./deltay;
  real_prec rdelta_z=1./deltaz;
  ULONG count_previous=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count_previous)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count_previous+=delta[i];
  ULONG counter_a=0;
  for(ULONG ip=0;ip<xp.size();++ip)
    {
      real_prec xpos=xp[ip];
      real_prec ypos=yp[ip];
      real_prec zpos=zp[ip];
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      int xc = static_cast<ULONG>(floor((xpos-min1)*rdelta_x)); // indices of the cell of the particle
      int yc = static_cast<ULONG>(floor((ypos-min2)*rdelta_y));
      int zc = static_cast<ULONG>(floor((zpos-min3)*rdelta_z));
      xc = static_cast<ULONG>(fmod(real_prec(xc),real_prec(N1)));
      yc = static_cast<ULONG>(fmod(real_prec(yc),real_prec(N2)));
      zc = static_cast<ULONG>(fmod(real_prec(zc),real_prec(N3)));
      real_prec xx  = deltax*(static_cast<real_prec>(xc)+0.5);
      real_prec yy  = deltay*(static_cast<real_prec>(yc)+0.5);
      real_prec zz  = deltaz*(static_cast<real_prec>(zc)+0.5);
      real_prec xxf = deltax*(static_cast<real_prec>(xc)+1.5);
      real_prec yyf = deltay*(static_cast<real_prec>(yc)+1.5);
      real_prec zzf = deltaz*(static_cast<real_prec>(zc)+1.5);
      real_prec xxff = deltax*(static_cast<real_prec>(xc)+2.5);
      real_prec yyff = deltay*(static_cast<real_prec>(yc)+2.5);
      real_prec zzff = deltaz*(static_cast<real_prec>(zc)+2.5);
      real_prec xxbb = deltax*(static_cast<real_prec>(xc)-1.5);
      real_prec yybb = deltay*(static_cast<real_prec>(yc)-1.5);
      real_prec zzbb = deltaz*(static_cast<real_prec>(zc)-1.5);
      real_prec xxb = deltax*(static_cast<real_prec>(xc)-0.5);
      real_prec yyb = deltay*(static_cast<real_prec>(yc)-0.5);
      real_prec zzb = deltaz*(static_cast<real_prec>(zc)-0.5);
      int Xb=(xc==0 ? N1: xc);
      int Yb=(yc==0 ? N2: yc);
      int Zb=(zc==0 ? N3: zc);
      int Xf=(xc==N1-1 ? -1: xc);
      int Yf=(yc==N2-1 ? -1: yc);
      int Zf=(zc==N3-1 ? -1: zc);
      int Xbb=(xc==0 ? N1: xc);
      int Ybb=(yc==0 ? N2: yc);
      int Zbb=(zc==0 ? N3: zc);
      if(xc!=0)
        Xbb=(xc==1 ? N1+1: xc);
      if(yc!=0)
       Ybb=(yc==1 ? N2+1: yc);
      if(zc!=0)
       Zbb=(zc==1 ? N3+1: zc);
      int Xff=(xc==N1-1 ? -1: xc);
      int Yff=(yc==N2-1 ? -1: yc);
      int Zff=(zc==N3-1 ? -1: zc);
      if(xc!=N1-1)
         Xff=(xc==N1-2 ? -2: xc);
      if(yc!=N2-1)
        Yff=(yc==N2-2 ? -2: yc);
      if(zc!=N3-1)
        Zff=(zc==N3-2 ? -2: zc);
      vector<int> i_idx = {Xbb-2, Xb-1, xc, Xf+1, Xff+2};
      vector<int> j_idx = {Ybb-2, Yb-1, yc, Yf+1, Xff+2};
      vector<int> k_idx = {Zbb-2, Zb-1, zc, Zf+1, Zff+2};
      vector<real_prec> MAS_xx=
        {
          MAS_PCS((xxbb - xpos)*rdelta_x),
          MAS_PCS((xxb  - xpos)*rdelta_x),
          MAS_PCS((xx   - xpos)*rdelta_x),
          MAS_PCS((xxf  - xpos)*rdelta_x),
          MAS_PCS((xxff - xpos)*rdelta_x)
      };
      vector<real_prec> MAS_yy =
        {
          MAS_PCS((yybb - ypos)*rdelta_y),
          MAS_PCS((yyb  - ypos)*rdelta_y),
          MAS_PCS((yy   - ypos)*rdelta_y),
          MAS_PCS((yyf  - ypos)*rdelta_y),
          MAS_PCS((yyff - ypos)*rdelta_y)
        };
      vector<real_prec> MAS_zz =
        {
          MAS_PCS((zzbb - zpos)*rdelta_z),
          MAS_PCS((zzb  - zpos)*rdelta_z),
          MAS_PCS((zz   - zpos)*rdelta_z),
          MAS_PCS((zzf  - zpos)*rdelta_z),
          MAS_PCS((zzff - zpos)*rdelta_z)
        };
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
      for(int ih=0;ih<MAX_MAS_DEG_PSC;++ih)
        for(int jh=0;jh<MAX_MAS_DEG_PSC;++jh)
          for(int kh=0;kh<MAX_MAS_DEG_PSC;++kh)
#ifdef USE_OMP
#pragma omp atomic update
#endif
            delta[index_3d(i_idx[ih],j_idx[jh], k_idx[kh], N3, N2) ] += MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
    counter_a++;
  }
  ULONG count=0;
#pragma omp parallel for reduction(+:count)
  for(ULONG i=0;i<delta.size();++i)
    count+=delta[i];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid = "<<count<<endl;
}
// **********************************************************************************************************************
void getDensity_TSC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec deltax, real_prec deltay, real_prec deltaz, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, vector<real_prec>&delta)
{

#ifndef  _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif

  real_prec rdelta_x=1./deltax;
  real_prec rdelta_y=1./deltay;
  real_prec rdelta_z=1./deltaz;
  ULONG count_previous=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count_previous)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count_previous+=delta[i];

  ULONG counter_a=0;
  for(ULONG ip=0;ip<xp.size();++ip)
    {
      real_prec xpos=xp[ip];
      real_prec ypos=yp[ip];
      real_prec zpos=zp[ip];
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      int xc = static_cast<ULONG>(floor((xpos-min1)*rdelta_x)); // indices of the cell of the particle
      int yc = static_cast<ULONG>(floor((ypos-min2)*rdelta_y));
      int zc = static_cast<ULONG>(floor((zpos-min3)*rdelta_z));
      xc = static_cast<ULONG>(fmod(real_prec(xc),real_prec(N1)));
      yc = static_cast<ULONG>(fmod(real_prec(yc),real_prec(N2)));
      zc = static_cast<ULONG>(fmod(real_prec(zc),real_prec(N3)));
      real_prec xx  = deltax*(static_cast<real_prec>(xc)+0.5);
      real_prec yy  = deltay*(static_cast<real_prec>(yc)+0.5);
      real_prec zz  = deltaz*(static_cast<real_prec>(zc)+0.5);
      real_prec xxf = deltax*(static_cast<real_prec>(xc)+1.5);
      real_prec yyf = deltay*(static_cast<real_prec>(yc)+1.5);
      real_prec zzf = deltaz*(static_cast<real_prec>(zc)+1.5);
      real_prec xxb = deltax*(static_cast<real_prec>(xc)-0.5);
      real_prec yyb = deltay*(static_cast<real_prec>(yc)-0.5);
      real_prec zzb = deltaz*(static_cast<real_prec>(zc)-0.5);
      int Xb=(xc==0 ? N1: xc);
      int Yb=(yc==0 ? N2: yc);
      int Zb=(zc==0 ? N3: zc);
      int Xf=(xc==N1-1 ? -1: xc);
      int Yf=(yc==N2-1 ? -1: yc);
      int Zf=(zc==N3-1 ? -1: zc);
      vector<int>i_idx{Xb-1,xc,Xf+1};
      vector<int>j_idx{Yb-1,yc,Yf+1};
      vector<int>k_idx{Zb-1,zc,Zf+1};
      vector<real_prec>MAS_xx{MAS_TSC((xxb- xpos)*rdelta_x),MAS_TSC((xx - xpos)*rdelta_x),MAS_TSC((xxf- xpos)*rdelta_x)};
      vector<real_prec>MAS_yy{MAS_TSC((yyb- ypos)*rdelta_y),MAS_TSC((yy - ypos)*rdelta_y),MAS_TSC((yyf- ypos)*rdelta_y)};
      vector<real_prec>MAS_zz{MAS_TSC((zzb- zpos)*rdelta_z),MAS_TSC((zz - zpos)*rdelta_z),MAS_TSC((zzf- zpos)*rdelta_z)};
#ifdef USE_OMP
#pragma omp parallel for collapse(3)
#endif
      for(int ih=0;ih<MAX_MAS_DEG;++ih)
        for(int jh=0;jh<MAX_MAS_DEG;++jh)
          for(int kh=0;kh<MAX_MAS_DEG;++kh)
#ifdef USE_OMP
#pragma omp atomic update
#endif
            delta[index_3d(i_idx[ih],j_idx[jh], k_idx[kh], N3, N2)] += MAS_xx[ih]*MAS_yy[jh]*MAS_zz[kh];
    counter_a++;
  }
  ULONG count=0;
#pragma omp parallel for reduction(+:count)
  for(ULONG i=0;i<delta.size();++i)
    count+=delta[i];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid = "<<count<<endl;
}
// ***************************************************************************************************************************************************************

void getDensity_CIC(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec> &zp, vector<real_prec>&delta)
{
// This function has been originally provided by F-S Kitaura. fskitaura@iac.es
#ifndef  _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#ifdef USE_OMP
#pragma omp parallel for
#endif
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;
#endif


  ULONG count_previous=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count_previous)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count_previous+=delta[i];

  ULONG counter_a=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:counter_a)
#endif
  for (ULONG ig=0; ig<xp.size(); ++ig)
    {   
      real_prec xpos=xp[ig]-num_0_5*d1;
      real_prec ypos=yp[ig]-num_0_5*d2;
      real_prec zpos=zp[ig]-num_0_5*d3;
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      ULONG i = static_cast<ULONG>(floor((xpos-min1)/d1));
      ULONG j = static_cast<ULONG>(floor((ypos-min2)/d2));
      ULONG k = static_cast<ULONG>(floor((zpos-min3)/d3));
      i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
      k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
      ULONG ii = static_cast<ULONG>(fmod(real_prec(i+1),real_prec(N1)));
      ULONG jj = static_cast<ULONG>(fmod(real_prec(j+1),real_prec(N2)));
      ULONG kk = static_cast<ULONG>(fmod(real_prec(k+1),real_prec(N3)));
      real_prec xc = static_cast<real_prec>(i);
      real_prec yc = static_cast<real_prec>(j);
      real_prec zc = static_cast<real_prec>(k);
      real_prec dx = (xpos-min1)/d1-xc;
      real_prec dy = (ypos-min2)/d2-yc;
      real_prec dz = (zpos-min3)/d3-zc;
      real_prec tx = num_1-dx;
      real_prec ty = num_1-dy;
      real_prec tz = num_1-dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(j+N2*i)]    += tx*ty*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(j+N2*ii)]   += dx*ty*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(jj+N2*i)]   += tx*dy*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(j+N2*i)]   += tx*ty*dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[k+N3*(jj+N2*ii)]  += dx*dy*tz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(j+N2*ii)]  += dx*ty*dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(jj+N2*i)]  += tx*dy*dz;
#ifdef USE_OMP
#pragma omp atomic update
#endif
      delta[kk+N3*(jj+N2*ii)] += dx*dy*dz;
      counter_a++;
    }
  ULONG count=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count)
#endif
  for(ULONG i=0;i<delta.size();++i)
    count+=delta[i];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid = "<<count<<endl;
}
// **********************************************************************************************************************
// **********************************************************************************************************************
void getDensity_NGP(ULONG N1, ULONG N2, ULONG N3,real_prec L1, real_prec L2, real_prec L3, real_prec d1, real_prec d2, real_prec d3, real_prec min1, real_prec min2, real_prec min3,const vector<real_prec>&xp, const vector<real_prec>&yp, const vector<real_prec>&zp, vector<real_prec>&delta)
{
#ifndef _GET_INTERPOLATED_FIELDS_FROM_SEVERAL_BIN_FILES_
#pragma omp parallel for
  for (ULONG i=0;i<delta.size(); i++)
    delta[i]= 0.;  //-1 if we want to calculate overdensity
#endif

  ULONG count_previous=0;
#pragma omp parallel for reduction(+:count_previous)
  for(ULONG ip=0;ip<delta.size();++ip)
    count_previous+=delta[ip];
  ULONG counter_a=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:counter_a)
#endif
  for (ULONG ig=0; ig<xp.size(); ++ig)
    {
      real_prec xpos=xp[ig];
      real_prec ypos=yp[ig];
      real_prec zpos=zp[ig];
      if(xpos< 0)xpos+=L1;
      if(ypos< 0)ypos+=L2;
      if(zpos< 0)zpos+=L3;
      if(xpos>=L1)xpos-=L1;
      if(ypos>=L2)ypos-=L2;
      if(zpos>=L3)zpos-=L3;
      ULONG i = static_cast<ULONG>(floor((xpos-min1)/d1)); // indices of the cell of the particle
      ULONG j = static_cast<ULONG>(floor((ypos-min2)/d2));
      ULONG k = static_cast<ULONG>(floor((zpos-min3)/d3));
      i = static_cast<ULONG>(fmod(real_prec(i),real_prec(N1)));
      j = static_cast<ULONG>(fmod(real_prec(j),real_prec(N2)));
      k = static_cast<ULONG>(fmod(real_prec(k),real_prec(N3)));
#pragma omp atomic update
      delta[k+N3*j+N3*N2*i]++;
      counter_a++;
    }

  ULONG count=0;
#ifdef USE_OMP
#pragma omp parallel for reduction(+:count)
#endif
  for(ULONG ip=0;ip<delta.size();++ip)
    count+=delta[ip];
  cout<<"Number of objects = "<<counter_a<<endl;
  cout<<"Number of objects assigned to the grid for the current sub-volume = "<<count-count_previous<<endl;
  cout<<"Cumulative number of objects assigned to the grid  = "<<count<<endl;
}
// ***************************************************************************************************************************************************************
void read_bin_file(string input_file,  int Nbytes_header, int Nbytes_data, int &Nparts, vector<float>&prop)
{
  cout<<"Reading input file "<<input_file<<endl;
  ifstream input;
  input.open(input_file.c_str(), ios::binary| ios::in);
  if(!input)
    {
      cout<<RED<<"File not found. Code exits here."<<RESET<<endl;
      exit(0);
    }
  input.read((char *)&Nparts,Nbytes_header);
  while(!input.eof())
    {
      float data_cat;
      input.read((char *)&data_cat,Nbytes_data);
      prop.push_back(data_cat);
    }
  input.close();
}



// ***************************************************************************************************************************************************************
void write_to_binary(string FNAME, vector<real_prec>&out)
{
  cout<<"Writting in binary file "<<FNAME<<endl;
  ofstream outf(FNAME,ios::out | ios::binary);
  for(ULONG i=0;i<out.size();++i)
    {
      float new_out=static_cast<float>(out[i]);
      outf.write((char *)&new_out,sizeof(float));
    }
  outf.close();
}


// ***************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************
// ********************************************MAIN MAIN MAIN MAIN ***********************************************************************************************
// ***************************************************************************************************************************************************************
// ***************************************************************************************************************************************************************
int main(int argc, char *argv[]){

  if(argc==1){
    cout<<RED<<"Error: code expects 4 parameters. Only "<<argc-1<<" provided."<<RESET<<endl;

    exit(1);
  }

  char temp;

  while((temp =  getopt(argc, argv, "hr:")) != -1)
    {
      if(temp=='h')
	{
          cout<<endl;
          cout<<"This code reads binary files with SLICS IC and intepolate to a mesh."<<endl;
          cout<<"balaguera@iac.es 2020"<<endl;
          cout<<endl;
          cout<<"Input parameter are:"<<endl;
          cout<<"1 -> <Path>: Path to the IC folder, e.g, /global/cscratch1/sd/jharno/DESI/IC/"<<endl;
          cout<<"2 -> <nIC>: Number x of the SLAC realization in path <Path> with the prefix LOS (e.g 1002, 986 etc )"<<endl;
          cout<<"3 -> <MAS>: Interpolation scheme: 0 (NGP), 1 (CIC), 2 (TSC), 3 (PCS)"<<endl;
          cout<<"4 -> <Nres>: Number of cells per dimention (e.g. 256 for a 256³ mesh)"<<endl;
          cout<<"5 -> <Outputdir>: Name of output directory "<<endl;
          cout<<"The code reads the 64 files located in the path <Path>/LOS<nIC>/xv##.ic"<<endl;
          cout<<"where ## is from 0 to 63 (64 sub-volumes used in the generation of IC)"<<endl;
          cout<<"The output is written in a single-precision binary file (C order)"<<endl;
          cout<<"at Outputdir/SLICS_IC_LOS<nIC>_Nres<Nres>_MAS<MAS>.dat"<<endl;
          cout<<endl;
      }

      else if(temp=='r')
        {
	  if(argc<5){
            cout<<RED<<"Error: code expects 5 parameters. Only "<<argc-1<<" provided."<<RESET<<endl;
	    exit(0);
	  }
          string path_ic = argv[2];
          int nIC  =atoi(argv[3]);
          int MAS = atoi(argv[4]);   // Mass assignment scheme: 0 for NGP, 1 for CIC
          int Nres= atoi(argv[5]);  // Number of cell/dimention
          string output_path = argv[6];
          if(Nres<=0)
	    {
	      cout<<RED<<"Error: Mesh resolution (Nres) must be > 0."<<RESET<<endl;
	      exit(0);
	    }
	  if(MAS>3)
	    {
	      cout<<RED<<"Error: requested mass assignment scheme not valid."<<RESET<<endl;
	      cout<<RED<<"Available: 0 (NGP), 1 (CIC), 2 (TSC), 3(PCS)"<<RESET<<endl;
	      exit(0);
	    }

	  string sMAS;
	  switch(MAS){case(0):sMAS="NGP";break; case(1):sMAS="CIC";break; case(2):sMAS="TSC";break; case(3):sMAS="PCS";break;}
	  cout<<"Getting IC on a mesh with size "<<Nres<<"³ using Mass Assignment Scheme "<<sMAS<<endl;
          string path=path_ic+"LOS"+to_string(nIC)+"/";

	  ULONG N_sub_boxes = pow(Nres_sim,3);
	  real_prec Lside_sub = Lside_sim/static_cast<real_prec>(Nres_sim);
	  real_prec delta_sim = static_cast<real_prec>(Lside/Lside_sim);
	  cout<<"Length of sub boxes in sim-units = "<<Lside_sub<<endl;
	  int Nbytes_header=sizeof(int);
	  int Nbytes_data=sizeof(float);


	  ULONG Ngrid=Nres*Nres*Nres;
	  vector<real_prec>density_field(Ngrid,0);


	  int i_sub_box=0;

	  for(i_sub_box=0; i_sub_box<N_sub_boxes;  ++ i_sub_box)
	    {   
	      string bin_file = path+"xv"+to_string(i_sub_box)+".ic";
	      vector<float>prop;
	      int Nparts;
	      cout<<endl;
	      read_bin_file(bin_file,Nbytes_header, Nbytes_data, Nparts, prop);
	      vector<real_prec>xcoord(Nparts,0), ycoord(Nparts,0), zcoord(Nparts,0);
	      cout<<"Number of particles "<<Nparts<<endl;
	      int ixb,iyb,izb;
	      index2coords(Nres_sim, i_sub_box,ixb,iyb,izb);
	      cout<<"Sub-box index = "<<i_sub_box<<". Sub-box coords: "<<ixb<<" "<<iyb<<" "<<izb<<endl;
#pragma omp parallel for
	      for(ULONG i=0;i< Nparts;++i)
		{
		  xcoord[i]=(static_cast<real_prec>(prop[column_x+i*Ncolumns])+static_cast<real_prec>(ixb*Lside_sub))*delta_sim;
		  ycoord[i]=(static_cast<real_prec>(prop[column_y+i*Ncolumns])+static_cast<real_prec>(iyb*Lside_sub))*delta_sim;
		  zcoord[i]=(static_cast<real_prec>(prop[column_z+i*Ncolumns])+static_cast<real_prec>(izb*Lside_sub))*delta_sim;
		}
	      prop.clear();prop.shrink_to_fit();

	      real_prec delta=Lside/static_cast<real_prec>(Nres);

	      cout<<"Interpolating ... "<<endl;
	      if(0==MAS)
		getDensity_NGP(Nres,Nres,Nres,Lside,Lside,Lside,delta,delta,delta,0,0,0,xcoord,ycoord,zcoord,density_field);
	      else if (1==MAS)
		getDensity_CIC(Nres,Nres,Nres,Lside,Lside,Lside,delta,delta,delta,0,0,0,xcoord,ycoord,zcoord,density_field);
	      else if (2==MAS)
		getDensity_TSC(Nres,Nres,Nres,Lside,Lside,Lside,delta,delta,delta,0,0,0,xcoord,ycoord,zcoord,density_field);
	      else if (3==MAS)
		getDensity_PCS(Nres,Nres,Nres,Lside,Lside,Lside,delta,delta,delta,0,0,0,xcoord,ycoord,zcoord,density_field);
	    }

          string output_den=output_path+"/SLICS_IC_LOS"+to_string(nIC)+"_Nres"+to_string(Nres)+"_MAS"+to_string(MAS)+".dat";
	  write_to_binary(output_den,density_field);
	}
    }
}
















