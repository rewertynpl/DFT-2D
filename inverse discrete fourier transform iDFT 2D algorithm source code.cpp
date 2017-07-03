//inverse discrete fourier transform iDFT 2D algorithm source code
//odwrotna dyskretna transformacja fouriera idft 2d algorytm kod źródłowy
//author marcin matysek (r)ewertyn.PL

 #include <iostream>
#include "conio.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>


using namespace std;

//complex number method:
void fun_inverse_discrete_fourier_transform_DFT_2D_method1(int N,double table[][12][12]);
void fun_discrete_fourier_transform_DFT_2D_method1(int N,double table[][12][12]);

//complex number method:
void fun_inverse_discrete_fourier_transform_DFT_2D_method2(int N,double table[][12][12]);
void fun_discrete_fourier_transform_DFT_2D_method2(int N,double table[][12][12]);

//other method not complex number signal need to be normal number
void fun_inverse_discrete_fourier_transform_DFT_2D_method3(int N,double table[][12][12]);
void fun_discrete_fourier_transform_DFT_2D_method3(int N,double table[][12][12]);


//inverse_method1 works only witch discrete _method1
//inverse_method2 works only witch discrete _method2
//inverse_method3 works only witch discrete _method3

static double diffclock(clock_t clock2,clock_t clock1)
{

    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}

int main()
{
    //zał N=okres sygnału w tablicy tab[] wtedy rodzielczość = 1 Hz

    int N=12;
    double time2;
    //tab[0]][N][N]=re
    //tab[1]][N][N]=im
    double tab[2][12][12]={{
    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-1.630526192,-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.735499212,-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192},

    {-0.923879533,0.70664666,0.996551596,0.923879533,1.659378744,1.369473808,-0.923879533,
    -2.29335334,-0.735499212,0.923879533,-0.072672064,-1.630526192}
    }};

//tab[0]][N][N]=re
//tab[1]][N][N]=im

    tab[1][0][0]=5.1234;//im number
    tab[1][1][0]=9.1234;//im number
    tab[1][2][1]=-3.1234;//im number

    cout<<" "<<tab[1][0][0]<<endl;
    cout<<" "<<tab[1][1][0]<<endl;
    cout<<" "<<tab[1][2][1]<<endl;

cout<<"signal g(x):"<<endl<<endl;
    for(int i=0;i<N;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[0][i][j]*100)/100<<"  ";
    }
    cout<<endl;
    }
    cout<<endl;

    clock_t start = clock();
    fun_discrete_fourier_transform_DFT_2D_method2(N,tab);
    time2=diffclock( start, clock() );

    cout<<" time="<<time2/1000<<endl;
    cout<<"transformation 2D re"<<endl<<endl;
    for(int i=0;i<N;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[0][i][j]*100)/100<<"  ";
    }
    cout<<endl;
    }
    cout<<endl;
 cout<<"transformation 2D im"<<endl<<endl;
    for(int i=0;i<N;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[1][i][j]*100)/100<<"  ";
    }
    cout<<endl;
    }
    cout<<endl;
    fun_inverse_discrete_fourier_transform_DFT_2D_method2(N,tab);

 cout<<"inverse transformation 2D"<<endl<<endl;
    for(int i=0;i<N;i++)
    {
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab[0][i][j]*100)/100<<"  ";
    }
    cout<<endl;
    }
    cout<<" "<<round(tab[1][0][0]*100)/100<<endl;
    cout<<" "<<round(tab[1][1][0]*100)/100<<endl;
    cout<<" "<<round(tab[1][2][1]*100)/100<<endl;
    cout<<endl;

    system("pause");
    return 0;
}

void fun_discrete_fourier_transform_DFT_2D_method1(int N,double table[][12][12])
{
     int i=N,j=N,k=N,m=N;
    const double pi=3.141592653589793238462;

    double*** table3 = new double**[4];
    //double*** table = new double**[4];//do 3 wymiarowej
    for( int i = 0; i < 4; ++i)
    {
      //table[i] = new double*[N];//do 3 wymiarowej
      table3[i] = new double*[N];

      for( int j = 0; j < N; ++j)
      {
          table3[i][j]= new double[N];
      }

    }
for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table3[i][j][k]=0;
    }
    }
}
for (int k=0;k<N;k++)
{
    for(int l=0;l<N;l++)
    {
for (int m=0;m<N;m++)
{
    for(int n=0;n<N;n++)
    {
        //nr 1 complex number method: other combinations are possible but you need other combination in inverse method too
          table3[0][k][l]=table3[0][k][l]+table[0][m][n]*cos((k*m+n*l)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]-table[0][m][n]*sin((k*m+n*l)*2*pi/(float)N);
        table3[0][k][l]=table3[0][k][l]-table[1][m][n]*sin((k*m+n*l)*2*pi/(float)N)*-1;//im*im
          table3[1][k][l]=table3[1][k][l]+table[1][m][n]*cos((k*m+n*l)*2*pi/(float)N);
    }}}
}

for(int i=0;i<N;i++)
{

    for(int j=0;j<N;j++)
    {
        //nr 1;
      table[0][i][j] =(table3[0][i][j]);
      table[1][i][j] =(table3[1][i][j]);

      //nr 2
      //table[0][i][j] =(table3[0][i][j]+table3[1][i][j]+table3[2][i][j]+table3[3][i][j]);
     // table[1][i][j] = 0;
    }
}

   for( int i = 0; i < 4; ++i)
   {
    for( int j = 0; j < N; ++j)
        {
       delete[] table3[i][j];
        }
       delete[] table3[i];
   }
    delete[] table3;

}
void fun_inverse_discrete_fourier_transform_DFT_2D_method1(int N,double table[][12][12])
{
     int i=N,j=N,k=N,m=N;
    const double pi=3.141592653589793238462;

    double*** table3 = new double**[4];
    //double*** table = new double**[4];//do 3 wymiarowej
    for( int i = 0; i < 4; ++i)
    {
      //table[i] = new double*[N];//do 3 wymiarowej
      table3[i] = new double*[N];

      for( int j = 0; j < N; ++j)
      {
          table3[i][j]= new double[N];
      }

    }
for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table3[i][j][k]=0;
    }
    }
}
for (int k=0;k<N;k++)
{
    for(int l=0;l<N;l++)
    {
for (int m=0;m<N;m++)
{
    for(int n=0;n<N;n++)
    {
        // nr 1 complex number method: other combinations are possible but you need other combination in discrete method too
          table3[0][k][l]=table3[0][k][l]+table[0][m][n]*cos((k*m+n*l)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]+table[0][m][n]*sin((k*m+n*l)*2*pi/(float)N);
        table3[0][k][l]=table3[0][k][l]+table[1][m][n]*sin((k*m+n*l)*2*pi/(float)N)*-1;//im*im
          table3[1][k][l]=table3[1][k][l]+table[1][m][n]*cos((k*m+n*l)*2*pi/(float)N);
    }}}
}

for(int i=0;i<N;i++)
{

    for(int j=0;j<N;j++)
    {
      table[0][i][j] =(table3[0][i][j])/(N*N);
      table[1][i][j] =(table3[1][i][j])/(N*N);
    }
}

   for( int i = 0; i < 4; ++i)
   {
    for( int j = 0; j < N; ++j)
        {
       delete[] table3[i][j];
        }
       delete[] table3[i];
   }
    delete[] table3;

}


void fun_discrete_fourier_transform_DFT_2D_method2(int N,double table[][12][12])
{
     int i=N,j=N,k=N,m=N;
    const double pi=3.141592653589793238462;

    double*** table3 = new double**[4];
    //double*** table = new double**[4];//do 3 wymiarowej
    for( int i = 0; i < 4; ++i)
    {
      //table[i] = new double*[N];//do 3 wymiarowej
      table3[i] = new double*[N];

      for( int j = 0; j < N; ++j)
      {
          table3[i][j]= new double[N];
      }

    }

    double table4[2][12][12]={};


for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table3[i][j][k]=0;
    }
    }
}
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table4[i][j][k]=0;
    }
    }
}

for (int m=0;m<N;m++)
{
    for(int l=0;l<N;l++)
    {
    for(int n=0;n<N;n++)
    {//complex number method:
        //nr 1 complex number method: other combinations are possible but you need other combination in inverse method too
          table4[0][m][l]=table4[0][m][l]+table[0][m][n]*cos((n*l)*2*pi/(float)N);
        table4[1][m][l]=table4[1][m][l]-table[0][m][n]*sin((n*l)*2*pi/(float)N);
        table4[0][m][l]=table4[0][m][l]-table[1][m][n]*sin((n*l)*2*pi/(float)N)*-1;//im*im;
          table4[1][m][l]=table4[1][m][l]+table[1][m][n]*cos((n*l)*2*pi/(float)N);
}}}

 for (int k=0;k<N;k++)
{
    for(int l=0;l<N;l++)
    {
    for(int m=0;m<N;m++)
    {//complex number method:
        table3[0][k][l]=table3[0][k][l]+table4[0][m][l]*cos((k*m)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]-table4[0][m][l]*sin((k*m)*2*pi/(float)N);
        table3[0][k][l]=table3[0][k][l]-table4[1][m][l]*sin((k*m)*2*pi/(float)N)*-1;//im*im;
          table3[1][k][l]=table3[1][k][l]+table4[1][m][l]*cos((k*m)*2*pi/(float)N);
    }}}

for(int i=0;i<N;i++)
{

    for(int j=0;j<N;j++)
    {

      table[0][i][j] =(table3[0][i][j]);
      table[1][i][j] =table3[1][i][j];
    }
}

   for( int i = 0; i < 4; ++i)
   {
    for( int j = 0; j < N; ++j)
        {
       delete[] table3[i][j];
        }
       delete[] table3[i];
   }
    delete[] table3;

}
void fun_inverse_discrete_fourier_transform_DFT_2D_method2(int N,double table[][12][12])
{
     int i=N,j=N,k=N,m=N;
    const double pi=3.141592653589793238462;

    double*** table3 = new double**[4];
    //double*** table = new double**[4];//do 3 wymiarowej
    for( int i = 0; i < 4; ++i)
    {
      //table[i] = new double*[N];//do 3 wymiarowej
      table3[i] = new double*[N];

      for( int j = 0; j < N; ++j)
      {
          table3[i][j]= new double[N];
      }

    }
    double table4[2][12][12]={};


for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table3[i][j][k]=0;
    }
    }
}
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table4[i][j][k]=0;
    }
    }
}

for (int m=0;m<N;m++)
{
    for(int l=0;l<N;l++)
    {
    for(int n=0;n<N;n++)
    {
        //complex number method:
        // nr 1 complex number method: other combinations are possible but you need other combination in discrete method too
          table4[0][m][l]=table4[0][m][l]+table[0][m][n]*cos((n*l)*2*pi/(float)N);
        table4[1][m][l]=table4[1][m][l]+table[0][m][n]*sin((n*l)*2*pi/(float)N);
        table4[0][m][l]=table4[0][m][l]+table[1][m][n]*sin((n*l)*2*pi/(float)N)*-1;//im*im;
          table4[1][m][l]=table4[1][m][l]+table[1][m][n]*cos((n*l)*2*pi/(float)N);
}}}

 for (int k=0;k<N;k++)
{
    for(int l=0;l<N;l++)
    {
    for(int m=0;m<N;m++)
    {
        //complex number method:
        table3[0][k][l]=table3[0][k][l]+table4[0][m][l]*cos((k*m)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]+table4[0][m][l]*sin((k*m)*2*pi/(float)N);
        table3[0][k][l]=table3[0][k][l]+table4[1][m][l]*sin((k*m)*2*pi/(float)N)*-1;//im*im;
          table3[1][k][l]=table3[1][k][l]+table4[1][m][l]*cos((k*m)*2*pi/(float)N);
    }}}

for(int i=0;i<N;i++)
{

    for(int j=0;j<N;j++)
    {
      table[0][i][j] =(table3[0][i][j])/(N*N);
      table[1][i][j] =table3[1][i][j]/(N*N);
    }
}

   for( int i = 0; i < 4; ++i)
   {
    for( int j = 0; j < N; ++j)
        {
       delete[] table3[i][j];
        }
       delete[] table3[i];
   }
    delete[] table3;

}

void fun_discrete_fourier_transform_DFT_2D_method3(int N,double table[][12][12])
{
     int i=N,j=N,k=N,m=N;
    const double pi=3.141592653589793238462;

    double*** table3 = new double**[4];
    //double*** table = new double**[4];//do 3 wymiarowej
    for( int i = 0; i < 4; ++i)
    {
      //table[i] = new double*[N];//do 3 wymiarowej
      table3[i] = new double*[N];

      for( int j = 0; j < N; ++j)
      {
          table3[i][j]= new double[N];
      }

    }

    double table4[2][12][12]={};


for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table3[i][j][k]=0;
    }
    }
}
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table4[i][j][k]=0;
    }
    }
}

for (int m=0;m<N;m++)
{
    for(int l=0;l<N;l++)
    {
    for(int n=0;n<N;n++)
    {

          table4[0][m][l]=table4[0][m][l]+table[0][m][n]*cos((n*l)*2*pi/(float)N);
        table4[1][m][l]=table4[1][m][l]+table[0][m][n]*sin((n*l)*2*pi/(float)N);
        table4[0][m][l]=table4[0][m][l]+table[1][m][n]*sin((n*l)*2*pi/(float)N)*-1;//im*im
          table4[1][m][l]=table4[1][m][l]+table[1][m][n]*cos((n*l)*2*pi/(float)N);
}}}

 for (int k=0;k<N;k++)
{
    for(int l=0;l<N;l++)
    {
    for(int m=0;m<N;m++)
    {
        table3[0][k][l]=table3[0][k][l]+table4[0][l][m]*cos((k*m)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]+table4[0][l][m]*sin((k*m)*2*pi/(float)N);
        table3[0][k][l]=table3[0][k][l]+table4[1][l][m]*sin((k*m)*2*pi/(float)N)*-1;//im*im
          table3[1][k][l]=table3[1][k][l]+table4[1][l][m]*cos((k*m)*2*pi/(float)N);
    }}}

for(int i=0;i<N;i++)
{

    for(int j=0;j<N;j++)
    {
        //nr 1;
      //table[0][i][j] =(table3[0][i][j]+table3[3][i][j]);
      //table[1][i][j] =(table3[1][i][j]+table3[2][i][j]);

      //nr 2
      table[0][i][j] =(table3[0][i][j]+table3[1][i][j]+table3[2][i][j]+table3[3][i][j]);
      table[1][i][j] =0;
    }
}

   for( int i = 0; i < 4; ++i)
   {
    for( int j = 0; j < N; ++j)
        {
       delete[] table3[i][j];
        }
       delete[] table3[i];
   }
    delete[] table3;

}
void fun_inverse_discrete_fourier_transform_DFT_2D_method3(int N,double table[][12][12])
{
     int i=N,j=N,k=N,m=N;
    const double pi=3.141592653589793238462;

    double*** table3 = new double**[4];
    //double*** table = new double**[4];//do 3 wymiarowej
    for( int i = 0; i < 4; ++i)
    {
      //table[i] = new double*[N];//do 3 wymiarowej
      table3[i] = new double*[N];

      for( int j = 0; j < N; ++j)
      {
          table3[i][j]= new double[N];
      }

    }
    double table4[2][12][12]={};


for (int i=0;i<4;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table3[i][j][k]=0;
    }
    }
}
for (int i=0;i<2;i++)
{
    for(int j=0;j<N;j++)
    {    for(int k=0;k<N;k++)
    {
       table4[i][j][k]=0;
    }
    }
}

for (int m=0;m<N;m++)
{
    for(int l=0;l<N;l++)
    {
    for(int n=0;n<N;n++)
    {

          table4[0][m][l]=table4[0][m][l]+table[0][m][n]*cos((n*l)*2*pi/(float)N);
        table4[1][m][l]=table4[1][m][l]+table[0][m][n]*sin((n*l)*2*pi/(float)N);
}}}

 for (int k=0;k<N;k++)
{
    for(int l=0;l<N;l++)
    {
    for(int m=0;m<N;m++)
    {
        table3[0][k][l]=table3[0][k][l]+table4[0][l][m]*cos((k*m)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]+table4[0][l][m]*sin((k*m)*2*pi/(float)N);
          table3[0][k][l]=table3[0][k][l]+table4[1][l][m]*cos((k*m)*2*pi/(float)N);
        table3[1][k][l]=table3[1][k][l]+table4[1][l][m]*sin((k*m)*2*pi/(float)N);
    }}}

for(int i=0;i<N;i++)
{

    for(int j=0;j<N;j++)
    {
      table[0][i][j] =(table3[0][i][j]+table3[1][i][j]+table3[2][i][j]+table3[3][i][j])/(N*N);
    }
}

   for( int i = 0; i < 4; ++i)
   {
    for( int j = 0; j < N; ++j)
        {
       delete[] table3[i][j];
        }
       delete[] table3[i];
   }
    delete[] table3;

}



//http://inverse-fast-fourier-transform-fft.blogspot.com/

