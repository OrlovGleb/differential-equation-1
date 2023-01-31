#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <math.h>
const static double  pi=3.14159265358979323846;


//разность строк
void raznost(int m, double *A, int q, int p, double z) 
{

    int j;
    for (j = 0; j < m; j++) {
        A[m * q + j] = A[m * q + j] - (A[m * p + j]) * z;
    }
}
//алгоритм Жордана-Гаусса
void f(int n, int m, double *A, double *X)
{
    int i;
    int j;
    int k;
    double b;
    double max;
    double q;
    for (k = 0; k < n; k++) {
        max = fabs(A[m * k + k]);

        for (i = k; i < n; i++) {
            if (A[m * i + k] > max) {
                max = A[m * i + k];
                for (j = 0; j < m; j++) {
                    q = A[m * i + j];
                    A[m * i + j] = A[m * k + j];
                    A[m * k + j] = q;
                }
            }
        }

        b = A[m * k + k]; //записываем в b коэффицент при первом символе строки
      
        for (j = 0; j < m; j++) {
            A[m * k + j] = (A[m * k + j]) / (b);
        }
        for (i = 0; i < n; i++) {
            if (i != k) {
                raznost(m, A, i, k, A[m * i + k]);
            } //вычитаем строки
        }
    }
    for (i = 0; i < n; i++) {
        X[i] = A[m * i + m - 1]; //заполнение массива с ответами
    }
}


int main(void)
{

    int j;
    int i;
    int n;
    int a1;
    int m;
    int t;
    double err=0;
    double h;
    scanf("%d", &n);
    h=pow(n,-1);

    m = n + 1;
    double *A;
    A = (double *)malloc((m) * (m+1) * sizeof(double)); //общая матрица
    double X[m+1]; //матрица для ответа
    double X1[m+1]; //матрица иксов
    double R[m+1];  //вспомогательный массив
    double Q[m+1];  //вспомогательный массив
    double P[m+1];  //вспомогательный массив
    double D[m+1];  //вспомогательный массив
    double D1[m+1];  //вспомогательный массив для аппроксимации
    scanf("%d", &t); //t-флаг, отвечающий за выбор решения
    scanf("%d", &a1);
    //зануляем матрицу
    for (i = 0; i < n+1; i++) {
        for (j = 0; j < n+2; j++) {
            A[i*(n+2)+j]=0;
        }     
    }
    //заполняем матрицу иксов
    for (i = 0; i < n+1; i++) {
        X1[i] = i*h;
        
    }
    //заполняем вспомогательные матрицы
    for (i = 0; i < n+1; i++) {
        Q[i] = 1+3*(h*h*cos(X1[i]))/2;
        P[i] = h*h*h-3-2*(h*h*cos(X1[i]));
        R[i] = 3+(h*h*cos(X1[i]))/2;
        if (t==0)
        D[i] = h*h*h*exp(X1[i]*X1[i]);
        else D[i] = h*h*h*(a1*sin(pi*X1[i]/2)*sin(pi*X1[i]/2)
        +(a1*pi/2)*sin(pi*X1[i])*cos(X1[i])-(a1*pi*pi*pi/2)*sin(pi*X1[i]));      
    }
    //задаем значения на границах
    A[0]=1;
    A[n+1]=0;
    A[n+2]=-3;
    A[n+3]=4;
    A[n+4]=-1;
    A[n+2+n+1]=0;
    A[n*(n+2)+n]=1;
    A[n*(n+2)+n+1]=a1;
    //заполняем оставшуюся часть матрицы
    for (i=2; i<n; i++){
        A[(n+2)*i+i-2]=-1;
        A[(n+2)*i+i-1]=R[i];
        A[(n+2)*i+i]=P[i];
        A[(n+2)*i+i+1]=Q[i];
        A[(n+2)*i+n+1]=D[i]; 
    }

    //решаем слу
    f(n+1, m+1, A, X);
    //находим ошибку
    if (t==1){
     for (i = 0; i < n+1; i++){
      if(fabs(X[i]-(a1*sin(pi*X1[i]/2)*sin(pi*X1[i]/2)))>err) 
      err=fabs(X[i]-(a1*sin(pi*X1[i]/2)*sin(pi*X1[i]/2))) ;
    }
    }
    printf("%lf",err);

    //открываем файл и записываем в него ответ        
    FILE *f;
    if (t==0)
    f=fopen("xxx.txt", "w");
    else f=fopen("xxx1.txt", "w");
    for (i = 0; i < n+1; i++) {
        fprintf(f,"%.3lf", X1[i]);
        fprintf(f," ");
        fprintf(f,"%.3lf", X[i]);
        fprintf(f,"\n");
    }
    //теперь решаем систему для другого равенства, 
    //чтобы аппроксимировать результат, имея частное решение
    

    free(A);
    return 0;
}