#include<stdio.h>
#include<stdlib.h>
#inlclude"mat.c"
void Gauss(mat_t U.vec_t b,int m,int n)
{
    if(m!=n)
    {
        printf("输入矩阵非方阵！\n");
        exit(-1);
    }
    else
 {
    int i,j;
    mat_t M;
    mat_t G;
    mat_t I;
    mat_new I;
    scalar(&I,m.n,1);
    mat_new(&G,(U.m),(U.n));
    mat_new(&M,(U.m),(U.n)+1);
    for(i=0;i<(U.m);i++)
    {
        for(j=0;j<(U.n);j++)
        {
            M.mat[i][j]=U.mat[i][j];
        }
      M.mat[i][(U.n)+1]=b[i];
    }
    mat_t M';
    mat_new(&M',(U,m),(U.n));
    ,
