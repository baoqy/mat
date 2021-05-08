#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
void Gauss(mat_t A,vec_t b,int m,int n)
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
    mat_t Mx;
    mat_t Gx;
    mat_t L;
    mat_t Lx;
    mat_t Ly;
    mat_t U;
    mat_new(&U,m,n);
    mat_new(&Ly,m,n);
    mat_new(&L,m,n);
    mat_new(&Lx,m,n);
    mat_new(&Mx,m,n+1);
    mat_new(&Gx,m,n);
    mat_new(&I,m,n);
    mat_scalar(&I,m,n,1);
    mat_new(&G,m,n);
    mat_new(&M,m,(n+1));
    for(i=0;i<(A.m);i++)          //计算矩阵U的增广矩阵M
    {
        for(j=0;j<(A.n);j++)
        {
            M.mat[i][j]=A.mat[i][j];
        }
      M.mat[i][A.n]=b.vec[i];
    }
    printf("方程的增广矩阵为:\n");
    mat_print(M);
    vec_t r;
    vec_new(&r,(A.n));
    vec_t v;
    vec_new(&v,(A.m));
    int k;
    for(k=0;k<((A.n)-1);k++)
    {                                 //取第k个元素为1，其余均为0的行向量r
      for(i=0;i<(A.n);i++)
      {
        if(i==k)
        {
            r.vec[i]=1;
        }
        else
        {
            r.vec[i]=0;
        }
      }
                                     //取增广矩阵第j列的下半部分v  
        for(i=0;i<(A.n);i++)
        {
            if(i>k)
            {
                v.vec[i]=M.mat[i][k];
            }
            else
            {
                v.vec[i]=0;
            }
        }                        
    for(i=0;i<(A.m);i++)
    {
        for(j=0;j<(A.n);j++)
        {
           L.mat[i][j]= G.mat[i][j]=(v.vec[i])*(r.vec[j])/(M.mat[k][k]);
            
        }
}
  mat_sub(Gx,I,G);                   //计算LU分解的L
  mat_add(Lx,I,G);
  if(k>0)
{
    mat_mul(L,Ly,Lx);
}
if(k<((U.n)-1))
{
    mat_mul(Mx,Gx,M);
    printf("上三角增广矩阵为:\n");
    for(i=0;i<(A.m);i++)                 //由最后的上三角增广矩阵得到LU分解的U
    {
        for(j=0;j<(A.n);j++)
        {
            U.mat[i][j]=Mx.mat[i][j];
        }
    }
    mat_print(Mx);
}
    mat_clone(M,Mx);  //将计算中间过程的增广矩阵和算子迭代
    mat_clone(G,Gx);
    mat_clone(Ly,Lx);
}
printf("所用LU分解为L=:\n");
mat_print(L);
printf("所用LU分解U=:\n");
mat_print(U);
}
}
 int main()
{
    int m,n;
    printf("输入矩阵行数m与列数n：\n");
    scanf("%d%d",&m,&n);
    mat_t U;
    mat_new(&U,m,n);
    mat_set(U,m,n);
    vec_t b;
    vec_new(&b,m);
    vec_set(b,m);
    Gauss(U,b,m,n);
    mat_free(U);
    free(b.vec);
    return 0;
 }
