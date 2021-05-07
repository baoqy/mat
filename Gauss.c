#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
void Gauss(mat_t U,vec_t b,int m,int n)
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
    mat_new(&Mx,m,n+1);
    mat_new(&Gx,m,n);
    mat_new(&I,m,n);
    mat_scalar(&I,m,n,1);
    mat_new(&G,m,n);
    mat_new(&M,m,(n+1));
    for(i=0;i<(U.m);i++)          //计算矩阵U的增广矩阵M
    {
        for(j=0;j<(U.n);j++)
        {
            M.mat[i][j]=U.mat[i][j];
        }
      M.mat[i][U.n]=b.vec[i];
    }
    printf("增广矩阵为:\n");
    mat_print(M);
    vec_t r;
    vec_new(&r,(U.n));
    vec_t v;
    vec_new(&v,(U.m));
    int k;
    for(k=0;k<((U.n)-1);k++)
    {                                 //取第k个元素为1，其余均为0的行向量r
      for(i=0;i<(U.n);i++)
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
    printf("行向量为:\n");
    vec_print(r);
                                     //取增广矩阵第j列的下半部分v  
        for(i=0;i<(U.n);i++)
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
  printf("列向量下部分为:\n");
  vec_print(v);
    for(i=0;i<(U.m);i++)
    {
        for(j=0;j<(U.n);j++)
        {
            G.mat[i][j]=(v.vec[i])*(r.vec[j])/(M.mat[k][k]);
        }
}
  printf("列向量与行向量之积为:\n");//计算算子G
  mat_print(G);
  printf("算子为:\n");
  mat_sub(Gx,I,G);
  mat_print(Gx);
if(k<((U.n)-1))
{
    mat_mul(Mx,Gx,M);
    printf("上三角增广矩阵为:\n");
    mat_print(Mx);
}
else if(k==((U.n)-1))
{
    exit(-1);
}
    mat_clone(M,Mx);  //将计算中间过程的增广矩阵和算子迭代
    mat_clone(G,Gx);
}
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
