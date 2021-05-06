
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
typedef struct
{
    int m;
    int n;
    double **mat;
}mat_t;
/*申请动态内存*/
void mat_new(mat_t *M,int m,int n)
{
    int i,j;
    M->mat=(double**)malloc(m*sizeof(double*));
        for(i=0;i<m;i++)
        {
            M->mat[i]=(double*)malloc(n*sizeof(double));
        }
            M->m=m;
            M->n=n;
        for(i=0;i<(M->m);i++)
        {
            for(j=0;j<(M->n);j++)
            {
                M->mat[i][j]=0;
            }                             //初始化零矩阵//
}
}
/*矩阵赋值*/
void mat_set(mat_t M,int m,int n)
{
   int i,j;
   int w;
   printf("输入矩阵是/否为数量矩阵(1/2)?:\n");
   scanf("%d",&w);
   if(w==2)
   {
       printf("输入矩阵所有元素是/否均相等(1/2)?:\n");
       scanf("%d",&w);
     if(w==2)
       {
           printf("请输入矩阵:\n");                 //非数量矩阵且并非所有元素均相等矩阵输入
             for(i=0;i<(M.m);i++)
               {
                  for(j=0;j<(M.n);j++)
                    {
                      scanf("%lf",&(M.mat[i][j]));
                    }
                }
        }
       else                                              //所有元素均相同矩阵输入
       {
           double h;
           printf("矩阵所有元素均为:\n");
           scanf("%lf",&h);
           for(i=0;i<(M.m);i++)
           {
               for(j=0;j<(M.n);j++)
               {
                   M.mat[i][j]=h;
               }
           }
       }
   }
   else                                                  //数量矩阵输入
   {
      double k;
       printf("输入数量矩阵kI,k=:\n");
       scanf("%lf",&k);
       for(i=0;i<(M.m);i++)
       {
           for(j=0;j<(M.n);j++)
           {
               if(i==j)
               {
                   M.mat[i][j]=k;
               }
               else
                   M.mat[i][j]=0;
           }
       }
    }
}
/*矩阵减法*/
void mat_sub(mat_t R,mat_t A,mat_t B)
{
        if((A.m)==(B.m)&&(A.n)==(B.n))
                {
                    int i,j;
                    R.m=A.m;
                    R.n=A.n;
                    for(i=0;i<(A.m);i++)
                    {
                        for(j=0;j<(A.n);j++)
                        { 
                            (R.mat[i][j])=(A.mat[i][j])-(B.mat[i][j]);
                        }
                    }
                }
            else
                        exit(-1);
}
/*矩阵加法*/
void mat_add(mat_t R,mat_t A,mat_t B)
{
        if((A.m)==(B.m)&&(A.n)==(B.n))
                {
                    int i,j;
                    R.m=A.m;
                    R.n=A.n;
                    for(i=0;i<(A.m);i++)
                    {
                        for(j=0;j<(A.n);j++)
                        {
                            (R.mat[i][j])=(A.mat[i][j])+(B.mat[i][j]);
                        }
                    }
                }
            else
                exit(-1);
}
/*矩阵乘法*/
void mat_mul(mat_t R,mat_t A,mat_t B)//矩阵乘法
{
    if((A.n)==(B.m))
    {
        int i,j,k;
        int sum=0;
        for(i=0;i<(A.m);i++)            
        {
            for(j=0;j<(B.n);j++)                  //(AB)ij=Sigma[a(ik)*b（kj)]
            {
                for(k=0;k<(A.n);k++)
                {
                    sum+=(A.mat[i][k])*(B.mat[k][j]);
                    R.mat[i][j]=sum;
                }
                sum=0;                     //R的每个元素计算完后要使sum归零，否则sum继续保持为上一轮计算的值  
            }
        }
    }
        else
        {
            printf("矩阵格式不匹配，无法相乘\n");
            exit(-1);
        }
}
/*矩阵乘法匹配判断*/
int mat_mul_judge(mat_t A,mat_t B)
{
        if((A.n)==(B.m))
        {
            return 0;
        }
        else
        {
            printf("矩阵格式不匹配，无法相乘\n");
            exit(-1);   
        }
}
/*矩阵转置*/
void mat_transpose(mat_t R,mat_t M)
{
    int i,j;
    R.m=M.n;
    R.n=M.m;
    for(i=0;i<(M.n);i++)
    {
        for(j=0;j<(M.m);j++)
        {
            R.mat[i][j]=M.mat[j][i];
        }
    }
}
/*矩阵定标输出*/
double mat_scaler(mat_t M,int i,int j)//输出矩阵的（i,j)元素
{
    int m,n;
    m=i-1;
    n=j-1;
    printf("所求元素为%lf\n",M.mat[m][n]);
}
/*矩阵元素替换*/
void mat_replace(mat_t M,int i,int j,double c)//将矩阵的（i,j)元素替换为c
{
        M.mat[i-1][j-1]=c;
}
//向量定义
typedef struct 
{
    int n;
    double *vec;
}vec_t;
void vec_new(vec_t *M,int n)
{
    int i;
    M->vec=(double*)malloc(n*sizeof(double));
    M->n=n;
    for(i=0;i<(M->n);i++)
    {
        M->vec[i]=0;
    }
}
//向量赋值
void vec_set(vec_t M,int n)
{
    int i;
    printf("输入列向量:\n");
    for(i=0;i<(M.n);i++)
    {
        scanf("%lf",&(M.vec[i]));
    }
}
//向量打印
void vec_print(vec_t M)      
{
    int i;
    printf("列向量为:\n");
    for(i=0;i<(M.n);i++)
    {
        printf("%lf ",M.vec[i]);
    }
    printf("\n");
}
//输出向量第i个元素
double vec_element(vec_t M,int i) 
{
    printf("列向量的第%d个元素为:\n",i);
    printf("%lf\n",M.vec[i]);
}
//替换向量第i个元素
void vec_replace(vec_t M,int i,double c)
{
    M.vec[i]=c;
}
/*数量矩阵输入*/
void mat_scalar(mat_t *K,int m,int n,int k)
{
    if(m==n)
    {
        int i,j,k;
        for(i=0;i<(K->n);i++)
        {
            for(j=0;j<(K->n);j++)
                K->mat[i][j]=k;
        }
    }
    else
        exit(-1);
}
//1范数
double norm_1(vec_t M)
{
    int i;
    double sum=0;
    for(i=0;i<(M.n);i++)
    {
       if((M.vec[i])<0)
        {
            M.vec[i]=0-M.vec[i];
        }
        sum+=M.vec[i];
    }
    return sum;
}
//无穷范数
double norm_infinite(vec_t M)
{
    int i;
    for(i=0;i<((M.n)-1);i++)
    {
        if((M.vec[i])<0)
        {
            M.vec[i]=-M.vec[i];
        }
    }
    double max=M.vec[0];
    for(i=0;i<(M.n);i++)
    {
        if((M.vec[i])<(M.vec[i+1]))
        {
            max=M.vec[i+1];
        }
    }
    return max;
}
//2范数
double norm_2(vec_t M)
{
    int i;
    double sum=0;
    for(i=0;i<(M.n);i++)
    {
        sum+=(M.vec[i])*(M.vec[i]);
    }
return sqrt(sum);
}
//F范数
double norm_F(mat_t M)
{
    int i,j;
    double sum=0;
    for(i=0;i<(M.m);i++)
    {
        for(j=0;j<(M.n);j++)
        {
            sum+=(M.mat[i][j])*(M.mat[i][j]);
        }
    }
    return sqrt(sum);
}
