#include <math.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct {
  int m;
  int n;
  double **mat;
} mat_t;

/*申请动态内存*/
mat_t new_mat(int m, int n) 
{
  int i, j;
  mat_t M;
  M.mat = (double **)malloc(m * sizeof(double *));
  for (i = 0; i < m; i++) 
  {
    M.mat[i] = (double *)malloc(n * sizeof(double));
  }
  M.m = m;
  M.n = n;
  for (i = 0; i < (M.m); i++)
  {
    for (j = 0; j < (M.n); j++) 
    {
      M.mat[i][j] = 0;
    } //初始化零矩阵//
  }
  return M;
}
/*矩阵释放*/
void free_mat(mat_t M) {
  int i;
  for (i = 0; i < (M.m); i++) {
    free(M.mat[i]);
  }
  free(M.mat);
}
/*矩阵打印*/
void mat_print(mat_t M) {
  int i, j;
  for (i = 0; i < (M.m); i++) {
    for (j = 0; j < (M.n); j++) {
      printf("%lf ", M.mat[i][j]);
    }
    printf("\n");
  }
}
/*矩阵复制*/
void mat_copy(mat_t A,mat_t B) //将已知矩阵B复制到一初始A中
{
  int i, j;
  for (i = 0; i < (B.m); i++) 
  {
    for (j = 0; j < (B.n); j++) 
    { A.mat[i][j] = B.mat[i][j];}
  }
}

// 矩阵克隆
mat_t new_mat_clone(mat_t A)//创建一新矩阵A*，并将A复制到A*中
{
    mat_t Ac=new_mat(A.m,A.n);
    mat_copy(Ac,A);
}
/*矩阵赋值*/
void mat_set_all(mat_t M) 
{
  int i, j;
  printf("请输入矩阵:\n"); 
  for (i = 0; i < (M.m); i++) 
  {
    for (j = 0; j < (M.n); j++) 
    {scanf("%lf", &(M.mat[i][j]));}
  }
 }
/*矩阵减法*/
double mat_sub(mat_t R, mat_t A, mat_t B) 
{
 if(((A.m)==(B.m))&&((R.m)==(A.m))&&((A.n)==(B.n))&&((R.n)==(A.n))) 
  {
    int i, j;
    for (i = 0; i < (A.m); i++) {
      for (j = 0; j < (A.n); j++) {
        (R.mat[i][j]) = (A.mat[i][j]) - (B.mat[i][j]);
      }
    }
  } else 
  {exit(-1);}
}
/*矩阵加法*/
double mat_add(mat_t R, mat_t A, mat_t B) 
{
  if( (((R.m)==(A.m))&&((R.m))==(B.m))&&((R.n)==(A.n))&&((R.n)==(B.n)) )
  {
    int i, j;
    for (i = 0; i < (A.m); i++) 
    {
      for (j = 0; j < (A.n); j++) 
      {(R.mat[i][j]) = (A.mat[i][j]) + (B.mat[i][j]); }
    }
  }
  else {exit(-1);}
}

/*矩阵乘法匹配判断*/
int mat_can_mul(mat_t A, mat_t B) 
{
  if ((A.n) == (B.m)) 
  {return 0;}
  else
  {
    printf("矩阵格式不匹配，无法相乘\n");
    exit(-1);
  }
}

/*矩阵乘法*/
void mat_mul(mat_t R, mat_t A, mat_t B) //矩阵乘法
{
  mat_can_mul(A, B);
  // TODO: R.m == A.m && R.n == B.n
 if(R.m==A.m&&R.n==B.n)
 {
  int i, j, k;
  int sum = 0;
  for (i = 0; i < (A.m); i++) 
  {
    for (j = 0; j < (B.n); j++) 
    {                                //(AB)ij=Sigma[a(ik)*b（kj)]
      sum = 0;
      for (k = 0; k < (A.n); k++) 
      {sum += (A.mat[i][k]) * (B.mat[k][j]);}
      R.mat[i][j] = sum;
    }
  }
}
else
{
    printf("储存结果矩阵格式不匹配");
    exit(-1);
}
}
/*矩阵转置*/
void mat_transpose(mat_t R, mat_t M)
{
  int i, j;
  // TODO: check size match
 if(R.m==M.n&&R.n==M.m)
 {
  for (i = 0; i < (M.n); i++)
  {
    for (j = 0; j < (M.m); j++) 
    {R.mat[i][j] = M.mat[j][i];}
  }}
 else
{
    printf("储存结果矩阵格式不匹配\n");
    exit(-1);
}}
/*矩阵数乘*/
// TODO: 数乘 （mat_t R,mat_t A,double a); R = a*A
void mat_scaler(mat_t R,mat_t A,double a)
{
    int i,j;
    for(i=0;i<(A.m);i++)
    {
        for(j=0;j<(A.n);j++)
        {R.mat[i][j]=a*A.mat[i][j];}
    }}
void mat_set(mat_t M, int i, int j, double c) //将矩阵的（i,j)元素替换为c
{
  M.mat[i - 1][j - 1] = c;
}
double mat_get(mat_t M, int i, int j) 
{
    if(0<i<=M.m&&0<j<=M.n)
  // TODO: make sure i \in [1,M.m], j \in [1,M.n]
    {return M.mat[i - 1][j - 1];}
    else
    {
        printf("提取数据越界");
        exit(-1);
    }
}
//数量矩阵
mat_t mat_scalar(mat_t M,double a)
{
   int m=M.m,n=M.n;
    if(m==n)
    {
        int i,j;
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(i==j)
                {M.mat[i][j]=a;}
                else
                {M.mat[i][j]=0;}
        }}}
return M;
}
mat_t new_mat_vec(int m)//定义m个元素的列向量M
{ return new_mat(m,1); }
mat_t new_mat_row(int n)//定义n个元素的行向量M
{return new_mat(1,n);}
mat_t new_vec_get(mat_t M,int j)//获取矩阵M的第j(从0开始）列为新的列向量
{
    int i;
    mat_t R=new_mat_vec(M.m);
    for(i=0;i<M.m;i++)
    {R.mat[i][0]=M.mat[i][j];}
    return R;
}
//将矩阵M所有元素取绝对值
mat_t mat_absolute(mat_t M)
{
    int i,j; 
    for(i=0;i<M.m;i++)
    {
        for(j=0;j<M.n;j++)
       {
        if(M.mat[i][j]<0)
        {M.mat[i][j]=0-M.mat[i][j];}
       }
     }
    return M;
}
//求矩阵M的最大元素max
double max(mat_t M)
{
    mat_t A=new_vec_get(M,0);//将列向量A初始化为矩阵M的第一（0）列
    int i,j;
    double max;
    for(i=0;i<M.m;i++) //求出矩阵M每一行的最大值构成列向量A
    {
        for(j=0;j<M.n-1;j++)
        {
            if(M.mat[i][j]<M.mat[i][j+1])
            {A.mat[i][0]=M.mat[i][j+1];}
            else
            {A.mat[i][0]=M.mat[i][j+1];}
        }
    }
    for(i=0;i<M.m-1;i++)//求出列向量A中的最大值max
    {
    if(A.mat[i][0]<A.mat[i+1][0])
    {max=A.mat[i+1][0];}
    else
    {max=A.mat[i][0];}
    }
    return max;
}
// TODO: mat_t
// 向量1范数
double norm_vec_1(mat_t M) 
{
  int i;
  double sum = 0;
  mat_t A=mat_absolute(M);
  for(i=0;i<A.m;i++)
  {sum += A.mat[0][i];}
  return sum;
}
//向量无穷范数
double norm_vec_infinite(mat_t M)
{
  int i;
  mat_t A=mat_absolute(M);
  return max(A);
}
// 向量2范数
double norm_vec_2(mat_t M) 
{
  int i;
  double sum = 0;
  for (i = 0; i<(M.n); i++) 
  {sum += (M.mat[i][0]) * (M.mat[i][0]);}
  return sqrt(sum);
}
// F范数
double norm_F(mat_t M) 
{
  int i, j;
  double sum = 0;
  for (i = 0; i < (M.m); i++) 
  {
    for (j = 0; j < (M.n); j++) 
    {sum += (M.mat[i][j]) * (M.mat[i][j]);}
  }
  return sqrt(sum);
}
//矩阵1范数
double norm_mat_1(mat_t M)
{
    int i,j;
    mat_t R=mat_absolute(M);
    mat_t A=new_mat_row(M.n);
    for(j=0;j<R.n;j++)
    {
        for(i=0;i<R.m;i++)
        {A.mat[0][j]+=R.mat[i][j];}
    }
return max(A);
}
//矩阵无穷范数
double norm_mat_infinite(mat_t M)
{
    mat_t A=new_mat(M.m,M.n);
    mat_transpose(A,M);
    return norm_mat_1(A);
}

