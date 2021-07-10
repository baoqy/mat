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
void mat_print(mat_t M) 
{
  int i, j;
  for (i = 0; i < (M.m); i++)
  {
    for (j = 0; j < (M.n); j++) 
    {
      printf("%lf\t", M.mat[i][j]);
    }
    printf("\n\n");
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
    return Ac;
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
  double sum;
  for (i = 0; i < (A.m); i++) 
  {
    for (j = 0; j < (B.n); j++) 
    {    
        sum=0;                           //(AB)ij=Sigma[a(ik)*b（kj)]
      for (k = 0; k < (A.n); k++) 
      {sum =sum+ (A.mat[i][k]) * (B.mat[k][j]);}
      R.mat[i][j]=sum;
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
    printf("转置储存结果矩阵格式不匹配\n");
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
  M.mat[i][j] = c;
}
//获取矩阵M的（i,j)元素
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
mat_t mat_scalar(int n,double a)
{
   mat_t A=new_mat(n,n);
        int i,j;
        for(i=0;i<n;i++)
        {
            for(j=0;j<n;j++)
            {
                if(i==j)
                {A.mat[i][j]=a;}
                else
                {A.mat[i][j]=0;}
        }}
return A;
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
mat_t new_row_get(mat_t M,int i)//获取矩阵M的第i(从0开始）行为新的行向量
{
    int j;
    mat_t R=new_mat_row(M.n);
    for(j=0;j<M.n;j++)
    {R.mat[0][j]=M.mat[i][j];}
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
double mat_max(mat_t M)
{
    int i,j,k;
    mat_t A=new_mat_vec((M.m)*(M.n));
    for(i=0;i<M.m;i++)//将就很M转化成列向量A
    {
      for(j=0;j<M.n;j++)
      {
          k=i*(M.n)+j;
          A.mat[k][0]=M.mat[i][j];
      }
    }
    double max=A.mat[0][0];
    for(i=0;i<(M.m)*(M.n)-1;i++)//求出列向量A中的最大值max
    {
      if(max<A.mat[i+1][0])
      {max=A.mat[i+1][0];}
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
  {sum += A.mat[i][0];}
  return sum;
}
//向量无穷范数
double norm_vec_infinite(mat_t M)
{
  mat_t A=mat_absolute(M);
  return mat_max(M);
}
// 向量2范数
double norm_vec_2(mat_t M) 
{
  int i;
  double sum = 0;
  for (i = 0; i<(M.m); i++) 
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
   return mat_max(A);
}
//矩阵无穷范数
double norm_mat_infinite(mat_t M)
{
    mat_t A=new_mat(M.m,M.n);
    mat_transpose(A,M);
    return norm_mat_1(A);
}
//解上三角方程
mat_t mat_U_solve(mat_t A,mat_t b)//Ax=b,A为上三角方阵且对角线元素均非0
{
   int m=A.m,n=A.n;
   int j,k;
   double sum;
   mat_t x=new_mat_vec(m);
   x.mat[m-1][0]=(b.mat[m-1][0])/(A.mat[m-1][n-1]);
       for(j=(m-2);j>=0;j--)
       { 
           sum=0;
       for(k=(m-2);k>=j;k--)
       { sum+=(x.mat[k+1][0])*A.mat[j][k+1];}
       x.mat[j][0]=(b.mat[j][0]-sum)/(A.mat[j][j]);
       }
   return x;
}
//解下三角方程
mat_t mat_L_solve(mat_t A,mat_t b)//A为下三角矩阵且对角线元素非0
{
    int m=A.m,n=A.n;
    mat_t At=new_mat(m,n);
    mat_t bt=new_mat_vec(m);
    mat_t Ax=new_mat(m,n);
    int i,j;
    for(i=0;i<m;i++)
        for(j=0;j<n;j++)
        {
            Ax.mat[i][j]=A.mat[i][n-1-j];
        }

    for(i=0;i<m;i++)
    {
        bt.mat[i][0]=b.mat[n-1-i][0];
        for(j=0;j<n;j++)
        {
            At.mat[i][j]=Ax.mat[n-1-i][j];
        }
    }
    mat_t y=(mat_U_solve(At,bt));
    mat_t x=new_mat_vec(m);
    for(i=0;i<m;i++)
{
    x.mat[i][0]=y.mat[m-1-i][0];
}
return x;;
}
//向量b的Givens矩阵
mat_t Givens(mat_t b,int i,int j)//b为m行列向量
{
    int n=b.m;
    mat_t T=mat_scalar(n,1);
    double c,s;
    c=(b.mat[i][0])/(sqrt(pow(b.mat[i][0],2)+pow(b.mat[j][0],2)));
    s=(b.mat[j][0])/(sqrt(pow(b.mat[i][0],2)+pow(b.mat[j][0],2)));
    mat_set(T,i,i,c);
    mat_set(T,i,j,s);
    mat_set(T,j,i,-s);
    mat_set(T,j,j,c);
    return T;
}
//求y=Hx的Householder矩阵H，其中x与y的2范数相等且x≠y
mat_t Householder(mat_t x,mat_t y)
{
    double a;
    int m=x.m;
    mat_t I=mat_scalar(m,1);
    mat_t M=new_mat_vec(m);
    mat_sub(M,x,y);
    a=pow(norm_vec_2(M),2);
    mat_t H=new_mat(m,m);
    mat_t K=new_mat(m,m);
    mat_t T=new_mat_row(m);
    mat_transpose(T,M);
    mat_mul(K,M,T);
    mat_t K1=new_mat(m,m);
    mat_scaler(K1,K,2/a);
    mat_sub(H,I,K1);
    return H;
}
//向量内积（数量积）
double dot_product(mat_t v1,mat_t v2)
{
    int i;
    int m=v1.m;
    double sum=0;
    for(i=0;i<m;i++)
    {
        sum+=v1.mat[i][0]*v2.mat[i][0];
    }
    return sum;
}
//用以在计算过程中判断精度是否达到要求（a），若达到要求则停止计算并且返回此时的误差b;
double solution_judge(mat_t A,mat_t x,mat_t b,double a)
{
    int m=A.m,n=x.n;
    int i,j;
    mat_t R=new_mat(m,n);
    mat_t R1=new_mat(m,n);
    mat_mul(R,A,x);
    mat_sub(R1,R,b);
    double c=norm_vec_2(R1);
    if(c<a)
    {
        return c;
    }
    else
    {
       ;
       
    }
}
/*int main()
{
    int m,n;
    printf("行数与列数:\n");
    scanf("%d%d",&m,&n);
    mat_t A=new_mat(m,n);
    mat_set_all(A);
    mat_t x=new_mat_vec(m);
    mat_set_all(x);
    mat_t b=new_mat_vec(m);
    mat_set_all(b);
   solution_judge(A,x,b);
  free_mat(A);
    free_mat(x);
    free_mat(b);
    return 0;
}*/
