#include <math.h>
#include <stdio.h>
#include <stdlib.h>
typedef struct {
  int m;
  int n;
  double **mat;
} mat_t;

/*申请动态内存*/ new_ *free_ *mat_t new_mat(int m, int n) {
  int i, j;
  mat_t M;
  M.mat = (double **)malloc(m * sizeof(double *));
  for (i = 0; i < m; i++) {
    M.mat[i] = (double *)malloc(n * sizeof(double));
  }
  M.m = m;
  M.n = n;
  for (i = 0; i < (M.m); i++) {
    for (j = 0; j < (M.n); j++) {
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
/*矩阵克隆*/
// TODO: copy
void mat_copy(mat_t A, mat_t B) //将已知矩阵B复制到一初始A中
{
  int i, j;
  for (i = 0; i < (B.m); i++) {
    for (j = 0; j < (B.n); j++) {
      A.mat[i][j] = B.mat[i][j];
    }
  }
}

// TODO:
mat_t new_mat_clone(mat_t A);
/*矩阵赋值*/
void mat_set_all(mat_t M, int m, int n) {
  int i, j;
  /*   int w;
     printf("输入矩阵是/否为数量矩阵(1/2)?:\n");
     scanf("%d",&w);
     if(w==2)
     {
         printf("输入矩阵所有元素是/否均相等(1/2)?:\n");
         scanf("%d",&w);
       if(w==2)
         {*/
  printf("请输入矩阵:\n"); //非数量矩阵且并非所有元素均相等矩阵输入
  for (i = 0; i < (M.m); i++) {
    for (j = 0; j < (M.n); j++) {
      scanf("%lf", &(M.mat[i][j]));
    }
  }
  /* }
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
}*/
}
/*矩阵减法*/
void mat_sub(mat_t R, mat_t A, mat_t B) {
  if ((A.m) == (B.m) && (A.n) == (B.n)) {
    int i, j;
    R.m = A.m;
    R.n = A.n;
    for (i = 0; i < (A.m); i++) {
      for (j = 0; j < (A.n); j++) {
        (R.mat[i][j]) = (A.mat[i][j]) - (B.mat[i][j]);
      }
    }
  } else {
    exit(-1);
  }
}
/*矩阵加法*/
void mat_add(mat_t R, mat_t A, mat_t B) {
  if ((A.m) == (B.m) && (A.n) == (B.n)) {
    int i, j;
    R.m = A.m;
    R.n = A.n;
    for (i = 0; i < (A.m); i++) {
      for (j = 0; j < (A.n); j++) {
        (R.mat[i][j]) = (A.mat[i][j]) + (B.mat[i][j]);
      }
    }
  } else {
    exit(-1);
  }
}

/*矩阵乘法匹配判断*/
int mat_can_mul(mat_t A, mat_t B) {
  if ((A.n) == (B.m)) {
    return 0;
  } else {
    printf("矩阵格式不匹配，无法相乘\n");
    exit(-1);
  }
}

/*矩阵乘法*/
void mat_mul(mat_t R, mat_t A, mat_t B) //矩阵乘法
{
  mat_can_mul(A, B);
  // TODO: R.m == A.m && R.n == B.n
  int i, j, k;
  int sum = 0;
  for (i = 0; i < (A.m); i++) {
    for (j = 0; j < (B.n); j++) { //(AB)ij=Sigma[a(ik)*b（kj)]
      sum = 0;
      for (k = 0; k < (A.n); k++) {
        sum += (A.mat[i][k]) * (B.mat[k][j]);
      }
      R.mat[i][j] = sum;
    }
  }
}

/*矩阵转置*/
void mat_transpose(mat_t R, mat_t M) {
  int i, j;
  // TODO: check size match
  R.m = M.n;
  R.n = M.m;

  for (i = 0; i < (M.n); i++) {
    for (j = 0; j < (M.m); j++) {
      R.mat[i][j] = M.mat[j][i];
    }
  }
}
/*矩阵定标输出*/
// TODO: 数乘 （mat_t R,mat_t A,double a); R = a*A
double mat_scaler(mat_t M, int i, int j) //输出矩阵的（i,j)元素
{
  int m, n;
  m = i - 1;
  n = j - 1;
  printf("所求元素为%lf\n", M.mat[m][n]);
}
/*矩阵元素替换*/
void mat_set(mat_t M, int i, int j, double c) //将矩阵的（i,j)元素替换为c
{
  M.mat[i - 1][j - 1] = c;
}
double mat_get(mat_t M, int i, int j) {
  // TODO: make sure i \in [1,M.m], j \in [1,M.n]
  return M.mat[i - 1][j - 1];
}

mat_t new_mat_vec(int n) { return mat_new(n, 1); }

// TODO: mat_t
// 1范数
double norm_vec_1(vec_t M) {
  int i;
  double sum = 0;
  for (i = 0; i < (M.n); i++) {
    if ((M.vec[i]) < 0) {
      M.vec[i] = 0 - M.vec[i];
    }
    sum += M.vec[i];
  }
  return sum;
}
//无穷范数
double norm_infinite(vec_t M) {
  int i;
  for (i = 0; i < ((M.n) - 1); i++) {
    if ((M.vec[i]) < 0) {
      M.vec[i] = -M.vec[i];
    }
  }
  double max = M.vec[0];
  for (i = 0; i < (M.n); i++) {
    if ((M.vec[i]) < (M.vec[i + 1])) {
      max = M.vec[i + 1];
    }
  }
  return max;
}
// 2范数
double norm_2(vec_t M) {
  int i;
  double sum = 0;
  for (i = 0; i < (M.n); i++) {
    sum += (M.vec[i]) * (M.vec[i]);
  }
  return sqrt(sum);
}
// F范数
double norm_F(mat_t M) {
  int i, j;
  double sum = 0;
  for (i = 0; i < (M.m); i++) {
    for (j = 0; j < (M.n); j++) {
      sum += (M.mat[i][j]) * (M.mat[i][j]);
    }
  }
  return sqrt(sum);
}
