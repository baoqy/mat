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
//修改矩阵的第i，j元素为a；
mat_t mat_set_one(mat_t M,int i,int j,double a)
{
    M.mat[i][j]=a;
    return M;
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
//输出对A做LU分解的L；
mat_t L_LU(mat_t A) 
{
    int m=A.m,n=A.n;
    int i, j;
    mat_t Ax=new_mat(m,n);
    mat_t G=new_mat(m,n);
    mat_t I=mat_scalar(n,1);
   
    mat_t Gx=new_mat(m,n);
    mat_t L=new_mat(m,n);
    mat_t Lx=new_mat(m,n);
    mat_t Ly=new_mat(m,n);
    mat_t U=new_mat(m,n);
    mat_t R=new_mat(m,n);

    int k;
    for(k=0;k<n-1;k++)
    {//取单位阵I的第K列
       mat_t r=new_mat_row(n);
       r.mat[0][k]=1;
      //取增广矩阵第k列的下半部分v
      mat_t v=new_vec_get(A,k);
      for (i = 0; i <n; i++) 
      {
        if (i <=k) 
        {v.mat[i][0] = 0;}
      }
      mat_mul(R,v,r);
      mat_scaler(G,R,1/A.mat[k][k]);
      mat_sub(Gx,I,G);
     mat_add(Lx,I,G);
      if (k > 0) 
      {mat_mul(L, Ly, Lx);}
       mat_mul(Ax, Gx, A);
      mat_copy(U,Ax);//计算LU分解的U
      
      mat_copy(A,Ax); //将计算中间过程的增广矩阵和算子迭代
      mat_copy(G,Gx);
      mat_copy(Ly,Lx);
    }
    
    return L;
  }
mat_t U_LU(mat_t A) 
{
    int m=A.m,n=A.n;
    int i, j;
    mat_t Ax=new_mat(m,n);
    mat_t G=new_mat(m,n);
    mat_t I=mat_scalar(n,1);
   
    mat_t Gx=new_mat(m,n);
    mat_t L=new_mat(m,n);
    mat_t Lx=new_mat(m,n);
    mat_t Ly=new_mat(m,n);
    mat_t U=new_mat(m,n);
    mat_t R=new_mat(m,n);

    int k;
    for(k=0;k<n-1;k++)
    {//取单位阵I的第K列
       mat_t r=new_mat_row(n);
       r.mat[0][k]=1;
      //取增广矩阵第k列的下半部分v
      mat_t v=new_vec_get(A,k);
      for (i = 0; i <n; i++) 
      {
        if (i <=k) 
        {v.mat[i][0] = 0;}
      }
      mat_mul(R,v,r);
      mat_scaler(G,R,1/A.mat[k][k]);
      mat_sub(Gx,I,G);
     mat_add(Lx,I,G);
      if (k > 0) 
      {mat_mul(L, Ly, Lx);}
       mat_mul(Ax, Gx, A);
      mat_copy(U,Ax);//计算LU分解的U
      
      mat_copy(A,Ax); //将计算中间过程的增广矩阵和算子迭代
      mat_copy(G,Gx);
      mat_copy(Ly,Lx);
    }
    
    return U;
  }


double Hermite(mat_t x,mat_t y,mat_t d_y ,double c)//用Hermite插值公式计算新坐标为c时对应的函数值,d_y存放已知点的导数值
{
double d_Li(mat_t v,int k)//计算Hermite插值公式中用到的Li（xi）的导数,v的元素两两不相等
{
    int i;
    int n=v.m;
    double a=0;
    for(i=0;i<n;i++)
    {
        if(i==k)
        {;}
        else
        {
            a+=1/(v.mat[k][0]-v.mat[i][0]);
        }
    }
    return a;
}
double Li(mat_t v,int k,double c)//计算Hermite（同时也是Lagrange）插值公式的Li(x）(i=k时)
{
    int i;
    int n=v.m;
    double a1=1,a2=1;
    for(i=0;i<n;i++)
    {
        if(i==k)
        {;}
        else
        {
           a1*=(c-v.mat[i][0]);
           a2*=(v.mat[k][0]-v.mat[i][0]);
        }
    }
return a1/a2;
 }
    int n=x.m;
    int i,j,k;
    double a1=0,a2=0;
    for(j=0;j<n;j++)
    {
        a1+=y.mat[j][0]*pow(Li(x,j,c),2)*(1-2*d_Li(x,j)*(c-x.mat[j][0]));
        a2+=d_y.mat[j][0]*pow(Li(x,j,c),2)*(c-x.mat[j][0]);
    }
    return a1+a2;
}


//向量v1和v2分别是点组的x和y构成的向量，求其前k(从0计数,故对于n个元素的向量，k最大为n-1)个元素的差商
double d_quotient(mat_t v1,mat_t v2,double k)
{
    int i,j;
    double a=1,b=0;
    for(j=0;j<=k;j++)
    {
        a=1;         //在上次使用完a之后，对a再次初始化为1，防止上次运行的结果a代入下次循环
        for(i=0;i<=k;i++)
        {
           if(i==j)  //当i==j时，执行空指令，跳过相乘运算
           {
               ;
           }
           else
           {
              a*=(v1.mat[j][0]-v1.mat[i][0]);//当i！=j时，正常运算
           }
        }
 
      b+=(v2.mat[j][0]/a);
     
    }
     return b;
}
 

//Newton法求在新坐标x=c时的函数值
double Newton(mat_t x, mat_t y,double c)
{
    int i,j;
    int n=x.m;
    double c1=1,c2=y.mat[0][0];
    for(i=0;i<n-1;i++)//由于在使用d_quotient()函数求差商时使用了i+1，这里规定i<n-1，防止数据访问越界
  {
        c1=1;          //在上次使用玩c1后，对c1重新初始化为1
        for(j=0;j<=i;j++)
        {
            c1*=(c-x.mat[j][0]);
        }
        c2+=c1*d_quotient(x,y,i+1);
    }
    return c2;
}


//共轭梯度法(要求A正定)
//TODO 精度控制以及迭代次数由用户指定
mat_t C_G(mat_t A,mat_t b)
{
    int m=A.m;
    int i,j,k;
    double a,a1,a2,b1;
    mat_t x=new_mat_vec(m);
    mat_t R=new_mat_vec(m);
    mat_t R1=new_mat_vec(m);
    mat_t r=new_mat_vec(m);
    mat_t rk=new_mat_vec(m);
    mat_t p=new_mat_vec(m);
    mat_t pk=new_mat_vec(m);
    mat_t xk=new_mat_vec(m);

    mat_mul(R,A,x);
    mat_sub(r,b,R);
    mat_copy(p,r);

    for(k=0;k<8;k++)
    {
        a1=dot_product(r,p);
        mat_mul(R,A,p);
        a2=dot_product(R,p);
        a=a1/a2;

        mat_scaler(R,p,a);
        mat_add(xk,x,R);

        mat_mul(R,A,p);
        mat_scaler(R1,R,a);
        mat_sub(rk,r,R1);

        a1=dot_product(rk,R);
        a2=dot_product(p,R);
        b1=-a1/a2;

        mat_scaler(R,p,b1);
        mat_add(pk,rk,R);

        mat_copy(r,rk);
        mat_copy(p,pk);
        mat_copy(x,xk);
        a=solution_judge(A,x,b,0.001);//Ax-b二范数精确度0.001
        if(a<0.001)
        {
            return x;
            break;
        }
    }

}

//TODO 迭代精度控制和迭代次数由用户指定
//Jacobi迭代
mat_t Jacobi(mat_t A,mat_t b)
{
    int m=A.m,n=A.n;
    int i,j,k;
    mat_t X=new_mat_vec(m);
    for(i=0;i<m;i++)
    {
        X.mat[i][0]=b.mat[i][0]/A.mat[i][i];
    }
    mat_t Xk=new_mat_vec(m);
    mat_t v=new_mat_row(m);
    for(k=0;k<=7;k++)
    {
        for(i=0;i<n;i++)
        {
           double sum1=0;
          for(j=0;j<=i-1;j++)
          {
              sum1+=A.mat[i][j]*X.mat[j][0];
          }

          double sum2=0;
          for(j=i+1;j<n;j++)
          {
              sum2+=A.mat[i][j]*X.mat[j][0];
          }

          Xk.mat[i][0]=-(1/A.mat[i][i])*(sum1+sum2-b.mat[i][0]);
         }
        mat_copy(X,Xk);
        mat_transpose(v,X);
      //  printf("第%d次运算结果为:\n",k+1);
       // mat_print(v);
    }
    return v;
}


//用QR分解解线性方程组Ax=b（用Givens矩阵实现）
mat_t QR(mat_t A,mat_t b)//n阶方阵A
{
    int n=A.n,m=A.m;
    int i,j;
    mat_t v=new_mat_vec(m);
    mat_t T=new_mat(n,n);
    mat_t bt=new_mat_vec(m);
    mat_t Ax=new_mat(n,n);
     for(j=0;j<n-1;j++)
     {
         for(i=j;i<n-1;i++)
         {
           mat_t v=new_vec_get(A,j);
           mat_t T=Givens(v,j,i+1);
           mat_mul(Ax,T,A);
          mat_mul(bt,T,b);
           mat_copy(b,bt);
           mat_copy(A,Ax);
         }
     }
    
    return mat_U_solve(A,bt);
}

//对矩阵A做Cholesky分解
mat_t Cholesky(mat_t A)
{
    int i,j,k,h;
    int m=A.m,n=A.n;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(A.mat[i][j]!=A.mat[j][i])
            {
                printf("矩阵非对称！\n");
                exit(-1);
            }                                  //疑问：正定与否如何判定？
        }
    }
    mat_t L=new_mat(m,n);
    L.mat[0][0]=sqrt(A.mat[0][0]);
        for(j=0;j<n;j++)
        {
            double sum1=0;
            for(k=0;k<=j-1;k++)
            {   
                sum1+=pow(L.mat[j][k],2);
            }
            L.mat[j][j]=sqrt(A.mat[j][j]-sum1);
            for(i=j+1;i<n;i++)
           {
            double sum2=0;
            for(h=0;h<=j-1;h++)
            {
                sum2+=(L.mat[i][h]*L.mat[j][h]);
            }
            L.mat[i][j]=(A.mat[i][j]-sum2)/L.mat[j][j]; 
        }
        }
        printf("对称正定矩阵的Cholesky分解L为:\n");
    mat_print(L);
    free_mat(L);
}



//Gasuss_Seidel迭代
//TODO 精度和迭代次数由用户指定
mat_t G_S(mat_t A,mat_t b)
{
    int m=A.m,n=A.n;
    int i,j,k;
    mat_t X=new_mat_vec(m);
    for(i=0;i<m;i++)
    {
        X.mat[i][0]=b.mat[i][0]/A.mat[i][i];
    }
    mat_t Xk=new_mat_clone(X);
    mat_t v=new_mat_row(m);
    for(k=0;k<=7;k++)
    {
        for(i=0;i<n;i++)
        {
           double sum1=0;
          for(j=0;j<=i-1;j++)
          {
              sum1+=A.mat[i][j]*Xk.mat[j][0];
          }

          double sum2=0;
          for(j=i+1;j<n;j++)
          {
              sum2+=A.mat[i][j]*X.mat[j][0];
          }

          Xk.mat[i][0]=-(1/A.mat[i][i])*(sum1+sum2-b.mat[i][0]);
         }
        mat_copy(X,Xk);
        mat_transpose(v,X);
        printf("第%d次运算结果为:\n",k+1);
        mat_print(v);
    }
}

double Li(mat_t v,int k,double c)//计算Larange插值公式里当x=c时的Li(x)值
{
    int i;
    int n=v.m;
    double a1=1,a2=1;
    for(i=0;i<n;i++)
    {
        if(k!=n-1)        //当k=n-1时，i=k+1会使数据访问越界，因此对k=n-1情形做判断
        {
         if(i==k)
         {  i=k+1; }
        }
        a1*=(c-v.mat[i][0]);
       if(i!=k)
       {
        a2*=(v.mat[k][0]-v.mat[i][0]);
       }
       else
           ;
    }
   return a1/a2;

    
}


//Lagrange插值公式计算当x=c时pn(x)的值
//mat_t x,mat_ty分别为n个点的x与y坐标按顺序生成的向量
double Lagrange(mat_t x,mat_t y,double c)
{
    int i;
    int n=x.m;
    double sum=0;
    for(i=0;i<n;i++)
    {
        sum+=( Li(x,i,c)*y.mat[i][0]);
    }
return sum;
}



//最速下降法求解方程Ax=b(要求A为正定矩阵！）;
//TODO 精度和迭代次数由用户指定
mat_t S_D(mat_t A,mat_t b)
{
   int m=A.m,n=A.n;
   int i,j,k;
   double a,a1,a2;
   mat_t R=new_mat_vec(m);
   mat_t R1=new_mat_vec(m);
   mat_t r=new_mat_vec(m);
   mat_t rk=new_mat_vec(m);
   mat_t x=new_mat_vec(m);
   mat_t xk=new_mat_clone(x);
for(k=0;k<=100;k++)           //进行100次迭代
{
   mat_mul(R,A,x);
   mat_sub(r,b,R);
   a1=dot_product(r,r);       //dot_product为计算两向量的内积函数
   mat_mul(R,A,r);
   a2=dot_product(R,r);
   a=a1/a2;
   mat_scaler(R,r,a);
   mat_add(xk,x,R);
   mat_copy(x,xk);
   mat_mul(R,A,r);
   mat_scaler(R1,R,a);
   mat_sub(rk,r,R1);
   mat_copy(r,rk);
}
  a=solution_judge(A,x,b,0.001);
   if(a<0.001)
  { 
      return x;
  }
 else
{
    printf("100次迭代后Ax-b二范数误差大于0.001");
}
}


//A为非亏损的n阶实矩阵,乘幂法求A得按模最大特征值,精确度要求为a
//TODO 精度和迭代次数由用户指定
double power(mat_t A,double a)
{
    int m=A.m;
    int k;
    double max,maxk;
    mat_t y=new_mat_vec(m);
    //取向量元素全为1的向量z；
    mat_t z=new_mat_vec(m);
    for(k=0;k<m;k++)
    {
        z.mat[k][0]=1;
    }
    mat_t zk=new_mat_vec(m);
    for(k=0;k<m;k++)
    {
        z.mat[k][0]=1;
    }
    //计算部分
    for(k=0;k<100;k++)  //设定最大计算次数为100次
    {
        mat_mul(y,A,z);
        mat_t Y=mat_absolute(y);
        max=mat_max(Y);
    //精度判断
        if(k>0)
        {
         if((max-maxk)<a&&(max-maxk)>(-a))
        {
            return max;
            break;
        }
        }

        mat_scaler(zk,y,1/max);
        mat_copy(z,zk);
        maxk=max;
    }
}


//Gauss消元法解发出那方程组Ax=b
mat_t Gauss(mat_t A, mat_t b) 
{
  if (A.m !=A.n) 
  {
    printf("输入矩阵非方阵！\n");
    exit(-1);
  } 
  else 
  {
    int m=A.m,n=A.n;
    int i, j;
    mat_t M=new_mat(m,n+1);
    mat_t G=new_mat(m,n);
    mat_t I=mat_scalar(n,1);
    mat_t Mx=new_mat(m,n+1);
    mat_t Gx=new_mat(m,n);
    mat_t L=new_mat(m,n);
    mat_t Lx=new_mat(m,n);
    mat_t Ly=new_mat(m,n);
    mat_t U=new_mat(m,n);
    mat_t R=new_mat(m,n);
    for (i = 0; i < m; i++) //计算矩阵U的增广矩阵M
    {
      for (j = 0; j < n; j++) 
      {M.mat[i][j] = A.mat[i][j];}
      M.mat[i][n] = b.mat[i][0];
    }
 
    int k;
    for(k=0;k<n-1;k++)
    {//取单位阵I的第K列
       mat_t r=new_mat_row(n);
       r.mat[0][k]=1;
      //取增广矩阵第k列的下半部分v
      mat_t v=new_vec_get(M,k);
      for (i = 0; i <n; i++) 
      {
        if (i <=k) 
        {v.mat[i][0] = 0;}
      }
      mat_mul(R,v,r);
      mat_scaler(G,R,1/M.mat[k][k]);
      mat_sub(Gx,I,G);
     mat_add(Lx,I,G);
      if (k > 0) 
      {mat_mul(L, Ly, Lx);}
       mat_mul(Mx, Gx, M);
      mat_copy(U,Mx);//计算LU分解的U
      
      mat_copy(M,Mx); //将计算中间过程的增广矩阵和算子迭代
      mat_copy(G,Gx);
      mat_copy(Ly,Lx);
    }
 /*   printf("所用LU分解为L=:\n");
    mat_print(L);
    printf("所用LU分解U=:\n");
    mat_print(U);
    printf("方程的解向量为:\n");*/
    b=new_vec_get(Mx,n);
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {A.mat[i][j]=Mx.mat[i][j];}
    }
    return mat_U_solve(A,b);
  }
}


//列主元消去法解方程Ax=b;
mat_t PLU(mat_t A,mat_t b)
{
   int m=A.m,n=A.n;
   int i,j,k,r;
   mat_t G=new_mat(m,n);
   mat_t Gk=new_mat(m,n);
   mat_t Ak=new_mat(m,n);
   mat_t R=new_mat(m,n);
   mat_t I=mat_scalar(m,1);
   mat_t bk=new_mat_vec(m);
       for(k=0;k<n-1;k++)
   {
      // 对第k列做主元选取
      mat_t Ac=new_mat_clone(A);
       for(i=k+1;i<m;i++)
       {
           if(Ac.mat[i][k]<0)              
           {Ac.mat[i][k]=0-Ac.mat[i][k];}
           if(Ac.mat[i][k]>Ac.mat[k][k])
           {
               Ac.mat[k][k]=Ac.mat[i][k];
               r=i;
           }
       }
      // 交换主元行
      mat_t A1=new_mat_clone(A);
      mat_t b1=new_mat_clone(b);
       for(j=0;j<n;j++)
       {
           A.mat[k][j]=A.mat[r][j];
           b.mat[k][0]=b.mat[r][0];
       }
       
       for(j=0;j<n;j++)
       {
           A.mat[r][j]=A1.mat[k][j];
           b.mat[r][0]=b1.mat[k][0];
       }
       //消元计算
     mat_t r=new_mat_row(n);
     r.mat[0][k]=1;
     mat_t v=new_vec_get(A,k);
     for(i=0;i<=k;i++)
     {
         v.mat[i][0]=0;
     }
     mat_mul(R,v,r);
     mat_scaler(G,R,1/A.mat[k][k]);
     mat_sub(Gk,I,G);
     mat_mul(Ak,Gk,A);mat_mul(bk,Gk,b);
     mat_copy(A,Ak);mat_copy(b,bk);
}
return mat_U_solve(A,b);
}


//追赶法解三对角方程组
mat_t Crout(mat_t A,mat_t b)
{
    int i,j,k;
    int n=A.m;
    double sum;
    mat_t L=new_mat(n,n);
    mat_t U=new_mat(n,n);
    for(k=0;k<n;k++)
    {
        L.mat[k][0]=A.mat[k][0];
        U.mat[0][k]=(A.mat[0][k])/(L.mat[0][0]);//先求的L的第一列和U的第一行
    }
    for(i=1;i<n;i++)
    {
        sum=0;
        for(j=0;j<i;j++)
        {sum+=(L.mat[i][j])*(U.mat[j][i]);}
         L.mat[i][i]=A.mat[i][i]-sum;     //计算L.mat[k][k]用于计算U时作为除数    
      for(k=i;k<n;k++)
      {
         U.mat[i][k]=(A.mat[i][k]-sum)/(L.mat[i][i]);
         //计算L的下一列   
         sum=0;
         for(j=0;j<i;j++)
         {sum+=(L.mat[k][j])*(U.mat[j][i]);}
          L.mat[k][i]=A.mat[k][i]-sum;
      }
    }
  mat_t y= mat_L_solve(L,b);
  mat_t x=mat_U_solve(U,y);
  return x;
}


//反幂法非奇异方阵A的按模最小特征值
//TODO 精度和迭代次数由用户指定
double in_power(mat_t A)
{
    int i,k;
    int n=A.m;
    double mk;
    mat_t y=new_mat_vec(n);
    mat_t yk=new_mat_vec(n);
    mat_t z=new_mat_vec(n);
    mat_t zk=new_mat_vec(n);
    for(i=0;i<n;i++)
    {
        z.mat[i][0]=1;
    }
    mat_t L=L_LU(A);
    mat_t U=U_LU(A);
    for(k=0;k<100;k++)//进行100次迭代计算
    {
        mat_t y=mat_L_solve(L,z);
        mat_t yk=mat_U_solve(U,y);
        mat_copy(y,mat_absolute(yk));
        mk=mat_max(y);
        mat_scaler(zk,yk,1/mk);
        mat_copy(z,zk);
    }
    return mk;
}


//TODO 加上区间转换以及对应的Gauss点和系数变换
//两点Gauss数值积分
double Gauss_2point(double f1,double f2)//f1和f2分别是被积函数在-sqrt(1/3)和sqrt(1/3)处的值
{
    double intgral;
    intgral=f1+f2;
    return intgral;
}


//1D有限元选取的基函数
double Psi(double x,int n,double h,double left,double right)
{
    double a;
    if(n==1)
    {
        if(left<=x&&x<=right)
        {
        a=(right-x)/h;
        return a;
        }
        else
        {
            return 0;
        }
    
    }
    if(n==2)
    {
        if(left<=x&&x<=right)
        {
        a=(x-left)/h;
        return a;
        }
        else
        {
            return 0;
        }

    }
}

double Phi_der(double x,int n,int j,double h,double left,double right)
{
    if(j==0)
    {
        if(left<=x&&x<=right)
        {
            return -1/h;
        }
        else
        {
            return 0;
        }
    }
    if(j==n)
    {
        if(right-h<=x&&x<=right)
        {
            return 1/h;
        }
        else
        {
            return 0;
        }
    }
    if(1<=j&&j<=n-1)
    {
        if(left<=x&&x<=right)
        {
            return 1/h;
        }
        if(right<=x&&x<=right+h)
        {
            return -1/h;
        }
        else
        {
            return 0;
        }
    }
}


double Phi(double x,int j,int n,double h,double left,double right)
{
 
    if(j==0)
    {
        if(left<=x&&x<=right)
        {
            return (right-x)/h;
        }
        else
        {
          return 0;
        }
    }
    if(j==n)
    {
        if(right-h<=x&&x<=right)
        {
           return 1;
        }
        else
        {
           return 0;
        }
    }
    if(0<=j&&j<n-1)
    {
        if(left<=x&&x<right)
        {
          return (x-left)/h;
          
        }
        if(right<=x&&x<=right+h)
        {
           return  (right+h-x)/h;
        }
        else
        {
         return 0;
        }
    }
    if(j==n-1)
    {
        if(left<=x&&x<right)
        {
            return (x-left)/h;
        }
        if(x==right)
        {
            return 0;
        }
    }
}


//1D有限元选取的局部线性基函数的一阶导数
double Psi_der(int n,double h)
{

    if(n==1)
    {
        return -1/h;
       
    }
    if(n==2)
    {
        return 1/h;
      
    }
}


//一维有限元算例1的函数f
double f(double x)
{
    double a;
    a=-exp(x)*(cos(x)-2*sin(x)-x*cos(x)-x*sin(x));
    return a;
 }


