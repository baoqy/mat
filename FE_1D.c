#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"mat.c"
void FE_1D(double left,double right,int n,double x)
{
    int i,j,k;
    double h=(right-left)/n;
       //A的封装
       int Nb=2;//在一个单元里选取的基函数有两段；
       double x1,x2;
       mat_t A=new_mat(n+1,n+1);
       mat_t S=new_mat(Nb,Nb);
       mat_t Ac=new_mat(n+1,n+1);
       for(k=0;k<n+1;k++)
       {
           for(i=0;i<Nb;i++)
           {
               for(j=0;j<Nb;j++)
               {
                   x1=(left+k*h+left+(k+1)*h)/2+h*(-1/sqrt(3))/2;//Gauss点转换,采取两点
                   x2=(left+k*h+left+(k+1)*h)/2+h*(1/sqrt(3))/2;
                   S.mat[j][i]=(h/2)*Gauss_2point(exp(x1)*Psi_der(i,h)*Psi_der(j,h),exp(x2)*Psi_der(i,h)*Psi_der(j,h));
         //     printf("x1=%lf,  x2=%lf,  Psi_der(%d,h)=%lf,  Psi_der(%d,h)=%lf,  S.mat[%d][%d]=%lf\n",x1,x2,i,Psi_der(i,h),j,Psi_der(j,h),j,i,S.mat[j][i]);
               }
           }
           
      
for(i=k;i<=n-1;i++)
{
    A.mat[i][i]+=S.mat[0][0];
    A.mat[i][i+1]+=S.mat[0][1];
    A.mat[i+1][i]+=S.mat[1][0];
    A.mat[i+1][i+1]+=S.mat[1][1];
       }
}
for(i=0;i<=n;i++)
{
    A.mat[0][i]=0;
    A.mat[n][i]=0;
}
mat_set_one(A,0,0,1);
mat_set_one(A,n,n,1);
printf("A=\n");
       mat_print(A);
//b矩阵的封装
mat_t b=new_mat_vec(n+1);
mat_t R=new_mat_vec(n+1);
for(i=0;i<=n;i++)
{
  /* for(j=0;j<Nb;j++)
   {
        x1=(left+i*h+left+(i+1)*h)/2+h*(-1/sqrt(3))/2;
        x2=(left+i*h+left+(i+1)*h)/2+h*(1/sqrt(3))/2;
        R.mat[j][0]=(h/2)*Gauss_2point(f(x1)*Psi(x1,j,h,left+i*h,left+(i+1)*h),f(x2)*Psi(x2,j,h,left+i*h,left+(i+1)*h));
      //  printf("R.mat[%d][0]=%lf\n",j,R.mat[j][0]);
      }
       b.mat[i][0]+=R.mat[0][0];
       b.mat[i+1][0]+=R.mat[1][0];
 
}*/
x1=(left+i*h+left+(i+1)*h)/2+h*(-1/sqrt(3))/2;
x2=(left+i*h+left+(i+1)*h)/2+h*(1/sqrt(3))/2;
b.mat[i][0]=(h/2)*Gauss_2point(f(x1)*Phi(x1,i,n,h,left+i*h,left+(i+1)*h),f(x2)*Phi(x2,i,n,h,left+i*h,left+(i+1)*h));
}
printf("b=\n");
mat_set_one(b,0,0,0);
mat_set_one(b,4,0,cos(1));
mat_print(b);
//得出系数矩阵
mat_t v=new_mat_vec(n+1);
mat_copy(v,Crout(A,b));
printf("解为\n");
mat_print(v);
//误差估计
//构造插值
double a=0;
for(j=0;j<=n;j++)
{   
           a+=v.mat[j][0]*Phi(x,j,n,h,left+j*h,left+(j+1)*h);
            printf("%lf=%lf*%lf\n",a,v.mat[j][0],Phi(x,j,n,h,left+j*h,left+(j+1)*h));
        
/*    if(j==n)
    {
        a+=v.mat[j][0]*Phi(x,j,n,h,right-h,right);
        printf("%lf=%lf*%lf\n",a,v.mat[j][0],Phi(x,j,n,h,right-h,right));
    }*/


}}
int main()
{
    double left=0,right=1;
    int n;
    double x;
    printf("单元数n=\n");
    scanf("%d",&n);
    printf("测试点x=:\n");
    scanf("%lf",&x);
    FE_1D(left,right,n,x);
    return 0;
}
