#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
double d_quotient(mat_t v1,mat_t v2,double k)//向量v1和v2分别是点组的x和y构成的向量，求其前k(从0计数,故对于n个元素的向量，k最大为n-1)个元素的差商
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
 

double Newton(mat_t x, mat_t y,double c)//Newton法求在新坐标x=c时的函数值
{
    int i,j;
    int n=x.m;
    double c1=1,c2=y.mat[0][0];
    for(i=0;i<n-1;i++)//由于在42行使用d_quotient()函数求差商时使用了i+1，这里规定i<n-1，防止数据访问越界
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
int main()
{
    int n;
    printf("点组个数n；\n");
   scanf("%d",&n);
   mat_t x=new_mat_vec(n);
   mat_t y=new_mat_vec(n);
   mat_set_all(x);
   mat_set_all(y);
   double c;
   printf("新坐标x=\n");
   scanf("%lf",&c);
   printf("新坐标对应函数值为%lf\n",Newton(x,y,c));
  free_mat(x);
  free_mat(y);
  return 0;
}
