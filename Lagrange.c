#include<stdio.h>
#include<stdlib.h>
#include"mat.c"
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


double Lagrange(mat_t x,mat_t y,double c)//Lagrange插值公式计算当x=c时pn(x)的值
                                          //mat_t x,mat_ty分别为n个点的x与y坐标按顺序生成的向量
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

int main()
{
    int n;
    double c;
    printf("点组的个数n:\n");
    scanf("%d",&n);
    mat_t x=new_mat_vec(n);
    mat_t y=new_mat_vec(n);
    mat_set_all(x);
    mat_set_all(y);
    printf("计算点的横坐标为:\n");
    scanf("%lf",&c);
    printf("计算结果为%lf\n",Lagrange(x,y,c));
    free_mat(x);
    free_mat(y);
    return 0;
}
