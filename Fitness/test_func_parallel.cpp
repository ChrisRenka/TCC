/*
  CEC13 Test Function Suite 
  Jane Jing Liang (email: liangjing@zzu.edu.cn) 
  Last Modified on 14th Feb. 2013
*/


//#include <WINDOWS.H>      
#include"test_func_parallel.hpp"

extern double *OShift,*M,**y,**z,*x_bound;;
extern int ini_flag,n_flag,func_flag;

void setupCEC2013(int nx, int func_num, int nThreads){
    
    int cf_num=10, i;
    int scanRet;
    
    if (ini_flag==1)
	{
		if ((n_flag!=nx)||(func_flag!=func_num))
		{
			ini_flag=0;
		}
	}

	if (ini_flag==0)
	{
		FILE *fpt;
		char FileName[80];
		free(M);
		free(OShift);
		free(y);
		free(z);
		free(x_bound);
		//y=(double *)malloc(sizeof(double)  *  nx);
		//z=(double *)malloc(sizeof(double)  *  nx);
        y = (double**)malloc(nThreads*sizeof(double*));
        z = (double**)malloc(nThreads*sizeof(double*));
        for(i=0;i<nThreads;i++){
            y[i]=(double *)malloc(sizeof(double)  *  nx);
            z[i]=(double *)malloc(sizeof(double)  *  nx);
        }
		x_bound=(double *)malloc(sizeof(double)  *  nx);
		for (i=0; i<nx; i++)
			x_bound[i]=100.0;

		if (!(nx==2||nx==5||nx==10||nx==20||nx==30||nx==40||nx==50||nx==60||nx==70||nx==80||nx==90||nx==100))
		{
			printf("\nError: Test functions are only defined for D=2,5,10,20,30,40,50,60,70,80,90,100.\n");
		}
		
		sprintf(FileName, "../Fitness/input_data/M_D%d.txt", nx);
		fpt = fopen(FileName,"r");
		if (fpt==NULL)
		{
		    printf("\n Error: Cannot open input file for reading \n");
		}

		M=(double*)malloc(cf_num*nx*nx*sizeof(double));
		if (M==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for (i=0; i<cf_num*nx*nx; i++)
		{
				scanRet = fscanf(fpt,"%lf",&M[i]);
		}
		fclose(fpt);
		

		fpt=fopen("../Fitness/input_data/shift_data.txt","r");
		if (fpt==NULL)
		{
			printf("\n Error: Cannot open input file for reading \n");
		}
		OShift=(double *)malloc(nx*cf_num*sizeof(double));
		if (OShift==NULL)
			printf("\nError: there is insufficient memory available!\n");
		for(i=0;i<cf_num*nx;i++)
		{
				scanRet = fscanf(fpt,"%lf",&OShift[i]);
		}
		fclose(fpt);

		n_flag=nx;
		func_flag=func_num;
		ini_flag=1;
		//printf("Function has been initialized!\n");
	}
}

void test_func(double *x, double *f, int nx, int mx,int func_num, int id)
{
    
    int i;

	for (i = 0; i < mx; i++)
	{
		switch(func_num)
		{
		case 1:	
			sphere_func(&x[i*nx],&f[i],nx,OShift,M,0,id);
			f[i]+=-1400.0;
			break;
		case 2:	
			ellips_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-1300.0;
			break;
		case 3:	
			bent_cigar_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-1200.0;
			break;
		case 4:	
			discus_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-1100.0;
			break;
		case 5:
			dif_powers_func(&x[i*nx],&f[i],nx,OShift,M,0,id);
			f[i]+=-1000.0;
			break;
		case 6:
			rosenbrock_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-900.0;
			break;
		case 7:	
			schaffer_F7_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-800.0;
			break;
		case 8:	
			ackley_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-700.0;
			break;
		case 9:	
			weierstrass_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-600.0;
			break;
		case 10:	
			griewank_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-500.0;
			break;
		case 11:	
			rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,0,id);
			f[i]+=-400.0;
			break;
		case 12:	
			rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-300.0;
			break;
		case 13:	
			step_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=-200.0;
			break;
		case 14:	
			schwefel_func(&x[i*nx],&f[i],nx,OShift,M,0,id);
			f[i]+=-100.0;
			break;
		case 15:	
			schwefel_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=100.0;
			break;
		case 16:	
			katsuura_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=200.0;
			break;
		case 17:	
			bi_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,0,id);
			f[i]+=300.0;
			break;
		case 18:	
			bi_rastrigin_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=400.0;
			break;
		case 19:	
			grie_rosen_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=500.0;
			break;
		case 20:	
			escaffer6_func(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=600.0;
			break;
		case 21:	
			cf01(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=700.0;
			break;
		case 22:	
			cf02(&x[i*nx],&f[i],nx,OShift,M,0,id);
			f[i]+=800.0;
			break;
		case 23:	
			cf03(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=900.0;
			break;
		case 24:	
			cf04(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=1000.0;
			break;
		case 25:	
			cf05(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=1100.0;
			break;
		case 26:
			cf06(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=1200.0;
			break;
		case 27:
			cf07(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=1300.0;
			break;
		case 28:
			cf08(&x[i*nx],&f[i],nx,OShift,M,1,id);
			f[i]+=1400.0;
			break;
		default:
			printf("\nError: There are only 28 test functions in this test suite!\n");
			f[i] = 0.0;
			break;
		}
		
	}


}

void sphere_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Sphere */
{
	int i;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
	f[0] = 0.0;
    for (i=0; i<nx; i++)
    {
        f[0] += z[id][i]*z[id][i];
    }
}

void ellips_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Ellipsoidal */
{
    int i;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
    oszfunc (z[id], y[id], nx);
	f[0] = 0.0;
    for (i=0; i<nx; i++)
    {
        f[0] += pow(10.0,6.0*i/(nx-1))*y[id][i]*y[id][i];
    }
}

void bent_cigar_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Bent_Cigar */
{
    int i;
	double beta=0.5;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
    asyfunc (z[id], y[id], nx,beta);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	f[0] = z[id][0]*z[id][0];
    for (i=1; i<nx; i++)
    {
        f[0] += pow(10.0,6.0)*z[id][i]*z[id][i];
    }
}

void discus_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Discus */
{
    int i;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
    oszfunc (z[id], y[id], nx);

	f[0] = pow(10.0,6.0)*y[id][0]*y[id][0];
    for (i=1; i<nx; i++)
    {
        f[0] += y[id][i]*y[id][i];
    }
}

void dif_powers_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Different Powers */
{
	int i;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
	f[0] = 0.0;
    for (i=0; i<nx; i++)
    {
        f[0] += pow(fabs(z[id][i]),2+4*i/(nx-1));
    }
	f[0]=pow(f[0],0.5);
}


void rosenbrock_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Rosenbrock's */
{
    int i;
	double tmp1,tmp2;
	shiftfunc(x, y[id], nx, Os);//shift
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]=y[id][i]*2.048/100;
    }
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);//rotate
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
	for (i=0; i<nx; i++)//shift to orgin
    {
        z[id][i]=z[id][i]+1;
    }

    f[0] = 0.0;
    for (i=0; i<nx-1; i++)
    {
		tmp1=z[id][i]*z[id][i]-z[id][i+1];
		tmp2=z[id][i]-1.0;
        f[0] += 100.0*tmp1*tmp1 +tmp2*tmp2;
    }
}

void schaffer_F7_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Schwefel's 1.2  */
{
    int i;
	double tmp;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];
	asyfunc (z[id], y[id], nx, 0.5);
	for (i=0; i<nx; i++)
		z[id][i] = y[id][i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

	for (i=0; i<nx-1; i++)
		z[id][i]=pow(y[id][i]*y[id][i]+y[id][i+1]*y[id][i+1],0.5);
    f[0] = 0.0;
    for (i=0; i<nx-1; i++)
    {
	  tmp=sin(50.0*pow(z[id][i],0.2));
      f[0] += pow(z[id][i],0.5)+pow(z[id][i],0.5)*tmp*tmp ;
    }
	f[0] = f[0]*f[0]/(nx-1)/(nx-1);
}

void ackley_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Ackley's  */
{
    int i;
    double sum1, sum2;

	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	asyfunc (z[id], y[id], nx, 0.5);
	for (i=0; i<nx; i++)
		z[id][i] = y[id][i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

    sum1 = 0.0;
    sum2 = 0.0;
    for (i=0; i<nx; i++)
    {
        sum1 += y[id][i]*y[id][i];
        sum2 += cos(2.0*PI*y[id][i]);
    }
    sum1 = -0.2*sqrt(sum1/nx);
    sum2 /= nx;
    f[0] =  E - 20.0*exp(sum1) - exp(sum2) +20.0;
}


void weierstrass_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Weierstrass's  */
{
    int i,j,k_max;
    double sum,sum2, a, b;

	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]=y[id][i]*0.5/100;
    }
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	asyfunc (z[id], y[id], nx, 0.5);
	for (i=0; i<nx; i++)
		z[id][i] = y[id][i]*pow(10.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

    a = 0.5;
    b = 3.0;
    k_max = 20;
    f[0] = 0.0;
    for (i=0; i<nx; i++)
    {
        sum = 0.0;
		sum2 = 0.0;
        for (j=0; j<=k_max; j++)
        {
            sum += pow(a,j)*cos(2.0*PI*pow(b,j)*(y[id][i]+0.5));
			sum2 += pow(a,j)*cos(2.0*PI*pow(b,j)*0.5);
        }
        f[0] += sum;
    }
	f[0] -= nx*sum2;
}


void griewank_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Griewank's  */
{
    int i;
    double s, p;

	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]=y[id][i]*600.0/100.0;
    }
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	for (i=0; i<nx; i++)
		z[id][i] = z[id][i]*pow(100.0,1.0*i/(nx-1)/2.0);


    s = 0.0;
    p = 1.0;
    for (i=0; i<nx; i++)
    {
        s += z[id][i]*z[id][i];
        p *= cos(z[id][i]/sqrt(1.0+i));
    }
    f[0] = 1.0 + s/4000.0 - p;
}

void rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Rastrigin's  */
{
    int i;
	double alpha=10.0,beta=0.2;
	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]=y[id][i]*5.12/100;
    }

	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

    oszfunc (z[id], y[id], nx);
    asyfunc (y[id], z[id], nx, beta);

	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

	for (i=0; i<nx; i++)
	{
		y[id][i]*=pow(alpha,1.0*i/(nx-1)/2);
	}

	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

    f[0] = 0.0;
    for (i=0; i<nx; i++)
    {
        f[0] += (z[id][i]*z[id][i] - 10.0*cos(2.0*PI*z[id][i]) + 10.0);
    }
}

void step_rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Noncontinuous Rastrigin's  */
{
    int i;
	double alpha=10.0,beta=0.2;
	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]=y[id][i]*5.12/100;
    }

	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

    for (i=0; i<nx; i++)
	{
		if (fabs(z[id][i])>0.5)
		z[id][i]=floor(2*z[id][i]+0.5)/2;
	}

    oszfunc (z[id], y[id], nx);
    asyfunc (y[id], z[id], nx, beta);

	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

	for (i=0; i<nx; i++)
	{
		y[id][i]*=pow(alpha,1.0*i/(nx-1)/2);
	}

	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

    f[0] = 0.0;
    for (i=0; i<nx; i++)
    {
        f[0] += (z[id][i]*z[id][i] - 10.0*cos(2.0*PI*z[id][i]) + 10.0);
    }
}

void schwefel_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Schwefel's  */
{
    int i;
	double tmp;
	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]*=1000/100;
    }
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	for (i=0; i<nx; i++)
		y[id][i] = z[id][i]*pow(10.0,1.0*i/(nx-1)/2.0);

	for (i=0; i<nx; i++)
		z[id][i] = y[id][i]+4.209687462275036e+002;
	
    f[0]=0;
    for (i=0; i<nx; i++)
	{
		if (z[id][i]>500)
		{
			f[0]-=(500.0-fmod(z[id][i],500))*sin(pow(500.0-fmod(z[id][i],500),0.5));
			tmp=(z[id][i]-500.0)/100;
			f[0]+= tmp*tmp/nx;
		}
		else if (z[id][i]<-500)
		{
			f[0]-=(-500.0+fmod(fabs(z[id][i]),500))*sin(pow(500.0-fmod(fabs(z[id][i]),500),0.5));
			tmp=(z[id][i]+500.0)/100;
			f[0]+= tmp*tmp/nx;
		}
		else
			f[0]-=z[id][i]*sin(pow(fabs(z[id][i]),0.5));
    }
    f[0]=4.189828872724338e+002*nx+f[0];
}

void katsuura_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Katsuura  */
{
    int i,j;
	double temp,tmp1,tmp2,tmp3;
	tmp3=pow(1.0*nx,1.2);
	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]*=5.0/100.0;
    }
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	for (i=0; i<nx; i++)
		z[id][i] *=pow(100.0,1.0*i/(nx-1)/2.0);

	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

    f[0]=1.0;
    for (i=0; i<nx; i++)
	{
		temp=0.0;
		for (j=1; j<=32; j++)
		{
			tmp1=pow(2.0,j);
			tmp2=tmp1*y[id][i];
			temp += fabs(tmp2-floor(tmp2+0.5))/tmp1;
		}
		f[0] *= pow(1.0+(i+1)*temp,10.0/tmp3);
    }
	tmp1=10.0/nx/nx;
    f[0]=f[0]*tmp1-tmp1;

}

void bi_rastrigin_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Lunacek Bi_rastrigin Function */
{
    int i;
	double mu0=2.5,d=1.0,s,mu1,tmp,tmp1,tmp2;
	double *tmpx;
	tmpx=(double *)malloc(sizeof(double)  *  nx);
	s=1.0-1.0/(2.0*pow(nx+20.0,0.5)-8.2);
	mu1=-pow((mu0*mu0-d)/s,0.5);

	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]*=10.0/100.0;
    }

	for (i = 0; i < nx; i++)
    {
		tmpx[i]=2*y[id][i];
        if (Os[i] < 0.)
            tmpx[i] *= -1.;
    }

	for (i=0; i<nx; i++)
	{
		z[id][i]=tmpx[i];
		tmpx[i] += mu0;
	}
	if (r_flag==1)
	rotatefunc(z[id], y[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		y[id][i]=z[id][i];

	for (i=0; i<nx; i++)
		y[id][i] *=pow(100.0,1.0*i/(nx-1)/2.0);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

    tmp1=0.0;tmp2=0.0;
    for (i=0; i<nx; i++)
	{
		tmp = tmpx[i]-mu0;
		tmp1 += tmp*tmp;
		tmp = tmpx[i]-mu1;
		tmp2 += tmp*tmp;
    }
	tmp2 *= s;
	tmp2 += d*nx;
	tmp=0;
	for (i=0; i<nx; i++)
	{
		tmp+=cos(2.0*PI*z[id][i]);
    }
	
	if(tmp1<tmp2)
		f[0] = tmp1;
	else
		f[0] = tmp2;
	f[0] += 10.0*(nx-tmp);
	free(tmpx);
}

void grie_rosen_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Griewank-Rosenbrock  */
{
    int i;
    double temp,tmp1,tmp2;

	shiftfunc(x, y[id], nx, Os);
	for (i=0; i<nx; i++)//shrink to the orginal search range
    {
        y[id][i]=y[id][i]*5/100;
    }
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	for (i=0; i<nx; i++)//shift to orgin
    {
        z[id][i]=y[id][i]+1;
    }

    f[0]=0.0;
    for (i=0; i<nx-1; i++)
    {
		tmp1 = z[id][i]*z[id][i]-z[id][i+1];
		tmp2 = z[id][i]-1.0;
        temp = 100.0*tmp1*tmp1 + tmp2*tmp2;
         f[0] += (temp*temp)/4000.0 - cos(temp) + 1.0;
    }
	tmp1 = z[id][nx-1]*z[id][nx-1]-z[id][0];
	tmp2 = z[id][nx-1]-1.0;
    temp = 100.0*tmp1*tmp1 + tmp2*tmp2;;
     f[0] += (temp*temp)/4000.0 - cos(temp) + 1.0 ;
}


void escaffer6_func (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Expanded Scaffer��s F6  */
{
    int i;
    double temp1, temp2;
	shiftfunc(x, y[id], nx, Os);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, Mr);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

	asyfunc (z[id], y[id], nx, 0.5);
	if (r_flag==1)
	rotatefunc(y[id], z[id], nx, &Mr[nx*nx]);
	else
    for (i=0; i<nx; i++)
		z[id][i]=y[id][i];

    f[0] = 0.0;
    for (i=0; i<nx-1; i++)
    {
        temp1 = sin(sqrt(z[id][i]*z[id][i]+z[id][i+1]*z[id][i+1]));
		temp1 =temp1*temp1;
        temp2 = 1.0 + 0.001*(z[id][i]*z[id][i]+z[id][i+1]*z[id][i+1]);
        f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
    }
    temp1 = sin(sqrt(z[id][nx-1]*z[id][nx-1]+z[id][0]*z[id][0]));
	temp1 =temp1*temp1;
    temp2 = 1.0 + 0.001*(z[id][nx-1]*z[id][nx-1]+z[id][0]*z[id][0]);
    f[0] += 0.5 + (temp1-0.5)/(temp2*temp2);
}


void cf01 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 1 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10, 20, 30, 40, 50};
	double bias[5] = {0, 100, 200, 300, 400};
	
	i=0;
	rosenbrock_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/1e+4;
	i=1;
	dif_powers_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/1e+10;
	i=2;
	bent_cigar_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/1e+30;
	i=3;
	discus_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/1e+10;
	i=4;
	sphere_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],0, id);
	fit[i]=10000*fit[i]/1e+5;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf02 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 2 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	for(i=0;i<cf_num;i++)
	{
		schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	}
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf03 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 3 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	for(i=0;i<cf_num;i++)
	{
		schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	}
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf04 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 4 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {20,20,20};
	double bias[3] = {0, 100, 200};
	i=0;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/400;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf05 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 4 */
{
	int i,cf_num=3;
	double fit[3];
	double delta[3] = {10,30,50};
	double bias[3] = {0, 100, 200};
	i=0;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/400;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf06 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 6 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,10,10};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/4e+3;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/1e+3;
	i=2;
	ellips_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/1e+10;
	i=3;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/400;
	i=4;
	griewank_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=1000*fit[i]/100;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);

}

void cf07 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 7 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,10,10,20,20};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	griewank_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/100;
	i=1;
	rastrigin_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/1e+3;
	i=2;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	weierstrass_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/400;
	i=4;
	sphere_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],0, id);
	fit[i]=10000*fit[i]/1e+5;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void cf08 (double *x, double *f, int nx, double *Os,double *Mr,int r_flag,int id) /* Composition Function 8 */
{
	int i,cf_num=5;
	double fit[5];
	double delta[5] = {10,20,30,40,50};
	double bias[5] = {0, 100, 200, 300, 400};
	i=0;
	grie_rosen_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/4e+3;
	i=1;
	schaffer_F7_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/4e+6;
	i=2;
	schwefel_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/4e+3;
	i=3;
	escaffer6_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],r_flag, id);
	fit[i]=10000*fit[i]/2e+7;
	i=4;
	sphere_func(x,&fit[i],nx,&Os[i*nx],&Mr[i*nx*nx],0, id);
	fit[i]=10000*fit[i]/1e+5;
	cf_cal(x, f, nx, Os, delta,bias,fit,cf_num);
}

void shiftfunc (double *x, double *xshift, int nx,double *Os)
{
	int i;
    for (i=0; i<nx; i++)
    {
        xshift[i]=x[i]-Os[i];
    }
}

void rotatefunc (double *x, double *xrot, int nx,double *Mr)
{
	int i,j;
    for (i=0; i<nx; i++)
    {
        xrot[i]=0;
			for (j=0; j<nx; j++)
			{
				xrot[i]=xrot[i]+x[j]*Mr[i*nx+j];
			}
    }
}

void asyfunc (double *x, double *xasy, int nx, double beta)
{
	int i;
    for (i=0; i<nx; i++)
    {
		if (x[i]>0)
        xasy[i]=pow(x[i],1.0+beta*i/(nx-1)*pow(x[i],0.5));
    }
}

void oszfunc (double *x, double *xosz, int nx)
{
	int i,sx;
	double c1,c2,xx;
    for (i=0; i<nx; i++)
    {
		if (i==0||i==nx-1)
        {
			if (x[i]!=0)
				xx=log(fabs(x[i]));
			if (x[i]>0)
			{	
				c1=10;
				c2=7.9;
			}
			else
			{
				c1=5.5;
				c2=3.1;
			}	
			if (x[i]>0)
				sx=1;
			else if (x[i]==0)
				sx=0;
			else
				sx=-1;
			xosz[i]=sx*exp(xx+0.049*(sin(c1*xx)+sin(c2*xx)));
		}
		else
			xosz[i]=x[i];
    }
}


void cf_cal(double *x, double *f, int nx, double *Os,double * delta,double * bias,double * fit, int cf_num)
{
	int i,j;
	double *w;
	double w_max=0,w_sum=0;
	w=(double *)malloc(cf_num * sizeof(double));
	for (i=0; i<cf_num; i++)
	{
		fit[i]+=bias[i];
		w[i]=0;
		for (j=0; j<nx; j++)
		{
			w[i]+=pow(x[j]-Os[i*nx+j],2.0);
		}
		if (w[i]!=0)
			w[i]=pow(1.0/w[i],0.5)*exp(-w[i]/2.0/nx/pow(delta[i],2.0));
		else
			w[i]=INF;
		if (w[i]>w_max)
			w_max=w[i];
	}

	for (i=0; i<cf_num; i++)
	{
		w_sum=w_sum+w[i];
	}
	if(w_max==0)
	{
		for (i=0; i<cf_num; i++)
			w[i]=1;
		w_sum=cf_num;
	}
	f[0] = 0.0;
    for (i=0; i<cf_num; i++)
    {
		f[0]=f[0]+w[i]/w_sum*fit[i];
    }
	free(w);
}

