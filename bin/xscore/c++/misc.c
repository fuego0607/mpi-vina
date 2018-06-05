# include "xtool.h"

float Distance(const float a[3], const float b[3])
{
	double d,tmpx,tmpy,tmpz;

	tmpx=(a[0]-b[0])*(a[0]-b[0]);
	tmpy=(a[1]-b[1])*(a[1]-b[1]);
	tmpz=(a[2]-b[2])*(a[2]-b[2]);

	d=sqrt(tmpx+tmpy+tmpz);

	return (float)d;
}

float Distance2(const float a[3], const float b[3])
{
	double d,tmpx,tmpy,tmpz;

	tmpx=(a[0]-b[0])*(a[0]-b[0]);
	tmpy=(a[1]-b[1])*(a[1]-b[1]);
	tmpz=(a[2]-b[2])*(a[2]-b[2]);

	d=tmpx+tmpy+tmpz;

	return (float)d;
}

void Translate_Point(const float start[3], const float move[3], float end[3])
{
	end[0]=start[0]+move[0];
	end[1]=start[1]+move[1];
	end[2]=start[2]+move[2];

	return;
}

void Rotate_Point(float theta, float axis[3], const float origin[3], 
		  const float start[3], float end[3])
{
        double mat[3][3],p[3];
        double a,b;
        int i,j;

	Normal_Vector(axis);

        a=sin((double)(theta*(PI)/180.0));
        b=cos((double)(theta*(PI)/180.0));

        for(i=0;i<=2;i++) end[i]=start[i]-origin[i];

        mat[0][0]=axis[0]*axis[0]+(1-axis[0]*axis[0])*b;
        mat[1][0]=axis[0]*axis[1]*(1-b)-axis[2]*a;
        mat[2][0]=axis[0]*axis[2]*(1-b)+axis[1]*a;
        mat[0][1]=axis[0]*axis[1]*(1-b)+axis[2]*a;
        mat[1][1]=axis[1]*axis[1]+(1-axis[1]*axis[1])*b;
        mat[2][1]=axis[1]*axis[2]*(1-b)-axis[0]*a;
        mat[0][2]=axis[0]*axis[2]*(1-b)-axis[1]*a;
        mat[1][2]=axis[1]*axis[2]*(1-b)+axis[0]*a;
        mat[2][2]=axis[2]*axis[2]+(1-axis[2]*axis[2])*b;

        for(i=0;i<=2;i++)
                {
                 p[i]=0.0;
                 for(j=0;j<=2;j++) p[i]+=(end[j]*mat[j][i]);
                }

        for(i=0;i<=2;i++) end[i]=p[i]+origin[i];

	return;
}

void Cross_Multiply(const float v1[3], const float v2[3], float result[3])
{
        result[0]=v1[1]*v2[2]-v2[1]*v1[2];
        result[1]=v1[2]*v2[0]-v2[2]*v1[0];
        result[2]=v1[0]*v2[1]-v2[0]*v1[1];

        return;
}

float Dot_Multiply(const float v1[3], const float v2[3])
{
	float result;

	result=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];

	return result;
}

float Triple_Multiply(const float v1[3], const float v2[3], const float v3[3])
{
	float result,tmp[3];

	// v1[]*v2[]

	tmp[0]=v1[1]*v2[2]-v2[1]*v1[2];
        tmp[1]=v1[2]*v2[0]-v2[2]*v1[0];
        tmp[2]=v1[0]*v2[1]-v2[0]*v1[1];

	// tmp[].v3[]

	result=tmp[0]*v3[0]+tmp[1]*v3[1]+tmp[2]*v3[2];

	return result;
}

float Angle_of_Two_Vectors(const float v1[3], const float v2[3])
{
        double angle;
        double l1,l2,tmp1,tmp2;

        l1=sqrt(v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]);
        l2=sqrt(v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);

        tmp1=v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2];
	tmp2=l1*l2;

        angle=acos(tmp1/tmp2); angle=angle/(PI)*180.0;

        return (float)angle;  // return angle in degree, 0-180
}

float Angle(const float a[3], const float b[3], const float c[3])
{
	int i;
	float angle,v1[3],v2[3];

	for(i=0;i<3;i++) {v1[i]=b[i]-a[i]; v2[i]=b[i]-c[i];}

	angle=Angle_of_Two_Vectors(v1,v2);

	return angle;  // the a-b-c angle in degree 0-180
}

float Torsion_Angle(const float p1[3], const float p2[3], 
		    const float p3[3], const float p4[3])
{
        int i;
        float angle;
        float v12[3],v23[3],v34[3],v1[3],v2[3],v3[3];

        for(i=0;i<=2;i++)
                {
                 v12[i]=p2[i]-p1[i];
                 v23[i]=p3[i]-p2[i];
                 v34[i]=p4[i]-p3[i];
                }

        Cross_Multiply(v12,v23,v1);
        Cross_Multiply(v23,v34,v2);

        angle=Angle_of_Two_Vectors(v1,v2);

        Cross_Multiply(v1,v2,v3);

	if(Angle_of_Two_Vectors(v23,v3)<=90.0) return angle;  
	else return -angle;

	// returned angle is in degree, -180~180
}

float Angle_of_Two_Planes(const float p1[3], const float p2[3],
                          const float p3[3], const float q1[3],
                          const float q2[3], const float q3[3])
{
	double angle,tmp1,tmp2;
	double A1,B1,C1,A2,B2,C2;

	A1=(p2[1]-p1[1])*(p3[2]-p1[2])-(p3[1]-p1[1])*(p2[2]-p1[2]);
	B1=(p3[0]-p1[0])*(p2[2]-p1[2])-(p2[0]-p1[0])*(p3[2]-p1[2]);
	C1=(p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]);

	A2=(q2[1]-q1[1])*(q3[2]-q1[2])-(q3[1]-q1[1])*(q2[2]-q1[2]);
        B2=(q3[0]-q1[0])*(q2[2]-q1[2])-(q2[0]-q1[0])*(q3[2]-q1[2]);
        C2=(q2[0]-q1[0])*(q3[1]-q1[1])-(q3[0]-q1[0])*(q2[1]-q1[1]);

	tmp1=A1*A2+B1*B2+C1*C2;
	tmp2=sqrt(A1*A1+B1*B1+C1*C1)*sqrt(A2*A2+B2*B2+C2*C2);

	angle=acos(tmp1/tmp2); angle=angle/(PI)*180.0;

	if(angle>90.0) angle=fabs(180.0-angle);

        return (float)angle;  // returned angle is in degree
}

// point[] is a point out of the plane(p1,p2,p3) 
// while centroid[] is a point on the plane 
float Angle_of_Point_and_Plane(const float point[3],
                               const float centroid[3],
                               const float p1[3],
                               const float p2[3],
                               const float p3[3])
{
	float angle;
	float v1[3],v2[3];
	float A,B,C;

	// determine the direction of the plane first

	A=(p2[1]-p1[1])*(p3[2]-p1[2])-(p3[1]-p1[1])*(p2[2]-p1[2]);
	B=(p3[0]-p1[0])*(p2[2]-p1[2])-(p2[0]-p1[0])*(p3[2]-p1[2]);
	C=(p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]);

	// use the centroid as the origin

	v1[0]=A; v1[1]=B; v1[2]=C;

	v2[0]=point[0]-centroid[0];
	v2[1]=point[1]-centroid[1];
	v2[2]=point[2]-centroid[2];

	angle=Angle_of_Two_Vectors(v1,v2);

	if(fabs(angle)>90.000) angle=180.0-angle;

	return (float)angle;  // return angle in degree, 0-90
}

void Unit_Vector(const float start[3], const float end[3], float result[3])
{
	float d,tmpx,tmpy,tmpz;

	tmpx=(end[0]-start[0])*(end[0]-start[0]);
	tmpy=(end[1]-start[1])*(end[1]-start[1]);
	tmpz=(end[2]-start[2])*(end[2]-start[2]);

	d=sqrt((double)(tmpx+tmpy+tmpz));

	result[0]=(end[0]-start[0])/d;
	result[1]=(end[1]-start[1])/d;
	result[2]=(end[2]-start[2])/d;

	return;
}

void Normal_Vector(float v[3])
{
	float d;

	d=sqrt((double)(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));

	v[0]/=d; v[1]/=d; v[2]/=d;

	return;
}

float Distance_From_Point_To_Line(const float origin[3], 
				  const float axis[3],
				  const float point[3])
{
	int i;
	float vec[3];
	float dt,d2,dist;

	for(i=0;i<3;i++) vec[i]=point[i]-origin[i];

	dt=vec[0]*axis[0]+vec[1]*axis[1]+vec[2]*axis[2];

	d2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]-dt*dt;

	if(d2<=0.000) d2=0.000;
	else dist=sqrt(d2);

	return dist;
}

void Set_Random_Number()
{
	long seed;
	struct tm *ptr;
	time_t lt;

	lt=time(NULL); ptr=localtime(&lt);

	seed=10000*ptr->tm_hour+100*ptr->tm_min+ptr->tm_sec;

	srand48(seed);	// Initialize drand48()

	return;
}

char *Get_Time()  // get local time in ascii format
{
	extern char *timestring;    
        struct tm *ptr;
        time_t lt;

        lt=time(NULL); ptr=localtime(&lt);
	timestring=asctime(ptr);

	Chomp(timestring); return timestring;
}

char *Chomp(char *string)  // remove the new-line at the end of a string
{
	int i,len;

	len=strlen(string); if(len<1) return NULL;

	for(i=len-1;i>=0;i--)
                {
                 if(string[i]!='\n') continue;
                 else {string[i]='\0'; break;}
                }

	return string;
}

void Get_Sp3_Coordinates(const float v1[3], const float v2[3], float v3[3], float v4[3])
// v1 and v2 are known vectors while v3 and v4 are desired vectors.
// v1, v2, v3, and v4 are all unified vectors.
{
	int i;
	float anchor[3], axis[3];

	for(i=0;i<3;i++) 
		{
		 axis[i]=-(v1[i]+v2[i]);
		 anchor[i]=0.000;
		 v3[i]=v4[i]=0.000;
		}

	Rotate_Point(90.0,axis,anchor,v1,v3);
	Rotate_Point(90.0,axis,anchor,v2,v4);

	for(i=0;i<3;i++)
		{
		 v3[i]=-v3[i];
		 v4[i]=-v4[i];
		}

	return;
}

void Get_Sp2_Coordinates(const float v1[3], const float v2[3], float v3[3])
// v1 and v2 are known vectors and v3 is the desired vector.
// v1, v2, and v3 are all unified vectors.
{
	int i;

	for(i=0;i<3;i++)
		{
		 v3[i]=-(v1[i]+v2[i]);
		}

	Normal_Vector(v3);

	return;
}

int Pow(int x, int y)
{
	int i,result;

	if(y<=0) return 1;
	
	result=1;

	for(i=1;i<=y;i++) result*=x;

	return result;
}

int Blank_Line_Check(const char *line)
{
	int i,mark;
	int len=0;

	len=strlen(line);
	if(len<1) return TRUE;	// blank line!

	mark=0;

	for(i=0;i<len;i++)
		{
		 if(line[i]==' ') continue;
		 else if(line[i]=='\t') continue;
		 else if(line[i]=='\n') break;
		 else mark++;
		}

	if(mark==0) return TRUE;	// blank line
	else return FALSE;		// not a blank line
}	

void Memory_Allocation_Error(char *position)
{
	printf("\n");
	printf("Memory allocation error!\n");
	printf("This is usually caused by memory overflow.\n");
	printf("Error happened at %s\n", position);
	exit(1);
}

void Open_File_Error(const char *filename)
{
	printf("\n");
	printf("Error: cannot open the file %s\n", filename);
	printf("Please make sure it exists.\n");
	exit(1);
}

void Read_File_Error(const char *filename)
{
	printf("\n");
	printf("Error: something wrong with %s\n", filename);
	printf("It may not have the correct format.\n");
	printf("Please check this file and try again.\n");
	exit(1);
}

void PDB_Format_Error(const char *filename)
{
        printf("\n");
        printf("Error: %s lacks necessary information.\n", filename);
        printf("This file may not be in PDB format.\n");
        printf("Please check it and try again.\n");
        exit(1);
}

void Mol2_Format_Error(const char *filename)
{
	printf("\n");
	printf("Error: %s lacks necessary information.\n", filename);
	printf("This file may not be in Mol2 format.\n");
	printf("Please check it and try again.\n");
	exit(1);
}

void Lig_Format_Error(const char *filename)
{
	printf("\n");
	printf("Error: %s lacks necessary information.\n", filename);
	printf("This file may not be in Lig format.\n");
	printf("Please check it and try again.\n");
	exit(1);
}

int Check_Mol2_File(const char *filename)
// return the number of molecules in the given file.
// if it is not a valid Mol2 file, return 0.
{
	FILE *fp;
	int count;
	char buf[256],head[256];

        if((fp=fopen(filename,"r"))==NULL) Open_File_Error(filename);

        count=0;

        for(;;)
                {
                 if(fgets(buf,256,fp)==NULL) break; 
                 else {strcpy(head,""); sscanf(buf,"%s",head);}

                 if(strcmp(head,"@<TRIPOS>MOLECULE")) continue;
                 else count++; 
                }

	fclose(fp);

	return count;
}

void Int_To_Char(int number, char string[], int digit)
// the difference between this function and 'sprintf' is:
// this function fills '0' on the vacant digits
{
	char tmp_string[80];
	int i,j,tmp,len;

	if(number<0) number=-number; 	// 'number' should not be negative 
	if(digit<=0) digit=3;		// default value

	strcpy(tmp_string,"");
	sprintf(tmp_string,"%d", number);
	len=strlen(tmp_string);

	if(len==digit) 
		{
		 strcpy(string,tmp_string); 
		}
	else
		{
		 tmp=digit-len;

		 for(i=0;i<tmp;i++) string[i]='0';	// fill in the zeros	

		 for(j=0,i=tmp;i<=digit;i++,j++) 
			{
			 string[i]=tmp_string[j];
			}
		}

	return;
}

float Calculate_RMSD(const Molecule &mol1, const Molecule &mol2)
{
	int i,num,count;
	float rmsd;

	if(mol1.num_atom!=mol2.num_atom) return -1.000;	// sth wrong 
	else num=mol1.num_atom;

	count=0; rmsd=0.000;

	for(i=0;i<num;i++)
		{
		 if(mol1.atom[i].type[0]=='H') continue;
		 else if(mol2.atom[i].type[0]=='H') continue;

		 rmsd+=Distance2(mol1.atom[i].coor,mol2.atom[i].coor);
		 count++;
		}

	if(count==0) return -1.000;	// sth wrong

	rmsd/=count; rmsd=sqrt(rmsd);

	return rmsd;
}

void Sort_Array(int n, double array[], char *order)
{
	int i,j;
	double temp;

	if(n<=0) 
		{
		 puts("Wrong usage in Sort_Array()");
		 exit(1);
		}

	if(n<=1) return; // unnecessary to rank

	if(!strcasecmp(order,"decreasing"))
		{
		 for(i=0;i<n-1;i++)
		 for(j=i+1;j<n;j++)
			{
			 if(array[i]>=array[j]) continue;
			 else {SWAP(array[i],array[j]);}
			}
		}
	else 	// increasing
		{
		 for(i=0;i<n-1;i++)
		 for(j=i+1;j<n;j++)
			{
			 if(array[i]<=array[j]) continue;
			 else {SWAP(array[i],array[j]);}
			}
		}

	return;
}

void Simplex_Minimization(int ndim, int mpts, float **p, float y[], float ftol,
		          float (*funk)(int, float []), int &nfunk)
/* use the down-hill simplex method to minimize a multi-dimensional function 
   inputs are: ndim, number of the independent variables;
               mpts, number of the initial points;
               p[mpts][ndim], values of the variables for the initial points;
               y[mpts], values of the function for the initial points;
               ftol, the convergence tolerance to be met by the minimization;
	       *funk, the energy function;
   outputs are: new points all within ftol of a minima, as in p[][] and y[];
               best point in slot 0;
               nfunk: number of function evaluations taken */ 
{
	static int NMAX=1000;
	int i,j,ihi,ilo,inhi;
	float rtol,sum,ysave,ytry,temp;
	float *psum=NULL;

	psum=new float[ndim];
	if(psum==NULL) Memory_Allocation_Error(); 

	nfunk=0;

	// get psum, an assembly of all the points

	for(j=0;j<ndim;j++)
		{
		 for(sum=0.0,i=0;i<mpts;i++) sum+=p[i][j];
		 psum[j]=sum;
		}

	for(;;)
       {
	ilo=0;

	// first, determine the highest(worst), next highest, and lowest(best)
	// ihi, the highest; inhi, the next highest; ilo, the lowest 

	ihi=y[0]>y[1]?(inhi=1,0):(inhi=0,1);

	for(i=0;i<mpts;i++)
		{
		 if(y[i]<=y[ilo]) ilo=i;
		 if(y[i]>y[ihi]) {inhi=ihi; ihi=i;}
		 else if((y[i]>y[inhi])&&(i!=ihi)) inhi=i;
		 else continue;
		}

	rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

	// compute the fractional range from highest to lowest and 
	// return if satisfactory

/*
	printf("Energy evaluation %d: ", nfunk);
	printf("best energy = %f ", y[ilo]);
	printf("second highest energy = %f ", y[inhi]);
	printf("highest energy = %f\n", y[ihi]);
*/

	if(rtol<ftol)	// if returning, put best point in slot 0 
		{
		 temp=y[0]; y[0]=y[ilo]; y[ilo]=temp;
		 for(i=0;i<ndim;i++) 
			{temp=p[0][i]; p[0][i]=p[ilo][i]; p[ilo][i]=temp;}
		 break;
		}

	if(nfunk>=NMAX) break;	// maximal iteration exceeded

	// begin a new iteration. First extrapolate by a factor of -1
	// through the face of the simplex across from the high point
	// i.e. reflect the simplex from the high point

	ytry=Simplex_Generation(ndim,mpts,p,y,psum,funk,ihi,-1.0);
	nfunk++;

	if(ytry<=y[ilo])
		{
		 // gives a result better than the best point, so try an
		 // additional extrapolation by a factor of 2

		 ytry=Simplex_Generation(ndim,mpts,p,y,psum,funk,ihi,2.0);
		 nfunk++;
		}
	else if(ytry>=y[inhi])
		{
		 // the reflected point is worse than the second-highest,
		 // so look for an intermediate lower point, i.e. do a one-
		 // dimensional contraction

		 ysave=y[ihi];

		 ytry=Simplex_Generation(ndim,mpts,p,y,psum,funk,ihi,0.5);
		 nfunk++;

		 if(ytry>=ysave)
			{
			 // cannot seem to get rid of that high point. Better
			 // contract around the lowest point.

			 for(i=0;i<mpts;i++)
				{
				 if(i!=ilo)
				{
				 for(j=0;j<ndim;j++)
					{
					 p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
					}
				 y[i]=(*funk)(ndim,psum);
				 nfunk++;	// keep track of function evaluation
				}
				}

			 // get new psum

			 for(j=0;j<ndim;j++)
                		{
                 		 for(sum=0.0,i=0;i<mpts;i++) sum+=p[i][j];
                 		 psum[j]=sum;
                		}
			}
		}
	else continue; 
       }

	if(psum) delete [] psum;

	return;
}

float Simplex_Generation(int ndim, int mpts, float **p, float y[], float psum[],
                         float (*funk)(int, float []), int ihi, float factor)
/* Extrapolate by a factor through the face of the simplex across from the high
   point, tries it, and replaces the high point if the new point is better. */
{
	int i;
	float ytry;
	float *ptry=NULL;

	ptry=new float[ndim];
	if(ptry==NULL) Memory_Allocation_Error();

	for(i=0;i<ndim;i++)
		{
		 ptry[i]=(1.0-factor)*(psum[i]-p[ihi][i])/(mpts-1)+factor*p[ihi][i];
		}
	
	ytry=(*funk)(ndim,ptry);   // evaluate the function at the trial point

	if(ytry<y[ihi])	// if it is better than the highest, then replace it 
		{
		 y[ihi]=ytry;
		 for(i=0;i<ndim;i++)
			{
			 psum[i]+=(ptry[i]-p[ihi][i]);
			 p[ihi][i]=ptry[i];
			}
		}

	if(ptry) delete [] ptry;

	return ytry;
}

void Powell_Minimization(int n, float p[], float ftol,
			 int &iter, float &fret, float (*func)(int, float[]))
/* Minimization of a function func() of n variables.
   On input:
   p[1,..,n] is an inital starting point;
   ftol is the fractional tolerance in the function value such that failure
   to decrease by more than this amount on one iteration signals doneness;
   On output:
   p[1..n] is set to the best point found; 
   fret is the returned function value at p; 
   iter is the number of iterations taken. */ 
{
	int i,j,ibig;
	float del,fp,fptt,t,*pt,*ptt,*xit;
	int ITMAX=100;

/*
   	xi[1..n][1..n] is the initial matrix, whose columns contain the initial
   	set of directions (usually the n unit vectors);
   	On output, xi[1..n][1..n] is the then-current direction set;
*/
	float **xi; xi=new float*[n]; for(i=0;i<n;i++) xi[i]=new float[n];	

	for(i=0;i<n;i++) 
	for(j=0;j<n;j++) 
		{
		 if(i==j) xi[i][j]=1.000;
		 else xi[i][j]=0.000;
		}

	pt=new float[n]; ptt=new float[n]; xit=new float[n];
	fret=(*func)(n,p);

	for(j=0;j<n;j++) pt[j]=p[j];	// save the initial point

	for(iter=1;iter<=ITMAX;iter++)
	{
	 fp=fret; ibig=0; 

	 del=0.0;  // will be the biggest function decrease 

	 for(i=0;i<n;i++)	// loop over all directions in the set
		{
		 for(j=0;j<n;j++) xit[j]=xi[j][i];  // copy the direction
		 fptt=fret;
		 linmin(p,xit,n,fret,func);	// minimize along it
		 if((fptt-fret)>del)
			{
			 del=fptt-fret; ibig=i;
			}
		}

	 if(2.0*(fp-fret)<=ftol*(fabs(fp)+fabs(fret))+TINY)
		{
		 break;	 // termination criterion fulfilled 
		}

	 for(j=0;j<n;j++)
		{
		 ptt[j]=2.0*p[j]-pt[j];
		 xit[j]=p[j]-pt[j];
		 pt[j]=p[j];
		}

	 fptt=(*func)(n,ptt);  // function value at extrapolated point

	 if(fptt<fp)
		{
		 t=2.0*(fp-2.0*fret+fptt)*SQR(fp-fret-del)-del*SQR(fp-fptt);
		 if(t<0.0)
			{
			 linmin(p,xit,n,fret,func);
			 for(j=0;j<n;j++)
				{
				 xi[j][ibig]=xi[j][n];
				 xi[j][n]=xit[j];
				}
			}
		}

	 // printf("Iteration %d: e = %6.1f\n", iter, fret);
	}

	for(i=0;i<n;i++) delete [] xi[i]; delete [] xi;
	delete [] xit; delete [] ptt; delete [] pt; return;
}

int ncom;
float *pcom,*xicom,(*nrfunc)(int, float []);
			
void linmin(float p[], float xi[], int n, float &fret, 
	    float (*func)(int, float []))
/* given an n-dimensional point p[1..n] and an n-dimensional
   direction xi[1..n], moves and resets p to where the function
   func() takes on a minimum along the direction xi from p,
   and replaces xi by the actual vector displacement that p
   was moved. also returns as fret the value of func() at the
   returned location p. this is actually all accomplished by
   calling the routines mnbrak() and brent().*/
{
	extern int ncom;
	extern float *pcom,*xicom,(*nrfunc)(int, float []);
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;	pcom=new float[n]; xicom=new float[n];
	nrfunc=func;  // define the global variables

	for(j=0;j<n;j++)
		{
		 pcom[j]=p[j]; xicom[j]=xi[j];
		}

	ax=0.0; xx=1.0;	  // initial guess for brackets

	mnbrak(ax,xx,bx,fa,fx,fb,f1dim);

	fret=brent(ax,xx,bx,f1dim,0.01,xmin);

	for(j=0;j<n;j++)
		{
		 xi[j]*=xmin;
		 p[j]+=xi[j];
		}

	 delete [] xicom; delete [] pcom; return;
}

float f1dim(float x)  // this function must accompany linmin
{
	extern int ncom;
	extern float *pcom,*xicom,(*nrfunc)(int, float []);
	int j;
	float f,*xt;

	xt=new float[ncom];

	for(j=0;j<ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(ncom,xt);
	
	delete [] xt; return f;
}

void mnbrak(float &ax, float &bx, float &cx, float &fa, float &fb, float &fc,
            float (*func)(float))
/* given a function func(), and given distinct initial points ax and bx, 
   this routine searches in the downhill direction (defined by the function
   as evaluated at the initial points) and returns new points ax, bx, cx
   that bracket a minimum of the function. Also returned are the function
   value at the three points, fa, fb, and fc. */
{
	float GOLD=1.618034;
	float GLIMIT=100.0;
	float ulim,u,r,q,fu,dum;

	fa=(*func)(ax);
	fb=(*func)(bx);

	if(fb>fa)	// switch roles of a and b to go downhill from a to b
	{
	 SHIFT(dum,ax,bx,dum);
	 SHIFT(dum,fb,fa,dum);
	}

	cx=bx+GOLD*(bx-ax);   // first guess for c
	fc=(*func)(cx);

	while(fb>fc)
	{
	 r=(bx-ax)*(fb-fc);
	 q=(bx-cx)*(fb-fa);
	 u=bx-((bx-cx)*q-(bx-ax)*r)/(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
	 ulim=bx+GLIMIT*(cx-bx);

	 if((bx-u)*(u-cx)>0.0)  // parabolic u is between b and c, give a try
		{
		 fu=(*func)(u);
		 if(fu<fc)   // got a minimum between b and c
			{
			 ax=bx; bx=u; fa=fb; fb=fu;
			 return;
			}
		 else if(fu>fb)  // got a minimum between a and u
			{
			 cx=u; fc=fu;
			 return;
			}
		 u=cx+GOLD*(cx-bx);
		 fu=(*func)(u);
		}
	 else if((cx-u)*(u-ulim)>0.0)
		{
		 fu=(*func)(u);
		 if(fu<fc)
			{
			 SHIFT(bx,cx,u,cx+GOLD*(cx-bx));
			 SHIFT(fb,fc,fu,(*func)(u));
			}
		}
	 else if((u-ulim)*(ulim-cx)>=0.0)
		{
		 u=ulim;
		 fu=(*func)(u);
		}
	 else
		{
		 u=cx+GOLD*(cx-bx);
		 fu=(*func)(u);
		}

	 SHIFT(ax,bx,cx,u);	// eliminate the oldest point and continue 
	 SHIFT(fa,fb,fc,fu);
	}
}

float brent(float ax, float bx, float cx, float (*f)(float), 
	    float tol, float &xmin)
/* given a function f, and given a bracketing triplet of ax,bx,cx (such
   that bx is between ax and cx, and f(bx) is less than both f(ax) and
   f(cx)), this routine isolates the minimum to a fractional precision
   of about tol using Brent's method. the abscissa of the minimum is
   returned as xmin, and the routine itself returns the minimum 
   function value.
*/
{
	float CGOLD=0.3819660;
	int ITMAX=100;
	int iter;
	float a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	float e=0.0;

	a=(ax<cx?ax:cx);	// a and b must be in ascending order
	b=(ax>cx?ax:cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);	// initialization

	for(iter=1;iter<=ITMAX;iter++)	// the main loop
		{
		 xm=0.5*(a+b);
		 tol2=2.0*(tol1=tol*fabs(x)+TINY);
		 if(fabs(x-xm)<=(tol2-0.5*(b-a)))  // test for termination
			{
			 xmin=x; return fx;
			}
		 if(fabs(e)>tol1)
			{
			 r=(x-w)*(fx-fv);  // construct a trial parabolic fit
			 q=(x-v)*(fx-fw);
			 p=(x-v)*q-(x-w)*r;
			 q=2.0*(q-r);
			 if(q>0.0) p=-p;
			 q=fabs(q);
			 etemp=e;
			 e=d;
			 if(fabs(p)>=fabs(0.5*q*etemp)||p<=q*(a-x)||p>=q*(b-x))
				{
				 // the above conditions determine the 
				 // acceptability of the parabolic fit.
				 // then take the golden section step into
				 // the larger of the two segments.
				 d=CGOLD*(e=(x>=xm?a-x:b-x));
				}
			 else
				{
				 d=p/q;
				 u=x+d;
				 if(u-a<tol2||b-u<tol2) d=SIGN(tol1,xm-x);
				}
			}
		 else
			{
			 d=CGOLD*(e=(x>=xm?a-x:b-x));
			}
		 u=(fabs(d)>=tol1?x+d:x+SIGN(tol1,d));
		 fu=(*f)(u);  // this is the one function evaluation per round
		 if(fu<=fx)
			{
			 if(u>=x) a=x; else b=x;
			 SHIFT(v,w,x,u);
			 SHIFT(fv,fw,fx,fu);
			}
		 else
			{
			 if(u<x) a=u; else b=u;
			 if(fu<=fw||w==x)
				{
				 v=w;
				 w=u;
				 fv=fw;
				 fw=fu;
				}
			 else if(fu<=fv||v==x||v==w)
				{
				 v=u;
				 fv=fu;
				}
			}
		}

	xmin=x; return fx;
}

