/*******************************************************************************
Singular value decomposition program, svdcmp, from "Numerical Recipes in C"
(Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling,
and B.P. Flannery
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdarg.h>

#define NR_END 1
#define FREE_ARG char*
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))

double **dmatrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;
	/* allocate pointers to rows */
	m= malloc((nrow+NR_END)*sizeof(double*));
	m += NR_END;
	m -= nrl;
	/* allocate rows and set pointers to them */
	m[nrl]= malloc((nrow*ncol+NR_END)*sizeof(double));
	m[nrl] += NR_END;
	m[nrl] -= ncl;
	for(i=nrl+1;i<=nrh;i++) {
		m[i]=m[i-1]+ncol;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

void
printdmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
	int i,j;
	for (i=nrl;i<=nrh;i++) {
		printf("[ ");
		for (j=ncl;j<=nch;j++) {
			printf("%lf ", m[i][j]);
		}
		printf("]\n");
	}
}


double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;
	v=malloc((nh-nl+1+NR_END)*sizeof(double));
	return v-nl+NR_END;
}

void free_dvector(double *v, int nl, int nh)
/* free a double vector allocated with dvector() */
{
	free(v+nl-NR_END);
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
        free( m[nrl]+ncl-NR_END);
        free(m+nrl-NR_END);
}

double pythag(double a, double b)
/* compute (a2 + b2)^1/2 without destructive underflow or overflow */
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

void 
printdvector(double *a, int nl , int nh, double xstart, double dx)
{
	int i;
	for (i = nl; i < nh; i++) {
		printf("%lf\t%lf\n",xstart, a[i]);
		xstart += dx;
	} 
	
}

/******************************************************************************/
void svdcmp(double **a, int m, int n, double w[], double **v)
/*******************************************************************************
Given a matrix a[1..m][1..n], this routine computes its singular value
decomposition, A = U.W.VT.  The matrix U replaces a on output.  The diagonal
matrix of singular values W is output as a vector w[1..n].  The matrix V (not
the transpose VT) is output as v[1..n][1..n].
*******************************************************************************/
{
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=dvector(1,n);
	g=scale=anorm=0.0; /* Householder reduction to bidiagonal form */
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l;k<=n;k++) rv1[k]=a[i][k]/h;
				for (j=l;j<=m;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[j][k]*a[i][k];
					for (k=l;k<=n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l;k<=n;k++) a[i][k] *= scale;
			}
		}
		anorm = DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n;i>=1;i--) { /* Accumulation of right-hand transformations. */
		if (i < n) {
			if (g) {
				for (j=l;j<=n;j++) /* Double division to avoid possible underflow. */
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=l;k<=n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<=n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n);i>=1;i--) { /* Accumulation of left-hand transformations. */
		l=i+1;
		g=w[i];
		for (j=l;j<=n;j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<=n;j++) {
				for (s=0.0,k=l;k<=m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<=m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<=m;j++) a[j][i] *= g;
		} else for (j=i;j<=m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n;k>=1;k--) { /* Diagonalization of the bidiagonal form. */
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=1;l--) { /* Test for splitting. */
				nm=l-1; /* Note that rv1[1] is always zero. */
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0; /* Cancellation of rv1[l], if l > 1. */
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { /* Convergence. */
				if (z < 0.0) { /* Singular value is made nonnegative. */
					w[k] = -z;
					for (j=1;j<=n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 30) printf("no convergence in 30 svdcmp iterations");
			x=w[l]; /* Shift from bottom 2-by-2 minor. */
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0; /* Next QR transformation: */
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z; /* Rotation can be arbitrary if z = 0. */
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1;jj<=m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);
}

/* remember! before calling zero out small w vals */
void
svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[])
{
	/* solves A*x = B */
	int jj, j, i;
	double s, *tmp;
	
	tmp = dvector(1, n);
	for(j = 1; j <= n; j++) {
		s = 0.0;
		if(w[j]) {
			for (i = 1; i <= m; i++) {
				s += u[i][j]*b[i];
			}
			s /=w [j];
		}
		tmp[j] = s;
	}
	for(j = 1; j <= n; j++) {
		s = 0.0;
		for( jj = 1; jj <= n; jj++) {
			s += v[j][jj]*tmp[jj];
		}
		x[j] = s;
	}
	free_dvector(tmp, 1, n);
}

struct xval {
	double x;
	double val;
};

void 
usage() 
{
	fprintf(stderr, "./a.out file\n file contains phi values for a Mach number\n");
}

void
printxvalarray(struct xval *a, size_t sz)
{
	size_t i = 0;
	for(i = 0; i < sz; i++) {
		printf("%lf\t%lf\n", a[i].x, a[i].val);
	}
}

int
xval_sortbyx(const void *a, const void *b)
{
	struct xval xva, xvb;
	xva = *((struct xval*)a);
	xvb = *((struct xval*)b);
	if (xva.x == xvb.x) 
		return 0;
	if (xva.x < xvb.x)
		return -1;
	return 1;
	
}

/* 
 * maps a function of double to double onto the vals of a xval array. 
 * caller must free. returns NULL if failure.
 */

struct xval*
mapfuncxvals(struct xval *ar, size_t sz, double (*fp)(double)) 
{
	struct xval *nar;
	size_t i;
	if ((nar = malloc(sz * sizeof(struct xval))) == NULL) {
		fprintf(stderr, "failed to malloc in mapfuncxvals\n");
		return NULL;
	}
	for(i=0; i < sz; i++) {
		nar[i].val = (*fp)(ar[i].val);
		nar[i].x = ar[i].x;
	}
	return nar;
}

/* takes var args as doubles and makes the ndiag matrix have those args.
 * returns -1 on fail
 */
int
ndiag(double **m, int nrl, int nrh, int ncl, int nch, int n, ...)
{
	int nargs=0, bufsiz=10, fail=0, i=0, j;
	double arg;
	double *args, *argstemp;
	va_list ap;
	if((args = malloc(bufsiz * sizeof(double))) == NULL) {
		fprintf(stderr, "failed to malloc in ndiag\n");
		return -1;
	}
	va_start(ap, n);
	while(i < n) {
		arg = va_arg(ap, double);
		args[nargs] = arg;
		if(++nargs == bufsiz-1) { /* realloc */
			if((argstemp = reallocarray(args, 2*bufsiz, sizeof(double))) == NULL) {
				fprintf(stderr, "failed to reallocarray in ndiag\n");
				fail=1;
				goto exit;
			}	
			args = argstemp;
		}
		i++;
	}
	va_end(ap);
	if(nargs > (nch - ncl)) {
		fprintf(stderr, "too many arguments (%d) for matrix that has %d columns\n", nargs, (nch-ncl));
		fail=1;
		goto exit;
	}
	printf("nargs:%d nargs/2 is :%d\n",nargs, nargs/2);
	argstemp = args+nargs/2;
	for (i=nrl; i<=nrh; i++) {
		for (j = ncl; j<=nch+nargs; j++) {
			if (argstemp < args+nargs) {
				m[i][j] = *argstemp++;
			} else {
				if(i < nargs/2) {
					argstemp = args + nargs/2 - i;
				} else {
					argstemp = args;
				}
				if (i > nargs/2) {
					ncl++;
				}
				break;
			}
		}
	}

exit:
	free(args);
	return (fail==1) ? - 1 : 0;
}

/* 
 * maps a function centrally by using the nearby points. it does not operate on 
 * the endpoints. it has 4 args to be able to see the x values too.
 * caller must free. returns NULL if failure.
 */
struct xval*
mapcentralfuncxvals(struct xval *ar, size_t sz, double (*fp) (double, double, double, double))
{
	struct xval *nar;
	size_t i;
	if (sz < 3) { // it's a noop jack.
		return NULL;
	}
	if ((nar = malloc(sz * sizeof(struct xval))) == NULL) {
		fprintf(stderr, "failed to malloc in mapfuncxvals\n");
		return NULL;
	}
	for(i=1; i < sz-1; i++) {
		nar[i].val = (*fp)(ar[i-1].val, ar[i+1].val, ar[i-1].x, ar[i+1].x);
		nar[i].x = ar[i].x;
	}
	nar[0].x = ar[0].x;
	nar[sz-1].x = ar[sz-1].x;
	return nar;
}

double
centralderiv(double val1, double val2, double x1, double x2)
{
	return (val2 - val1)/(x2 - x1);
}

double
densityfromphi(double a)
{
	return 1/sqrt(1-(a/1.1)); //XXX: hardcoded M
}

double
velocityfromphi(double a)
{
	return 1.1 - sqrt((1.1 * 1.1) - 2*a);
}

double
phifromdensity(double a)
{
	return a/(1+ (0.1*0.1)); /* XXX: hardcoded k = 2 */
}

void
simloop(double *phivec, double *nvec, double *uvec, double **A, double **U, double *resdiagar, int imin, int imax, int tsteps, double minx, double maxx) 
{
	int numpoints, i, t;
	double xrange, dx, dt, tstart=0;
	double *phistart, *nstart, *ustart, *phitemp, *ntemp, *utemp;
	numpoints = imax - imin;
	xrange = maxx - minx;
	dx = xrange / (numpoints - 1);
	dt = 0.5 * dx;
	printf("SIMLOOP> dx:%lf dt:%lf\n", dx, dt);
	/* allocating temporary and working vectors */
	phistart = dvector(imin, imax);
	nstart = dvector(imin, imax);
	ustart = dvector (imin, imax);
	phitemp = dvector(imin, imax);
	ntemp = dvector(imin, imax);
	utemp = dvector(imin, imax);
	/* copy values from input into start */
	for(i=1; i<= imax; i++) {
		phistart[i] = phivec[i];
		nstart[i] = nvec[i];
		ustart[i] = uvec[i];
	}
	for(t=0; t<tsteps; t++) {
		printf("#phi tstep:%d simtime:%lf\n", t, tstart);
		printdvector(phistart, imin, imax, minx, dx);
		printf("#n tstep:%d simtime:%lf\n", t, tstart);
		printdvector(nstart, imin, imax, minx, dx);
		printf("#u tstep:%d simtime:%lf\n", t, tstart);
		printdvector(ustart, imin, imax, minx, dx);

		for(i=1; i<= imax; i++) { /* copy to temp space for calc */
			phitemp[i] = phistart[i];
			ntemp[i] = nstart[i];
			utemp[i] = ustart[i];
		}
		for(i=1; i<= imax; i++) {
			if(i < imax-1) { /*going forward so don't disturb the last elem */
				ntemp[i] = nstart[i] - (dt/dx) * (nstart[i+1]*ustart[i+1] - nstart[i]*ustart[i]);
				utemp[i] = ustart[i] - (dt/(2*dx)) * ((ustart[i+1] * ustart[i+1]) - (ustart[i] * ustart[i]) + (2 * phistart[i+1]) - (2 * phistart[i]));
			}
		}
		// now with n and u calculated , get phi
		svbksb(A, resdiagar, U, imax, imax, nstart, phitemp);
		for(i=1; i<= imax; i++) { /* copy back for next iter */
			phistart[i] = phitemp[i];
			nstart[i] = ntemp[i];
			ustart[i] = utemp[i];
		}
		/* fill in kvec to validate results. on a random point */		
		tstart += dt;
	}

	/*deallocating working vectors*/
	free_dvector(phistart, imin, imax);
	free_dvector(nstart, imin, imax);
	free_dvector(ustart, imin, imax);
	free_dvector(phitemp, imin, imax);
	free_dvector(ntemp, imin, imax);
	free_dvector(utemp, imin, imax);
}

struct xval *
genpulse(double a, double f, double k, double xmin, double dx, int numpoints)
{
	int i=0;
	struct xval *pulsedata;
	if ((pulsedata = malloc(numpoints * sizeof(struct xval))) == NULL) {
		fprintf(stderr, "failed to allocate pulsedata\n");
		return NULL;
	}
	for (i=0; i<numpoints; i++) {
		pulsedata[i].x = xmin;
		pulsedata[i].val = a+f*sin(k*xmin);	
		xmin += dx;
	}
	return pulsedata;
}

struct xval *
gensechpulse(double a, double b, double xmin, double xmax, int numpoints)
{
	int i=0;
	double xrange, dx;
	struct xval *pulsedata;
	xrange = xmax - xmin;
	dx = xrange / (numpoints - 1);
	if((pulsedata = malloc(numpoints * sizeof(struct xval))) == NULL) {
		fprintf(stderr, "failed to allocate pulsedata\n");
		return NULL;
	}
	for (i=0; i<numpoints; i++) {
		pulsedata[i].x = xmin;
		pulsedata[i].val = a*(1.0/cosh(b*xmin));
		xmin += dx;
	}
	return pulsedata;
}

int 
main(int argc, char *argv[])
{
	struct xval *phidata, *phitemp, *ndata, *udata;
	size_t bufsiz = 20, valsread = 0;
	FILE *fp;
	int fscanfret, qsortret, ndiagret, i;
	struct xval tmpval;
	double **A, **U;
	double *resdiagar, *phivector, *nvector, *uvector, dx;
	/*if (argc != 2) {
		usage();
		return EXIT_FAILURE;
	}
	if ((phidata = malloc(bufsiz*sizeof(struct xval))) == NULL) {
		return EXIT_FAILURE;
	}
	if ((fp = fopen(argv[1], "r")) == NULL) {
		fprintf(stderr, "failed to open %s\n", argv[1]);
		return EXIT_FAILURE;
	}
	while ((fscanfret = fscanf(fp, "%lf    %lf", &tmpval.x, &tmpval.val)) != EOF) {
		phidata[valsread].x = tmpval.x;
		phidata[valsread].val = tmpval.val;
		if (valsread == (bufsiz - 1)) {
			fprintf(stderr, "reallocating at element:%d\n", valsread);
			if ((phitemp = reallocarray(phidata, 2 * bufsiz, sizeof(struct xval))) == NULL) {
				fprintf(stderr, "failed to realloc %d bytes\n", 2 * bufsiz);
				return EXIT_FAILURE;
			}
			bufsiz *= 2;
			phidata = phitemp;
		}
		valsread++;
	}
	fclose(fp);
        */
	/*printf("read %d xvals from file %s\n", valsread, argv[1]);*/
	/*qsort(phidata, valsread, sizeof(struct xval), &xval_sortbyx);*/
	/*printxvalarray(phidata, valsread);*/
	ndata = gensechpulse(1,0.1,-40,40,100);
	valsread = 100;
	dx = ndata[1].x - ndata[0].x;
	printf("dx is :%lf\n", dx);
	phidata = mapfuncxvals(ndata, valsread, &phifromdensity);
	udata = mapfuncxvals(phidata, valsread, &velocityfromphi);
	phivector = dvector(1, valsread);
	nvector = dvector(1, valsread);
	uvector = dvector(1, valsread);
	for(i=1;i<=valsread;i++) {
		//phivector[i] = phidata[i-1].val;
		//uvector[i] = udata[i-1].val;
		phivector[i] = 0;
		uvector[i] = 1.1;
		//nvector[i] = (dx * dx)*(1.0 - ndata[i-1].val);
		nvector[i] = ndata[i-1].val;
	}
	A = dmatrix(1, valsread, 1, valsread);
	if ((ndiagret = ndiag(A, 1, valsread, 1, valsread, 3, 1.0, -2.0*(1.0+(dx * dx)), 1.0)) == -1 ) {
		fprintf(stderr, "error in ndiag\n");
		return EXIT_FAILURE;
	}
	printf("finished ndiag\n");
	/* do PBC */
	A[1][valsread] = 1.0;
	A[valsread][1] = 1.0;
	/*printdmatrix(A, 1, 10, 1, 10);*/
	U = dmatrix(1, valsread, 1, valsread);
	resdiagar = dvector(1,valsread);
	svdcmp(A, valsread, valsread, resdiagar, U);
	/*printf("after decomp\n");*/
	/*printdmatrix(A, 1, 10, 1, 10);*/
	/*printf("diags:\n");*/
	for(i=1;i<valsread;i++) { /* according to NR we have to zero out the small wjs. so i zero less than 1 */
		if (resdiagar[i] < 1) {
			resdiagar[i] = 0.0;
		}
		/*printf(" %lf ", resdiagar[i]);*/
	}
	simloop(phivector, nvector, uvector, A, U, resdiagar, 1, valsread, 100, -40, 40);
	free_dvector(resdiagar, 1, valsread);
	free_dvector(nvector, 1, valsread);
	free_dvector(phivector, 1, valsread);
	free_dmatrix(A, 1, valsread, 1, valsread);
	free_dmatrix(U, 1, valsread, 1, valsread);
	free(ndata);
	free(phidata);
	free(udata);
	return EXIT_SUCCESS;
}
