// Lasciate ogne speranza, voi ch’intrate
//
// Numerical stability is a joke to those scary people :C

#include <stdio.h>
#include <math.h>

#define RECURS_LEVEL_MAX  10
#define N_MAX             10

// pointer to function of one variable
typedef double (*Func1vFp)(double);

// example functions of one variable
double f_poly(double x) { // polynomial a[0] + a[1]x + ... + a[n]x^n
    double res = 2*pow(x,5)-4*pow(x,4)+3.5*x*x+1.35*x-6.25;
    return res;
}

double f_rat(double x) {
    double denom = (x-0.5)*(x-0.5)+0.01;
    return 1.0/denom;
}

double f_exp(double x) {
    double res = 2*x*exp(-1.5*x)-1;
    return res;
}

double f_trig(double x) {
    return x*tan(x);
}

// Quadratures

// rectangle rule, leftpoint
double quad_rect_left(Func1vFp f1, double a, double b, int n) {
    double dx = (b-a)/n;
    double sum = 0.0;
    for(int i = 0;i<n;i++){
        double x = a+dx*i;
        sum+=f1(x);
    }
    return dx*sum;
}

// rectangle rule, rightpoint
double quad_rect_right(Func1vFp f1, double a, double b, int n) {
    double dx = (b-a)/n;
    double sum = 0.0;
    for(int i = 0;i<n;i++){
        double x = a+dx*(i+1);
        sum+=f1(x);
    }
    return dx*sum;
}

// rectangle rule, midpoint
double quad_rect_mid(Func1vFp f1, double a, double b, int n) {
	double dx = (b-a)/n;
    double sum = 0.0;
    for(int i = 0;i<n;i++){
        double x = a+dx*i+dx/2;
        sum+=f1(x);
    }
    return dx*sum;
}

// trapezoidal rule
double quad_trap(Func1vFp func, double a, double b, int n) {
    double dx = (b-a)/n;
    double sum = func(a)+func(b);
    for(int i = 1;i<n;i++){
        sum+=func(a+dx*i)*2.0;
    }
    return dx*sum/2.0;
}

// Simpson's rule
double quad_simpson(Func1vFp f, double a, double b, int n) {
    double dx = (b-a)/n;
    double sum = f(a)-f(b);

    for(int i = 0;i<n;i++){
        sum+= 2*f(a+dx*(i+1))+4*f(a+dx*i+dx/2);
    }

    return sum*dx/6;
}

// pointer to quadrature function
typedef double (*QuadratureFp)(Func1vFp,double, double,int);

// array of pointers to integrand functions
Func1vFp func_tab[] = { f_poly, f_rat, f_trig, f_exp };

// array of pointers to quadrature functions
QuadratureFp quad_tab[] = {
	quad_rect_left, quad_rect_right, quad_rect_mid, quad_trap, quad_simpson };

// calls 'quad_no' quadrature function for 'fun_no' integrand function
// on interval [a, b] and n subintervals
double quad_select(int fun_no, int quad_no, double a, double b, int n) {
    return (quad_tab[quad_no](func_tab[fun_no],a,b,n));
}

// adaptive algorithm
double recurs(Func1vFp f, double a, double b, double S, double delta, QuadratureFp quad, int level) {
    double c = (a+b)/2;
    double S_1 = quad(f,a,c,1);
    double S_2 = quad(f,c,b,1);
    if(fabs((S_1+S_2)-S)<delta){
        return S_1+S_2;
    }
    level++;
    if(level>RECURS_LEVEL_MAX){
        return NAN;
    }
    return recurs(f,a,c,S_1,delta/2,quad,level)+recurs(f,c,b,S_2,delta/2,quad,level);
}

// initialization for adaptive algorithm
double init_recurs(Func1vFp f, double a, double b, double delta, QuadratureFp quad) {
    int level = 0;
    double S = quad(f,a,b,1);
    return recurs(f,a,b,S,delta,quad,level);
}

// double integrals

// pointer to function of two variables
typedef double (*Func2vFp)(double,double);

double func2v_2(double x, double y) {
	return 2 - x*x - y*y*y;
}

// sample functions to define the normal domain

double lower_bound2(double x) {
	return 0.7*exp(-2*x*x);
}
double upper_bound2(double x) {
	return sin(10*x);
}

// rectangle rule (leftpoint) - double integral over rectangular domain
double dbl_integr(Func2vFp f, double x1, double x2, int nx, double y1, double y2, int ny) {
    double dx = (x2-x1)/nx;
    double dy = (y2-y1)/ny;
    double sum = 0.0;

    for(int y = 0;y<ny;y++){
        for(int x = 0;x<nx;x++){
            sum+=f((x1+dx*x),(y1+y*dy));
        }
    }
    return sum*dx*dy;
}

// rectangle rule (midpoint) - double integral over normal domain with respect to the x-axis
double dbl_integr_normal_1(Func2vFp f, double x1, double x2, int nx, double hy,
						   Func1vFp fg, Func1vFp fh) {
    double sum = 0.0;
    double dx = (x2-x1)/nx;
    for(int i = 0;i<nx;i++){
        double x = x1+i*dx+dx/2;
        double g = fg(x);
        double h = fh(x);
        int ny = (int)round((h-g)/hy);
        double dy = (h-g)/ny;
        for(int j = 0;j<ny;j++){
            double y = g+dy*j+dy/2;
            sum+=f(x,y)*dy;
        }
    }
    return sum*dx;
}

// rectangle rule (leftpoint) - double integral over multiple normal
// domains with respect to the x-axis
double dbl_integr_normal_n(Func2vFp f, double x1, double x2, int nx, double y1, double y2,
		int ny, Func1vFp fg, Func1vFp fh)  {
    double dx = (x2-x1)/nx;
    double dy = (y2-y1)/ny;
    double sum = 0.0;
    double x = x1;
    for(int i = 0;i<nx;i++){
        double gx = fg(x);
        double hx = fh(x);
        double up = fmin(y2,hx);
        double low = fmax(y1,gx);

        int new_ny = ceil((up-low)/dy);
        double new_dy = (up-low)/(double)new_ny;
        double y = fmax(y1,gx);
        for(int j = 0;j<new_ny;j++){
            sum+=f(x,y)*new_dy;
            y+=new_dy;
        }
        x+=dx;
    }
    return sum*dx;
}

// multiple quadratures

typedef double (*FuncNvFp)(const double*, int);
typedef int (*BoundNvFp)(const double*, int);

// sample function of three variables
double func3v(const double v[], int n) {
	return v[0] - v[1] + 2*v[2];
}

// sample predicate in 3D
int bound3v(const double v[], int n) { // a cylinder
	return v[0] > 0 && v[0] < 0.5 && v[1]*v[1] + (v[2]-1)*(v[2]-1) < 1;
}

// sample function of n variables
double funcNv(const double v[], int n) {
	double fv = 1.;
	for (int i = 1; i < n; ++i) {
		fv += sin(i*v[i]);
	}
	return fv;
}
// sample n-dimensional predicate (n-dim hypersphere)
int boundNv(const double v[], int n) {
	double r = 0.0;
	for (int i = 0; i < n; ++i) r += (v[i]-1)*(v[i]-1);
	return r <= 1.;
}

// multiple integrals over a cuboid with predicate (if boundary != NULL)
// rectangular rule (rightpoint)
double trpl_quad_rect(FuncNvFp f, double variable_lim[][2], const int tn[], BoundNvFp boundary) {
    double dx = (variable_lim[0][1]-variable_lim[0][0])/tn[0];
    double dy = (variable_lim[1][1]-variable_lim[1][0])/tn[1];
    double dz = (variable_lim[2][1]-variable_lim[2][0])/tn[2];
    double sum=0.0;
    double x = variable_lim[0][0]+dx;
    for(int i = 0;i<tn[0];i++){
        double y = variable_lim[1][0]+dy;
        for(int j = 0;j<tn[1];j++){
            double z = variable_lim[2][0]+dz;
            for(int k = 0;k<tn[2];k++){
                double v[]={x,y,z};

                if(boundary==NULL||boundary(v,3)==1) {
                    sum+=f(v,3);
                }
                z=z+dz;
            }
            y=y+dy;
        }
        x=x+dx;
    }
    return sum*dx*dy*dz;
}

// multiple integrals over a n-dim hypercuboid with predicate (if boundary != NULL)
// rectangular rule (midpoint)
void recur_quad_rect_mid(double *p_sum, FuncNvFp f, int variable_no, double t_variable[],
         double variable_lim[][2], const int tn[], int level, BoundNvFp boundary) {
    double h = (variable_lim[level][1]-variable_lim[level][0])/tn[level];
    double x = variable_lim[level][0]+h/2;
    double result = 0.0;

    for(int i = 0;i<tn[level];i++){
        t_variable[level]=x;
        if(level < variable_no-1){
            recur_quad_rect_mid(&result,f,variable_no,t_variable,variable_lim,tn,level+1,boundary);
        }
        else if(boundary==NULL || boundary(t_variable,variable_no)==1){
            result += f(t_variable,variable_no);
        }
        x+=h;
    }
    *p_sum += result * h;
}

int main(void) {
	int to_do, n, nx, ny, integrand_fun_no, method_no, flag;
	int no_funcs = sizeof(func_tab) / sizeof(Func1vFp);
	int no_quads = sizeof(quad_tab) / sizeof(QuadratureFp);
	double a, b, x1, x2, y1, y2, hy, sum, delta;
	double t_variable[N_MAX], variable_lim[N_MAX][2];
	int tn[N_MAX];

	scanf("%d", &to_do);
	switch (to_do) {
		case 1: // loop over quadratures and integrands
			scanf("%lf %lf %d", &a, &b, &n);
			for(int q = 0; q < no_quads; ++q) {
				for(int f = 0; f < no_funcs; ++f) {
					printf("%.5f ",quad_select(f, q, a, b, n));
				}
				printf("\n");
			}
			break;

		case 2: // adaptive algorithm
			scanf("%d %d",&integrand_fun_no,&method_no);
			scanf("%lf %lf %lf", &a, &b, &delta);
			printf("%.5f\n", init_recurs(func_tab[integrand_fun_no],a,b,delta,quad_tab[method_no]));
			break;
		case 3: // double integral over a rectangle
			scanf("%lf %lf %d", &x1, &x2, &nx);
			scanf("%lf %lf %d", &y1, &y2, &ny);
			printf("%.5f\n", dbl_integr(func2v_2, x1, x2, nx, y1, y2, ny));
			break;
		case 4: // double integral over normal domain
			scanf("%lf %lf %d", &x1, &x2, &nx);
			scanf("%lf", &hy);
			printf("%.5f\n", dbl_integr_normal_1(func2v_2, x1, x2, nx, hy, lower_bound2, upper_bound2));
			break;
		case 5: // double integral over multiple normal domains
			scanf("%lf %lf %d", &x1, &x2, &nx);
			scanf("%lf %lf %d", &y1, &y2, &ny);
			printf("%.5f\n",dbl_integr_normal_n(func2v_2, x1, x2, nx, y1, y2, ny, lower_bound2, upper_bound2));
			break;
		case 6: // triple integral over a cuboid
			scanf("%lf %lf %d", &variable_lim[0][0], &variable_lim[0][1], tn);
			scanf("%lf %lf %d", &variable_lim[1][0], &variable_lim[1][1], tn+1);
			scanf("%lf %lf %d", &variable_lim[2][0], &variable_lim[2][1], tn+2);
			scanf("%d", &flag);
			printf("%.5f\n", trpl_quad_rect(func3v, variable_lim, tn, flag ? bound3v : NULL));
			break;
		case 7: // multiple integral over hyper-cuboid
			scanf("%d", &n);
			if(n > N_MAX) break;
			for(int i = 0; i < n; ++i) {
				scanf("%lf %lf %d", &variable_lim[i][0], &variable_lim[i][1], tn+i);
			}
			scanf("%d", &flag);
			sum = 0.;
			recur_quad_rect_mid(&sum, funcNv, n, t_variable, variable_lim, tn, 0, flag ? boundNv : NULL);
			printf("%.5f\n", sum);
			break;
		default:
			printf("Nothing to do for %d\n", to_do);
			break;
	}

	return 0;
}
