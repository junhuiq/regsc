/*********************************************************************
 * bcd.c
 * This program implements the BCD (Block-Coordinate Descent) algorithm.
 *     bcd(y,x,z,n,p,q,lambda,weight,XTol,maxIter,b,gamma);
 *  Output: 
 *      b: an n-by-p matrix
 *		gamma: a q-by-1 vector
 *  Inputs:
 *		y: lhs variable, an n-by-1 vector
 *      x: rhs variables, an n-by-p matrix
 *      z: rhs variables, an n-by-q matrix
 *		n: length of y (sample size)
 *		p: number of variables in x 
 *		q: number of variables in z
 *      lambda: tuning parameter, scalar
 *      weight: weighting vector, (n-1)-by-1
 *      XTol: error tolerance, scalar
 *      maxIter: max num of iterations, scalar integer
 *  Subroutines needed:
 *      useful.c useful.h
 *
 *  Junhui Qian
 *  Shanghai Jiao Tong University
 *  Dec 22, 2018
 *  junhuiq@gmail.com
 *
 ********************************************************************/
//#include "mex.h"   
//#include <matrix.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include "useful.h"


double func(double gam, const double *R, int p, const double lam, const double *g){
    int i;
    double *A, *q, sum;
    
	if(gam<0)
		return 1e+6;

	A   = (double *)malloc(p*p*sizeof(double));
	q   = (double *)malloc(p*sizeof(double));

    for(i=0;i<p*p;i++)
        A[i] = gam*R[i];
    for(i=0;i<p;i++)
        A[i*p+i] = A[i*p+i]+lam*lam/2;
    for(i=0;i<p;i++)
        q[i] = g[i];
    gaussj(A, p, q, 1);
    sum=0;
    for(i=0;i<p;i++)
        sum = sum + g[i]*q[i];
    free(A);
    free(q);
    return gam*(1-0.5*sum);
}

void bcd(const double *y, const double *x, const double *z, const int n, const int p, const int q, const double lambda, const double *weight, const double XTol, const int maxIter, double *b, double *gamma)
{
    double *A,*vp,*vp2,*err,*Rxx,*Sxx,*Rxy,*Sxy,*b0,*g; 
    double *beta,*beta1,*delta,*s1,*s2,*b1m,*bm,*v1,*v2,*v3,*np;
    double *Szz,*Sxz,*Rxz,*vq,*A1,*A2,*nq;
    int i,j,k,m,t,exitflag;
    double dTmp,lam,ax,bx,cx,fa,fb,fc,xmin,errn;
    
    ax = 1e-16; bx = 1e+6;
    
	A   = (double *)malloc(p*p*sizeof(double));
	vp  = (double *)malloc(p*sizeof(double));
	vp2   = (double *)malloc(p*sizeof(double));
	err = (double *)malloc((n*p)*sizeof(double));
	Rxx = (double *)malloc(n*p*p*sizeof(double));
	Sxx = (double *)malloc(n*p*p*sizeof(double));
	Rxy = (double *)malloc(n*p*sizeof(double));
	Sxy = (double *)malloc(n*p*sizeof(double));
	b0  = (double *)malloc(n*p*sizeof(double));
	g   = (double *)malloc(p*sizeof(double));

    
	beta = (double *)malloc(n*p*sizeof(double));
	beta1= (double *)malloc((n-1)*p*sizeof(double));
	delta= (double *)malloc(n*sizeof(double));
	s1   = (double *)malloc(n*sizeof(double));
	s2   = (double *)malloc(n*sizeof(double));
	b1m  = (double *)malloc((n-1)*p*sizeof(double));
	bm   = (double *)malloc(n*p*sizeof(double));
	v1   = (double *)malloc(p*sizeof(double));
	v2   = (double *)malloc(p*sizeof(double));
	v3   = (double *)malloc(p*sizeof(double));
	np   = (double *)malloc(n*p*sizeof(double));

    
    if(q>0){
		Szz = (double *)malloc(q*q*sizeof(double));
		Sxz = (double *)malloc(n*p*q*sizeof(double));
		Rxz = (double *)malloc(n*p*q*sizeof(double));
		vq = (double *)malloc(q*sizeof(double));
		A1 = (double *)malloc(p*q*sizeof(double));
		A2 = (double *)malloc(q*q*sizeof(double));
		nq = (double *)malloc(n*q*sizeof(double));
	}
    
    for(i=0;i<n;i++){
        getrow(x,n,p,1,i,v1);

        fillv(v2,p,y[i]);
        dotmult(v1,1,p,v2);
        insert_row(Rxy,n,p,1,i,v2);
        
        mtimes(v1,p,1,v1,p,A);
        insert_row(Rxx,n,p,p,i,A);
    }
    invshift(Rxx,n,p*p,Sxx); cumsum(Sxx,n,p*p,Rxx); invshift(Rxx,n,p*p,Sxx);
    invshift(Rxy,n,p,Sxy); cumsum(Sxy,n,p,Rxy); invshift(Rxy,n,p,Sxy);
    mtrans(x,n,p,np); // np = x'
    if(q>0){
        for(i=0;i<n;i++){
            getrow(x,n,p,1,i,v1); 
            getrow(z,n,q,1,i,vq); 
            mtimes(v1,p,1,vq,q,A1); // A1=x_i*z_i'
            insert_row(Rxz,n,q,p,i,A1);
        }
    	invshift(Rxz,n,p*q,Sxz); cumsum(Sxz,n,p*q,Rxz); invshift(Rxz,n,p*q,Sxz);
        mtrans(z,n,q,nq);   //nq = z';
        mtimes(nq,q,n,z,q,Szz); // Szz=z'*z;
    }
//    mexPrintf("%f,%f,%f,%f\n",Sxz[0],Sxz[n],Sxz[2*n],Sxz[3*n]);
	fillv(b0,n*p,0);
    b0[0]=0.01;
	fillv(b,n*p,0);
    errn=1e+10; 
    i = 0;
    while(errn>XTol && i<maxIter){
        i=i+1;
        for(j=1;j<=n;j++){
            if(j==1){
                /*
            if q > 0
                beta = cumsum(d1);
                gamma = Szz\(z'*(y-sum(x.*beta,2)));
                s2 = z*gamma;
            else
                s2 = zeros(n,1);
            end */
                if(q>0){
                    cumsum(b,n,p,beta);
                    dotmult(x,n,p,beta);
                    summ(beta,n,p,2,delta);
                    for(k=0;k<n;k++)
                        delta[k]=y[k]-delta[k];
                    mtimes(nq,q,n,delta,1,vq);//(z'*(y-sum(x.*beta,2)))
					for(m=0;m<q*q;m++)
						A2[m]=Szz[m];
                    gaussj(A2, q, vq, 1);
                    for(m=0;m<q;m++)
                        gamma[m]=vq[m];
                    mtimesv(z,n,q,vq,s2);
                }else{
                    fillv(s2,n,0);
                }
                /*
            beta1 = cumsum(d1(2:n,:));
            s1 = [0;sum(x(2:n,:).*beta1,2)];
            d1(1,:) = (reshape(Sxx(1,:),p,p) \ (x'*(y-s1-s2)))';
                 */
                for(k=0;k<n-1;k++){
                    for(m=0;m<p;m++){
                        b1m[m*(n-1)+k]=b[m*n+k+1];
                    }
                }
                cumsum(b1m,n-1,p,beta1);
                fillzero(s1,n);
                for(k=1;k<n;k++){
                    for(m=0;m<p;m++){
                        s1[k]=s1[k]+x[m*n+k]*beta1[m*(n-1)+k-1];
                    }
                }
                vcopy(y,n,delta);
                vsub(delta,n,s1);
                vsub(delta,n,s2);
                mtimes(np,p,n,delta,1,vp);
                getrow(Sxx,n,p,p,0,A);
                gaussj(A,p,vp,1);// (reshape(Sxx(1,:),p,p)\(x'*(y-s1-s2)))';
                for(k=0;k<p;k++)
                    b[k*n] = vp[k];
            }else{
                if(q>0){
                    getrow(Sxz,n,q,p,j-1,A1);
                    mtimesv(A1,p,q,gamma,v3);
                }else{
                    fillzero(v3,p);
                }
           /* if q > 0
                beta = cumsum([d1(1:nn-1,:); d0(nn:n,:)]);
                gamma = Szz\(z'*(y-sum(x.*beta,2)));
                s3 = reshape(Sxz(nn,:),p,q)*gamma;
            else
                s3 = zeros(p,1);
            end
*/
                sumr(b,n,p,0,j-1,vp2);
                getrow(x,n,p,1,j-1,vp);
				for(m=0;m<p;m++)
					v1[m]=dotprod(vp,p,vp2)*vp[m];

                fillv(v2,p,0);
                for(t=j+1;t<=n;t++){
                    sumr(b,n,p,0,j-1,vp2);
                    sumr(b0,n,p,j,t-j,vp);
                    vadd(vp2,p,vp);
                    
                    getrow(x,n,p,1,t-1,vp);
					for(m=0;m<p;m++)
						v2[m]=v2[m]+dotprod(vp,p,vp2)*vp[m];
                }
                getrow(Sxy,n,p,1,j-1,g);
                for(m=0;m<p;m++)
                    g[m]=v1[m]+v2[m]+v3[m]-g[m];
                getrow(Sxx,n,p,p,j-1,A);
                /*
            s1 =  x(nn,:)'*x(nn,:)*(sum(d1(1:nn-1,:),1))';
            s2 = zeros(p,1);
            for t=nn+1:n
                s2 = s2 + x(t,:)'*x(t,:)*(sum(d1(1:nn-1,:),1) + sum(d0(nn+1:t,:),1))';
            end
            g = -(Sxy(nn,:)'-s1-s2-s3);
*/
                /* Equivalent Matlab code:
                 * if norm(g)<=lambda*alpha*weight(nn-1)    d1(nn,:) = 0; */
                dTmp = norm2(g,p);
                if(dTmp <= lambda*weight[j-2]){
                    for(k=0;k<p;k++)
                        b[k*n+j-1]=0;
                }else{
                /* Equivalent Matlab code:
               *   R = reshape(Sxx(nn,:),p,p)+lambda*(1-alpha)/2*eye(p);
                 * [gam,fval,exitflag] = fminbnd(@myfun2,0,1e10,[],R,lambda*alpha*weight(nn-1),g); */
                    lam = lambda*weight[j-2];
                    mnbrak(&ax,&bx,&cx,&fa,&fb,&fc,A,p,lam,g,func);
                    exitflag=brent(ax,bx,cx,A,p,lam,g,func,XTol,&xmin);
                /* Equivalent Matlab code:
                 * d1(nn,:) = -gam*((gam*reshape(Rn(nn,:),p,p) + lambda^2/2*eye(p))\g(nn,:)')'; */    
                    if (exitflag==1){
                        for(k=0;k<p*p;k++)
                            A[k] = xmin*A[k];
                        for(k=0;k<p;k++)
                            A[k*p+k] = A[k*p+k]+lam*lam/2;
                        gaussj(A,p,g,1);
                        for(k=0;k<p;k++)
                            b[k*n+j-1]=-xmin*g[k];
                    }else{
                        for(k=0;k<p;k++)
                            b[k*n+j-1]=0;
                    }
                }
            }
        }
        for(k=0;k<n*p;k++)
            err[k]=b[k]-b0[k];
        errn = norm2(err,n*p)/norm2(b0,n*p);
        vcopy(b,n*p,b0);
//        mexPrintf("%d, %f\n",i,errn);
    }

    free(A);
    free(vp);
    free(vp2);
    free(err);
    free(Rxx);
    free(Sxx);
    free(Rxy);
    free(Sxy);
    free(b0);
    free(g);
    
    free(beta);
    free(beta1);
    free(delta);
    free(s1);
    free(s2);
    free(b1m);
    free(bm);
    free(v1);
    free(v2);
    free(v3);
    free(np);

	
	if(q>0){
        free(Szz);
        free(Sxz);
		free(Rxz);
        free(vq);
        free(A1);
        free(A2);
        free(nq);
    }
}

void bcd_C(const double *y, const double *x, const double *z, const int *n, const int *p, const int *q, const double *lambda, const double *weight, const double *XTol, const int *maxIter, double *b, double *gamma)
{
  bcd(y, x, z, *n, *p, *q, *lambda, weight, *XTol, *maxIter, b, gamma);
}
  
