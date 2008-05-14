#include <stdio.h>
#include <math.h>
#include <string.h>
#include "randlib.h"
#include "utils.h"
#include "vecutils.h"

#define min(a,b) ((a) <= (b) ? (a) : (b))
#define MXNCAT 100L
#define MAXHAP  MXNCAT
#define MAXROOK MXNCAT

typedef float  frvec[MAXROOK];
typedef float  fhvec[MAXHAP];
typedef int    irvec[MAXROOK];
typedef int    ihvec[MAXHAP];
typedef float  fhrmat[MAXHAP][MAXROOK];
typedef float  frhmat[MAXROOK][MAXHAP];
typedef long   ihrmat[MAXHAP][MAXROOK];
typedef int    irhmat[MAXROOK][MAXHAP];

/* 25 July 2000: revamp memory allocation etc. to avoid double-allocating:
 * single, etc. etc. */
extern void rdirich(float *prob, int ncat, float *results);

void gibbswrap(int *H, int *R, double *a, int *startiter, int *maxiter,
	       int * poolsamp, int * rooksamp, int *startfval, int *thin,
	       double *dprior, double *results,
	       int *outfile,
 	       char **outfn,
/* 	       char **randphrase, */
	       int *ranseed,
	       int *rptiter);
	  
int gibbs(int H, int R, float a, long startiter, long maxiter,
	  int * poolsamp, irhmat rooksamp, int startfval, int thin,
	  frvec fprior, double *resmat,
	  int outfile, char *outfn, int from_R, int rptiter);

void BUGSout(double * g, int H, int R, long tot, char *fn);

int main(int argc, char ** argv) {
  int H, R;
  float a;
  long startiter, maxiter;
  int rptiter=1000;
  int poolsamp[MAXHAP];
  irhmat rooksamp;
  int startfval;
  int thin;
  frvec fprior;
  double *resmat;

  int h,r;

  static char phrase[100];  /* random-number seed */
  static long is1,is2;

  /* no sanity checks: should check maxiter<startiter, (maxiter-startiter)/thin is
   * an integer, startfval <= R, etc., etc. */
  printf("Enter phrase (no spaces) for PRNG seed:\n");
  scanf("%s",phrase);
  phrtsd(phrase,&is1,&is2);
  setall(is1,is2);
  printf("Enter # haps, # rooks:\n");
  scanf("%d %d",&H,&R);
  /*   rooksamp=i2matrix(R,H); */
  printf("Enter hap samples from pooled pop (%d)\n",H);
  for (h=0; h<H; h++)
    scanf("%d",poolsamp+h);
  for (r=0; r<R; r++) {
    printf("Hap samples from rookery %d:\n",r+1);
    for (h=0; h<H; h++)
      scanf("%d",rooksamp[r]+h);
  }
  printf("startfval (0 for equal contribs, 1<=n<=R for biased contrib):\n");
  scanf("%d",&startfval);
  printf("burn-in, total, thinning factor:\n");
  scanf("%ld %ld %d",&startiter,&maxiter,&thin);
  a=1.0;
  gibbs(H,R,a,startiter,maxiter,poolsamp,rooksamp,startfval,thin,fprior,resmat,
	1,"turtle-gibbs",FALSE,rptiter);
  exit(0);
}

void gibbswrap(int *H, int *R, double *a, int *startiter, int *maxiter,
	       int * poolsamp, int * rooksamp, int *startfval, int *thin,
	       double *dprior, double *results,
	       int  *outfile,
 	       char **outfn,
/* 	       char **randphrase, */
	       int *ranseed,
	       int *rptiter) {

  /* get everything in double format from R */

  frvec fprior;
  irhmat rookmat;
  int h,r;
  int npts,nvars;
  static long is1,is2=234234;
  /* char *outfn="gibbsout"; */

  /* seed random-number generator */
  /*   phrtsd(*randphrase,&is1,i&s2); */
  is1 = (long)ranseed;
  setall(is1,is2);

  /*   fprintf(stderr,"PRNG seed: (%s) %ld %ld\n",*randphrase,is1,is2); */
  /*  if (*outfile==1) {
   * fprintf(stderr,"Printing output to file %s.{BOT,FRQ}\n",outfn[0]);
   *   }
   */

  for (h=0; h<(*H); h++) 
    for (r=0; r<(*R); r++)
      rookmat[r][h] = rooksamp[(*H)*r+h];

  /*   if (dprior[0]<0) */
  /*     fprior=NULL; */
  /*   else { */
  /*     fprintf(stderr,"allocating fprior\n"); */
  /*     fprior=farray(*R); */
    for (r=0; r<(*R); r++)
      fprior[r]=(float)dprior[r];
    /*   } */

  gibbs(*H,*R,(float)*a,(long)*startiter,(long)*maxiter,
	poolsamp,rookmat,*startfval,*thin,fprior,results,
	*outfile,*outfn,TRUE,*rptiter);

  npts =(*maxiter-*startiter)/(*thin);
  nvars = (*R)+(*H)*(*R);
  /*   for (i=0; i<npts; i++) */
  /*     for (j=0; j<nvars; j++) */
  /*       results[i+j*npts]=(double)(resmat[i][j]); */
      
  /* fprintf(stderr,"%f %f %f %f\n",
   * resmat[0][0],resmat[0][1],resmat[1][0],resmat[1][1]); */
}

int gibbs(int H, int R, float a, long startiter, long maxiter,
	  int * poolsamp, irhmat rooksamp, int startfval, int thin,
	  frvec fprior, double * results,
	  int outfile, char * outfn, int from_R, int rptiter) {

  static float genmul_kluge=0.9999F;
  int genmul_flag=0;
  float tmptot;

  int r,h,w0,w1;
  float sum;
  float totpool;
  float harmn;
  long it;
  char tmpbuf[1000];
  FILE * botfile, * frqfile;
  int print_out=FALSE;
  int save_out=TRUE;
  
  long npts;

  frhmat rookprior,  /* prior (Dirichlet) parameters for rookery haplotype freqs */
    rookfreqval,     /* current Gibbs sample of rookery haplotype freqs */
    rooktot2;        /* total Gibbs sample+prior+real sample rook-hap estimates */
  irvec rooktot;     /* total sample size for each rookery */
/*   frvec fval;         *//* current Gibbs sample for rookery contributions */
  float fval[MAXROOK]; /* test of funny display in gdb */
  fhvec ybar;        /* average overall haplotype freqs in rookery samples */
  fhvec poolfreq;    /* current computed Gibbs sample for pool haplotype freqs */
  fhrmat w;          /* temporary (transposed) Gibbs prob. that indiv with hap H comes from rookery R */
  ihrmat tmpmat;     /* Gibbs sample of numbers of each hap from each rookery */
  frvec rookcontrib, /* total (Gibbs) contribution of each rookery */
    rooktotp;        /* Gibbs contr + prior for rookery contribs */


  if (outfile>0) {
     print_out=TRUE;
     save_out=FALSE;
     /* open output files */
     strncpy(tmpbuf,outfn,100);
     strncat(tmpbuf,".bot",4);
     botfile = fopen(tmpbuf,"w"); 
     strncpy(tmpbuf,outfn,100);
     strncat(tmpbuf,".frq",4);
     frqfile = fopen(tmpbuf,"w"); 
  }
     
  if (H>MAXHAP) {
    fprintf(stderr,"# haplotypes (%d) exceeds maximum (%ld)\n",H,MAXHAP);
    return(1);
  }
  if (R>MAXROOK) {
    fprintf(stderr,"# rookeries (%d) exceeds maximum (%ld)\n",R,MAXROOK);
    return(1);
  }
  for (r=0; r<R; r++) {
    for (rooktot[r]=0, h=0; h<H; h++)
      rooktot[r] += rooksamp[r][h];
    if (rooktot[r]==0) {
      fprintf(stderr,"Can't do Gibbs with all-missing loci ...\n");
      return(2);
    }
  }
  for (h=0, totpool=0.0; h<H; h++)
    totpool += poolsamp[h];
  /*  calculate prior according to Pella and Masuda from harmonic mean:
      a parameter scales the strength of the prior */
  for (r=0, sum=0; r<R; r++) {
    sum += (float)1/rooktot[r];
  }
  harmn = 1/(sum/R);
  for (h=0; h<H; h++) {
    for (r=0,ybar[h]=0.0; r<R; r++)
      ybar[h] += (float)rooksamp[r][h]/rooktot[r];
    ybar[h] /= R;
  }
  for(h=0; h<H; h++)
    for (r=0; r<R; r++)
      rookprior[r][h] =  a*sqrt(harmn)*ybar[h];
/*   if (fprior==NULL) { */
/*     fprintf(stderr,"allocating fprior\n"); */
/*     fprior = farray(R); */
    /** default prior for contributions is EQUAL contrib from all rooks **/
  if (fprior[0]<0)
    for (r=0; r<R; r++)
      fprior[r]=(float)1/R;
  /* } */
  /* allocate results matrix if necessary */
  npts = (maxiter-startiter)/thin;
  if (!from_R && outfile<=0) {
    results = darray(npts*(R+H*R));
    fprintf(stderr,"Allocating space: from_R=%d\n",from_R);
  }
    /*     results = fmatrix((maxiter-startiter)/thin,R+H*R); */
    /*  dimnames(results) <- list(NULL,
	c(paste("contrib",dimnames(rooksamp)$rookery,sep="."),
	outer(1:H,1:R,function(x,y)paste("rookhap",dimnames(rooksamp)$rookery[y],
	dimnames(rooksamp)$haplotype[x],sep="."))))
    */
    /*   # set initial rookery freqs */
  for (r=0; r<R; r++)
    rdirich(rookprior[r],H,rookfreqval[r]);
  if (startfval<0)  /* ## use random start */
    rdirich(fprior,R,fval);
  else if (startfval==0)  /* ## equal-contribution start */
    for (r=0; r<R; r++)
      fval[r]=(float)1/R;
    else if (startfval<=R) {  /* ## start with 95% in one rookery, the rest evenly divided */
      for (r=0; r<R; r++) {
	if (r==startfval)
	  fval[r]=0.95;
	else
	  fval[r]=0.05/(R-1);
      }
    }
    else {
      fprintf(stderr,"startfval must be between 0 and R\n");
      return(3);
    }
    /*   ## pool contribs (f): vector, length R
	 ## rook haplotype freqs (h): matrix, H rows (haplotypes) x R cols (rookeries)
	 ## pool freqs (pool.freq): h %*% f, vector, length R
	 ## "val" indicates realized (Gibbs-sampler) value as opposed to Dirichlet params
    */
    for (it=0; it<maxiter; it++) {
      if (rptiter>0 && it % rptiter==0)
	fprintf(stderr,"it. %ld\n",it);
      for (h=0; h<H; h++) poolfreq[h]=0.0;
      for (h=0; h<H; h++)
	for (r=0; r<R; r++)
	  poolfreq[h]+=rookfreqval[r][h]*fval[r];
      /*  ## probability that an individual with hap H (row) comes from rookery R (column);
	  ## use R's "columns first" rule to calculate */
      for (h=0; h<H; h++)
	for (r=0; r<R; r++)
	  w[h][r] = rookfreqval[r][h]*fval[r]/poolfreq[h];
      /** genmul KLUGE: complains if sum(prob) > 0.99999 **/
      /* genmul_flag=0;
	 for (h=0; h<H; h++) {
	for (r=0,tmptot=0; r<R; r++)
	  tmptot += w[h][r];
	if (tmptot>genmul_kluge) {
	  genmul_flag=1;
	  if (tmptot>1.01) 
	    fprintf(stderr,"weird gibbs sum>1 (%f): adjusting down by %f\n",
		    tmptot,genmul_kluge/tmptot);
	  for (r=0; r<R; r++)
	    w[h][r] *= genmul_kluge/tmptot;
	}
	}*/
	  
      /*     ## take multinomial samples of each type ... */
      for (h=0; h<H; h++)
	/* 	genmul(poolsamp[h],&(w[h][0]),H,&(tmpmat[h][0])); */
	/*** DUMB DUMB DUMB DUMB! fixed bug -- was "H" instead of "R" for
	 *** number of categories  ... ***/
	genmul(poolsamp[h],&(w[h][0]),R,&(tmpmat[h][0]));
	/*     ## get posteriors for p (pool contribs, f) and Q (rook freqs, rookfreq) */
	/*     ## supposing we're doing things the easy way: */
	/*     ## posterior of p =  (pool sample plus any priors if desired) */
	/*     rookcontrib <- apply(tmpmat,2,sum) */
      for (r=0; r<R; r++) rookcontrib[r]=0.0;
      for (h=0; h<H; h++)
	for (r=0; r<R; r++)
	  rookcontrib[r] += tmpmat[h][r];
      for (r=0; r<R; r++)
	rooktotp[r] = rookcontrib[r]+fprior[r];
      rdirich(rooktotp,R,fval);
      /*     ## posterior of Q = (rookery sample + pool sample (known) + priors) */
      for (r=0; r<R; r++) {
	for (h=0; h<H; h++)
	  rooktot2[r][h] = rooksamp[r][h]+tmpmat[h][r]+rookprior[r][h];
	rdirich(rooktot2[r],H,rookfreqval[r]);
      }
      if (it>=startiter) {
	w0 = it-startiter;
	if (w0 % thin == 0) {
	  w1 = w0 / thin;
          if (print_out) {
  	     for (r=0; r<R; r++)
	        fprintf(botfile,"%1.10g ",fval[r]);
             fprintf(botfile,"\n");
	      for (r=0; r<R; r++)
	        for (h=0; h<H; h++)
	           fprintf(frqfile,"%1.10g ",rookfreqval[r][h]);
             fprintf(frqfile,"\n");
          } /* print output */
          if (save_out) {
  	     for (r=0; r<R; r++)
	        results[w1+r*npts] = (double)fval[r];
	      for (r=0; r<R; r++)
	        for (h=0; h<H; h++)
	           results[w1+(R+r*H+h)*npts] = (double)rookfreqval[r][h];
	  } /* save output */
	}  /* if thin */
      }  /* if beyond burn-in */
    } /* main Gibbs iteration loop */
    /* if (outfile) {
       fprintf(stderr,"check BUGS output file -- not tested since changes\n");
       BUGSout(results,H,R,(maxiter-startiter)/thin,outfn);
       } */
    if (print_out) {
      fclose(botfile);
      fclose(frqfile);
    }
    return(0);
}


void BUGSout(double * g, int H, int R, long tot, char *fn) {
  /*   ## take a gibbs run (nxp dataframe) and output to files in BUGS/CODA-compatible format */
  FILE *ind, *out;
  int p;
  int i;
  long j;

  ind=tryfile(fn,".ind",1,1);
  out=tryfile(fn,".out",1,2);

  p = R+H*R; /* number of columns */

  for (i=0; i<p; i++) {
    fprintf(ind,"var%d %ld %ld\n",i+1,i*tot+1,(i+1)*tot);
    for (j=0; j<tot; j++)
      fprintf(out,"%ld %f\n",j+1,g[j+i*tot]);
  }
  return;
}
