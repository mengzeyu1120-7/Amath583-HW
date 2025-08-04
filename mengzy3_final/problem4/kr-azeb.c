// kenneth.roche@pnl.gov; k8r@uw.edu 
// az=b 
// input: dim(a)
// output: accept / reject based on residual norm ... inf-norm
// with Andrew, Jan, Will ... tinkering with numerics a little
// theta: cc -c -DSEQ -DSANITY kr-azeb.c ; cc -o xazeb-seq-sanity kr-azeb.o

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <complex.h>

#include <mpi.h>

// if I want to use the MKL builds then ...
#define MKL_Complex16 double complex
#include <mkl.h>
#include <mkl_scalapack.h>
#include <mkl_pblas.h>
//#include <scalapack.h>

#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>

//internal clocks for caliper dtime 
static double  t_gettimeofday;
static clock_t t_clock, t_clock1;
static struct  timeval s;

void b_t(void)
{
  gettimeofday( &s, NULL);
  t_clock = clock();
  t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec;
}

double e_t(int type)
{
  switch (type) {
  case 0:
    t_clock1 = clock();
    gettimeofday( &s, NULL);
    t_clock = t_clock1 - t_clock;
    t_gettimeofday = s.tv_sec + 1e-6 * s.tv_usec - t_gettimeofday;
    return t_gettimeofday;
  case 1:
    return t_gettimeofday;
  case 2:
    return t_clock / (double)CLOCKS_PER_SEC;
  }
  return t_gettimeofday;
}

int b_rts ( int argc , char ** argv ) 
{
  int i ;
  FILE * fp ;
  char fn[ FILENAME_MAX ] ; 
  /* time stamping and clock synchronization */
  struct tm * tp ;
  time_t t ;
  
  /* initialize the time_t structure */
  t = time( NULL ) ; tp = localtime( &t ) ;
  
  /* name and open the log file -no PATH listed */
  sprintf( fn , "roche-run-log-%d.log" , 1900 + tp->tm_year ) ;
  
  if ( ( fp = fopen( fn , "a" ) ) == NULL ) { /* note the lack of a PATH to the FILE */
    fprintf( stderr , "error: cannot fopen() file %s\n", fn ) ;
    return -1 ;
  }
  
  // time stamp run log
  fprintf( fp , "\nstart_time(%d/%d::%d:%d:%d) [" , 1 + tp->tm_mon , tp->tm_mday , tp->tm_hour , tp->tm_min , tp->tm_sec ) ;
  
  i = 0 ;
  while ( i < argc ) 
    {
      fprintf( fp , "%s " , argv[ i ] ) ;
      i++ ;
    }
  fprintf( fp , "] " ) ;
  
  /* close the log file */
  if ( fclose( fp ) == EOF ) 
    { /* returns 0 on success */
      fprintf( stderr , "error: cannot fclose() file '%s' .\n" , fn ) ;
      return( EXIT_FAILURE ) ;
    }
  return 1 ;
}

int e_rts ( void ) 
{
  FILE * fp ;
  char fn[ FILENAME_MAX ] ; 
  /* time stamping and clock synchronization */
  struct tm * tp ;
  time_t t ;
  
  /* initialize the time_t structure */
  t = time( NULL ) ; tp = localtime( &t ) ;
  
  /* name and open the log file -no PATH listed */
  sprintf( fn , "roche-run-log-%d.log" , 1900 + tp->tm_year ) ;
  
  if ( ( fp = fopen( fn , "a" ) ) == NULL ) 
    { /* note the lack of a PATH to the FILE */
      fprintf( stderr , "error: cannot fopen() file %s\n" , fn ) ;
      return( EXIT_FAILURE ) ;
    }
  
  t = time( NULL ) ;
  tp = localtime( &t ) ;
  fprintf( fp , "end_time(%d/%d::%d:%d:%d)\n" , 1 + tp->tm_mon , tp->tm_mday , tp->tm_hour , tp->tm_min , tp->tm_sec ) ;
  
  /* close the file */
  if ( fclose( fp ) == EOF ) 
    fprintf( stderr , "error: cannot fclose() file %s\n" , fn ) ;
  
  return 1 ;
}

void mepsilon( double * rv ) 
{ /* double precision machine epsilon */
  double dpval , dpval_ , dptmp , dpone ; 
  int i ; /* count the powers of 2 enumerated */
  
  i = 0 ;
  dpone = dptmp = 1. ;
  dpval = 0. ;
  while ( ( dpval - dpone ) != 0. ) 
    { /* the machine is asked to work within its own limitations */
      i++ ; /* bump the current power of two */
      dpval_ = dpval ; /* store the last value which has meaning 
			  in case you want to do something with it */
      dptmp /= 2. ; /* form a smaller and smaller number */
      dpval = dpone + dptmp ; /* add this number to one as a basis 
				 for comparison to the number one */
      if ( i > 64 ) //note, 1 / 2^64 is O(1.e-19) 
	break ; 
    }
#ifdef VERBOSE 
  printf( "stop iteration %d , estimate of machine 'single' precision ~ 1 / (2^%d)\n" , i - 1 , i - 1 ) ; 
#endif
  *rv = 2. * dptmp ;   
  return ; 
}

//for allocating linear array indices to processes in a block cyclic manner 
void get_blk_cyc( int ip , int npdim , int ma , int mblk , int * nele ) 
{ //virtual process(0,0) owns the first element(s)
  int my_ele, nele_, np;
  int srcproc, extra;

  srcproc = 0;
  //virtual process(0,0) owns the first element(s) 
  my_ele = ( npdim + ip - srcproc ) % npdim;
  nele_ = ma / mblk;
  np = ( nele_ / npdim ) * mblk;
  extra = nele_ % npdim;
  if ( my_ele < extra ) np += mblk; 
  else if ( my_ele == extra ) np += ma % mblk;
  *nele = np;
}

//for determining how many elements / data portion of the symbolic matrix the process at (ip,iq) of the pXq virtual rectangular process grid will manage in memory
void get_blk_cyc( int ip , int npdim , int ma , int mblk , int * nele ) ;
void get_mem_req_blk_cyc ( int ip , int iq , int np , int nq , int ma , int na , int mblk , int nblk , int * nip , int * niq ) 
{ // irows are traversed first in the loop of this map - this is known as a column major ordering 
  if ( ip >= 0 && iq >= 0 )
    {
      get_blk_cyc( ip , np , ma , mblk , nip ); //get rows
      get_blk_cyc( iq , nq , na , nblk , niq ); //get columns ...
#ifdef VERBOSE
      printf( "ip,iq = %d,%d\t nip, niq = %d %d\n" , ip , iq , *nip , *niq );
#endif
    }
  else 
    {
      *nip = -1;
      *niq = -1;
    }
}

void get_mem_req_blk_cyc ( int ip , int iq , int np , int nq , int ma , int na , int mblk , int nblk , int * nip , int * niq ) ;
long int fast_pmat_herm( int ip , int iq , int p , int q , double complex * a , int na , int mb )
{
  /*`
    subroutine ZLAGHE(integer N, integer K, double precision dimension( * ) D, complex*16 dimension( lda, * ) A, integer LDA, integer dimension( 4 ) ISEED, complex*16 dimension( * ) WORK, integer INFO ); ... returns hermitian matrix 	

    ZLANHE(),  returns the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value  of a complex hermitian matrix A.
  */
  
  long int nelements = 0L ;
  int gi , gj , li , lj , lit , ljt , iptmp , iqtmp , iptmpt , iqtmpt ;
  int nip , niq ;
  int seed = 11 ; 
  double complex ztmp = 0. + I * 0. ; 

  get_mem_req_blk_cyc ( ip , iq , p , q , na , na , mb , mb , &nip , &niq ) ;
  srand( seed ) ;
  for ( gj = 0 ; gj < na ; gj++ )
    { /* left to right over the global matrix columns */ 
      iqtmp = ( int ) floor( ( double ) ( gj / mb ) ) % q ;
      iptmpt = ( int ) floor( ( double ) ( gj / mb ) ) % p ;
      for ( gi = gj ; gi < na ; gi++ )
	{ /* lower triangular part explicitly generated including the diagonal */ 
	  iptmp = ( int ) floor( ( double ) ( gi / mb ) ) % p ;
	  iqtmpt = ( int ) floor( ( double ) ( gi / mb ) ) % q ;
	  ztmp = ( 0.5 - ( double ) rand() / ( double ) RAND_MAX ) + I * ( 0.5 - ( double ) rand() / ( double ) RAND_MAX ) ; 
	  if ( gi == gj )
	    ztmp = creal( ztmp ) ;
	  if ( ip == iptmp && iq == iqtmp )
	    { /* the lower triangular part */
	      nelements++ ;
	      li = ( int ) floor( ( double ) ( floor( ( double ) ( gi / mb ) ) / p ) ) * mb + gi % mb ;
	      lj = ( int ) floor( ( double ) ( floor( ( double ) ( gj / mb ) ) / q ) ) * mb + gj % mb ;
	      *( a + li + lj * nip ) = ztmp ;
	    }
	  if ( ip == iptmpt && iq == iqtmpt && gi != gj )
	    { 
	      nelements++ ;
	      lit = ( int ) floor( ( double ) ( floor( ( double ) ( gj / mb ) ) / p ) ) * mb + gj % mb ;
	      ljt = ( int ) floor( ( double ) ( floor( ( double ) ( gi / mb ) ) / q ) ) * mb + gi % mb ;
	      *( a + lit + ljt * nip ) = creal( ztmp ) + -1. * I * cimag( ztmp ) ; /* conjugate the ij term into the ji term */
	    }
	}
    }
  return nelements ; 
}

double cmatrix_infnorm( int ma, int na , double complex * a )
{ // inf_norm(A) is the largest row-sum of |A|
  // ma, na are row and column dimensions of A respectively
  int i, j;
  double dtmp, inf_norm_a;

  inf_norm_a = 0.;
  for ( i = 0 ; i < ma ; i++ )
    {
      dtmp = 0.;
      for ( j = 0 ; j < na ; j++ )
	dtmp+= cabs(a[i+j*ma]);
      if ( (dtmp-inf_norm_a) > 0.)
	inf_norm_a = dtmp;
    }
  return inf_norm_a;
}

void get_mem_req_blk_cyc ( int ip , int iq , int np , int nq , int ma , int na , int mblk , int nblk , int * nip , int * niq ) ;
void pcmatrix_infnorm ( int ip , int iq , int p , int q , double complex * a , int ma , int na , int mb , MPI_Comm comm , double *infn, double *onen )
{
  int Ione = 1 ;
  int iam , np ;
  int i , j , gi , gj , li , lj , nip , niq ;
  double norm ; 
  double * sv , * rv ;

  MPI_Comm_size( comm , &np ) ;
  MPI_Comm_rank( comm , &iam ) ;
  get_mem_req_blk_cyc ( ip , iq , p , q , ma , na , mb , mb , &nip , &niq ) ;
  
  // parallel computation of inf-norm of matrix A
  //pick the largest row sum, max(1<=j<=na) Sum_j |a_ij| 
  if ( ( sv = malloc( ma * sizeof( double ) ) ) == NULL ) fprintf( stderr, "error: cannot malloc sv_\n" ) ;
  if ( iam == 0 ) if ( ( rv = malloc( ma * sizeof( double ) ) ) == NULL ) fprintf( stderr, "error: cannot malloc rv_\n" ) ; //only root needs this
  for ( i = 0 ; i < ma ; i++ ) sv[ i ] = 0. ;
  if ( iam == 0 ) for ( i = 0 ; i < ma ; i++ ) rv[ i ] = 0. + I * 0. ; //zero the receive buffer
  for ( li = 0 ; li < nip ; li++ )
    {
      gi = ip * mb + ( int ) floor( ( double ) ( li / mb ) ) * p * mb + li % mb ; 
      //this is the local contribution to the global row sum
      for ( lj = 0 ; lj < niq ; lj++ ) sv[ gi ] += cabs( a[ li + lj * nip ] ) ;
    }
  MPI_Reduce( sv , rv , ma , MPI_DOUBLE , MPI_SUM , 0 , comm ) ;
  if ( iam == 0 ) 
    { //find the largest one 
      norm = 0. ;
      for ( i = 0 ; i < ma ; i++ )
	if ( rv[ i ] - norm > 0. )
	  norm = rv[ i ] ;
    }
  MPI_Bcast( &norm , Ione , MPI_DOUBLE , 0 , comm ) ;
  *infn = norm; 
  if ( iam == 0 ) free( rv ) ; 
  free( sv ) ;

  /* parallel computation of 1-norm of matrix A */
  //pick the largest col sum, max(1<=j<=na) Sum_j |a_ij| 
  if ( ( sv = malloc( na * sizeof( double ) ) ) == NULL )
    fprintf( stderr, "error: cannot malloc sv_\n" ) ;
  if ( iam == 0 )
    if ( ( rv = malloc( na * sizeof( double ) ) ) == NULL )
      fprintf( stderr, "error: cannot malloc rv_\n" ) ;
  for ( j = 0 ; j < na ; j++ ) sv[ j ] = 0. ;

  for ( lj = 0 ; lj < niq ; lj++ ) 
    {
      gj = iq * mb + ( int ) floor( ( double ) ( lj / mb ) ) * q * mb + lj % mb ; 
      for ( li = 0 ; li < nip ; li++ )
	sv[gj] += cabs( a[ li + lj * nip ] ) ;
    }
  MPI_Reduce( sv , rv , na , MPI_DOUBLE , MPI_SUM , 0 , comm ) ;
  if ( iam == 0 ) 
    { //find the largest one 
      norm = 0. ;
      for ( j = 0 ; j < na ; j++ )
	if ( rv[ j ] > norm )
	  norm = rv[ j ] ;
    }
  MPI_Bcast( &norm , Ione , MPI_DOUBLE , 0 , comm ) ;
  *onen = norm; 
  if ( iam == 0 ) free( rv ) ; 
  free( sv ) ;
  
  return;
}

//zgesv_( &ma , &nrhs , a , &lda , ipiv , z , &ldb , &info ) ;
void lu_solver(int n, double complex * a, double complex *bz )
{ //no pivoting
  //a is over written with LU
  //b is over written with z
  int i,j,k;

  //forward reduction 
  for (k=0;k<n-1;k++)
    for (i=k+1;i<n;i++)
      {
	a[i+n*k] /= a[k+n*k]; //l(i,k) = a(i,k) / a(k,k) 
	for (j=k+1;j<n;j++)
	  a[i+n*j] -= a[i+n*k] * a[k+n*j];
	bz[i] -= a[i+n*k]*bz[k];
      }

  //back substitution
  for (k=n-1;k>=0;k--)
    {
      for (i=k+1;i<n;i++)
	bz[k] -= a[k+n*i]*bz[i];
      bz[k] /= a[k+n*k];  
    }
  
  return;
}

//pivoting assumes that the maximum element in each row and column is the same order of magnitude
void plu_solver(int n, double complex * a, double complex *bz )
{ //with pivoting this time
  //a is over written with LU
  //b is over written with z
  int i,j,k,l;
  int ipiv;
  double zmax, ztmp;

  //forward reduction 
  for (k=0;k<n-1;k++)
    {
      //check for pivot
      zmax = 0.;
      for (l=k;l<n;l++)
	{//check the rows for the largest value
	  ztmp = cabs(a[l+n*k]);
	  if (ztmp-zmax>0.)
	    {
	      zmax = ztmp;
	      ipiv = l;
	    }
	}
      if (ipiv!=k)
	{//execute the swap
	}

      for (i=k+1;i<n;i++)
	{
	  a[i+n*k] /= a[k+n*k]; //l(i,k) = a(i,k) / a(k,k) 
	  for (j=k+1;j<n;j++)
	    a[i+n*j] -= a[i+n*k] * a[k+n*j];
	  bz[i] -= a[i+n*k]*bz[k];
	}
    }
  //form inverse of the pivot
  //back substitution
  for (k=n-1;k>=0;k--)
    {
      for (i=k+1;i<n;i++)
	bz[k] -= a[k+n*i]*bz[i];
      bz[k] /= a[k+n*k];  
    }
  
  return;
}

void hqr_solver(int n, double complex * a, double complex *bz )
{ //Householder QR solve of Az=b
  //A is overwritten with R
  //b is over written with z
  //QRz=b; Rz = Q*b; z = R^-1Q*b
  //Q*b is computed implicitly, Q is not formed explicitly

  int i,j,k;
  double complex *x, *vk;
  double complex sgnz; //sgn(z) = 1./cabs(z) * z 
  double xnorm, vknorm;
  
  assert(x=malloc(sizeof(double complex)*n));
  assert(vk=malloc(sizeof(double complex)*n));
  
  //left to right over columns of A
  for (k=0;k<n-1;k++) //last column is free
    {
      //copy to x, row k to n in column k of A
      memcpy( x , a+k*n+k , (n-k) * sizeof(double complex) ); 
      memcpy( vk , x , (n-k) * sizeof(double complex) ); 

      xnorm = 0.; vknorm = 0.;
      for (i=0;i<n-k;i++) //get the 2-norm of x
	vknorm += cabs(x[i])*cabs(x[i]);
      xnorm = sqrt(vknorm);

      sgnz = x[0] / cabs(x[0]);
      vk[0] += sgnz * xnorm; //assume + and redirect the first element
      if (cabs(x[0]-sgnz * xnorm) > cabs(vk[0]))
	vk[0] = x[0] - sgnz * xnorm; //correction 

      // normalize vk
      xnorm = 0.; vknorm = 0.;
      for (i=0;i<n-k;i++) //get the 2-norm of vk
	xnorm += cabs(vk[i])*cabs(vk[i]);
      vknorm = sqrt(xnorm);
      for (i=0;i<n-k;i++) 
	vk[i] /= vknorm; 

      //compose vk*A from (I - 2 vk vk*)A(k:n-1,k:n-1)
      for (i=0;i<n-k;i++) x[i]=0.+I*0.; //use x as a buffer 
      for (j=k;j<n;j++)
	for (i=k;i<n;i++)
	  x[j-k] += conj(vk[i-k]) * a[j*n+i];  
      //converting A to R here (upper triangular)
      for (j=k;j<n;j++)
	for (i=k;i<n;i++)
	  a[i+j*n]-=2.*vk[i-k]*x[j-k];      
#ifdef VERBOSE
      //should be zero below the diagonal
      for (i=0;i<n;i++)
	fprintf(stdout,"\t[k=%d] a[%d,%d]=%f + I %f\n",k,i,k,creal(a[i+k*n]),cimag(a[i+k*n]));
      fprintf(stdout,"\n");
#endif

      //implicit calculation of Q*b
      //compose vk*b from (I - 2 vk vk*)b(k:n-1)
      sgnz = 0. + I * 0.; //use this as a buffer again
      for (i=k;i<n;i++)
	sgnz += conj(vk[i-k]) * bz[i];  
      for (i=k;i<n;i++)
	bz[i]-=2.*vk[i-k]*sgnz;
    
    }

  //back substitution
  for (k=n-1;k>=0;k--)
    {
      for (i=k+1;i<n;i++)
	bz[k] -= a[k+n*i]*bz[i];
      bz[k] /= a[k+n*k];  
    }

  free(x);
  free(vk);
  
  return;
}

void Cblacs_pinfo( int * , int * ) ;
void Cblacs_setup( int * , int * ) ;
void Cblacs_get( int , int , int * ) ;
void Cblacs_gridinit( int * , char * , int , int ) ;
void Cblacs_gridinfo( int , int * , int * , int * , int * ) ;
void Cblacs_exit( int ) ;

void Cdgesd2d( int , int , int , double * , int , int , int ) ;
void Cdgerv2d( int , int , int , double * , int , int , int ) ;
void Cigebs2d( int ictxt , char * scope , char * top , int m , int n , int * A , int lda ) ; 
void Cigebr2d( int ictxt , char * scope , char * top , int m , int n , int * A , int lda , int rsrc , int csrc ) ;

//extern void zgesv( int * , int * , double complex * , int * , int * , double complex * , int * , int * ) ;
//void pzgesv( int * , int * , double complex * , int * , int * , int * , int * , double complex * , int * , int * , int * , int * ) ;
//void pzgemv( char * , int * , int * , double complex * , double complex * , int * , int * , int * , double complex * , int * , int * , int * , int * , double complex * , double complex * , int * , int * , int * , int * ) ;
/*
  void pzgesv_( int * , int * , double complex * , int * , int * , int * , int * , double complex * , int * , int * , int * , int * ) ;
  void pzgemv_( char * , int * , int * , double complex * , double complex * , int * , int * , int * , double complex * , int * , int * , int * , int * , double complex * , double complex * , int * , int * , int * , int * ) ;
*/

void b_t(void);
double e_t(int);
int b_rts (int ,char **);
int e_rts (void);
void mepsilon( double *);
void get_blk_cyc(int ,int ,int ,int ,int *);
void get_mem_req_blk_cyc (int ,int ,int ,int ,int ,int ,int ,int ,int * ,int *);
double cmatrix_infnorm( int ma, int na , double complex * a );
void pcmatrix_infnorm ( int ip , int iq , int p , int q , double complex * a , int ma , int na , int mb , MPI_Comm comm , double *infn, double *onen );
void lu_solver(int n, double complex * a, double complex *bz );
void hqr_solver(int n, double complex * a, double complex *bz );
void plu_solver(int n, double complex * a, double complex *bz );

int main( int argc , char ** argv ) 
{

  /* for the kernel */
  int i , j , li , lj , gi , gj ;
  int * ipiv ; /* the pivot array */
  int info , maxlocr , * Ibuf ; 
  double meps , dtmp , * sv , * rv ; 
  double a_infnorm_p , z_norm_p , residual_p , err_p ; /* set of MPI processes executes the problem in parallel */
  double complex ztmp , * a , * z , * b , * a_bak , *sz ;
  double complex zone = 1. + I * 0. ; double complex zmone = -1. + I * 0. ;

  int np , p , q ; /* np ~ p q , np := number of processes , p := number of process rows , q := number of process columns */
  int ip , iq , nip , niq , iptmp , iqtmp ; /* id (ip,iq) in (p,q) rectangular , virtual process grid owns (nip,niq) elements */
  int ma , na , mb , nb ; /* matrix dimensions [ma,na][na,nb]+[ma,nb] , block size */
  int iam , npmpi ; 
  int DESCA[ 9 ] , DESCB[ 9 ] ; /* array descriptors */
  char *b_order, *scope ;
  int ione = 1 , mone = -1 , zero = 0 ;
  int iam_blacs , ictxt , nprocs_blacs ;

  /* just recreate the same data sets for sequential and parallel cases */
  int isd_a = 17 ;
  int isd_b = 37 ;  

  // internal timing terms 
  double wall , cpu ;

  /* initialize the MPI communicator MPI_COMM_WORLD */
  MPI_Init( &argc , &argv ) ;
  MPI_Comm_size( MPI_COMM_WORLD , &npmpi ) ;
  MPI_Comm_rank( MPI_COMM_WORLD , &iam ) ;

  if ( iam == 0 ) /* let root timestamp the job */
    if ( b_rts( argc , argv ) < 0 ) 
      printf( "error: cannot update run log file\n" ) ;
  
  /* parse the command line */
  if ( argc < 6 )  
    { 
      fprintf( stderr , "usage: %s p q mb nb ma\n" , argv[ 0 ] ) ;
      return( EXIT_FAILURE ) ;
    }
  p = atoi( argv[ 1 ] ) ;
  q = atoi( argv[ 2 ] ) ;
  np = p * q ;
  if ( np < npmpi ) return( EXIT_FAILURE ) ;
  mb = atoi( argv[ 3 ] ) ;
  nb = atoi( argv[ 4 ] ) ;
  ma = atoi( argv[ 5 ] ) ;
  na = ma ; //work with well determined linear systems - square cases only 

  // compare zgesv_() to pzgesv_() on common small problems
  // see my note on long double versus double precision
  mepsilon(&meps); //calibrate the machine numerics

#ifdef SEQ
#ifdef SANITY
  double complex * z_s_cp, * z_s_cp_ ;
  if ( ( z_s_cp = malloc( sizeof( double complex ) * na ) ) == NULL ) return( EXIT_SUCCESS ) ;
  if ( ( z_s_cp_ = malloc( sizeof( double complex ) * na ) ) == NULL ) return( EXIT_SUCCESS ) ;
#endif
  double a_infnorm_s , z_norm_s , residual_s , err_s ; /* single MPI process executes the problem */
  int lda , ldb , nrhs ; 
  if ( iam == 0 ) //there are two versions ... with and without pivoting
    { //state the problem and solve it with a single MPI process using LAPACK
      printf("eps_machine = %e\n",meps);

      if ( ( a = malloc( sizeof( double complex ) * ma * na ) ) == NULL ) return( EXIT_SUCCESS ) ;
      if ( ( a_bak = malloc( sizeof( double complex ) * ma * na ) ) == NULL ) return( EXIT_SUCCESS ) ;
    
      /* generate matrix elements -no fixing the diagonal */
      srand( isd_a ) ;
      for ( i = 0 ; i < ma * na ; i++ )
	a[ i ] = .5 - ( double ) rand() / ( double ) RAND_MAX + I * ( .5 - ( double ) rand() / ( double ) RAND_MAX ) ;
#ifdef ILLCON
      for ( j = 0 ; j < na ; j++ )
        a[ma-1+j*ma] = creal(a[j*ma])+meps + I*(cimag(a[j*ma])-meps); 
#endif
      memcpy( a_bak , a , ma * na * sizeof(double complex) ); //A is overwritten with its factors, so make a copy

#ifdef SHOWA 
      //extract columns and print to screen 
      for ( j = 0 ; j < na ; j++ )
	{
	  for ( i = 0 ; i < ma ; i++ )
	    printf( "[%d]As[%d] = (%g, %g)\n" , j , i , creal( a[ i+j*ma ] ) , cimag( a[ i+j*ma ] ) ) ;
	  printf( "\n\n" ) ;
	}
#endif 

      //estimate the norm of A
      a_infnorm_s = cmatrix_infnorm( ma , na , a ) ;
      printf("(seq):\tN=dim(A)=\t%d\t||A||_inf=\t%f\n",na,a_infnorm_s); //square matrix

      //assuming now ma = na ... square matrices only
      if ( ( b = malloc( sizeof( double complex ) * ma ) ) == NULL ) return( EXIT_SUCCESS ) ;
      if ( ( z = malloc( sizeof( double complex ) * na ) ) == NULL ) return( EXIT_SUCCESS ) ;
  
      srand( isd_b ) ;
      for ( i = 0 ; i < ma ; i++ ) 
	b[ i ] = 0.5 - ( double ) rand() / ( double ) RAND_MAX + I * ( 0.5 - ( double ) rand() / ( double ) RAND_MAX ) ;
      /* we will keep b as is for the residual check and send z to 
	 the solver as the rhs since the rhs is 
	 overwritten on return w/ solution vector z from Az=b */
      memcpy( z, b , ma * sizeof(double complex) );
      
      b_t();
      //let's try the reference algorithm ... 
      lu_solver(na,a,z);
      wall = e_t(0); cpu  = e_t(2);
      fprintf(stdout,"(seq-LU):\tN=\t\t%d\tnp=\t%d ( %d x %d )\twall_t: %8.6f\tcpu_t: %8.6f\n",na,ione*ione,ione,ione,wall,cpu);

#ifdef SANITY
      memcpy( z_s_cp , z , na * sizeof(double complex) );
#endif

      // compute the induced inf-norm (az-b) ~ the residual
      residual_s  = 0. ;
      for ( i = 0 ; i < ma ; i++ ) 
	{ 
	  ztmp = 0. + I * 0. ;
	  for ( j = 0 ; j < na ; j++ ) ztmp += a_bak[ i + j * na ] * z[ j ] ; 
	  ztmp -= b[ i ] ; // ztmp ~ az - b 
	  if ( cabs( ztmp ) > residual_s ) residual_s = cabs( ztmp ) ;
	}
      
      // inf-norm( z ) 
      z_norm_s = 0. ;
      for ( i = 0 ; i < ma ; i++ ) if ( cabs( z[ i ] ) > z_norm_s ) z_norm_s = cabs( z[ i ] ) ;
      
      // error estimate
      err_s = residual_s / ( a_infnorm_s * z_norm_s * meps ) ; // O( na ) or smaller 
      fprintf(stdout, "(seq-LU)\tn= %d\t|A|= %e\t|z|= %e\t|res|= %e\tMeps= %e\terror= %e\n" , na , a_infnorm_s , z_norm_s , residual_s , meps , err_s ) ;

      // and AGAIN, AGAIN
      // this time by way of (H)QR -Householder reflectors
      memcpy( a, a_bak , ma * na * sizeof(double complex) ); //get A again
      memcpy( z, b , ma * sizeof(double complex) );
      b_t();
      hqr_solver(na, a, z );
      //void hqr_solver(int n, double complex * a, double complex *bz )
      wall = e_t(0); cpu  = e_t(2);
      fprintf(stdout,"(seq-HQR):\tN=\t\t%d\tnp=\t%d ( %d x %d )\twall_t: %8.6f\tcpu_t: %8.6f\n",na,ione*ione,ione,ione,wall,cpu);

      // compute the induced inf-norm (az-b) ~ the residual
      residual_s  = 0. ;
      for ( i = 0 ; i < ma ; i++ ) 
	{ 
	  ztmp = 0. + I * 0. ;
	  for ( j = 0 ; j < na ; j++ ) ztmp += a_bak[ i + j * na ] * z[ j ] ; 
	  ztmp -= b[ i ] ; // ztmp ~ az - b 
	  if ( cabs( ztmp ) > residual_s ) residual_s = cabs( ztmp ) ;
	}
      
      // inf-norm( z ) 
      z_norm_s = 0. ;
      for ( i = 0 ; i < ma ; i++ ) if ( cabs( z[ i ] ) > z_norm_s ) z_norm_s = cabs( z[ i ] ) ;
      
      // error estimate
      err_s = residual_s / ( a_infnorm_s * z_norm_s * meps ) ; // O( na ) or smaller 
      fprintf(stdout, "(seq-HQR)\tn= %d\t|A|= %e\t|z|= %e\t|res|= %e\tMeps= %e\terror= %e\n" , na , a_infnorm_s , z_norm_s , residual_s , meps , err_s ) ;

      // and AGAIN, AGAIN
      // now with pivoting ... 
      // just use LAPACK
      memcpy( a, a_bak , ma * na * sizeof(double complex) ); //get A again
      memcpy( z, b , ma * sizeof(double complex) ); //get b again

      // initialize some of the parameters for the call to LAPACK routine */
      if ( ( ipiv = malloc( ma * sizeof( int ) ) ) == NULL ) printf( "error: cannot malloc ipiv\n" ) ;
      /* 
	 ipiv is empty on entry and on return holds row indices indicating
	 row ( i ) of a was swapped with row ( ipiv[ i ] - 1 ) of a ( C convention ) 
	 .NOT. row( ipiv[ i ] ) as suggested on p.237 of the LAPACK manual 
      */

      nrhs = 1 ;
      lda = ldb = ma ; 
      
      b_t();

      //void zgesv(const long long *, const long long *, MKL_Complex16 *, const long long *, long long *, MKL_Complex16 *, const long long *, long long *)
      //info = LAPACKE_zgesv (LAPACK_COL_MAJOR , ma , nrhs , a , lda , ipiv , z , ldb );
      zgesv( &ma , &nrhs , a , &lda , ipiv , z , &ldb , &info ) ;
      
      wall = e_t(0); cpu  = e_t(2);
      if ( info < 0 ) printf( "error: routine zgesv()\t INFO=%d\n" , info ) ;
      fprintf(stdout,"(seq-PLU):\tN=\t\t%d\tnp=\t%d ( %d x %d )\twall_t: %8.6f\tcpu_t: %8.6f\n",na,ione*ione,ione,ione,wall,cpu);

      
#ifdef SANITY
      memcpy( z_s_cp_ , z , na * sizeof(double complex) );
#endif
      

      // compute the induced inf-norm (az-b) ~ the residual
      residual_s  = 0. ;
      for ( i = 0 ; i < ma ; i++ ) 
	{ 
	  ztmp = 0. + I * 0. ;
	  for ( j = 0 ; j < na ; j++ ) ztmp += a_bak[ i + j * na ] * z[ j ] ; 
	  ztmp -= b[ i ] ; // ztmp ~ az - b 
	  if ( cabs( ztmp ) > residual_s ) residual_s = cabs( ztmp ) ;
	}
      
      // inf-norm( z ) 
      z_norm_s = 0. ;
      for ( i = 0 ; i < ma ; i++ ) if ( cabs( z[ i ] ) > z_norm_s ) z_norm_s = cabs( z[ i ] ) ;
      
      // error estimate
      err_s = residual_s / ( a_infnorm_s * z_norm_s * meps ) ; // O( na ) or smaller 
      fprintf(stdout, "(seq-PLU)\tn= %d\t|A|= %e\t|z|= %e\t|res|= %e\tMeps= %e\terror= %e\n" , na , a_infnorm_s , z_norm_s , residual_s , meps , err_s ) ;

      free ( a ) ; free( a_bak );
      free ( z ) ;
      free( b ) ; 
      free( ipiv ) ;
    } /* end root study of the problem */
  
  MPI_Barrier( MPI_COMM_WORLD ) ;
#endif

  // distributed memory version of az=b , i.e. in parallel 

  /* initialize the BLACS grid */
  b_order = "R" ;
  scope = "All" ;
  Cblacs_pinfo( &iam_blacs , &nprocs_blacs ) ;
  if ( nprocs_blacs < 1 ) 
    Cblacs_setup( &iam_blacs , &nprocs_blacs ) ;
  Cblacs_get( mone , zero , &ictxt ) ; 
  Cblacs_gridinit( &ictxt , b_order , p , q ) ;  /* 'Row-Major' */
  Cblacs_gridinfo( ictxt , &p , &q , &ip , &iq ) ; /* ip,iq: the process row,column id */

  get_mem_req_blk_cyc( ip , iq , p , q , ma , na , mb , nb , &nip , &niq ) ;

  /* the array descriptors */

  /* A -the original matrix elements */
  DESCA[ 0 ] = 1 ;
  DESCA[ 1 ] = ictxt ;
  DESCA[ 2 ] = ma ;
  DESCA[ 3 ] = na ;
  DESCA[ 4 ] = mb ;
  DESCA[ 5 ] = nb ;
  DESCA[ 6 ] = 0 ; /* again, C vs F conventions ???? */
  DESCA[ 7 ] = 0 ; /* again, C vs F conventions ???? */
  DESCA[ 8 ] = nip ;

  /* allocate memory for the data and work arrays */
  assert( a = malloc( nip * niq * sizeof( double complex ) ) ) ;
  assert( a_bak = malloc( nip * niq * sizeof( double complex ) ) ) ;
  assert( sz = malloc( ma * sizeof( double complex ) ) ) ; //complex send buffer 

  // create and distribute same matrix elements -not very clever but OK 
  // to create the data -have to visit all intermediate values generated by the rng
  srand( isd_a ) ;
  for ( gj = 0 ; gj < na ; gj++ )
    {
      iqtmp = ( int ) floor( ( double ) ( gj / nb ) ) % q ;
      for ( gi = 0 ; gi < ma ; gi++ ) 
	{
	  iptmp = ( int ) floor( ( double ) ( gi / mb ) ) % p ;
	  ztmp = .5 - ( double ) rand() / ( double ) RAND_MAX + I * ( .5 - ( double ) rand() / ( double ) RAND_MAX ) ;
	  if ( ip == iptmp && iq == iqtmp )
	    { //in this case, get local offset, and copy
	      li = ( int ) floor( ( double ) ( floor( ( double ) ( gi / mb ) ) / p ) ) * mb + gi % mb ;
	      lj = ( int ) floor( ( double ) ( floor( ( double ) ( gj / nb ) ) / q ) ) * nb + gj % nb ;
	      a[ li + lj * nip ] = ztmp ;
	    }
	}
    }
  memcpy( a_bak , a , nip * niq * sizeof(double complex) ); //A is overwritten with its factors, so make a copy

#ifdef SHOWA 
  int k ;
  double complex *rz;
  assert( rz = malloc( ma * sizeof( double complex ) ) ) ;
  
  //extract columns and print to screen 
  for ( k = 0 ; k < na ; k ++ ) 
    { /* loop over the global column vectors */
      for ( i = 0 ; i < ma ; i++ )
	rz[ i ] = sz[ i ] = 0. + I * 0. ;
      for ( lj = 0 ; lj < niq ; lj++ )
	{ /* extract a single vector */
	  gj = iq * nb + ( int ) floor( ( double ) ( lj / nb ) ) * q * nb + lj % nb ; 
	  if ( gj == k ) 
	    {
	      for ( li = 0 ; li < nip ; li++ )
		{ /* form contribution to single global vector */
		  gi = ip * mb + ( int ) floor( ( double ) ( li / mb ) ) * p * mb + li % mb ; 
		  sz[ gi ] = a[ li + lj * nip ] ; 
		}
	    }
	}
      MPI_Barrier( MPI_COMM_WORLD ) ;

      /* form single vector by reduction and then broadcast it */
      MPI_Reduce( sz , rz , ma , MPI_DOUBLE_COMPLEX , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
      if ( iam == 0 )
	for ( i = 0 ; i < ma ; i++ )
	  printf( "[%d]A[%d] = (%g, %g)\n" , k , i , creal( rz[ i ] ) , cimag( rz[ i ] ) ) ;
      
      if ( iam == 0 )
	printf( "\n\n" ) ;

      MPI_Barrier( MPI_COMM_WORLD ) ;
    }

  MPI_Barrier( MPI_COMM_WORLD ) ;
  free( rz );
#endif 

  //estimate the norm of A
  a_infnorm_p = 0. ;
  pcmatrix_infnorm( ip, iq, p, q, a, ma, na, mb, MPI_COMM_WORLD, &a_infnorm_p , &dtmp );
  if (iam==0) fprintf(stdout, "(par):\tN=dim(A)=\t%d\t||A||_inf=\t%f\n",na,a_infnorm_p); //square matrix
  
  /* the vectors */
  DESCB[ 0 ] = 1 ;
  DESCB[ 1 ] = ictxt ;
  DESCB[ 2 ] = ma ;
  DESCB[ 3 ] = 1 ; /* nrhs */
  DESCB[ 4 ] = mb ;
  DESCB[ 5 ] = 1 ;
  DESCB[ 6 ] = 0 ; /* row source in the pgrid */
  DESCB[ 7 ] = 0 ; /* column source in the pgrid */
  DESCB[ 8 ] = nip ;

  // work with the same rhs as the LAPACK routine 
  get_mem_req_blk_cyc ( ip , iq , p , q , na , 1 , nb , nb , &nip , &niq ) ;
  if ( ( b = malloc( sizeof( double complex ) * nip * 1 ) ) == NULL ) return( EXIT_FAILURE ) ;
  if ( ( z = malloc( sizeof( double complex ) * nip * 1 ) ) == NULL ) return( EXIT_FAILURE ) ;
      
  srand( isd_b ) ;
  if ( iq == 0 )
    for ( gi = 0 ; gi < ma ; gi++ ) 
      {
	iptmp = ( int ) floor( ( double ) ( gi / mb ) ) % p ;
	ztmp = .5 - ( double ) rand() / ( double ) RAND_MAX + I * ( .5 - ( double ) rand() / ( double ) RAND_MAX ) ;
	if ( ip == iptmp )
	  { // get local offset, and copy
	    li = ( int ) floor( ( double ) ( floor( ( double ) ( gi / mb ) ) / p ) ) * mb + gi % mb ;
	    //li = ( ( gi / mb ) / p ) * mb + gi % mb ;
	    b[ li ] = ztmp ;
	  }
      }
  memcpy( z , b , nip * sizeof(double complex) );
  //if ( iq != 0 ) for ( i = 0 ; i < nip ; i++ ) b[ i ] = z[ i ] = 0. + I * 0. ;
  
  /* for the solve : get some memory for the local storage needs here */
  if ( ( Ibuf = malloc( ( nprocs_blacs ) * ( sizeof *Ibuf ) ) ) == NULL ) 
    {
      fprintf( stderr , "error: cannot malloc() for Ibuf\n" ) ;
      MPI_Finalize() ;
      return( EXIT_SUCCESS ) ;
    }
  maxlocr = 0 ; /* each node needs this guy for ipiv */
  MPI_Gather( &nip , ione , MPI_INT , Ibuf , ione , MPI_INT , 0 , MPI_COMM_WORLD ) ; 
  MPI_Bcast( Ibuf , nprocs_blacs , MPI_INT , 0 , MPI_COMM_WORLD ) ;
  for ( i = 0 ; i < nprocs_blacs ; i++ ) /* make sure we select the largest one */
    if ( Ibuf[ i ] > maxlocr ) maxlocr = Ibuf[ i ] ;
  free( Ibuf ) ;
  
  // get the memory for the pivot array and the local work arrays 
  if ( ( ipiv = malloc( ( maxlocr + nb ) * ( sizeof *ipiv ) ) ) == NULL ) 
    {
      fprintf(stderr,"error: cannot malloc() for ipiv\n");
      MPI_Finalize() ;
      return( EXIT_SUCCESS ) ;
    }

  if (iam == 0) b_t();
  //pzgesv_( &na , &ione , a , &ione , &ione , DESCA , ipiv , z , &ione , &ione , DESCB , &info ) ;
  pzgesv( &na , &ione , a , &ione , &ione , DESCA , ipiv , z , &ione , &ione , DESCB , &info ) ;
  if (iam == 0)
    {
      wall = e_t(0);
      cpu  = e_t(2);
      fprintf(stdout,"(par-PLU):\tN=\t\t%d\tnp=\t%d ( %d x %d )\tblk: %d\twall_t: %8.6f\tcpu_t: %8.6f\n",na,p*q,p,q,nb,wall,cpu);
    }
#ifdef SANITY
#ifdef VERBOSE
  double complex * rz ;
  assert( rz = malloc( ma * sizeof( double complex ) ) ) ; //complex send buffer 

  //kr - remove IPIV verbose
  /*
    if (iam == 0)
    for(i=0;i<maxlocr+nb;i++)
    printf("ipiv[%d] = %d\t(row %d of ipiv exchanged with ipiv[%d]=%d)\n",i,ipiv[i],i,i,ipiv[i]);
  */

  // get global (z) to compare to sequential (z_s_cp) from LAPACK
  for ( i = 0 ; i < ma ; i++ )
    {
      sz[ i ] = 0. + I * 0. ;
      if ( iam == 0 )
        rz[ i ] = 0. + I * 0. ;
    }
  if ( iq == 0 )
    for ( li = 0 ; li < nip ; li++ )
      {
	gi = ip * mb + ( li / mb ) * p * mb + li % mb ;
	sz[ gi ] = z[ li ];
      }
  //just use reduce ... lazy
  MPI_Reduce( sz , rz , ma , MPI_DOUBLE_COMPLEX , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
  if ( iam == 0 )
    { //compare the answers numerically -parallel one is now in rz
      printf("(no pivot: N=\t%d\n",ma);
      for ( i = 0 ; i < ma ; i++ )
	{
	  if ( (creal(z_s_cp[i])-creal(rz[i])>meps)|| (cimag(z_s_cp[i])-cimag(rz[i])>meps) )
	    printf("S[%d]=\t( %2.15f , %2.15f )\tP[%d]=\t( %2.15f , %2.15f )\t[%d]\t= ( %2.15f , %2.15f )\n",i,creal(z_s_cp[i]),cimag(z_s_cp[i]),i,creal(rz[i]),cimag(rz[i]),i,creal(z_s_cp[i])-creal(rz[i]),cimag(z_s_cp[i])-cimag(rz[i]));
	}
      free(z_s_cp);
      printf("(pivot: N=\t%d\n",ma);
      for ( i = 0 ; i < ma ; i++ )
	{
	  if ( (creal(z_s_cp_[i])-creal(rz[i])>meps)|| (cimag(z_s_cp_[i])-cimag(rz[i])>meps) )
	    printf("S[%d]=\t( %2.15f , %2.15f )\tP[%d]=\t( %2.15f , %2.15f )\t[%d]\t= ( %2.15f , %2.15f )\n",i,creal(z_s_cp_[i]),cimag(z_s_cp_[i]),i,creal(rz[i]),cimag(rz[i]),i,creal(z_s_cp_[i])-creal(rz[i]),cimag(z_s_cp_[i])-cimag(rz[i]));
	}
      free(z_s_cp_);
    }
  free( rz );
  MPI_Barrier( MPI_COMM_WORLD ) ;
#endif
#endif

  // form the residual vector 
  //void pzgemv_( TRANS, M, N, ALPHA, A, IA, JA, DESCA, X, IX, JX, DESCX, INCX, BETA, Y, IY, JY, DESCY, INCY )
  pzgemv( "N" , &ma , &na , &zone , a_bak , &ione , &ione , DESCA , z , &ione , &ione , DESCB , &ione , &zmone , b , &ione , &ione , DESCB , &ione ) ;
  
  // parallel computation of infnorm of vector residual 
  get_mem_req_blk_cyc( ip , iq , p , q , ma , na , mb , nb , &nip , &niq ) ;
  if ( ( sv = malloc( ma * sizeof( double ) ) ) == NULL ) fprintf( stderr, "error: cannot malloc sv_\n" ) ;
  if ( iam == 0 ) if ( ( rv = malloc( ma * sizeof( double ) ) ) == NULL ) fprintf( stderr, "error: cannot malloc rv_\n" ) ; //only root needs this
  for ( i = 0 ; i < ma ; i++ )
    {
      sv[ i ] = 0. ;
      if ( iam == 0 )
        rv[ i ] = 0. ;
    }
  if ( iq == 0 )
    for ( li = 0 ; li < nip ; li++ )
      {
	gi = ip * mb + ( li / mb ) * p * mb + li % mb ;
	sv[ gi ] = cabs( b[ li ] ) ;
      }
  /* now reduce this to root */
  MPI_Reduce( sv , rv , ma , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
  if( iam == 0 )
    { /* find the largest one */
      residual_p = 0. ;
      for ( i = 0 ; i < ma ; i++ )
        {
          if ( rv[ i ] > residual_p )
            residual_p = rv[ i ] ;
        }
    }
  // inform the other processes of the result
  MPI_Bcast( &residual_p , ione , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;

  // inf-norm(z)
  for ( i = 0 ; i < ma ; i++ )
    {
      sv[ i ] = 0. ;
      if ( iam == 0 )
        rv[ i ] = 0. ;
    }
  if ( iq == 0 )
    for ( li = 0 ; li < nip ; li++ )
      {
	gi = ip * mb + ( li / mb ) * p * mb + li % mb ;
	sv[ gi ] = cabs(z[li]);
      }
  /* now reduce this to root */
  MPI_Reduce( sv , rv , ma , MPI_DOUBLE , MPI_SUM , 0 , MPI_COMM_WORLD ) ;
  if( iam == 0 )
    { /* find the largest one */
      z_norm_p = 0. ;
      for ( i = 0 ; i < ma ; i++ )
        {
          if ( rv[ i ] > z_norm_p )
            z_norm_p = rv[ i ] ;
        }
    }
  /* ... and inform the other processes of the result */
  MPI_Bcast( &z_norm_p , ione , MPI_DOUBLE , 0 , MPI_COMM_WORLD ) ;

  free( sv ) ; 
  if( iam == 0 ) 
    free( rv ) ;

  /* estimate error */
  
  if ( iam == 0 )
    {
      err_p = residual_p / ( a_infnorm_p * z_norm_p * meps ) ; /* should be O( na ) */
      printf( "(par-PLU)\tn= %d\t|A|= %e\t|z|= %e\t|res|= %e\tMeps= %e\terror= %e\n" , na , a_infnorm_p , z_norm_p , residual_p , meps , err_p ) ;
    }

  /* pzasum_( &ma , ztmp , b , &ione , &ione , DESCB , &ione ) ; */
  free( a ) ; free( z ) ; free( b ) ; free( ipiv ) ; free( sz ) ;

  if ( iam == 0 )
    if ( e_rts( ) < 0 ) 
      printf( "error: cannot update run log file\n" ) ;
  
  MPI_Barrier( MPI_COMM_WORLD ) ;

#ifdef VERBOSE
  fprintf(stderr,"...[%d][%d]clean exit\n" , ip , iq ) ;
#endif

  MPI_Finalize( ) ;
  return( EXIT_SUCCESS ) ;
}

