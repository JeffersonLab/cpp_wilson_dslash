#include <cpp_dslash_parscalar.h>

#include <cstdlib>
using namespace std;

#include <cache.h>

#include <cpp_dslash_parscalar_utils_64bit.h>

#include <tables_parscalar.h>
#include <dispatch_parscalar.h>

#include <cpp_dslash_parscalar_types.h>
using namespace CPlusPlusWilsonDslash::DslashParscalar64BitTypes;

#include <shift_table_parscalar.h> 

namespace CPlusPlusWilsonDslash {

  namespace DslashParscalar64Bit {
  
void decomp_plus(size_t lo,size_t hi, int id, const void *ptr) /*need to fix decomp_minus */
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;


  HalfSpinor *s3;
  FourSpinor* sp ALIGN;

  for (ix1=low;ix1<high;ix1++) {
    int thissite = tab->siteTable(ix1);

    sp=&psi[thissite];
    /******************************* direction +0 *********************************/	   
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_minus(*sp, *s3);
    /******************************* direction +1 *********************************/
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_minus(*sp, *s3);
    
    /******************************* direction +2 *********************************/
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_minus(*sp, *s3);

    /******************************* direction +3 *********************************/
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
    decomp_gamma3_minus(*sp, *s3);

  }
}


/* the basic operations in this routine include loading a spinor, doing 
 * the spin projection, and multiplying the halfspinor by the appropriate 
 * gauge field, and saving the resulting halfspinor to a lattice temporary */

/* need gauge fields on opposite cb */
void decomp_hvv_plus(size_t lo,size_t hi, int id, const void *ptr)
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix *um ALIGN;

  HalfSpinor *s3 ALIGN;

  FourSpinor *sm ALIGN; 

  /*  printf("ID: %d low=%d high=%d: DecompHvvPlus\n", id, low, high);*/

  for (ix1=low;ix1<high;ix1++) 
  {
    int thissite = tab->siteTable(ix1);

    /* Spinor to project*/
    sm=&psi[thissite];

    /******************************* direction -1 *********************************/
    um=&gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_plus(*sm, *um, *s3);
    
    /******************************* direction -1 *********************************/
    um=&gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_plus(*sm, *um, *s3);

    /******************************* direction -2 *********************************/
    um=&gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_plus(*sm, *um, *s3);

    /******************************* direction -3 *********************************/
    um=&gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
    decomp_hvv_gamma3_plus(*sm, *um, *s3);
  }
}
/***************end of decomp_hvv****************/


/* the basic operations in this routine include loading the halfspinor 
 * from memory, multiplying it by the appropriate gauge field, doing the 
 * spin reconstruction, and summing over directions, and saving the partial 
 * sum over directions */

void mvv_recons_plus(size_t lo,size_t hi, int id, const void *ptr)
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr;
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();
 
  int cb = a->cb; 

  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;
  GaugeMatrix *up ALIGN;
  HalfSpinor *s3 ALIGN, *s4 ALIGN;
  FourSpinor part_sum ALIGN, *result ALIGN;

  /* printf("ID: %d low=%d high=%d: MvvReconsPlus\n", id, low, high); */
 


  for (ix1=low;ix1<high;ix1++) {
    int thissite=tab->siteTable(ix1);
    result=&psi[thissite];

    up=&gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_plus(*s3, *up, part_sum);    

    up=&gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_plus_add(*s3, *up, part_sum);    

    up=&gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
    mvv_recons_gamma2_plus_add(*s3, *up, part_sum);    

    up=&gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
    mvv_recons_gamma3_plus_add_store(*s3, *up, part_sum, *result);    

  }
}


   
/* this routine takes the partial sum from mvv_recons() and loops 
 * over the output spin components, 2 at a time doing a sum over directions 
 * for each set, accumulating in xmm0-2 and loading the halfspinor 
 * temporaries into xmm3-5 */

void recons_plus(size_t lo,size_t hi, int id, const void *ptr )	
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  HalfSpinor *hs0 ALIGN , *hs1 ALIGN, *hs2 ALIGN, *hs3 ALIGN;
  FourSpinor *rn ALIGN;

  /*  printf("ID: %d low=%d high=%d: ReconsPlus\n", id, low, high); */
  for (ix1=low;ix1<high;ix1++) 
  {
    int thissite = tab->siteTable(ix1);
    rn=&psi[thissite];

    /* first spin component of result */
    hs0 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,0);
    hs1 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,1);	  
    hs2 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
    hs3 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,3);
    recons_4dir_plus(*hs0, *hs1, *hs2, *hs3, *rn);
  }
}


/*************** now for isign corresponding to -1  ****************************************/

void decomp_minus(size_t lo,size_t hi, int id, const void *ptr ) /*need to fix decomp_minus */
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  HalfSpinor *s3 ALIGN; 
  FourSpinor *sp ALIGN;
 
  /*  printf("ID: %d low=%d high=%d: DecompMinus\n", id, low, high); */
  for (ix1=low;ix1<high;ix1++) {
    int thissite = tab->siteTable(ix1);

    sp=&psi[thissite];

    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,0);
    decomp_gamma0_plus(*sp, *s3);

    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,1);
    decomp_gamma1_plus(*sp, *s3);
    
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,2);
    decomp_gamma2_plus(*sp, *s3);
    
    s3 = tab->halfspinorBufferOffset(DECOMP_SCATTER,ix1,3);
    decomp_gamma3_plus(*sp, *s3); 
  
  }
}


/* need gauge fields on opposite cb */
void decomp_hvv_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1 = 0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  GaugeMatrix *um ALIGN;
  HalfSpinor *s3 ALIGN;
  FourSpinor  *sm ALIGN;

  /* printf("ID: %d low=%d high=%d: DecompHvvMinus\n", id, low, high); */

  for (ix1=low;ix1<high;ix1++) 
  {
    int thissite = tab->siteTable(ix1);
    sm=&psi[thissite];

    um=&gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,0);
    decomp_hvv_gamma0_minus(*sm, *um, *s3);	   

    um=&gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,1);
    decomp_hvv_gamma1_minus(*sm, *um, *s3);	   

    um=&gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,2);
    decomp_hvv_gamma2_minus(*sm, *um, *s3);	   
    
    um=&gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(DECOMP_HVV_SCATTER,ix1,3);
    decomp_hvv_gamma3_minus(*sm, *um, *s3);	   
  }
}


void mvv_recons_minus(size_t lo,size_t hi, int id, const void *ptr )
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 

  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();
 
  int cb = a->cb; 
  int ix1;

  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  GaugeMatrix (*gauge_field)[4] = (GaugeMatrix(*)[4])a->u;
  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  GaugeMatrix *up ALIGN;

  HalfSpinor *s3 ALIGN;

  FourSpinor rs ALIGN;
  FourSpinor *rn ALIGN;

  /*  printf("ID: %d low=%d high=%d: MvvReconsMinus\n", id, low, high); */

  for (ix1=low;ix1<high;ix1++) {
    int thissite = tab->siteTable(ix1);
    rn=&psi[thissite];

    up = &gauge_field[thissite][0];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,0);
    mvv_recons_gamma0_minus(*s3, *up, rs);


    up = &gauge_field[thissite][1];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,1);
    mvv_recons_gamma1_minus_add(*s3, *up, rs);


    up = &gauge_field[thissite][2];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,2);
     mvv_recons_gamma2_minus_add(*s3, *up, rs);

    up = &gauge_field[thissite][3];
    s3 = tab->halfspinorBufferOffset(RECONS_MVV_GATHER,ix1,3);
    mvv_recons_gamma3_minus_add_store(*s3, *up, rs,*rn);

  }

}
/******************end of mvv_recons*************************/


void recons_minus(size_t lo,size_t hi, int id, const void *ptr )	
{
  ThreadWorkerArgs *a =(ThreadWorkerArgs *)ptr; 
  ShiftTable<HalfSpinor>* tab = (ShiftTable<HalfSpinor>*)a->s;
  int subgrid_vol_cb = tab->subgridVolCB();

  int cb = a->cb; 
  int ix1=0;
  const int low = cb*subgrid_vol_cb + lo; 
  const int high = cb*subgrid_vol_cb + hi;

  FourSpinor *psi = (FourSpinor*)a->spinor;
  HalfSpinor *chi = (HalfSpinor*)a->half_spinor;

  HalfSpinor *hs0 ALIGN;
  HalfSpinor *hs1 ALIGN;
  HalfSpinor *hs2 ALIGN;
  HalfSpinor *hs3 ALIGN;   
  FourSpinor  *s1 ALIGN,  *rn ALIGN;

  for (ix1=low; ix1<high; ix1++)  {
    int thissite = tab->siteTable(ix1);
    rn=&psi[thissite];

    /* first spin component of result */
    hs0 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,0);
    hs1 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,1);
    hs2 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,2);
    hs3 = tab->halfspinorBufferOffset(RECONS_GATHER,ix1,3);
    recons_4dir_minus(*hs0, *hs1, *hs2, *hs3, *rn);
    /* end of loop */
  }
}
/*****************end of isign corresponding to -1 **********************/

}


  void Dslash<double>::operator()(double* res, 
				    double* psi, 
				    double* u, 
				    int isign,
				    int cb)
    {
      HalfSpinor* chi1 = tab->getChi1();
      HalfSpinor* chi2 = tab->getChi2();
      int subgrid_vol_cb = s_tab->subgridVolCB();
      
      
      if(isign==1) {
	
#ifndef SSEDSLASH_4D_NOCOMMS
	tab->startReceives();
#endif
	
#ifndef SSEDSLASH_4D_NOCOMPUTE
	dispatchToThreads(DslashParscalar64Bit::decomp_plus,
			  (void *)psi,
			  (void *)chi1,
			  (void *)u,
			  (void *)s_tab,
			  cb,
			  subgrid_vol_cb);
#endif
	
#ifndef SSEDSLASH_4D_NOCOMMS
	tab->startSendForward(); 
#endif
	
#ifndef SSEDSLASH_4D_NOCOMPUTE
	dispatchToThreads(DslashParscalar64Bit::decomp_hvv_plus,
			  (void*)psi,
			  (void*)chi2,
			  (void*)u,
			  (void*)s_tab,
			  cb,
			  subgrid_vol_cb);
	
#endif	
	
#ifndef SSEDSLASH_4D_NOCOMMS
	tab->finishSendForward();
	tab->finishReceiveFromBack();
	tab->startSendBack();
#endif   // NOCOMMS
	
#ifndef SSEDSLASH_4D_NOCOMPUTE
	dispatchToThreads(DslashParscalar64Bit::mvv_recons_plus,
			  (void*)res,
			  (void *)chi1,
			  (void *)u,
			  (void *)s_tab,
			  1-cb,
			  subgrid_vol_cb);
#endif
	
#ifndef SSEDSLASH_4D_NOCOMMS
	tab->finishSendBack();
	tab->finishReceiveFromForward();    
#endif  // NOCOMMS
	
	
#ifndef SSEDSLASH_4D_NOCOMPUTE
	dispatchToThreads(DslashParscalar64Bit::recons_plus,
			  (void*)res, 
			  (void*)chi2,
			  (void*)u,	
			  (void*)s_tab,
			  1-cb,
			  subgrid_vol_cb);
#endif
	
      }		
      
      if(isign==-1) 
	{
	  
	  
#ifndef SSEDSLASH_4D_NOCOMMS
	  tab->startReceives();
#endif // NOCOMMS
	  
	  
#ifndef SSEDSLASH_4D_NOCOMPUTE
	  dispatchToThreads(DslashParscalar64Bit::decomp_minus,
			    (void*)psi,
			    (void *)chi1,
			    (void *)u,
			    (void *)s_tab,
			    cb,
			    subgrid_vol_cb);
#endif
	  
#ifndef SSEDSLASH_4D_NOCOMMS
	  tab->startSendForward();
#endif
	  
#ifndef SSEDSLASH_4D_NOCOMPUTE
	  dispatchToThreads(DslashParscalar64Bit::decomp_hvv_minus,
			    (void*)psi,
			    (void*)chi2,
			    (void*)u,
			    (void*)s_tab,
			    cb,
			    subgrid_vol_cb);
#endif
	  
#ifndef SSEDSLASH_4D_NOCOMMS
	  tab->finishSendForward();
	  tab->finishReceiveFromBack();
	  tab->startSendBack();
#endif
	  
#ifndef SSEDSLASH_4D_NOCOMPUTE
	  dispatchToThreads(DslashParscalar64Bit::mvv_recons_minus,
			    (void*)res,
			    (void *)chi1,
			    (void *)u,
			    (void *)s_tab,
			    1-cb,
			    subgrid_vol_cb);
#endif
	  
#ifndef SSEDSLASH_4D_NOCOMMS
	  tab->finishSendBack();
	  tab->finishReceiveFromForward();
#endif // #ifndef NOCOMMS
	  
#ifndef SSEDSLASH_4D_NOCOMPUTE
	  dispatchToThreads(DslashParscalar64Bit::recons_minus,
			    (void*)res, 
			    (void *)chi2,
			    (void *)u,	
			    (void *)s_tab,
			    1-cb,
			    subgrid_vol_cb);
#endif
	}		
    }
    



  


  /* INITIALIZE ROUTINE */
  /* Constructor */
  Dslash<double>::Dslash(const int latt_size[],      
		   void (*getSiteCoords)(int coord[], int node, int linearsite),
		   int (*getLinearSiteIndex)(const int coord[]),
		   int (*getNodeNumber)(const int coord[])
		   ) 
    {
      /* Get the dimensions of the machine */
      const int *machine_size = QMP_get_logical_dimensions();
      
      /* Check we are in 4D */
      if (QMP_get_logical_number_of_dimensions() != 4) {
	QMP_error("init_sse_su3dslash: number of logical dimensions does not match problem");
	QMP_abort(1);
      }
      
      /* Check problem size in 4D */
      for(int mu=0; mu < 4; mu++)  {
	if ( latt_size[mu] % 2 != 0 ) {
	  fprintf(stderr,"This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the GLOBAL requirement latt_size[%d]=%d\n", 
		  mu, latt_size[mu]);
	  
	  exit(1);
	}
      }
      
      /* Check x-checkerboarding */
      if ( (latt_size[0] / machine_size[0]) % 2 != 0 ) {
	fprintf(stderr,"This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even. In addition LOCAL dimension 0 (x) has to be even ,  Your lattice does not meet the LOCAL requirement\n");
	QMP_abort(1);
      }
      
      /* Check requisite number of even local sizes */
      int num_even=0;
      int sx = latt_size[0]/machine_size[0]; if( sx%2 == 0 ) num_even++;
      int sy = latt_size[1]/machine_size[1]; if( sy%2 == 0 ) num_even++;
      int sz = latt_size[2]/machine_size[2]; if( sz%2 == 0 ) num_even++;
      int st = latt_size[3]/machine_size[3]; if( st%2 == 0 ) num_even++;
      
      if( num_even < 2 ) { 
	fprintf(stderr, "Need at least 2 subdimensions to be even");
	QMP_abort(1);
      }
      
      int subgrid[4] = {sx,sy,sz,st};
      tab = new DslashTables<HalfSpinor,4>(subgrid);
      s_tab = new ShiftTable<HalfSpinor>(subgrid, 
					 tab->getChi1(), 
					 tab->getChi2(), 
					 (HalfSpinor*(*)[4])(tab->getRecvBufptr()), 
					 (HalfSpinor*(*)[4])(tab->getSendBufptr()), 
					 getSiteCoords,
					 getLinearSiteIndex,
					 getNodeNumber);
      
    }
    
    /* Destructor */
    Dslash<double>::~Dslash() 
    {
      delete s_tab;
      delete tab;
    }


} // Namespace
