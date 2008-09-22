#include <shift_table_scalar.h>

#include <iostream>
#include <cstdlib>
#include <cstddef>
using namespace std;

namespace CPlusPlusWilsonDslash { 

  ShiftTable::ShiftTable(const int latt_size[],  
			 void (*getSiteCoords)(int coord[], int node, int linearsite),
			 int (*getLinearSiteIndex)(const int coord[]),
			 int (*nodeNum)(const int coord[])
			 ) : Nd(4)
  {

    for(int mu=0; mu < 4; mu++)  {
      if ( latt_size[mu] % 2 != 0 ) {
	cerr << "This is a Dslash with checkerboarding in 4 dimensions. Each GLOBAL dimension must be even,  Your lattice is not like this: latt_size[" <<
	  mu <<"]="<<latt_size[mu]<< endl; 
	exit(1);
      }

      /* Copy this internally */
      tot_size[mu] = latt_size[mu];
    }

    /* Compute total volume and total checkerboarded volume */
    total_vol = tot_size[0];
    for(int mu=1; mu < Nd; mu++) { 
      total_vol *= tot_size[mu];
    }
    total_vol_cb = total_vol/2;
      
    xshift_table = (int *)malloc(4*total_vol*2*sizeof(int)+Cache::CacheLineSize);
      
    if ( xshift_table == 0x0 ) {
      cerr << "Could not allocate xshift table" << endl;
      exit(1);
    }

    ptrdiff_t pad = 0;
    if ( (ptrdiff_t)xshift_table % Cache::CacheLineSize != 0 ) {
      pad=Cache::CacheLineSize-((ptrdiff_t)xshift_table % Cache::CacheLineSize);
    }
    shift_table = (int *)((char *)xshift_table + pad);


    
    xsite_table = (int *)malloc(total_vol*sizeof(int)+Cache::CacheLineSize);
    
    if ( xsite_table == 0x0 ) {
      cerr << "Could not allocate site table " << endl;
      exit(1);
    }

    pad = 0;
    if ( (ptrdiff_t)xsite_table % Cache::CacheLineSize != 0 ) {
      pad=(ptrdiff_t)Cache::CacheLineSize-((ptrdiff_t)xsite_table % Cache::CacheLineSize);
    }
    site_table = (int *)((char *)xsite_table + pad);    
    
    for(int p=0; p < 2; p++) { 	    
      for(int t=0; t < tot_size[3]; t++) { 
	for(int z=0; z < tot_size[2]; z++) {
	  for(int y=0; y < tot_size[1]; y++) {     
	    for(int x=0; x < tot_size[0]/2; x++) {
	      int coord[4];
	      
	      coord[0] = 2*x+p;
	      coord[1] = y;
	      coord[2] = z; 
	      coord[3] = t;
	      
	      /* Get the site and N-parity of the chosen victim */
	      int qdp_index = getLinearSiteIndex(coord); /* get the lexico index */
	      int my_index = myLinearSiteIndex4D(coord);
	      
	      /* Add lexico site into site_table, for current cb3 and linear */
	      /* Map (cb3, linear) -> lexico */ 
	      site_table[ my_index ] = qdp_index;
	      
	    }
	  }
	}
      }
    }
    
    /* Get the offsets needed for neighbour comm. */
    /* soffsets(position,direction,isign,cb)   */ 
    /*  where  isign    = +1 : plus direction */
    /*                  =  0 : negative direction */
    /*         cb       =  0 : even lattice (includes origin) */
    /*                  = +1 : odd lattice (does not include origin) */
    /* the offsets cotain the current site, i.e the neighbour for site i  */
    /* is  shift_table(i,dir,cb,mu) and NOT  i + soffset(..)    */
    
    /* Loop over directions and sites, building up shift tables */
    for(int cb=0; cb < 2; cb++) {
      for(int site = 0; site < total_vol_cb; ++site) { 
        int fcoord[4], bcoord[4], coord[4];
	int blinear, flinear;
	
	int my_index = cb*total_vol_cb + site;
	
	int qdp_index = site_table[ my_index ];
	
	getSiteCoords(coord, 0, qdp_index); 
	
	for(int dir=0; dir < 4; dir++) {
	  
	  /* Backwards displacement*/
	  offs(bcoord, coord, dir, -1);
	  blinear = getLinearSiteIndex(bcoord);
	  
	  /* Forward displacement */
	  offs(fcoord, coord, dir, +1);
	  flinear = getLinearSiteIndex(fcoord);
	  
	  
	  /* Gather */
	  shift_table[dir+Nd*my_index ] = blinear;
	  shift_table[dir+Nd*(my_index+total_vol)] = flinear;
	}
      }
    }
    
  }



  void ShiftTable::mySiteCoords4D(int gcoords[], int node, int linearsite)
  {
    int mu;
    int subgrid_cb_nrow[4];
    int tmp_coord[4];
    int cb,cbb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = tot_size[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

    /* Base coordinate single processor: 0,0,0,0 always */
    for(mu=0; mu < 4; mu++) { 
      gcoords[mu] = 0;
    }
    
    cb=linearsite/total_vol_cb;

    crtesn4d(linearsite % total_vol_cb, subgrid_cb_nrow, tmp_coord);

    // Add on position within the node
    // NOTE: the cb for the x-coord is not yet determined
    gcoords[0] += 2*tmp_coord[0];
    for(mu=1; mu < 4; ++mu) {
      gcoords[mu] += tmp_coord[mu];
    }

    cbb = cb;
    for(mu=1; mu < 4; ++mu) {
      cbb += gcoords[mu];
    }
    gcoords[0] += (cbb & 1);
  }

  int ShiftTable::myLinearSiteIndex4D(const int gcoords[]) 
  {
    int mu;
    int subgrid_cb_nrow[4];
    int subgrid_cb_coord[4];
    int cb;

    for(mu=0; mu < 4; mu++) { 
      subgrid_cb_nrow[mu] = tot_size[mu];
    }
    subgrid_cb_nrow[0] /=2;  /* Checkerboarding */

    cb=0;
    for(mu=0; mu < Nd; ++mu) { 
      cb += gcoords[mu];
    }
    cb &=1;
    
    subgrid_cb_coord[0] = (gcoords[0]/2)% subgrid_cb_nrow[0];
    for(mu=1; mu < 4; mu++) { 
      subgrid_cb_coord[mu] = gcoords[mu] % subgrid_cb_nrow[mu];
    }

    return localSite4d(subgrid_cb_coord, subgrid_cb_nrow) + cb*total_vol_cb;
  }

}
