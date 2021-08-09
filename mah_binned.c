/* 

Calculates the mass accretion history along the main branch (i.e. for
direct progenitors) for consistent-trees merger trees.

TODO 

    + Implement total 

*/

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <math.h>
#include "read_tree.h"

int get_bin(float m, float *mb, int nm) {
  // Get the index of the bin [0, nm) (bin edges defined by mb) in
  // which m belongs. Return -1 if unbinnable.
  //
  // ii = get_bin(m, mb, nm);
  
  int i, j;

  j = -1;
  for (i=0; i<nm; i++) {
    if ((m > mb[i]) && (m <= mb[i+1])) {
      j = i;
      break;
    }
  }
  
  return j;
}


int main(void) {
  int i, ii, jj, ic, c, nh, np, nl, ns, check, ihnp;
  int j = 0, cc=-1;
  long iid, iid_next;
  float lmmin, lmmax, mmin, mmax, dlogm, dm, m0, m02;
  int nm, im;
  struct halo_list *hl, *ihl;
  struct halo ih, ihn;
  struct halo *ihp;

  FILE *fptr1, *fptr2, *fptr3, *fptr4;

  read_tree("tree_0_0_0.dat");
  // printf("%"PRId64" halos found in tree_0_0_0.dat!\n", all_halos.num_halos);
  
  nh = all_halos.num_halos;
  nl = halo_tree.num_lists;
  ns = halo_tree.num_scales;

  // Set up log-spaced mass bins
  mmin = 1e8;   // Msol/h
  mmax = 1e12;  // Msol/h
  lmmin = log10(mmin);
  lmmax = log10(mmax);
  nm = 8;
  dlogm = (lmmax - lmmin) / (float)nm;
  float mb[nm+1];
  float mba[nm];

  printf("dlogm %f", dlogm); //, dm);
  
  // Bin edges
  for (i=0; i<nm+1; i++) {
    mb[i] = pow(10, lmmin + (i * dlogm));
    printf("bin %d %f\n", i, mb[i]);
  }
  mb[nm+1] = mmax;

  // Bin centres
  for (i=0; i<nm; i++) {
    mba[i] = 0.5 * (mb[i] + mb[i+1]);
  }

  // n(z), M(z), M(z)^2
  int nz[nm][nl];
  float mz[nm][nl], mz2[nm][nl];

  memset(nz, 0, nm * nl * sizeof(int));
  memset(mz, 0, nm * nl * sizeof(float));
  memset(mz2, 0, nm * nl * sizeof(float));
    
  
  
  printf("Found %d haloes, %d lists, %d scales \n", nh, nl, ns);
  
  long root_ids[nh];  // tree root IDs
  float aexp[nl], z[nl]; // scale factors 
  // float mvir[nl]; // virial masses
  
  
  // Now we have read in the haloes, we want to do something with them

  // Look at all the haloes with more than two progenitors
  
  /* for (i = 0; i < nh; ++i) { */
  /*   np = all_halos.halos[i].num_prog; */
  /*   if (np > 2) { */
  /* 	printf("number of progenitors %d \n", np); */
  /*     } */
  /* } */

  
  for (i=0; i<nh; i++) {
    ih = all_halos.halos[i];
    iid = ih.tree_root_id;
      
    // printf("%ld \n", ih.tree_root_id);
    // Check if we have this id
    check = 1;
    for (jj=0; jj<j; jj++) {
      if (root_ids[jj] == iid)
	check = 0;
    }
    // If we don't already have the ID, store it
    if (check) {
      root_ids[j] = iid;
      j += 1;
    }
  }

  printf("%d unique roots found \n", j);

  
  // Walk back up the tree
  for (i=0; i<j; i++) {
    iid = root_ids[i];
    printf("---- working on root halo %d \n", iid);
    // The root halo will be in the latest list ihl =
    // &(halo_tree.halo_lists[0]);

    // Why bother using the individual halo lists? Instead can we use
    // the all_halos halo_list?
    ihl = &(all_halos);
    
    // printf("---- selected halo list at scale %f \n", ihl->scale);

    // For some reason this throws a segfault...
    // ih = *lookup_halo_in_list(ihl, iid);  // *: struct *halo -> struct halo
    // Search using a loop instead
    for (ii=0; ii<nh; j++) {
      ih = all_halos.halos[i];
      if (ih.id == iid)
	break;
    }
    
    ihnp = ih.num_prog;
    iid_next = ih.depth_first_id;
    printf("-------- halo %d has %d progenitors \n", iid, ihnp);
    printf("-------- depth-first id %d \n", iid_next);

    // Store properties
    c = 0;
    if (c > cc) {
      aexp[cc+1] = ih.scale;
      z[cc+1] = 1./aexp[cc+1] - 1.;  // TODO check these are unique
      cc++;
    }
    /* mvir[c] = ih.mvir; */
    im = get_bin(ih.mvir, mb, nm);  // mass bin
    m0 = ih.mvir;
    m02 = ih.mvir * ih.mvir;
    nz[im][c] += 1;
    mz[im][c] += 1.;
    mz2[im][c] += 1;
    c += 1;
    
    while (ihnp > 0) {
      ihn = *ih.prog;
      ih = ihn;
      iid = ih.id;
      ihnp = ih.num_prog;
      printf("-------- halo %d has %d progenitors \n", iid, ihnp);

      if (c > cc) {
	aexp[cc+1] = ih.scale;
	z[cc+1] = 1./aexp[cc+1] - 1.;  // TODO check these are unique
	cc++;
      }

      nz[im][c] += 1;
      mz[im][c] += ih.mvir / m0;
      mz2[im][c] += ih.mvir * ih.mvir / m02;

      // Store properties
      /* aexp[c] = ih.scale; */
      /* mvir[c] = ih.mvir; */
      c += 1;
    }

    /* // Now move up the tree */
    /* while (ihnp > 0) { */
    /*   printf("-------- working on halo %d \n", iid_next); */
    /*   ih = *lookup_halo_in_list(ihl, iid_next); */
    /*   ihnp = ih.num_prog; */
    /*   iid_next = ih.depth_first_id; */
    /*   printf("-------- halo %d has %d progenitors \n", iid, ihnp); */
    /*   printf("-------- depth-first id %d \n", ih.depth_first_id); */
    /* } */
    
  }

  // Open the file pointer and prepare the file for writing. Overwrite
  // it if it exists.
  fptr1 = fopen("./out_nz.txt", "w");
  fptr2 = fopen("./out_mz.txt", "w");
  fptr3 = fopen("./out_sz.txt", "w");
  fptr4 = fopen("./out_mba.txt", "w");
  //The header will contain the number of roots and scales. Don't use
  // /n so the first # will write on the same line

  for (im=0; im<nm; im++) {
    fprintf(fptr4, "%f ", mba[im]);
    for (ic=0; ic<nl; ic++) {
      mz[im][ic] /= (float)nz[im][ic];   // <M(z)>
      mz2[im][ic] /= (float)nz[im][ic];  // <M(z)^2>
      fprintf(fptr1, "%d ", nz[im][ic]);
      fprintf(fptr2, "%f ", mz[im][ic]);
      fprintf(fptr3, "%f ", sqrt(mz2[im][ic] -
				 (mz[im][ic] * mz[im][ic])));  // sigma 
    }
    fprintf(fptr1, "\n");
    fprintf(fptr2, "\n");
    fprintf(fptr3, "\n");
  }
  fclose(fptr1);
  fclose(fptr2);
  fclose(fptr3);

  fptr1 = fopen("./out_z.txt", "w");
  for (ic=0; ic<nl; ic++) {
    fprintf(fptr1, "%f ", z[ic]);
  }
  fclose(fptr1);
  // Now store the unique haloes in their own array. Do we need to do
  // this? We already have the length of the nonzero elements j.
  /* long uroot_ids[j]; */
  /* for (i=0, i<j, i++){ */
  /*   uroot_ids[i] = root_ids[j];  */
  /* } */

  
  // Look at the halo tree
  /* nt = halo_tree.num_lists; */
  /* printf("Found %d lists \n", nt); */
  /* for (i=0; i<nt; i++) { */
  /*   hl = &(halo_tree.halo_lists[i]); */
  /*   for (j=0; j<hl->num_halos; ++j) { */
  /*     printf("%d %d\n", i, hl->halos[j].descid);//, hl->halos[j].index); */
  /*   } */
  /* } */
  // print("Found %d haloes in those lists \n", halo_tree.num_halos);
  return 0;
}
