/* 

Calculates the mass accretion history along the main branch (i.e. for
direct progenitors) for consistent-trees merger trees.

TODO 

    + Implement total 

*/

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "read_tree.h"

int main(void) {
  int i, ii, jj, ic, c, nh, np, nl, ns, check, ihnp;
  int j = 0;
  long iid, iid_next;
  struct halo_list *hl, *ihl;
  struct halo ih, ihn;
  struct halo *ihp;

  FILE *fptr;
  
  read_tree("tree_0_0_0.dat");
  // printf("%"PRId64" halos found in tree_0_0_0.dat!\n", all_halos.num_halos);
  
  nh = all_halos.num_halos;
  nl = halo_tree.num_lists;
  ns = halo_tree.num_scales;

  printf("Found %d haloes, %d lists, %d scales \n", nh, nl, ns);
  
  long root_ids[nh];  // tree root IDs
  float aexp[nl]; // scale factors 
  float mvir[nl]; // virial masses
  
  
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

  // Open the file pointer and prepare the file for writing. Overwrite
  // it if it exists.
  fptr = fopen("./out.txt", "w");
  //The header will contain the number of roots and scales. Don't use
  // /n so the first # will write on the same line
  fprintf(fptr, "# %d %d ", j, nl);
  
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
    aexp[c] = ih.scale;
    mvir[c] = ih.mvir;
    c += 1;
    
    while (ihnp > 0) {
      ihn = *ih.prog;
      ih = ihn;
      iid = ih.id;
      ihnp = ih.num_prog;
      printf("-------- halo %d has %d progenitors \n", iid, ihnp);

      // Store properties
      aexp[c] = ih.scale;
      mvir[c] = ih.mvir;
      c += 1;
    }

    fprintf(fptr, "#\n"); // print # before each entry
    for (ic=0; ic<c;ic++) {
      fprintf(fptr, "%f %f\n", aexp[ic], mvir[ic]);
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

  fclose(fptr);
  
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
