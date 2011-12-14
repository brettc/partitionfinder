/*

PhyML:  a program that  computes maximum likelihood phyLOGenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "spr.h"
#include "utilities.h"
#include "lk.h"
#include "optimiz.h"
#include "bionj.h"
#include "models.h"
#include "free.h"
#include "help.h"
#include "simu.h"
#include "eigen.h"
#include "pars.h"
#include "alrt.h"

#ifdef MPI
#include "mpi_boot.h"
#endif

#ifdef PHYML

int main(int argc, char **argv)
{
  
  calign *cdata;
  option *io;
  t_tree *tree;
  int n_otu, num_data_set;
  int num_tree,tree_line_number,num_rand_tree;
  matrix *mat;
  model *mod;
  time_t t_beg,t_end;
  phydbl best_lnL,most_likely_size,tree_size;
  int r_seed;
  char *most_likely_tree=NULL;

  
#ifdef MPI
  int rc;
  rc = MPI_Init(&argc,&argv);
  if (rc != MPI_SUCCESS) {
    PhyML_Printf("\n. Error starting MPI program. Terminating.\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  MPI_Comm_size(MPI_COMM_WORLD,&Global_numTask);
  MPI_Comm_rank(MPI_COMM_WORLD,&Global_myRank);
#endif

#ifdef QUIET
  setvbuf(stdout,NULL,_IOFBF,2048);
#endif


  tree             = NULL;
  mod              = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;

  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  srand(r_seed);
  io->r_seed = r_seed;


  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  mat = NULL;
  tree_line_number = 0;


  if((io->n_data_sets > 1) && (io->n_trees > 1))
    {
      io->n_data_sets = MIN(io->n_trees,io->n_data_sets);
      io->n_trees     = MIN(io->n_trees,io->n_data_sets);
    }

  For(num_data_set,io->n_data_sets)
    {
      n_otu = 0;
      best_lnL = UNLIKELY;
      Get_Seq(io);
      Make_Model_Complete(io->mod);
      Set_Model_Name(io->mod);
      Print_Settings(io);
      mod = io->mod;
        
      
      if(io->data)
	{
	  if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
	  cdata = Compact_Data(io->data,io);

	  Free_Seq(io->data,cdata->n_otu);
	  
	  if(cdata) Check_Ambiguities(cdata,io->datatype,io->mod->state_len);
	  else
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
		    if(!io->quiet) PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);


		  Init_Model(cdata,mod,io);



		  switch(io->in_tree)
		    {
		    case 0 : case 1 : { tree = Dist_And_BioNJ(cdata,mod,io);    break; }
		    case 2 :          { tree = Read_User_Tree(cdata,mod,io); break; }
		    }


		  if(!tree) continue;

		  time(&t_beg);
		  time(&(tree->t_beg));
      
		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = cdata;
		  tree->both_sides  = YES;
		  tree->n_pattern   = tree->data->crunch_len;

		  if(mod->s_opt->random_input_tree) Random_Tree(tree);

		  if((!num_data_set) && (!num_tree) && (!num_rand_tree)) Check_Memory_Amount(tree);		  

		  Prepare_Tree_For_Lk(tree);

		  if(io->in_tree == 1) Spr_Pars(tree);
		 
		  if(io->do_alias_subpatt)
		    {
		      tree->update_alias_subpatt = YES;
		      Lk(tree);
		      tree->update_alias_subpatt = NO;
		    }
		  
		  if(tree->mod->s_opt->opt_topo)
		    {
		      if(tree->mod->s_opt->topo_search      == NNI_MOVE) Simu_Loop(tree);
		      else if(tree->mod->s_opt->topo_search == SPR_MOVE) Speed_Spr_Loop(tree);
		      else                                               Best_Of_NNI_And_SPR(tree);

		      if(tree->n_root) Add_Root(tree->t_edges[0],tree);
		    }
		  else
		    {
		      if(tree->mod->s_opt->opt_subst_param || 
			 tree->mod->s_opt->opt_bl)                       Round_Optimize(tree,tree->data,ROUND_MAX);
		      else                                               Lk(tree);
		    }

		  tree->both_sides = 1;
		  Lk(tree);
		  Pars(tree);
		  Get_Tree_Size(tree);
		  PhyML_Printf("\n. Log likelihood of the current tree: %f.\n",tree->c_lnL);

		  Br_Len_Involving_Invar(tree);
		  Rescale_Br_Len_Multiplier_Tree(tree);

		  if(!tree->n_root) Get_Best_Root_Position(tree);

		  /* Print the tree estimated using the current random (or BioNJ) starting tree */
		  if(io->mod->s_opt->n_rand_starts > 1)
		    {
		      Print_Tree(io->fp_out_trees,tree);
		      fflush(NULL);
		    }

		  /* Record the most likely tree in a string of characters */
		  if(tree->c_lnL > best_lnL)
		    {
		      best_lnL = tree->c_lnL;
		      if(most_likely_tree) Free(most_likely_tree);
		      most_likely_tree = Write_Tree(tree);
		      most_likely_size = Get_Tree_Size(tree);
		    }

/* 		  JF(tree); */

		  time(&t_end);
		  Print_Fp_Out(io->fp_out_stats,t_beg,t_end,tree,
			       io,num_data_set+1,
			       (tree->mod->s_opt->n_rand_starts > 1)?
			       (num_rand_tree):(num_tree));
		  
		  if(tree->io->print_site_lnl) Print_Site_Lk(tree,io->fp_out_lk);

		  /* Start from BioNJ tree */
		  if((num_rand_tree == io->mod->s_opt->n_rand_starts-1) && (tree->mod->s_opt->random_input_tree))
		    {
		      /* Do one more iteration in the loop, but don't randomize the tree */
		      num_rand_tree--;
		      tree->mod->s_opt->random_input_tree = 0;
		    }
		  
		  Free_Spr_List(tree);
		  Free_One_Spr(tree->best_spr);
		  if(tree->mat) Free_Mat(tree->mat);
		  Free_Triplet(tree->triplet_struct);
		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}


	      /* Launch bootstrap analysis */
	      if(mod->bootstrap) 
		{
		  if(!io->quiet) PhyML_Printf("\n. Launch bootstrap analysis on the most likely tree...\n");

                  #ifdef MPI
		  MPI_Bcast (most_likely_tree, strlen(most_likely_tree)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
		  if(!io->quiet)  PhyML_Printf("\n. The bootstrap analysis will use %d CPUs.\n",Global_numTask);
		  #endif

		  most_likely_tree = Bootstrap_From_String(most_likely_tree,cdata,mod,io);
		}
	      else if(io->ratio_test) 
		{
		  /* Launch aLRT */
		  if(!io->quiet) PhyML_Printf("\n. Compute aLRT branch supports on the most likely tree...\n");
		  most_likely_tree = aLRT_From_String(most_likely_tree,cdata,mod,io);
		}

	      /* Print the most likely tree in the output file */
	      if(!io->quiet) PhyML_Printf("\n. Printing the most likely tree in file '%s'...\n", Basename(io->out_tree_file));
	      if(io->n_data_sets == 1) rewind(io->fp_out_tree);
	      PhyML_Fprintf(io->fp_out_tree,"%s\n",most_likely_tree);
	      

	      if(io->n_trees > 1 && io->n_data_sets > 1) break;
	    }
	  Free_Cseq(cdata);
	}
      else
	{
	  PhyML_Printf("\n. No data was found.\n");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      Free_Model_Complete(mod);
    }
  
  if(most_likely_tree) Free(most_likely_tree);

  if(mod->s_opt->n_rand_starts > 1) PhyML_Printf("\n. Best log likelihood: %f\n",best_lnL);

  Free_Optimiz(mod->s_opt);
  Free_Model_Basic(mod);
  Free_Custom_Model(mod);
  M4_Free_M4_Model(mod->m4mod);
  Free(mod);

  if(io->fp_in_align)  fclose(io->fp_in_align);
  if(io->fp_in_tree)   fclose(io->fp_in_tree);
  if(io->fp_out_lk)    fclose(io->fp_out_lk);
  if(io->fp_out_tree)  fclose(io->fp_out_tree);
  if(io->fp_out_trees) fclose(io->fp_out_trees);
  if(io->fp_out_stats) fclose(io->fp_out_stats);

  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}


#elif(M4)
#include "m4.h"
int main(int argc, char **argv)
{
  M4_main(argc, argv);
  return 1;
}

#elif(PART)
#include "mg.h"
int main(int argc, char **argv)
{
  PART_main(argc, argv);
  return 1;
}

#elif(PHYTIME)
#include "times.h"
int main(int argc, char **argv)
{
  TIMES_main(argc, argv);
  return 1;
}

#elif(PHYCONT)
#include "continuous.h"
int main(int argc, char **argv)
{
  CONT_main(argc, argv);
  return 1;
}

#elif(RF)
int main(int argc, char **argv)
{
  t_tree *tree1, *tree2;
  FILE *fp_tree1, *fp_tree2;
  int i,j;

  fp_tree1 = (FILE *)fopen(argv[1],"r");
  fp_tree2 = (FILE *)fopen(argv[2],"r");

  tree1 = Read_Tree_File(fp_tree1);
  tree2 = Read_Tree_File(fp_tree2);


  For(i,tree1->n_otu)
    {
      For(j,tree2->n_otu)
	{
	  if(!strcmp(tree1->noeud[i]->name,tree2->noeud[j]->name))
	    {
	      break;
	    }
	}

      if(j == tree2->n_otu)
	{
	  Prune_Subtree(tree1->noeud[i]->v[0],tree1->noeud[i],NULL,NULL,tree1);
	}
    }

  For(i,tree2->n_otu)
    {
      For(j,tree1->n_otu)
	{
	  if(!strcmp(tree2->noeud[i]->name,tree1->noeud[j]->name))
	    {
	      break;
	    }
	}

      if(j == tree1->n_otu)
	{
	  Prune_Subtree(tree2->noeud[i]->v[0],tree2->noeud[i],NULL,NULL,tree2);
	}
    }
  
  PhyML_Printf("%s\n",Write_Tree(tree1));
  PhyML_Printf("%s\n",Write_Tree(tree2));



/*   t_tree *tree1, *tree2; */
/*   FILE *fp_tree1, *fp_tree2; */
/*   int i,j,rf,n_edges,n_common,bip_size; */
/*   phydbl thresh; */
/*   t_edge *b; */


/*   fp_tree1 = (FILE *)fopen(argv[1],"r"); */
/*   fp_tree2 = (FILE *)fopen(argv[2],"r"); */
/*   thresh = (phydbl)atof(argv[3]); */

/*   tree1 = Read_Tree_File(fp_tree1); */
/*   tree2 = Read_Tree_File(fp_tree2); */

/*   Get_Rid_Of_Prefix('_',tree1); */

/* /\*   Find_Common_Tips(tree1,tree2); *\/ */

/*   Alloc_Bip(tree1); */
/*   Alloc_Bip(tree2); */

/*   Get_Bip(tree1->noeud[0],tree1->noeud[0]->v[0],tree1); */
/*   Get_Bip(tree2->noeud[0],tree2->noeud[0]->v[0],tree2); */
  
/* /\*   PhyML_Printf("\n. rf=%f\n",Compare_Bip_On_Existing_Edges(thresh,tree1,tree2)); *\/ */
/*   For(i,2*tree1->n_otu-3) tree1->t_edges[i]->bip_score = 0; */
/*   For(i,2*tree2->n_otu-3) tree2->t_edges[i]->bip_score = 0; */
  
/*   rf = 0; */
/*   n_edges = 0; */

/*   /\* First tree *\/ */
/*   For(i,2*tree1->n_otu-3)  */
/*     { */
/*       /\* Consider the branch only if the corresponding bipartition has size > 1 *\/ */
/*       b = tree1->t_edges[i]; */
/*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */
	  
/*       if(bip_size > 1) */
/* 	{ */
/* 	  /\* with non-zero length *\/ */
/* 	  if(tree1->t_edges[i]->l > thresh)   */
/* 	    { */
/* 	      n_edges++; */
/* 	      /\* This t_edge is not found in tree2 *\/ */
/* 	      if(!tree1->t_edges[i]->bip_score) rf++; ; */
/* 	    } */
/* 	} */
/*     } */


/*   /\* Second tree *\/ */
/*   For(i,2*tree2->n_otu-3)  */
/*     { */
/*       b = tree2->t_edges[i]; */
/*       bip_size = MIN(b->left->bip_size[b->l_r],b->rght->bip_size[b->r_l]); */

/*       if(bip_size > 1) */
/* 	{ */
/* 	  if(tree2->t_edges[i]->l > thresh)   */
/* 	    { */
/* 	      n_edges++; */
/* 	      /\* This t_edge is not found in tree1 *\/ */
/* 	      if(!tree2->t_edges[i]->bip_score) rf++; ; */
/* 	    } */
/* 	} */
/*     } */

/*   if(!n_edges) */
/*     { */
/*       Exit("\n. No comparable internal edges were found.\n"); */
/*     } */
/*   else */
/*     { */
/*       PhyML_Printf("\n. Robinson and Foulds distance: %f.",(double)rf/(n_edges)); */
/* /\*       PhyML_Printf("\n. %d internal edges were processed (%d in the first tree, %d in the second).\n",n_edges,n_edges_t1,n_edges-n_edges_t1); *\/ */
/*       PhyML_Printf("\n"); */
/*     } */

  return 1;
}

#elif(TIPORDER)
#include "tiporder.h"
int main(int argc, char **argv)
{
  TIPO_main(argc, argv);
  return 1;
}

#endif

