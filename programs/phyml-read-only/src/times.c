/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/


/* Routines for molecular clock trees and molecular dating */


#include "times.h"

/*********************************************************/
#ifdef PHYTIME

int TIMES_main(int argc, char **argv)
{
  align **data;
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
  char *most_likely_tree;
  int i;
  int user_lk_approx;
  
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
  data             = NULL;
  most_likely_tree = NULL;
  best_lnL         = UNLIKELY;
  most_likely_size = -1.0;
  tree_size        = -1.0;

  io = (option *)Get_Input(argc,argv);
  r_seed = (io->r_seed < 0)?(time(NULL)):(io->r_seed);
  io->r_seed = r_seed;
  /* !!!!!!!!!!!!!!!!!!!!!!!! */
/*   r_seed = 1289246338; */
/*   r_seed = 1289266727; */
/*   r_seed = 1289422815; */
/*   r_seed = 1289443891; */
/*   r_seed = 1290652518; */
/*   r_seed = 1292195490; */
  /* r_seed =  1298284669; */
  /* r_seed  = 1298403366; */
  /* r_seed = 1298509108; */
  /* sys = system("sleep 5s"); */
  /* r_seed = 1299649586; */
  /* r_seed = 1302160422; */
  /* r_seed = 1302576741; */
  /* r_seed = 1302588678; */
  /* r_seed = 1303247709; */
  /* r_seed =  1303970631; */
  /* r_seed = 1304059976; */
  /* r_seed = 1306315195; */
  /* r_seed = 1308263660; */

  srand(r_seed); rand();
  PhyML_Printf("\n. Seed: %d\n",r_seed);
  PhyML_Printf("\n. Pid: %d\n",getpid());
  Make_Model_Complete(io->mod);
  mod = io->mod;
  if(io->in_tree == 2) Test_Multiple_Data_Set_Format(io);
  else io->n_trees = 1;

  io->colalias = 1;  /* Do not compress sites if you're using Evolve function */

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
      data = Get_Seq(io);

      if(data)
	{
	  if(io->n_data_sets > 1) PhyML_Printf("\n. Data set [#%d]\n",num_data_set+1);
	  PhyML_Printf("\n. Compressing sequences...\n");
	  cdata = Compact_Data(data,io);
	  Free_Seq(data,cdata->n_otu);
	  Check_Ambiguities(cdata,io->mod->io->datatype,io->mod->state_len);

	  for(num_tree=(io->n_trees == 1)?(0):(num_data_set);num_tree < io->n_trees;num_tree++)
	    {
	      if(!io->mod->s_opt->random_input_tree) io->mod->s_opt->n_rand_starts = 1;

	      For(num_rand_tree,io->mod->s_opt->n_rand_starts)
		{
		  if((io->mod->s_opt->random_input_tree) && (io->mod->s_opt->topo_search != NNI_MOVE))
		    PhyML_Printf("\n. [Random start %3d/%3d]\n",num_rand_tree+1,io->mod->s_opt->n_rand_starts);


		  Init_Model(cdata,mod,io);

		  if(io->mod->use_m4mod) M4_Init_Model(mod->m4mod,cdata,mod);

		  /* A BioNJ tree is built here */
		  if(!io->in_tree) tree = Dist_And_BioNJ(cdata,mod,io);
		  /* A user-given tree is used here instead of BioNJ */
		  else             tree = Read_User_Tree(cdata,mod,io);

		  if(!tree) continue;

		  if(!tree->n_root) 
		    {
		      PhyML_Printf("\n. Sorry, PhyTime requires a rooted tree as input.");
		      Exit("\n");      
		    }

		  time(&t_beg);
		  time(&(tree->t_beg));

		  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
		  RATES_Init_Rate_Struct(tree->rates,io->rates,tree->n_otu);

		  Update_Ancestors(tree->n_root,tree->n_root->v[0],tree);
		  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);		  


		  RATES_Fill_Lca_Table(tree);

		  tree->mod         = mod;
		  tree->io          = io;
		  tree->data        = cdata;
		  tree->both_sides  = 1;
		  tree->n_pattern   = tree->data->crunch_len/tree->mod->state_len;

		  Prepare_Tree_For_Lk(tree);

		  /* Read node age priors */
		  Read_Clade_Priors(io->clade_list_file,tree);

		  /* Set upper and lower bounds for all node ages */
		  TIMES_Set_All_Node_Priors(tree);


		  /* Count the number of time slices */
		  TIMES_Get_Number_Of_Time_Slices(tree);
		  
		  /* Work with log of branch lengths? */
		  if(tree->mod->log_l == YES) Log_Br_Len(tree);
		  
		  /* Force the exact likelihood score */
		  user_lk_approx = tree->io->lk_approx;
		  tree->io->lk_approx = EXACT;

		  /* MLE for branch lengths */
		  PhyML_Printf("\n");
		  Round_Optimize(tree,tree->data,ROUND_MAX);

		  /* Set vector of mean branch lengths for the Normal approximation
		     of the likelihood */
		  RATES_Set_Mean_L(tree);

		  /* Estimate the matrix of covariance for the Normal approximation of
		     the likelihood */
		  PhyML_Printf("\n");
		  PhyML_Printf("\n. Computing Hessian...");
		  tree->rates->bl_from_rt = 0;
		  Free(tree->rates->cov_l);
		  tree->rates->cov_l = Hessian_Seo(tree);
		  /* tree->rates->cov_l = Hessian_Log(tree); */
		  For(i,(2*tree->n_otu-3)*(2*tree->n_otu-3)) tree->rates->cov_l[i] *= -1.0;
		  if(!Iter_Matinv(tree->rates->cov_l,2*tree->n_otu-3,2*tree->n_otu-3,YES)) Exit("\n");
		  tree->rates->covdet = Matrix_Det(tree->rates->cov_l,2*tree->n_otu-3,YES);
		  For(i,(2*tree->n_otu-3)*(2*tree->n_otu-3)) tree->rates->invcov[i] = tree->rates->cov_l[i];
		  if(!Iter_Matinv(tree->rates->invcov,2*tree->n_otu-3,2*tree->n_otu-3,YES)) Exit("\n");
		  tree->rates->grad_l = Gradient(tree);

		  /* Pre-calculation of conditional variances to speed up calculations */
		  RATES_Bl_To_Ml(tree);
		  RATES_Get_Conditional_Variances(tree);
		  RATES_Get_All_Reg_Coeff(tree);
		  RATES_Get_Trip_Conditional_Variances(tree);
		  RATES_Get_All_Trip_Reg_Coeff(tree);

		  Lk(tree);
		  PhyML_Printf("\n");
		  PhyML_Printf("\n. p(data|model) [exact ] ~ %.2f",tree->c_lnL);

		  tree->io->lk_approx = NORMAL;
		  For(i,2*tree->n_otu-3) tree->rates->u_cur_l[i] = tree->rates->mean_l[i] ;
		  tree->c_lnL = Lk(tree);
		  PhyML_Printf("\n. p(data|model) [approx] ~ %.2f",tree->c_lnL);

		  tree->io->lk_approx = user_lk_approx;

		  tree->rates->model = io->rates->model;		  
		  PhyML_Printf("\n. Selected model '%s'",(io->rates->model == THORNE)?("thorne"):("guindon"));
		  if(tree->rates->model == GUINDON) tree->mod->gamma_mgf_bl = YES;
		  
		  tree->rates->bl_from_rt = YES;

		  time(&t_beg);
		  tree->mcmc = MCMC_Make_MCMC_Struct();
		  MCMC_Copy_MCMC_Struct(tree->io->mcmc,tree->mcmc,"phytime");
		  MCMC_Complete_MCMC(tree->mcmc,tree);
		  tree->mcmc->randomize = YES;
		  tree->mcmc->is_burnin = NO;
		  MCMC(tree);
		  MCMC_Close_MCMC(tree->mcmc);
		  MCMC_Free_MCMC(tree->mcmc);
		  PhyML_Printf("\n");


/* 		  tree->mcmc = MCMC_Make_MCMC_Struct(); */
/* 		  MCMC_Copy_MCMC_Struct(tree->io->mcmc,tree->mcmc,"burnin"); */
/* 		  MCMC_Complete_MCMC(tree->mcmc,tree); */
/* 		  tree->mcmc->adjust_tuning = YES; */
/* 		  tree->mcmc->is_burnin     = YES; */
/* 		  tree->mcmc->chain_len = tree->io->mcmc->chain_len_burnin; */
/* 		  MCMC(tree); */
/* 		  MCMC_Close_MCMC(tree->mcmc); */

		  
/* 		  new_mcmc = MCMC_Make_MCMC_Struct(tree); */
/* 		  MCMC_Complete_MCMC(new_mcmc,tree); */
/* 		  MCMC_Copy_MCMC_Struct(tree->mcmc,new_mcmc,"phytime"); */
/* 		  MCMC_Free_MCMC(tree->mcmc); */

/* 		  tree->mcmc                  = new_mcmc; */
/* 		  tree->mcmc->chain_len       = tree->io->mcmc->chain_len; */
/* 		  tree->mcmc->randomize       = NO; */
/* 		  tree->mcmc->adjust_tuning   = NO; */
/* 		  tree->mcmc->is_burnin       = NO; */

/* 		  time(&t_beg); */
/* 		  MCMC(tree); */
/* 		  MCMC_Close_MCMC(tree->mcmc); */
/* 		  MCMC_Free_MCMC(tree->mcmc); */
/* 		  PhyML_Printf("\n"); */


		  Free_Tree_Pars(tree);
		  Free_Tree_Lk(tree);
		  Free_Tree(tree);
		}

	      break;
	    }
	  Free_Cseq(cdata);
	}
    }

  Free_Model(mod);

  if(io->fp_in_align)   fclose(io->fp_in_align);
  if(io->fp_in_tree)    fclose(io->fp_in_tree);
  if(io->fp_out_lk)     fclose(io->fp_out_lk);
  if(io->fp_out_tree)   fclose(io->fp_out_tree);
  if(io->fp_out_trees)  fclose(io->fp_out_trees);
  if(io->fp_out_stats)  fclose(io->fp_out_stats);

  Free(most_likely_tree);
  Free_Input(io);

  time(&t_end);
  Print_Time_Info(t_beg,t_end);

#ifdef MPI
  MPI_Finalize();
#endif

  return 0;
}

#endif

/*********************************************************/

void TIMES_Least_Square_Node_Times(t_edge *e_root, t_tree *tree)
{

  /* Solve A.x = b, where x are the t_node time estimated
     under the least square criterion.

     A is a n x n matrix, with n being the number of
     nodes in a rooted tree (i.e. 2*n_otu-1).

   */

  phydbl *A, *b, *x;
  int n;
  int i,j;
  t_node *root;
  
    
  n = 2*tree->n_otu-1;
  
  A = (phydbl *)mCalloc(n*n,sizeof(phydbl));
  b = (phydbl *)mCalloc(n,  sizeof(phydbl));
  x = (phydbl *)mCalloc(n,  sizeof(phydbl));
  
  
  if(!tree->n_root && e_root) Add_Root(e_root,tree);
  else if(!e_root)            Add_Root(tree->t_edges[0],tree);
  
  root = tree->n_root;

  TIMES_Least_Square_Node_Times_Pre(root,root->v[0],A,b,n,tree);
  TIMES_Least_Square_Node_Times_Pre(root,root->v[1],A,b,n,tree);
  
  b[root->num] = tree->e_root->l/2.;
  
  A[root->num * n + root->num]       = 1.0;
  A[root->num * n + root->v[0]->num] = -.5;
  A[root->num * n + root->v[1]->num] = -.5;
    
  if(!Matinv(A, n, n,YES))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");      
    }

  For(i,n) x[i] = .0;
  For(i,n) For(j,n) x[i] += A[i*n+j] * b[j];

  For(i,n-1) { tree->rates->nd_t[tree->noeud[i]->num] = x[i]; }
  tree->rates->nd_t[root->num] = x[n-1];
  tree->n_root->l[0] = tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[0]->num];
  tree->n_root->l[1] = tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num];


  /* Rescale the t_node times such that the time at the root
     is -100. This constraint implies that the clock rate
     is fixed to the actual tree length divided by the tree
     length measured in term of differences of t_node times */

  phydbl scale_f,time_tree_length,tree_length;

  scale_f = -100./tree->rates->nd_t[root->num];
  For(i,2*tree->n_otu-1) tree->rates->nd_t[i] *= scale_f;
  For(i,2*tree->n_otu-1) if(tree->rates->nd_t[i] > .0) tree->rates->nd_t[i] = .0;


  time_tree_length = 0.0;
  For(i,2*tree->n_otu-3)
    if(tree->t_edges[i] != tree->e_root)
      time_tree_length +=
	FABS(tree->rates->nd_t[tree->t_edges[i]->left->num] -
	     tree->rates->nd_t[tree->t_edges[i]->rght->num]);
  time_tree_length += FABS(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[0]->num]);
  time_tree_length += FABS(tree->rates->nd_t[root->num] - tree->rates->nd_t[root->v[1]->num]);
  
  tree_length = 0.0;
  For(i,2*tree->n_otu-3) tree_length += tree->t_edges[i]->l;

  tree->rates->clock_r = tree_length / time_tree_length;

  Free(A);
  Free(b);
  Free(x);

}

/*********************************************************/

void TIMES_Least_Square_Node_Times_Pre(t_node *a, t_node *d, phydbl *A, phydbl *b, int n, t_tree *tree)
{
  if(d->tax)
    {
      A[d->num * n + d->num] = 1.;
      
      /* Set the time stamp at tip nodes to 0.0 */
/*       PhyML_Printf("\n. Tip t_node date set to 0"); */
      b[d->num] = 0.0;
      return;
    }
  else
    {
      int i;
      
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  TIMES_Least_Square_Node_Times_Pre(d,d->v[i],A,b,n,tree);
      
      A[d->num * n + d->num] = 1.;
      b[d->num] = .0;
      For(i,3)
	{
	  A[d->num * n + d->v[i]->num] = -1./3.;
	  if(d->v[i] != a) b[d->num] += d->b[i]->l;
	  else             b[d->num] -= d->b[i]->l;
	}
      b[d->num] /= 3.;
    }
}

/*********************************************************/

/* Adjust t_node times in order to have correct time stamp ranking with
 respect to the tree topology */

void TIMES_Adjust_Node_Times(t_tree *tree)
{
  TIMES_Adjust_Node_Times_Pre(tree->n_root->v[0],tree->n_root->v[1],tree);
  TIMES_Adjust_Node_Times_Pre(tree->n_root->v[1],tree->n_root->v[0],tree);

  if(tree->rates->nd_t[tree->n_root->num] > MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
						tree->rates->nd_t[tree->n_root->v[1]->num]))
    {
      tree->rates->nd_t[tree->n_root->num] = MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
						 tree->rates->nd_t[tree->n_root->v[1]->num]);
    }
}

/*********************************************************/

void TIMES_Adjust_Node_Times_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl min_height;

      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    TIMES_Adjust_Node_Times_Pre(d,d->v[i],tree);
	  }

      min_height = 0.0;
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      if(tree->rates->nd_t[d->v[i]->num] < min_height)
		{
		  min_height = tree->rates->nd_t[d->v[i]->num];
		}
	    }
	}

      if(tree->rates->nd_t[d->num] > min_height) tree->rates->nd_t[d->num] = min_height;

      if(tree->rates->nd_t[d->num] < -100.) tree->rates->nd_t[d->num] = -100.;

    }
}

/*********************************************************/

  /* Multiply each time stamp at each internal 
     t_node by  'tree->time_stamp_mult'.
   */

void TIMES_Mult_Time_Stamps(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) tree->rates->nd_t[tree->noeud[i]->num] *= FABS(tree->mod->s_opt->tree_size_mult);
  tree->rates->nd_t[tree->n_root->num] *= FABS(tree->mod->s_opt->tree_size_mult);
}

/*********************************************************/

void TIMES_Print_Node_Times(t_node *a, t_node *d, t_tree *tree)
{
  t_edge *b;
  int i;
  
  b = NULL;
  For(i,3) if((d->v[i]) && (d->v[i] == a)) {b = d->b[i]; break;}

  PhyML_Printf("\n. (%3d %3d) a->t = %12f d->t = %12f (#=%12f) b->l = %12f [%12f;%12f]",
	       a->num,d->num,
	       tree->rates->nd_t[a->num],
	       tree->rates->nd_t[d->num],
	       tree->rates->nd_t[a->num]-tree->rates->nd_t[d->num],
	       (b)?(b->l):(-1.0),
	       tree->rates->t_prior_min[d->num],
	       tree->rates->t_prior_max[d->num]);
  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  TIMES_Print_Node_Times(d,d->v[i],tree);
    }
}

/*********************************************************/

void TIMES_Set_All_Node_Priors(t_tree *tree)
{
  int i;
  phydbl min_prior;

  /* Set all t_prior_max values */
  TIMES_Set_All_Node_Priors_Bottom_Up(tree->n_root,tree->n_root->v[0],tree);
  TIMES_Set_All_Node_Priors_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);

  tree->rates->t_prior_max[tree->n_root->num] = 
    MIN(tree->rates->t_prior_max[tree->n_root->num],
	MIN(tree->rates->t_prior_max[tree->n_root->v[0]->num],
	    tree->rates->t_prior_max[tree->n_root->v[1]->num]));


  /* Set all t_prior_min values */
  if(!tree->rates->t_has_prior[tree->n_root->num])
    {
      min_prior = 1.E+10;
      For(i,2*tree->n_otu-2)
	{
	  if(tree->rates->t_has_prior[i])
	    {
	      if(tree->rates->t_prior_min[i] < min_prior)
		min_prior = tree->rates->t_prior_min[i];
	    }
	}
      tree->rates->t_prior_min[tree->n_root->num] = 2.0 * min_prior;
      /* tree->rates->t_prior_min[tree->n_root->num] = 10. * min_prior; */
    }

  TIMES_Set_All_Node_Priors_Top_Down(tree->n_root,tree->n_root->v[0],tree);
  TIMES_Set_All_Node_Priors_Top_Down(tree->n_root,tree->n_root->v[1],tree);

  Get_Node_Ranks(tree);
  TIMES_Set_Floor(tree);
}

/*********************************************************/

void TIMES_Set_All_Node_Priors_Bottom_Up(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl t_sup;
  
  if(d->tax) return;
  else 
    {
      t_node *v1, *v2; /* the two sons of d */

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      TIMES_Set_All_Node_Priors_Bottom_Up(d,d->v[i],tree);	      
	    }
	}
      
      v1 = v2 = NULL;
      For(i,3) if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	{
	  if(!v1) v1 = d->v[i]; 
	  else    v2 = d->v[i];
	}
      
      if(tree->rates->t_has_prior[d->num])
	{
	  t_sup = MIN(tree->rates->t_prior_max[d->num],
		      MIN(tree->rates->t_prior_max[v1->num],
			  tree->rates->t_prior_max[v2->num]));

	  tree->rates->t_prior_max[d->num] = t_sup;

	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      PhyML_Printf("\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
	      PhyML_Printf("\n. Inconsistency in the prior settings detected at t_node %d",d->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->rates->t_prior_max[d->num] = 
	    MIN(tree->rates->t_prior_max[v1->num],
		tree->rates->t_prior_max[v2->num]);
	}
    }
}

/*********************************************************/

void TIMES_Set_All_Node_Priors_Top_Down(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      
      if(tree->rates->t_has_prior[d->num])
	{
	  tree->rates->t_prior_min[d->num] = MIN(tree->rates->t_prior_min[d->num],tree->rates->t_prior_min[a->num]);
	  
	  if(tree->rates->t_prior_max[d->num] < tree->rates->t_prior_min[d->num])
	    {
	      PhyML_Printf("\n. prior_min=%f prior_max=%f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
	      PhyML_Printf("\n. Inconsistency in the prior settings detected at t_node %d",d->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Warn_And_Exit("\n");
	    }
	}
      else
	{
	  tree->rates->t_prior_min[d->num] = tree->rates->t_prior_min[a->num];
	}
            
      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      TIMES_Set_All_Node_Priors_Top_Down(d,d->v[i],tree);
	    }
	}
    }
}

/*********************************************************/

void TIMES_Set_Floor(t_tree *tree)
{
  TIMES_Set_Floor_Post(tree->n_root,tree->n_root->v[0],tree);
  TIMES_Set_Floor_Post(tree->n_root,tree->n_root->v[1],tree);
  tree->rates->t_floor[tree->n_root->num] = MIN(tree->rates->t_floor[tree->n_root->v[0]->num],
						tree->rates->t_floor[tree->n_root->v[1]->num]);
}

/*********************************************************/

void TIMES_Set_Floor_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax)
    {
      tree->rates->t_floor[d->num] = tree->rates->nd_t[d->num];
      d->rank_max = d->rank;
      return;
    }
  else
    {
      int i;
      t_node *v1,*v2;

      v1 = v2 = NULL;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Set_Floor_Post(d,d->v[i],tree);
	      if(!v1) v1 = d->v[i];
	      else    v2 = d->v[i];
	    }
	}
      tree->rates->t_floor[d->num] = MIN(tree->rates->t_floor[v1->num],
					 tree->rates->t_floor[v2->num]);

      if(tree->rates->t_floor[v1->num] < tree->rates->t_floor[v2->num])
	{
	  d->rank_max = v1->rank_max;
	}
      else if(tree->rates->t_floor[v2->num] < tree->rates->t_floor[v1->num])
	{
	  d->rank_max = v2->rank_max;
	}
      else
	{
	  d->rank_max = MAX(v1->rank_max,v2->rank_max);
	}
    }
}

/*********************************************************/
/* Does it work for serial samples? */
phydbl TIMES_Log_Conditional_Uniform_Density(t_tree *tree)
{
  phydbl min,max;
  phydbl dens;
  int i;

  min = tree->rates->nd_t[tree->n_root->num];

  dens = 0.0;
  For(i,2*tree->n_otu-1)
    {
      if((tree->noeud[i]->tax == NO) && (tree->noeud[i] != tree->n_root))
	{
	  max = tree->rates->t_floor[i];

	  dens += LOG(Dorder_Unif(tree->rates->nd_t[i],
				  tree->noeud[i]->rank-1,
				  tree->noeud[i]->rank_max-2,
				  min,max));
	}
    }
  return dens;
}

/*********************************************************/

/* Logarithm of Formula (9) in Rannala and Yang (1996) */
/* Not sure that is works for serial samples */
phydbl TIMES_Log_Yule(t_tree *tree)
{
  phydbl sumti,density,lambda;
  int n,i;

  sumti = 0.0;
  for(i=tree->n_otu;i<2*tree->n_otu-1;i++) sumti += tree->rates->nd_t[i];
  sumti -= tree->rates->nd_t[i-1];

  lambda = tree->rates->birth_rate;
  n = tree->n_otu;
  
  density = 
    (n-1.)*LOG(2.) + 
    (n-2.)*LOG(lambda) - 
    lambda*sumti - 
    Factln(n) - 
    LOG(n-1.) - 
    (n-2.)*LOG(1.-EXP(-lambda));
  
  return density;
}

/*********************************************************/

phydbl TIMES_Lk_Times(t_tree *tree)
{
  phydbl condlogdens;

  condlogdens = 0.0;

  /* TIMES_Lk_Times_Trav(tree->n_root,tree->n_root->v[0], */
  /* 		      tree->rates->nd_t[tree->n_root->num], */
  /* 		      tree->rates->t_floor[tree->n_root->v[0]->num],&condlogdens,tree); */

  /* TIMES_Lk_Times_Trav(tree->n_root,tree->n_root->v[1], */
  /* 		      tree->rates->nd_t[tree->n_root->num], */
  /* 		      tree->rates->t_floor[tree->n_root->v[1]->num],&condlogdens,tree); */


  /* Count the number of internal node (minus the root node) in the deepest slice */
  /* n_nodes = 0; */
  /* For(i,2*tree->n_otu-2) */
  /*   { */
  /*     if(tree->noeud[i]->tax == NO && tree->noeud[i] != tree->n_root) */
  /* 	{ */
  /* 	  if(!(tree->rates->nd_t[i] > tree->rates->t_floor[tree->n_root->num])) */
  /* 	    { */
  /* 	      n_nodes++; */
  /* 	    } */
  /* 	} */
  /*   } */
   
  /* condlogdens = -(phydbl)(n_nodes)*LOG(tree->rates->t_floor[tree->n_root->num] - tree->rates->nd_t[tree->n_root->num]); */

  condlogdens = TIMES_Lk_Uniform_Core(tree);

  tree->rates->c_lnL_times = condlogdens;

  return(condlogdens);
}

/*********************************************************/

void TIMES_Lk_Times_Trav(t_node *a, t_node *d, phydbl lim_inf, phydbl lim_sup, phydbl *logdens, t_tree *tree)
{
  int i;
  
  if(!d->tax)
    {
      /* if(tree->rates->nd_t[d->num] > lim_sup) */
      /* 	{ */
      /* 	  lim_inf = lim_sup; */
      /* 	  lim_sup = 0.0; */
      /* 	  For(i,2*tree->n_otu-2) */
      /* 	    if((tree->rates->t_floor[i] < lim_sup) && (tree->rates->t_floor[i] > tree->rates->nd_t[d->num])) */
      /* 	      lim_sup = tree->rates->t_floor[i]; */
      /* 	} */
      
      /* if(tree->rates->nd_t[d->num] < lim_inf || tree->rates->nd_t[d->num] > lim_sup) */
      /* 	{ */
      /* 	  PhyML_Printf("\n. nd_t = %f lim_inf = %f lim_sup = %f",tree->rates->nd_t[d->num],lim_inf,lim_sup); */
      /* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
      /* 	  Exit("\n"); */
      /* 	} */
  
      lim_inf = tree->rates->nd_t[tree->n_root->num];
      lim_sup = tree->rates->t_floor[d->num];
      
      *logdens = *logdens + LOG(lim_sup - lim_inf);   
    }
  
  if(d->tax == YES) return;
  else
    {      
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Lk_Times_Trav(d,d->v[i],lim_inf,lim_sup,logdens,tree);
	    }
	}
    }
}

/*********************************************************/

phydbl TIMES_Log_Number_Of_Ranked_Labelled_Histories(t_node *root, int per_slice, t_tree *tree)
{
  int i;
  phydbl logn;
  t_node *v1,*v2;
  int dir1r,dir2r;
  int n1,n2;
  
  TIMES_Update_Curr_Slice(tree);

  logn = .0;
  v1 = v2 = NULL;
  if(root == tree->n_root)
    {
      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(root,root->v[0],per_slice,&logn,tree);
      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(root,root->v[1],per_slice,&logn,tree);
      v1 = root->v[0];
      v2 = root->v[1];
    }
  else
    {
      For(i,3)
	{
	  if(root->v[i] != root->anc && root->b[i] != tree->e_root)
	    {
	      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(root,root->v[i],per_slice,&logn,tree);
	      if(!v1) v1 = root->v[i];
	      else    v2 = root->v[i];
	    }
	}
    }

  dir1r = dir2r = -1;
  For(i,3) if(v1->v[i] == root || v1->b[i] == tree->e_root) { dir1r = i; break; }
  For(i,3) if(v2->v[i] == root || v2->b[i] == tree->e_root) { dir2r = i; break; }

  if(per_slice == NO)
    {
      n1 = tree->rates->n_tips_below[v1->num];
      n2 = tree->rates->n_tips_below[v2->num];
    }
  else
    {
      if(tree->rates->curr_slice[v1->num] == tree->rates->curr_slice[root->num])
	n1 = tree->rates->n_tips_below[v1->num];
      else
	n1 = 1;
      
      if(tree->rates->curr_slice[v2->num] == tree->rates->curr_slice[root->num])
	n2 = tree->rates->n_tips_below[v2->num];
      else
	n2 = 1;
    }

  tree->rates->n_tips_below[root->num] = n1+n2;

  logn += Factln(n1+n2-2) - Factln(n1-1) - Factln(n2-1);  

  return(logn);
}

/*********************************************************/

void TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(t_node *a, t_node *d, int per_slice, phydbl *logn, t_tree *tree)
{
  if(d->tax == YES) 
    {
      tree->rates->n_tips_below[d->num] = 1;
      return;
    }
  else
    {
      int i,n1,n2;
      t_node *v1, *v2;
      int dir1d,dir2d;

      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Log_Number_Of_Ranked_Labelled_Histories_Post(d,d->v[i],per_slice,logn,tree);
	    }
	}

      v1 = v2 = NULL;
      dir1d = dir2d = -1;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(v1 == NULL) {v1 = d->v[i];}
	      else           {v2 = d->v[i];}
	    }
	}

      For(i,3) if(v1->v[i] == d) {dir1d = i; break;}
      For(i,3) if(v2->v[i] == d) {dir2d = i; break;}

      if(per_slice == NO)
	{
	  n1 = tree->rates->n_tips_below[v1->num];
	  n2 = tree->rates->n_tips_below[v2->num];
	}
      else
	{
	  if(tree->rates->curr_slice[v1->num] == tree->rates->curr_slice[d->num])
	    n1 = tree->rates->n_tips_below[v1->num];
	  else
	    n1 = 1;

	  if(tree->rates->curr_slice[v2->num] == tree->rates->curr_slice[d->num])
	    n2 = tree->rates->n_tips_below[v2->num];
	  else
	    n2 = 1;
	}

      tree->rates->n_tips_below[d->num] = n1+n2;

      (*logn) += Factln(n1+n2-2) - Factln(n1-1) - Factln(n2-1);
    }
}

/*********************************************************/

void TIMES_Update_Curr_Slice(t_tree *tree)
{
  int i,j;

  For(i,2*tree->n_otu-1)
    {
      For(j,tree->rates->n_time_slices)
	{
	  if(!(tree->rates->nd_t[i] > tree->rates->time_slice_lims[j])) break;
	}
      tree->rates->curr_slice[i] = j;

      /* PhyML_Printf("\n. Node %3d [%12f] is in slice %3d.",i,tree->rates->nd_t[i],j); */
    }
}

/*********************************************************/

phydbl TIMES_Lk_Uniform_Core(t_tree *tree)
{  
  phydbl logn;

  logn = TIMES_Log_Number_Of_Ranked_Labelled_Histories(tree->n_root,YES,tree);


  tree->rates->c_lnL_times = 0.0;
  TIMES_Lk_Uniform_Post(tree->n_root,tree->n_root->v[0],tree);
  TIMES_Lk_Uniform_Post(tree->n_root,tree->n_root->v[1],tree);

  /* printf("\n. ^ %f %f %f", */
  /* 	 (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.), */
  /* 	 LOG(tree->rates->time_slice_lims[tree->rates->curr_slice[tree->n_root->num]] - */
  /* 	     tree->rates->nd_t[tree->n_root->num]), */
  /* 	 (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.) * */
  /* 	 LOG(tree->rates->time_slice_lims[tree->rates->curr_slice[tree->n_root->num]] - */
  /* 	     tree->rates->nd_t[tree->n_root->num])); */
	
  
  tree->rates->c_lnL_times += 
    Factln(tree->rates->n_tips_below[tree->n_root->num]-2.) -
    (phydbl)(tree->rates->n_tips_below[tree->n_root->num]-2.) *
    LOG(tree->rates->time_slice_lims[tree->rates->curr_slice[tree->n_root->num]] -
	tree->rates->nd_t[tree->n_root->num]);
  
  tree->rates->c_lnL_times -= logn;
  
  return(tree->rates->c_lnL_times);
}

/*********************************************************/

void TIMES_Get_Number_Of_Time_Slices(t_tree *tree)
{
  int i;

  tree->rates->n_time_slices=0;
  TIMES_Get_Number_Of_Time_Slices_Post(tree->n_root,tree->n_root->v[0],tree);
  TIMES_Get_Number_Of_Time_Slices_Post(tree->n_root,tree->n_root->v[1],tree);
  Qksort(tree->rates->time_slice_lims,NULL,0,tree->rates->n_time_slices-1);

  if(tree->rates->n_time_slices > 1)
    {
      PhyML_Printf("\n");
      PhyML_Printf("\n. Sequence were collected at %d different time points.",tree->rates->n_time_slices);
      For(i,tree->rates->n_time_slices) printf("\n+ [%3d] time point @ %12f ",i+1,tree->rates->time_slice_lims[i]);
    }
}

/*********************************************************/

void TIMES_Get_Number_Of_Time_Slices_Post(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax == YES)
    {
      For(i,tree->rates->n_time_slices) 
	if(Are_Equal(tree->rates->t_floor[d->num],tree->rates->time_slice_lims[i],1.E-6)) break;

      if(i == tree->rates->n_time_slices) 
	{
	  tree->rates->time_slice_lims[i] = tree->rates->t_floor[d->num];
	  tree->rates->n_time_slices++;
	}
      return;
    }
  else
    {
      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  TIMES_Get_Number_Of_Time_Slices_Post(d,d->v[i],tree);
    }
}

/*********************************************************/

void TIMES_Get_N_Slice_Spans(t_tree *tree)
{
  int i,j;

  For(i,2*tree->n_otu-2)
    {
      if(tree->noeud[i]->tax == NO)
	{
	  For(j,tree->rates->n_time_slices)
	    {
	      if(Are_Equal(tree->rates->t_floor[i],tree->rates->time_slice_lims[j],1.E-6))
		{
		  tree->rates->n_time_slice_spans[i] = j+1;
		  /* PhyML_Printf("\n. Node %3d spans %3d slices [%12f].", */
		  /* 	       i+1, */
		  /* 	       tree->rates->n_slice_spans[i], */
		  /* 	       tree->rates->t_floor[i]); */
		  break;
		}
	    }
	}
    }
}

/*********************************************************/

void TIMES_Lk_Uniform_Post(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax == YES) return;
  else
    {
      int i;

      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      TIMES_Lk_Uniform_Post(d,d->v[i],tree);
	    }
	}
      
      if(tree->rates->curr_slice[a->num] != tree->rates->curr_slice[d->num])
	{
	  tree->rates->c_lnL_times += 
	    Factln(tree->rates->n_tips_below[d->num]-1.) - 
	    (phydbl)(tree->rates->n_tips_below[d->num]-1.) *
	    LOG(tree->rates->time_slice_lims[tree->rates->curr_slice[d->num]] -
		tree->rates->nd_t[d->num]);
	}
    }
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
