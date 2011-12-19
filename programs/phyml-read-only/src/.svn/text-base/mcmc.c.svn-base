/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homoLOGous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/



#include "mcmc.h"

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC(t_tree *tree)
{
  int move;
  phydbl u;
  int first,secod;
  int i;

  RATES_Set_Clock_And_Nu_Max(tree);
  RATES_Set_Birth_Rate_Boundaries(tree);

  if(tree->mcmc->randomize == YES)
    {
      MCMC_Randomize_Birth(tree);
      MCMC_Randomize_Nu(tree);
      MCMC_Randomize_Node_Times(tree);
      MCMC_Sim_Rate(tree->n_root,tree->n_root->v[0],tree);
      MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree);
      MCMC_Randomize_Node_Rates(tree);
      MCMC_Randomize_Clock_Rate(tree); /* Clock Rate must be the last parameter to be randomized */
      MCMC_Randomize_Rate_Across_Sites(tree);
      MCMC_Randomize_Kappa(tree);
      MCMC_Randomize_Covarion_Rates(tree);
      MCMC_Randomize_Covarion_Switch(tree);
    }
  else
    {
      MCMC_Read_Param_Vals(tree);
    }

  tree->mod->update_eigen = YES;

  MCMC_Initialize_Param_Val(tree->mcmc,tree);
  
  Update_Ancestors(tree->n_root,tree->n_root->v[0],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  RATES_Lk_Rates(tree);
  TIMES_Lk_Times(tree);
  tree->both_sides = NO;
  if(tree->mcmc->use_data) Lk(tree);
  else tree->c_lnL = 0.0;
  tree->mod->update_eigen = NO;
  MCMC_Print_Param(tree->mcmc,tree);

       
  //////////////////
  if(tree->io->mutmap == YES)
    {
      int j;
      char *s,*t;
      FILE *fp;
      
      Make_Ancestral_Seq(tree);
      Make_MutMap(tree);
 
      For(j,tree->n_otu)
	{
	  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	  strcpy(s,tree->t_nodes[j]->name);
	  tree->t_nodes[j]->name = s;
	}
      
      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      t = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      
      tree->write_tax_names = YES;
      For(i,tree->n_pattern)
	{
	  For(j,tree->n_otu)
	    {
	      strcpy(t,tree->t_nodes[j]->name);
	      s[0]=tree->data->c_seq[j]->state[i];
	      s[1]='\0';
	      strcat(s,"--");
	      sprintf(s+strlen(s),"%.0f",tree->rates->nd_t[j]);
	      /* strcat(s,tree->t_nodes[j]->name); */
	      strcpy(tree->t_nodes[j]->name,s);
	    }
	  
	  strcpy(s,"rosettatree.");
	  sprintf(s+strlen(s),"%d",i);
	  fp = fopen(s,"w");
	  s = Write_Tree(tree,NO);
	  PhyML_Fprintf(fp,"%s",s);
	  fclose(fp);
	  
	  For(j,tree->n_otu) strcpy(tree->t_nodes[j]->name,t);
	}
      Free(s);
      Free(t);
    }


  first = 0;
  secod = 1;
  do
    {
      /* if(tree->mcmc->ess[tree->mcmc->num_move_tree_height] > 100 &&  */
      /* 	 tree->mcmc->ess[tree->mcmc->num_move_nu] > 100          && */
      /* 	 tree->mcmc->ess[tree->mcmc->num_move_clock_r] > 100     && */
      /* 	 tree->mcmc->run > 1000) */
      /* 	{ */
      /* 	  FILE *fp; */
      /* 	  char *s; */

      /* 	  s = (char *)mCalloc(100,sizeof(char)); */

      /* 	  sprintf(s,"simul_par.%d",getpid()); */
      /* 	  fclose(tree->mcmc->out_fp_stats); */
      /* 	  tree->mcmc->out_fp_stats = fopen(s,"w"); */
      /* 	  tree->mcmc->run = 0; */
      /* 	  tree->mcmc->nd_t_digits = 4; */
      /* 	  MCMC_Print_Param(tree->mcmc,tree); */

      /* 	  RATES_Update_Cur_Bl(tree); */
      /* 	  printf("\n. %s",Write_Tree(tree,NO)); */
      /* 	  Evolve(tree->data,tree->mod,tree); */

      /* 	  sprintf(s,"simul_seq.%d",getpid()); */
      /* 	  fp = fopen(s,"w"); */
      /* 	  Print_CSeq(fp,NO,tree->data); */
      /* 	  fflush(NULL); */
      /* 	  fclose(fp); */
      /* 	  Free(s); */

      /* 	  Exit("\n"); */
      /* 	} */


      if(tree->mcmc->ess[tree->mcmc->num_move_tree_height] > 100 &&
      	 tree->mcmc->ess[tree->mcmc->num_move_nu] > 100          &&
      	 tree->mcmc->ess[tree->mcmc->num_move_clock_r] > 100     &&
      	 tree->mcmc->run > 1000)
      	{
      	  FILE *fp;
      	  char *s,*t;

      	  s = (char *)mCalloc(100,sizeof(char));
	  
      	  t = strrchr(tree->io->in_align_file,'.');
      	  sprintf(s,"res%s",t);
      	  fp = fopen(s,"w");
      	  fclose(tree->mcmc->out_fp_stats);
      	  tree->mcmc->out_fp_stats = fopen(s,"w");
      	  tree->mcmc->run = 0;
      	  MCMC_Print_Param(tree->mcmc,tree);
      	  fclose(fp);
      	  Free(s);
      	  Exit("\n");
      	}

      /* if(tree->mcmc->run > 50000) */
      /* 	{ */
      /* 	  int i; */
      /* 	  For(i,tree->n_otu) */
      /* 	    printf("\n %15f %15f %15f",tree->rates->nd_t[i],EXP(tree->rates->br_r[i]),tree->rates->br_r[i]); */
      /* 	  Exit("\n"); */
      /* 	} */

      u = Uni();

      For(move,tree->mcmc->n_moves) if(tree->mcmc->move_weight[move] > u) break;
      
      if(u < .5) { first = 0; secod = 1; }
      else       { first = 1; secod = 0; }
      

      phydbl cur_lnL;
      cur_lnL = tree->c_lnL;

      /* printf("\n. [%15s] %15f ",tree->mcmc->move_name[move],tree->c_lnL); */

      /* Clock rate */
      if(!strcmp(tree->mcmc->move_name[move],"clock"))
      	{
      	  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
	  
	  MCMC_Single_Param_Generic(&(tree->rates->clock_r),tree->rates->min_clock,tree->rates->max_clock,move,
      	  			    NULL,&(tree->c_lnL),
      	  			    NULL,Wrap_Lk,tree->mcmc->move_type[move],NO,NULL,tree,NULL);
      	}

      /* Nu */
      else if(!strcmp(tree->mcmc->move_name[move],"nu"))
      	{
      	  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
	  MCMC_Nu(tree);
      	}

      /* Tree height */
      else if(!strcmp(tree->mcmc->move_name[move],"tree_height"))
      	{
      	  MCMC_Tree_Height(tree);
      	}

      /* Subtree height */
      else if(!strcmp(tree->mcmc->move_name[move],"subtree_height"))
      	{
      	  MCMC_Subtree_Height(tree);
      	}

      /* Subtree rates */
      else if(!strcmp(tree->mcmc->move_name[move],"subtree_rates"))
      	{
      	  MCMC_Subtree_Rates(tree);
      	}

      /* Birth rate */
      else if(!strcmp(tree->mcmc->move_name[move],"birth_rate"))
      	{
      	  MCMC_Birth_Rate(tree);
      	}


      /* Swing rates */
      else if(!strcmp(tree->mcmc->move_name[move],"tree_rates"))
      	{
      	  MCMC_Tree_Rates(tree);
      	}

      else if(!strcmp(tree->mcmc->move_name[move],"updown_nu_cr"))
      	{
      	  MCMC_Updown_Nu_Cr(tree);
      	}

      else if(!strcmp(tree->mcmc->move_name[move],"updown_t_cr"))
      	{
      	  MCMC_Updown_T_Cr(tree);
      	}

      /* Ts/tv ratio */
      else if(!strcmp(tree->mcmc->move_name[move],"kappa"))
      	{
      	  MCMC_Kappa(tree);
	}

      /* Gamma shape parameter */
      else if(!strcmp(tree->mcmc->move_name[move],"ras"))
      	{
      	  MCMC_Rate_Across_Sites(tree);
	}

      /* Covarion model parameters */
      else if(!strcmp(tree->mcmc->move_name[move],"cov_rates"))
      	{
      	  MCMC_Covarion_Rates(tree);
	}

      /* Covarion model parameters */
      else if(!strcmp(tree->mcmc->move_name[move],"cov_switch"))
      	{
      	  MCMC_Covarion_Switch(tree);
	}


      /* Times */
      else if(!strcmp(tree->mcmc->move_name[move],"time"))
      	{
      	  tree->both_sides = YES;
      	  if(tree->mcmc->use_data == YES) Lk(tree);
      	  tree->both_sides = NO;


      	  if(tree->mcmc->is == NO || tree->rates->model_log_rates == YES)
      	    {
	      MCMC_One_Time(tree->n_root,tree->n_root->v[first],YES,tree);
	      MCMC_One_Time(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	  else
      	    {
	      /* MCMC_One_Time(tree->n_root,tree->n_root->v[first],YES,tree); */
	      /* MCMC_One_Time(tree->n_root,tree->n_root->v[secod],YES,tree); */
      	      RATES_Posterior_One_Time(tree->n_root,tree->n_root->v[first],YES,tree);
      	      RATES_Posterior_One_Time(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	}
      
      /* Node Rates */

      else if(!strcmp(tree->mcmc->move_name[move],"nd_rate"))
      	{
      	  MCMC_One_Node_Rate(tree->n_root,tree->n_root->v[first],YES,tree);
      	  MCMC_One_Node_Rate(tree->n_root,tree->n_root->v[secod],YES,tree);
      	}

      /* Edge Rates */
      else if(!strcmp(tree->mcmc->move_name[move],"br_rate"))
      	{

      	  tree->both_sides = YES;
      	  if(tree->mcmc->use_data == YES) Lk(tree);
      	  tree->both_sides = NO;

      	  if(tree->mcmc->is == NO)
      	    {
      	      /* MCMC_Slice_One_Rate(tree->n_root,tree->n_root->v[first],YES,tree); */
      	      /* MCMC_Slice_One_Rate(tree->n_root,tree->n_root->v[secod],YES,tree); */

	      MCMC_One_Rate(tree->n_root,tree->n_root->v[first],YES,tree);
	      MCMC_One_Rate(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }
      	  else
      	    {
      	      RATES_Posterior_One_Rate(tree->n_root,tree->n_root->v[first],YES,tree);
      	      RATES_Posterior_One_Rate(tree->n_root,tree->n_root->v[secod],YES,tree);
      	    }

	  /* MCMC_Sim_Rate(tree->n_root,tree->n_root->v[0],tree); */
	  /* MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree); */
      	  /* if(tree->mcmc->use_data == YES) Lk(tree); */
	  /* RATES_Lk_Rates(tree); */
      	}

      tree->mcmc->run++;
      MCMC_Get_Acc_Rates(tree->mcmc);

      MCMC_Print_Param(tree->mcmc,tree);
      MCMC_Print_Param_Stdin(tree->mcmc,tree);

      if(tree->io->mutmap == YES)
	{
	  if(!(tree->mcmc->run%tree->mcmc->sample_interval)) 
	    {
	      Sample_Ancestral_Seq(YES,!tree->mcmc->use_data,tree);
	      
	      phydbl sum = 0.0;
	      int edge,site,mut;
	      char *s;
	      FILE *fp;
	      
	      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));
	      
	      strcpy(s,tree->mcmc->io->in_align_file);
	      strcat(s,"_");
	      strcat(s,tree->mcmc->out_filename);
	      strcat(s,".mutmap");
	      fp = fopen(s,"w");
	      
	      Free(s);
	      
	      For(i,(2*tree->n_otu-3)*(tree->n_pattern)*6) sum += tree->mutmap[i];
	      PhyML_Fprintf(fp,"edge\t site\t mut\t count");
	      For(i,(2*tree->n_otu-3)*(tree->n_pattern)*6) 
		{
		  Get_Mutmap_Coord(i,&edge,&site,&mut,tree);
		  PhyML_Fprintf(fp,"\n%4d\t %4d\t %4d\t %10f",edge,site,mut,(phydbl)tree->mutmap[i]/sum);
		}
	      
	      fclose(fp);	      
	    }
	}

      (void)signal(SIGINT,MCMC_Terminate); 
    }
  while(tree->mcmc->run < tree->mcmc->chain_len);


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Single_Param_Generic(phydbl *val, 
			       phydbl lim_inf, 
			       phydbl lim_sup, 
			       int move_num,
			       phydbl *lnPrior,
			       phydbl *lnLike,
			       phydbl (*prior_func)(t_edge *,t_tree *,supert_tree *), 
			       phydbl (*like_func)(t_edge *,t_tree *,supert_tree *),
			       int move_type,
			       int _log, /* _log == YES: the model describes the distribution of log(val) but the move applies to val. Need a correction factor */
			       t_edge *branch, t_tree *tree, supert_tree *stree)
{
  phydbl cur_val,new_val,new_lnLike,new_lnPrior,cur_lnLike,cur_lnPrior;
  phydbl u,alpha,ratio;
  phydbl K,mult;
  phydbl new_lnval, cur_lnval;
 
  Record_Br_Len(tree);

  cur_val       = *val;
  new_val       = -1.0;
  ratio         =  0.0;
  mult          = -1.0;
  K             = tree->mcmc->tune_move[move_num];
  cur_lnval     = LOG(*val);
  new_lnval     = cur_lnval;

  if(lnLike)
    {
      cur_lnLike = *lnLike;
      new_lnLike = *lnLike;
    }
  else
    {
      cur_lnLike = 0.0;
      new_lnLike = 0.0;
    }
  
  if(lnPrior)
    {
      cur_lnPrior = *lnPrior;
      new_lnPrior = *lnPrior;
    }
  else
    {
      cur_lnPrior = 0.0;
      new_lnPrior = 0.0;
    }

  MCMC_Make_Move(&cur_val,&new_val,lim_inf,lim_sup,&ratio,K,move_type);

  if(new_val < lim_sup && new_val > lim_inf)
    {
      *val = new_val;
  
      RATES_Update_Cur_Bl(tree);
            
      if(_log == YES) ratio += (cur_lnval - new_lnval);
      
      if(prior_func) /* Prior ratio */
	{ 
	  new_lnPrior = (*prior_func)(branch,tree,stree); 
	  ratio += (new_lnPrior - cur_lnPrior); 
	}
      
      if(like_func && tree->mcmc->use_data == YES)  /* Likelihood ratio */
	{ 
	  new_lnLike  = (*like_func)(branch,tree,stree);  
	  ratio += (new_lnLike - cur_lnLike);  
	}
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      if(u > alpha) /* Reject */
	{
	  *val    = cur_val;
	  new_val = cur_val;
	  if(lnPrior) *lnPrior = cur_lnPrior;
	  if(lnLike)  *lnLike  = cur_lnLike;
	  Restore_Br_Len(tree);
	  if(tree->mod->update_eigen) Update_Eigen(tree->mod);
	}
      else /* Accept */
	{
	  tree->mcmc->acc_move[move_num]++;
	  if(lnPrior) *lnPrior = new_lnPrior;
	  if(lnLike)  *lnLike  = new_lnLike;
	}
    }
      
  tree->mcmc->run_move[move_num]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Update_Effective_Sample_Size(int move_num, t_mcmc *mcmc, t_tree *tree)
{
  phydbl new_val,cur_val;

  mcmc->ess_run[move_num]++;

  new_val = mcmc->new_param_val[move_num];
  cur_val = mcmc->old_param_val[move_num];

  if(mcmc->ess_run[move_num] == 1)
    {
      mcmc->first_val[move_num] = cur_val;
      mcmc->sum_val[move_num]   = cur_val;
      mcmc->sum_valsq[move_num] = POW(cur_val,2);
      return;
    }

  mcmc->sum_val[move_num]            += new_val;
  mcmc->sum_valsq[move_num]          += POW(new_val,2);
  mcmc->sum_curval_nextval[move_num] += cur_val * new_val;
  	
  mcmc->ess[move_num] = 
    Effective_Sample_Size(mcmc->first_val[move_num],
			  new_val,
			  mcmc->sum_val[move_num],
			  mcmc->sum_valsq[move_num],
			  mcmc->sum_curval_nextval[move_num],
			  mcmc->ess_run[move_num]);
	
  mcmc->old_param_val[move_num] = new_val;

  if(move_num == mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu)
    {
      /* FILE *fp; */
      /* fp = fopen("out","a"); */
      /* fprintf(fp,"%f\n",new_val); */
      /* fclose(fp); */
      /* printf("\n. first=%G last=%G sum=%f sum_valsq=%f sum_cur_next=%f run=%d ess=%f", */
      /* 	     mcmc->first_val[move_num],new_val, */
      /* 	     mcmc->sum_val[move_num], */
      /* 	     mcmc->sum_valsq[move_num], */
      /* 	     mcmc->sum_curval_nextval[move_num], */
      /* 	     mcmc->run_move[move_num]+1, */
      /* 	     mcmc->ess[move_num] */
      /* 	     ); */
    }
      
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Sample_Joint_Rates_Prior(t_tree *tree)
{
  int i,dim;
  phydbl T;
  phydbl *r,*t,*lambda;
  phydbl *min_r,*max_r;
  phydbl k;

  dim    = 2*tree->n_otu-2;
  lambda = tree->rates->_2n_vect1;
  min_r  = tree->rates->_2n_vect2;
  max_r  = tree->rates->_2n_vect3;
  r      = tree->rates->br_r;
  t      = tree->rates->nd_t;

  For(i,dim) tree->rates->mean_r[i] = 1.0;

  RATES_Fill_Lca_Table(tree);
  RATES_Covariance_Mu(tree);

  T = .0;
  For(i,dim) T += (t[tree->t_nodes[i]->num] - t[tree->t_nodes[i]->anc->num]);
  For(i,dim) lambda[i] = (t[tree->t_nodes[i]->num] - t[tree->t_nodes[i]->anc->num])/T;
  For(i,dim) r[i] = 1.0;
  For(i,dim) min_r[i] = tree->rates->min_rate;
  For(i,dim) max_r[i] = tree->rates->max_rate;

  k = 1.; /* We want \sum_i lambda[i] r[i] = 1 */

  Rnorm_Multid_Trunc_Constraint(tree->rates->mean_r, 
				tree->rates->cov_r, 
				min_r,max_r, 
				lambda,
				k, 
				r,
				dim);
}
  
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  t_edge *b;
  int i;
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data, new_lnL_rate, cur_lnL_rate;
  phydbl ratio, alpha;
  phydbl new_mu, cur_mu;
  phydbl r_min, r_max;
  phydbl K;
  t_edge *b1,*b2,*b3;
  t_node *v2,*v3;
  phydbl nu,cr;
  
  if(tree->rates->model == STRICTCLOCK) return;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else
    For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }
   
  cur_mu       = tree->rates->br_r[d->num];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  r_min        = tree->rates->min_rate;
  r_max        = tree->rates->max_rate;
  K            = tree->mcmc->tune_move[tree->mcmc->num_move_br_r+d->num];
  ratio        = 0.0;
  nu           = tree->rates->nu;
  cr           = tree->rates->clock_r;
  
  Record_Br_Len(tree);
  
  u = Uni();
  
  /* MCMC_Make_Move(&cur_mu,&new_mu,r_min,r_max,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_br_r+d->num]); */

  phydbl dt,sd,mean;
  int err;
  dt     = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  sd     = SQRT(tree->rates->nu * dt);
  mean   = tree->rates->br_r[a->num] - 0.5*sd*sd;
  new_mu = Rnorm_Trunc(mean,sd,r_min,r_max,&err);
  ratio  = Log_Dnorm_Trunc(cur_mu,mean,sd,r_min,r_max,&err) - Log_Dnorm_Trunc(new_mu,mean,sd,r_min,r_max,&err);

  if(new_mu > r_min && new_mu < r_max)
    {      
      tree->rates->br_r[d->num] = new_mu;
            
      v2 = v3 = NULL;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    if(!v2) { v2 = d->v[i]; }
	    else    { v3 = d->v[i]; }
	  }
      
      
      b1 = NULL;
      if(a == tree->n_root) b1 = tree->e_root;
      else For(i,3) if(d->v[i] == a) { b1 = d->b[i]; break; }
      
      b2 = b3 = NULL;
      if(!d->tax)
	{
	  For(i,3)
	    if((d->v[i] != a) && (d->b[i] != tree->e_root))
	      {
		if(!b2) { b2 = d->b[i]; }
		else    { b3 = d->b[i]; }
	      }
	}
      
      tree->rates->br_do_updt[d->num] = YES;
      if(!d->tax)
      	{
      	  tree->rates->br_do_updt[v2->num] = YES;
      	  tree->rates->br_do_updt[v3->num] = YES;
      	}
      
      RATES_Update_Cur_Bl(tree);
      
      /* printf("\n. r0=%f r1=%f cr=%f mean=%f var=%f nu=%f dt=%f", */
      /* 	 r0,r1,tree->rates->clock_r,b1->gamma_prior_mean,b1->gamma_prior_var,nu,t1-t0); */
      
      
      if(tree->mcmc->use_data)
	{
	  if(tree->io->lk_approx == EXACT)
	    {
	      Update_PMat_At_Given_Edge(b1,tree);
	      if(!d->tax)
		{
		  Update_PMat_At_Given_Edge(b2,tree);
		  Update_PMat_At_Given_Edge(b3,tree);
		}
	      Update_P_Lk(tree,b1,d);
	    }
	  new_lnL_data = Lk_At_Given_Edge(b1,tree);
	  
	  /* tree->both_sides = NO; */
	  /* new_lnL_data = Lk(tree); */
	}
      
      new_lnL_rate = RATES_Lk_Rates(tree);
      
      /* Likelihood ratio */
      if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

      /* Prior ratio */
      ratio += (new_lnL_rate - cur_lnL_rate);
      
      
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      
      if(u > alpha) /* Reject */
	{
	  tree->rates->br_r[d->num] = cur_mu;	  
	  tree->c_lnL        = cur_lnL_data;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  
	  Restore_Br_Len(tree);
	  
	  if(tree->mcmc->use_data && tree->io->lk_approx == EXACT)
	    {
	      Update_PMat_At_Given_Edge(b1,tree);
	      if(!d->tax)
		{
		  Update_PMat_At_Given_Edge(b2,tree);
		  Update_PMat_At_Given_Edge(b3,tree);
		}
	      Update_P_Lk(tree,b1,d);
	    }
	  
	  /* tree->both_sides = YES; */
	  /* new_lnL_data = Lk(tree); */
	  /* tree->both_sides = NO; */
	  
	}
      else
	{
	  tree->mcmc->acc_move[tree->mcmc->num_move_br_r+d->num]++;
	}
    }
  tree->mcmc->run_move[tree->mcmc->num_move_br_r+d->num]++;
  
  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,d->b[i],d);
		/* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
		MCMC_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,b,d);
      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_One_Node_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  t_edge *b;
  int i;

  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else
    For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }

  /* Only the LOG_RANDWALK move seems to work here. Change with caution then. */
  tree->rates->br_do_updt[d->num] = YES;
  MCMC_Single_Param_Generic(&(tree->rates->nd_r[d->num]),
			    tree->rates->min_rate,
			    tree->rates->max_rate,
			    tree->mcmc->num_move_nd_r+d->num,
			    &(tree->rates->c_lnL_rates),NULL,
			    Wrap_Lk_Rates,NULL,
			    tree->mcmc->move_type[tree->mcmc->num_move_nd_r+d->num],
			    NO,NULL,tree,NULL);


  Update_PMat_At_Given_Edge(b,tree);

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		MCMC_One_Node_Rate(d,d->v[i],YES,tree);
	      }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_One_Time(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl u;
  phydbl t_min,t_max;
  phydbl t1_cur, t1_new;
  phydbl cur_lnL_data, new_lnL_data;
  phydbl cur_lnL_rate, new_lnL_rate;
  phydbl cur_lnL_time, new_lnL_time;
  phydbl ratio,alpha;
  t_edge *b1,*b2,*b3;
  int    i;
  phydbl t0,t2,t3;
  t_node *v2,*v3;
  phydbl K;
  int move_num;
  
  
  if(d->tax) return; /* Won't change time at tip */
  if(d == tree->n_root) return; /* Won't change time at root. Use tree height move instead */

  if(FABS(tree->rates->t_prior_min[d->num] - tree->rates->t_prior_max[d->num]) < 1.E-10) return;

  Record_Br_Len(tree);
  RATES_Record_Rates(tree);
  RATES_Record_Times(tree);
  
  move_num       = d->num-tree->n_otu+tree->mcmc->num_move_nd_t;
  K              = tree->mcmc->tune_move[move_num];
  cur_lnL_data   = tree->c_lnL;
  cur_lnL_rate   = tree->rates->c_lnL_rates;
  t1_cur         = tree->rates->nd_t[d->num];
  new_lnL_data   = cur_lnL_data;
  new_lnL_rate   = cur_lnL_rate;
  ratio          = 0.0;
  cur_lnL_time   = tree->rates->c_lnL_times;
  new_lnL_time   = cur_lnL_time;

  v2 = v3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!v2) { v2 = d->v[i]; }
	else    { v3 = d->v[i]; }
      }


  b1 = NULL;
  if(a == tree->n_root) b1 = tree->e_root;
  else For(i,3) if(d->v[i] == a) { b1 = d->b[i]; break; }

  b2 = b3 = NULL;
  For(i,3)
    if((d->v[i] != a) && (d->b[i] != tree->e_root))
      {
	if(!b2) { b2 = d->b[i]; }
	else    { b3 = d->b[i]; }
      }
  
  t0 = (a)?(tree->rates->nd_t[a->num]):(tree->rates->t_prior_min[d->num]);
  t2 = tree->rates->nd_t[v2->num];
  t3 = tree->rates->nd_t[v3->num];

  t_min = MAX(t0,tree->rates->t_prior_min[d->num]);
  t_max = MIN(MIN(t2,t3),tree->rates->t_prior_max[d->num]);

  t_min += tree->rates->min_dt;
  t_max -= tree->rates->min_dt;

  if(t_min > t_max) 
    {
      PhyML_Printf("\n. t_min = %f t_max = %f",t_min,t_max);
      PhyML_Printf("\n. prior_min = %f prior_max = %f",tree->rates->t_prior_min[d->num],tree->rates->t_prior_max[d->num]);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      /* Exit("\n"); */
    }

  MCMC_Make_Move(&t1_cur,&t1_new,t_min,t_max,&ratio,K,tree->mcmc->move_type[move_num]);

  if(t1_new > t_min && t1_new < t_max) 
    {
      tree->rates->nd_t[d->num] = t1_new;

      /* Update branch lengths */
      tree->rates->br_do_updt[d->num]  = YES;
      tree->rates->br_do_updt[v2->num] = YES;
      tree->rates->br_do_updt[v3->num] = YES;
      RATES_Update_Cur_Bl(tree);

      if(tree->mcmc->use_data)
	{
	  if(tree->io->lk_approx == EXACT)
	    {
	      Update_PMat_At_Given_Edge(b1,tree);
	      Update_PMat_At_Given_Edge(b2,tree);
	      Update_PMat_At_Given_Edge(b3,tree);
	      Update_P_Lk(tree,b1,d);
	    }
	  new_lnL_data = Lk_At_Given_Edge(b1,tree);

	  /* /\* !!!!!!!!!!!!!!!!!!!!1 *\/ */
	  /* if(FABS(Lk(tree) - new_lnL_data) > 1.E-5) */
	  /*   { */
	  /*     PhyML_Printf("\n. b1->l=%f b2->l=%f b3->l=%f",b1->l,b2->l,b3->l); */
	  /*     PhyML_Printf("\n. a->num=%d d->num=%d root=%d (%d %d)",a->num,d->num,a==tree->n_root,tree->n_root->v[0]->num,tree->n_root->v[1]->num); */
	  /*     PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max); */
	  /*     PhyML_Printf("\n. t1_new=%f t1_cur=%f",t1_new,t1_cur); */
	  /*     PhyML_Printf("\n. %f %f",cur_lnL_data,tree->c_lnL); */
	  /*     PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
	  /*     Warn_And_Exit(""); */
	  /*   } */
	}

      new_lnL_rate = RATES_Lk_Rates(tree);
      new_lnL_time = TIMES_Lk_Times(tree);

      if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);
      ratio += (new_lnL_rate - cur_lnL_rate);
      ratio += (new_lnL_time - cur_lnL_time);

      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      u = Uni();
	           
      if(u > alpha) /* Reject */
	{
	  tree->rates->nd_t[d->num] = t1_cur;
	  tree->c_lnL              = cur_lnL_data;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->rates->c_lnL_times = cur_lnL_time;
	  Restore_Br_Len(tree);
	  RATES_Reset_Rates(tree);
	  RATES_Reset_Times(tree);

	  if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) 
	    {
	      Update_PMat_At_Given_Edge(b1,tree);
	      Update_PMat_At_Given_Edge(b2,tree);
	      Update_PMat_At_Given_Edge(b3,tree);
	      Update_P_Lk(tree,b1,d);
	    }
	  /* printf("\n. R new_lnL_data = %f cur_lnL_data = %f t_new=%f t_cur=%f tmin=%f tmax=%f", */
	  /* 	 new_lnL_data,cur_lnL_data,t1_new,t1_cur,t_min,t_max); */
	}
      else
	{
	  /* printf("\n. A new_lnL_data = %f cur_lnL_data = %f t_new=%f t_cur=%f tmin=%f tmax=%f", */
	  /* 	 new_lnL_data,cur_lnL_data,t1_new,t1_cur,t_min,t_max); */
	  tree->mcmc->acc_move[move_num]++;
	}
      
      if(t1_new < t0)
	{
	  t1_new = t0+1.E-4;
	  PhyML_Printf("\n");
	  PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
	  PhyML_Printf("\n. t0 = %f t1_new = %f",t0,t1_new);
	  PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      if(t1_new > MIN(t2,t3))
	{
	  PhyML_Printf("\n");
	  PhyML_Printf("\n. a is root -> %s",(a == tree->n_root)?("YES"):("NO"));
	  PhyML_Printf("\n. t0 = %f t1_new = %f t1 = %f t2 = %f t3 = %f MIN(t2,t3)=%f",t0,t1_new,t1_cur,t2,t3,MIN(t2,t3));
	  PhyML_Printf("\n. t_min=%f t_max=%f",t_min,t_max);
	  PhyML_Printf("\n. (t1-t0)=%f (t2-t1)=%f",t1_cur-t0,t2-t1_cur);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
      
      if(isnan(t1_new))
	{
	  PhyML_Printf("\n. run=%d",tree->mcmc->run);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  /*       Exit("\n"); */
	}
    }
  
  tree->mcmc->run_move[move_num]++;

  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	For(i,3)
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,d->b[i],d);
	      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
	      MCMC_One_Time(d,d->v[i],YES,tree);
	    }
      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,b1,d);
      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) {tree->both_sides = YES; Lk(tree); } */
    }	    
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Tree_Height(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;
  phydbl cur_height, new_height;
  phydbl min_t,max_t;

  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_tree_height];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_height   = tree->rates->nd_t[tree->n_root->num];
  cur_lnL_time = tree->rates->c_lnL_times;
  min_t        = tree->rates->t_prior_min[tree->n_root->num];
  max_t        = tree->rates->t_prior_max[tree->n_root->num];

  u = Uni();
  mult = EXP(K*(u-0.5));

 /* WARNING: It must not be floor = tree->rates->t_prior_max[tree->n_root->num]; 
     floor is the maximum value a node height can take when one ignores the 
     calibration nodes, i.e., floor is set by the height of the tips
  */
  /* floor = 0.0; */
  floor = tree->rates->t_floor[tree->n_root->num];
  /* floor = tree->rates->t_prior_max[tree->n_root->num]; */

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);
  new_height = tree->rates->nd_t[tree->n_root->num];
  
  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
  	  Restore_Br_Len(tree);
  	  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
	  return;
  	}
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  /* ratio += LOG(mult); */
  ratio += (phydbl)(n_nodes)*LOG(mult);

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);


  /* printf("\n. [%f %f] [%f %f] [%f %f] [%f]",  */
  /* 	 new_lnL_rate , cur_lnL_rate, */
  /* 	 new_lnL_time , cur_lnL_time, */
  /* 	 cur_height,new_height,mult); */

  /* phydbl diff = (new_lnL_time - cur_lnL_time)+(phydbl)(n_nodes-1.)*LOG(mult); */
  /* printf("\n. diff_lk = %12f -(n-1)*log(mult) = %12f n_nodes=%d [%12f] log(mult)=%f", */
  /* 	 (new_lnL_time - cur_lnL_time), */
  /* 	 -(phydbl)(n_nodes-1.)*LOG(mult), */
  /* 	 n_nodes, */
  /* 	 diff, */
  /* 	 LOG(mult)); */

  /* !!!!!!!!!!!! */
  /* ratio += LOG(Dexp(FABS(new_height-floor),1./10.) / Dexp(FABS(cur_height-floor),1./10.)); */
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_height]++;
      tree->mcmc->acc_move[tree->mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_height]++;
  tree->mcmc->run_move[tree->mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Updown_T_Cr(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl floor;
  int n_nodes;

  if(FABS(tree->rates->t_prior_max[tree->n_root->num] - tree->rates->t_prior_min[tree->n_root->num]) < 1.E-10) return;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K            = tree->mcmc->tune_move[tree->mcmc->num_move_updown_t_cr];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  ratio        = 0.0;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  cur_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = EXP(K*(u-0.5));

  floor = 0.0;

  Scale_Subtree_Height(tree->n_root,mult,floor,&n_nodes,tree);

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  Restore_Br_Len(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
  	  return;
  	}
    }

  if(RATES_Check_Node_Times(tree))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->rates->clock_r /= mult;
  if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
    {
      tree->rates->clock_r *= mult;
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
      return;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio is actually mult^(n) when changing the absolute
     node heights. When considering the relative heights, this ratio combined
     to the Jacobian for the change of variable ends up to being equal to mult. 
  */
  ratio += (n_nodes - 1)*LOG(mult);
  /* ratio += -LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);

  /* !!!!!!!!!!!!1 */
  /* ratio += LOG(Dexp(FABS(new_height-floor),1./10.) / Dexp(FABS(cur_height-floor),1./10.)); */
  
  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      tree->rates->clock_r *= mult;
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_updown_t_cr]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_updown_t_cr]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Subtree_Height(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_time,new_lnL_time;
  phydbl cur_height,new_height;
  phydbl floor;
  int target;
  int n_nodes;

  RATES_Record_Times(tree);
  Record_Br_Len(tree);

  K = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_height];
  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;
  ratio        = 0.0;
  cur_lnL_time = tree->rates->c_lnL_times;

  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);

  cur_height = tree->rates->nd_t[target];
  floor = tree->rates->t_floor[target];

  if(tree->t_nodes[target] == tree->n_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!Scale_Subtree_Height(tree->t_nodes[target],mult,floor,&n_nodes,tree))
    {
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
      return;
    }

  new_height = tree->rates->nd_t[target];

  For(i,2*tree->n_otu-1)
    {
      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
  	 tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
  	{
  	  RATES_Reset_Times(tree);
	  Restore_Br_Len(tree);
	  tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;
  	  return;
  	}
    }

     
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  new_lnL_rate = RATES_Lk_Rates(tree);
  new_lnL_time = TIMES_Lk_Times(tree);

  /* The Hastings ratio here is mult^(n_nodes) and the ratio of the prior joint densities
     of the modified node heigths given the unchanged one is 1. This is different from the 
     case where all the nodes, including the root node, are scaled. 
  */
  ratio += (phydbl)(n_nodes)*LOG(mult);

  /* Likelihood ratio */
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  /* Prior ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  ratio += (new_lnL_time - cur_lnL_time);


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->rates->c_lnL_times = cur_lnL_time;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_subtree_height]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_subtree_height]++;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Tree_Rates(t_tree *tree)
{
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_data,new_lnL_data;
  int n_nodes;
  phydbl init_clock;
  
  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_tree_rates];
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = tree->rates->c_lnL_rates;
  init_clock    = tree->rates->clock_r;
  ratio         = 0.0;

  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */
    
  /* Multiply branch rates (or add to log of rates) */
  if(!Scale_Subtree_Rates(tree->n_root,mult,&n_nodes,tree))
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;
      return;
    }

  if(n_nodes != 2*tree->n_otu-2)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  /* Divide clock_r */
  tree->rates->clock_r /= mult;
  if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
    {
      tree->rates->clock_r = init_clock;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;
      return;
    }

/*   if(tree->rates->model == GUINDON) */
/*     { */
      int i;
      For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
      RATES_Update_Cur_Bl(tree);
      if(tree->mcmc->use_data) new_lnL_data = Lk(tree);
/*     } */

  new_lnL_rate = RATES_Lk_Rates(tree);


  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  ratio += (+(2*tree->n_otu-2)-1)*LOG(mult);
  /* ratio += (+(2*tree->n_otu-2)-1-2)*LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */

  /* If modelling log of rates instead of rates */
  if(tree->rates->model_log_rates == YES) ratio -= (2*tree->n_otu-2)*LOG(mult);

  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);

  /* Likelihood density ratio */
  ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
/*   printf("\n. cur_lnL=%f new_lnL=%f ratio=%f mult=%f %f [%f %f]", */
/*   	 cur_lnL_rate,new_lnL_rate,ratio,mult,(+(2*tree->n_otu-2)-1-2)*LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)),new_lnL_data,cur_lnL_data); */

  if(u > alpha)
    {
/*       PhyML_Printf("\n. Reject mult=%f",mult); */
      tree->rates->clock_r = init_clock;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_data;
    }
  else
    {
/*       PhyML_Printf("\n. Accept mult=%f",mult); */
      tree->mcmc->acc_move[tree->mcmc->num_move_tree_rates]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_tree_rates]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Subtree_Rates(t_tree *tree)
{  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_data,new_lnL_data;
  int target;
  int n_nodes;
  
  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_subtree_rates];
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = cur_lnL_rate;
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = cur_lnL_data;
  ratio         = 0.0;

  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = Rgamma(1./K,K); */

  target = Rand_Int(tree->n_otu,2*tree->n_otu-3);

  /* Multiply branch rates */
  if(!Scale_Subtree_Rates(tree->t_nodes[target],mult,&n_nodes,tree))
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_subtree_rates]++;
      return;
    }

  new_lnL_rate = RATES_Lk_Rates(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  ratio += (+n_nodes)*LOG(mult);
  /* ratio += (n_nodes-2)*LOG(mult) + LOG(Dgamma(1./mult,1./K,K)/Dgamma(mult,1./K,K)); */
  
  /* If modelling log of rates instead of rates */
  if(tree->rates->model_log_rates == YES) ratio -= (n_nodes)*LOG(mult);

  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);

  /* Likelihood density ratio */
  ratio += (new_lnL_data - cur_lnL_data);


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_subtree_rates]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_subtree_rates]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Swing(t_tree *tree)
{
  int i;
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_data,new_lnL_data;
  phydbl cur_lnL_rate,new_lnL_rate;

  if(tree->rates->model == STRICTCLOCK) return;

  RATES_Record_Times(tree);
  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = 3.;
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = cur_lnL_data;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = cur_lnL_rate;
  ratio         = 0.0;

  u = Uni();
  /* mult = EXP(K*(u-0.5)); */
  mult = u*(K - 1./K) + 1./K;


  For(i,2*tree->n_otu-1)
    {
      if(tree->t_nodes[i]->tax == NO) 
	{
	  tree->rates->nd_t[i] *= mult;
	}

      if(tree->rates->nd_t[i] > tree->rates->t_prior_max[i] ||
         tree->rates->nd_t[i] < tree->rates->t_prior_min[i])
        {
          RATES_Reset_Times(tree);
          Restore_Br_Len(tree);
          return;
        }
    }

  For(i,2*tree->n_otu-2)
    {
      tree->rates->br_r[i] /= mult;

      if(tree->rates->br_r[i] > tree->rates->max_rate ||
         tree->rates->br_r[i] < tree->rates->min_rate)
        {
          RATES_Reset_Times(tree);
          RATES_Reset_Rates(tree);
          Restore_Br_Len(tree);
          return;
        }
    }

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio += (-(tree->n_otu-1.)-2.)*LOG(mult);
  ratio += (new_lnL_rate - cur_lnL_rate);
  if(tree->mcmc->use_data) ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();

  if(u > alpha)
    {
      RATES_Reset_Times(tree);
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
      tree->rates->c_lnL_rates = cur_lnL_rate;
      /* printf("\n. Reject %8f",mult); */
    }
  else
    {
      /* printf("\n. Accept %8f",mult); */
    }


}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Updown_Nu_Cr(t_tree *tree)
{
  phydbl K,mult,u,alpha,ratio;
  phydbl cur_lnL_rate,new_lnL_rate;
  phydbl cur_lnL_data,new_lnL_data;
  int i;
  int n_nodes;

  RATES_Record_Rates(tree);
  Record_Br_Len(tree);

  K             = tree->mcmc->tune_move[tree->mcmc->num_move_updown_nu_cr];
  cur_lnL_data  = tree->c_lnL;
  new_lnL_data  = tree->c_lnL;
  cur_lnL_rate  = tree->rates->c_lnL_rates;
  new_lnL_rate  = tree->rates->c_lnL_rates;
  n_nodes       = 0;

  u = Uni();
  mult = EXP(K*(u-0.5));

  /* Multiply branch rates */
  /* if(!Scale_Subtree_Rates(tree->n_root,mult,&n_nodes,tree)) */
  /*   { */
  /*     RATES_Reset_Rates(tree); */
  /*     Restore_Br_Len(tree); */
  /*     tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++; */
  /*     return; */
  /*   } */

  tree->rates->clock_r /= mult;
  if(tree->rates->clock_r < tree->rates->min_clock || tree->rates->clock_r > tree->rates->max_clock)
    {
      tree->rates->clock_r *= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
      return;
    }

  tree->rates->nu *= mult;
  if(tree->rates->nu < tree->rates->min_nu || tree->rates->nu > tree->rates->max_nu)
    {
      tree->rates->nu      /= mult;
      tree->rates->clock_r *= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
      return;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
  RATES_Update_Cur_Bl(tree);
  if(tree->mcmc->use_data) new_lnL_data = Lk(tree);

  new_lnL_rate = RATES_Lk_Rates(tree);

  ratio = 0.0;
  /* Proposal ratio: 2n-2=> number of multiplications, 1=>number of divisions */
  /* ratio += n_nodes*LOG(mult); /\* (1-1)*LOG(mult); *\/ */
  ratio += 0.0*LOG(mult); /* (1-1)*LOG(mult); */
  /* Prior density ratio */
  ratio += (new_lnL_rate - cur_lnL_rate);
  /* Likelihood density ratio */
  ratio += (new_lnL_data - cur_lnL_data);

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  u = Uni();
  
  if(u > alpha)
    {
      tree->rates->clock_r *= mult;
      tree->rates->nu /= mult;
      RATES_Reset_Rates(tree);
      Restore_Br_Len(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
      tree->c_lnL        = cur_lnL_data;
    }
  else
    {
      tree->mcmc->acc_move[tree->mcmc->num_move_updown_nu_cr]++;
    }

  tree->mcmc->run_move[tree->mcmc->num_move_updown_nu_cr]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Param_Stdin(t_mcmc *mcmc, t_tree *tree)
{
  time_t cur_time;
  phydbl min;
  int i;
  time(&cur_time);
  
  min = MDBL_MAX;
  For(i,tree->n_otu-1)
    {
      /*       printf("\n. %d %f %f %f %f",i, */
      /* 	     tree->mcmc->new_param_val[tree->mcmc->num_move_nd_t+i], */
      /* 	     tree->mcmc->old_param_val[tree->mcmc->num_move_nd_t+i], */
      /* 	     tree->mcmc->ess[tree->mcmc->num_move_nd_t+i], */
      /* 	     tree->mcmc->sum_val[tree->mcmc->num_move_nd_t+i]); */

      if(tree->mcmc->ess[tree->mcmc->num_move_nd_t+i] < min)
	min = tree->mcmc->ess[tree->mcmc->num_move_nd_t+i];
    }


  if(mcmc->run == 1)
    {
      PhyML_Printf("\n\n");
      PhyML_Printf("\t%7s","Run");
      PhyML_Printf("\t%5s","Time");
      PhyML_Printf("\t%10s","Likelihood");
      PhyML_Printf("\t%10s","Prior");
      PhyML_Printf("\t%17s","SubstRate[ESS]");
      PhyML_Printf("\t%15s","TreeHeight[ESS]");    
      if(tree->rates->model == THORNE || tree->rates->model == GUINDON) PhyML_Printf("\t%14s","AutoCor[ESS]");    
      else PhyML_Printf("\t%14s","RateVar[ESS]");
      PhyML_Printf("\t%13s","BirthR[ESS]");
      PhyML_Printf("\t%7s","MinESS");    
    }

  if((cur_time - mcmc->t_last_print) >  mcmc->print_every)
    {
      mcmc->t_last_print = cur_time;
      PhyML_Printf("\n");
      PhyML_Printf("\t%7d",tree->mcmc->run);
      PhyML_Printf("\t%5d",(int)(cur_time-mcmc->t_beg));
      PhyML_Printf("\t%10.2f",tree->c_lnL);
      PhyML_Printf("\t%10.2f",tree->rates->c_lnL_rates+tree->rates->c_lnL_times);
      PhyML_Printf("\t%12.6f[%3.0f]",RATES_Average_Substitution_Rate(tree),tree->mcmc->ess[tree->mcmc->num_move_clock_r]);
      /* PhyML_Printf("\t%12.6f[%4.0f]",tree->rates->clock_r,tree->mcmc->ess[tree->mcmc->num_move_clock_r]); */
      PhyML_Printf("\t%10.1f[%3.0f]",tree->rates->nd_t[tree->n_root->num],tree->mcmc->ess[tree->mcmc->num_move_nd_t+tree->n_root->num-tree->n_otu]);
      PhyML_Printf("\t%9.4G[%3.0f]",tree->rates->nu,tree->mcmc->ess[tree->mcmc->num_move_nu]);
      PhyML_Printf("\t%8.4G[%3.0f]",tree->rates->birth_rate,tree->mcmc->ess[tree->mcmc->num_move_birth_rate]);
      PhyML_Printf("\t%7.0f",min);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Param(t_mcmc *mcmc, t_tree *tree)
{
  int i;
  FILE *fp;
  char *s;
  int orig_approx;
  phydbl orig_lnL;
  char *s_tree;

  if(tree->mcmc->run > mcmc->chain_len) return;

  s = (char *)mCalloc(100,sizeof(char));

  fp = mcmc->out_fp_stats;
  
/*   if(tree->mcmc->run == 0) */
/*     { */
/*       PhyML_Fprintf(stdout," ["); */
/*       fflush(NULL); */
/*     } */
  
/*   if(!(mcmc->run%(mcmc->chain_len/10))) */
/*     { */
/*       PhyML_Fprintf(stdout,"."); */
/*       fflush(NULL); */
/*     } */

/*   if(tree->mcmc->run == mcmc->chain_len) */
/*     { */
/*       PhyML_Fprintf(stdout,"]"); */
/*       fflush(NULL); */
/*     } */


/*   MCMC_Print_Means(mcmc,tree); */
/*   MCMC_Print_Last(mcmc,tree); */



  if(!(mcmc->run%mcmc->sample_interval)) 
    {

      MCMC_Copy_To_New_Param_Val(tree->mcmc,tree);

      For(i,tree->mcmc->n_moves)
	{
	  if((tree->mcmc->acc_rate[i] > .1) && (tree->mcmc->start_ess[i] == NO)) tree->mcmc->start_ess[i] = YES;
	  if(tree->mcmc->start_ess[i] == YES) MCMC_Update_Effective_Sample_Size(i,tree->mcmc,tree);
	  if(tree->mcmc->run > (int) tree->mcmc->chain_len * 0.1) tree->mcmc->adjust_tuning[i] = NO;
	}

      if(tree->mcmc->run == 0)
	{
	  time(&(mcmc->t_beg));
	  time(&(mcmc->t_last_print));

	  PhyML_Fprintf(fp,"# Random seed: %d",tree->io->r_seed);
	  PhyML_Fprintf(fp,"\n");
	  PhyML_Fprintf(fp,"Run\t");
/* 	  PhyML_Fprintf(fp,"Time\t"); */
	  /* PhyML_Fprintf(fp,"MeanRate\t"); */
/* 	  PhyML_Fprintf(fp,"NormFact\t"); */

	  /* For(i,mcmc->n_moves) */
	  /*   { */
	  /*     strcpy(s,"Acc."); */
	  /*     PhyML_Fprintf(fp,"%s%d\t",strcat(s,mcmc->move_name[i]),i); */
	  /*   } */

	  /* For(i,mcmc->n_moves) */
	  /*   { */
	  /*     strcpy(s,"Tune."); */
	  /*     PhyML_Fprintf(fp,"%s%d\t",strcat(s,mcmc->move_name[i]),i); */
	  /*   } */

	  /* For(i,mcmc->n_moves) */
	  /*   { */
	  /*     strcpy(s,"Run."); */
	  /*     PhyML_Fprintf(fp,"%s\t",strcat(s,mcmc->move_name[i])); */
	  /*   } */

	  PhyML_Fprintf(fp,"LnLike[Exact]\t");
	  PhyML_Fprintf(fp,"LnLike[Approx]\t");
	  PhyML_Fprintf(fp,"LnRates\t");
	  PhyML_Fprintf(fp,"LnTimes\t");
	  PhyML_Fprintf(fp,"LnPosterior\t");
	  PhyML_Fprintf(fp,"ClockRate\t");
	  PhyML_Fprintf(fp,"EvolRate\t");
	  PhyML_Fprintf(fp,"Nu\t");
	  PhyML_Fprintf(fp,"BirthRate\t");
	  PhyML_Fprintf(fp,"TsTv\t");


	  if(tree->mod->n_catg > 1)
	    {
	      if(tree->mod->free_mixt_rates == NO) PhyML_Fprintf(fp,"Alpha\t");
	      else
		{
		  For(i,tree->mod->n_catg) PhyML_Fprintf(fp,"p%d\t",i);
		  For(i,tree->mod->n_catg) PhyML_Fprintf(fp,"r%d\t",i);
		}
	    }


	  if(tree->mod->m4mod->n_h > 1 && tree->mod->use_m4mod == YES)
	    {
	      For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"cov_p%d\t",i);
	      For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"cov_r%d\t",i);
	      PhyML_Fprintf(fp,"cov_switch\t");
	    }


	  if(fp != stdout)
	    {
	      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	  	{
	  	  if(tree->t_nodes[i] == tree->n_root->v[0])
		    PhyML_Fprintf(fp,"T%d%s\t",i,"[LeftRoot]");
		  else if(tree->t_nodes[i] == tree->n_root->v[1])
		    PhyML_Fprintf(fp,"T%d%s\t",i,"[RightRoot]");
		  else if(tree->t_nodes[i] == tree->n_root)
		    PhyML_Fprintf(fp,"T%d%s\t",i,"[Root]");
		  else
		    PhyML_Fprintf(fp,"T%d[%d]\t",i,tree->t_nodes[i]->anc->num);
	  	}
	    }


/* 	  if(fp != stdout) */
/* 	    { */
/* 	      For(i,2*tree->n_otu-1) */
/* 	  	{ */
/* 	  	  if(tree->t_nodes[i] == tree->n_root->v[0]) */
/* 	  	    PhyML_Fprintf(fp,"R%d[LeftRoot]\t",i); */
/* 	  	  else if(tree->t_nodes[i] == tree->n_root->v[1]) */
/* 	  	    PhyML_Fprintf(fp,"R%d[RightRoot]\t",i); */
/* 	  	  else if(tree->t_nodes[i] != tree->n_root) */
/* 	  	    PhyML_Fprintf(fp," R%d[%d]\t",i,tree->t_nodes[i]->anc->num); */
/* 		  else */
/* 	  	    PhyML_Fprintf(fp," R%d[Root]\t",i); */

/* /\* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); *\/ */
/* 	  	} */
/* 	    } */


	  if(fp != stdout)
	    {
	      For(i,2*tree->n_otu-2)
	  	{
	  	  if(tree->t_nodes[i] == tree->n_root->v[0])
	  	    PhyML_Fprintf(fp,"B%d[LeftRoot]\t",i);
	  	  else if(tree->t_nodes[i] == tree->n_root->v[1])
	  	    PhyML_Fprintf(fp,"B%d[RightRoot]\t",i);
	  	  else
	  	    PhyML_Fprintf(fp," B%d[%d]\t",i,tree->t_nodes[i]->anc->num);

/* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); */
	  	}
	    }
	  

      /* 	  if(fp != stdout) */
      /* 	    { */
      /* 	      For(i,2*tree->n_otu-2) */
      /* 	  	{ */
      /* 	  	  if(tree->t_nodes[i] == tree->n_root->v[0]) */
      /* 	  	    PhyML_Fprintf(fp,"G%d[LeftRoot]\t",i); */
      /* 	  	  else if(tree->t_nodes[i] == tree->n_root->v[1]) */
      /* 	  	    PhyML_Fprintf(fp,"G%d[RightRoot]\t",i); */
      /* 	  	  else */
      /* 	  	    PhyML_Fprintf(fp," G%d[%f]\t",i,tree->rates->ml_l[i]); */


      /* /\* 		  PhyML_Fprintf(fp," R%d[%f]\t",i,tree->rates->mean_l[i]); *\/ */
      /* 	  	} */
      /* 	    } */


	  if(fp != stdout)
	    {
	      For(i,2*tree->n_otu-3)
		{
		  if(tree->t_edges[i] == tree->e_root)
		    PhyML_Fprintf(fp,"*L[%f]%d\t",i,tree->rates->u_ml_l[i]);
		  else
		    PhyML_Fprintf(fp," L[%f]%d\t",i,tree->rates->u_ml_l[i]);
		}
	    }


	  PhyML_Fprintf(mcmc->out_fp_trees,"#NEXUS\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"BEGIN TREES;\n");
	  PhyML_Fprintf(mcmc->out_fp_trees,"\tTRANSLATE\n");
	  For(i,tree->n_otu-1) PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s,\n",tree->t_nodes[i]->num,tree->t_nodes[i]->name);
	  PhyML_Fprintf(mcmc->out_fp_trees,"\t%3d\t%s;\n",tree->t_nodes[i]->num,tree->t_nodes[i]->name);
	  tree->write_tax_names = NO;
	}

      PhyML_Fprintf(fp,"\n");

      PhyML_Fprintf(fp,"%6d\t",tree->mcmc->run);

/*       time(&mcmc->t_cur); */
/*       PhyML_Fprintf(fp,"%6d\t",(int)(mcmc->t_cur-mcmc->t_beg)); */
      
/*       RATES_Update_Cur_Bl(tree); */
/*       PhyML_Fprintf(fp,"%f\t",RATES_Check_Mean_Rates(tree)); */

/*       PhyML_Fprintf(fp,"%f\t",tree->rates->norm_fact); */
      /* For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%f\t",tree->mcmc->acc_rate[i]); */
      /* For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%f\t",(phydbl)(tree->mcmc->tune_move[i])); */
/*       For(i,tree->mcmc->n_moves) PhyML_Fprintf(fp,"%d\t",(int)(tree->mcmc->run_move[i])); */

      orig_approx = tree->io->lk_approx;
      orig_lnL = tree->c_lnL;
      tree->io->lk_approx = EXACT;
      if(tree->mcmc->use_data)  Lk(tree);  else tree->c_lnL = 0.0;
      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);
      tree->io->lk_approx = NORMAL;
      tree->c_lnL = 0.0;
      if(tree->mcmc->use_data)  Lk(tree);  else tree->c_lnL = 0.0;
      PhyML_Fprintf(fp,"%.1f\t",tree->c_lnL);
      tree->io->lk_approx = orig_approx;
      tree->c_lnL = orig_lnL;

/*       PhyML_Fprintf(fp,"0\t0\t"); */


      PhyML_Fprintf(fp,"%G\t",tree->rates->c_lnL_rates);
      PhyML_Fprintf(fp,"%G\t",tree->rates->c_lnL_times);
      PhyML_Fprintf(fp,"%G\t",tree->c_lnL+tree->rates->c_lnL_rates+tree->rates->c_lnL_times);
      PhyML_Fprintf(fp,"%G\t",tree->rates->clock_r);
      PhyML_Fprintf(fp,"%G\t",RATES_Average_Substitution_Rate(tree));
      PhyML_Fprintf(fp,"%G\t",tree->rates->nu);
      PhyML_Fprintf(fp,"%G\t",tree->rates->birth_rate);
      PhyML_Fprintf(fp,"%G\t",tree->mod->kappa->v);
      
      if(tree->mod->n_catg > 1)
	{
	  if(tree->mod->free_mixt_rates == NO)
	    PhyML_Fprintf(fp,"%G\t",tree->mod->alpha);
	  else
	    {
	      For(i,tree->mod->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->gamma_r_proba->v[i]);
	      For(i,tree->mod->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->gamma_rr->v[i]);
	      /* For(i,tree->mod->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->gamma_r_proba_unscaled[i]); */
	      /* For(i,tree->mod->n_catg) PhyML_Fprintf(fp,"%G\t",tree->mod->gamma_rr_unscaled[i]); */
	    }
	}

      if(tree->mod->m4mod->n_h > 1 && tree->mod->use_m4mod == YES)
	{
	  For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->h_fq[i]);
	  For(i,tree->mod->m4mod->n_h) PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->multipl[i]);
	  PhyML_Fprintf(fp,"%G\t",tree->mod->m4mod->delta);
	}


      char *format = (char *)mCalloc(100,sizeof(char));
      sprintf(format,"%%.%df\t",tree->mcmc->nd_t_digits);
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,format,tree->rates->nd_t[i]);
      Free(format);

      /* for(i=0;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%.4f\t",LOG(tree->rates->nd_r[i])); */

      for(i=0;i<2*tree->n_otu-2;i++) 
	PhyML_Fprintf(fp,"%.4f\t",
		      (tree->rates->model_log_rates == YES)?
		      EXP(tree->rates->br_r[i]):
		      tree->rates->br_r[i]);

      /* fp_pred = fopen("predict.txt","a");       */
      /* for(i=0;i<2*tree->n_otu-2;i++)  */
      /* 	PhyML_Fprintf(fp_pred,"B%d\t%12f\t%12f\t%4d\n",i,EXP(tree->rates->br_r[i]),tree->rates->nd_t[i],tree->rates->has_survived[i]); */
      /* fclose(fp_pred); */


      /* phydbl p,sd,mean; */
      
      /* for(i=0;i<2*tree->n_otu-2;i++) */
      /* 	{ */
      /* 	  sd = tree->rates->nu * (tree->rates->nd_t[i] - tree->rates->nd_t[tree->t_nodes[i]->anc->num]); */
      /* 	  mean = tree->rates->br_r[tree->t_nodes[i]->anc->num] - .5*sd*sd; */
      /* 	  p = Pnorm(tree->rates->br_r[i],mean,sd); */
      /* 	  PhyML_Fprintf(fp,"%f\t",p); */
      /* 	  tree->rates->mean_r[i] += p; */
      /* 	} */

      /* for(i=0;i<2*tree->n_otu-2;i++) PhyML_Fprintf(fp,"%.4f\t",tree->rates->cur_gamma_prior_mean[i]); */
      /* if(fp != stdout) for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(fp,"%G\t",tree->rates->t_prior[i]); */
/*       For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",EXP(tree->t_edges[i]->l)); */


      RATES_Update_Cur_Bl(tree);
      For(i,2*tree->n_otu-3) PhyML_Fprintf(fp,"%f\t",tree->t_edges[i]->l);


      For(i,2*tree->n_otu-2) tree->rates->mean_r[i] = EXP(tree->rates->br_r[i]);
      For(i,2*tree->n_otu-1) tree->rates->mean_t[i] = tree->rates->nd_t[i];
      
      /* Branch_Lengths_To_Time_Lengths(tree); */
      /* char *s; */
      /* s = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
      /* strcpy(s,mcmc->io->in_align_file); */
      /* strcat(s,"_"); */
      /* strcat(s,mcmc->out_filename); */
      /* strcat(s,".ps"); */
      /* DR_Draw_Tree(s,tree); */
      /* Free(s);	 */

      /* FILE *fp; */
      /* int j; */
      /* t_node *d, *v1, *v2; */
      /* int n1, n2; */
      /* phydbl r1, r2; */

      /* s = (char *)mCalloc(T_MAX_NAME,sizeof(char)); */
      /* strcpy(s,mcmc->io->in_align_file); */
      /* strcat(s,"_"); */
      /* strcat(s,mcmc->out_filename); */
      /* strcat(s,".corr"); */
      /* fp = fopen(s,"w"); */

      /* n1 = n2 = 0; */
      /* r1 = r2 = 0.f; */
      /* For(i,2*tree->n_otu-3) */
      /* 	{ */
      /* 	  if(tree->t_nodes[i]->tax == NO) */
      /* 	    { */
      /* 	      d = tree->t_nodes[i]; */
      /* 	      v1 = v2 = NULL; */
      /* 	      For(j,3) */
      /* 		{ */
      /* 		  if(d->v[j] != d->anc && d->b[j] != tree->e_root) */
      /* 		    { */
      /* 		      if(!v1) v1 = d->v[j]; */
      /* 		      else    v2 = d->v[j]; */
      /* 		    } */
      /* 		} */

      /* 	      For(j,3) */
      /* 		{ */
      /* 		  if(v1->v[j] && v1->v[j] == d) */
      /* 		    { */
      /* 		      n1 = v1->bip_size[j]; */
      /* 		      /\* r1 = RATES_Get_Mean_Rate_In_Subtree(v1,tree); *\/ */
      /* 		      r1 = EXP(tree->rates->br_r[v1->num]); */
      /* 		      break; */
      /* 		    } */
      /* 		} */


      /* 	      For(j,3) */
      /* 		{ */
      /* 		  if(v2->v[j] && v2->v[j] == d) */
      /* 		    { */
      /* 		      n2 = v2->bip_size[j]; */
      /* 		      /\* r2 = RATES_Get_Mean_Rate_In_Subtree(v2,tree); *\/ */
      /* 		      r2 = EXP(tree->rates->br_r[v2->num]); */
      /* 		      break; */
      /* 		    } */
      /* 		} */

      /* 	      fprintf(fp,"\n%4d %4d %15f %15f",n1,n2,r1,r2); */
      /* 	    } */
      /* 	} */

      /* fclose(fp); */
      
      // TREES
      Branch_Lengths_To_Time_Lengths(tree);
      tree->bl_ndigits = 3;
      s_tree = Write_Tree(tree,NO);
      tree->bl_ndigits = 7;
      PhyML_Fprintf(mcmc->out_fp_trees,"TREE %8d [%f] = [&R] %s\n",mcmc->run,tree->c_lnL,s_tree);
      Free(s_tree);
    }

  if(tree->mcmc->run == mcmc->chain_len) PhyML_Fprintf(mcmc->out_fp_trees,"END;\n");

  fflush(NULL);   

  Free(s);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Means(t_mcmc *mcmc, t_tree *tree)
{

  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      int i;
      char *s;

      s = (char *)mCalloc(T_MAX_FILE,sizeof(char));

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".means");
      
      fclose(mcmc->out_fp_means);

      mcmc->out_fp_means = fopen(s,"w");
      
      PhyML_Fprintf(mcmc->out_fp_means,"#");
      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(mcmc->out_fp_means,"T%d\t",i);	  

      PhyML_Fprintf(mcmc->out_fp_means,"\n");      

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) tree->rates->t_mean[i] *= (phydbl)(mcmc->run / mcmc->sample_interval);

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	{
	  tree->rates->t_mean[i] += tree->rates->nd_t[i];
	  tree->rates->t_mean[i] /= (phydbl)(mcmc->run / mcmc->sample_interval + 1);

/* 	  PhyML_Fprintf(tree->mcmc->out_fp_means,"%d\t",tree->mcmc->run / tree->mcmc->sample_interval);	   */
	  PhyML_Fprintf(tree->mcmc->out_fp_means,"%.1f\t",tree->rates->t_mean[i]);
	}

      PhyML_Fprintf(tree->mcmc->out_fp_means,"\n");
      fflush(NULL);

      Free(s);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Print_Last(t_mcmc *mcmc, t_tree *tree)
{

  if(!(mcmc->run%mcmc->sample_interval)) 
    {
      int i;
      char *s;

      s = (char *)mCalloc(T_MAX_FILE,sizeof(char));

      strcpy(s,tree->mcmc->out_filename);
      strcat(s,".lasts");
      
      fclose(mcmc->out_fp_last);

      mcmc->out_fp_last = fopen(s,"w");

/*       rewind(mcmc->out_fp_last); */

      PhyML_Fprintf(mcmc->out_fp_last,"#");
      PhyML_Fprintf(tree->mcmc->out_fp_last,"Time\t");

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++)
	PhyML_Fprintf(tree->mcmc->out_fp_last,"T%d\t",i);

      PhyML_Fprintf(tree->mcmc->out_fp_last,"\n");

      if(mcmc->run)
	{
	  time(&(mcmc->t_cur));
	  PhyML_Fprintf(tree->mcmc->out_fp_last,"%d\t",(int)(mcmc->t_cur-mcmc->t_beg));
/* 	  PhyML_Fprintf(tree->mcmc->out_fp_last,"%d\t",(int)(mcmc->t_beg)); */
	}

      for(i=tree->n_otu;i<2*tree->n_otu-1;i++) PhyML_Fprintf(tree->mcmc->out_fp_last,"%.1f\t",tree->rates->nd_t[i]);

      PhyML_Fprintf(tree->mcmc->out_fp_last,"\n");
      fflush(NULL);

      Free(s);
  }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

t_mcmc *MCMC_Make_MCMC_Struct()
{
  t_mcmc *mcmc;

  mcmc               = (t_mcmc *)mCalloc(1,sizeof(t_mcmc));
  mcmc->out_filename = (char *)mCalloc(T_MAX_FILE,sizeof(char));

  return(mcmc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Free_MCMC(t_mcmc *mcmc)
{
  int i;
  
  Free(mcmc->move_type);
  Free(mcmc->adjust_tuning);
  Free(mcmc->out_filename);
  Free(mcmc->move_weight);
  Free(mcmc->acc_move);
  Free(mcmc->run_move);
  Free(mcmc->prev_acc_move);
  Free(mcmc->prev_run_move);
  Free(mcmc->acc_rate);
  Free(mcmc->tune_move);
  For(i,mcmc->n_moves) Free(mcmc->move_name[i]);
  Free(mcmc->move_name);
  Free(mcmc->ess_run);
  Free(mcmc->start_ess);
  Free(mcmc->ess);
  Free(mcmc->sum_val);
  Free(mcmc->sum_valsq);
  Free(mcmc->sum_curval_nextval);
  Free(mcmc->first_val);
  Free(mcmc->old_param_val);
  Free(mcmc->new_param_val);
  Free(mcmc);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Pause(t_mcmc *mcmc)
{
  char choice;
  char *s;
  int len;

  s = (char *)mCalloc(100,sizeof(char));

  
  if(!(mcmc->run%mcmc->chain_len) && (mcmc->is_burnin == NO))
    {
      PhyML_Printf("\n. Do you wish to stop the analysis [N/y] ");
      if(!scanf("%c",&choice)) Exit("\n");
      if(choice == '\n') choice = 'N';
      else getchar(); /* \n */
      
      Uppercase(&choice);
	  
      switch(choice)
	{
	case 'N': 
	  {	    
	    len = 1E+4;
	    PhyML_Printf("\n. How many extra generations is required [default: 1E+4] ");
	    Getstring_Stdin(s);
	    if(s[0] == '\0') len = 1E+4; 
	    else len = (int)atof(s); 

	    if(len < 0)
	      {
		PhyML_Printf("\n. The value entered must be an integer greater than 0.\n");
		Exit("\n");
	      }	    
	    mcmc->chain_len += len;
	    break;
	  }
	      
	case 'Y': 
	  {
	    PhyML_Printf("\n. Ok. Done.\n");
	    Exit("\n");
	    break;
	  }
	  
	default: 
	  {
	    PhyML_Printf("\n. Please enter 'Y' or 'N'.\n");
	    Exit("\n");
	  }
	}
    }

  Free(s);

}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Terminate()
{
  char choice;
  PhyML_Printf("\n\n. Do you really want to terminate [Y/n]: ");
  if(!scanf("%c",&choice)) Exit("\n");
  if(choice == '\n') choice = 'Y';
  else getchar(); /* \n */      
  Uppercase(&choice);	  
  if(choice == 'Y') raise(SIGTERM);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Init_MCMC_Struct(char *filename, option *io, t_mcmc *mcmc)
{
  int pid;

  mcmc->io               = io;
  mcmc->is               = NO;
  mcmc->use_data         = YES;
  mcmc->run              = 0;
  mcmc->sample_interval  = 1E+3;
  mcmc->chain_len        = 1E+6;
  mcmc->chain_len_burnin = 1E+5;
  mcmc->randomize        = YES;
  mcmc->norm_freq        = 1E+3;
  mcmc->max_tune         = 1.E+20;
  mcmc->min_tune         = 1.E-10;
  mcmc->print_every      = 2;
  mcmc->is_burnin        = NO;
  mcmc->nd_t_digits      = 1;

  if(filename)
    {
      char *s;

      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      strcpy(mcmc->out_filename,filename);
      pid = getpid();
      sprintf(mcmc->out_filename+strlen(mcmc->out_filename),".%d",pid);

      strcpy(s,mcmc->io->in_align_file);
      strcat(s,"_");
      strcat(s,mcmc->out_filename);
      mcmc->out_fp_stats = fopen(s,"w");

      strcpy(s,mcmc->io->in_align_file);
      strcat(s,"_");
      strcat(s,mcmc->out_filename);
      strcat(s,".trees");
      mcmc->out_fp_trees = fopen(s,"w");

      strcpy(s,mcmc->io->in_align_file);
      strcat(s,"_");
      strcat(s,mcmc->out_filename);
      strcat(s,".constree");
      mcmc->out_fp_constree = fopen(s,"w");

/*       strcpy(s,tree->mcmc->out_filename); */
/*       strcat(s,".means"); */
/*       tree->mcmc->out_fp_means = fopen(s,"w"); */

/*       strcpy(s,tree->mcmc->out_filename); */
/*       strcat(s,".lasts"); */
/*       tree->mcmc->out_fp_last  = fopen(s,"w"); */

      Free(s);
    }
  else 
    {
      mcmc->out_fp_stats = stderr;
      mcmc->out_fp_trees = stderr;
      /* tree->mcmc->out_fp_means = stderr; */
      /* tree->mcmc->out_fp_last  = stderr; */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Copy_MCMC_Struct(t_mcmc *ori, t_mcmc *cpy, char *filename)
{
  int pid;
  int i;
  
  cpy->use_data           = ori->use_data        ;
  cpy->sample_interval    = ori->sample_interval ;
  cpy->chain_len          = ori->chain_len       ;
  cpy->randomize          = ori->randomize       ;
  cpy->norm_freq          = ori->norm_freq       ;
  cpy->n_moves            = ori->n_moves         ;
  cpy->max_tune           = ori->max_tune        ;
  cpy->min_tune           = ori->min_tune        ;
  cpy->print_every        = ori->print_every     ;
  cpy->is_burnin          = ori->is_burnin       ;
  cpy->is                 = ori->is              ;
  cpy->io                 = ori->io              ;
  cpy->in_fp_par          = ori->in_fp_par       ;
  cpy->nd_t_digits        = ori->nd_t_digits     ;

  For(i,cpy->n_moves) 
    {
      cpy->old_param_val[i]      = ori->old_param_val[i];
      cpy->new_param_val[i]      = ori->new_param_val[i];
      cpy->start_ess[i]          = ori->start_ess[i];
      cpy->ess_run[i]            = ori->ess_run[i];
      cpy->ess[i]                = ori->ess[i];
      cpy->sum_val[i]            = ori->sum_val[i];
      cpy->sum_valsq[i]          = ori->sum_valsq[i];
      cpy->first_val[i]          = ori->first_val[i];
      cpy->sum_curval_nextval[i] = ori->sum_curval_nextval[i];
      cpy->move_weight[i]        = ori->move_weight[i];
      cpy->run_move[i]           = ori->run_move[i];
      cpy->acc_move[i]           = ori->acc_move[i];
      cpy->prev_run_move[i]      = ori->prev_run_move[i];
      cpy->prev_acc_move[i]      = ori->prev_acc_move[i];
      cpy->acc_rate[i]           = ori->acc_rate[i];
      cpy->tune_move[i]          = ori->tune_move[i];
      strcpy(cpy->move_name[i],ori->move_name[i]);
      cpy->adjust_tuning[i]      = ori->adjust_tuning[i];
    }


  
  if(filename) 
    {
      char *s;

      s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

      strcpy(cpy->out_filename,filename);
      pid = getpid();
      sprintf(cpy->out_filename+strlen(cpy->out_filename),".%d",pid);

      strcpy(s,cpy->io->in_align_file);
      strcat(s,"_");
      strcat(s,cpy->out_filename);
      cpy->out_fp_stats = fopen(s,"w");

      strcpy(s,cpy->io->in_align_file);
      strcat(s,"_");
      strcpy(s,cpy->out_filename);
      strcat(s,".trees");
      cpy->out_fp_trees = fopen(s,"w");

      strcpy(s,cpy->io->in_align_file);
      strcat(s,"_");
      strcpy(s,cpy->out_filename);
      strcat(s,".constree");
      cpy->out_fp_constree = fopen(s,"w");

      Free(s);
    }
  else 
    {
      cpy->out_fp_stats = stderr;
      cpy->out_fp_trees = stderr;
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Close_MCMC(t_mcmc *mcmc)
{
  fclose(mcmc->out_fp_trees);
  fclose(mcmc->out_fp_stats);
  fclose(mcmc->out_fp_constree);
  /* fclose(mcmc->out_fp_means); */
  /* fclose(mcmc->out_fp_last); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Kappa(t_tree *tree)
{
  tree->mod->kappa->v = Uni()*5.;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Rate_Across_Sites(t_tree *tree)
{
  if(tree->mod->n_catg == 1) return;
  
  if(tree->mod->free_mixt_rates == YES)
    {
      int i;
      For(i,tree->mod->n_catg)
	{
	  tree->mod->gamma_r_proba_unscaled->v[i] = Uni()*100.;
	  tree->mod->gamma_rr_unscaled->v[i] = (phydbl)i+1.; /* Do not randomize those as their ordering matter */
	}
    }
  else
    {
      tree->mod->alpha->v = Uni()*5.;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Covarion_Rates(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;
 
  if(tree->mod->m4mod->n_h == 1) return;
  
  int i;

  For(i,tree->mod->m4mod->n_h)
    {
      tree->mod->m4mod->multipl_unscaled[i] = (phydbl)i+1.;
      tree->mod->m4mod->h_fq_unscaled[i] = Uni()*(100.-0.01) + 0.01;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Covarion_Switch(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;

  tree->mod->m4mod->delta = Uni()*(10.-0.01)+0.01;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Branch_Lengths(t_tree *tree)
{
  int i;

  if(tree->mod->log_l == NO)
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l = Rexp(10.);
  else
    For(i,2*tree->n_otu-3) tree->t_edges[i]->l = -4* Uni();
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Rates(t_tree *tree)
{

  int i,err;
  phydbl mean_r, var_r;
  phydbl min_r, max_r;

  mean_r = 1.0;
  var_r  = 0.5;
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;

  For(i,2*tree->n_otu-2)
    if(tree->t_nodes[i] != tree->n_root)
      tree->rates->nd_r[i] = Rnorm_Trunc(mean_r,SQRT(var_r),min_r,max_r,&err);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Rates(t_tree *tree)
{

  /* Should be called once t_node times have been determined */

  int i;
  phydbl u;
  phydbl r_min,r_max;

  For(i,2*tree->n_otu-2) tree->rates->br_r[i] = 1.0;

  r_min = 0.9;
  r_max = 1.1;
  u     = 0.0;

/*   For(i,2*tree->n_otu-2) */
/*     { */
/*       u = Uni(); */
/*       u = u * (r_max-r_min) + r_min; */
/*       tree->rates->br_r[i] = u; */

/*       if(tree->rates->br_r[i] < tree->rates->min_rate) tree->rates->br_r[i] = tree->rates->min_rate;  */
/*       if(tree->rates->br_r[i] > tree->rates->max_rate) tree->rates->br_r[i] = tree->rates->max_rate;  */
/*     } */

  MCMC_Randomize_Rates_Pre(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Rates_Pre(tree->n_root,tree->n_root->v[1],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Rates_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  phydbl mean_r, var_r;
  phydbl min_r, max_r;
  int err;

/*   mean_r = tree->rates->br_r[a->num]; */
/*   var_r  = tree->rates->nu * (tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]); */
  mean_r = 1.0;
  var_r  = 0.5;
  min_r  = tree->rates->min_rate;
  max_r  = tree->rates->max_rate;
  
  tree->rates->br_r[d->num] = Rnorm_Trunc(mean_r,SQRT(var_r),min_r,max_r,&err);

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Randomize_Rates_Pre(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Randomize_Birth(t_tree *tree)
{
  phydbl min_b,max_b;
  phydbl u;

  min_b = tree->rates->birth_rate_min;
  max_b = tree->rates->birth_rate_max;
  
  u = Uni();
  tree->rates->birth_rate = (max_b - min_b) * u + min_b;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Randomize_Nu(t_tree *tree)
{
  phydbl min_nu,max_nu;
  phydbl u;

  /* It is preferable to start with small values of nu 
     as if is difficult for the MCMC sampler to sample
     equal rates on edge (i.e., molecular clock) since
     such combination of rate lies on the boundary of 
     the space of all edge rate combination. We give
     here a bit of help to the sampler by considering 
     starting points close to the molecular clock
     constraint.
  */
  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu/10.;
  
  u = Uni();
  tree->rates->nu = (max_nu - min_nu) * u + min_nu;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Clock_Rate(t_tree *tree)
{
  phydbl u;
  u = Uni();
  tree->rates->clock_r = u * (tree->rates->max_clock - tree->rates->min_clock) + tree->rates->min_clock;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Alpha(t_tree *tree)
{
  phydbl u;

  u = Uni();
  tree->rates->alpha = u*6.0+1.0;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Times(t_tree *tree)
{
  phydbl t_sup, t_inf;
  phydbl u;
  int iter;
  int i;
  phydbl dt,min_dt;
  int min_node;

  t_inf = tree->rates->t_prior_min[tree->n_root->num];
  t_sup = tree->rates->t_prior_max[tree->n_root->num];

  u = Uni();
  u *= (t_sup - t_inf);
  u += t_inf;
  
  tree->rates->nd_t[tree->n_root->num] = u;

  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Randomize_Node_Times_Top_Down(tree->n_root,tree->n_root->v[1],tree);
  
  min_node = -1;
  iter = 0;
  do
    {
      min_dt = MDBL_MAX;
      For(i,2*tree->n_otu-2) 
	{
	  dt = tree->rates->nd_t[i] - tree->rates->nd_t[tree->t_nodes[i]->anc->num];
	  if(dt < min_dt)
	    {
	      min_dt = dt;
	      min_node = i;
	    }
	}

      if(min_dt > .1 * FABS(tree->rates->nd_t[tree->n_root->num])/(phydbl)(tree->n_otu-1)) break;

      RATES_Record_Times(tree);
      For(i,2*tree->n_otu-1) 
	{
	  if(tree->t_nodes[i]->tax == NO) 
	    tree->rates->nd_t[i] -= 0.1*FABS(tree->rates->nd_t[tree->n_root->num])/(phydbl)(tree->n_otu-1);

	  if(tree->rates->nd_t[i] < tree->rates->t_prior_min[i] || 
	     tree->rates->nd_t[i] > tree->rates->t_prior_max[i])
	    {
	      RATES_Reset_Times(tree);
	      break;
	    }
	}

      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[0],tree);
      MCMC_Randomize_Node_Times_Bottom_Up(tree->n_root,tree->n_root->v[1],tree);
      
      iter++;
    }
  while(iter < 1000);

  if(iter == 1000)
    {      
      PhyML_Printf("\n. min_dt = %f",min_dt);
      PhyML_Printf("\n. min->t=%f min->anc->t=%f",tree->rates->nd_t[min_node],tree->rates->nd_t[tree->t_nodes[min_node]->anc->num]);
      PhyML_Printf("\n. up=%f down=%f",tree->rates->t_prior_min[min_node],tree->rates->t_floor[tree->t_nodes[min_node]->anc->num]);
      PhyML_Printf("\n. min_node = %d",min_node);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }


/*   PhyML_Printf("\n. Needed %d iterations to randomize node heights.",iter); */
/*   TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[0],tree); */
/*   TIMES_Print_Node_Times(tree->n_root,tree->n_root->v[1],tree); */

  if(RATES_Check_Node_Times(tree))
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Times_Bottom_Up(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      phydbl u;
      phydbl t_inf, t_sup;
      t_node *v1, *v2;


      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Bottom_Up(d,d->v[i],tree);
	    }
	}

      v1 = v2 = NULL;
      For(i,3)
      	{
      	  if(d->v[i] != a && d->b[i] != tree->e_root)
      	    {
      	      if(!v1) v1 = d->v[i];
      	      else    v2 = d->v[i];
      	    }
      	}

      t_sup = MIN(tree->rates->nd_t[v1->num],tree->rates->nd_t[v2->num]);
      t_inf = tree->rates->nd_t[a->num];

      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
 
      if(u > tree->rates->t_prior_min[d->num] && u < tree->rates->t_prior_max[d->num])
	tree->rates->nd_t[d->num] = u;
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Randomize_Node_Times_Top_Down(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;      
      phydbl u;
      phydbl t_inf, t_sup;

      t_inf = MAX(tree->rates->nd_t[a->num],tree->rates->t_prior_min[d->num]);
      t_sup = tree->rates->t_prior_max[d->num];

      u = Uni();
      u *= (t_sup - t_inf);
      u += t_inf;
      
      tree->rates->nd_t[d->num] = u;

      For(i,3)
	{
	  if((d->v[i] != a) && (d->b[i] != tree->e_root))
	    {
	      MCMC_Randomize_Node_Times_Top_Down(d,d->v[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Get_Acc_Rates(t_mcmc *mcmc)
{
  int i;
  phydbl eps;
  int lag;


  if(mcmc->run < (int)(0.01*mcmc->chain_len)) lag = 100;
  else lag = 100;

  eps = 1.E-6;

  For(i,mcmc->n_moves)
    {
      if(mcmc->run_move[i] - mcmc->prev_run_move[i] > lag)
	{
	  mcmc->acc_rate[i] = 
	    (phydbl)(mcmc->acc_move[i] - mcmc->prev_acc_move[i] + eps) / 
	    (phydbl)(mcmc->run_move[i] - mcmc->prev_run_move[i] + eps) ;


	  mcmc->prev_run_move[i] = mcmc->run_move[i];
	  mcmc->prev_acc_move[i] = mcmc->acc_move[i];
	  
	  MCMC_Adjust_Tuning_Parameter(i,mcmc);
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Adjust_Tuning_Parameter(int move, t_mcmc *mcmc)
{
  if(mcmc->adjust_tuning[move])
    {
      phydbl scale;
      phydbl rate;
      phydbl rate_inf,rate_sup;
      
      if(mcmc->run < (int)(0.01*mcmc->chain_len)) scale = 1.5;
      else scale = 1.2;

      if(!strcmp(mcmc->move_name[move],"tree_height"))
	{
	  rate_inf = 0.3;
	  rate_sup = 0.3;
	}
      else if(!strcmp(mcmc->move_name[move],"subtree_height"))
	{
	  rate_inf = 0.2;
	  rate_sup = 0.2;
	}
      /* if(!strcmp(mcmc->move_name[move],"tree_rates")) */
      /* 	{ */
      /* 	  rate_inf = 0.05; */
      /* 	  rate_sup = 0.05; */
      /* 	} */
      else
	{
	  rate_inf = 0.3;
	  rate_sup = 0.3;
	}

      /* PhyML_Printf("\n. %s acc=%d run=%d tune=%f", */
      /* 		   mcmc->move_name[move], */
      /* 		   mcmc->acc_move[move], */
      /* 		   mcmc->run_move[move], */
      /* 		   mcmc->tune_move[move]); */

      rate = mcmc->acc_rate[move];
      
      if(rate < rate_inf)      
	{
	  mcmc->tune_move[move] /= scale;
	}

      else if(rate > rate_sup) 
	{
	  mcmc->tune_move[move] *= scale;
	}

      if(mcmc->tune_move[move] > mcmc->max_tune) mcmc->tune_move[move] = mcmc->max_tune;
      if(mcmc->tune_move[move] < mcmc->min_tune) mcmc->tune_move[move] = mcmc->min_tune;
	  
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_One_Length(t_edge *b, t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data;
  phydbl ratio, alpha;
  phydbl new_l, cur_l;
  phydbl K,mult;


  cur_l        = b->l;
  cur_lnL_data = tree->c_lnL;
  K            = 0.1;
  
  u = Uni();
  mult = EXP(K*(u-0.5));
  /* mult = u*(K-1./K)+1./K; */
  new_l = cur_l * mult;

  if(new_l < tree->mod->l_min || new_l > tree->mod->l_max) return;

  b->l = new_l;
  new_lnL_data = Lk_At_Given_Edge(b,tree);
/*   tree->both_sides = NO; */
/*   new_lnL_data = Lk(tree); */


  ratio =
    (new_lnL_data - cur_lnL_data) +
    (LOG(mult));


  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      b->l = cur_l;
      Update_PMat_At_Given_Edge(b,tree);
      tree->c_lnL = cur_lnL_data;
    }

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Scale_Br_Lens(t_tree *tree)
{
  phydbl u;
  phydbl new_lnL_data, cur_lnL_data;
  phydbl ratio, alpha;
  phydbl K,mult;
  int i;

  Record_Br_Len(tree);

  cur_lnL_data = tree->c_lnL;
  K            = 1.2;
  
  u = Uni();
  mult = u*(K-1./K)+1./K;

  For(i,2*tree->n_otu-3) 
    {
      tree->t_edges[i]->l *= mult;
      if(tree->t_edges[i]->l < tree->mod->l_min || 
	 tree->t_edges[i]->l > tree->mod->l_max) return;
    }

  tree->both_sides = NO;
  new_lnL_data = Lk(tree);

  ratio =
    (new_lnL_data - cur_lnL_data) +
    (2*tree->n_otu-5) * (LOG(mult));

  ratio = EXP(ratio);
  alpha = MIN(1.,ratio);
  
  u = Uni();
  
  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(tree);
      tree->c_lnL = cur_lnL_data;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Br_Lens(t_tree *tree)
{
  MCMC_Br_Lens_Pre(tree->t_nodes[0],
  		   tree->t_nodes[0]->v[0],
  		   tree->t_nodes[0]->b[0],tree);

  /* int i; */
  /* For(i,2*tree->n_otu-3) */
  /*   { */
  /*     MCMC_One_Length(tree->t_edges[Rand_Int(0,2*tree->n_otu-4)],acc,run,tree); */
  /*   } */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Br_Lens_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;


  if(a == tree->n_root || d == tree->n_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  MCMC_One_Length(b,tree);
  if(d->tax) return;
  else 
    {
      For(i,3) 
	if(d->v[i] != a)
	  {
	    Update_P_Lk(tree,d->b[i],d);
	    MCMC_Br_Lens_Pre(d,d->v[i],d->b[i],tree);
	  }
      Update_P_Lk(tree,b,d);
    }  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Kappa(t_tree *tree)
{
  int change,i;

  change = NO;

  tree->mod->update_eigen = YES;
		
  if(tree->io->lk_approx == NORMAL)
    {
      tree->io->lk_approx = EXACT;
      if(tree->mcmc->use_data == YES) Lk(tree);
      change = YES;
    }
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = NO;
  MCMC_Single_Param_Generic(&(tree->mod->kappa->v),0.,15.,tree->mcmc->num_move_kappa,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_kappa],NO,NULL,tree,NULL);
	  
  if(change == YES)
    {
      tree->io->lk_approx = NORMAL;
      if(tree->mcmc->use_data == YES) Lk(tree);
    }
  
  tree->mod->update_eigen = NO;  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Rate_Across_Sites(t_tree *tree)
{
  if(tree->mod->n_catg == 1) return;

  if(tree->mod->free_mixt_rates == YES)
    {
      MCMC_Free_Mixt_Rate(tree);
    }
  else
    {
      MCMC_Alpha(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Alpha(t_tree *tree)
{
  int i;
  
  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = NO;
  MCMC_Single_Param_Generic(&(tree->mod->alpha->v),0.,100.,tree->mcmc->num_move_ras,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_ras],NO,NULL,tree,NULL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Free_Mixt_Rate(t_tree *tree)
{
  int i, class;
  phydbl u;
  phydbl min,max;

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = NO;

  class = Rand_Int(0,tree->mod->n_catg-1);

  min = 0.01;
  max = +100.;
  u = Uni();
  if(u < .5)
    {
      if(!class)
      	{
      	  min = 0.01;
      	  max = tree->mod->gamma_rr_unscaled->v[1];
      	}
      else if(class == tree->mod->n_catg-1)
      	{
      	  min = tree->mod->gamma_rr_unscaled->v[tree->mod->n_catg-2];
      	  max = +100.;
      	}
      else
      	{
      	  min = MIN(tree->mod->gamma_rr_unscaled->v[class-1],tree->mod->gamma_rr_unscaled->v[class+1]);
      	  max = MAX(tree->mod->gamma_rr_unscaled->v[class-1],tree->mod->gamma_rr_unscaled->v[class+1]);
      	}
            
      /* int i; */
      /* printf("\n %3d (%15f %15f)",class,min,max); */
      /* For(i,tree->mod->n_catg) */
      /* 	printf("%15f ",tree->mod->gamma_rr_unscaled[i]); */

      MCMC_Single_Param_Generic(&(tree->mod->gamma_rr_unscaled->v[class]),min,max,tree->mcmc->num_move_ras+class+tree->mod->n_catg,
				NULL,&(tree->c_lnL),
				NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_ras+class+tree->mod->n_catg],NO,NULL,tree,NULL);

    }
  else
    {
      MCMC_Single_Param_Generic(&(tree->mod->gamma_r_proba_unscaled->v[class]),0.01,+100.,tree->mcmc->num_move_ras+class,
				NULL,&(tree->c_lnL),
				NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_ras+class],NO,NULL,tree,NULL);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Covarion_Rates(t_tree *tree)
{
  int i, class;
  phydbl u;
  phydbl min,max;

  if(tree->mod->use_m4mod == NO) return;

  tree->mod->update_eigen = YES;

  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;

  class = Rand_Int(0,tree->mod->m4mod->n_h-1);


  min = 0.01;
  max = 100.;
  u = Uni();
  if(u < .5)
    {
      if(!class)
	{
	  min = 0.01;
	  max = tree->mod->m4mod->multipl_unscaled[1];
	}
      else if(class == tree->mod->m4mod->n_h-1)
	{
	  min = tree->mod->m4mod->multipl_unscaled[tree->mod->m4mod->n_h-2];
	  max = +100.;
	}
      else
	{
	  min = MIN(tree->mod->m4mod->multipl_unscaled[class-1],tree->mod->m4mod->multipl_unscaled[class+1]);
	  max = MAX(tree->mod->m4mod->multipl_unscaled[class-1],tree->mod->m4mod->multipl_unscaled[class+1]);
	}

      MCMC_Single_Param_Generic(&(tree->mod->m4mod->multipl_unscaled[class]),min,max,tree->mcmc->num_move_cov_rates+class+tree->mod->m4mod->n_h,
				NULL,&(tree->c_lnL),
				NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_rates+class+tree->mod->m4mod->n_h],NO,NULL,tree,NULL);
    }
  else
    {
      MCMC_Single_Param_Generic(&(tree->mod->m4mod->h_fq_unscaled[class]),0.01,+100.,tree->mcmc->num_move_cov_rates+class,
				NULL,&(tree->c_lnL),
				NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_rates+class],NO,NULL,tree,NULL);
    }

  tree->mod->update_eigen = NO;

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Covarion_Switch(t_tree *tree)
{
  if(tree->mod->use_m4mod == NO) return;

  tree->mod->update_eigen = YES;
  MCMC_Single_Param_Generic(&(tree->mod->m4mod->delta),0.01,+100.,tree->mcmc->num_move_cov_switch,
			    NULL,&(tree->c_lnL),
			    NULL,Wrap_Lk,tree->mcmc->move_type[tree->mcmc->num_move_cov_switch],NO,NULL,tree,NULL); 
  tree->mod->update_eigen = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Birth_Rate(t_tree *tree)
{
  MCMC_Single_Param_Generic(&(tree->rates->birth_rate),
			    tree->rates->birth_rate_min,
			    tree->rates->birth_rate_max,
			    tree->mcmc->num_move_birth_rate,
			    &(tree->rates->c_lnL_times),NULL,
			    Wrap_Lk_Times,NULL,tree->mcmc->move_type[tree->mcmc->num_move_birth_rate],NO,NULL,tree,NULL); 
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Nu(t_tree *tree)
{
  phydbl cur_nu,new_nu,cur_lnL_rate,new_lnL_rate;
  phydbl u,alpha,ratio;
  phydbl min_nu,max_nu;
  phydbl K;
  phydbl cur_lnL_data, new_lnL_data;
  int i;

  Record_Br_Len(tree);

  cur_nu        = -1.0;
  new_nu        = -1.0;
  ratio         =  0.0;

  K = tree->mcmc->tune_move[tree->mcmc->num_move_nu];

  cur_lnL_rate = tree->rates->c_lnL_rates;
  new_lnL_rate = tree->rates->c_lnL_rates;

  cur_lnL_data = tree->c_lnL;
  new_lnL_data = tree->c_lnL;
  
  cur_nu       = tree->rates->nu;

  min_nu = tree->rates->min_nu;
  max_nu = tree->rates->max_nu;

  MCMC_Make_Move(&cur_nu,&new_nu,min_nu,max_nu,&ratio,K,tree->mcmc->move_type[tree->mcmc->num_move_nu]);
  
  if(new_nu < max_nu && new_nu > min_nu) 
    {
      tree->rates->nu = new_nu;
      
      new_lnL_rate = RATES_Lk_Rates(tree);      

      if((tree->rates->model == GUINDON) && (tree->mcmc->use_data))
	{
	  For(i,2*tree->n_otu-2) tree->rates->br_do_updt[i] = YES;
	  new_lnL_data = Lk(tree);
	}
      
      ratio +=
      	(new_lnL_rate - cur_lnL_rate);

      ratio +=
	(new_lnL_data - cur_lnL_data);


      /* !!!!!!!!!!!!!!!! */
      /* Modelling exp(nu) and making move on nu */
      /* ratio += (new_nu - cur_nu); */
	 
      /* Exponential prior on nu */
      /* ratio += LOG(Dexp(new_nu,10.) / Dexp(cur_nu,10.)); */

 
      ratio = EXP(ratio);
      alpha = MIN(1.,ratio);
      
      u = Uni();
      if(u > alpha) /* Reject */
	{
	  tree->rates->nu    = cur_nu;
	  tree->rates->c_lnL_rates = cur_lnL_rate;
	  tree->c_lnL        = cur_lnL_data;
	  Restore_Br_Len(tree);
	}
      else
	{
	  tree->mcmc->acc_move[tree->mcmc->num_move_nu]++;
	}
    }
  tree->mcmc->run_move[tree->mcmc->num_move_nu]++;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_All_Rates(t_tree *tree)
{
  phydbl cur_lnL_data, new_lnL_data, cur_lnL_rate;
  phydbl u, ratio, alpha;

  new_lnL_data = tree->c_lnL;
  cur_lnL_data = tree->c_lnL;
  cur_lnL_rate = tree->rates->c_lnL_rates;
  ratio        = 0.0;
  
  Record_Br_Len(tree);
  RATES_Record_Rates(tree);
  
  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[0],tree);
  MCMC_Sim_Rate(tree->n_root,tree->n_root->v[1],tree);

  new_lnL_data = Lk(tree);
  
  ratio += (new_lnL_data - cur_lnL_data);
  ratio = EXP(ratio);

  alpha = MIN(1.,ratio);
  u = Uni();

  if(u > alpha) /* Reject */
    {
      Restore_Br_Len(tree);
      RATES_Reset_Rates(tree);
      tree->rates->c_lnL_rates = cur_lnL_rate;
    }
  else
    {
      tree->rates->c_lnL_rates = RATES_Lk_Rates(tree);;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Only works when simulating from prior */
void MCMC_Sim_Rate(t_node *a, t_node *d, t_tree *tree)
{
  int err;
  phydbl mean,sd,br_r_a,dt_d;

  br_r_a = tree->rates->br_r[a->num];
  dt_d   = tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num];
  sd     = SQRT(dt_d*tree->rates->nu);
  
  if(tree->rates->model_log_rates == YES)
    {
      mean = br_r_a - .5*sd*sd;
    }
  else
    {
      mean = br_r_a;
    }
      
  if(tree->rates->model == STRICTCLOCK) tree->rates->br_r[d->num] = 1.0;
  else
    tree->rates->br_r[d->num] = Rnorm_Trunc(mean,sd,tree->rates->min_rate,tree->rates->max_rate,&err);

  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  MCMC_Sim_Rate(d,d->v[i],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Complete_MCMC(t_mcmc *mcmc, t_tree *tree)
{
  int i;
  phydbl sum;

  mcmc->io      = tree->io;
  mcmc->n_moves = 0;

  mcmc->num_move_nd_r = mcmc->n_moves;
  mcmc->n_moves += 2*tree->n_otu-1;

  mcmc->num_move_br_r = mcmc->n_moves;
  mcmc->n_moves += 2*tree->n_otu-2;

  mcmc->num_move_nd_t = mcmc->n_moves;
  mcmc->n_moves += tree->n_otu-1;

  mcmc->num_move_nu             = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_clock_r        = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_height    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_height = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_kappa          = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_tree_rates     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_subtree_rates  = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_updown_nu_cr   = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_ras            = mcmc->n_moves; mcmc->n_moves += 2*tree->mod->n_catg;
  mcmc->num_move_updown_t_cr    = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_cov_rates      = mcmc->n_moves; mcmc->n_moves += 2*tree->mod->m4mod->n_h;
  mcmc->num_move_cov_switch     = mcmc->n_moves; mcmc->n_moves += 1;
  mcmc->num_move_birth_rate     = mcmc->n_moves; mcmc->n_moves += 1;

  mcmc->run_move           = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->acc_move           = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->prev_run_move      = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->prev_acc_move      = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->acc_rate           = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->move_weight        = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->move_type          = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  
  /* TO DO: instead of n_moves here we should have something like n_param */
  mcmc->sum_val            = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->first_val          = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->sum_valsq          = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->sum_curval_nextval = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->ess                = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->ess_run            = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->start_ess          = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->new_param_val      = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->old_param_val      = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));
  mcmc->adjust_tuning      = (int *)mCalloc(mcmc->n_moves,sizeof(int));
  mcmc->tune_move          = (phydbl *)mCalloc(mcmc->n_moves,sizeof(phydbl));

  mcmc->move_name = (char **)mCalloc(mcmc->n_moves,sizeof(char *));
  For(i,mcmc->n_moves) mcmc->move_name[i] = (char *)mCalloc(50,sizeof(char));

  For(i,mcmc->n_moves) mcmc->adjust_tuning[i] = YES;

  for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) strcpy(mcmc->move_name[i],"br_rate");  
  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) strcpy(mcmc->move_name[i],"nd_rate");
  for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   strcpy(mcmc->move_name[i],"time");
  strcpy(mcmc->move_name[mcmc->num_move_nu],"nu");
  strcpy(mcmc->move_name[mcmc->num_move_clock_r],"clock");
  strcpy(mcmc->move_name[mcmc->num_move_tree_height],"tree_height");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_height],"subtree_height");
  strcpy(mcmc->move_name[mcmc->num_move_kappa],"kappa");
  strcpy(mcmc->move_name[mcmc->num_move_tree_rates],"tree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_subtree_rates],"subtree_rates");
  strcpy(mcmc->move_name[mcmc->num_move_updown_nu_cr],"updown_nu_cr");
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+2*tree->mod->n_catg;i++) 
    strcpy(mcmc->move_name[i],"ras");  
  strcpy(mcmc->move_name[mcmc->num_move_updown_t_cr],"updown_t_cr");
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+2*tree->mod->m4mod->n_h;i++) 
    strcpy(mcmc->move_name[i],"cov_rates");  
  strcpy(mcmc->move_name[mcmc->num_move_cov_switch],"cov_switch");
  strcpy(mcmc->move_name[mcmc->num_move_birth_rate],"birth_rate");
  
  if(tree->rates->model_log_rates == YES)
    for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_type[i] = MCMC_MOVE_RANDWALK_NORMAL;  
  else
    for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;  

  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_type[i] = MCMC_MOVE_RANDWALK_UNIFORM;
  /* for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE; */
  mcmc->move_type[mcmc->num_move_nu] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_clock_r] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_height] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_kappa] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_tree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_subtree_rates] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_updown_nu_cr] = MCMC_MOVE_RANDWALK_NORMAL;
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+2*tree->mod->n_catg;i++) mcmc->move_type[i] = MCMC_MOVE_RANDWALK_NORMAL;  
  mcmc->move_type[mcmc->num_move_updown_t_cr] = MCMC_MOVE_SCALE_THORNE;
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+2*tree->mod->m4mod->n_h;i++) mcmc->move_type[i] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_cov_switch] = MCMC_MOVE_SCALE_THORNE;
  mcmc->move_type[mcmc->num_move_birth_rate] = MCMC_MOVE_SCALE_THORNE;

  /* We start with small tuning parameter values in order to have inflated ESS
     for clock_r */
  For(i,mcmc->n_moves)
    {
      switch(mcmc->move_type[i])
	{
	case MCMC_MOVE_RANDWALK_NORMAL:
	  {
	    mcmc->tune_move[i] = 1.E-1;
	    break;
	  }
	case MCMC_MOVE_RANDWALK_UNIFORM:
	  {
	    mcmc->tune_move[i] = 2.0;
	    break;
	  }
	case MCMC_MOVE_SCALE_GAMMA:
	  {
	    mcmc->tune_move[i] = 1.0;
	    break;
	  }
	case MCMC_MOVE_SCALE_THORNE:
	  {
	    mcmc->tune_move[i] = 1.0;
	    break;
	  }
	}
    }
  
  for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_weight[i] = (phydbl)(1./(2.*tree->n_otu-2)); /* Rates */
  for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_weight[i] = 0.0; /* Node rates */
  for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_weight[i] = (phydbl)(1./(tree->n_otu-1));  /* Times */
  mcmc->move_weight[mcmc->num_move_clock_r]         = 1.0;
  mcmc->move_weight[mcmc->num_move_tree_height]     = 10.0;
  mcmc->move_weight[mcmc->num_move_subtree_height]  = 0.0;
  mcmc->move_weight[mcmc->num_move_nu]              = 2.0;
  mcmc->move_weight[mcmc->num_move_kappa]           = 0.5;
  mcmc->move_weight[mcmc->num_move_tree_rates]      = 1.0;
  mcmc->move_weight[mcmc->num_move_subtree_rates]   = 0.5;
  mcmc->move_weight[mcmc->num_move_updown_nu_cr]    = 0.0;
  for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+2*tree->mod->n_catg;i++) mcmc->move_weight[i] = 0.5*(1./(phydbl)tree->mod->n_catg);
  mcmc->move_weight[mcmc->num_move_updown_t_cr]     = 0.0; /* Does not seem to work well (does not give uniform prior on root height
  							      when sampling from prior) */
  for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+2*tree->mod->m4mod->n_h;i++) mcmc->move_weight[i] = 0.5*(1./(phydbl)tree->mod->m4mod->n_h);
  mcmc->move_weight[mcmc->num_move_cov_switch]      = 1.0;
  mcmc->move_weight[mcmc->num_move_birth_rate]      = 1.0;


  /* for(i=mcmc->num_move_br_r;i<mcmc->num_move_br_r+2*tree->n_otu-2;i++) mcmc->move_weight[i] = 0.0; /\* Rates *\/ */
  /* for(i=mcmc->num_move_nd_r;i<mcmc->num_move_nd_r+2*tree->n_otu-1;i++) mcmc->move_weight[i] = 0.0; /\* Node rates *\/ */
  /* for(i=mcmc->num_move_nd_t;i<mcmc->num_move_nd_t+tree->n_otu-1;i++)   mcmc->move_weight[i] = (phydbl)(1./(tree->n_otu-1));  /\* Times *\/ */
  /* mcmc->move_weight[mcmc->num_move_clock_r]         = 1.0; */
  /* mcmc->move_weight[mcmc->num_move_tree_height]     = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_subtree_height]  = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_nu]              = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_kappa]           = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_tree_rates]      = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_subtree_rates]   = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_updown_nu_cr]    = 0.0; */
  /* for(i=mcmc->num_move_ras;i<mcmc->num_move_ras+2*tree->mod->n_catg;i++) mcmc->move_weight[i] = 0.5*(1./(phydbl)tree->mod->n_catg);; */
  /* mcmc->move_weight[mcmc->num_move_updown_t_cr]     = 0.0; /\* Does not seem to work well (does not give uniform prior on root height */
  /* 							      when sampling from prior) *\/ */
  /* for(i=mcmc->num_move_cov_rates;i<mcmc->num_move_cov_rates+2*tree->mod->m4mod->n_h;i++) mcmc->move_weight[i] = 0.0; */
  /* mcmc->move_weight[mcmc->num_move_cov_switch]      = 0.0; */



  sum = 0.0;
  For(i,mcmc->n_moves) sum += mcmc->move_weight[i];
  For(i,mcmc->n_moves) mcmc->move_weight[i] /= sum;
  for(i=1;i<mcmc->n_moves;i++) mcmc->move_weight[i] += mcmc->move_weight[i-1];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Initialize_Param_Val(t_mcmc *mcmc, t_tree *tree)
{
  int i;

  mcmc->old_param_val[mcmc->num_move_nu]          = tree->rates->nu;
  mcmc->old_param_val[mcmc->num_move_clock_r]     = tree->rates->clock_r;
  mcmc->old_param_val[mcmc->num_move_tree_height] = tree->rates->nd_t[tree->n_root->num];
  mcmc->old_param_val[mcmc->num_move_kappa]       = tree->mod->kappa->v;

  For(i,2*tree->n_otu-2)
    mcmc->old_param_val[mcmc->num_move_br_r+i] = tree->rates->br_r[i];
  
  For(i,tree->n_otu-1)
    mcmc->old_param_val[mcmc->num_move_nd_t+i] = tree->rates->nd_t[tree->n_otu+i];

  For(i,2*tree->n_otu-1)
    mcmc->old_param_val[mcmc->num_move_nd_r+i] = tree->rates->nd_r[i];

  For(i,mcmc->n_moves) tree->mcmc->new_param_val[i] = tree->mcmc->old_param_val[i];
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Copy_To_New_Param_Val(t_mcmc *mcmc, t_tree *tree)
{
  int i;

  mcmc->new_param_val[mcmc->num_move_nu]          = tree->rates->nu;
  mcmc->new_param_val[mcmc->num_move_clock_r]     = tree->rates->clock_r;
  mcmc->new_param_val[mcmc->num_move_tree_height] = tree->rates->nd_t[tree->n_root->num];
  mcmc->new_param_val[mcmc->num_move_kappa]       = tree->mod->kappa->v;
  mcmc->new_param_val[mcmc->num_move_birth_rate]  = tree->rates->birth_rate;

  For(i,2*tree->n_otu-2)
    mcmc->new_param_val[mcmc->num_move_br_r+i] = tree->rates->br_r[i];
  
  For(i,tree->n_otu-1)
    mcmc->new_param_val[mcmc->num_move_nd_t+i] = tree->rates->nd_t[tree->n_otu+i];

  For(i,2*tree->n_otu-1)
    mcmc->new_param_val[mcmc->num_move_nd_r+i] = tree->rates->nd_r[i];
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Slice_One_Rate(t_node *a, t_node *d, int traversal, t_tree *tree)
{
  phydbl L,R; /* Left and Right limits of the slice */
  phydbl w; /* window width */
  phydbl u;
  phydbl x0,x1;
  phydbl logy;
  t_edge *b;
  int i;



  b = NULL;
  if(a == tree->n_root) b = tree->e_root;
  else For(i,3) if(d->v[i] == a) { b = d->b[i]; break; }
  
  w = 0.05;
  /* w = 10.; */

  x0 = tree->rates->br_r[d->num];
  logy = tree->c_lnL+tree->rates->c_lnL_rates - Rexp(1.);
  
  u = Uni();

  L = x0 - w*u;  
  R = L + w;
  
  do
    {
      tree->rates->br_r[d->num] = L;
      tree->rates->br_do_updt[d->num] = YES;
      RATES_Update_Cur_Bl(tree);
      Lk_At_Given_Edge(b,tree);
      RATES_Lk_Rates(tree);
      if(L < tree->rates->min_rate) { L = tree->rates->min_rate - w; break;}
      L = L - w;
    }
  while(tree->c_lnL + tree->rates->c_lnL_rates > logy);
  L = L + w;


  do
    {
      tree->rates->br_r[d->num] = R;
      tree->rates->br_do_updt[d->num] = YES;
      RATES_Update_Cur_Bl(tree);
      Lk_At_Given_Edge(b,tree);
      RATES_Lk_Rates(tree);
      if(R > tree->rates->max_rate) { R = tree->rates->max_rate + w; break;}
      R = R + w;
    }
  while(tree->c_lnL + tree->rates->c_lnL_rates > logy);
  R = R - w;


  for(;;)
    {
      u = Uni();
      x1 = L + u*(R-L);

      tree->rates->br_r[d->num] = x1;
      tree->rates->br_do_updt[d->num] = YES;
      RATES_Update_Cur_Bl(tree);
      Lk_At_Given_Edge(b,tree);
      RATES_Lk_Rates(tree);
      
      if(tree->c_lnL + tree->rates->c_lnL_rates > logy) break;
      
      if(x1 < x0) L = x1;
      else        R = x1;
    }


  if(traversal == YES)
    {
      if(d->tax == YES) return;
      else
	{
	  For(i,3)
	    if(d->v[i] != a && d->b[i] != tree->e_root)
	      {
		if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,d->b[i],d);
		/* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) { tree->both_sides = YES; Lk(tree); } */
		MCMC_Slice_One_Rate(d,d->v[i],YES,tree);
	      }
	}
      if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) Update_P_Lk(tree,b,d);
      /* if(tree->io->lk_approx == EXACT && tree->mcmc->use_data) { tree->both_sides = YES; Lk(tree); } */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void MCMC_Make_Move(phydbl *cur, phydbl *new, phydbl inf, phydbl sup, phydbl *loghr, phydbl tune, int move_type)
{
  phydbl u;

  switch(move_type)
    {

    case MCMC_MOVE_RANDWALK_NORMAL:
      {
	*new   = *cur + Rnorm(0.0,tune);
	/* Do not use reflection here */
	*loghr = 0.0;

	break;
      }

    case MCMC_MOVE_RANDWALK_UNIFORM:
      {
	u = Uni();
	/* *new   = u * (2.*tune) + (*cur) - tune; */
	/* *new   = Reflect(*new,inf,sup); */
	*new   = u*(sup-inf)+inf;
	*loghr = 0.0;
	break;
      }

    case MCMC_MOVE_SCALE_THORNE:
      {
	u = Uni();
	*new   = (*cur) * EXP(tune*(u-.5));
	*loghr = LOG((*new)/(*cur));
	break;
      }

    case MCMC_MOVE_SCALE_GAMMA:
      {
	phydbl r;
	*new = (*cur) * Rgamma(1./tune,tune);
	r = (*new)/(*cur);
	*loghr = -LOG(r) + LOG(Dgamma(1./r,1./tune,tune)/Dgamma(r,1./tune,tune));
	break;
      }

    default : 
      {
	PhyML_Printf("\n. Move not implemented");
	Exit("");
	break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void MCMC_Read_Param_Vals(t_tree *tree)
{
  char *token;
  int sizemax;
  FILE *in_fp;
  phydbl val;
  int i,v;

  in_fp = tree->mcmc->in_fp_par;
  

  sizemax = T_MAX_LINE;

  token = (char *)mCalloc(sizemax,sizeof(char));
  
  if(fgets(token,sizemax,in_fp)); // Skip first line
  if(fgets(token,sizemax,in_fp)); // Skip second
  
  v=fscanf(in_fp,"%lf\t",&val); // Run
  /* PhyML_Printf("\n. Run = %d",(int)val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnLike[Exact]
  /* PhyML_Printf("\n. LnLike = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnLike[Approx]
  /* PhyML_Printf("\n. LnLike = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnPriorRate
  /* PhyML_Printf("\n. LnPrior = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnPriorTime
  /* PhyML_Printf("\n. LnPrior = %f",val); */
  v=fscanf(in_fp,"%lf\t",&val); // LnPosterior
  /* PhyML_Printf("\n. LnPost = %f",val); */

  v=fscanf(in_fp,"%lf\t",&val); // ClockRate
  tree->rates->clock_r = val;
  /* PhyML_Printf("\n. Clock = %f",val); */

  v=fscanf(in_fp,"%lf\t",&val); // EvolRate

  v=fscanf(in_fp,"%lf\t",&val); // Nu
  tree->rates->nu = val;
  /* PhyML_Printf("\n. Nu = %f",val); */

  v=fscanf(in_fp,"%lf\t",&val); // Birth rate
  tree->rates->birth_rate = val;

  v=fscanf(in_fp,"%lf\t",&val); // TsTv
  tree->mod->kappa->v = val;
  /* PhyML_Printf("\n. TsTv = %f",val); */

  For(i,tree->n_otu-1)
    {
      v=fscanf(in_fp,"%lf\t",&val); // Node heights
      tree->rates->nd_t[i+tree->n_otu] = val;
    }

  For(i,2*tree->n_otu-2)
    {
      v=fscanf(in_fp,"%lf\t",&val); // Edge average rates
      tree->rates->br_r[i] = LOG(val);
    }
  
  
  Free(token);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

