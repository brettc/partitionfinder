/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "lk.h"


/* int    LIM_SCALE; */
/* phydbl LIM_SCALE_VAL; */
/* phydbl BIG; */
/* phydbl SMALL; */

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Tips_At_One_Site_Nucleotides_Float(char state, int pos, phydbl *p_lk)
{
  switch(state)
    {
    case 'A' : p_lk[pos+0]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'C' : p_lk[pos+1]=1.; p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'G' : p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'T' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'U' : p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+0]=.0;
      break;
    case 'M' : p_lk[pos+0]=p_lk[pos+1]=1.; p_lk[pos+2]=p_lk[pos+3]=.0;
      break;
    case 'R' : p_lk[pos+0]=p_lk[pos+2]=1.; p_lk[pos+1]=p_lk[pos+3]=.0;
      break;
    case 'W' : p_lk[pos+0]=p_lk[pos+3]=1.; p_lk[pos+1]=p_lk[pos+2]=.0;
      break;
    case 'S' : p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+0]=p_lk[pos+3]=.0;
      break;
    case 'Y' : p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+2]=.0;
      break;
    case 'K' : p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=p_lk[pos+1]=.0;
      break;
    case 'B' : p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+0]=.0;
      break;
    case 'D' : p_lk[pos+0]=p_lk[pos+2]=p_lk[pos+3]=1.; p_lk[pos+1]=.0;
      break;
    case 'H' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+3]=1.; p_lk[pos+2]=.0;
      break;
    case 'V' : p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=1.; p_lk[pos+3]=.0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_lk[pos+0]=p_lk[pos+1]=p_lk[pos+2]=p_lk[pos+3]=1.;break;
    default :
      {
	PhyML_Printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Tips_At_One_Site_Nucleotides_Int(char state, int pos, short int *p_pars)
{
  switch(state)
    {
    case 'A' : p_pars[pos+0]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'C' : p_pars[pos+1]=1; p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'G' : p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'T' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'U' : p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+0]=0;
      break;
    case 'M' : p_pars[pos+0]=p_pars[pos+1]=1; p_pars[pos+2]=p_pars[pos+3]=0;
      break;
    case 'R' : p_pars[pos+0]=p_pars[pos+2]=1; p_pars[pos+1]=p_pars[pos+3]=0;
      break;
    case 'W' : p_pars[pos+0]=p_pars[pos+3]=1; p_pars[pos+1]=p_pars[pos+2]=0;
      break;
    case 'S' : p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+0]=p_pars[pos+3]=0;
      break;
    case 'Y' : p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+2]=0;
      break;
    case 'K' : p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=p_pars[pos+1]=0;
      break;
    case 'B' : p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+0]=0;
      break;
    case 'D' : p_pars[pos+0]=p_pars[pos+2]=p_pars[pos+3]=1; p_pars[pos+1]=0;
      break;
    case 'H' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+3]=1; p_pars[pos+2]=0;
      break;
    case 'V' : p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=1; p_pars[pos+3]=0;
      break;
    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
      p_pars[pos+0]=p_pars[pos+1]=p_pars[pos+2]=p_pars[pos+3]=1;break;
    default :
      {
	PhyML_Printf("\n. Unknown character state : %c\n",state);
	Exit("\n. Init failed (check the data type)\n");
	break;
      }
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Tips_At_One_Site_AA_Float(char aa, int pos, phydbl *p_lk)
{
  int i;

  For(i,20) p_lk[pos+i] = .0;

  switch(aa){
  case 'A' : p_lk[pos+0]= 1.; break;/* Alanine */
  case 'R' : p_lk[pos+1]= 1.; break;/* Arginine */
  case 'N' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'D' : p_lk[pos+3]= 1.; break;/* Aspartic acid */
  case 'C' : p_lk[pos+4]= 1.; break;/* Cysteine */
  case 'Q' : p_lk[pos+5]= 1.; break;/* Glutamine */
  case 'E' : p_lk[pos+6]= 1.; break;/* Glutamic acid */
  case 'G' : p_lk[pos+7]= 1.; break;/* Glycine */
  case 'H' : p_lk[pos+8]= 1.; break;/* Histidine */
  case 'I' : p_lk[pos+9]= 1.; break;/* Isoleucine */
  case 'L' : p_lk[pos+10]=1.; break;/* Leucine */
  case 'K' : p_lk[pos+11]=1.; break;/* Lysine */
  case 'M' : p_lk[pos+12]=1.; break;/* Methionine */
  case 'F' : p_lk[pos+13]=1.; break;/* Phenylalanin */
  case 'P' : p_lk[pos+14]=1.; break;/* Proline */
  case 'S' : p_lk[pos+15]=1.; break;/* Serine */
  case 'T' : p_lk[pos+16]=1.; break;/* Threonine */
  case 'W' : p_lk[pos+17]=1.; break;/* Tryptophan */
  case 'Y' : p_lk[pos+18]=1.; break;/* Tyrosine */
  case 'V' : p_lk[pos+19]=1.; break;/* Valine */

  case 'B' : p_lk[pos+2]= 1.; break;/* Asparagine */
  case 'Z' : p_lk[pos+5]= 1.; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_lk[pos+i] = 1.; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Tips_At_One_Site_AA_Int(char aa, int pos, short int *p_pars)
{
  int i;

  For(i,20) p_pars[pos+i] = .0;

  switch(aa){
  case 'A' : p_pars[pos+0]  = 1; break;/* Alanine */
  case 'R' : p_pars[pos+1]  = 1; break;/* Arginine */
  case 'N' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'D' : p_pars[pos+3]  = 1; break;/* Aspartic acid */
  case 'C' : p_pars[pos+4]  = 1; break;/* Cysteine */
  case 'Q' : p_pars[pos+5]  = 1; break;/* Glutamine */
  case 'E' : p_pars[pos+6]  = 1; break;/* Glutamic acid */
  case 'G' : p_pars[pos+7]  = 1; break;/* Glycine */
  case 'H' : p_pars[pos+8]  = 1; break;/* Histidine */
  case 'I' : p_pars[pos+9]  = 1; break;/* Isoleucine */
  case 'L' : p_pars[pos+10] = 1; break;/* Leucine */
  case 'K' : p_pars[pos+11] = 1; break;/* Lysine */
  case 'M' : p_pars[pos+12] = 1; break;/* Methionine */
  case 'F' : p_pars[pos+13] = 1; break;/* Phenylalanin */
  case 'P' : p_pars[pos+14] = 1; break;/* Proline */
  case 'S' : p_pars[pos+15] = 1; break;/* Serine */
  case 'T' : p_pars[pos+16] = 1; break;/* Threonine */
  case 'W' : p_pars[pos+17] = 1; break;/* Tryptophan */
  case 'Y' : p_pars[pos+18] = 1; break;/* Tyrosine */
  case 'V' : p_pars[pos+19] = 1; break;/* Valine */

  case 'B' : p_pars[pos+2]  = 1; break;/* Asparagine */
  case 'Z' : p_pars[pos+5]  = 1; break;/* Glutamine */

  case 'X' : case '?' : case '-' : For(i,20) p_pars[pos+i] = 1; break;
  default :
    {
      PhyML_Printf("\n. Unknown character state : %c\n",aa);
      Exit("\n. Init failed (check the data type)\n");
      break;
    }
  }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Tips_At_One_Site_Generic_Float(char *state, int ns, int state_len, int pos, phydbl *p_lk)
{
  int i;
  int state_int;

  For(i,ns) p_lk[pos+i] = 0.;

  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_lk[pos+i] = 1.;
  else
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
	{
	  PhyML_Printf("\n. state='%c'",state);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(state_int > ns)
	{
	  PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");	  
	}
      p_lk[pos+state_int] = 1.;
      /*       PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d ns: %d pos: %d",__FILE__,__LINE__,state,state_int,state_len,ns,pos); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_Tips_At_One_Site_Generic_Int(char *state, int ns, int state_len, int pos, short int *p_pars)
{
  int i;
  int state_int;

  For(i,ns) p_pars[pos+i] = 0;
  
  if(Is_Ambigu(state,GENERIC,state_len)) For(i,ns) p_pars[pos+i] = 1;
  else 
    {
      char format[6];
      sprintf(format,"%%%dd",state_len);
      if(!sscanf(state,format,&state_int))
	{
	  PhyML_Printf("\n. state='%c'",state);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      if(state_int > ns)
	{
	  PhyML_Printf("\n. %s %d cstate: %.2s istate: %d state_len: %d.\n",__FILE__,__LINE__,state,state_int,state_len);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");	  
	}
      p_pars[pos+state_int] = 1;
/*       PhyML_Printf("\n* %s %d cstate: %.2s istate: %d state_len: %d ns: %d pos: %d",__FILE__,__LINE__,state,state_int,state_len,ns,pos); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Get_All_Partial_Lk_Scale(t_tree *tree, t_edge *b_fcus, t_node *a, t_node *d)
{
  Update_P_Lk(tree,b_fcus,d);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Post_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i,dir;

  dir = -1;
  
  if(d->tax) 
    {
      Get_All_Partial_Lk_Scale(tree,d->b[0],a,d);
      return;
    }
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Post_Order_Lk(d,d->v[i],tree);
	  else dir = i;
	}      
      Get_All_Partial_Lk_Scale(tree,d->b[dir],a,d);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Pre_Order_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) 
    {
      Get_All_Partial_Lk_Scale(tree,d->b[0],a,d);
      return;
    }
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Get_All_Partial_Lk_Scale(tree,d->b[i],d->v[i],d);
	      Pre_Order_Lk(d,d->v[i],tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk(t_tree *tree)
{
  int br;
  int n_patterns;

  tree->old_lnL = tree->c_lnL;

#ifdef PHYTIME
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
#endif

  Check_Br_Len_Bounds(tree);
  
  if(tree->rates && tree->io->lk_approx == NORMAL)
    {
      tree->c_lnL = Lk_Normal_Approx(tree);
      return tree->c_lnL;
    }

  n_patterns = tree->n_pattern;

  Set_Model_Parameters(tree->mod);
  
  For(br,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->t_edges[br],tree);

  Post_Order_Lk(tree->t_nodes[0],tree->t_nodes[0]->v[0],tree);
  if(tree->both_sides)
    Pre_Order_Lk(tree->t_nodes[0],
		 tree->t_nodes[0]->v[0],
		 tree);

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;
  For(tree->curr_site,n_patterns)
    if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(tree->t_nodes[0]->b[0],tree);

/*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); */

/*   tree->c_lnL = .0; */
/*   For(tree->curr_site,n_patterns) */
/*     { */
/*       tree->c_lnL += tree->c_lnL_sorted[tree->curr_site]; */
/*     } */

    
  Adjust_Min_Diff_Lk(tree);

  tree->c_lnL_mixt = tree->c_lnL;

  if(tree->nextree) 
    {
      Copy_Edge_Lengths(tree->nextree,tree);
      Lk(tree->nextree);
    }
  else Lk_LastFirst(tree);

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  int n_patterns;
  
  n_patterns = tree->n_pattern;

#ifdef PHYTIME
  if((tree->rates) && (tree->rates->bl_from_rt)) RATES_Update_Cur_Bl(tree);
#endif
  
  Check_Br_Len_Bounds(tree);

  if(tree->rates && tree->io->lk_approx == NORMAL)
    {
      tree->c_lnL = Lk_Normal_Approx(tree);
      return tree->c_lnL;
    }

  Update_PMat_At_Given_Edge(b_fcus,tree);

  if(b_fcus->left->tax)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->c_lnL             = .0;
  tree->sum_min_sum_scale = .0;
  For(tree->curr_site,n_patterns)
    {
      if(tree->data->wght[tree->curr_site] > SMALL) Lk_Core(b_fcus,tree);
    }

/*   Qksort(tree->c_lnL_sorted,NULL,0,n_patterns-1); */

/*   tree->c_lnL = .0; */
/*   For(tree->curr_site,n_patterns) */
/*     { */
/*       tree->c_lnL += tree->c_lnL_sorted[tree->curr_site]; */
/*     } */


  Adjust_Min_Diff_Lk(tree);

  tree->c_lnL_mixt = tree->c_lnL;
  
  if(tree->nextree) 
    {
      tree->nextree->t_edges[b_fcus->num]->l = b_fcus->l;
      Lk_At_Given_Edge(tree->nextree->t_edges[b_fcus->num],tree->nextree);
    }
  else Lk_LastFirst(tree);

  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

/* Core of the likelihood calculcation. Assume that the partial likelihoods on both
   sides of t_edge *b are up-to-date. Calculate the log-likelihood at one site.
*/

#ifndef USE_OLD_LK

phydbl Lk_Core(t_edge *b, t_tree *tree)
{
  phydbl log_site_lk;
  phydbl site_lk_cat, site_lk;
  int sum_scale_left, sum_scale_rght;
  int fact_sum_scale;
  phydbl max_sum_scale;
  phydbl sum;
  int ambiguity_check,state;
  int catg,ns,k,l,site;
  int dim1,dim2,dim3;
  int *sum_scale_left_cat,*sum_scale_rght_cat;
  phydbl multiplier;
  int exponent, piecewise_exponent;
  phydbl tmp;
  phydbl logbig;
  phydbl inv_site_lk;


  logbig = LOG((phydbl)BIG);

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  log_site_lk     = .0;
  site_lk         = .0;
  site_lk_cat     = .0;
  ambiguity_check = -1;
  state           = -1;
  site            = tree->curr_site;
  ns              = tree->mod->ns;

  if((b->rght->tax) && (!tree->mod->s_opt->greedy))
    {
      ambiguity_check = tree->data->c_seq[b->rght->num]->is_ambigu[site];
      if(!ambiguity_check) state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
    }

  sum_scale_left = .0;
  sum_scale_rght = .0;
  
  if(tree->mod->use_m4mod) ambiguity_check = 1;

  
  sum_scale_left_cat = b->sum_scale_left_cat;
  sum_scale_rght_cat = b->sum_scale_rght_cat;


  /* Actual likelihood calculation */
  /* For all classes of rates */
  For(catg,tree->mod->n_catg)
    {
      site_lk_cat = .0;

      /* b is an external edge */
      if((b->rght->tax) && (!tree->mod->s_opt->greedy))
        {
          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
          if(!ambiguity_check)
            {
              sum = .0;
              For(l,ns)
                {
                  sum +=
                    b->Pij_rr[catg*dim3+state*dim2+l] *
                    b->p_lk_left[site*dim1+catg*dim2+l];
		  
		  /* if(tree->curr_site == 0) printf("\n. %f %f %d",  */
		  /* 				  b->Pij_rr[catg*dim3+state*dim2+l],                     */
		  /* 				  b->p_lk_left[site*dim1+catg*dim2+l],b->left->num); */
                }

              site_lk_cat += sum * tree->mod->pi->v[state];	
	    }
          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
          else
            {
              For(k,ns)
                {
                  sum = .0;
                  if(b->p_lk_tip_r[site*dim2+k] > .0)
                    {
                      For(l,ns)
                        {
                          sum +=
                            b->Pij_rr[catg*dim3+k*dim2+l] *
                            b->p_lk_left[site*dim1+catg*dim2+l];
                        }

                      site_lk_cat +=
                        sum *
                        tree->mod->pi->v[k] *
                        b->p_lk_tip_r[site*dim2+k];
		    }
                }
            }
        }
      /* b is an internal edge: ns x ns terms to consider */
      else
        {
          For(k,ns)
            {
              sum = .0;
              if(b->p_lk_rght[site*dim1+catg*dim2+k] > .0)
                {
                  For(l,ns)
                    {
                      sum +=
                        b->Pij_rr[catg*dim3+k*dim2+l] *
                        b->p_lk_left[site*dim1+catg*dim2+l];
                    }

                  site_lk_cat +=
                    sum *
                    tree->mod->pi->v[k] *
                    b->p_lk_rght[site*dim1+catg*dim2+k];
		}
            }
	}
      tree->site_lk_cat[catg] = site_lk_cat;
    }
    
  max_sum_scale =  (phydbl)BIG;
  For(catg,tree->mod->n_catg)
    {
      sum_scale_left_cat[catg] =
        (b->sum_scale_left)?
        (b->sum_scale_left[catg*tree->n_pattern+site]):
        (0.0);
      
      sum_scale_rght_cat[catg] =
        (b->sum_scale_rght)?
        (b->sum_scale_rght[catg*tree->n_pattern+site]):
        (0.0);

      sum = sum_scale_left_cat[catg] + sum_scale_rght_cat[catg];

      if(sum < .0)
	{
	  PhyML_Printf("\n. sum = %G",sum);
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Warn_And_Exit("\n");
	}

      tmp = sum + (logbig - LOG(tree->site_lk_cat[catg]))/(phydbl)LOG2;
      if(tmp < max_sum_scale) max_sum_scale = tmp; /* min of the maxs */
    }

/*   fact_sum_scale = (int)((max_sum_scale + min_sum_scale) / 2); */

  fact_sum_scale = (int)(max_sum_scale / 2);

  /* Apply scaling factors */
  For(catg,tree->mod->n_catg)
    {
      exponent = -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale;

      site_lk_cat = tree->site_lk_cat[catg];
     
      if(exponent > 0)
	{
	  do
	    {
	      piecewise_exponent = MIN(exponent,63);
	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	      site_lk_cat *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	}
      else
	{
	  do
	    {
	      piecewise_exponent = MAX(exponent,-63);
	      multiplier = 1. / (phydbl)((unsigned long long)(1) << -piecewise_exponent);
	      site_lk_cat *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	}

      if(isinf(site_lk_cat))
	{
	  PhyML_Printf("\n+ site=%4d cat=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
		       tree->curr_site,
		       catg,
		       tree->site_lk_cat[catg],
		       sum_scale_left_cat[catg]+sum_scale_rght_cat[catg],
		       max_sum_scale,
		       fact_sum_scale,
		       -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale,
		       (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale));

	  site_lk_cat = BIG / 10;
	}

      
      if(site_lk_cat < SMALL)
	{
	  if(tree->mod->n_catg == 1)
	    {
	      PhyML_Printf("\n- site=%4d cat=%4d site_lk_cat=%G sum_scale=%d max=%G fact=%d expo=%d dbl=%G",
			   tree->curr_site,
			   catg,
			   tree->site_lk_cat[catg],
			   sum_scale_left_cat[catg]+sum_scale_rght_cat[catg],
			   max_sum_scale,
			   fact_sum_scale,
			   -(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale,
			   (double)tree->site_lk_cat[catg] * pow(2.,-(sum_scale_left_cat[catg]+sum_scale_rght_cat[catg])+fact_sum_scale));

	      Exit("\n");
	    }

	  site_lk_cat = .0;
	}

      tree->site_lk_cat[catg] = site_lk_cat;
    }

  site_lk = .0;
  For(catg,tree->mod->n_catg)
    {
      site_lk += tree->site_lk_cat[catg] * tree->mod->gamma_r_proba->v[catg];
    }

  inv_site_lk = 0.;
  
  /* The substitution model does include invariable sites */
  if(tree->mod->invar)
    {
      /* The site is invariant */
      if(tree->data->invar[site] > -0.5)
        {
	  /* Multiply P(D|r=0) by 2^(fact_sum_scale) */

	  inv_site_lk = tree->mod->pi->v[tree->data->invar[site]];
	  exponent = fact_sum_scale;
	  do
	    {
	      piecewise_exponent = MIN(exponent,63);
	      multiplier = (phydbl)((unsigned long long)(1) << piecewise_exponent);
	      inv_site_lk *= multiplier;
	      exponent -= piecewise_exponent;
	    }
	  while(exponent != 0);
	  inv_site_lk *= tree->mod->pinvar->v;
	  /* Update the value of site_lk */
	  if(isinf(inv_site_lk)) // P(D|r=0) >> P(D|r>0) => assume P(D) = P(D|r=0)P(r=0)
	    {
	      PhyML_Printf("\n== Numerical precision issue alert.");
	      PhyML_Printf("\n== P(D|r=0)=%G P(D|r>0)=%G",
			   tree->mod->pi->v[tree->data->invar[site]],
			   EXP(LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale));	      
	      site_lk = (1. - tree->mod->pinvar->v)*tree->mod->pi->v[tree->data->invar[site]];
	      fact_sum_scale = 0.0;
	    }
	  else
	    {
	      site_lk *= (1. - tree->mod->pinvar->v);
	      site_lk += inv_site_lk;	      
	    }
	}
      else
        {
          /* Same formula as above with P(D | subs, rate = 0) = 0 */
	  site_lk *= (1. - tree->mod->pinvar->v);
        }
    }

  log_site_lk = LOG(site_lk) - (phydbl)LOG2 * fact_sum_scale;

  For(catg,tree->mod->n_catg) tree->log_site_lk_cat[catg][site] = LOG(tree->site_lk_cat[catg]) - (phydbl)LOG2 * fact_sum_scale;
  
  if(isinf(log_site_lk) || isnan(log_site_lk))
    {
      PhyML_Printf("\n. Site = %d",site);
      PhyML_Printf("\n. invar = %d",tree->data->invar[site]);
      PhyML_Printf("\n. scale_left = %d scale_rght = %d",sum_scale_left,sum_scale_rght);
      PhyML_Printf("\n. inv_site_lk = %f",inv_site_lk);
      PhyML_Printf("\n. Lk = %G LOG(Lk) = %f < %G",site_lk,log_site_lk,-BIG);
      For(catg,tree->mod->n_catg) PhyML_Printf("\n. rr=%f p=%f",tree->mod->gamma_rr->v[catg],tree->mod->gamma_r_proba->v[catg]);
      PhyML_Printf("\n. pinv = %G",tree->mod->pinvar->v);
      PhyML_Printf("\n. bl mult = %G",tree->mod->br_len_multiplier->v);

      /* int i; */
      /* For(i,2*tree->n_otu-3) */
      /* 	{ */
      /* 	  PhyML_Printf("\n. b%3d->l = %f %f [%G] %f %f %f %f [%s]",i, */
      /* 		       tree->t_edges[i]->l, */
      /* 		       tree->t_edges[i]->gamma_prior_mean, */
      /* 		       tree->t_edges[i]->gamma_prior_var, */
      /* 		       tree->rates->nd_t[tree->t_edges[i]->left->num], */
      /* 		       tree->rates->nd_t[tree->t_edges[i]->rght->num], */
      /* 		       tree->rates->br_r[tree->t_edges[i]->left->num], */
      /* 		       tree->rates->br_r[tree->t_edges[i]->rght->num], */
      /* 		       tree->t_edges[i] == tree->e_root ? "YES" : "NO"); */
      /* 	  fflush(NULL); */
		       
      /* 	} */

      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }
  
  tree->cur_site_lk[site] = log_site_lk;

  /* Multiply log likelihood by the number of times this site pattern is found in the data */
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;

  tree->c_lnL += tree->data->wght[site]*log_site_lk;
/*   tree->sum_min_sum_scale += (int)tree->data->wght[site]*min_sum_scale; */
  
  return log_site_lk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/

void Update_P_Lk(t_tree *tree, t_edge *b, t_node *d)
{

  if((tree->io->do_alias_subpatt == YES) && 
     (tree->update_alias_subpatt == YES)) 
    Alias_One_Subpatt((d==b->left)?(b->rght):(b->left),d,tree);

  if(d->tax) return;
  
  if(tree->mod->use_m4mod == NO)
    {
      if(tree->io->datatype == NT)
	{
	  Update_P_Lk_Nucl(tree,b,d);
	}
      else if(tree->io->datatype == AA)
	{
	  Update_P_Lk_AA(tree,b,d);
	}
    }
  else
    {
      Update_P_Lk_Generic(tree,b,d);
    }

  if(tree->nextree) Update_P_Lk(tree->nextree,
				tree->nextree->t_edges[b->num],
				tree->nextree->t_nodes[d->num]);
  
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Update_P_Lk_Generic(t_tree *tree, t_edge *b, t_node *d)
{
/*
           |
	   |<- b
	   |
	   d
          / \
       	 /   \
       	/     \
	n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i,j;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int NsNg, Ns, NsNs;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;
  int *p_lk_loc;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;
 
  NsNg = tree->mod->n_catg * tree->mod->ns;
  Ns = tree->mod->ns;
  NsNs = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;
  
  /* TO DO: Might be worth keeping these directions in memory instead of
     calculating them every time... */
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);

/*   if((dir1 == -1) || (dir2 == -1)) */
/*     { */
/*       PhyML_Printf("\n. d = %d",d->num); */
/*       PhyML_Printf("\n. d->v[0] = %d, d->v[1] = %d, d->v[2] = %d",d->v[0]->num,d->v[1]->num,d->v[2]->num); */
/*       PhyML_Printf("\n. d->b[0] = %d, d->b[1] = %d, d->b[2] = %d",d->b[0]->num,d->b[1]->num,d->b[2]->num); */
/*       PhyML_Printf("\n. d->num = %d dir1 = %d dir2 = %d",d->num,dir1,dir2); */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Exit(""); */
/*     } */
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_left;
      p_lk_loc = b->p_lk_loc_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_rght;
      p_lk_loc = b->p_lk_loc_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
      sum_scale_v1 = d->b[dir1]->sum_scale_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
      sum_scale_v1 = d->b[dir1]->sum_scale_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
      sum_scale_v2 = d->b[dir2]->sum_scale_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
      sum_scale_v2 = d->b[dir2]->sum_scale_left;
    }
  
  /* Change probability matrices on the two pendant edges */
  Pij1 = d->b[dir1]->Pij_rr;
  Pij2 = d->b[dir2]->Pij_rr;
  
  /* For every site in the alignment */
  For(site,n_patterns)
    {
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;
      
      if(!tree->mod->s_opt->greedy)
	{
	  /* n_v1 and n_v2 are tip nodes */
	  if(n_v1->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
	      if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*Ns,tree);
	    }
	      
	  if(n_v2->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
	      if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*Ns,tree);
	    }
	}
      
      if(tree->mod->use_m4mod)
	{
	  ambiguity_check_v1 = YES;
	  ambiguity_check_v2 = YES;
	}

      if(p_lk_loc[site] < site)
	{
	  Copy_P_Lk(p_lk,p_lk_loc[site],site,tree);
	  Copy_Scale(sum_scale,p_lk_loc[site],site,tree);
	}
      else
	{
	  /* For all the rate classes */
	  For(catg,tree->mod->n_catg)
	    {
	      smallest_p_lk  =  BIG;
	      
	      /* For all the states at node d */
	      For(i,tree->mod->ns)
		{
		  p1_lk1 = .0;
		  
		  /* n_v1 is a tip */
		  if((n_v1->tax) && (!tree->mod->s_opt->greedy))
		    {
		      if(ambiguity_check_v1 == NO)
			{
			  /* For the (non-ambiguous) state at node n_v1 */
			  p1_lk1 = Pij1[catg*NsNs+i*Ns+state_v1];
			}
		      else
			{
			  /* For all the states at node n_v1 */
			  For(j,tree->mod->ns)
			    {
			      p1_lk1 += Pij1[catg*NsNs+i*Ns+j] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*Ns+j];
			    }
			}
		    }
		  /* n_v1 is an internal node */
		  else
		    {
		      /* For the states at node n_v1 */
		      For(j,tree->mod->ns)
			{
			  p1_lk1 += Pij1[catg*NsNs+i*Ns+j] * p_lk_v1[site*NsNg+catg*Ns+j];
			}
		    }
		  
		  p2_lk2 = .0;
		  
		  /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
		  /* n_v2 is a tip */
		  if((n_v2->tax) && (!tree->mod->s_opt->greedy))
		    {
		      if(ambiguity_check_v2 == NO)
			{
			  /* For the (non-ambiguous) state at node n_v2 */
			  p2_lk2 = Pij2[catg*NsNs+i*Ns+state_v2];
			}
		      else
			{
			  /* For all the states at node n_v2 */
			  For(j,tree->mod->ns)
			    {
			      p2_lk2 += Pij2[catg*NsNs+i*Ns+j] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*Ns+j];
			    }
			}
		    }
		  /* n_v2 is an internal node */
		  else
		    {
		      /* For all the states at node n_v2 */
		      For(j,tree->mod->ns)
			{
			  p2_lk2 += Pij2[catg*NsNs+i*Ns+j] * p_lk_v2[site*NsNg+catg*Ns+j];
			}
		    }
		  
		  p_lk[site*NsNg+catg*Ns+i] = p1_lk1 * p2_lk2;
		  
		  /* 	      PhyML_Printf("\n+ %G",p_lk[site*NsNg+catg*Ns+i]); */
		  
		  if(p_lk[site*NsNg+catg*Ns+i] < smallest_p_lk) smallest_p_lk = p_lk[site*NsNg+catg*Ns+i] ;
		}
	      
	      /* Current scaling values at that site */
	      sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
	      sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);
	      
	      sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
	      
	      /* 	  PhyML_Printf("\n@ %d",sum_scale[catg*n_patterns+site]); */
	      
	      /* Scaling */
	      if(smallest_p_lk < p_lk_lim_inf)
		{
		  curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
		  curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);
		  
		  /* 	      if(fabs(curr_scaler_pow) > 63 || fabs(curr_scaler_pow) > 63) */
		  /* 		{ */
		  /* 		  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk); */
		  /* 		  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow); */
		  /* 		  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__); */
		  /* 		  Warn_And_Exit("\n"); */
		  /* 		} */
		  
		  sum_scale[catg*n_patterns+site] += curr_scaler_pow;
		  
		  do
		    {
		      piecewise_scaler_pow = MIN(curr_scaler_pow,63);
		      curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
		      For(i,tree->mod->ns)
			{
			  p_lk[site*NsNg+catg*Ns+i] *= curr_scaler;
			  
			  if(p_lk[site*NsNg+catg*Ns+i] > BIG)
			    {
			      PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
			      PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
			      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			      Warn_And_Exit("\n");
			    }
			}
		      curr_scaler_pow -= piecewise_scaler_pow;
		    }
		  while(curr_scaler_pow != 0);
		}
	    }
	}
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_P_Lk_Nucl(t_tree *tree, t_edge *b, t_node *d)
{
/*
           |
	   |<- b
	   |
	   d
          / \
       	 /   \
       	/     \
	n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;
  phydbl p0,p1,p2,p3;
  int *p_lk_loc;


  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;
 
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;
  
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);

  if((dir1 == -1) || (dir2 == -1))
    {
      PhyML_Printf("\n. d = %d",d->num);
      PhyML_Printf("\n. d->v[0] = %d, d->v[1] = %d, d->v[2] = %d",d->v[0]->num,d->v[1]->num,d->v[2]->num);
      PhyML_Printf("\n. d->b[0] = %d, d->b[1] = %d, d->b[2] = %d",d->b[0]->num,d->b[1]->num,d->b[2]->num);
      PhyML_Printf("\n. d->num = %d dir1 = %d dir2 = %d",d->num,dir1,dir2);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Exit("");
    }
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];


  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_left;
      p_lk_loc = b->p_lk_loc_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_rght;
      p_lk_loc = b->p_lk_loc_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
      sum_scale_v1 = d->b[dir1]->sum_scale_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
      sum_scale_v1 = d->b[dir1]->sum_scale_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
      sum_scale_v2 = d->b[dir2]->sum_scale_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
      sum_scale_v2 = d->b[dir2]->sum_scale_left;
    }
  
  /* Change probability matrices on the two pendant edges */
  Pij1 = d->b[dir1]->Pij_rr;
  Pij2 = d->b[dir2]->Pij_rr;
  
  /* For every site in the alignment */
  For(site,n_patterns)
    {
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;
      
      if(!tree->mod->s_opt->greedy)
	{
	  /* n_v1 and n_v2 are tip nodes */
	  if(n_v1->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
	      if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
	      
	  if(n_v2->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
	      if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
	}
      
      if(tree->mod->use_m4mod)
	{
	  ambiguity_check_v1 = YES;
	  ambiguity_check_v2 = YES;
	}

      if(p_lk_loc[site] < site)
	{
	  Copy_P_Lk(p_lk,p_lk_loc[site],site,tree);
	  Copy_Scale(sum_scale,p_lk_loc[site],site,tree);
	}
      else
	{
	  /* For all the rate classes */
	  For(catg,tree->mod->n_catg)
	    {
	      smallest_p_lk  =  BIG;
	      
	      /* For all the state at node d */
	      For(i,tree->mod->ns)
		{
		  p1_lk1 = .0;
		  
		  /* n_v1 is a tip */
		  if((n_v1->tax) && (!tree->mod->s_opt->greedy))
		    {
		      if(ambiguity_check_v1 == NO)
			{
			  /* For the (non-ambiguous) state at node n_v1 */
			  p1_lk1 = Pij1[catg*dim3+i*dim2+state_v1];
			}
		      else
			{
			  /* For all the states at node n_v1 */
			  p0=Pij1[catg*dim3+i*dim2+0] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+0];
			  p1=Pij1[catg*dim3+i*dim2+1] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+1];
			  p2=Pij1[catg*dim3+i*dim2+2] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+2];
			  p3=Pij1[catg*dim3+i*dim2+3] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+3];
			  p1_lk1 = p0+p1+p2+p3;

			  if(isnan(p1_lk1))
			    {
			      PhyML_Printf("\n. p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3);
			      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			      Warn_And_Exit("\n");		  
			    }
			}
		    }
		  /* n_v1 is an internal node */
		  else
		    {
		      /* For the states at node n_v1 */
		      /* 		  For(j,tree->mod->ns) */
		      /* 		    { */
		      /* 		      p1_lk1 += Pij1[catg*dim3+i*dim2+j] * p_lk_v1[site*dim1+catg*dim2+j]; */
		      /* 		    } */
		      
		      p0=Pij1[catg*dim3+i*dim2+0] * p_lk_v1[site*dim1+catg*dim2+0];
		      p1=Pij1[catg*dim3+i*dim2+1] * p_lk_v1[site*dim1+catg*dim2+1];
		      p2=Pij1[catg*dim3+i*dim2+2] * p_lk_v1[site*dim1+catg*dim2+2];
		      p3=Pij1[catg*dim3+i*dim2+3] * p_lk_v1[site*dim1+catg*dim2+3];
		      p1_lk1 = p0+p1+p2+p3;
		      if(isnan(p1_lk1))
			{
			  PhyML_Printf("\n. p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3);
			  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			  Warn_And_Exit("\n");		  
			}
		    }
		  
		  p2_lk2 = .0;
		  
		  /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
		  /* n_v2 is a tip */
		  if((n_v2->tax) && (!tree->mod->s_opt->greedy))
		    {
		      if(ambiguity_check_v2 == NO)
			{
			  /* For the (non-ambiguous) state at node n_v2 */
			  p2_lk2 = Pij2[catg*dim3+i*dim2+state_v2];
			  if(isnan(p2_lk2))
			    {
			      PhyML_Printf("\n. Tree %d",tree->tree_num);
			      PhyML_Printf("\n. catg=%d dim3=%d dim2=%d i=%d state_v2=%d",catg,dim3,dim2,i,state_v2);
			      PhyML_Printf("\n. Pij2[0] = %f",Pij2[0]);
			      PhyML_Printf("\n. q[0]=%f",tree->mod->eigen->q[0]);
			      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			      Warn_And_Exit("\n");		  
			    }
			}
		      else
			{
			  /* For all the states at node n_v2 */
			  p0=Pij2[catg*dim3+i*dim2+0] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+0];
			  p1=Pij2[catg*dim3+i*dim2+1] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+1];
			  p2=Pij2[catg*dim3+i*dim2+2] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+2];
			  p3=Pij2[catg*dim3+i*dim2+3] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+3];
			  p2_lk2 = p0+p1+p2+p3;
			  if(isnan(p2_lk2))
			    {
			      PhyML_Printf("\n. p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3);
			      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			      Warn_And_Exit("\n");		  
			    }
			}
		    }
		  /* n_v2 is an internal node */
		  else
		    {
		      /* For all the states at node n_v2 */
		      p0=Pij2[catg*dim3+i*dim2+0] * p_lk_v2[site*dim1+catg*dim2+0];
		      p1=Pij2[catg*dim3+i*dim2+1] * p_lk_v2[site*dim1+catg*dim2+1];
		      p2=Pij2[catg*dim3+i*dim2+2] * p_lk_v2[site*dim1+catg*dim2+2];
		      p3=Pij2[catg*dim3+i*dim2+3] * p_lk_v2[site*dim1+catg*dim2+3];
		      p2_lk2 = p0+p1+p2+p3;
		      if(isnan(p2_lk2))
			{
			  PhyML_Printf("\n. p0=%f p1=%f p2=%f p3=%f",p0,p1,p2,p3);
			  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			  Warn_And_Exit("\n");		  
			}
		    }
		  
		  p_lk[site*dim1+catg*dim2+i] = p1_lk1 * p2_lk2;	    
		  
		  
		  if(isnan(p2_lk2))
		    {
		      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
		      Warn_And_Exit("\n");		  
		    }

		  if(p_lk[site*dim1+catg*dim2+i] < smallest_p_lk) smallest_p_lk = p_lk[site*dim1+catg*dim2+i] ; 
		}
	      
	      

	      /* Current scaling values at that site */
	      sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
	      sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);
	      
	      sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
	      
	      /* Scaling */
	      if(smallest_p_lk < p_lk_lim_inf)
		{
		  curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
		  curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);
		  
		  /* 	      if(fabs(curr_scaler_pow) > 63 || fabs(curr_scaler_pow) > 63) */
		  /* 		{ */
		  /* 		  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk); */
		  /* 		  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow); */
		  /* 		  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__); */
		  /* 		  Warn_And_Exit("\n"); */
		  /* 		} */
		  
		  sum_scale[catg*n_patterns+site] += curr_scaler_pow;
		  
		  do
		    {
		      piecewise_scaler_pow = MIN(curr_scaler_pow,63);
		      curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
		      For(i,tree->mod->ns)
			{
			  p_lk[site*dim1+catg*dim2+i] *= curr_scaler;
			  
			  if(p_lk[site*dim1+catg*dim2+i] > BIG)
			    {
			      PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
			      PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
			      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			      Warn_And_Exit("\n");
			    }
			}
		      curr_scaler_pow -= piecewise_scaler_pow;
		    }
		  while(curr_scaler_pow != 0);
		}
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_P_Lk_AA(t_tree *tree, t_edge *b, t_node *d)
{
/*
           |
	   |<- b
	   |
	   d
          / \
       	 /   \
       	/     \
	n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  phydbl *Pij1,*Pij2;
  int *sum_scale, *sum_scale_v1, *sum_scale_v2;
  int sum_scale_v1_val, sum_scale_v2_val;
  int i;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  short int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;
  phydbl curr_scaler;
  int curr_scaler_pow, piecewise_scaler_pow;
  phydbl p_lk_lim_inf;
  phydbl smallest_p_lk;
  phydbl p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19;
  int *p_lk_loc;

  p_lk_lim_inf = (phydbl)P_LK_LIM_INF;
 
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = NO;
  sum_scale_v1_val = sum_scale_v2_val = 0;
  p1_lk1 = p2_lk2 = .0;

  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;
  
  /* TO DO: Might be worth keeping these directions in memory instead of
     calculating them every time... */
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);

/*   if((dir1 == -1) || (dir2 == -1)) */
/*     { */
/*       PhyML_Printf("\n. d = %d",d->num); */
/*       PhyML_Printf("\n. d->v[0] = %d, d->v[1] = %d, d->v[2] = %d",d->v[0]->num,d->v[1]->num,d->v[2]->num); */
/*       PhyML_Printf("\n. d->b[0] = %d, d->b[1] = %d, d->b[2] = %d",d->b[0]->num,d->b[1]->num,d->b[2]->num); */
/*       PhyML_Printf("\n. d->num = %d dir1 = %d dir2 = %d",d->num,dir1,dir2); */
/*       PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*       Exit(""); */
/*     } */
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_left;
      p_lk_loc = b->p_lk_loc_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_rght;
      p_lk_loc = b->p_lk_loc_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
      sum_scale_v1 = d->b[dir1]->sum_scale_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
      sum_scale_v1 = d->b[dir1]->sum_scale_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
      sum_scale_v2 = d->b[dir2]->sum_scale_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
      sum_scale_v2 = d->b[dir2]->sum_scale_left;
    }
  
  /* Change probability matrices on the two pendant edges */
  Pij1 = d->b[dir1]->Pij_rr;
  Pij2 = d->b[dir2]->Pij_rr;
  
  /* For every site in the alignment */
  For(site,n_patterns)
    {
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = NO;
      
      if(!tree->mod->s_opt->greedy)
	{
	  /* n_v1 and n_v2 are tip nodes */
	  if(n_v1->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
	      if(ambiguity_check_v1 == NO) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
	      
	  if(n_v2->tax)
	    {
	      /* Is the state at this tip ambiguous? */
	      ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
	      if(ambiguity_check_v2 == NO) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dim2,tree);
	    }
	}
      
      if(tree->mod->use_m4mod)
	{
	  ambiguity_check_v1 = YES;
	  ambiguity_check_v2 = YES;
	}

      if(p_lk_loc[site] < site)
	{
	  Copy_P_Lk(p_lk,p_lk_loc[site],site,tree);
	  Copy_Scale(sum_scale,p_lk_loc[site],site,tree);
	}
      else
	{
	  /* For all the rate classes */
	  For(catg,tree->mod->n_catg)
	    {
	      smallest_p_lk  =  BIG;
	      
	      /* For all the state at node d */
	      For(i,tree->mod->ns)
		{
		  p1_lk1 = .0;
		  
		  /* n_v1 is a tip */
		  if((n_v1->tax) && (!tree->mod->s_opt->greedy))
		    {
		      if(ambiguity_check_v1 == NO)
			{
			  /* For the (non-ambiguous) state at node n_v1 */
			  p1_lk1 = Pij1[catg*dim3+i*dim2+state_v1];
			}
		      else
			{
			  /* For all the states at node n_v1 */
			  p0  = Pij1[catg*dim3+i*dim2+ 0] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 0];
			  p1  = Pij1[catg*dim3+i*dim2+ 1] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 1];
			  p2  = Pij1[catg*dim3+i*dim2+ 2] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 2];
			  p3  = Pij1[catg*dim3+i*dim2+ 3] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 3];
			  p4  = Pij1[catg*dim3+i*dim2+ 4] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 4];
			  p5  = Pij1[catg*dim3+i*dim2+ 5] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 5];
			  p6  = Pij1[catg*dim3+i*dim2+ 6] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 6];
			  p7  = Pij1[catg*dim3+i*dim2+ 7] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 7];
			  p8  = Pij1[catg*dim3+i*dim2+ 8] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 8];
			  p9  = Pij1[catg*dim3+i*dim2+ 9] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+ 9];
			  p10 = Pij1[catg*dim3+i*dim2+10] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+10];
			  p11 = Pij1[catg*dim3+i*dim2+11] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+11];
			  p12 = Pij1[catg*dim3+i*dim2+12] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+12];
			  p13 = Pij1[catg*dim3+i*dim2+13] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+13];
			  p14 = Pij1[catg*dim3+i*dim2+14] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+14];
			  p15 = Pij1[catg*dim3+i*dim2+15] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+15];
			  p16 = Pij1[catg*dim3+i*dim2+16] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+16];
			  p17 = Pij1[catg*dim3+i*dim2+17] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+17];
			  p18 = Pij1[catg*dim3+i*dim2+18] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+18];
			  p19 = Pij1[catg*dim3+i*dim2+19] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+19];
			  p1_lk1 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
			}
		    }
		  /* n_v1 is an internal node */
		  else
		    {
		      /* For the states at node n_v1 */
		      p0  = Pij1[catg*dim3+i*dim2+ 0] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 0];
		      p1  = Pij1[catg*dim3+i*dim2+ 1] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 1];
		      p2  = Pij1[catg*dim3+i*dim2+ 2] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 2];
		      p3  = Pij1[catg*dim3+i*dim2+ 3] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 3];
		      p4  = Pij1[catg*dim3+i*dim2+ 4] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 4];
		      p5  = Pij1[catg*dim3+i*dim2+ 5] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 5];
		      p6  = Pij1[catg*dim3+i*dim2+ 6] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 6];
		      p7  = Pij1[catg*dim3+i*dim2+ 7] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 7];
		      p8  = Pij1[catg*dim3+i*dim2+ 8] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 8];
		      p9  = Pij1[catg*dim3+i*dim2+ 9] * (phydbl)p_lk_v1[site*dim1+catg*dim2+ 9];
		      p10 = Pij1[catg*dim3+i*dim2+10] * (phydbl)p_lk_v1[site*dim1+catg*dim2+10];
		      p11 = Pij1[catg*dim3+i*dim2+11] * (phydbl)p_lk_v1[site*dim1+catg*dim2+11];
		      p12 = Pij1[catg*dim3+i*dim2+12] * (phydbl)p_lk_v1[site*dim1+catg*dim2+12];
		      p13 = Pij1[catg*dim3+i*dim2+13] * (phydbl)p_lk_v1[site*dim1+catg*dim2+13];
		      p14 = Pij1[catg*dim3+i*dim2+14] * (phydbl)p_lk_v1[site*dim1+catg*dim2+14];
		      p15 = Pij1[catg*dim3+i*dim2+15] * (phydbl)p_lk_v1[site*dim1+catg*dim2+15];
		      p16 = Pij1[catg*dim3+i*dim2+16] * (phydbl)p_lk_v1[site*dim1+catg*dim2+16];
		      p17 = Pij1[catg*dim3+i*dim2+17] * (phydbl)p_lk_v1[site*dim1+catg*dim2+17];
		      p18 = Pij1[catg*dim3+i*dim2+18] * (phydbl)p_lk_v1[site*dim1+catg*dim2+18];
		      p19 = Pij1[catg*dim3+i*dim2+19] * (phydbl)p_lk_v1[site*dim1+catg*dim2+19];
		      p1_lk1 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
		    }
		  
		  p2_lk2 = .0;
		  
		  /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
		  /* n_v2 is a tip */
		  if((n_v2->tax) && (!tree->mod->s_opt->greedy))
		    {
		      if(ambiguity_check_v2 == NO)
			{
			  /* For the (non-ambiguous) state at node n_v2 */
			  p2_lk2 = Pij2[catg*dim3+i*dim2+state_v2];
			}
		      else
			{
			  /* For all the states at node n_v2 */
			  p0  = Pij2[catg*dim3+i*dim2+ 0] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 0];
			  p1  = Pij2[catg*dim3+i*dim2+ 1] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 1];
			  p2  = Pij2[catg*dim3+i*dim2+ 2] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 2];
			  p3  = Pij2[catg*dim3+i*dim2+ 3] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 3];
			  p4  = Pij2[catg*dim3+i*dim2+ 4] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 4];
			  p5  = Pij2[catg*dim3+i*dim2+ 5] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 5];
			  p6  = Pij2[catg*dim3+i*dim2+ 6] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 6];
			  p7  = Pij2[catg*dim3+i*dim2+ 7] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 7];
			  p8  = Pij2[catg*dim3+i*dim2+ 8] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 8];
			  p9  = Pij2[catg*dim3+i*dim2+ 9] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+ 9];
			  p10 = Pij2[catg*dim3+i*dim2+10] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+10];
			  p11 = Pij2[catg*dim3+i*dim2+11] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+11];
			  p12 = Pij2[catg*dim3+i*dim2+12] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+12];
			  p13 = Pij2[catg*dim3+i*dim2+13] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+13];
			  p14 = Pij2[catg*dim3+i*dim2+14] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+14];
			  p15 = Pij2[catg*dim3+i*dim2+15] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+15];
			  p16 = Pij2[catg*dim3+i*dim2+16] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+16];
			  p17 = Pij2[catg*dim3+i*dim2+17] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+17];
			  p18 = Pij2[catg*dim3+i*dim2+18] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+18];
			  p19 = Pij2[catg*dim3+i*dim2+19] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+19];
			  p2_lk2 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
			  
			}
		    }
		  /* n_v2 is an internal node */
		  else
		    {
		      /* For all the states at node n_v2 */
		      p0  = Pij2[catg*dim3+i*dim2+ 0] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 0];
		      p1  = Pij2[catg*dim3+i*dim2+ 1] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 1];
		      p2  = Pij2[catg*dim3+i*dim2+ 2] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 2];
		      p3  = Pij2[catg*dim3+i*dim2+ 3] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 3];
		      p4  = Pij2[catg*dim3+i*dim2+ 4] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 4];
		      p5  = Pij2[catg*dim3+i*dim2+ 5] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 5];
		      p6  = Pij2[catg*dim3+i*dim2+ 6] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 6];
		      p7  = Pij2[catg*dim3+i*dim2+ 7] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 7];
		      p8  = Pij2[catg*dim3+i*dim2+ 8] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 8];
		      p9  = Pij2[catg*dim3+i*dim2+ 9] * (phydbl)p_lk_v2[site*dim1+catg*dim2+ 9];
		      p10 = Pij2[catg*dim3+i*dim2+10] * (phydbl)p_lk_v2[site*dim1+catg*dim2+10];
		      p11 = Pij2[catg*dim3+i*dim2+11] * (phydbl)p_lk_v2[site*dim1+catg*dim2+11];
		      p12 = Pij2[catg*dim3+i*dim2+12] * (phydbl)p_lk_v2[site*dim1+catg*dim2+12];
		      p13 = Pij2[catg*dim3+i*dim2+13] * (phydbl)p_lk_v2[site*dim1+catg*dim2+13];
		      p14 = Pij2[catg*dim3+i*dim2+14] * (phydbl)p_lk_v2[site*dim1+catg*dim2+14];
		      p15 = Pij2[catg*dim3+i*dim2+15] * (phydbl)p_lk_v2[site*dim1+catg*dim2+15];
		      p16 = Pij2[catg*dim3+i*dim2+16] * (phydbl)p_lk_v2[site*dim1+catg*dim2+16];
		      p17 = Pij2[catg*dim3+i*dim2+17] * (phydbl)p_lk_v2[site*dim1+catg*dim2+17];
		      p18 = Pij2[catg*dim3+i*dim2+18] * (phydbl)p_lk_v2[site*dim1+catg*dim2+18];
		      p19 = Pij2[catg*dim3+i*dim2+19] * (phydbl)p_lk_v2[site*dim1+catg*dim2+19];
		      p2_lk2 = p0+p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p13+p14+p15+p16+p17+p18+p19;
		      
		    }
		  
		  p_lk[site*dim1+catg*dim2+i] = p1_lk1 * p2_lk2;	    
		  
		  if(p_lk[site*dim1+catg*dim2+i] < smallest_p_lk) smallest_p_lk = p_lk[site*dim1+catg*dim2+i] ; 
		}
	      
	      /* Current scaling values at that site */
	      sum_scale_v1_val = (sum_scale_v1)?(sum_scale_v1[catg*n_patterns+site]):(0);
	      sum_scale_v2_val = (sum_scale_v2)?(sum_scale_v2[catg*n_patterns+site]):(0);
	      
	      sum_scale[catg*n_patterns+site] = sum_scale_v1_val + sum_scale_v2_val;
	      
	      /* Scaling */
	      if(smallest_p_lk < p_lk_lim_inf)
		{
		  curr_scaler_pow = (int)(LOG(p_lk_lim_inf)-LOG(smallest_p_lk))/LOG2;
		  curr_scaler     = (phydbl)((unsigned long long)(1) << curr_scaler_pow);
		  
		  /* 	      if(fabs(curr_scaler_pow) > 63 || fabs(curr_scaler_pow) > 63) */
		  /* 		{ */
		  /* 		  PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk); */
		  /* 		  PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow); */
		  /* 		  PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__); */
		  /* 		  Warn_And_Exit("\n"); */
		  /* 		} */
		  
		  sum_scale[catg*n_patterns+site] += curr_scaler_pow;
		  
		  do
		    {
		      piecewise_scaler_pow = MIN(curr_scaler_pow,63);
		      curr_scaler = (phydbl)((unsigned long long)(1) << piecewise_scaler_pow);
		      For(i,tree->mod->ns)
			{
			  p_lk[site*dim1+catg*dim2+i] *= curr_scaler;
			  
			  if(p_lk[site*dim1+catg*dim2+i] > BIG)
			    {
			      PhyML_Printf("\n. p_lk_lim_inf = %G smallest_p_lk = %G",p_lk_lim_inf,smallest_p_lk);
			      PhyML_Printf("\n. curr_scaler_pow = %d",curr_scaler_pow);
			      PhyML_Printf("\n. Err in file %s at line %d.",__FILE__,__LINE__);
			      Warn_And_Exit("\n");
			    }
			}
		      curr_scaler_pow -= piecewise_scaler_pow;
		    }
		  while(curr_scaler_pow != 0);
		}
	    }
	}
    }
}

#endif

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Return_Abs_Lk(t_tree *tree)
{
  Lk(tree);
  return FABS(tree->c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


matrix *ML_Dist(calign *data, model *mod)
{
  int i,j,k,l;
  phydbl init;
  int n_catg;
  phydbl d_max,sum;
  matrix *mat;
  calign *twodata,*tmpdata;
  int state0, state1,len;
  phydbl *F;
  eigen *eigen_struct;

  tmpdata         = (calign *)mCalloc(1,sizeof(calign));
  tmpdata->c_seq  = (align **)mCalloc(2,sizeof(align *));
  tmpdata->b_frq  = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  tmpdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  F               = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl ));
  eigen_struct    = (eigen *)Make_Eigen_Struct(mod->ns);

  tmpdata->n_otu  = 2;

  tmpdata->crunch_len = data->crunch_len;
  tmpdata->init_len   = data->init_len;


  mat = NULL;
  if(mod->io->datatype == NT) mat = (mod->whichmodel < 10)?(K80_dist(data,1E+6)):(JC69_Dist(data,mod));
  else if(mod->io->datatype == AA) mat = JC69_Dist(data,mod);
  else if(mod->io->datatype == GENERIC) mat = JC69_Dist(data,mod);
   

  For(i,mod->n_catg) /* Don't use the discrete gamma distribution */
    {
      mod->gamma_rr->v[i]      = 1.0;
      mod->gamma_r_proba->v[i] = 1.0;
    }

  n_catg = mod->n_catg;
  mod->n_catg = 1;

  For(j,data->n_otu-1)
    {
      tmpdata->c_seq[0]       = data->c_seq[j];
      tmpdata->c_seq[0]->name = data->c_seq[j]->name;
      tmpdata->wght           = data->wght;

      for(k=j+1;k<data->n_otu;k++)
	{

	  tmpdata->c_seq[1]       = data->c_seq[k];
	  tmpdata->c_seq[1]->name = data->c_seq[k]->name;

	  twodata = Compact_Cdata(tmpdata,mod->io);

	  Check_Ambiguities(twodata,mod->io->datatype,mod->io->mod->state_len);
	  
	  Hide_Ambiguities(twodata);
	  
	  init = mat->dist[j][k];
	  
	  if((init > DIST_MAX-SMALL) || (init < .0)) init = 0.1;
	  	  
	  d_max = init;
  
	  For(i,mod->ns*mod->ns) F[i]=.0;
	  len = 0;
 	  For(l,twodata->c_seq[0]->len)
	    {
	      state0 = Assign_State(twodata->c_seq[0]->state+l*mod->io->mod->state_len,mod->io->datatype,mod->io->mod->state_len);
	      state1 = Assign_State(twodata->c_seq[1]->state+l*mod->io->mod->state_len,mod->io->datatype,mod->io->mod->state_len);

	      if((state0 > -1) && (state1 > -1))
		{
		  F[mod->ns*state0+state1] += twodata->wght[l];
		  len += (int)twodata->wght[l];
		}
	    }
	  if(len > .0) {For(i,mod->ns*mod->ns) F[i] /= (phydbl)len;}
	  	  
	  sum = 0.;
	  For(i,mod->ns*mod->ns) sum += F[i];
	  	      
	      
	  if(sum < .001) d_max = -1.;
	  else if((sum > 1. - .001) && (sum < 1. + .001)) Opt_Dist_F(&(d_max),F,mod);
	  else
	    {
	      PhyML_Printf("\n. sum = %f\n",sum);
	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	      Exit("");
	    }
	  	  	  
	  if(d_max >= DIST_MAX) d_max = DIST_MAX;
	  
	  
	  /* Do not correct for dist < BL_MIN, otherwise Fill_Missing_Dist
	   *  will not be called
	   */
	  
	  mat->dist[j][k] = d_max;
	  mat->dist[k][j] = mat->dist[j][k];
	  Free_Cseq(twodata);
	}
    }

  mod->n_catg = n_catg;
    
  
  Free(tmpdata->ambigu);
  Free(tmpdata->b_frq);
  Free(tmpdata->c_seq);
  free(tmpdata);
  Free_Eigen(eigen_struct);
  Free(F);

  return mat;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk_Given_Two_Seq(calign *data, int numseq1, int numseq2, phydbl dist, model *mod, phydbl *loglk)
{
  align *seq1,*seq2;
  phydbl site_lk,log_site_lk;
  int i,j,k,l;
/*   phydbl **p_lk_l,**p_lk_r; */
  phydbl *p_lk_l,*p_lk_r;
  phydbl len;
  int dim1,dim2;

  dim1 = mod->ns;
  dim2 = mod->ns * mod->ns;

  DiscreteGamma(mod->gamma_r_proba->v, mod->gamma_rr->v, mod->alpha->v,
		mod->alpha->v,mod->n_catg,mod->gamma_median);

  seq1 = data->c_seq[numseq1];
  seq2 = data->c_seq[numseq2];


  p_lk_l = (phydbl *)mCalloc(data->c_seq[0]->len * mod->ns,sizeof(phydbl));
  p_lk_r = (phydbl *)mCalloc(data->c_seq[0]->len * mod->ns,sizeof(phydbl));


  For(i,mod->n_catg)
    {
      len = dist*mod->gamma_rr->v[i];
      if(len < mod->l_min) len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,dim2*i,mod->Pij_rr->v);
    }

  if(mod->io->datatype == NT)
    {
      For(i,data->c_seq[0]->len)
	{
	  Init_Tips_At_One_Site_Nucleotides_Float(seq1->state[i],i*mod->ns,p_lk_l);
	  Init_Tips_At_One_Site_Nucleotides_Float(seq2->state[i],i*mod->ns,p_lk_r);
	}
    }
  else if(mod->io->datatype == AA)
    {
      For(i,data->c_seq[0]->len)
	{
	  Init_Tips_At_One_Site_AA_Float(seq1->state[i],i*mod->ns,p_lk_l);
	  Init_Tips_At_One_Site_AA_Float(seq2->state[i],i*mod->ns,p_lk_r);
	}
    }
  else
    {
      PhyML_Printf("\n. Not implemented yet...");
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }


  site_lk = .0;
  *loglk = 0;

  For(i,data->c_seq[0]->len)
    {
      if(data->wght[i])
	{
	  site_lk = log_site_lk = .0;
	  if(!data->ambigu[i])
	    {
	      For(k,mod->ns) {if(p_lk_l[i*mod->ns+k] > .0001) break;}
	      For(l,mod->ns) {if(p_lk_r[i*mod->ns+l] > .0001) break;}
	      For(j,mod->n_catg)
		{
		  site_lk +=
		    mod->gamma_r_proba->v[j] *
		    mod->pi->v[k] *
		    p_lk_l[i*dim1+k] *
		    mod->Pij_rr->v[j*dim2+k*dim1+l] *
		    p_lk_r[i*dim1+l];
		}
	    }
	  else
	    {
	      For(j,mod->n_catg)
		{
		  For(k,mod->ns) /*sort sum terms ? No global effect*/
		    {
		      For(l,mod->ns)
			{
			  site_lk +=
			    mod->gamma_r_proba->v[j] *
			    mod->pi->v[k] *
			    p_lk_l[i*dim1+k] *
			    mod->Pij_rr->v[j*dim2+k*dim1+l] *
			    p_lk_r[i*dim1+l];
			}
		    }
		}
	    }

	  if(site_lk <= .0)
	    {
	      PhyML_Printf("'%c' '%c'\n",seq1->state[i],seq2->state[i]);
	      Exit("\n. Err: site lk <= 0\n");
	    }

	  log_site_lk += (phydbl)LOG(site_lk);

	  *loglk += data->wght[i] * log_site_lk;/* sort sum terms ? No global effect*/
	}
    }

/*   For(i,data->c_seq[0]->len) */
/*     { */
/*       Free(p_lk_l[i]); */
/*       Free(p_lk_r[i]); */
/*     } */

  Free(p_lk_l); Free(p_lk_r);
  return *loglk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Unconstraint_Lk(t_tree *tree)
{
  int i;

  tree->unconstraint_lk = .0;

  For(i,tree->data->crunch_len)
    {
      tree->unconstraint_lk +=
	tree->data->wght[i]*(phydbl)LOG(tree->data->wght[i]);
    }
  tree->unconstraint_lk -=
    tree->data->init_len*(phydbl)LOG(tree->data->init_len);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Make_Tree_4_Lk(t_tree *tree, calign *cdata, int n_site)
{
  int i;

  tree->c_lnL_sorted = (phydbl *)mCalloc(tree->n_pattern, sizeof(phydbl));
  tree->cur_site_lk  = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
  tree->old_site_lk  = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));
  tree->site_lk_cat  = (phydbl *)mCalloc(tree->mod->n_catg,sizeof(phydbl));  
  tree->log_site_lk_cat      = (phydbl **)mCalloc(tree->mod->n_catg,sizeof(phydbl *));
  For(i,tree->mod->n_catg)
    tree->log_site_lk_cat[i] = (phydbl *)mCalloc(cdata->crunch_len,sizeof(phydbl));


  tree->log_lks_aLRT = (phydbl **)mCalloc(3,sizeof(phydbl *));
  For(i,3) tree->log_lks_aLRT[i] = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  For(i,2*tree->n_otu-3)
    {
      Make_Edge_Lk(tree->t_edges[i],tree);
      Make_Edge_NNI(tree->t_edges[i]);
    }

  For(i,2*tree->n_otu-2) Make_Node_Lk(tree->t_nodes[i]);

  if(tree->mod->s_opt->greedy) 
    {
      Init_P_Lk_Tips_Double(tree);
    }
  else 
    {
      Init_P_Lk_Tips_Int(tree);
    }
  
  Init_P_Lk_Loc(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_P_Lk_Tips_Double(t_tree *tree)
{
  int curr_site,i,j,k,dim1,dim2;
  
  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;


  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if (tree->io->datatype == NT)
	    Init_Tips_At_One_Site_Nucleotides_Float(tree->data->c_seq[i]->state[curr_site],
						    curr_site*dim1+0*dim2,
						    tree->t_nodes[i]->b[0]->p_lk_rght);

	  else if(tree->io->datatype == AA)
	    Init_Tips_At_One_Site_AA_Float(tree->data->c_seq[i]->state[curr_site],
					   curr_site*dim1+0*dim2,
					   tree->t_nodes[i]->b[0]->p_lk_rght);

	  else if(tree->io->datatype == GENERIC)
	    Init_Tips_At_One_Site_Generic_Float(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						curr_site*dim1+0*dim2,
						tree->t_nodes[i]->b[0]->p_lk_rght);

	  for(j=1;j<tree->mod->n_catg;j++)
	    {
	      For(k,tree->mod->ns)
		{
		  tree->t_nodes[i]->b[0]->p_lk_rght[curr_site*dim1+j*dim2+k] = 
		    tree->t_nodes[i]->b[0]->p_lk_rght[curr_site*dim1+0*dim2+k];
		}
	    }
	}
    }
  
  if(tree->mod->use_m4mod) 
    {
      M4_Init_P_Lk_Tips_Double(tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_P_Lk_Tips_Int(t_tree *tree)
{
  int curr_site,i,dim1;

  dim1 = tree->mod->ns;

  For(curr_site,tree->data->crunch_len)
    {
      For(i,tree->n_otu)
	{
	  if(tree->io->datatype == NT)
	    Init_Tips_At_One_Site_Nucleotides_Int(tree->data->c_seq[i]->state[curr_site],
						  curr_site*dim1,
						  tree->t_nodes[i]->b[0]->p_lk_tip_r);

	  else if(tree->io->datatype == AA)
	    Init_Tips_At_One_Site_AA_Int(tree->data->c_seq[i]->state[curr_site],
					 curr_site*dim1,					   
					 tree->t_nodes[i]->b[0]->p_lk_tip_r);

	  else if(tree->io->datatype == GENERIC)
	    {
	      Init_Tips_At_One_Site_Generic_Int(tree->data->c_seq[i]->state+curr_site*tree->mod->state_len,
						tree->mod->ns,
						tree->mod->state_len,
						curr_site*dim1,
						tree->t_nodes[i]->b[0]->p_lk_tip_r);
	    }
	}
    }
   if(tree->mod->use_m4mod) 
     {
       M4_Init_P_Lk_Tips_Int(tree);
     }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Update_PMat_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  int i;
  phydbl len;
  phydbl l_min, l_max;

  l_min = tree->mod->l_min;
  l_max = tree->mod->l_max;

  if(tree->mod->gamma_mgf_bl == YES)
    {
      phydbl shape,scale,mean,var;

      mean = b_fcus->gamma_prior_mean;
      var  = b_fcus->gamma_prior_var;
      
      if(mean < l_min) mean = l_min;
      if(mean > l_max) mean = l_max;

      shape = mean*mean/var;
      scale = var/mean;

      For(i,tree->mod->n_catg) 
	{
	  PMat_MGF_Gamma(b_fcus->Pij_rr+tree->mod->ns*tree->mod->ns*i,
			 shape,scale,tree->mod->gamma_rr->v[i],tree->mod);
	}
    }
  else
    {
      len = -1.0;
      
            
      if(tree->mod->log_l == YES) b_fcus->l = EXP(b_fcus->l);
      
      For(i,tree->mod->n_catg)
	{
	  if(b_fcus->has_zero_br_len == YES) 
	    {
	      len = -1.0;
	    }
	  else
	    {
	      len = b_fcus->l*tree->mod->gamma_rr->v[i];	  
	      if(len < l_min)      len = l_min;
	      else if(len > l_max) len = l_max;
	    }

	  len *= tree->mod->br_len_multiplier->v;

	  PMat(len,tree->mod,tree->mod->ns*tree->mod->ns*i,b_fcus->Pij_rr);
	}

      if(tree->mod->log_l == YES) b_fcus->l = LOG(b_fcus->l);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* void Update_P_Lk_On_A_Path(t_node *a, t_node *d, t_edge *b, t_node *target_one_side, t_node *target_other_side, t_tree *tree) */
/* { */


/*   /\* */
/*                 \               / */
/* 	   target\___________b_/ */
/* 		 /  \	\  \   \ */
/* 		/    \	 \  \	\ */

/*     target is the root of the subtree at which we want */
/*     the likelihood to be updated */
/*   *\/ */



/* /\*   PhyML_Printf("Update_p_lk on (%d %d) at %d (target=%d %d)\n", *\/ */
/* /\* 	 b->left->num, *\/ */
/* /\* 	 b->rght->num, *\/ */
/* /\* 	 a->num, *\/ */
/* /\* 	 target_one_side->num, *\/ */
/* /\* 	 target_other_side->num); *\/ */

/*   Update_P_Lk(tree,b,a); */
/*   if((a == target_one_side) && (d == target_other_side))  */
/*     return; */
/*   else */
/*     { */
/*       Update_P_Lk_On_A_Path(d, */
/* 			    d->v[tree->t_dir[d->num][target_other_side->num]], */
/* 			    d->b[tree->t_dir[d->num][target_other_side->num]], */
/* 			    target_one_side, */
/* 			    target_other_side, */
/* 			    tree);  */
/*     } */
/* } */

void Update_P_Lk_Along_A_Path(t_node **path, int path_length, t_tree *tree)
{
  int i,j;

  For(i,path_length-1)
    {
      For(j,3)
	if(path[i]->v[j] == path[i+1])
	  {
	    if(path[i] == path[i]->b[j]->left)
	      {
		Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->left);		    
	      }
	    else if(path[i] == path[i]->b[j]->rght)
	      {
		Update_P_Lk(tree,path[i]->b[j],path[i]->b[j]->rght);
	      }
	    else
	      {
		PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		Exit("");
	      }
	    break;
	  }      
#ifdef DEBUG
      if(j == 3)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
	  Exit("");
	}
#endif
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk_Dist(phydbl *F, phydbl dist, model *mod)
{
  int i,j,k;
  phydbl lnL,len;
  int dim1,dim2;

  if(mod->log_l == YES) dist = EXP(dist);

  For(k,mod->n_catg)
    {
      len = dist*mod->gamma_rr->v[k];
      if(len < mod->l_min)      len = mod->l_min;
      else if(len > mod->l_max) len = mod->l_max;
      PMat(len,mod,mod->ns*mod->ns*k,mod->Pij_rr->v);
    }
  
  dim1 = mod->ns*mod->ns;
  dim2 = mod->ns;
  lnL = .0;

/*   For(i,mod->ns) */
/*     { */
/*       For(j,mod->ns) */
/* 	{ */
/* 	  For(k,mod->n_catg) */
/* 	    { */
/*  	      lnL += */
/* 		F[dim1*k+dim2*i+j] * */
/* 		LOG(mod->pi[i] * mod->Pij_rr[dim1*k+dim2*i+j]); */
/* 	    } */
/* 	} */
/*     } */

  For(i,mod->ns-1)
    {
      for(j=i+1;j<mod->ns;j++)
	{
	  For(k,mod->n_catg)
	    {
 	      lnL +=
		(F[dim1*k+dim2*i+j] + F[dim1*k+dim2*j+i])*
		LOG(mod->pi->v[i] * mod->Pij_rr->v[dim1*k+dim2*i+j]);
	    }
	}
    }

  For(i,mod->ns) For(k,mod->n_catg) lnL += F[dim1*k+dim2*i+i]* LOG(mod->pi->v[i] * mod->Pij_rr->v[dim1*k+dim2*i+i]);

  return lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Update_Lk_At_Given_Edge(t_edge *b_fcus, t_tree *tree)
{
  Update_P_Lk(tree,b_fcus,b_fcus->left);
  Update_P_Lk(tree,b_fcus,b_fcus->rght);
  tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Print_Lk_Given_Edge_Recurr(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  PhyML_Printf("\n___ Edge %3d (left=%3d rght=%3d) lnL=%f",
	 b->num,
	 b->left->num,
	 b->rght->num,
	 Lk_At_Given_Edge(b,tree));

  if(d->tax) return;
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a)
	  Print_Lk_Given_Edge_Recurr(d,d->v[i],d->b[i],tree);
    }
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

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Core of the likelihood calculcation. Assume that the partial likelihoods on both
   sides of t_edge *b are up-to-date. Calculate the log-likelihood at one site.
*/
#ifdef USE_OLD_LK
phydbl Lk_Core(t_edge *b, t_tree *tree)
{
  phydbl log_site_lk, site_lk, site_lk_cat;
  phydbl scale_left, scale_rght;
  phydbl sum;
  int ambiguity_check,state;
  int catg,ns,k,l,site;
  int dim1,dim2,dim3;


  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  log_site_lk     = .0;
  site_lk         = .0;
  site_lk_cat     = .0;
  ambiguity_check = -1;
  state           = -1;
  site            = tree->curr_site;
  ns              = tree->mod->ns;

  tree->old_site_lk[site] = tree->cur_site_lk[site];

  /* Get the vectors of partial likelihood scaling */
  scale_left =
    (b->sum_scale_left)?
    (b->sum_scale_left[site]):
    (0.0);
  
  scale_rght =
    (b->sum_scale_rght)?
    (b->sum_scale_rght[site]):
    (0.0);

  if((b->rght->tax) && (!tree->mod->s_opt->greedy))
    {
      ambiguity_check = tree->data->c_seq[b->rght->num]->is_ambigu[site];
      if(!ambiguity_check) state = Get_State_From_P_Pars(b->p_lk_tip_r,site*dim2,tree);
    }

  if(tree->mod->use_m4mod) ambiguity_check = 1;
  
  /* For all classes of rates */
  For(catg,tree->mod->n_catg)
    {
      site_lk_cat = .0;

      /* b is an external edge */
      if((b->rght->tax) && (!tree->mod->s_opt->greedy))
        {
          /* If the character observed at the tip is NOT ambiguous: ns x 1 terms to consider */
          if(!ambiguity_check)
            {
              sum = .0;
              For(l,ns)
                {
                  sum +=
                    b->Pij_rr[catg*dim3+state*dim2+l] *
                    (phydbl)b->p_lk_left[site*dim1+catg*dim2+l];
                }
              site_lk_cat += sum * tree->mod->pi->v[state];
            }
          /* If the character observed at the tip is ambiguous: ns x ns terms to consider */
          else
            {
              For(k,ns)
                {
                  sum = .0;
                  if(b->p_lk_tip_r[site*dim2+k] > .0)
                    {
                      For(l,ns)
                        {
                          sum +=
                            b->Pij_rr[catg*dim3+k*dim2+l] *
                            (phydbl)b->p_lk_left[site*dim1+catg*dim2+l];
                        }
                      site_lk_cat +=
                        sum *
                        tree->mod->pi->v[k] *
                        (phydbl)b->p_lk_tip_r[site*dim2+k];
                    }
                }
            }
        }
      /* b is an internal edge: ns x ns terms to consider */
      else
        {
          For(k,ns)
            {
              sum = .0;
              if(b->p_lk_rght[site*dim1+catg*dim2+k] > .0)
                {
                  For(l,ns)
                    {
                      sum +=
                        b->Pij_rr[catg*dim3+k*dim2+l] *
                        (phydbl)b->p_lk_left[site*dim1+catg*dim2+l];
                    }
                  site_lk_cat +=
                    sum *
                    tree->mod->pi->v[k] *
                    (phydbl)b->p_lk_rght[site*dim1+catg*dim2+k];
                }
            }
        }

      tree->log_site_lk_cat[catg][site] = site_lk_cat;
      site_lk += site_lk_cat * tree->mod->gamma_r_proba->v[catg];
    }

  /* site_lk may be too small ? */
  if(site_lk < 1.E-300) site_lk = 1.E-300;

  /* The substitution model does not consider invariable sites */
  if(!tree->mod->invar)
    {
      log_site_lk = (phydbl)LOG(site_lk) + (phydbl)scale_left + (phydbl)scale_rght;
    }
  else
    {
      /* The site is invariant */
      if((phydbl)tree->data->invar[site] > -0.5)
        {
          if((scale_left + scale_rght > 0.0) || (scale_left + scale_rght < 0.0))
            site_lk *= (phydbl)EXP(scale_left + scale_rght);
          
          log_site_lk =
            (phydbl)LOG(site_lk*(1.0-tree->mod->pinvar->v) + tree->mod->pinvar->v*tree->mod->pi->v[tree->data->invar[site]]);
          /* LOG(P(D)) = LOG(P(D | subst. rate > 0) * P(subst. rate > 0) + P(D | subst. rate = 0) * P(subst. rate = 0)) */
        }
      else
        {
          log_site_lk = (phydbl)LOG(site_lk*(1.0-tree->mod->pinvar->v)) + (phydbl)scale_left + (phydbl)scale_rght;
          /* Same formula as above with P(D | subs, rate = 0) = 0 */
        }
    }
  
  if(log_site_lk < -BIG) Warn_And_Exit("\nlog_site_lk < -BIG\n");

  For(catg,tree->mod->n_catg)
    tree->log_site_lk_cat[catg][site] =
    LOG(tree->log_site_lk_cat[catg][site]) +
    scale_left +
    scale_rght;
  
  tree->cur_site_lk[site] = log_site_lk;
  /* Multiply log likelihood by the number of times this site pattern is found oin the data */
  tree->c_lnL_sorted[site] = tree->data->wght[site]*log_site_lk;
  tree->c_lnL += tree->c_lnL_sorted[site];
  return log_site_lk;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


/* Update partial likelihood on edge b on the side of b where
   node d lies.
*/
void Update_P_Lk(t_tree *tree, t_edge *b, t_node *d)
{
/*
           |
           |<- b
           |
           d
          / \
         /   \
        /     \
        n_v1   n_v2
*/
  t_node *n_v1, *n_v2;
  phydbl p1_lk1,p2_lk2;
  phydbl *p_lk,*p_lk_v1,*p_lk_v2;
  phydbl *Pij1,*Pij2;
  phydbl max_p_lk;
  phydbl *sum_scale;
  phydbl scale_v1, scale_v2;
  int i,j;
  int catg,site;
  int dir1,dir2;
  int n_patterns;
  int ambiguity_check_v1,ambiguity_check_v2;
  int state_v1,state_v2;
  int dim1, dim2, dim3;

  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;
  dim3 = tree->mod->ns * tree->mod->ns;

  state_v1 = state_v2 = -1;
  ambiguity_check_v1 = ambiguity_check_v2 = -1;
  scale_v1 = scale_v2 = 0.0;
  p1_lk1 = p2_lk2 = .0;


  if(d->tax)
    {
      PhyML_Printf("\n. t_node %d is a leaf...",d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Warn_And_Exit("\n");
    }

  n_patterns = tree->n_pattern;
  
  /* TO DO: Might be worth keeping these directions in memory instead of
     calculating them every time... */
  dir1=dir2=-1;
  For(i,3) if(d->b[i] != b) (dir1<0)?(dir1=i):(dir2=i);

  if((dir1 == -1) || (dir2 == -1))
    {
      PhyML_Printf("\n. d = %d",d->num);
      PhyML_Printf("\n. d->v[0] = %d, d->v[1] = %d, d->v[2] = %d",d->v[0]->num,d->v[1]->num,d->v[2]->num);
      PhyML_Printf("\n. d->b[0] = %d, d->b[1] = %d, d->b[2] = %d",d->b[0]->num,d->b[1]->num,d->b[2]->num);
      PhyML_Printf("\n. d->num = %d dir1 = %d dir2 = %d",d->num,dir1,dir2);
      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
      Exit("");
    }
  
  n_v1 = d->v[dir1];
  n_v2 = d->v[dir2];

  /* Get the partial likelihood vectors on edge b and the two pendant
     edges (i.e., the two other edges connected to d) */
  if(d == b->left)
    {
      p_lk = b->p_lk_left;
      sum_scale = b->sum_scale_left;
    }
  else
    {
      p_lk = b->p_lk_rght;
      sum_scale = b->sum_scale_rght;
    }
      
  if(d == d->b[dir1]->left)
    {
      p_lk_v1 = d->b[dir1]->p_lk_rght;
    }
  else
    {
      p_lk_v1 = d->b[dir1]->p_lk_left;
    }
  
  if(d == d->b[dir2]->left)
    {
      p_lk_v2 = d->b[dir2]->p_lk_rght;
    }
  else
    {
      p_lk_v2 = d->b[dir2]->p_lk_left;
    }
  
  /* Change probability matrices on the two pendant edges */
  Pij1 = d->b[dir1]->Pij_rr;
  Pij2 = d->b[dir2]->Pij_rr;
  
  /* For every site in the alignment */
  For(site,n_patterns)
    {
      /* Current scaling values at that site */
/*       scale_v1 = (sum_scale_v1)?(sum_scale_v1[site]):(0.0); */
/*       scale_v2 = (sum_scale_v2)?(sum_scale_v2[site]):(0.0);  */
/*       sum_scale[site] = scale_v1 + scale_v2; */

      max_p_lk = -BIG;
      state_v1 = state_v2 = -1;
      ambiguity_check_v1 = ambiguity_check_v2 = -1;
      
      if(!tree->mod->s_opt->greedy)
        {
          /* n_v1 and n_v2 are tip nodes */
          if(n_v1->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v1 = tree->data->c_seq[n_v1->num]->is_ambigu[site];
              if(!ambiguity_check_v1) state_v1 = Get_State_From_P_Pars(n_v1->b[0]->p_lk_tip_r,site*dim2,tree);
            }
              
          if(n_v2->tax)
            {
              /* Is the state at this tip ambiguous? */
              ambiguity_check_v2 = tree->data->c_seq[n_v2->num]->is_ambigu[site];
              if(!ambiguity_check_v2) state_v2 = Get_State_From_P_Pars(n_v2->b[0]->p_lk_tip_r,site*dim2,tree);
            }
        }
      
      if(tree->mod->use_m4mod)
        {
          ambiguity_check_v1 = 1;
          ambiguity_check_v2 = 1;
        }

      /* For all the rate classes */
      For(catg,tree->mod->n_catg)
        {
          /* For all the state at node d */
          For(i,tree->mod->ns)
            {
              p1_lk1 = .0;
              
              /* n_v1 is a tip */
              if((n_v1->tax) && (!tree->mod->s_opt->greedy))
                {
                  if(!ambiguity_check_v1)
                    {
                      /* For the (non-ambiguous) state at node n_v1 */
                      p1_lk1 = Pij1[catg*dim3+i*dim2+state_v1];
                    }
                  else
                    {
                      /* For all the states at node n_v1 */
                      For(j,tree->mod->ns)
                        {
                          p1_lk1 += Pij1[catg*dim3+i*dim2+j] * (phydbl)n_v1->b[0]->p_lk_tip_r[site*dim2+j];
                        }
                    }
                }
              /* n_v1 is an internal node */
              else
                {
                  /* For the states at node n_v1 */
                  For(j,tree->mod->ns)
                    {
                      p1_lk1 += Pij1[catg*dim3+i*dim2+j] * (phydbl)p_lk_v1[site*dim1+catg*dim2+j];
                    }
                }
              
              p2_lk2 = .0;
              
              /* We do exactly the same as for node n_v1 but for node n_v2 this time.*/
              /* n_v2 is a tip */
              if((n_v2->tax) && (!tree->mod->s_opt->greedy))
                {
                  if(!ambiguity_check_v2)
                    {
                      /* For the (non-ambiguous) state at node n_v2 */
                      p2_lk2 = Pij2[catg*dim3+i*dim2+state_v2];
                    }
                  else
                    {
                      /* For all the states at node n_v2 */
                      For(j,tree->mod->ns)
                        {
                          p2_lk2 += Pij2[catg*dim3+i*dim2+j] * (phydbl)n_v2->b[0]->p_lk_tip_r[site*dim2+j];
                        }
                    }
                }
              /* n_v2 is an internal node */
              else
                {
                  /* For all the states at node n_v2 */
                  For(j,tree->mod->ns)
                    {
                      p2_lk2 += Pij2[catg*dim3+i*dim2+j] * (phydbl)p_lk_v2[site*dim1+catg*dim2+j];
                    }
                }
              
              p_lk[site*dim1+catg*dim2+i] = (phydbl)(p1_lk1 * p2_lk2);
              
              /* Work out the maximum value of the partial likelihoods at node d */
              if(p_lk[site*dim1+catg*dim2+i] > max_p_lk) max_p_lk = p_lk[site*dim1+catg*dim2+i];
            }
        }
      
      /* Scaling */
      if((max_p_lk < LIM_SCALE_VAL) || (max_p_lk > (1./LIM_SCALE_VAL)))
        {
          /* For each rate class */
          For(catg,tree->mod->n_catg)
            {
              /* For each state at node d */
              For(i,tree->mod->ns)
                {
                  /* Divide the corresponding partial likelihood by the maximum
                     value of the partial likelihoods calculated over all rate
                     classes and states. */
                  p_lk[site*dim1+catg*dim2+i] /= max_p_lk;
                  
/*                if((p_lk[site][catg][i] > BIG) || (p_lk[site][catg][i] < SMALL)) */
/*                  { */
/*                    PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/*                    PhyML_Printf("\n. p_lk[%3d][%2d][%3d] = %G max_p_lk = %G",site,catg,i,p_lk[site][catg][i],max_p_lk); */
/*                    PhyML_Printf("\n. alpha=%f pinv=%f",tree->mod->alpha->v,tree->mod->pinvar->v); */
/*                    For(i,tree->mod->n_catg) PhyML_Printf("\n. rr[%2d] = %G",i,tree->mod->rr[i]); */
/*                    PhyML_Printf("\n. d->b[dir1]->l = %f, d->b[dir2]->l = %f",d->b[dir1]->l,d->b[dir2]->l); */
/*                    PhyML_Printf("\n. d->v[dir1]->num = %d, d->v[dir2]->num = %d",d->v[dir1]->num,d->v[dir2]->num); */
/*                    if(d->v[dir1]->tax) */
/*                      { */
/*                        PhyML_Printf("\n. Character observed at d->v[dir1] = %d",state_v1); */
/*                      } */
/*                    if(d->v[dir2]->tax) */
/*                      { */
/*                        PhyML_Printf("\n. Character observed at d->v[dir2] = %d",state_v2); */
/*                      } */
/*                    Warn_And_Exit("\n. Numerical precision problem ! (send me an e-mail : s.guindon@auckland.ac.nz)\n"); */
/*                  } */
                }
            }
          sum_scale[site] += (phydbl)LOG(max_p_lk);
        }
    }
}

#endif 

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_Subpatt(t_tree *tree)
{
  
  if(tree->n_root)
    {
      Alias_Subpatt_Post(tree->n_root,tree->n_root->v[0],tree);
      Alias_Subpatt_Post(tree->n_root,tree->n_root->v[1],tree);
    }
  else
    {
      Alias_Subpatt_Post(tree->t_nodes[0],tree->t_nodes[0]->v[0],tree);
      /* if(tree->both_sides)  */
      Alias_Subpatt_Pre(tree->t_nodes[0],tree->t_nodes[0]->v[0],tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_One_Subpatt(t_node *a, t_node *d, t_tree *tree)
{
  int i,j;
  int *patt_id_v1, *patt_id_v2, *patt_id_d;
  int *p_lk_loc_d, *p_lk_loc_v1, *p_lk_loc_v2;
  t_node *v1, *v2;
  t_edge *b0, *b1, *b2;
  int curr_patt_id_v1, curr_patt_id_v2;
  int curr_p_lk_loc_v1, curr_p_lk_loc_v2;
  int num_subpatt;

  b0 = b1 = b2 = NULL;

  if(d->tax) 
    {
      patt_id_d  = (d == d->b[0]->left)?(d->b[0]->patt_id_left):(d->b[0]->patt_id_rght);
      p_lk_loc_d = (d == d->b[0]->left)?(d->b[0]->p_lk_loc_left):(d->b[0]->p_lk_loc_rght);

      For(i,tree->n_pattern)
	{
	  For(j,tree->n_pattern)
	    {
	      if(patt_id_d[i] == patt_id_d[j])
		{
		  p_lk_loc_d[i] = j;
		  break;
		}
	      if(j > i)
		{
		  PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	    }
	}
      return;
    }
  else
    {
      v1 = v2 = NULL;
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      if(!v1) { v1=d->v[i]; b1=d->b[i];}
	      else    { v2=d->v[i]; b2=d->b[i];}
	    }
	  else
	    {
	      b0 = d->b[i];
	    }
	}
      

      patt_id_v1  = (v1 == b1->left)?(b1->patt_id_left):(b1->patt_id_rght);
      patt_id_v2  = (v2 == b2->left)?(b2->patt_id_left):(b2->patt_id_rght);
      patt_id_d   = (d  == b0->left)?(b0->patt_id_left):(b0->patt_id_rght);  
      p_lk_loc_d  = (d  == b0->left)?(b0->p_lk_loc_left):(b0->p_lk_loc_rght);
      p_lk_loc_v1 = (v1 == b1->left)?(b1->p_lk_loc_left):(b1->p_lk_loc_rght);
      p_lk_loc_v2 = (v2 == b2->left)?(b2->p_lk_loc_left):(b2->p_lk_loc_rght);
      
      num_subpatt = 0;
      For(i,tree->n_pattern)
	{
	  curr_patt_id_v1  = patt_id_v1[i];
	  curr_patt_id_v2  = patt_id_v2[i];
	  curr_p_lk_loc_v1 = p_lk_loc_v1[i];
	  curr_p_lk_loc_v2 = p_lk_loc_v2[i];

	  p_lk_loc_d[i] = i;

	  if((curr_p_lk_loc_v1 == i) || (curr_p_lk_loc_v2 == i))
	    {
	      p_lk_loc_d[i] = i;
	      patt_id_d[i] = num_subpatt;
	      num_subpatt++;
	    }
	  else
	    if(curr_p_lk_loc_v1 == curr_p_lk_loc_v2)
	      {
		p_lk_loc_d[i] = curr_p_lk_loc_v1;
		patt_id_d[i] = patt_id_d[curr_p_lk_loc_v1];
	      }
	    else
	      {
		for(j=MAX(curr_p_lk_loc_v1,curr_p_lk_loc_v2);j<tree->n_pattern;j++)
		  {
		    if((patt_id_v1[j] == curr_patt_id_v1) && 
		       (patt_id_v2[j] == curr_patt_id_v2))
		      {
			p_lk_loc_d[i] = j;
			
			if(j == i)
			  {
			    patt_id_d[i] = num_subpatt;
			    num_subpatt++;
			  }
			else patt_id_d[i] = patt_id_d[j]; 
			break;
		      }
		    if(j > i)
		      {
			PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
			Warn_And_Exit("");
		      }
		  }
	      }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_Subpatt_Post(t_node *a, t_node *d, t_tree *tree)
{

  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      Alias_Subpatt_Post(d,d->v[i],tree);	      
	    }
	}
      Alias_One_Subpatt(a, d, tree);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Alias_Subpatt_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      Alias_One_Subpatt(d->v[i],d,tree);
	      Alias_Subpatt_Pre(d,d->v[i],tree);	      
	    }      
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Copy_P_Lk(phydbl *p_lk, int site_from, int site_to, t_tree *tree)
{
  int i,j;
  int dim1,dim2;


  dim1 = tree->mod->n_catg * tree->mod->ns;
  dim2 = tree->mod->ns;

/*   PhyML_Printf("\n# %d %d",site_to,site_from); */
      
  For(i,tree->mod->ns) For(j,tree->mod->n_catg)
    {
      p_lk[site_to*dim1+j*dim2+i] = p_lk[site_from*dim1+j*dim2+i];
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Copy_Scale(int *scale, int site_from, int site_to, t_tree *tree)
{
  int i;

  For(i,tree->mod->n_catg) 
    {
      scale[i*tree->n_pattern+site_to] = scale[i*tree->n_pattern+site_from];
/*       PhyML_Printf("\n. %d",scale[i*tree->n_pattern+site_to]); */
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Init_P_Lk_Loc(t_tree *tree)
{
  int i,j;
  t_node *d;
  int *patt_id_d;

  For(i,2*tree->n_otu-3)
    {
      For(j,tree->n_pattern)
	{
	  tree->t_edges[i]->p_lk_loc_left[j] = j;
	  tree->t_edges[i]->p_lk_loc_rght[j] = j;
	}
    }

  For(i,tree->n_otu)
    {
      d = tree->t_nodes[i];
      patt_id_d = (d == d->b[0]->left)?(d->b[0]->patt_id_left):(d->b[0]->patt_id_rght);
      For(j,tree->n_pattern)
	{
	  patt_id_d[j] = (int)tree->data->c_seq[d->num]->state[j];
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Lk_Normal_Approx(t_tree *tree)
{
  phydbl lnL;
  int i;
  int dim;
  phydbl first_order;

  dim = 2*tree->n_otu-3;
  
  lnL = Dnorm_Multi_Given_InvCov_Det(tree->rates->u_cur_l,
				     tree->rates->mean_l,
				     tree->rates->invcov,
				     tree->rates->covdet,
				     2*tree->n_otu-3,YES);

  first_order = 0.0;
  For(i,dim) first_order += (tree->rates->u_cur_l[i] - tree->rates->mean_l[i]) * tree->rates->grad_l[i];
  
  lnL += first_order;
  
  /* printf("\n"); */
  /* For(i,dim) printf("%f\t",tree->rates->u_cur_l[i]); */
  /* printf("\n. Lk=%f %f",lnL,tree->mod->l_min); */


/*   err = NO; */
/*   dim = 2*tree->n_otu-3; */
/*   For(i,dim) */
/*     { */
/*       if((tree->rates->mean_l[i] / tree->mod->l_min < 1.1) &&  */
/* 	 (tree->rates->mean_l[i] / tree->mod->l_min > 0.9)) */
/* 	{	   */
/* 	  lnL -= Log_Dnorm(tree->t_edges[i]->l,tree->rates->mean_l[i],SQRT(tree->rates->cov_l[i*dim+i]),&err); */
/* 	  if(err) */
/* 	    { */
/* 	      PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__); */
/* 	      Warn_And_Exit(""); */
/* 	    } */
/* 	  lambda = 1./SQRT(tree->rates->cov_l[i*dim+i]); */
/* 	  l = (tree->mod->log_l == YES)?(EXP(tree->t_edges[i]->l)):(tree->t_edges[i]->l); */
/* 	  lnL += LOG(Dexp(l,lambda)); */
/* /\* 	  printf("\n. lambda = %f",lambda); *\/ */
/* 	} */
/*     } */

  return(lnL);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Part_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return PART_Lk_At_Given_Edge(b,stree);;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Part_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  return PART_Lk(stree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk(tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Geo_Lk(t_edge *b, t_tree *tree, supert_tree *stree)
{
  TIPO_Lk(tree);
  return tree->geo_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  Lk_At_Given_Edge(b,tree);
  return tree->c_lnL;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Lk_Rates(t_edge *b, t_tree *tree, supert_tree *stree)
{
  RATES_Lk_Rates(tree);
  return tree->rates->c_lnL_rates;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Wrap_Lk_Times(t_edge *b, t_tree *tree, supert_tree *stree)
{
  TIMES_Lk_Times(tree);
  return tree->rates->c_lnL_times;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



phydbl Wrap_Lk_Linreg(t_edge *b, t_tree *tree, supert_tree *stree)
{
  /* RATES_Lk_Linreg(tree); */
  return -1.;
  /* return tree->rates->c_lnL_linreg; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


phydbl Wrap_Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree, supert_tree *stree)
{
  phydbl diff;
  diff = Diff_Lk_Norm_At_Given_Edge(b,tree);
  return(-diff);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Sample_Ancestral_Seq(int mutmap, int fromprior, t_tree *tree)
{
  int rate_cat;
  int i,j,k,l;
  phydbl *probs;
  phydbl sum;
  phydbl u;
  int n_mut;
  FILE *fp;
  phydbl *muttime;
  int *muttype;
  char *s;
  int *ordering;

  probs = (phydbl *)mCalloc(tree->mod->n_catg,sizeof(phydbl));
  
  muttime = (phydbl *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
			   sizeof(phydbl));

  muttype = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
			   sizeof(int));

  ordering = (int *)mCalloc((2*tree->n_otu-3)*5, // At most 5 mutations per branch on average
			    sizeof(int));

  s = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  if(fromprior == YES)
    {
      /* Update P(D_x|X=i) for each state i and node X */
      tree->both_sides = YES;
      Lk(tree);
    }

  For(i,tree->n_pattern)
    {
      /* Sample the rate class from its posterior density */
      For(j,tree->mod->n_catg) 
	{
	  if(fromprior == NO)
	    probs[j] = 
	      EXP(tree->log_site_lk_cat[j][i])*
	      tree->mod->gamma_r_proba->v[j];
	  else
	    probs[j] = tree->mod->gamma_r_proba->v[j];
	}


      /* Scale probas. */
      sum = .0;
      For(j,tree->mod->n_catg) sum += probs[j];	  
      For(j,tree->mod->n_catg) probs[j]/=sum;

      /* CDF */
      for(j=1;j<tree->mod->n_catg;j++) probs[j] += probs[j-1];

      /* Sample rate */
      u = Uni();
      rate_cat = -1;
      For(j,tree->mod->n_catg) 
	if(probs[j] > u) 
	  { 
	    rate_cat = j; 
	    break; 
	  }

      n_mut = 0;
      Sample_Ancestral_Seq_Pre(tree->t_nodes[0],tree->t_nodes[0]->v[0],tree->t_nodes[0]->b[0],
			       i,rate_cat,
			       muttype,muttime,&n_mut,
			       mutmap,fromprior,tree);


      For(j,n_mut) ordering[j] = 0;

      For(j,n_mut-1)
	{
	  for(k=j+1;k<n_mut;k++)
	    {
	      if(muttime[k] > muttime[j]) ordering[k]++;
	      else ordering[j]++;
	    }
	}
      
      strcpy(s,"rosetta.");
      sprintf(s+strlen(s),"%d",i);
      fp = fopen(s,"a");
      PhyML_Fprintf(fp,"\n-1 -1 -1.0 -1");
      
      For(j,n_mut)
	{
	  For(k,n_mut)
	    {
	      if(ordering[k] == j)
		{
		  For(l,tree->data->init_len) if(tree->data->sitepatt[l] == i) break;
		  PhyML_Fprintf(fp,"\n%4d %4d %12f %4d",j,muttype[k],muttime[k],l);
		  /* PhyML_Fprintf(fp,"\n%d",muttype[ordering[j]]); */
		  break;
		}
	    }
	}


      For(j,n_mut)
	{
	  muttype[j] = -2;
	  muttime[j] = +1.;
	}

      fclose(fp);
    }

  Free(s);
  Free(muttype);
  Free(muttime);
  Free(ordering);
  Free(probs);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Sample_Ancestral_Seq_Pre(t_node *a, t_node *d, t_edge *b, 
			      int site, int rate_cat, 
			      int *muttype, phydbl *muttime, int *n_mut, 
			      int mutmap, int fromprior, t_tree *tree)
{

  int i,j;
  int sa,sd;
  phydbl *Pij;
  phydbl *p_lk;
  int dim1, dim2, dim3;
  phydbl sum;
  phydbl u;
  char *c;
  phydbl *probs;

  probs = (phydbl *)mCalloc(tree->mod->ns,sizeof(phydbl));
      
  if(a->tax) 
    c = tree->data->c_seq[a->num]->state+site*tree->mod->state_len;
  else
    c = tree->anc_data->c_seq[a->num-tree->n_otu]->state+site*tree->mod->state_len;
  
  sa = Assign_State(c,
		    tree->mod->io->datatype,
		    tree->mod->state_len);
  
  if(sa == -1) // c is an indel
    {
      For(j,tree->mod->ns) probs[j] = tree->mod->pi->v[j];
      
      for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
      
      u = Uni();
      For(j,tree->mod->ns) 
	if(probs[j] > u) 
	  { 
	    sa = j; 
	    break; 
	  }
    }

  if(d->tax == NO) // Need to sample state at node d
    {

      dim1 = tree->mod->n_catg * tree->mod->ns;
      dim2 = tree->mod->ns;
      dim3 = tree->mod->ns * tree->mod->ns;
      sum  = 0.0;
      
      
      Pij  = b->Pij_rr;
      
      if(d == b->left)
	p_lk = b->p_lk_left;
      else
	p_lk = b->p_lk_rght;
            
      For(j,tree->mod->ns) probs[j] = 0.0;
      
      /* Formula (10) in Nielsen's Mutation Maping paper, e.g. */
      For(j,tree->mod->ns) 
	{
	  if(fromprior == NO)
	    probs[j] = 
	      p_lk[site*dim1+rate_cat*dim2+j] * 
	      Pij[rate_cat*dim3+sa*dim2+j];
	  else
	    probs[j] = Pij[rate_cat*dim3+sa*dim2+j];
	}
      

      /* if(site == 92) */
      /* 	printf("\n. l=%f pr[%f %f %f %f] Pij[%f %f %f %f] plk[%f %f %f %f]", */
      /* 	       b->l * tree->mod->gamma_rr->v[rate_cat], */
      /* 	       probs[0],probs[1],probs[2],probs[3], */
      /* 	       Pij[rate_cat*dim3+sa*dim2+0], */
      /* 	       Pij[rate_cat*dim3+sa*dim2+1], */
      /* 	       Pij[rate_cat*dim3+sa*dim2+2], */
      /* 	       Pij[rate_cat*dim3+sa*dim2+3], */
      /* 	       p_lk[site*dim1+rate_cat*dim2+0], */
      /* 	       p_lk[site*dim1+rate_cat*dim2+1], */
      /* 	       p_lk[site*dim1+rate_cat*dim2+2], */
      /* 	       p_lk[site*dim1+rate_cat*dim2+3]); */


      /* Scale the probabilities */
      sum = 0.0;
      For(j,tree->mod->ns) sum += probs[j];
      For(j,tree->mod->ns) probs[j] /= sum;
      
      /* CDF */
      for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
      
      /* Sample state according to their posterior probas. */
      sd = -1;
      u = Uni();
      For(j,tree->mod->ns) 
	if(probs[j] > u) 
	  { 
	    sd = j; 
	    break; 
	  }

      /* Assign state */
      tree->anc_data->c_seq[d->num-tree->n_otu]->state[site] = Reciproc_Assign_State(sd,tree->io->datatype);
    }
  else
    {
      c = tree->data->c_seq[d->num]->state+site*tree->mod->state_len;
        
      sd = Assign_State(c,
			tree->mod->io->datatype,
			tree->mod->state_len);

      if(sd == -1) // c is an indel
	{
	  For(j,tree->mod->ns) probs[j] = tree->mod->pi->v[j];
	  
	  for(j=1;j<tree->mod->ns;j++) probs[j] += probs[j-1];
	  
	  u = Uni();
	  For(j,tree->mod->ns) 
	    if(probs[j] > u) 
	      { 
		sd = j; 
		break; 
	      }
	}
    }

  /* if(site == 92) */
  /*   { */
  /*     printf("\n. sa=%d (%s,%d) sd=%d (%s,%d) b->l=%f", */
  /* 	     sa,a->tax?a->name:"",a->num,sd,d->tax?d->name:"",d->num,b->l); */
  /*   } */
	       
  if(mutmap == YES) Map_Mutations(a,d,sa,sd,b,site,rate_cat,muttype,muttime,n_mut,tree);
  
  Free(probs);

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Sample_Ancestral_Seq_Pre(d,d->v[i],d->b[i],site,rate_cat,muttype,muttime,n_mut,mutmap,fromprior,tree);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Map_Mutations(t_node *a, t_node *d, int sa, int sd, t_edge *b, int site, int rate_cat, int *muttype, phydbl *muttime, int *n_mut, t_tree *tree)
{
  int i,j;
  phydbl *probs,*all_probs;
  int slast; // Last state visited
  phydbl tlast, td;
  phydbl *Q;
  phydbl u,sum;
  int *mut; // Array of mutations
  int thismut;
  int n_mut_branch;
  phydbl br,cr,ta,gr;
  int n_iter;
  int first_mut;

  all_probs = (phydbl *)mCalloc(tree->mod->ns*tree->mod->ns,sizeof(phydbl));
  mut = (int *)mCalloc(tree->mod->ns*tree->mod->ns,sizeof(int));
  
  // Edge rate
  br = 		      
    (tree->rates->model_log_rates == YES)?
    EXP(tree->rates->br_r[d->num]):
    tree->rates->br_r[d->num];

  // Clock (i.e., overall) rate
  cr = tree->rates->clock_r;
  
  // Age of node a
  ta = tree->rates->nd_t[a->num];
  
  // Site relative rate
  gr = tree->mod->gamma_rr->v[rate_cat];


  // Rate matrix
  Q = tree->mod->qmat->v;
  
  // Length of the 'time' interval considered: product of the branch length by the 
  // current relative rate at that site (set when sampling ancestral sequences)
  td = b->l * gr;
  
  // Matrix of change probabilities
  For(i,tree->mod->ns)
    {
      // We only care about the non-diagonal elements here
      For(j,tree->mod->ns) all_probs[i*tree->mod->ns+j] = -Q[i*tree->mod->ns+j] / Q[i*tree->mod->ns+i];
      
      // Set the diagonal to 0
      all_probs[i*tree->mod->ns+i] = 0.0;
      
      // \sum_{j != i} -q_{ij}/q_{ii}
      sum = 0;
      For(j,tree->mod->ns) sum += all_probs[i*tree->mod->ns+j];
      
      // Diagonal: 1 - \sum_{j != i} -q_{ij}/q_{ii}
      all_probs[i*tree->mod->ns+i] = 1.-sum;
      
      // Get the cumulative probas
      for(j=1;j<tree->mod->ns;j++) all_probs[i*tree->mod->ns+j] += all_probs[i*tree->mod->ns+j-1];
    }
  
  For(i,tree->mod->ns*tree->mod->ns) mut[i] = 0;
  tlast = .0;
  slast = sa;
  probs = NULL;
  n_mut_branch = 0;
  n_iter = 0;
  first_mut = YES;

  do
    {
      if((sa != sd) && (first_mut == YES)) // ancestral and descendant states are distinct
	{
	  // Sample a time for the first mutation conditional on at least one mutation
	  // occurring (see formula A2 in Nielsen's Genetics paper (2001))
	  u = Uni();
	  tlast = -LOG(1. - u*(1.-EXP(Q[sa*tree->mod->ns+sa]*td)))/-Q[sa*tree->mod->ns+sa];
	}
      else
	{
	  // Sample a time for the next mutation
	  tlast = tlast + Rexp(-Q[slast*tree->mod->ns+slast]);
	}

      // Select the appropriate vector of change probabilities
      probs = all_probs+slast*tree->mod->ns;
      
      /* printf("\n. slast=%2d sd=%2d tlast=%12G td=%12G p=%12f rcat=%12f site=%4d", */
      /* 	 slast,sd,tlast,td,-Q[slast*tree->mod->ns+slast], */
      /* 	 tree->mod->gamma_rr->v[rate_cat],site); */
      
      // The time for the next mutation does not exceed the length
      // of the time interval -> sample a new mutation event
      if(tlast < td)
	{
	  first_mut = NO;

	  n_mut_branch++;

	  u = Uni();
	  For(i,tree->mod->ns)
	    if(probs[i] > u) 
	      {
		// Record mutation type
		mut[slast*tree->mod->ns+i]++;
		
		// Record mutation type in the site mutation array
		thismut = MIN(i,slast) * tree->mod->ns + MAX(i,slast) - (MIN(i,slast)+1+(int)POW(MIN(i,slast)+1,2))/2;
		muttype[(*n_mut)+n_mut_branch-1] = thismut;

		if((thismut > 5) || (thismut < 0)) 
		  {
		    PhyML_Printf("\n. thismut = %d",thismut);
		    PhyML_Printf("\n. Err in file %s at line %d\n\n",__FILE__,__LINE__);
		    Warn_And_Exit("");
		  }

		// Record time of mutation
		muttime[(*n_mut)+n_mut_branch-1] = ta + br*cr*gr;
		  
		// Update the last state
		slast = i; 		    
		break; 
	      }
	}
      else
	{
	  if(slast == sd) break;
	  else
	    {
	      // Restart from the beginning
	      For(i,tree->mod->ns*tree->mod->ns) mut[i] = 0;
	      For(i,n_mut_branch) muttype[(*n_mut)+n_mut_branch-1] = -2;
	      For(i,n_mut_branch) muttime[(*n_mut)+n_mut_branch-1] = +1.;
	      tlast = 0.0;
	      slast = sa;
	      n_mut_branch = 0;
	      first_mut = YES;
	      n_iter++;
	    }
	}
    }
  while(1);
  
  (*n_mut) += n_mut_branch;


  For(i,tree->mod->ns)
    {
      for(j=i+1;j<tree->mod->ns;j++)
	{
	  if(mut[i*tree->mod->ns+j] + mut[j*tree->mod->ns+i] > 0)
	    {
	      thismut = MIN(i,j) * tree->mod->ns + MAX(i,j) - (MIN(i,j)+1+(int)POW(MIN(i,j)+1,2))/2;
	      tree->mutmap[thismut*(tree->n_pattern)*(2*tree->n_otu-3) + b->num*(tree->n_pattern) + site]++;
	      /* if(site == 92) */
	      /* 	{ */
	      /* 	  printf("\nx sa=%d sd=%d mut=%d",sa,sd,thismut); */
	      /* 	} */
	    }
	}
    }
  
  Free(all_probs);
  Free(mut);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

phydbl Lk_LastFirst(t_tree *tree)
{
  phydbl c_lnL;
  t_tree *buff;
  phydbl sum;

  c_lnL = 0.0;

  // Normalise the mixture model class proba
  buff = tree;
  sum  = 0.0;
  do
    {      
      sum += buff->mod->unscaled_prob;
      buff   = buff->prevtree;      
    }while(buff);


  buff = tree;
  do
    {
      buff->mod->prob = buff->mod->unscaled_prob/sum;
      buff   = buff->prevtree;      
    }while(buff);



  // Calculate the likelihood
  buff = tree;
  do
    {
      c_lnL += buff->c_lnL_mixt * buff->mod->prob;
      /* printf("\n. lk tree %d: %f %f",buff->tree_num,buff->c_lnL_mixt,buff->mod->kappa->v); */
      buff   = buff->prevtree;
   }while(buff);

  buff = tree;
  do
    {
      buff->c_lnL = c_lnL; // Overrides the values of c_lnL;
      buff   = buff->prevtree;
    }while(buff);

  return(c_lnL);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

