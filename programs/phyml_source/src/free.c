/*

PHYML :  a program that  computes maximum likelihood  phyLOGenies from
DNA or AA homoLOGous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "free.h"


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_All_Nodes_Light(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-1) 
    {
      if(tree->t_nodes[i])
	{
	  Free_Node(tree->t_nodes[i]);
	}
    }
  Free(tree->t_nodes);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_All_Edges_Light(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-2) 
    if(tree->t_edges[i])
      {
	Free_Edge(tree->t_edges[i]);
      }
  Free(tree->t_edges);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Mat(matrix *mat)
{
  int i;

  For(i,mat->n_otu)
    {
      Free(mat->P[i]);
      Free(mat->Q[i]);
      Free(mat->dist[i]);
      Free(mat->name[i]);
    }

  Free(mat->P);
  Free(mat->Q);
  Free(mat->dist);
  Free(mat->name);
  Free(mat->tip_node);
      
  Free(mat->on_off);
  Free(mat);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Partial_Lk(phydbl *p_lk, int len, int n_catg)
{
  Free(p_lk);

/*   int i,j; */
/*   For(i,len) */
/*     { */
/*       For(j,n_catg) Free((*p_lk)[i][j]); */
/*       Free((*p_lk)[i]); */
/*     } */
/*   Free((*p_lk)); */
/*   (*p_lk) = NULL; */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Tree(t_tree *tree)
{
  Free(tree->t_dir);
  if(tree->short_l) Free(tree->short_l);
  if(tree->mutmap) Free(tree->mutmap);
  Free_Bip(tree);
  Free(tree->curr_path);
  Free_All_Edges_Light(tree);
  Free_All_Nodes_Light(tree);
  Free(tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Bip(t_tree *tree)
{
  int i,j;

  if(tree->has_bip)
    {
      For(i,2*tree->n_otu-2)
	{
	  Free(tree->t_nodes[i]->bip_size);
	  For(j,3) Free(tree->t_nodes[i]->bip_node[j]);
	  Free(tree->t_nodes[i]->bip_node);
	}
    }
  tree->has_bip = NO;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Edge_Labels(t_edge *b)
{
  int i;
  
  if(b->labels)
    {
      For(i,b->n_labels-(b->n_labels%BLOCK_LABELS)+BLOCK_LABELS) Free(b->labels[i]);
      Free(b->labels);
      b->labels = NULL;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Edge(t_edge *b)
{
  Free_Edge_Labels(b);
  Free(b);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Node(t_node *n)
{
  Free(n->b);
  Free(n->v);
  Free(n->l);
  Free(n->score);
  Free(n->s_ingrp);
  Free(n->s_outgrp);
  if(n->ori_name) { Free(n->ori_name); n->ori_name = NULL; }
  /* if(n->name)     { Free(n->name);     n->name     = NULL; }  */
  /* Don't do that: see Copy_Tax_Names_To_Tip_Labels       
     tree->t_nodes[i]->ori_name = tree->t_nodes[i]->name; */
  Free(n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Cseq(calign *data)
{
  int i;
  
  Free(data->invar);
  Free(data->wght);
  Free(data->ambigu);
  Free(data->b_frq);
  Free(data->sitepatt);
  For(i,data->n_otu)
    {
      Free(data->c_seq[i]->name);
      if(data->c_seq[i]->state) 
	{
	  Free(data->c_seq[i]->state);
	  if(data->c_seq[i]->is_ambigu) Free(data->c_seq[i]->is_ambigu);
	}
      Free(data->c_seq[i]);
    }
  Free(data->c_seq);
  Free(data);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Seq(align **d, int n_otu)
{
  int i;
  For(i,n_otu)
    {
      Free(d[i]->name);
      Free(d[i]->state);
      if(d[i]->is_ambigu) Free(d[i]->is_ambigu);
      Free(d[i]);
    }
  Free(d);
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_All(align **d, calign *cdata, t_tree *tree)
{
  Free_Cseq(cdata);
  Free_Seq(d,tree->n_otu);
  Free_Tree(tree);
}      

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_SubTree(t_edge *b_fcus, t_node *a, t_node *d, t_tree *tree)
{
  int i;

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      Free_SubTree(d->b[i],d,d->v[i],tree);
	      Free_Edge(d->b[i]);
	      Free_Node(d->v[i]);
	    }
	}
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Tree_Ins_Tar(t_tree *tree)
{
  return;
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Tree_Pars(t_tree *tree)
{
  int i;
  
  Free(tree->step_mat);
  Free(tree->site_pars);
  For(i,2*tree->n_otu-3)
    Free_Edge_Pars(tree->t_edges[i],tree);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Edge_Pars(t_edge *b, t_tree *tree)
{
/*   int i; */

  Free(b->pars_l);
  Free(b->pars_r);
  
/*   For(i,tree->data->crunch_len)  */
/*     { */
/*       Free(b->p_pars_l[i]); */
/*       Free(b->p_pars_r[i]); */
/*     } */
  
  Free(b->ui_l);
  Free(b->ui_r);
  Free(b->p_pars_l);
  Free(b->p_pars_r);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Tree_Lk(t_tree *tree)
{
  int i;
  t_edge *b;
  t_node *n;

  b = NULL;
  n = NULL;

  For(i,3) Free(tree->log_lks_aLRT[i]);
  Free(tree->log_lks_aLRT);

  Free(tree->c_lnL_sorted);
  Free(tree->cur_site_lk);
  Free(tree->old_site_lk);
  Free(tree->site_lk_cat);

  For(i,tree->mod->n_catg) Free(tree->log_site_lk_cat[i]);
  Free(tree->log_site_lk_cat);
				
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];      
      Free_Edge_Lk(tree,b);
    }
}


//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void Free_Node_Lk(t_node *n)
{
/*   Free(n->n_ex_nodes); */
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Edge_Lk(t_tree *tree, t_edge *b)
{
  Free(b->nni);

  Free(b->div_post_pred_left);
  Free(b->div_post_pred_rght);

  if(b->p_lk_left)
    {
      Free(b->p_lk_left);
      if(b->sum_scale_left) Free(b->sum_scale_left);
    }

  if(b->p_lk_tip_l) Free(b->p_lk_tip_l);


  if(b->p_lk_rght)
    {
      Free(b->p_lk_rght);
      if(b->sum_scale_rght) Free(b->sum_scale_rght);
    }
  
  if(b->p_lk_tip_r) Free(b->p_lk_tip_r);

  Free(b->sum_scale_left_cat);
  Free(b->sum_scale_rght_cat);

  Free(b->patt_id_left);
  Free(b->patt_id_rght);
  Free(b->p_lk_loc_left);
  Free(b->p_lk_loc_rght);

  Free(b->Pij_rr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Model_Complete(model *mod)
{
  Free(mod->pi->v);
  Free(mod->pi_unscaled->v);
  Free(mod->gamma_r_proba->v);
  Free(mod->gamma_r_proba_unscaled->v);
  Free(mod->gamma_rr->v);
  Free(mod->gamma_rr_unscaled->v);
  Free(mod->Pij_rr->v);
  Free(mod->qmat->v);
  Free(mod->qmat_buff->v);
  Free_Eigen(mod->eigen);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Model_Basic(model *mod)
{
  Free(mod->modelname);
  Free(mod->custom_mod_string);
  Free(mod->user_b_freq->v);
  Free(mod->user_b_freq);
  Free(mod->kappa);
  Free(mod->lambda);
  Free(mod->alpha);
  Free(mod->pinvar);
  Free(mod->alpha_old);
  Free(mod->kappa_old);
  Free(mod->lambda_old);
  Free(mod->pinvar_old);
  Free(mod->rr);
  Free(mod->rr_val);
  Free(mod->rr_num);
  Free(mod->n_rr_per_cat);
  Free(mod->Pij_rr);
  Free(mod->mr);
  Free(mod->qmat);
  Free(mod->qmat_buff);
  Free(mod->pi);
  Free(mod->pi_unscaled);
  Free(mod->gamma_r_proba);
  Free(mod->gamma_r_proba_unscaled);
  Free(mod->gamma_rr);
  Free(mod->gamma_rr_unscaled);
  Free(mod->br_len_multiplier);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Custom_Model(model *mod)
{
  if(mod->rr->v)
    {
      Free(mod->rr_num->v);
      Free(mod->rr->v);
      Free(mod->rr_val->v);
      Free(mod->n_rr_per_cat->v);
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Model(model *mod)
{
  Free_Custom_Model(mod);
  Free_Model_Complete(mod);
  Free_Model_Basic(mod);
  M4_Free_M4_Model(mod->m4mod);
  Free(mod);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free(void *p)
{
  free(p);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Input(option *io)
{
  int i;
  RATES_Free_Rates(io->rates);
  MCMC_Free_MCMC(io->mcmc);
  Free(io->in_align_file);
  Free(io->in_tree_file);
  Free(io->in_constraint_tree_file);
  Free(io->out_tree_file);
  Free(io->out_trees_file);
  Free(io->out_boot_tree_file);
  Free(io->out_boot_stats_file);
  Free(io->out_stats_file);
  Free(io->out_lk_file); 
  Free(io->out_ps_file);
  Free(io->out_trace_file);
  Free(io->aa_rate_mat_file);
  Free(io->nt_or_cd);
  Free(io->run_id_string);
  Free(io->clade_list_file);
  For(i,T_MAX_ALPHABET) Free(io->alphabet[i]);
  Free(io->alphabet);
  if(io->short_tax_names)
    {
      For(i,io->size_tax_names) 
	{
	  Free(io->short_tax_names[i]);
	  Free(io->long_tax_names[i]);
	}
      Free(io->long_tax_names);
      Free(io->short_tax_names);
    }
  Free_Tree_List(io->treelist);
  if(io->lon) Free(io->lon);
  if(io->lat) Free(io->lat);
  Free(io);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Tree_List(t_treelist *list)
{
  Free(list->tree);
  Free(list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_St(supert_tree *st)
{
  int i;

  For(i,2*st->tree->n_otu-3) 
    Free(st->tree->t_edges[i]->nni);

  For(i,st->n_part) Free(st->match_st_node_in_gt[i]);

  Free(st->match_st_node_in_gt);

  Free_Tree(st->tree);
  
  Free(st);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Eigen(eigen *eigen_struct)
{
  Free(eigen_struct->space_int);
  Free(eigen_struct->space);
  Free(eigen_struct->e_val);
  Free(eigen_struct->e_val_im);
  Free(eigen_struct->r_e_vect);
  Free(eigen_struct->r_e_vect_im);
  Free(eigen_struct->l_e_vect);
  Free(eigen_struct->q);
  Free(eigen_struct);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_One_Spr(spr *this_spr)
{
  Free(this_spr->path);
  Free(this_spr);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Spr_List(t_tree *tree)
{
  int i;

  For(i,tree->size_spr_list+1) Free_One_Spr(tree->spr_list[i]);
  Free(tree->spr_list);

}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////



void Free_Triplet(triplet *t)
{
  int i,j,k;

  Free(t->F_bc);
  Free(t->F_cd);
  Free(t->F_bd);
  Free(t->pi_bc);
  Free(t->pi_cd);
  Free(t->pi_bd);

  For(k,t->mod->n_catg) 
    {
      For(i,t->size) 
	{
	  For(j,t->size) Free(t->core[k][i][j]);  
	  Free(t->core[k][i]);
	}
      Free(t->core[k]);	  
    }
  Free(t->core);

  For(i,t->size) 
    {
      For(j,t->size) Free(t->p_one_site[i][j]);  
      Free(t->p_one_site[i]);
    }
  Free(t->p_one_site);

  For(i,t->size) 
    {
      For(j,t->size) Free(t->sum_p_one_site[i][j]);  
      Free(t->sum_p_one_site[i]);
    }
  Free(t->sum_p_one_site);

  Free_Eigen(t->eigen_struct);
  
  Free(t);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Actual_CSeq(calign *data)
{
  int i;
  For(i,data->n_otu)
    {
      Free(data->c_seq[i]->state);
      data->c_seq[i]->state = NULL;
    }
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Prefix_Tree(pnode *n, int size)
{
  int i;
  
  For(i,size)
    {
      if(n->next[i])
	{
	  Free_Prefix_Tree(n->next[i],size);
	}
    }
  Free_Pnode(n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Pnode(pnode *n)
{
  Free(n->next);
  Free(n);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Optimiz(optimiz *s_opt)
{
  Free(s_opt);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

void Free_Nexus(option *io)
{
  int i,j;
  
  For(i,N_MAX_NEX_COM)
    {
      For(j,io->nex_com_list[i]->nparm) Free_Nexus_Parm(io->nex_com_list[i]->parm[j]);
      Free(io->nex_com_list[i]->parm);
      Free(io->nex_com_list[i]->name);
      Free(io->nex_com_list[i]);      
    }
  Free(io->nex_com_list);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Nexus_Com(nexcom **com)
{
  int i;

  For(i,N_MAX_NEX_COM)
    {
      Free(com[i]->parm);
      Free(com[i]->name);
      Free(com[i]);
    }
  Free(com);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void Free_Nexus_Parm(nexparm *parm)
{
  Free(parm->value);
  Free(parm->name);
  Free(parm);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

