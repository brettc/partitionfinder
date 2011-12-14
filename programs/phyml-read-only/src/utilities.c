/*

PhyML:  a program that  computes maximum likelihood phylogenies from
DNA or AA homologous sequences.

Copyright (C) Stephane Guindon. Oct 2003 onward.

All parts of the source except where indicated are distributed under
the GNU public licence. See http://www.opensource.org for details.

*/

#include "utilities.h"


/*********************************************************/

/* void Make_All_Edges_Light(t_node *a, t_node *d, int *curr_num_edge) */
/* { */
/*   int i; */

/*   Make_Edge_Light(a,d,*curr_num_edge); */
/*   (*curr_num_edge)++; */
/*   if(d->tax) return; */
/*   else */
/*     { */
/*       For(i,3) */
/* 	{ */
/* 	  if(d->v[i] != a) */
/* 	    Make_All_Edges_Light(d,d->v[i],curr_num_edge); */
/* 	} */
/*     } */
/* } */

/*********************************************************/

void Make_All_Edges_Lk(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  For(i,3) if((a->v[i]) && (a->v[i] == d)) Make_Edge_Lk(a->b[i],tree);
  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Make_All_Edges_Lk(d,d->v[i],tree);
	}
    }
}

/*********************************************************/

t_tree *Read_Tree(char *s_tree)
{
  char **subs;
  int i,n_ext,n_int,n_otu;
  t_tree *tree;
  int degree;
  t_node *root_node;

  n_int = n_ext = 0;
  
  n_otu=0;
  For(i,(int)strlen(s_tree)) if(s_tree[i] == ',') n_otu++;
  n_otu+=1;

  tree = Make_Tree_From_Scratch(n_otu,NULL);

  subs = Sub_Trees(s_tree,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree == 2) 
    {
      /* Unroot_Tree(subs); */
      /* degree = 3; */
      /* root_node = tree->noeud[n_otu]; */
      root_node      = tree->noeud[2*n_otu-2];
      root_node->num = 2*n_otu-2;
      tree->n_root   = root_node;
      n_int         -= 1;
    }
  else
    {
      root_node      = tree->noeud[n_otu];
      root_node->num = n_otu;
      tree->n_root   = NULL;
   }
  
  root_node->tax = 0;

  tree->has_branch_lengths = 0;
  tree->num_curr_branch_available = 0;
  For(i,degree) R_rtree(s_tree,subs[i],root_node,tree,&n_int,&n_ext);
  
  if(tree->n_root)
    {
      tree->e_root = tree->t_edges[tree->num_curr_branch_available];
            
      For(i,3) if(tree->n_root->v[0]->v[i] == tree->n_root) { tree->n_root->v[0]->v[i] = tree->n_root->v[1]; break; }
      For(i,3) if(tree->n_root->v[1]->v[i] == tree->n_root) { tree->n_root->v[1]->v[i] = tree->n_root->v[0]; break; }

      Connect_One_Edge_To_Two_Nodes(tree->n_root->v[0],
      				    tree->n_root->v[1],
      				    tree->e_root,
      				    tree);

      tree->e_root->l = tree->n_root->l[0] + tree->n_root->l[1];
      if(tree->e_root->l > 0.0)
	tree->n_root_pos = tree->n_root->l[0] / tree->e_root->l;
      else
	tree->n_root_pos = .5;
    }
  
  For(i,NODE_DEG_MAX) Free(subs[i]);
  Free(subs);
  return tree;
}

/*********************************************************/


/*********************************************************/
/* 'a' in t_node a stands for ancestor. 'd' stands for descendant */ 
void R_rtree(char *s_tree_a, char *s_tree_d, t_node *a, t_tree *tree, int *n_int, int *n_ext)
{
  int i;
  t_node *d;
  int n_otu = tree->n_otu;

  if(strstr(s_tree_a," ")) 
    {
      PhyML_Printf("\n. [%s]",s_tree_a);
      Warn_And_Exit("\n. Err: the tree must not contain a ' ' character\n");
    }

  if(s_tree_d[0] == '(')
    {
      char **subs;
      int degree;

      (*n_int)+=1;
      d      = tree->noeud[n_otu+*n_int];
      d->num = n_otu+*n_int;
      d->tax = 0;
      
      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]);
      Read_Branch_Length(s_tree_d,s_tree_a,tree);
      
      For(i,3)
	{
	  if(!a->v[i])
	    {
	      a->v[i]=d;
	      d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l;
	      a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	      break;
	    }
	}
      d->v[0]=a;

      if(a != tree->n_root)
	{
	  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
	  tree->num_curr_branch_available++;
	}
      
      subs=Sub_Trees(s_tree_d,&degree);
      Clean_Multifurcation(subs,degree,2);
      R_rtree(s_tree_d,subs[0],d,tree,n_int,n_ext);
      R_rtree(s_tree_d,subs[1],d,tree,n_int,n_ext);
      For(i,NODE_DEG_MAX) Free(subs[i]);
      Free(subs);
    }

  else
    {
      int i;

      d      = tree->noeud[*n_ext];
      d->tax = 1;

      Read_Node_Name(d,s_tree_d,tree);
      Read_Branch_Label(s_tree_d,s_tree_a,tree->t_edges[tree->num_curr_branch_available]); 
      Read_Branch_Length(s_tree_d,s_tree_a,tree);
      
      For(i,3)
	{
	 if(!a->v[i])
	   {
	     a->v[i]=d;
	     d->l[0]=tree->t_edges[tree->num_curr_branch_available]->l;
	     a->l[i]=tree->t_edges[tree->num_curr_branch_available]->l;
	     break;
	   }
	}
      d->v[0]=a;

      if(a != tree->n_root)
	{
	  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
	  tree->num_curr_branch_available++;
	}

      d->num=*n_ext;
      (*n_ext)+=1;
    }
}

/*********************************************************/

void Read_Branch_Label(char *s_d, char *s_a, t_edge *b)
{
  char *sub_tp;
  char *p;
  int i,pos;

  sub_tp = (char *)mCalloc(3+(int)strlen(s_d)+1,sizeof(char));

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,"#");
  p = strstr(s_a,sub_tp);
  if(!p)
    {
      sub_tp[0] = ',';
      sub_tp[1] = '\0';
      strcat(sub_tp,s_d);
      strcat(sub_tp,"#");
      p = strstr(s_a,sub_tp);
    }

  i = 0;
  b->n_labels = 0;
  if(p)
    {
      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      b->n_labels++;
      
      pos = 0;
      do 
	{
	  b->labels[b->n_labels-1][pos] = p[i+strlen(s_d)+1];
	  i++;
	  pos++;
	  if(p[i+strlen(s_d)+1] == '#') 
	    { 
	      b->labels[b->n_labels-1][pos] = '\0';
	      b->n_labels++;
	      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
	      i++;
	      pos=0;
	    }
	}
      while((p[i+strlen(s_d)+1] != ':') && 
	    (p[i+strlen(s_d)+1] != ',') && 
	    (p[i+strlen(s_d)+1] != '('));

      b->labels[b->n_labels-1][pos] = '\0';
    }

  if(p)
    {
      if(b->n_labels == 1)
	PhyML_Printf("\n. Found label '%s' on t_edge %3d.",b->labels[0],b->num);
      else
	{
	  PhyML_Printf("\n. Found labels ");
	  For(i,b->n_labels) PhyML_Printf("'%s' ",b->labels[i]);
	  PhyML_Printf("on t_edge %3d.",b->num);
	}
    }

  Free(sub_tp);
}

/*********************************************************/

void Read_Branch_Length(char *s_d, char *s_a, t_tree *tree)
{
  char *sub_tp;
  char *p;
  t_edge *b;
  int i;

  b = tree->t_edges[tree->num_curr_branch_available];

/*   sub_tp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
  sub_tp = (char *)mCalloc(4+strlen(s_d)+1,sizeof(char));

  For(i,b->n_labels)
    {
      strcat(s_d,"#");
      strcat(s_d,b->labels[i]);
    }

  sub_tp[0] = '(';
  sub_tp[1] = '\0';
  strcat(sub_tp,s_d);
  strcat(sub_tp,":");
  p = strstr(s_a,sub_tp);
  if(!p)
    {
      sub_tp[0] = ',';
      sub_tp[1] = '\0';
      strcat(sub_tp,s_d);
      strcat(sub_tp,":");
      p = strstr(s_a,sub_tp);
    }

  if(p)
    {
      b->l = atof((char *)p+(int)strlen(sub_tp));
      tree->has_branch_lengths = 1;
    }      
  Free(sub_tp);
}

/*********************************************************/

void Read_Node_Name(t_node *d, char *s_tree_d, t_tree *tree)
{
  int i;
  

  if(!tree->t_edges[tree->num_curr_branch_available]->n_labels)
    {
      d->name = (char *)mCalloc(strlen(s_tree_d)+1,sizeof(char ));
      strcpy(d->name,s_tree_d);
    }
  else
    {
      i = 0;
      do
	{
	  d->name = (char *)realloc(d->name,(i+1)*sizeof(char ));
	  d->name[i] = s_tree_d[i];
	  i++;
	}
      while(s_tree_d[i] != '#');
      d->name[i] = '\0';
    }
  d->ori_name = d->name;

}
/*********************************************************/

void Unroot_Tree(char **subtrees)
{
  char **tmp_sub;
  int degree,i,j;

  PhyML_Printf("\n. Removing the root...\n");
  
  tmp_sub = Sub_Trees(subtrees[0],&degree);
  if(degree >= 2)
    {
      strcpy(subtrees[2],subtrees[1]);
      Clean_Multifurcation(tmp_sub,degree,2);
      For(j,2) strcpy(subtrees[j],tmp_sub[j]);
    }
  else
    {
      tmp_sub = Sub_Trees(subtrees[1],&degree);
      strcpy(subtrees[2],subtrees[0]);
      Clean_Multifurcation(tmp_sub,degree,2);
      For(j,2) strcpy(subtrees[j],tmp_sub[j]);
    }

  For(i,degree) Free(tmp_sub[i]);
  Free(tmp_sub);
}

/*********************************************************/

void Clean_Multifurcation(char **subtrees, int current_deg, int end_deg)
{

  if(current_deg <= end_deg) return;
  else
    {
      char *s_tmp;
      int i;

/*       s_tmp = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
      s_tmp = (char *)mCalloc(6+
			      (int)strlen(subtrees[0])+1+
			      (int)strlen(subtrees[1])+1,
			      sizeof(char));

      strcat(s_tmp,"(\0");
      strcat(s_tmp,subtrees[0]);
      strcat(s_tmp,",\0");
      strcat(s_tmp,subtrees[1]);
      strcat(s_tmp,")\0");
      Free(subtrees[0]);
      subtrees[0] = s_tmp;

      for(i=1;i<current_deg-1;i++) strcpy(subtrees[i],subtrees[i+1]);

      Clean_Multifurcation(subtrees,current_deg-1,end_deg);
    }
}

/*********************************************************/

char **Sub_Trees(char *tree, int *degree)
{
  char **subs;
  int posbeg,posend;
  int i;

  if(tree[0] != '(') {*degree = 1; return NULL;}

  subs=(char **)mCalloc(NODE_DEG_MAX,sizeof(char *));

  For(i,NODE_DEG_MAX) subs[i]=(char *)mCalloc(strlen(tree)+1,sizeof(char));

  
  posbeg=posend=1;
  (*degree)=0;
  do
    {
      posbeg = posend;
      if(tree[posend] != '(')
	{
	  while((tree[posend] != ',' ) &&
		(tree[posend] != ':' ) &&
		(tree[posend] != '#' ) &&
		(tree[posend] != ')' )) 
	    {
	      posend++ ;
	    }
	  posend -= 1;
	}
      else posend=Next_Par(tree,posend);

      while((tree[posend+1] != ',') &&
	    (tree[posend+1] != ':') &&
	    (tree[posend+1] != '#') &&
	    (tree[posend+1] != ')')) {posend++;}


      strncpy(subs[(*degree)],tree+posbeg,posend-posbeg+1);
/*       strcat(subs[(*degree)],"\0"); */
      subs[(*degree)][posend-posbeg+1]='\0'; /* Thanks to Jean-Baka Domelevo-Entfellner */

      posend += 1;
      while((tree[posend] != ',') &&
	    (tree[posend] != ')')) {posend++;}
      posend+=1;


      (*degree)++;
      if((*degree) == NODE_DEG_MAX)
	{
	  For(i,(*degree))
	    PhyML_Printf("\n. Subtree %d : %s\n",i+1,subs[i]);

	  PhyML_Printf("\n. The degree of a t_node cannot be greater than %d\n",NODE_DEG_MAX);
	  Warn_And_Exit("\n");
	}
    }
  while(tree[posend-1] != ')');

  return subs;
}


/*********************************************************/

int Next_Par(char *s, int pos)
{
  int curr;

  curr=pos+1;

  while(*(s+curr) != ')')
    {
      if(*(s+curr) == '(') curr=Next_Par(s,curr);
      curr++;
    }

  return curr;
}

/*********************************************************/

void Print_Tree(FILE *fp, t_tree *tree)
{
  char *s_tree;
  int i;

  s_tree = (char *)Write_Tree(tree);

  if(OUTPUT_TREE_FORMAT == NEWICK) PhyML_Fprintf(fp,"%s\n",s_tree);
  else if(OUTPUT_TREE_FORMAT == NEXUS)
    {
      PhyML_Fprintf(fp,"#NEXUS\n");
      PhyML_Fprintf(fp,"BEGIN TREES;\n");
      PhyML_Fprintf(fp,"\tTRANSLATE\n");
      For(i,tree->n_otu) PhyML_Fprintf(fp,"\t%3d\t%s,\n",i+1,tree->noeud[i]->name);
      PhyML_Fprintf(fp,"\tUTREE PAUP_1=\n");
      PhyML_Fprintf(fp,"%s\n",s_tree);
      PhyML_Fprintf(fp,"ENDBLOCK;");
    }
  Free(s_tree);
}

/*********************************************************/

char *Write_Tree(t_tree *tree)
{
  char *s;
  int i,available;
  int pos;

  s=(char *)mCalloc((int)T_MAX_NAME,sizeof(char));
  available = (int)T_MAX_NAME-1;

  
  s[0]='(';
  pos = 1;
  
  if(!tree->n_root)
    {
      i = 0;
      while((!tree->noeud[tree->n_otu+i]->v[0]) ||
	    (!tree->noeud[tree->n_otu+i]->v[1]) ||
	    (!tree->noeud[tree->n_otu+i]->v[2])) i++;
      
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[0],&available,&s,&pos,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[1],&available,&s,&pos,tree);
      R_wtree(tree->noeud[tree->n_otu+i],tree->noeud[tree->n_otu+i]->v[2],&available,&s,&pos,tree);
    }
  else
    {
      R_wtree(tree->n_root,tree->n_root->v[0],&available,&s,&pos,tree);
      R_wtree(tree->n_root,tree->n_root->v[1],&available,&s,&pos,tree);
    }

  s[pos-1]=')';
  s[pos]=';';
  s[pos+1]='\0';

  return s;
}

/*********************************************************/

void R_wtree(t_node *pere, t_node *fils, int *available, char **s_tree, int *pos, t_tree *tree)
{
  int i,p,ori_len;
  char *format;

  format = (char *)mCalloc(100,sizeof(char));

  sprintf(format,"%%.%df",tree->bl_ndigits);
  /* strcpy(format,"%f"); */

  p = -1;
  if(fils->tax)
    {
/*       printf("\n- Writing on %p",*s_tree); */
      ori_len = *pos;

      if(OUTPUT_TREE_FORMAT == NEWICK)
	{
	  if(tree->write_tax_names == YES)
	    {
	      if(tree->io && tree->io->long_tax_names) 
		{
		  strcat(*s_tree,tree->io->long_tax_names[fils->num]);
		  (*pos) += (int)strlen(tree->io->long_tax_names[fils->num]);
		}
	      else
		{
		  strcat(*s_tree,fils->name);
		  (*pos) += (int)strlen(fils->name);
		}	  
	    }
	  else if(tree->write_tax_names == NO)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%d",fils->num);
	    }
	}
      else if(OUTPUT_TREE_FORMAT == NEXUS)
	{
	  (*pos) += sprintf(*s_tree+*pos,"%d",fils->num+1);
	}
      else
	{
	  PhyML_Printf("\n. Unknown tree format.");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  PhyML_Printf("\n. s=%s\n",*s_tree);
	}

      if((fils->b) && (fils->b[0]) && (fils->b[0]->l > -1.))
	{
	  /* if(tree->print_labels) */
	  /*   { */
	  /*     if(fils->b[0]->n_labels < 10) */
	  /* 	For(i,fils->b[0]->n_labels)  */
	  /* 	  { */
	  /* 	    (*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[0]->labels[i]); */
	  /* 	  } */
	  /*     else */
	  /* 	{ */
	  /* 	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[0]->n_labels); */
	  /* 	} */
	  /*   } */

	  strcat(*s_tree,":");
	  (*pos)++;

#ifndef PHYTIME
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[0]->l);
	    }
	  else
	    {
	      if(pere == tree->n_root)
		{
		  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
		  (*pos) += sprintf(*s_tree+*pos,format,tree->e_root->l * root_pos);
		}
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,format,fils->b[0]->l);
		}
	    }		
#else
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[0]->l);
	    }
	  else
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,tree->rates->cur_l[fils->num]);
	    }
#endif

	  /* !!!!!!!!!!!!!!!!!!!!1 */
	  if(tree->print_labels)
	    {
	      if(fils->b[0]->n_labels < 10)
		For(i,fils->b[0]->n_labels) 
		  {
		    (*pos) += sprintf(*s_tree+*pos,"::%s",fils->b[0]->labels[i]);
		  }
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,"::%d_labels",fils->b[0]->n_labels);
		}
	    }

	}

      strcat(*s_tree,",");
      (*pos)++;

      (*available) = (*available) - (*pos - ori_len);

      if(2 < 1)
	{
	  PhyML_Printf("\n. s=%s\n",*s_tree);
	  PhyML_Printf("\n. len=%d\n",strlen(*s_tree));
	  PhyML_Printf("\n. The sequence names in your input file might be too long.");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(*available < (int)T_MAX_NAME/2)
	{
	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
	  (*available) = (int)T_MAX_NAME;
	}
/*       printf(" %s [%d,%d]",*s_tree,(int)strlen(*s_tree),*available); */
    }
  else
    {

      (*s_tree)[(*pos)]='(';
      (*s_tree)[(*pos)+1]='\0';
      (*pos)++;
      (*available)--;

      if(tree->n_root)
	{
	  For(i,3)
	    {
	      if((fils->v[i] != pere) && (fils->b[i] != tree->e_root))
		R_wtree(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
	}
      else
	{
	  For(i,3)
	    {
	      if(fils->v[i] != pere)
		R_wtree(fils,fils->v[i],available,s_tree,pos,tree);
	      else p=i;
	    }
	}

      ori_len = *pos;
      
      if(p < 0)
	{
	  PhyML_Printf("\n. fils=%p root=%p root->v[0]=%p root->v[1]=%p",fils,tree->n_root,tree->n_root->v[0],tree->n_root->v[1]);
	  PhyML_Printf("\n. tree->e_root=%p fils->b[0]=%p fils->b[1]=%p fils->b[2]=%p",tree->e_root,fils->b[0],fils->b[1],fils->b[2]);		       
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

/*       printf("\n+ Writing on %p",*s_tree); */
      (*s_tree)[(*pos)-1] = ')';
      (*s_tree)[(*pos)]   = '\0';

      if((fils->b) && (fils->b[p]->l > -1.))
	{
	  if(tree->print_boot_val)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%d",fils->b[p]->bip_score);
	    }
	  else if(tree->print_alrt_val)
	    {
	      (*pos) += sprintf(*s_tree+*pos,"%f",fils->b[p]->ratio_test);
	    }
	  
	  fflush(NULL);

	  /* if((tree->print_labels) && (fils->b[p]->labels != NULL)) */
	  /*   { */
	  /*     if(fils->b[p]->n_labels < 10) */
	  /* 	For(i,fils->b[p]->n_labels)  */
	  /* 	  { */
	  /* 	    (*pos) += sprintf(*s_tree+*pos,"#%s",fils->b[p]->labels[i]); */
	  /* 	  } */
	  /*     else */
	  /* 	{ */
	  /* 	  (*pos) += sprintf(*s_tree+*pos,"#%d_labels",fils->b[p]->n_labels); */
	  /* 	} */
	  /*   } */

	  strcat(*s_tree,":");
	  (*pos)++;

#ifndef PHYTIME
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[p]->l);
	    }
	  else
	    {
	      if(pere == tree->n_root)
		{
		  phydbl root_pos = (fils == tree->n_root->v[0])?(tree->n_root_pos):(1.-tree->n_root_pos);
		  (*pos) += sprintf(*s_tree+*pos,format,tree->e_root->l * root_pos);
		}
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,format,fils->b[p]->l);
		}
	    }
#else
	  if(!tree->n_root)
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,fils->b[p]->l);
	    }
	  else
	    {
	      (*pos) += sprintf(*s_tree+*pos,format,tree->rates->cur_l[fils->num]);
	    }
#endif

	  /* !!!!!!!!!!!!!!!!!!!!!!!!!!!!!1 */
	  if((tree->print_labels) && (fils->b[p]->labels != NULL))
	    {
	      if(fils->b[p]->n_labels < 10)
		For(i,fils->b[p]->n_labels) 
		  {
		    (*pos) += sprintf(*s_tree+*pos,"::%s",fils->b[p]->labels[i]);
		  }
	      else
		{
		  (*pos) += sprintf(*s_tree+*pos,"::%d_labels",fils->b[p]->n_labels);
		}
	    }

	}
      strcat(*s_tree,",");
      (*pos)++;
      (*available) = (*available) - (*pos - ori_len);

      if(*available < 0)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      if(*available < (int)T_MAX_NAME/2)
	{
	  (*s_tree) = (char *)mRealloc(*s_tree,*pos+(int)T_MAX_NAME,sizeof(char));
	  (*available) = (int)T_MAX_NAME;
	}
/*       printf(" %s [%d,%d]",*s_tree,(int)strlen(*s_tree),*available); */
    }

  Free(format);
}

/*********************************************************/

void Init_Tree(t_tree *tree, int n_otu)
{
  tree->n_otu                     = n_otu;
  tree->mat                       = NULL;
  tree->n_root                    = NULL;
  tree->e_root                    = NULL;
  tree->ps_tree                   = NULL;
  tree->short_l                   = NULL;

  tree->depth_curr_path           = 0;
  tree->has_bip                   = NO;
  tree->n_moves                   = 0;
  tree->n_improvements            = 0;
  tree->bl_from_node_stamps       = 0;
  tree->lock_topo                 = 0;
  tree->ps_page_number            = 0;
  tree->init_lnL                  = UNLIKELY;
  tree->best_lnL                  = UNLIKELY;
  tree->old_lnL                   = UNLIKELY;
  tree->c_lnL                     = UNLIKELY;
  tree->sum_min_sum_scale         = .0;
  tree->n_swap                    = 0;
  tree->best_pars                 = 1E+5;
  tree->n_pattern                 = -1;
  tree->n_root_pos                = -1.;
  tree->print_labels              = 1;
  tree->print_boot_val            = 0;
  tree->print_alrt_val            = 0;
  tree->num_curr_branch_available = 0;
  tree->tip_order_score           = .0;
  tree->write_tax_names           = YES;
  tree->update_alias_subpatt      = NO;
  tree->bl_ndigits                = 8;
  tree->n_short_l                 = 100;
  tree->norm_scale                = 0.0;
}

/*********************************************************/

void Make_New_Edge_Label(t_edge *b)
{
  int i;

  b->labels = (char **)realloc(b->labels,(b->n_labels+BLOCK_LABELS)*sizeof(char *));

  if(!b->labels)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      for(i=b->n_labels;i<b->n_labels+100;i++) b->labels[i] = (char *)mCalloc(T_MAX_LABEL,sizeof(char));
    }
}

/*********************************************************/

t_edge *Make_Edge_Light(t_node *a, t_node *d, int num)
{
  t_edge *b;

  b = (t_edge *)mCalloc(1,sizeof(t_edge));

  Init_Edge_Light(b,num);

  if(a && b)
    {
      b->left = a;  b->rght = d;
      if(a->tax) {b->rght = a; b->left = d;} /* root */
      /* a tip is necessary on the right side of the t_edge */

      (b->left == a)?
	(Make_Edge_Dirs(b,a,d,NULL)):
	(Make_Edge_Dirs(b,d,a,NULL));

      b->l             = a->l[b->l_r];
      if(a->tax) b->l  = a->l[b->r_l];
      b->l_old         = b->l;
    }
  else
    {
      b->left = NULL;
      b->rght = NULL;
    }

  return b;

}

/*********************************************************/

void Init_Edge_Light(t_edge *b, int num)
{
  b->num                  = num;
  b->bip_score            = 0;
  b->dist_btw_edges       = .0;
  b->topo_dist_btw_edges  = 0;
  b->has_zero_br_len      = NO;
  b->n_jumps              = 0;
  b->gamma_prior_mean     = 1.E-0;
  b->gamma_prior_var      = 1.E-1;

  b->p_lk_left            = NULL;
  b->p_lk_rght            = NULL;
  b->p_lk_loc_left        = NULL;
  b->p_lk_loc_rght        = NULL;
  b->Pij_rr               = NULL;
}

/*********************************************************/

void Init_Node_Light(t_node *n, int num)
{
  n->num                    = num;
  n->tax                    = -1;
  n->dist_to_root           = .0;
  n->common                 = 1;
  n->ext_node               = NULL;
  n->name                   = NULL;
  n->ori_name               = NULL;
  n->y_rank                 = 0.;
  n->y_rank_ori             = 0.;
  n->y_rank_max             = 0.;
  n->y_rank_min             = 0.;
  n->anc                    = NULL;
  n->rank                   = 0;
}

/*********************************************************/

void Make_Edge_Dirs(t_edge *b, t_node *a, t_node *d, t_tree *tree)
{
  int i;
  t_edge *e_root;

  e_root = (tree)?(tree->e_root):(NULL);
 
  if(a == b->rght)
    {
      PhyML_Printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(d == b->left)
    {
      PhyML_Printf("\n. a->num = %3d ; d->num = %3d",a->num,d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  b->l_r = b->r_l = -1;
  For(i,3)
    {
      if((a->v[i]) && ((a->v[i] == d) || (e_root && a->b[i] == e_root)))
	{
	  b->l_r  = i; /* we consider here that 'a' is on the left handside of 'b'*/
	  a->b[i] = b;
	}
      if((d->v[i]) && ((d->v[i] == a) || (e_root && d->b[i] == e_root)))
	{
	  b->r_l  = i; /* we consider here that 'd' is on the right handside of 'b'*/
	  d->b[i] = b;
	}
    }

  if(a->tax) {b->r_l = 0; For(i,3) if(d->v[i]==a) {b->l_r = i; break;}}

  b->l_v1 = b->l_v2 = b->r_v1 = b->r_v2 = -1;
  For(i,3)
    {
      if(b->left->v[i] != b->rght)
	{
	  if(b->l_v1 < 0) b->l_v1 = i;
	  else            b->l_v2 = i;
	}

      if(b->rght->v[i] != b->left)
	{
	  if(b->r_v1 < 0) b->r_v1 = i;
	  else            b->r_v2 = i;
	}
    }
}

/*********************************************************/


void Make_Edge_Pars(t_edge *b, t_tree *tree)
{
/*   int site; */

  b->pars_l = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->pars_r = (int *)mCalloc(tree->data->crunch_len,sizeof(int));


  b->ui_l = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));
  b->ui_r = (unsigned int *)mCalloc(tree->data->crunch_len,sizeof(unsigned int));


  b->p_pars_l = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));
  b->p_pars_r = (int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(int ));

}

/*********************************************************/

void Make_Edge_Lk(t_edge *b, t_tree *tree)
{
  int ns;
  
  ns = -1;

  ns = tree->mod->ns;

  b->l_old = b->l;

  b->div_post_pred_left = (short int *)mCalloc(ns,sizeof(short int));
  b->div_post_pred_rght = (short int *)mCalloc(ns,sizeof(short int));

  b->Pij_rr = (phydbl *)mCalloc(tree->mod->n_catg*tree->mod->ns*tree->mod->ns,sizeof(phydbl));
  
  b->sum_scale_left_cat = (int *)mCalloc(tree->mod->n_catg,sizeof(int));  
  b->sum_scale_rght_cat = (int *)mCalloc(tree->mod->n_catg,sizeof(int));

  if(!b->left->tax)
    b->sum_scale_left = (int *)mCalloc(tree->data->crunch_len*tree->mod->n_catg,sizeof(int));
  else
    b->sum_scale_left = NULL;
  
  if(!b->rght->tax)
    b->sum_scale_rght = (int *)mCalloc(tree->data->crunch_len*tree->mod->n_catg,sizeof(int));
  else
    b->sum_scale_rght = NULL;
  
  
  if((!b->left->tax) || (tree->mod->s_opt->greedy))
    {
      b->p_lk_left = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->n_catg*tree->mod->ns,sizeof(phydbl));
      b->p_lk_tip_l = NULL;
    }
  else if(b->left->tax)
    {
      b->p_lk_left   = NULL;      
      b->p_lk_tip_l  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int ));
    }  
  
  if((!b->rght->tax) || (tree->mod->s_opt->greedy))
    {
      b->p_lk_rght = (phydbl *)mCalloc(tree->data->crunch_len*tree->mod->n_catg*tree->mod->ns,sizeof(phydbl));
      b->p_lk_tip_r = NULL;
    }
  else if(b->rght->tax)
    {
      b->p_lk_rght = NULL;      
      b->p_lk_tip_r  = (short int *)mCalloc(tree->data->crunch_len*tree->mod->ns,sizeof(short int));
    }

  b->patt_id_left  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->patt_id_rght  = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->p_lk_loc_left = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
  b->p_lk_loc_rght = (int *)mCalloc(tree->data->crunch_len,sizeof(int));
}


/*********************************************************/

void Make_Edge_NNI(t_edge *b)
{
  b->nni    = Make_NNI();
  b->nni->b = b;
  b->nni->left = b->left;
  b->nni->rght = b->rght;
}

/*********************************************************/

nni *Make_NNI()
{
  nni *a_nni;
  a_nni = (nni *)mCalloc(1,sizeof(nni ));
  Init_NNI(a_nni);
  return a_nni;
}

/*********************************************************/

void Init_NNI(nni *a_nni)
{
  a_nni->left         = NULL;
  a_nni->rght         = NULL;
  a_nni->b            = NULL;
  a_nni->init_l       = -1.;
  a_nni->init_lk      = .0;
  a_nni->score        = +1.0;
  a_nni->best_l       = -1.;
  a_nni->swap_node_v1 = NULL;
  a_nni->swap_node_v2 = NULL;
  a_nni->swap_node_v3 = NULL;
  a_nni->swap_node_v4 = NULL;
  a_nni->lk0          = UNLIKELY;
  a_nni->lk1          = UNLIKELY;
  a_nni->lk2          = UNLIKELY;
  a_nni->l0           = -1.0;
  a_nni->l1           = -1.0;
  a_nni->l2           = -1.0;
}

/*********************************************************/

t_node *Make_Node_Light(int num)
{
  t_node *n;
  n           = (t_node *)mCalloc(1,sizeof(t_node));
  n->v        = (t_node **)mCalloc(3,sizeof(t_node *));
  n->l        = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->b        = (t_edge **)mCalloc(3,sizeof(t_edge *));
  n->score    = (phydbl *)mCalloc(3,sizeof(phydbl));
  n->s_ingrp  = (int *)mCalloc(3,sizeof(int));
  n->s_outgrp = (int *)mCalloc(3,sizeof(int));

  Init_Node_Light(n,num);

  return n;
}

/*********************************************************/


void Make_Node_Lk(t_node *n)
{
/*   n->n_ex_nodes = (int *)mCalloc(2,sizeof(int)); */
  return;
}

/*********************************************************/

void Detect_Align_File_Format(option *io)
{
  int c;
  fpos_t curr_pos;
  
  fgetpos(io->fp_in_align,&curr_pos);
  
  errno = 0;

  while((c=fgetc(io->fp_in_align)) != EOF)
    {
      if(errno) io->data_file_format = PHYLIP;
      else if(c == '#')
	{
	  char s[10],t[6]="NEXUS";
	  if(!fgets(s,6,io->fp_in_align))
	    {     
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(!strcmp(t,s)) 
	    {
	      fsetpos(io->fp_in_align,&curr_pos);
	      io->data_file_format = NEXUS;
	      return;
	    }
	}
    }
  
  fsetpos(io->fp_in_align,&curr_pos);
}

/*********************************************************/

void Detect_Tree_File_Format(option *io)
{
  int c;
  fpos_t curr_pos;
  
  fgetpos(io->fp_in_tree,&curr_pos);

  errno = 0;

  while((c=fgetc(io->fp_in_tree)) != EOF)
    {
      if(errno) 
	{
	  io->tree_file_format = PHYLIP;
	  PhyML_Printf("\n. Detected PHYLIP tree file format.");
	}
      else if(c == '#')
	{
	  char s[10],t[6]="NEXUS";
	  if(!fgets(s,6,io->fp_in_tree))
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(!strcmp(t,s))
	    {
	      fsetpos(io->fp_in_tree,&curr_pos);
	      io->tree_file_format = NEXUS;
	      PhyML_Printf("\n. Detected NEXUS tree file format.");
	      return;
	    }
	}
    }
  
  fsetpos(io->fp_in_tree,&curr_pos);
}

/*********************************************************/

align **Get_Seq(option *io)
{
  io->data = NULL;

  Detect_Align_File_Format(io);

  switch(io->data_file_format)
    {
    case PHYLIP: 
      {
	io->data = Get_Seq_Phylip(io);
	break;
      }
    case NEXUS:
      {
	io->nex_com_list = Make_Nexus_Com();
	Init_Nexus_Format(io->nex_com_list);
	Get_Nexus_Data(io->fp_in_align,io);
	Free_Nexus(io);
	break;
      }
    default:
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
	break;
      }
    }

  if(!io->data)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else
    {
      int i,j;
      char **buff;
      int *remove;
      int n_unkn,n_removed,pos;

      buff = (char **)mCalloc(io->n_otu,sizeof(char *));
      For(i,io->n_otu) buff[i] = (char *)mCalloc(io->data[0]->len,sizeof(char));
      remove = (int *)mCalloc(io->data[0]->len,sizeof(int));

      n_removed = 0;

      For(i,io->data[0]->len)
	{
	  For(j,io->n_otu)
	    {
	      if((io->data[j]->state[i] == '?') || (io->data[j]->state[i] == '-')) io->data[j]->state[i] = 'X';
	      if((io->datatype == NT) && (io->data[j]->state[i] == 'N')) io->data[j]->state[i] = 'X';
	      if(io->data[j]->state[i] == 'U') io->data[j]->state[i] = 'T';
	    }

	  n_unkn = 0;
	  For(j,io->n_otu) if(io->data[j]->state[i] == 'X') n_unkn++;

	  if(n_unkn == io->n_otu)
	    {
	      remove[i] = 1;
	      n_removed++;
	    }

	  For(j,io->n_otu) buff[j][i] = io->data[j]->state[i];
	}

      pos = 0;
      For(i,io->data[0]->len)
	{
/* 	  if(!remove[i]) */
/* 	    { */
	      For(j,io->n_otu) io->data[j]->state[pos] = buff[j][i];
	      pos++;
/* 	    } */
	}

      For(i,io->n_otu) io->data[i]->len = pos;
      For(i,io->n_otu) Free(buff[i]);
      Free(buff);
      Free(remove);
    }


  return io->data;
}

/*********************************************************/

/* align **Get_Seq_Nexus(option *io) */
/* { */
/*   char *s,*ori_s; */
/*   char *token; */
/*   int in_comment; */
/*   nexcom *curr_com; */
/*   nexparm *curr_parm; */
/*   int nxt_token_t,cur_token_t; */

/*   s = (char *)mCalloc(T_MAX_LINE,sizeof(char)); */
/*   token = (char *)mCalloc(T_MAX_TOKEN,sizeof(char)); */
      	  
/*   ori_s      = s; */
/*   in_comment = NO; */
/*   curr_com   = NULL; */
/*   curr_parm  = NULL; */
/*   nxt_token_t = NEXUS_COM;  */
/*   cur_token_t = -1;  */

/*   while(fgets(s,T_MAX_LINE,io->fp_in_align)) */
/*     {       */
/*       do */
/* 	{	   */
/* 	  Get_Token(&s,token);	   */

/* /\* 	  PhyML_Printf("\n. Token: '%s' next_token=%d cur_token=%d",token,nxt_token_t,cur_token_t); *\/ */

/* 	  if(token[0] == '\0') break; */

/* 	  if(token[0] == ';')  */
/* 	    { */
/* 	      curr_com   = NULL; */
/* 	      curr_parm  = NULL; */
/* 	      nxt_token_t = NEXUS_COM; */
/* 	      cur_token_t = -1; */
/* 	      break; /\* End of command *\/  */
/* 	    } */

/* 	  if(nxt_token_t == NEXUS_EQUAL)  */
/* 	    { */
/* 	      cur_token_t = NEXUS_VALUE; */
/* 	      nxt_token_t = NEXUS_PARM; */
/* 	      continue; */
/* 	    } */

/* 	  if((nxt_token_t == NEXUS_COM) && (cur_token_t != NEXUS_VALUE))  */
/* 	    { */
/* 	      Find_Nexus_Com(token,&curr_com,&curr_parm,io->nex_com_list); */
/* 	      if(curr_com)  */
/* 		{ */
/* 		  nxt_token_t = curr_com->nxt_token_t; */
/* 		  cur_token_t = curr_com->cur_token_t; */
/* 		} */
/* 	      if(cur_token_t != NEXUS_VALUE) continue; */
/* 	    } */

/* 	  if((nxt_token_t == NEXUS_PARM) && (cur_token_t != NEXUS_VALUE))  */
/* 	    { */
/* 	      Find_Nexus_Parm(token,&curr_parm,curr_com); */
/* 	      if(curr_parm)  */
/* 		{ */
/* 		  nxt_token_t = curr_parm->nxt_token_t; */
/* 		  cur_token_t = curr_parm->cur_token_t; */
/* 		} */
/* 	      if(cur_token_t != NEXUS_VALUE) continue; */
/* 	    } */

/* 	  if(cur_token_t == NEXUS_VALUE) */
/* 	    { */
/* 	      if((curr_parm->fp)(token,curr_parm,io))  /\* Read in parameter value *\/ */
/* 		{ */
/* 		  nxt_token_t = NEXUS_PARM; */
/* 		  cur_token_t = -1; */
/* 		} */
/* 	    } */
/* 	} */
/*       while(strlen(token) > 0); */
/*     } */

/*   Free(ori_s); */
/*   Free(token); */

/*   return io->data; */
/* } */

/* /\*********************************************************\/ */

void Get_Nexus_Data(FILE *fp, option *io)
{
  char *token;
  int in_comment;
  nexcom *curr_com;
  nexparm *curr_parm;
  int nxt_token_t,cur_token_t;

  token = (char *)mCalloc(T_MAX_TOKEN,sizeof(char));
      	  
  in_comment = NO;
  curr_com   = NULL;
  curr_parm  = NULL;
  nxt_token_t = NEXUS_COM; 
  cur_token_t = -1; 

  do
    {
      if(!Get_Token(fp,token)) break;

/*       PhyML_Printf("\n+ Token: '%s' next_token=%d cur_token=%d",token,nxt_token_t,cur_token_t); */

      if(token[0] == ';') 
	{
	  curr_com    = NULL;
	  curr_parm   = NULL;
	  nxt_token_t = NEXUS_COM;
	  cur_token_t = -1;
	}
      
      if(nxt_token_t == NEXUS_EQUAL) 
	{
	  cur_token_t = NEXUS_VALUE;
	  nxt_token_t = NEXUS_PARM;
	  continue;
	}
      
      if((nxt_token_t == NEXUS_COM) && (cur_token_t != NEXUS_VALUE)) 
	{
	  Find_Nexus_Com(token,&curr_com,&curr_parm,io->nex_com_list);
	  if(curr_com) 
	    {
	      nxt_token_t = curr_com->nxt_token_t;
	      cur_token_t = curr_com->cur_token_t;
	    }
	  if(cur_token_t != NEXUS_VALUE) continue;
	}
      
      if((nxt_token_t == NEXUS_PARM) && (cur_token_t != NEXUS_VALUE)) 
	{
	  Find_Nexus_Parm(token,&curr_parm,curr_com);
	  if(curr_parm) 
	    {
	      nxt_token_t = curr_parm->nxt_token_t;
	      cur_token_t = curr_parm->cur_token_t;
	    }
	  if(cur_token_t != NEXUS_VALUE) continue;
	}
      
      if(cur_token_t == NEXUS_VALUE)
	{
	  if((curr_parm->fp)(token,curr_parm,io))  /* Read in parameter value */
	    {
	      nxt_token_t = NEXUS_PARM;
	      cur_token_t = -1;
	    }
	}
    }
  while(strlen(token) > 0);
  
  Free(token);
}

/*********************************************************/

int Get_Token(FILE *fp, char *token)
{
  char c;
  
  c = ' ';
  while(c == ' ' || c == '\t' || c == '\n') 
    {
      c = fgetc(fp);
      if(c == EOF) return 0;
    }

  if(c == '"')
    {
      do
	{
	  *token = c;
	  token++;
	  c = fgetc(fp);
	  if(c == EOF) return 0;
	}
      while(c != '"');
      *token = c;
      c = fgetc(fp);
      if(c == EOF) return 0;
      *(token+1) = '\0';
      return 1;
    }

  if(c == '[')
    {
      Skip_Comment(fp);
      c = fgetc(fp);
      if(c == EOF) return 0;
      return 1;
    }

  if(c == '#')      { *token = c; token++; }
  else if(c == ';') { *token = c; token++; }
  else if(c == ',') { *token = c; token++; }
  else if(c == '.') { *token = c; token++; }
  else if(c == '=') { *token = c; token++; }
  else if(c == '(') { *token = c; token++; }
  else if(c == ')') { *token = c; token++; }
  else if(c == '{') { *token = c; token++; }
  else if(c == '}') { *token = c; token++; }
  else if(c == '?') { *token = c; token++; }
  else if(c == '-') { *token = c; token++; }
  else
    {
      while(isgraph(c) && c != ';' && c != '-' && c != ',')
	{
	  *(token++) = c;
	  c = fgetc(fp);
	  if(c == EOF) return 0;
	}

      fseek(fp,-1*sizeof(char),SEEK_CUR);

    }
  *token = '\0';
  return 1;
}


/* void Get_Token(char *line, char *token) */
/* { */
/*   while(**line == ' ' || **line == '\t') (*line)++; */


/*   if(**line == '"')  */
/*     { */
/*       do { *token = **line; (*line)++; token++; } while(**line != '"'); */
/*       *token = **line; */
/*       (*line)++; */
/*       *(token+1) = '\0'; */
/*       return; */
/*     } */

/*   if(**line == '[')  */
/*     { */
/*       int in_comment; */

/*       in_comment = 1; */
/*       do  */
/* 	{  */
/* 	  (*line)++;  */
/* 	  if(**line == '[')  */
/* 	    { */
/* 	      in_comment++; */
/* 	    } */
/* 	  else if(**line == ']') in_comment--;	   */
/* 	} */
/*       while(in_comment); */
/*       (*line)++; */
/*       return; */
/*     } */


/*   if(**line == '#')      {*token = **line; (*line)++; token++; } */
/*   else if(**line == ';') {*token = **line; (*line)++; token++; } */
/*   else if(**line == ',') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '.') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '=') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '(') {*token = **line; (*line)++; token++; } */
/*   else if(**line == ')') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '{') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '}') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '?') {*token = **line; (*line)++; token++; } */
/*   else if(**line == '-') {*token = **line; (*line)++; token++; } */
/*   else */
/*     { */
/*       while(isgraph(**line) && **line != ';' && **line != '=' && **line != ',')  */
/* 	{ */
/* 	  *(token++) = **line; */
/* 	  (*line)++;  */
/* 	} */
/*     } */
/*   *token = '\0'; */
/* } */

/*********************************************************/

nexcom **Make_Nexus_Com()
{
  nexcom **com;
  int i;

  com = (nexcom **)mCalloc(N_MAX_NEX_COM,sizeof(nexcom *));
  
  For(i,N_MAX_NEX_COM)
    {
      com[i]       = (nexcom *)mCalloc(1,sizeof(nexcom));
      com[i]->name = (char *)mCalloc(T_MAX_NEX_COM,sizeof(char));
      com[i]->parm = (nexparm **)mCalloc(N_MAX_NEX_PARM,sizeof(nexparm *));
    }

  return com;
}

/*********************************************************/

nexparm *Make_Nexus_Parm()
{
  nexparm *parm;

  parm        = (nexparm *)mCalloc(1,sizeof(nexparm));  
  parm->name  = (char *)mCalloc(T_MAX_TOKEN,sizeof(char ));
  parm->value = (char *)mCalloc(T_MAX_TOKEN,sizeof(char ));

  return parm;
}

/*********************************************************/

void Init_Nexus_Format(nexcom **com)
{

  /*****************************/

  strcpy(com[0]->name,"dimensions");
  com[0]->nparm = 2;
  com[0]->nxt_token_t = NEXUS_PARM;
  com[0]->cur_token_t = NEXUS_COM;

  com[0]->parm[0] = Make_Nexus_Parm();
  strcpy(com[0]->parm[0]->name,"ntax");
  com[0]->parm[0]->fp = Read_Nexus_Dimensions;
  com[0]->parm[0]->com = com[0];
  com[0]->parm[0]->nxt_token_t = NEXUS_EQUAL;
  com[0]->parm[0]->cur_token_t = NEXUS_PARM;

  com[0]->parm[1] = Make_Nexus_Parm();
  strcpy(com[0]->parm[1]->name,"nchar");
  com[0]->parm[1]->fp = Read_Nexus_Dimensions;
  com[0]->parm[1]->com = com[0];
  com[0]->parm[1]->nxt_token_t = NEXUS_EQUAL;
  com[0]->parm[1]->cur_token_t = NEXUS_PARM;

  /*****************************/

  strcpy(com[1]->name,"format");
  com[1]->nparm = 11;
  com[1]->nxt_token_t = NEXUS_PARM;
  com[1]->cur_token_t = NEXUS_COM;

  com[1]->parm[0] = Make_Nexus_Parm();
  strcpy(com[1]->parm[0]->name,"datatype");
  com[1]->parm[0]->fp = Read_Nexus_Format;
  com[1]->parm[0]->com = com[1];
  com[1]->parm[0]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[0]->cur_token_t = NEXUS_PARM;

  com[1]->parm[1] = Make_Nexus_Parm();
  strcpy(com[1]->parm[1]->name,"respectcase");
  com[1]->parm[1]->fp = Read_Nexus_Format;
  com[1]->parm[1]->com = com[1];
  com[1]->parm[1]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[1]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[2] = Make_Nexus_Parm();
  strcpy(com[1]->parm[2]->name,"missing");
  com[1]->parm[2]->fp = Read_Nexus_Format;
  com[1]->parm[2]->com = com[1];
  com[1]->parm[2]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[2]->cur_token_t = NEXUS_PARM;

  com[1]->parm[3] = Make_Nexus_Parm();
  strcpy(com[1]->parm[3]->name,"gap");
  com[1]->parm[3]->fp = Read_Nexus_Format;
  com[1]->parm[3]->com = com[1];
  com[1]->parm[3]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[3]->cur_token_t = NEXUS_PARM;

  com[1]->parm[4] = Make_Nexus_Parm();
  strcpy(com[1]->parm[4]->name,"symbols");
  com[1]->parm[4]->fp = Read_Nexus_Format;
  com[1]->parm[4]->com = com[1];
  com[1]->parm[4]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[4]->cur_token_t = NEXUS_PARM;

  com[1]->parm[5] = Make_Nexus_Parm();
  strcpy(com[1]->parm[5]->name,"equate");
  com[1]->parm[5]->fp = Read_Nexus_Format;
  com[1]->parm[5]->com = com[1];
  com[1]->parm[5]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[5]->cur_token_t = NEXUS_PARM;

  com[1]->parm[6] = Make_Nexus_Parm();
  strcpy(com[1]->parm[6]->name,"matchchar");
  com[1]->parm[6]->fp = Read_Nexus_Format;
  com[1]->parm[6]->com = com[1];
  com[1]->parm[6]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[6]->cur_token_t = NEXUS_PARM;

  com[1]->parm[7] = Make_Nexus_Parm();
  strcpy(com[1]->parm[7]->name,"transpose");
  com[1]->parm[7]->fp = Read_Nexus_Format;
  com[1]->parm[7]->com = com[1];
  com[1]->parm[7]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[7]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[8] = Make_Nexus_Parm();
  strcpy(com[1]->parm[8]->name,"interleave");
  com[1]->parm[8]->fp = Read_Nexus_Format;
  com[1]->parm[8]->com = com[1];
  com[1]->parm[8]->nxt_token_t = NEXUS_PARM;
  com[1]->parm[8]->cur_token_t = NEXUS_VALUE;

  com[1]->parm[9] = Make_Nexus_Parm();
  strcpy(com[1]->parm[9]->name,"items");
  com[1]->parm[9]->fp = Read_Nexus_Format;
  com[1]->parm[9]->com = com[1];
  com[1]->parm[9]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[9]->cur_token_t = NEXUS_PARM;

  com[1]->parm[10] = Make_Nexus_Parm();
  strcpy(com[1]->parm[10]->name,"statesformat");
  com[1]->parm[10]->fp = Read_Nexus_Format;
  com[1]->parm[10]->com = com[1];
  com[1]->parm[10]->nxt_token_t = NEXUS_EQUAL;
  com[1]->parm[10]->cur_token_t = NEXUS_PARM;

  /*****************************/

  strcpy(com[2]->name,"eliminate");
  com[2]->nparm = 0;
  com[2]->nxt_token_t = NEXUS_VALUE;
  com[2]->cur_token_t = NEXUS_COM;

  /*****************************/

  strcpy(com[3]->name,"taxlabels");
  com[3]->nparm = 0;
  com[3]->nxt_token_t = -1;
  com[3]->cur_token_t = -1;
 
 /*****************************/

  strcpy(com[4]->name,"charstatelabels");
  com[4]->nparm = 0;
  com[4]->nxt_token_t = -1;
  com[4]->cur_token_t = -1;

  /*****************************/

  strcpy(com[5]->name,"charlabels");
  com[5]->nparm = 0;
  com[5]->nxt_token_t = -1;
  com[5]->cur_token_t = -1;

  /*****************************/

  strcpy(com[6]->name,"statelabels");
  com[6]->nparm = 0;
  com[6]->nxt_token_t = -1;
  com[6]->cur_token_t = -1;

  /*****************************/

  strcpy(com[7]->name,"matrix");
  com[7]->nparm = 1;
  com[7]->nxt_token_t = NEXUS_COM;
  com[7]->cur_token_t = NEXUS_VALUE; /* This will allow us to skip directly 
					to the matrix reading function */

  com[7]->parm[0] = Make_Nexus_Parm();
  strcpy(com[7]->parm[0]->name,"matrix");
  com[7]->parm[0]->fp = Read_Nexus_Matrix;
  com[7]->parm[0]->com = com[7];
  com[7]->parm[0]->nxt_token_t = NEXUS_COM;
  com[7]->parm[0]->cur_token_t = -1; 

  /*****************************/

  strcpy(com[8]->name,"begin");
  com[8]->nparm = 3;

  com[8]->nxt_token_t = NEXUS_PARM;
  com[8]->cur_token_t = NEXUS_COM;

  com[8]->parm[0] = Make_Nexus_Parm();
  strcpy(com[8]->parm[0]->name,"data");
  com[8]->parm[0]->fp = Read_Nexus_Begin;
  com[8]->parm[0]->com = com[8];
  com[8]->parm[0]->nxt_token_t = NEXUS_COM;
  com[8]->parm[0]->cur_token_t = NEXUS_PARM;


  com[8]->parm[1] = Make_Nexus_Parm();
  strcpy(com[8]->parm[1]->name,"trees");
  com[8]->parm[1]->fp = Read_Nexus_Begin;
  com[8]->parm[1]->com = com[8];
  com[8]->parm[1]->nxt_token_t = NEXUS_COM;
  com[8]->parm[1]->cur_token_t = NEXUS_PARM;


  com[8]->parm[2] = Make_Nexus_Parm();
  strcpy(com[8]->parm[2]->name,"taxa");
  com[8]->parm[2]->fp = Read_Nexus_Taxa;
  com[8]->parm[2]->com = com[8];
  com[8]->parm[2]->nxt_token_t = NEXUS_COM;
  com[8]->parm[2]->cur_token_t = NEXUS_VALUE; 


  /*****************************/

  strcpy(com[9]->name,"end");
  com[9]->nparm = 0;
  com[9]->nxt_token_t = -1;
  com[9]->cur_token_t = -1;
 
  /*****************************/

  strcpy(com[10]->name,"translate");
  com[10]->nparm = 1;
  com[10]->nxt_token_t = NEXUS_COM;
  com[10]->cur_token_t = NEXUS_VALUE;

  com[10]->parm[0] = Make_Nexus_Parm();
  strcpy(com[10]->parm[0]->name,"translate");
  com[10]->parm[0]->fp = Read_Nexus_Translate;
  com[10]->parm[0]->com = com[10];
  com[10]->parm[0]->nxt_token_t = NEXUS_COM;
  com[10]->parm[0]->cur_token_t = -1; 

  /*****************************/

  strcpy(com[11]->name,"tree");
  com[11]->nparm = 1;
  com[11]->nxt_token_t = NEXUS_COM;
  com[11]->cur_token_t = NEXUS_VALUE;

  com[11]->parm[0] = Make_Nexus_Parm();
  strcpy(com[11]->parm[0]->name,"tree");
  com[11]->parm[0]->fp = Read_Nexus_Tree;
  com[11]->parm[0]->com = com[11];
  com[11]->parm[0]->nxt_token_t = -1;
  com[11]->parm[0]->cur_token_t = -1; 


  /*****************************/

}



/*********************************************************/

void Find_Nexus_Com(char *token, nexcom **found_com, nexparm **default_parm, nexcom **com_list)
{
  int i,j,tokenlen,ndiff;
    
  For(i,N_MAX_NEX_COM) 
    {
      tokenlen = strlen(token);
      ndiff = -1;
      if(tokenlen && (tokenlen == strlen(com_list[i]->name)))
	{
	  ndiff = 0;
	  For(j,tokenlen)
	    {
	      Lowercase(token+j);
	      Lowercase(com_list[i]->name+j);
	      if(token[j] != com_list[i]->name[j]) ndiff++;
	    }
	}
      if(!ndiff) { *found_com = com_list[i]; break; }
    }

  if(*found_com && (*found_com)->nparm) *default_parm = (*found_com)->parm[0];

  if(*found_com) PhyML_Printf("\n. Found command '%s'.\n",(*found_com)->name);
}

/*********************************************************/

void Find_Nexus_Parm(char *token, nexparm **found_parm, nexcom *curr_com)
{
  int i,j;
  int tokenlen;
  int ndiff;

  if(!curr_com)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  For(i,curr_com->nparm)
    {
      tokenlen = strlen(token);
      ndiff = -1;
      if(tokenlen == strlen(curr_com->parm[i]->name))
	{
	  ndiff = 0;
	  For(j,tokenlen)
	    {
	      Lowercase(token+j);
	      Lowercase(curr_com->parm[i]->name+j);
	      if(token[j] != curr_com->parm[i]->name[j]) ndiff++;
	    }
	}
      if(!ndiff) { *found_parm = curr_com->parm[i]; break; }
    }

  if(*found_parm) PhyML_Printf("\n. Found parameter '%s'.\n",(*found_parm)->name);
}

/*********************************************************/

int Read_Nexus_Taxa(char *token, nexparm *curr_parm, option *io)
{

  PhyML_Printf("\n. Skipping 'taxa' block");

  do
    {
      Get_Token(io->fp_in_tree,token);
      if(token[0] == ';') break;
    }while(strlen(token) > 0);
  
  fseek(io->fp_in_tree,-1*sizeof(char),SEEK_CUR);

  return 1;
}

/*********************************************************/

int Read_Nexus_Translate(char *token, nexparm *curr_parm, option *io)
{
  int tax_num;
  char *end;

  PhyML_Printf("\n. Reading 'translate' block");
  io->size_tax_names = 0;

  do
    {
      Get_Token(io->fp_in_tree,token);
      if(token[0] == ';') break;
      tax_num = (int)strtol(token,&end,10);
      if(*end =='\0' && token[0])
	{
	  io->size_tax_names++;

	  io->short_tax_names = (char **)realloc(io->short_tax_names,io->size_tax_names*sizeof(char *));
	  io->short_tax_names[io->size_tax_names-1] = (char *)mCalloc(strlen(token)+1,sizeof(char));
	  sprintf(io->short_tax_names[io->size_tax_names-1],"%d",tax_num);

	  Get_Token(io->fp_in_tree,token);

	  io->long_tax_names = (char **)realloc(io->long_tax_names,io->size_tax_names*sizeof(char *));
	  io->long_tax_names[io->size_tax_names-1] = (char *)mCalloc(strlen(token)+1,sizeof(char));
	  strcpy(io->long_tax_names[io->size_tax_names-1],token);

/* 	  printf("\n. Copying %s number %d",io->long_tax_names[io->size_long_tax_names-1],tax_num-1); */
	}
    }while(strlen(token) > 0);
  
  fseek(io->fp_in_tree,-1*sizeof(char),SEEK_CUR);

  return 1;
}

/*********************************************************/

int Read_Nexus_Matrix(char *token, nexparm *curr_parm, option *io)
{

  if(io->interleaved) io->data = Read_Seq_Interleaved(io);
  else                io->data = Read_Seq_Sequential(io);

  fseek(io->fp_in_align,-1*sizeof(char),SEEK_CUR);

  return 1;
}

/*********************************************************/

int Read_Nexus_Tree(char *token, nexparm *curr_parm, option *io)
{
  io->treelist->tree = (t_tree **)realloc(io->treelist->tree,(io->treelist->list_size+1)*sizeof(t_tree *));
  io->tree = Read_Tree_File_Phylip(io->fp_in_tree);
  if(!(io->treelist->list_size%10) && io->treelist->list_size > 1) 
    {
      PhyML_Printf("\n. Reading tree %d",io->treelist->list_size);
      if(io->tree->n_root) PhyML_Printf(" (that is a rooted tree)");
      else                 PhyML_Printf(" (that is an unrooted tree)");
    }
  io->treelist->tree[io->treelist->list_size] = io->tree;
  io->treelist->list_size++;
  fseek(io->fp_in_tree,-1*sizeof(char),SEEK_CUR);
  return 1;
}

/*********************************************************/

int Read_Nexus_Begin(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  if(!curr_parm)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!strcmp(curr_parm->name,"data") || !strcmp(curr_parm->name,"trees")) 
    PhyML_Printf("\n. Reading '%s' block.\n",curr_parm->value);
  else
    {
      PhyML_Printf("\n. The '%s' block type is not supported by PhyML. Sorry.\n",curr_parm->name);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  return 1;
}

/*********************************************************/

int Read_Nexus_Dimensions(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  if(!curr_parm)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  strcpy(curr_parm->value,token);

  if(!strcmp(curr_parm->name,"ntax"))
    {
      sscanf(curr_parm->value,"%d",&(io->n_otu));
    }

  if(!strcmp(curr_parm->name,"nchar"))
    {
      sscanf(curr_parm->value,"%d",&(io->init_len));
    }
  return 1;
}

/*********************************************************/

int Read_Nexus_Format(char *token, nexparm *curr_parm, option *io)
{
  int i;  
  
  if(token[0] == '=') return 0;
  
  if(!curr_parm)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  For(i,strlen(token)) Lowercase(token+i);

  strcpy(curr_parm->value,token);
    
  if(!strcmp(curr_parm->name,"datatype"))
    {
      if(!strcmp(curr_parm->value,"standard"))
	{
	  io->datatype = GENERIC;
	  io->mod->whichmodel = JC69;
	  io->mod->s_opt->opt_kappa  = NO;
	  io->mod->s_opt->opt_lambda = NO;
	  io->mod->ns = 2;
	  io->alphabet[0][0] = '0'; io->alphabet[0][1] = '\0';
	  io->alphabet[1][0] = '1'; io->alphabet[1][1] = '\0';
	}

      else if(!strcmp(curr_parm->value,"dna"))
	{
	  io->datatype = NT;
	  io->mod->ns = 4;
	}

      else if(!strcmp(curr_parm->value,"rna"))
	{
	  io->datatype = NT;
	  io->mod->ns = 4;
	}

      else if(!strcmp(curr_parm->value,"nucleotide"))
	{
	  io->datatype = NT;
	  io->mod->ns = 4;
	}

      else if(!strcmp(curr_parm->value,"protein"))
	{
	  io->datatype = AA;
	  io->mod->ns = 20;
	}
      
      else if(!strcmp(curr_parm->value,"continuous"))
	{
	  PhyML_Printf("\n. The 'continuous' format is not supported by PhyML. Sorry.\n");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
    }

  else if(!strcmp(curr_parm->name,"missing"))
    {
      PhyML_Printf("\n. The 'missing' subcommand is not supported by PhyML. Sorry.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  else if(!strcmp(curr_parm->name,"gap"))
    {
      PhyML_Printf("\n. The 'gap' subcommand is not supported by PhyML. Sorry.\n");
      PhyML_Printf("\n. But the characters 'X', '?' and '-' will be considered as indels by default.\n"); 
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  else if(!strcmp(curr_parm->name,"symbols"))
    {
      if(*token != '"' || *(token+strlen(token)-1) != '"')
	{
	  PhyML_Printf("\n. Symbols list is supposed to be displayed between quotation marks (e.g., \"ACTG\").\n");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}


      int i,has_spaces,state_len;

      i          = 0;
      has_spaces = 0;     
      token++; /* Get rid of the first '"' character */
      while(token[i] != '"')  { if(token[i] == ' ') { has_spaces = 1; break; } i++; }

      io->mod->ns = 0;
      if(!has_spaces)
	{
	  while(token[i] != '"') 
	    { 
	      io->alphabet[io->mod->ns][0] = token[i]; 
	      io->alphabet[io->mod->ns][1] = '\0'; 
	      io->mod->ns++;
	      i++;
	      if(io->mod->ns > T_MAX_ALPHABET)
		{
		  PhyML_Printf("\n. The alphabet cannot contain more than %d characters. Sorry.",T_MAX_ALPHABET);
		  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		  Warn_And_Exit("");
		}
	    }
	}
      else
	{
	  i = 0;
	  do
	    {
	      state_len = 0;
	      while(token[i] != ' ' && token[i] != '"') 
		{ 
		  io->alphabet[io->mod->ns][state_len] = token[i];
		  state_len++;
		  i++;
		  if(state_len > T_MAX_STATE)
		    {
		      PhyML_Printf("\n. A state cannot contain more than %d characters. Sorry.\n",T_MAX_STATE);
		      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
		      Warn_And_Exit("");
		    }
		}
	      
	      io->alphabet[io->mod->ns][state_len] = '\0';
	      io->mod->ns++;
	      if(token[i] != '"') i++;
	    }
	  while(token[i] != '"');
	}

      int len;
      len = strlen(io->alphabet[0]);
      For(i,io->mod->ns)
	{
	  if(strlen(io->alphabet[i]) != len)
	    {
	      PhyML_Printf("\n. All character states defined in the symbol list are supposed to have the same length.\n");
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
      io->mod->state_len = len;      

/*       For(i,io->mod->ns) PhyML_Printf("\n. '%s'",io->alphabet[i]); */
    }

  else if(!strcmp(curr_parm->name,"equate"))
    {
      PhyML_Printf("\n. PhyML does not recognize the command '%s' yet. Sorry.",curr_parm->name);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  else if(!strcmp(curr_parm->name,"matchchar"))
    {
      PhyML_Printf("\n. PhyML does not recognize the command '%s' yet. Sorry.",curr_parm->name);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  else if(!strcmp(curr_parm->name,"items"))
    {
      PhyML_Printf("\n. PhyML does not recognize the command '%s' yet. Sorry.",curr_parm->name);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  else if(!strcmp(curr_parm->name,"interleave"))
    {
      io->interleaved = 1;
    }

  return 1;
}

/*********************************************************/

int Read_Nexus_Eliminate(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n. 'Eliminate' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return 1;
}

/*********************************************************/

int Read_Nexus_Taxlabel(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n. 'Taxlabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return 1;
}

/*********************************************************/

int Read_Nexus_Charstatelabels(char *token, nexparm *curr_parm, option *io)
{

  if(token[0] == '=') return 0;

  PhyML_Printf("\n. 'CharStateLabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return 1;
}

/*********************************************************/

int Read_Nexus_Charlabels(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n. 'CharLabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return 1;
}

/*********************************************************/

int Read_Nexus_Statelabels(char *token, nexparm *curr_parm, option *io)
{
  if(token[0] == '=') return 0;

  PhyML_Printf("\n. 'StateLabels' command is not supported by PhyML. Sorry.");
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return 1;
}

/*********************************************************/

align **Get_Seq_Phylip(option *io)
{
  Read_Ntax_Len_Phylip(io->fp_in_align,&io->n_otu,&io->init_len);

  if(io->interleaved) io->data = Read_Seq_Interleaved(io);
  else                io->data = Read_Seq_Sequential(io);

  return io->data;
}

/*********************************************************/

void Read_Ntax_Len_Phylip(FILE *fp ,int *n_otu, int *n_tax)
{
  char *line;
  int readok;
  
  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));  

  readok = 0;
  do
    {
      if(fscanf(fp,"%s",line) == EOF)
	{
	  Free(line); 
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      else
	{
	  if(strcmp(line,"\n") && strcmp(line,"\r") && strcmp(line,"\t"))
	    {
	      sscanf(line,"%d",n_otu);
	      if(*n_otu <= 0) Warn_And_Exit("\n. The number of taxa cannot be negative.\n");

	      if(!fscanf(fp,"%s",line)) Exit("\n");
	      sscanf(line,"%d",n_tax);
	      if(*n_tax <= 0) Warn_And_Exit("\n. The sequence length cannot be negative.\n");
	      else readok = 1;
	    }
	}
    }while(!readok);
  
  Free(line);
}

/*********************************************************/

align **Read_Seq_Sequential(option *io)
{
  int i;
  char *line;
  align **data;
/*   char c; */
  char *format;


  format = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));

/*   while((c=fgetc(in))!='\n'); */
 /*  while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') && (c != '\t')); */

  For(i,io->n_otu)
    {
      data[i]        = (align *)mCalloc(1,sizeof(align));
      data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(io->init_len*io->mod->state_len+1,sizeof(char));

      data[i]->is_ambigu = NULL;
      data[i]->len = 0;

      sprintf(format, "%%%ds", T_MAX_NAME);

      if(!fscanf(io->fp_in_align,format,data[i]->name)) Exit("\n");

      Check_Sequence_Name(data[i]->name);

      while(data[i]->len < io->init_len * io->mod->state_len) Read_One_Line_Seq(&data,i,io->fp_in_align);

      if(data[i]->len != io->init_len * io->mod->state_len)
	{
	  PhyML_Printf("\n. Err: Problem with species %s's sequence (check the format).\n",data[i]->name);
	  PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, io->init_len * io->mod->state_len);
	  Warn_And_Exit("");
	}
    }

  For(i,io->n_otu) data[i]->state[data[i]->len] = '\0';
  
  Free(format);
  Free(line);
  
  return data;
}

/*********************************************************/

align **Read_Seq_Interleaved(option *io)
{
  int i,end,num_block;
  char *line;
  align **data;
/*   char c; */
  char *format;


  line   = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  format = (char *)mCalloc(T_MAX_NAME, sizeof(char));
  data   = (align **)mCalloc(io->n_otu,sizeof(align *));


/*   while(((c=fgetc(io->fp_in_align))!='\n') && (c != ' ') && (c != '\r') && (c != '\t')); */

  end = 0;
  For(i,io->n_otu)
    {

      data[i]        = (align *)mCalloc(1,sizeof(align));
      data[i]->name  = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      data[i]->state = (char *)mCalloc(io->init_len*io->mod->state_len+1,sizeof(char));

      data[i]->len       = 0;
      data[i]->is_ambigu = NULL;

      sprintf(format, "%%%ds", T_MAX_NAME);
/*       sprintf(format, "%%%ds", 10); */
      
      if(!fscanf(io->fp_in_align,format,data[i]->name)) Exit("\n");

      Check_Sequence_Name(data[i]->name);

      if(!Read_One_Line_Seq(&data,i,io->fp_in_align))
	{
	  end = 1;
	  if((i != io->n_otu) && (i != io->n_otu-1))
	    {
	      PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
	      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d\n",data[i]->len, io->init_len * io->mod->state_len);
	      Warn_And_Exit("");
	    }
	  break;
	}
    }
  
  if(data[0]->len == io->init_len * io->mod->state_len) end = 1;

/*   if(end) printf("\n. finished yet '%c'\n",fgetc(io->fp_in_align)); */
  if(!end)
    {

      end = 0;

      num_block = 1;
      do
	{
	  num_block++;

	  /* interblock */
	  if(!fgets(line,T_MAX_LINE,io->fp_in_align)) break;

	  if(line[0] != 13 && line[0] != 10)
	    {
	      PhyML_Printf("\n. One or more missing sequences in block %d.\n",num_block-1);
	      Warn_And_Exit("");
	    }

	  For(i,io->n_otu) if(data[i]->len != io->init_len * io->mod->state_len) break;

	  if(i == io->n_otu) break;

	  For(i,io->n_otu)
	    {
	      if(data[i]->len > io->init_len * io->mod->state_len)
		{
		  PhyML_Printf("\n. Observed length=%d expected length=%d.\n",data[i]->len,io->init_len * io->mod->state_len);
		  PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
		  Warn_And_Exit("");
		}
	      else if(!Read_One_Line_Seq(&data,i,io->fp_in_align))
		{
		  end = 1;
		  if((i != io->n_otu) && (i != io->n_otu-1))
		    {
		      PhyML_Printf("\n. Err: Problem with species %s's sequence.\n",data[i]->name);
		      PhyML_Printf("\n. Observed sequence length: %d, expected length: %d.\n",data[i]->len, io->init_len * io->mod->state_len);
		      Warn_And_Exit("");
		    }
		  break;
		}
	    }
	}while(!end);
    }

  For(i,io->n_otu) data[i]->state[data[i]->len] = '\0';

  For(i,io->n_otu)
    {
      if(data[i]->len != io->init_len * io->mod->state_len)
	{
	  PhyML_Printf("\n. Check sequence '%s' length (expected length: %d, observed length: %d) [OTU %d].\n",data[i]->name,io->init_len,data[i]->len,i+1);
	  Warn_And_Exit("");
	}
    }

  Free(format);
  Free(line);

  return data;
}

/*********************************************************/

int Read_One_Line_Seq(align ***data, int num_otu, FILE *in)
{
  char c = ' ';
  int nchar = 0;
  
  while(1)
    {
/*       if((c == EOF) || (c == '\n') || (c == '\r')) break; */
      
      if((c == 13) || (c == 10)) 
	{
/* 	  PhyML_Printf("[%d %d]\n",c,nchar); fflush(NULL); */
	  if(!nchar)
	    {
	      c=(char)fgetc(in);
	      continue;
	    }
	  else 
	    { 
/* 	      PhyML_Printf("break\n");  */
	      break; 
	    }
	}
      else if(c == EOF)
	{
/* 	  PhyML_Printf("EOL\n"); */
	  break;
	}
      else if((c == ' ') || (c == '\t') || (c == 32)) 
	{
/* 	  PhyML_Printf("[%d]",c); */
	  c=(char)fgetc(in); 
	  continue;
	}

      nchar++;
      Uppercase(&c);
      
      if(c == '.')
	{
	  c = (*data)[0]->state[(*data)[num_otu]->len];
	  if(!num_otu)
	    Warn_And_Exit("\n. Err: Symbol \".\" should not appear in the first sequence\n");
	}
      (*data)[num_otu]->state[(*data)[num_otu]->len]=c;
      (*data)[num_otu]->len++;
/*       PhyML_Printf("%c",c); */
      c = (char)fgetc(in);
      if(c == ';') break;
    }

  if(c == EOF) return 0;
  else return 1;  
}

/*********************************************************/

void Uppercase(char *ch)
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
   *ch = isupper((int)*ch) ? *ch : toupper((int)*ch);
}

/*********************************************************/

void Lowercase(char *ch)
{
  /* convert ch to upper case -- either ASCII or EBCDIC */
  *ch = isupper((int)*ch) ? tolower((int)*ch) : *ch;
}

/*********************************************************/

calign *Compact_Data(align **data, option *io)
{
  calign *cdata_tmp,*cdata;
  int i,j,k,site;
  int n_patt,which_patt,n_invar;
  char **sp_names;
  int n_otu, n_sites;
  pnode *proot;
  int compress;
  int n_ambigu,is_ambigu;

  n_otu      = io->n_otu;
  n_patt     = 0;
  which_patt = 0;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
    {
      sp_names[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(sp_names[i],data[i]->name);
    }

  cdata_tmp = Make_Cseq(n_otu,data[0]->len,io->mod->state_len,data[0]->len,sp_names);
  proot     = (pnode *)Create_Pnode(T_MAX_ALPHABET);
 
  For(i,n_otu) Free(sp_names[i]);
  Free(sp_names);

  if(data[0]->len%io->mod->state_len)
    {
      PhyML_Printf("\n. Sequence length is not a multiple of %d\n",io->mod->state_len);
      Warn_And_Exit("");
    }
  
  compress = io->colalias;
  n_ambigu = 0;
  is_ambigu = 0;

  if(!io->quiet && !compress) { PhyML_Printf("\n. WARNING: sequences are not compressed !\n");}

  Fors(site,data[0]->len,io->mod->state_len)
    {
      if(io->rm_ambigu)
	{
	  is_ambigu = 0;
	  For(j,n_otu) if(Is_Ambigu(data[j]->state+site,io->datatype,io->mod->state_len)) break;
	  if(j != n_otu)
	    {
	      is_ambigu = 1;
	      n_ambigu++;
	    }
	}

      if(!is_ambigu)
	{
	  if(compress)
	    {
	      which_patt = -1;
	      Traverse_Prefix_Tree(site,-1,&which_patt,&n_patt,data,io,proot);
	      if(which_patt == n_patt-1) /* New pattern found */
		{
		  n_patt--;
		  k=n_patt;
		}
	      else
		{
		  k = n_patt-10;
		}
	    }
	  else
	    {
	      k = n_patt;
	    }
	  
	  if(k == n_patt) /* add a new site pattern */
	    {
	      For(j,n_otu)
		Copy_One_State(data[j]->state+site,
			       cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len,
			       io->mod->state_len);
	      	      
	      For(i,n_otu)
		{
		  For(j,n_otu)
		    {
		      if(!(Are_Compatible(cdata_tmp->c_seq[i]->state+n_patt*io->mod->state_len,
					  cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len,
					  io->mod->state_len,
					  io->datatype))) break;
		    }
		  if(j != n_otu) break;
		}
	      
	      if((j == n_otu) && (i == n_otu)) /* all characters at that site are compatible with one another:
						  the site may be invariant */
		{
		  For(j,n_otu)
		    {
		      cdata_tmp->invar[n_patt] = Assign_State(cdata_tmp->c_seq[j]->state+n_patt*io->mod->state_len,
							      io->datatype,
							      io->mod->state_len);
		      if(cdata_tmp->invar[n_patt] > -1.) break; /* It is not actually (at least one state in the column is ambiguous) */
		    }
		}
	      else cdata_tmp->invar[n_patt] = -1;
	      
	      cdata_tmp->sitepatt[site] = n_patt;
	      cdata_tmp->wght[n_patt]  += 1;
	      n_patt                   += 1;
	    }
	  else
	    {
	      cdata_tmp->sitepatt[site]    = which_patt;
	      cdata_tmp->wght[which_patt] += 1;
	    }
	}
    }
  
  data[0]->len -= n_ambigu;
  
  cdata_tmp->init_len                   = data[0]->len;
  cdata_tmp->crunch_len                 = n_patt;
  For(i,n_otu) cdata_tmp->c_seq[i]->len = n_patt;
  
  if(!io->quiet) PhyML_Printf("\n. %d patterns found (out of a total of %d sites). \n",n_patt,data[0]->len);

  if((io->rm_ambigu) && (n_ambigu)) PhyML_Printf("\n. Removed %d columns of the alignment as they contain ambiguous characters (e.g., gaps) \n",n_ambigu);

/*   For(i,n_otu) */
/*     { */
/*       For(j,cdata_tmp->crunch_len) */
/* 	{ */
/* 	  printf("%c",cdata_tmp->c_seq[i]->state[j*io->mod->state_len+1]); */
/* 	} */
/*       printf("\n"); */
/*     } */

  n_invar=0;
  For(i,cdata_tmp->crunch_len) 
    {
      if(cdata_tmp->invar[i] > -1.) n_invar+=(int)cdata_tmp->wght[i];
    }

  if(!io->quiet) PhyML_Printf("\n. %d sites without polymorphism (%.2f%c).\n",n_invar,100.*(phydbl)n_invar/data[0]->len,'%');
  
  cdata_tmp->obs_pinvar = (phydbl)n_invar/data[0]->len;

  n_sites = 0;
  For(i,cdata_tmp->crunch_len) n_sites += cdata_tmp->wght[i];
  if(n_sites != data[0]->len / io->mod->state_len)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(io->datatype == NT)      Get_Base_Freqs(cdata_tmp);
  else if(io->datatype == AA) Get_AA_Freqs(cdata_tmp);
  else {/* Uniform state frequency distribution.*/}

  cdata = Copy_Cseq(cdata_tmp,io);
  
  Free_Cseq(cdata_tmp);
  Free_Prefix_Tree(proot,T_MAX_ALPHABET);

  return cdata;
}

/*********************************************************/

calign *Compact_Cdata(calign *data, option *io)
{
  calign *cdata;
  int i,j,k,site;
  int n_patt,which_patt;
  int n_otu;

  n_otu = data->n_otu;

  cdata         = (calign *)mCalloc(1,sizeof(calign));
  cdata->n_otu  = n_otu;
  cdata->c_seq  = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->wght   = (int *)mCalloc(data->crunch_len,sizeof(int));
  cdata->b_frq  = (phydbl *)mCalloc(io->mod->ns,sizeof(phydbl));
  cdata->ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
  cdata->invar  = (short int *)mCalloc(data->crunch_len,sizeof(short int));

  cdata->crunch_len = cdata->init_len = -1;
  For(j,n_otu)
    {
      cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
      cdata->c_seq[j]->name      = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(cdata->c_seq[j]->name,data->c_seq[j]->name);
      cdata->c_seq[j]->state     = (char *)mCalloc(data->crunch_len,sizeof(char));
      cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(data->crunch_len,sizeof(short int));
      cdata->c_seq[j]->state[0]  = data->c_seq[j]->state[0];
    }


  n_patt = which_patt =  0;

  For(site,data->crunch_len)
    {
      if(data->wght[site])
	{
	  For(k,n_patt)
	    {
	      For(j,n_otu)
		{
		  if(strncmp(cdata->c_seq[j]->state+k*io->mod->state_len,
			     data->c_seq[j]->state+site*io->mod->state_len,
			     io->mod->state_len))
		    break;
		}
	      
	      if(j == n_otu)
		{
		  which_patt = k;
		  break;
		}
	    }
	  
	  /*       /\* TO DO *\/ */
	  /*       k = n_patt; */
	  
	  if(k == n_patt)
	    {
	      For(j,n_otu) Copy_One_State(data->c_seq[j]->state+site*io->mod->state_len,
					  cdata->c_seq[j]->state+n_patt*io->mod->state_len,
					  io->mod->state_len);
	      
	      For(i,n_otu)
		{
		  For(j,n_otu)
		    {
		      if(!(Are_Compatible(cdata->c_seq[i]->state+n_patt*io->mod->state_len,
					  cdata->c_seq[j]->state+n_patt*io->mod->state_len,
					  io->mod->state_len,
					  io->datatype))) break;
		    }
		  if(j != n_otu) break;
		}
	      
	      if((j == n_otu) && (i == n_otu)) 
		{
		  For(j,n_otu)
		    {
		      cdata->invar[n_patt] = Assign_State(cdata->c_seq[j]->state+n_patt*io->mod->state_len,
							    io->datatype,
							    io->mod->state_len);
		      if(cdata->invar[n_patt] > -1.) break;
		    }
		}
	      else cdata->invar[n_patt] = -1;
	      
	      cdata->wght[n_patt] += data->wght[site];
	      n_patt+=1;
	    }
	  else cdata->wght[which_patt] += data->wght[site];
	  
	  /*       Print_Site(cdata,k,n_otu,"\n",io->stepsize); */
	}
    }
  
  cdata->init_len   = data->crunch_len;
  cdata->crunch_len = n_patt;
  For(i,n_otu) cdata->c_seq[i]->len = n_patt;

  if(io->datatype == NT)      Get_Base_Freqs(cdata);
  else if(io->datatype == AA) Get_AA_Freqs(cdata);
  else {/* Not implemented yet */}

  return cdata;
}

/*********************************************************/

void Traverse_Prefix_Tree(int site, int seqnum, int *patt_num, int *n_patt, align **data, option *io, pnode *n)
{
  int ret_val;

  ret_val = -1;

  if(seqnum == io->n_otu-1)
    {
      n->weight++;
      if(n->weight == 1)
	{
	  n->num = *n_patt;
	  (*n_patt) += 1;
	}
      (*patt_num) = n->num;
      return;
    }
  else
    {
      int next_state;

      next_state = -1;
      next_state = Assign_State_With_Ambiguity(data[seqnum+1]->state+site,
					       io->datatype,
					       io->mod->state_len);

      if(!n->next[next_state]) n->next[next_state] = Create_Pnode(T_MAX_ALPHABET);
      Traverse_Prefix_Tree(site,seqnum+1,patt_num,n_patt,data,io,n->next[next_state]);
    }
}

/*********************************************************/

pnode *Create_Pnode(int size)
{
  pnode *n;
  int i;

  n = (pnode *)mCalloc(1,sizeof(pnode ));
  n->next = (pnode **)mCalloc(size,sizeof(pnode *));
  For(i,size) n->next[i] = NULL;
  n->weight = 0;
  n->num = -1;
  return n;
}
/*********************************************************/
/*********************************************************/

void Get_Base_Freqs(calign *data)
{
  int i,j,k;
  phydbl A,C,G,T;
  phydbl fA,fC,fG,fT;
  int w;

  fA = fC = fG = fT = .25;

  For(k,8)
    {
      A = C = G = T = .0;
      For(i,data->n_otu)
	{
	  For(j,data->crunch_len)
	    {
	      w = data->wght[j];
	      if(w)
		{
		  switch(data->c_seq[i]->state[j])
		    {
		    case 'A' : A+=w;
		      break;
		    case 'C' : C+=w;
		      break;
		    case 'G' : G+=w;
		      break;
		    case 'T' : T+=w;
		      break;
		    case 'U' : T+=w;
		      break;
		    case 'M' : C+=w*fC/(fC+fA); A+=w*fA/(fA+fC);
		      break;
		    case 'R' : G+=w*fG/(fA+fG); A+=w*fA/(fA+fG);
		      break;
		    case 'W' : T+=w*fT/(fA+fT); A+=w*fA/(fA+fT);
		      break;
		    case 'S' : C+=w*fC/(fC+fG); G+=w*fG/(fC+fG);
		      break;
		    case 'Y' : C+=w*fC/(fC+fT); T+=w*fT/(fT+fC);
		      break;
		    case 'K' : G+=w*fG/(fG+fT); T+=w*fT/(fT+fG);
		      break;
		    case 'B' : C+=w*fC/(fC+fG+fT); G+=w*fG/(fC+fG+fT); T+=w*fT/(fC+fG+fT);
		      break;
		    case 'D' : A+=w*fA/(fA+fG+fT); G+=w*fG/(fA+fG+fT); T+=w*fT/(fA+fG+fT);
		      break;
		    case 'H' : A+=w*fA/(fA+fC+fT); C+=w*fC/(fA+fC+fT); T+=w*fT/(fA+fC+fT);
		      break;
		    case 'V' : A+=w*fA/(fA+fC+fG); C+=w*fC/(fA+fC+fG); G+=w*fG/(fA+fC+fG);
		      break;
		    case 'N' : case 'X' : case '?' : case 'O' : case '-' :
		      A+=w*fA; C+=w*fC; G+=w*fG; T+=w*fT; break;
		    default : break;
		    }
		}
	    }
	}
      fA = A/(A+C+G+T);
      fC = C/(A+C+G+T);
      fG = G/(A+C+G+T);
      fT = T/(A+C+G+T);
    }
  
  data->b_frq[0] = fA;
  data->b_frq[1] = fC;
  data->b_frq[2] = fG;
  data->b_frq[3] = fT;
}

/*********************************************************/

void Get_AA_Freqs(calign *data)
{
  int i,j,k;
  phydbl A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,Y;
  phydbl fA,fC,fD,fE,fF,fG,fH,fI,fK,fL,fM,fN,fP,fQ,fR,fS,fT,fV,fW,fY;
  int w;
  phydbl sum;
  
  fA = fC = fD = fE = fF = fG = fH = fI = fK = fL =
  fM = fN = fP = fQ = fR = fS = fT = fV = fW = fY = 1./20.;
  
  For(k,8)
    {
      A = C = D = E = F = G = H = I = K = L =
      M = N = P = Q = R = S = T = V = W = Y = .0;

      For(i,data->n_otu)
	{
	  For(j,data->crunch_len)
	    {
	      w = data->wght[j];
	      if(w)
		{
		  switch(data->c_seq[i]->state[j])
		    {
		    case 'A' : A+=w;		break;
		    case 'C' : C+=w;		break;
		    case 'D' : D+=w;		break;
		    case 'E' : E+=w;		break;
		    case 'F' : F+=w;		break;
		    case 'G' : G+=w;		break;
		    case 'H' : H+=w;		break;
		    case 'I' : I+=w;		break;
		    case 'K' : K+=w;		break;
		    case 'L' : L+=w;		break;
		    case 'M' : M+=w;		break;
		    case 'N' : N+=w;		break;
		    case 'P' : P+=w;		break;
		    case 'Q' : Q+=w;		break;
		    case 'R' : R+=w;		break;
		    case 'S' : S+=w;		break;
		    case 'T' : T+=w;		break;
		    case 'V' : V+=w;		break;
		    case 'W' : W+=w;		break;
		    case 'Y' : Y+=w;		break;
		    case 'Z' : Q+=w;		break;
		    case 'X' : case '?' : case 'O' : case '-' :
		      A+=w*fA;
		      C+=w*fC;
		      D+=w*fD;
		      E+=w*fE;
		      F+=w*fF;
		      G+=w*fG;
		      H+=w*fH;
		      I+=w*fI;
		      K+=w*fK;
		      L+=w*fL;
		      M+=w*fM;
		      N+=w*fN;
		      P+=w*fP;
		      Q+=w*fQ;
		      R+=w*fR;
		      S+=w*fS;
		      T+=w*fT;
		      V+=w*fV;
		      W+=w*fW;
		      Y+=w*fY;
		      break;
		    default : break;
		    }
		}
	    }
	}
      sum = (A+C+D+E+F+G+H+I+K+L+M+N+P+Q+R+S+T+V+W+Y);
      fA = A/sum;      fC = C/sum;      fD = D/sum;      fE = E/sum;
      fF = F/sum;      fG = G/sum;      fH = H/sum;      fI = I/sum;
      fK = K/sum;      fL = L/sum;      fM = M/sum;      fN = N/sum;
      fP = P/sum;      fQ = Q/sum;      fR = R/sum;      fS = S/sum;
      fT = T/sum;      fV = V/sum;      fW = W/sum;      fY = Y/sum;
    }

  data->b_frq[0]  = fA;  data->b_frq[1]  = fR;  data->b_frq[2]  = fN;  data->b_frq[3]  = fD;
  data->b_frq[4]  = fC;  data->b_frq[5]  = fQ;  data->b_frq[6]  = fE;  data->b_frq[7]  = fG;
  data->b_frq[8]  = fH;  data->b_frq[9]  = fI;  data->b_frq[10] = fL;  data->b_frq[11] = fK;
  data->b_frq[12] = fM;  data->b_frq[13] = fF;  data->b_frq[14] = fP;  data->b_frq[15] = fS;
  data->b_frq[16] = fT;  data->b_frq[17] = fW;  data->b_frq[18] = fY;  data->b_frq[19] = fV;

}

/*********************************************************/

t_tree *Read_Tree_File(option *io)
{
  t_tree *tree;

  Detect_Tree_File_Format(io);

  io->treelist->list_size = 0;

  switch(io->tree_file_format)
    {
    case PHYLIP: 
      {
	do
	  {
	    io->treelist->tree = (t_tree **)realloc(io->treelist->tree,(io->treelist->list_size+1)*sizeof(t_tree *));
	    io->tree = Read_Tree_File_Phylip(io->fp_in_tree);
	    if(!io->tree) break;
	    if(io->treelist->list_size > 1) PhyML_Printf("\n. Reading tree %d",io->treelist->list_size+1);
	    io->treelist->tree[io->treelist->list_size] = io->tree;
	    io->treelist->list_size++;
	  }while(io->tree);
	break;
      }
    case NEXUS:
      {
	io->nex_com_list = Make_Nexus_Com();
	Init_Nexus_Format(io->nex_com_list);
	Get_Nexus_Data(io->fp_in_tree,io);
	Free_Nexus(io);
	break;
      }
    default:
      {
	PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	Warn_And_Exit("");
	break;
      }
    }
  
  if(!io->long_tax_names)
    {
      int i;

      tree = io->treelist->tree[0];

      io->long_tax_names  = (char **)mCalloc(tree->n_otu,sizeof(char *));
      io->short_tax_names = (char **)mCalloc(tree->n_otu,sizeof(char *));

      For(i,tree->n_otu)
	{
	  io->long_tax_names[i] = (char *)mCalloc(strlen(tree->noeud[i]->name)+1,sizeof(char));
	  io->short_tax_names[i] = (char *)mCalloc(strlen(tree->noeud[i]->name)+1,sizeof(char));
	  strcpy(io->long_tax_names[i],tree->noeud[i]->name);
	  strcpy(io->short_tax_names[i],tree->noeud[i]->name);
	}
    }
  return NULL;
}

/*********************************************************/

t_tree *Read_Tree_File_Phylip(FILE *fp_input_tree)
{
  char *line;
  t_tree *tree;
  int i;
  char c;

 
  do
    {
      c=fgetc(fp_input_tree);
    }
  while((c != '(') && (c != EOF));
  
  if(c==EOF) return NULL;

  line = (char *)mCalloc(1,sizeof(char));
  
  i=0;
  for(;;)
    {
      if((c == ' ') || (c == '\n'))
	{
	  c=fgetc(fp_input_tree);
	  if(c == EOF || c == ';') break;
	  else continue;
	}
      
      if(c == '[')
	{
	  Skip_Comment(fp_input_tree);
	  c = fgetc(fp_input_tree);
	  if(c == EOF || c == ';') break;
	}

      line = (char *)mRealloc(line,i+2,sizeof(char));

      line[i]=c;
      i++;
      c=fgetc(fp_input_tree);
      if(c==EOF || c==';') break;
    }
  line[i] = '\0';

  tree = Read_Tree(line);
  Free(line);
  return tree;
}

/*********************************************************/

void Connect_Edges_To_Nodes_Recur(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  Connect_One_Edge_To_Two_Nodes(a,d,tree->t_edges[tree->num_curr_branch_available],tree);
  tree->num_curr_branch_available += 1;

  if(d->tax) return;
  else For(i,3) if(d->v[i] != a) Connect_Edges_To_Nodes_Recur(d,d->v[i],tree);
}

/*********************************************************/

void Connect_One_Edge_To_Two_Nodes(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i,dir_a_d,dir_d_a;

  dir_a_d = -1;
  For(i,3) if(a->v[i] == d) {dir_a_d = i; break;}

  dir_d_a = -1;
  For(i,3) if(d->v[i] == a) {dir_d_a = i; break;}

  if(dir_a_d == -1 || dir_d_a == -1)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  a->b[dir_a_d] = b;
  d->b[dir_d_a] = b;
  b->num        = tree->num_curr_branch_available;
  b->left       = a;
  b->rght       = d;
  if(a->tax) {b->rght = a; b->left = d;} /* root */
  /* a tip is necessary on the right hand side of the t_edge */
  
  (b->left == a)?
    (Make_Edge_Dirs(b,a,d,tree)):
    (Make_Edge_Dirs(b,d,a,tree));

  b->l                    = a->l[b->l_r];
  if(a->tax) b->l         = a->l[b->r_l];
  b->l_old                = b->l;
}

/*********************************************************/

void Update_Dirs(t_tree *tree)
{
  int i;
  int buff;
  t_edge *b;

  b = NULL;
  buff = -1;
  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      
      if((!b->left->tax) && (b->left->v[b->l_v1]->num < b->left->v[b->l_v2]->num))
	{
	  buff    = b->l_v1;
	  b->l_v1 = b->l_v2;
	  b->l_v2 = buff;
	}
      if((!b->rght->tax) && (b->rght->v[b->r_v1]->num < b->rght->v[b->r_v2]->num))
	{
	  buff    = b->r_v1;
	  b->r_v1 = b->r_v2;
	  b->r_v2 = buff;
	}
    }

}

/*********************************************************/

void Exit(char *message)
{
  fflush(NULL);
  PhyML_Fprintf(stderr,"%s",message);
  exit(1);
}

/*********************************************************/

void *mCalloc(int nb, size_t size)
{
  void *allocated;

  if((allocated = calloc((size_t)nb,size)) != NULL)
/*   if((allocated = malloc((size_t)nb*(size_t)size)) != NULL) */
    {
      return allocated;
    }
  else
    Warn_And_Exit("\n. Err: low memory\n");

  return NULL;
}

/*********************************************************/

void *mRealloc(void *p,int nb, size_t size)
{
  if((p = realloc(p,(size_t)nb*size)) != NULL)
	return p;
  else
    Warn_And_Exit("\n. Err: low memory\n");

  return NULL;
}

/*********************************************************/

/* t_tree *Make_Light_Tree_Struct(int n_otu) */
/* { */
/*   t_tree *tree; */
/*   int i; */

/*   tree          = (t_tree *)mCalloc(1,sizeof(t_tree )); */
/*   tree->t_edges = (t_edge **)mCalloc(2*n_otu-3,sizeof(t_edge *)); */
/*   tree->noeud   = (t_node **)mCalloc(2*n_otu-2,sizeof(t_node *)); */
/*   tree->n_otu   = n_otu; */

/*   For(i,2*n_otu-3) */
/*     tree->t_edges[i] = Make_Edge_Light(NULL,NULL,i); */

/*   For(i,2*n_otu-2) */
/*     tree->noeud[i] = Make_Node_Light(i); */

/*   return tree; */
/* } */

/*********************************************************/

int Sort_Phydbl_Decrease(const void *a, const void *b)
{
    if((*(phydbl *)(a)) >= (*(phydbl *)(b))) return -1;
    else return 1;
}

/*********************************************************/

void Qksort_Int(int *A, int *B, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    int tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
      if (A[uhi] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo]    = A[uhi];
	    A[uhi]    = tempEntry;

	    if(B)
	      {
		tempEntry = B[ulo];
		B[ulo]    = B[uhi];
		B[uhi]    = tempEntry;
	      }



	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;


		if(B)
		  {
		    tempEntry = B[ieq];
		    B[ieq] = B[ulo];
		    B[ulo] = tempEntry;
		  }


		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort_Int(A, B, ilo, ieq - 1);
    Qksort_Int(A, B, uhi + 1, ihi);

}

/*********************************************************/

/* Sort in ascending order. Elements in B (if provided) are also re-ordered according to the ordering of A  */
void Qksort(phydbl *A, phydbl *B, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl tempEntry;	// temporary entry used for swapping

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
      if (A[uhi] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo]    = A[uhi];
	    A[uhi]    = tempEntry;

	    if(B)
	      {
		tempEntry = B[ulo];
		B[ulo]    = B[uhi];
		B[uhi]    = tempEntry;
	      }



	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;


		if(B)
		  {
		    tempEntry = B[ieq];
		    B[ieq] = B[ulo];
		    B[ulo] = tempEntry;
		  }


		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort(A, B, ilo, ieq - 1);
    Qksort(A, B, uhi + 1, ihi);
}

/********************************************************/

void Qksort_Matrix(phydbl **A, int col, int ilo, int ihi)
{
    phydbl pivot;	// pivot value for partitioning array
    int ulo, uhi;	// indices at ends of unpartitioned region
    int ieq;		// least index of array entry with value equal to pivot
    phydbl *tempEntry;	// temporary entry used for swapping

    tempEntry = NULL;

    if (ilo >= ihi) {
	return;
    }
    // Select a pivot value.
    pivot = A[(ilo + ihi)/2][col];
    // Initialize ends of unpartitioned region and least index of entry
    // with value equal to pivot.
    ieq = ulo = ilo;
    uhi = ihi;
    // While the unpartitioned region is not empty, try to reduce its size.
    while (ulo <= uhi) {
	if (A[uhi][col] > pivot) {
	    // Here, we can reduce the size of the unpartitioned region and
	    // try again.
	    uhi--;
	} else {
	    // Here, A[uhi] <= pivot, so swap entries at indices ulo and
	    // uhi.
	    tempEntry = A[ulo];
	    A[ulo] = A[uhi];
	    A[uhi] = tempEntry;
	    // After the swap, A[ulo] <= pivot.
	    if (A[ulo][col] < pivot) {
		// Swap entries at indices ieq and ulo.
		tempEntry = A[ieq];
		A[ieq] = A[ulo];
		A[ulo] = tempEntry;
		// After the swap, A[ieq] < pivot, so we need to change
		// ieq.
		ieq++;
		// We also need to change ulo, but we also need to do
		// that when A[ulo] = pivot, so we do it after this if
		// statement.
	    }
	    // Once again, we can reduce the size of the unpartitioned
	    // region and try again.
	    ulo++;
	}
    }
    // Now, all entries from index ilo to ieq - 1 are less than the pivot
    // and all entries from index uhi to ihi + 1 are greater than the
    // pivot.  So we have two regions of the array that can be sorted
    // recursively to put all of the entries in order.
    Qksort_Matrix(A, col, ilo, ieq - 1);
    Qksort_Matrix(A, col, uhi + 1, ihi);
}

/********************************************************/

void Print_Site(calign *cdata, int num, int n_otu, char *sep, int stepsize)
{
  int i,j;
  For(i,n_otu)
    {
      PhyML_Printf("%s   ",cdata->c_seq[i]->name);
      For(j,stepsize)
	PhyML_Printf("%c",cdata->c_seq[i]->state[num+j]);
      PhyML_Printf("%s",sep);
    }
  PhyML_Fprintf(stderr,"%s",sep);
}

/*********************************************************/

void Print_Site_Lk(t_tree *tree, FILE *fp)
{
  int site;
  int catg;
  char *s;
  phydbl postmean;

  if(!tree->io->print_site_lnl)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(!tree->io->print_trace)
    {
      s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
      
      PhyML_Fprintf(fp,"Note : P(D|M) is the probability of site D given the model M (i.e., the site likelihood)\n");
      if(tree->mod->n_catg > 1 || tree->mod->invar)
	PhyML_Fprintf(fp,"P(D|M,rr[x]) is the probability of site D given the model M and the relative rate\nof evolution rr[x], where x is the class of rate to be considered.\nWe have P(D|M) = \\sum_x P(x) x P(D|M,rr[x]).\n");
      PhyML_Fprintf(fp,"\n\n");
      
      sprintf(s,"Site");
      PhyML_Fprintf(fp, "%-7s",s);
      
      sprintf(s,"P(D|M)");
      PhyML_Fprintf(fp,"%-16s",s);
      
      if(tree->mod->n_catg > 1)
	{
	  For(catg,tree->mod->n_catg)
	    {
	      sprintf(s,"P(D|M,rr[%d]=%5.4f)",catg+1,tree->mod->gamma_rr[catg]);
	      PhyML_Fprintf(fp,"%-22s",s);
	    }
	  
	  sprintf(s,"Posterior mean");
	  PhyML_Fprintf(fp,"%-22s",s);
	}
      
      
      if(tree->mod->invar)
	{
	  sprintf(s,"P(D|M,rr[0]=0)");
	  PhyML_Fprintf(fp,"%-16s",s);
	}
      PhyML_Fprintf(fp,"\n");
      
      For(site,tree->data->init_len)
	{
	  PhyML_Fprintf(fp,"%-7d",site+1);
	  PhyML_Fprintf(fp,"%-16g",(phydbl)EXP(tree->cur_site_lk[tree->data->sitepatt[site]]));      
	  if(tree->mod->n_catg > 1)
	    {
	      For(catg,tree->mod->n_catg)
		PhyML_Fprintf(fp,"%-22g",(phydbl)EXP(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]));

	      postmean = .0;
	      For(catg,tree->mod->n_catg) 
		postmean += 
		tree->mod->gamma_rr[catg] * 
		EXP(tree->log_site_lk_cat[catg][tree->data->sitepatt[site]]) * 
		tree->mod->gamma_r_proba[catg];
	      postmean /= EXP(tree->cur_site_lk[tree->data->sitepatt[site]]);

	      PhyML_Fprintf(fp,"%-22g",postmean);
	    }
	  if(tree->mod->invar)
	    {
	      if((phydbl)tree->data->invar[tree->data->sitepatt[site]] > -0.5)
		PhyML_Fprintf(fp,"%-16g",tree->mod->pi[tree->data->invar[tree->data->sitepatt[site]]]);
	      else
		PhyML_Fprintf(fp,"%-16g",0.0);
	    }
	  PhyML_Fprintf(fp,"\n");
	}
      Free(s);
    }
  else
    {
      For(site,tree->data->init_len)
	PhyML_Fprintf(fp,"%.2f\t",tree->cur_site_lk[tree->data->sitepatt[site]]);
      PhyML_Fprintf(fp,"\n");
    }
}


/*********************************************************/

void Print_Seq(align **data, int n_otu)
{
  int i,j;

  PhyML_Printf("%d\t%d\n",n_otu,data[0]->len);
  For(i,n_otu)
    {
      For(j,20)
	{
	  if(j<(int)strlen(data[i]->name))
	     putchar(data[i]->name[j]);
	  else putchar(' ');
	}
/*       PhyML_Printf("%10d  ",i); */
      For(j,data[i]->len)
	{
	  PhyML_Printf("%c",data[i]->state[j]);
	}
      PhyML_Printf("\n");
    }
}

/*********************************************************/

void Print_CSeq(FILE *fp, calign *cdata)
{
  int i,j;
  int n_otu;
  
  n_otu = cdata->n_otu;
  if(cdata->format == 0)
    {
      PhyML_Fprintf(fp,"%d\t%d\n",n_otu,cdata->init_len);
    }
  else
    {
      PhyML_Fprintf(fp,"#NEXUS\n");
      PhyML_Fprintf(fp,"begin data\n");
      PhyML_Fprintf(fp,"dimensions ntax=%d nchar=%d;\n",n_otu,cdata->init_len);
      PhyML_Fprintf(fp,"format sequential datatype=dna;\n");
      PhyML_Fprintf(fp,"matrix\n");
    }
  For(i,n_otu)
    {
      For(j,50)
	{
	  if(j<(int)strlen(cdata->c_seq[i]->name))
	    fputc(cdata->c_seq[i]->name[j],fp);
	  else fputc(' ',fp);
	}
      
      PhyML_Fprintf(fp,"%s",cdata->c_seq[i]->state);
      PhyML_Fprintf(fp,"\n");
    }
  PhyML_Fprintf(fp,"\n");

  if(cdata->format == 1)
    {
      PhyML_Fprintf(fp,";\n");
      PhyML_Fprintf(fp,"END;\n");
    }


/*   PhyML_Printf("\t"); */
/*   For(j,cdata->crunch_len) */
/*     PhyML_Printf("%.0f ",cdata->wght[j]); */
/*   PhyML_Printf("\n"); */
}

/*********************************************************/

void Order_Tree_Seq(t_tree *tree, align **data)
{
    int i,j,n_otu;
    align *buff;

    n_otu = tree->n_otu;

    For(i,n_otu)
      {
	For(j,n_otu)
	  {
	    if(!strcmp(tree->noeud[i]->name,data[j]->name))
	      break;
	  }
	buff = data[j];
	data[j] = data[i];
	data[i] = buff;
      }
}

/*********************************************************/

void Order_Tree_CSeq(t_tree *tree, calign *cdata)
{
    int i,j,n_otu_tree,n_otu_cdata;
    align *buff;

    n_otu_tree  = tree->n_otu;
    n_otu_cdata = cdata->n_otu;

    if(n_otu_tree != n_otu_cdata) 
      {
	PhyML_Printf("\n. Number of taxa in the tree: %d, number of sequences: %d.\n",n_otu_tree,n_otu_cdata);
	Warn_And_Exit("\n. The number of tips in the tree is not the same as the number of sequences\n");
      }

    For(i,MAX(n_otu_tree,n_otu_cdata))
      {
	For(j,MIN(n_otu_tree,n_otu_cdata))
	  {
	    if(!strcmp(tree->noeud[i]->name,cdata->c_seq[j]->name))
	      break;
	  }
	
	if(j==MIN(n_otu_tree,n_otu_cdata))
	  {
	    PhyML_Printf("\n. Err: %s is not found in sequence data set\n",tree->noeud[i]->name);
	    Warn_And_Exit("");
	  }
	
	buff            = cdata->c_seq[j];
	cdata->c_seq[j] = cdata->c_seq[i];
	cdata->c_seq[i] = buff;
      }
}

/*********************************************************/

matrix *Make_Mat(int n_otu)
{
  matrix *mat;
  int i;

  mat = (matrix *)mCalloc(1,sizeof(matrix));

  mat->n_otu = n_otu;

  mat->P        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->Q        = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->dist     = (phydbl **)mCalloc(n_otu,sizeof(phydbl *));
  mat->on_off   = (int *)mCalloc(n_otu,sizeof(int));
  mat->name     = (char **)mCalloc(n_otu,sizeof(char *));
  mat->tip_node = (t_node **)mCalloc(n_otu,sizeof(t_node *));


  For(i,n_otu)
    {
      mat->P[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->Q[i]    = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->dist[i] = (phydbl *)mCalloc(n_otu,sizeof(phydbl));
      mat->name[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));
    }

  return mat;
}

/*********************************************************/

void Init_Mat(matrix *mat, calign *data)
{
  int i;

  mat->n_otu = data->n_otu;
  mat->r = mat->n_otu;
  mat->curr_int = mat->n_otu;
  mat->method = 1;

  For(i,data->n_otu)
    {
      strcpy(mat->name[i],data->c_seq[i]->name);
      mat->on_off[i] = 1;
    }
}

/*********************************************************/

t_tree *Make_Tree_From_Scratch(int n_otu, calign *data)
{
  t_tree *tree;

  tree = Make_Tree(n_otu);
  Init_Tree(tree,n_otu);
  Make_All_Tree_Nodes(tree);
  Make_All_Tree_Edges(tree);
  Make_Tree_Path(tree);
  if(data)
    {
      Copy_Tax_Names_To_Tip_Labels(tree,data);
      tree->data = data;
    }
  return tree;
}

/*********************************************************/

t_tree *Make_Tree(int n_otu)
{
  t_tree *tree;
  tree = (t_tree *)mCalloc(1,sizeof(t_tree ));
  tree->t_dir = (short int *)mCalloc((2*n_otu-2)*(2*n_otu-2),sizeof(int));
  return tree;
}

/*********************************************************/

void Make_Tree_Path(t_tree *tree)
{
  tree->curr_path = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *));
}

/*********************************************************/

void Make_All_Tree_Nodes(t_tree *tree)
{
  int i;

  tree->noeud = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));

  For(i,2*tree->n_otu-1)
    {
      tree->noeud[i] = (t_node *)Make_Node_Light(i);
      if(i < tree->n_otu) tree->noeud[i]->tax = 1;
      else                tree->noeud[i]->tax = 0;
    }
}

/*********************************************************/

void Make_All_Tree_Edges(t_tree *tree)
{
  int i;
  tree->t_edges      = (t_edge **)mCalloc(2*tree->n_otu-2,sizeof(t_edge *));
  For(i,2*tree->n_otu-2) tree->t_edges[i] = (t_edge *)Make_Edge_Light(NULL,NULL,i);
}

/*********************************************************/

void Copy_Tax_Names_To_Tip_Labels(t_tree *tree, calign *data)
{
  int i;

  For(i,tree->n_otu)
    {
      tree->noeud[i]->name = (char *)mCalloc((int)strlen(data->c_seq[i]->name)+1,sizeof(char));
      tree->noeud[i]->ori_name = tree->noeud[i]->name;
      strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
      tree->noeud[i]->tax = 1;
      tree->noeud[i]->num = i;
    }
}

/*********************************************************/

void Print_Dist(matrix *mat)
{
  int i,j;

  For(i,mat->n_otu)
    {
      PhyML_Printf("%s ",mat->name[i]);

      For(j,mat->n_otu)
	PhyML_Printf("%9.6f ",mat->dist[i][j]);
      PhyML_Printf("\n");
    }
}

/*********************************************************/

void Print_Node(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  int dir;
  dir = -1;
  For(i,3) if(a->v[i] == d) {dir = i; break;}
  PhyML_Printf("Node nums = %3d %3d  (dir=%d);",a->num,d->num,dir);
  PhyML_Printf("Node names = '%s' '%s' ; ",a->name,d->name);
  For(i,3) if(a->v[i] == d)
    {
      if(a->b[i])
	{
	  PhyML_Printf("Branch num = %3d%c (%d %d) %f",
		       a->b[i]->num,a->b[i]==tree->e_root?'*':' ',a->b[i]->left->num,
		       a->b[i]->rght->num,a->b[i]->l);
	  if(a->b[i]->left->tax) PhyML_Printf(" WARNING LEFT->TAX!");
	  break;
	}
    }
  PhyML_Printf("\n");

  if(d->tax) return;
  else
    For(i,3)
      if(d->v[i] != a) Print_Node(d,d->v[i],tree);
}

/*********************************************************/

void Share_Lk_Struct(t_tree *t_full, t_tree *t_empt)
{
  int i,j,n_otu;
  t_edge *b_e,*b_f;
  t_node *n_e, *n_f;

  n_otu                   = t_full->n_otu;
  t_empt->n_root          = t_full->n_root;
  t_empt->e_root          = t_full->e_root;
  t_empt->c_lnL_sorted    = t_full->c_lnL_sorted;
  t_empt->log_site_lk_cat = t_full->log_site_lk_cat;
  t_empt->cur_site_lk     = t_full->cur_site_lk;
  t_empt->old_site_lk     = t_full->old_site_lk;
  t_empt->triplet_struct  = t_full->triplet_struct;
  t_empt->log_lks_aLRT    = t_full->log_lks_aLRT;
  t_empt->site_lk_cat     = t_full->site_lk_cat;

  For(i,2*n_otu-3)
    {
      b_f = t_full->t_edges[i];
      b_e = t_empt->t_edges[i];

      b_e->Pij_rr = b_f->Pij_rr;

      b_e->nni = b_f->nni;
    }


  for(i=n_otu;i<2*n_otu-2;i++)
    {
      n_f = t_full->noeud[i];
      n_e = t_empt->noeud[i];
            
      For(j,3)
	{
	  if(n_f->b[j]->left == n_f)
	    {
	      if(n_e->b[j]->left == n_e)
		{
		  n_e->b[j]->p_lk_left          = n_f->b[j]->p_lk_left;
		  n_e->b[j]->p_lk_loc_left      = n_f->b[j]->p_lk_loc_left;
		  n_e->b[j]->patt_id_left       = n_f->b[j]->patt_id_left;
		  n_e->b[j]->sum_scale_left     = n_f->b[j]->sum_scale_left;
		  n_e->b[j]->sum_scale_left_cat = n_f->b[j]->sum_scale_left_cat;
		  n_e->b[j]->p_lk_tip_l         = n_f->b[j]->p_lk_tip_l;
		}
	      else
		{
		  n_e->b[j]->p_lk_rght          = n_f->b[j]->p_lk_left;
		  n_e->b[j]->p_lk_loc_rght      = n_f->b[j]->p_lk_loc_left;
		  n_e->b[j]->patt_id_rght       = n_f->b[j]->patt_id_left;
		  n_e->b[j]->sum_scale_rght     = n_f->b[j]->sum_scale_left;
		  n_e->b[j]->sum_scale_rght_cat = n_f->b[j]->sum_scale_left_cat;
		  n_e->b[j]->p_lk_tip_r         = n_f->b[j]->p_lk_tip_l;
		}
	    }
	  else
	    {
	      if(n_e->b[j]->rght == n_e)
		{
		  n_e->b[j]->p_lk_rght          = n_f->b[j]->p_lk_rght;
		  n_e->b[j]->p_lk_loc_rght      = n_f->b[j]->p_lk_loc_rght;
		  n_e->b[j]->patt_id_rght       = n_f->b[j]->patt_id_rght;
		  n_e->b[j]->sum_scale_rght     = n_f->b[j]->sum_scale_rght;
		  n_e->b[j]->sum_scale_rght_cat = n_f->b[j]->sum_scale_rght_cat;
		  n_e->b[j]->p_lk_tip_r         = n_f->b[j]->p_lk_tip_r;
		}
	      else
		{
		  n_e->b[j]->p_lk_left          = n_f->b[j]->p_lk_rght;
		  n_e->b[j]->p_lk_loc_left      = n_f->b[j]->p_lk_loc_rght;
		  n_e->b[j]->patt_id_left       = n_f->b[j]->patt_id_rght;
		  n_e->b[j]->sum_scale_left     = n_f->b[j]->sum_scale_rght;
		  n_e->b[j]->sum_scale_left_cat = n_f->b[j]->sum_scale_rght_cat;
		  n_e->b[j]->p_lk_tip_l         = n_f->b[j]->p_lk_tip_r;
		}
	    }
	}
    }

  For(i,n_otu)
    {
      n_f = t_full->noeud[i];
      n_e = t_empt->noeud[i];

      if(n_f->b[0]->rght == n_f)
	{
	  n_e->b[0]->p_lk_rght          = n_f->b[0]->p_lk_rght;
	  n_e->b[0]->p_lk_loc_rght      = n_f->b[0]->p_lk_loc_rght;
	  n_e->b[0]->patt_id_rght       = n_f->b[0]->patt_id_rght;
	  n_e->b[0]->sum_scale_rght     = n_f->b[0]->sum_scale_rght;
	  n_e->b[0]->sum_scale_rght_cat = n_f->b[0]->sum_scale_rght_cat;
	  n_e->b[0]->p_lk_tip_r         = n_f->b[0]->p_lk_tip_r;
	}
      else
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
    }
}

/*********************************************************/

void Share_Spr_Struct(t_tree *t_full, t_tree *t_empt)
{
  t_empt->size_spr_list = t_full->size_spr_list;
  t_empt->spr_list      = t_full->spr_list;
  t_empt->best_spr      = t_full->best_spr;
}

/*********************************************************/

void Share_Pars_Struct(t_tree *t_full, t_tree *t_empt)
{
  int i;

  t_empt->site_pars = t_full->site_pars;
  t_empt->step_mat  = t_full->step_mat;

  For(i,2*t_full->n_otu-3)
    {
      t_empt->t_edges[i]->ui_l     = t_full->t_edges[i]->ui_l;
      t_empt->t_edges[i]->ui_r     = t_full->t_edges[i]->ui_r;

      t_empt->t_edges[i]->pars_l   = t_full->t_edges[i]->pars_l;
      t_empt->t_edges[i]->pars_r   = t_full->t_edges[i]->pars_r;

      t_empt->t_edges[i]->p_pars_l = t_full->t_edges[i]->p_pars_l;
      t_empt->t_edges[i]->p_pars_r = t_full->t_edges[i]->p_pars_r;
    }
}

/*********************************************************/

void Print_Mat(matrix *mat)
{
  int i,j;

  PhyML_Printf("%d",mat->n_otu);
  PhyML_Printf("\n");

  For(i,mat->n_otu)
    {
      For(j,13)
	{
	  if(j>=(int)strlen(mat->name[i])) putchar(' ');
	  else putchar(mat->name[i][j]);
	}

      For(j,mat->n_otu)
	{
	  char s[2]="-";
	  if(mat->dist[i][j] < .0)
	    PhyML_Printf("%12s",s);
	  else
	    PhyML_Printf("%12f",mat->dist[i][j]);
	}
      PhyML_Printf("\n");
    }
}

/*********************************************************/

int Sort_Edges_NNI_Score(t_tree *tree, t_edge **sorted_edges, int n_elem)
{
  int i,j;
  t_edge *buff;

  For(i,n_elem-1)
    {
      for(j=i+1;j<n_elem;j++)
	{
	  if(sorted_edges[j]->nni->score  < sorted_edges[i]->nni->score)
	    {
	      buff = sorted_edges[j];
	      sorted_edges[j] = sorted_edges[i];
	      sorted_edges[i] = buff;
	    }
	}
    }
  return 1;
}

/*********************************************************/

int Sort_Edges_Depth(t_tree *tree, t_edge **sorted_edges, int n_elem)
{
  int i,j;
  t_edge *buff;
  phydbl *depth,buff_depth;
  
  depth = (phydbl *)mCalloc(n_elem,sizeof(phydbl));

  For(i,n_elem) 
    depth[i] = 
    sorted_edges[i]->left->bip_size[sorted_edges[i]->l_r] * 
    sorted_edges[i]->rght->bip_size[sorted_edges[i]->r_l] ;


  For(i,n_elem-1)
    {
      for(j=i+1;j<n_elem;j++)
	{
	  if(depth[i] > depth[j])
	    {
	      buff = sorted_edges[i];
	      sorted_edges[i] = sorted_edges[j];
	      sorted_edges[j] = buff;

	      buff_depth = depth[i];
	      depth[i] = depth[j];
	      depth[j] = buff_depth;
	    }
	}
    }

  Free(depth);

  return 1;
}

/*********************************************************/

void NNI(t_tree *tree, t_edge *b_fcus, int do_swap)
{
  int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
  t_node *v1,*v2,*v3,*v4;
  phydbl lk0, lk1, lk2;
  phydbl lk0_init, lk1_init, lk2_init;
  phydbl bl_init;
  phydbl l0,l1,l2;
  phydbl l_infa, l_infb, l_max;
/*   phydbl lk_infa, lk_infb, lk_max; */
  phydbl lk_init;

  bl_init                = b_fcus->l;
  lk_init                = tree->c_lnL;

  b_fcus->nni->init_l    = b_fcus->l;
  b_fcus->nni->init_lk   = tree->c_lnL;;

  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;

  lk0 = lk1 = lk2        = UNLIKELY;
  v1 = v2 = v3 = v4      = NULL;

  l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

  l_r                    = b_fcus->l_r;
  r_l                    = b_fcus->r_l;

  v1                     = b_fcus->left->v[b_fcus->l_v1];
  v2                     = b_fcus->left->v[b_fcus->l_v2];
  v3                     = b_fcus->rght->v[b_fcus->r_v1];
  v4                     = b_fcus->rght->v[b_fcus->r_v2];

  Record_Br_Len(tree);

  if(v1->num < v2->num)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(v3->num < v4->num)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  l0 = l1 = l2 = -1.;

  
  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  tree->both_sides = 1;

  tree->update_alias_subpatt = YES;
  lk1_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  tree->update_alias_subpatt = NO;

  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = tree->mod->l_min;

  if(tree->mod->s_opt->fast_nni)
    {
      Fast_Br_Len(b_fcus,tree,1);
      lk1 = Lk_At_Given_Edge(b_fcus,tree);
    }
  else
    {
      lk1 = Br_Len_Brent(l_infa,l_max,l_infb,
			 tree->mod->s_opt->min_diff_lk_local,
			 b_fcus,tree,
			 tree->mod->s_opt->brent_it_max,
			 tree->mod->s_opt->quickdirty);
    }

  if(lk1 < lk1_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("%f %f %f %G\n",l_infa,l_max,l_infb,b_fcus->l);
      PhyML_Printf("%f -- %f \n",lk1_init,lk1);
      PhyML_Printf("\n. Err. in NNI (1)\n");
    }

  l1  = b_fcus->l;
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/


  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  Restore_Br_Len(tree);
  /* b_fcus->l = bl_init; */
  tree->both_sides = 1;

  tree->update_alias_subpatt = YES;
  lk2_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  tree->update_alias_subpatt = NO;

  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = tree->mod->l_min;

  if(tree->mod->s_opt->fast_nni)
    {
      Fast_Br_Len(b_fcus,tree,1);
      lk2 = Lk_At_Given_Edge(b_fcus,tree);
    }
  else
    {
      lk2 = Br_Len_Brent(l_infa,l_max,l_infb,
			 tree->mod->s_opt->min_diff_lk_local,
			 b_fcus,tree,
			 tree->mod->s_opt->brent_it_max,
			 tree->mod->s_opt->quickdirty);
    }

  if(lk2 < lk2_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("%f %f %f %G\n",l_infa,l_max,l_infb,b_fcus->l);
      PhyML_Printf("%f -- %f \n",lk2_init,lk2);
      PhyML_Printf("\n. Err. in NNI (2)\n");
   }

  l2  = b_fcus->l;
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/


  
  /***********/
  Restore_Br_Len(tree);
  /* b_fcus->l = bl_init; */
  tree->both_sides = 1;

  tree->update_alias_subpatt = YES;
  lk0_init = Update_Lk_At_Given_Edge(b_fcus,tree);
  tree->update_alias_subpatt = NO;

  if(FABS(lk0_init - lk_init) > tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. lk_init = %f; lk = %f diff = %f l = %G\n",
		   lk_init,
		   lk0_init,
		   lk_init-lk0_init,
		   b_fcus->l);
      PhyML_Printf("\n. Curr_lnL = %f\n",Lk(tree));
      Warn_And_Exit("\n. Err. in NNI (3)\n");
    }
  
  l_infa = 10.*b_fcus->l;
  l_max  = b_fcus->l;
  l_infb = tree->mod->l_min;
  
  if(tree->mod->s_opt->fast_nni)
    {
      Fast_Br_Len(b_fcus,tree,1);
      lk0 = Lk_At_Given_Edge(b_fcus,tree);
    }
  else
    {
      lk0 = Br_Len_Brent(l_infa,l_max,l_infb,
			 tree->mod->s_opt->min_diff_lk_local,
			 b_fcus,tree,
			 tree->mod->s_opt->brent_it_max,
			 tree->mod->s_opt->quickdirty);
    }
  
  if(lk0 < lk_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. %f %f %f %f",l_infa,l_max,l_infb,b_fcus->l);
      PhyML_Printf("\n. %f -- %f",lk0_init,lk0);
      PhyML_Printf("\n. Err. in NNI (3)\n");
      Warn_And_Exit("\n");
    }
  
  l0  = b_fcus->l;
  /***********/
  
  b_fcus->nni->lk0 = lk0;
  b_fcus->nni->lk1 = lk1;
  b_fcus->nni->lk2 = lk2;
  
  b_fcus->nni->l0  = l0;
  b_fcus->nni->l1  = l1;
  b_fcus->nni->l2  = l2;
  
  b_fcus->nni->score = lk0 - MAX(lk1,lk2);
  
  if((b_fcus->nni->score <  tree->mod->s_opt->min_diff_lk_local) &&
     (b_fcus->nni->score > -tree->mod->s_opt->min_diff_lk_local))
    {
      b_fcus->nni->score = .0;
      b_fcus->nni->lk1 = b_fcus->nni->lk0;
      b_fcus->nni->lk2 = b_fcus->nni->lk0;
    }
  
  if(lk0 > MAX(lk1,lk2))
    {
      b_fcus->nni->best_conf    = 0;
      b_fcus->nni->best_l       = l0;
      b_fcus->nni->swap_node_v1 = NULL;
      b_fcus->nni->swap_node_v2 = NULL;
      b_fcus->nni->swap_node_v3 = NULL;
      b_fcus->nni->swap_node_v4 = NULL;
    }
  else if(lk1 > MAX(lk0,lk2))
    {
      b_fcus->nni->best_conf    = 1;
      b_fcus->nni->best_l       = l1;
      b_fcus->nni->swap_node_v1 = v2;
      b_fcus->nni->swap_node_v2 = b_fcus->left;
      b_fcus->nni->swap_node_v3 = b_fcus->rght;
      b_fcus->nni->swap_node_v4 = v3;
    }
  else if(lk2 > MAX(lk0,lk1))
    {
      b_fcus->nni->best_conf    = 2;
      b_fcus->nni->best_l       = l2;
      b_fcus->nni->swap_node_v1 = v2;
      b_fcus->nni->swap_node_v2 = b_fcus->left;
      b_fcus->nni->swap_node_v3 = b_fcus->rght;
      b_fcus->nni->swap_node_v4 = v4;
    }
  else
    {
      b_fcus->nni->score        = +1.0;
      b_fcus->nni->best_conf    = 0;
      b_fcus->nni->best_l       = l0;
      b_fcus->nni->swap_node_v1 = NULL;
      b_fcus->nni->swap_node_v2 = NULL;
      b_fcus->nni->swap_node_v3 = NULL;
      b_fcus->nni->swap_node_v4 = NULL;
    }
  
  if((do_swap) && ((lk1 > lk0) || (lk2 > lk0)))
    {
      tree->n_swap++;
      PhyML_Printf("Swap t_edge %d -> %f\n",b_fcus->num,MAX(lk1,lk2));
      
      if(lk1 > lk2)
	{
	  tree->best_lnL = lk1;
	  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
	  b_fcus->l = l1;
	  tree->both_sides = 1;
	  Lk(tree);
	}
      else
	{
	  tree->best_lnL = lk2;
	  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
	  b_fcus->l = l2;
	  tree->both_sides = 1;
	  Lk(tree);
	}
    }
  else
    {
      /* b_fcus->l = bl_init; */
      Restore_Br_Len(tree);
      Update_PMat_At_Given_Edge(b_fcus,tree);
      tree->c_lnL = lk_init;
    }
}

/*********************************************************/

void NNI_Pars(t_tree *tree, t_edge *b_fcus, int do_swap)
{
  int l_r, r_l, l_v1, l_v2, r_v3, r_v4;
  t_node *v1,*v2,*v3,*v4;
  int pars0, pars1, pars2;
  int pars_init;

  pars_init              = tree->c_pars;
  b_fcus->nni->best_conf = 0;
  b_fcus->nni->score     = +1.0;

  pars0 = pars1 = pars2  = 0;
  v1 = v2 = v3 = v4      = NULL;

  l_r = r_l = l_v1 = l_v2 = r_v3 = r_v4 = -1;

  l_r                    = b_fcus->l_r;
  r_l                    = b_fcus->r_l;

  v1                     = b_fcus->left->v[b_fcus->l_v1];
  v2                     = b_fcus->left->v[b_fcus->l_v2];
  v3                     = b_fcus->rght->v[b_fcus->r_v1];
  v4                     = b_fcus->rght->v[b_fcus->r_v2];

  if(v1->num < v2->num)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(v3->num < v4->num)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  
  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v3,tree);
  tree->both_sides = 1;
  pars1 = Update_Pars_At_Given_Edge(b_fcus,tree);
  Swap(v3,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/

  /***********/
  Swap(v2,b_fcus->left,b_fcus->rght,v4,tree);
  tree->both_sides = 1;
  pars2 = Update_Pars_At_Given_Edge(b_fcus,tree);
  Swap(v4,b_fcus->left,b_fcus->rght,v2,tree);
  /***********/


  /***********/
   tree->both_sides = 1;
   pars0 = Update_Pars_At_Given_Edge(b_fcus,tree);
 
   if(pars0 != pars_init)
     {
       PhyML_Printf("\n. pars_init = %d; pars0 = %d\n",
	      pars_init,
	      pars0);
       Warn_And_Exit("\n. Err. in NNI (3)\n");
     }
   /***********/

   tree->c_pars = pars0;

   b_fcus->nni->score = MIN(pars1,pars2) - pars0;

   if(pars0 < MIN(pars1,pars2))
     {
       b_fcus->nni->best_conf    = 0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
      }
   else if(pars1 < MIN(pars0,pars2))
     {
       b_fcus->nni->best_conf    = 1;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v3;
     }
   else if(pars2 > MIN(pars0,pars1))
     {
       b_fcus->nni->best_conf    = 2;
       b_fcus->nni->swap_node_v1 = v2;
       b_fcus->nni->swap_node_v2 = b_fcus->left;
       b_fcus->nni->swap_node_v3 = b_fcus->rght;
       b_fcus->nni->swap_node_v4 = v4;
     }
   else
     {
       b_fcus->nni->score        = +1.0;
       b_fcus->nni->swap_node_v1 = NULL;
       b_fcus->nni->swap_node_v2 = NULL;
       b_fcus->nni->swap_node_v3 = NULL;
       b_fcus->nni->swap_node_v4 = NULL;
     }
}

/*********************************************************/

void Swap(t_node *a, t_node *b, t_node *c, t_node *d, t_tree *tree)
{
  int ab, ba, cd, dc, bc;
  int i;


  /* \             /d      \             /a
   *  \           /         \           /
   *   \b__...__c/    ->     \b__...__c/
   *   /         \	     /		\
   *  /           \	    /		 \
   * /a            \  	   /d             \ 
   *
   * nodes b and c are not necessarily on the same branch 
   */


#ifdef DEBUG
  if(!a || !b || !c || !d)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  ab = ba = cd = dc = bc = -1;

  For(i,3) if(a->v[i] == b) { ab = i; break; }
  For(i,3) if(b->v[i] == a) { ba = i; break; }
  For(i,3) if(c->v[i] == d) { cd = i; break; }
  For(i,3) if(d->v[i] == c) { dc = i; break; }
  For(i,3) if(b->v[i] == c) { bc = i; break; }

#ifdef DEBUG
  if(ab < 0 || ba < 0 || cd < 0 || dc < 0)
    {
      PhyML_Printf("\n. Nodes %d %d %d %d\n",a->num,b->num,c->num,d->num);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  a->v[ab] = c;
  d->v[dc] = b;
  b->v[ba] = d;
  c->v[cd] = a;
  b->b[ba] = d->b[dc];
  c->b[cd] = a->b[ab];

  (a->b[ab]->left == b)?
  (a->b[ab]->left = c):
  (a->b[ab]->rght = c);

  (d->b[dc]->left == c)?
  (d->b[dc]->left = b):
  (d->b[dc]->rght = b);

  For(i,3)
    {
      if(a->b[ab]->left->v[i] == a->b[ab]->rght) a->b[ab]->l_r = i;
      if(a->b[ab]->rght->v[i] == a->b[ab]->left) a->b[ab]->r_l = i;
      if(d->b[dc]->left->v[i] == d->b[dc]->rght) d->b[dc]->l_r = i;
      if(d->b[dc]->rght->v[i] == d->b[dc]->left) d->b[dc]->r_l = i;
    }

  a->b[ab]->l_v1 = a->b[ab]->l_v2 =
  a->b[ab]->r_v1 = a->b[ab]->r_v2 =
  d->b[dc]->l_v1 = d->b[dc]->l_v2 =
  d->b[dc]->r_v1 = d->b[dc]->r_v2 = -1;

  For(i,3)
    {
      if(i != a->b[ab]->l_r)
	{
	  if(a->b[ab]->l_v1 < 0) a->b[ab]->l_v1 = i;
	  else a->b[ab]->l_v2 = i;
	}
      if(i != a->b[ab]->r_l)
	{
	  if(a->b[ab]->r_v1 < 0) a->b[ab]->r_v1 = i;
	  else a->b[ab]->r_v2 = i;
	}
      if(i != d->b[dc]->l_r)
	{
	  if(d->b[dc]->l_v1 < 0) d->b[dc]->l_v1 = i;
	  else d->b[dc]->l_v2 = i;
	}
      if(i != d->b[dc]->r_l)
	{
	  if(d->b[dc]->r_v1 < 0) d->b[dc]->r_v1 = i;
	  else d->b[dc]->r_v2 = i;
	}
    }
  Update_Dirs(tree);
}

/*********************************************************/

void Update_All_Partial_Lk(t_edge *b_fcus, t_tree *tree)
{

  Update_SubTree_Partial_Lk(b_fcus->left->b[b_fcus->l_v1],
			    b_fcus->left,
			    b_fcus->left->v[b_fcus->l_v1],
			    tree);

  Update_SubTree_Partial_Lk(b_fcus->left->b[b_fcus->l_v2],
			    b_fcus->left,
			    b_fcus->left->v[b_fcus->l_v2],
			    tree);

  Update_SubTree_Partial_Lk(b_fcus->rght->b[b_fcus->r_v1],
			    b_fcus->rght,
			    b_fcus->rght->v[b_fcus->r_v1],
			    tree);

  Update_SubTree_Partial_Lk(b_fcus->rght->b[b_fcus->r_v2],
			    b_fcus->rght,
			    b_fcus->rght->v[b_fcus->r_v2],
			    tree);

  tree->c_lnL = Lk_At_Given_Edge(b_fcus,tree);
}

/*********************************************************/

void Update_SubTree_Partial_Lk(t_edge *b_fcus, t_node *a, t_node *d, t_tree *tree)
{
  int i;

  Update_P_Lk(tree,b_fcus,a);
  if(d->tax) return;
  else For(i,3) if(d->v[i] != a)
    Update_SubTree_Partial_Lk(d->b[i],d,d->v[i],tree);
}

/*********************************************************/

calign *Make_Cseq(int n_otu, int crunch_len, int state_len, int init_len, char **sp_names)
{
  calign *cdata;
  int j;

  cdata           = (calign *)mCalloc(1,sizeof(calign));
  cdata->n_otu    = n_otu;
  cdata->c_seq    = (align **)mCalloc(n_otu,sizeof(align *));
  cdata->b_frq    = (phydbl *)mCalloc(T_MAX_ALPHABET,sizeof(phydbl));
  cdata->wght     = (int *)mCalloc(crunch_len,sizeof(int));
  cdata->ambigu   = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->invar    = (short int *)mCalloc(crunch_len,sizeof(short int));
  cdata->sitepatt = (int *)mCalloc(init_len,sizeof(int ));
  cdata->format   = 0;

  cdata->crunch_len = crunch_len;
  cdata->init_len   = init_len;
  cdata->obs_pinvar = .0;

  For(j,n_otu)
    {
      cdata->c_seq[j]            = (align *)mCalloc(1,sizeof(align));
      cdata->c_seq[j]->name      = (char *)mCalloc((int)(strlen(sp_names[j])+1),sizeof(char));
      strcpy(cdata->c_seq[j]->name,sp_names[j]);
      cdata->c_seq[j]->state     = (char *)mCalloc(crunch_len*state_len,sizeof(char));
      cdata->c_seq[j]->is_ambigu = (short int *)mCalloc(crunch_len,sizeof(short int));
    }

  return cdata;
}

/*********************************************************/

t_treelist *Make_Treelist(int list_size)
{
  t_treelist *tlist;

  tlist = (t_treelist *)mCalloc(1,sizeof(t_treelist));
  tlist->list_size = list_size;
  tlist->tree = (t_tree **)mCalloc(list_size,sizeof(t_tree *));

  return tlist;
}


/*********************************************************/

void Copy_Seq_Names_To_Tip_Labels(t_tree *tree, calign *data)
{
  int i;
  For(i,tree->n_otu)
    {
      strcpy(tree->noeud[i]->name,data->c_seq[i]->name);
    }
}

/*********************************************************/

calign *Copy_Cseq(calign *ori, option *io)
{
  calign *new;
  int i,j,k,n_otu,c_len;
  char **sp_names;

  n_otu = ori->n_otu;
  c_len = ori->crunch_len;

  sp_names = (char **)mCalloc(n_otu,sizeof(char *));
  For(i,n_otu)
    {
      sp_names[i] = (char *)mCalloc(strlen(ori->c_seq[i]->name)+1,sizeof(char));
      strcpy(sp_names[i],ori->c_seq[i]->name);
    }

  new = Make_Cseq(n_otu,c_len+1,io->mod->state_len,ori->init_len,sp_names);

  new->obs_pinvar = ori->obs_pinvar;

  For(i,ori->init_len) new->sitepatt[i] = ori->sitepatt[i];

  For(j,ori->crunch_len)
    {
      For(i,ori->n_otu) 
	{
	  For(k,io->mod->state_len) 
	    {
	      new->c_seq[i]->state[j*io->mod->state_len+k] = 
		ori->c_seq[i]->state[j*io->mod->state_len+k];
	    }
	  new->c_seq[i]->is_ambigu[j] = ori->c_seq[i]->is_ambigu[j];
	}

      new->wght[j]   = ori->wght[j];
      new->ambigu[j] = ori->ambigu[j];
      new->invar[j]  = ori->invar[j];
    }

  For(i,ori->n_otu)
    {
      new->c_seq[i]->len = ori->c_seq[i]->len;
      strcpy(new->c_seq[i]->name,ori->c_seq[i]->name);
    }

  For(i,ori->n_otu) new->c_seq[i]->state[c_len*io->mod->state_len] = '\0';

  For(i,io->mod->ns) new->b_frq[i] = ori->b_frq[i];

  new->init_len           = ori->init_len;
  new->clean_len          = ori->clean_len;
  new->crunch_len         = ori->crunch_len;
  new->n_otu              = ori->n_otu;

  For(i,n_otu) Free(sp_names[i]);
  Free(sp_names);

  return new;
}

/*********************************************************/

optimiz *Make_Optimiz()
{
  optimiz *s_opt;
  s_opt = (optimiz *)mCalloc(1,sizeof(optimiz));
  return s_opt;
}

/*********************************************************/


int Filexists(char *filename)
{
  FILE *fp;
  fp =fopen(filename,"r");
  if (fp) {
    fclose(fp);
    return 1;
  } else
    return 0;
}

/*********************************************************/

FILE *Openfile(char *filename, int mode)
{
  /* mode = 0 -> read */
  /* mode = 1 -> write */
  /* mode = 2 -> append */

  FILE *fp;
  char *s;
  int open_test=0;

/*   s = (char *)mCalloc(T_MAX_FILE,sizeof(char)); */

/*   strcpy(s,filename); */

  s = filename;

  fp = NULL;

  switch(mode)
    {
    case 0 :
      {
	while(!(fp = (FILE *)fopen(s,"r")) && ++open_test<10)
	  {
	    PhyML_Printf("\n. Can't open file '%s', enter a new name : ",s);
	    Getstring_Stdin(s);
	  }
	break;
      }
    case 1 :
      {
	fp = (FILE *)fopen(s,"w");
	break;
      }
    case 2 :
      {
	fp = (FILE *)fopen(s,"a");
	break;
      }

    default : break;

    }

/*   Free(s); */

  return fp;
}

/*********************************************************/

void Print_Fp_Out(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *io, int n_data_set, int num_tree)
{
  char *s;
  div_t hour,min;
  int i;

/*   For(i,2*tree->n_otu-3) fprintf(fp_out,"\n. * Edge %3d: %f",i,tree->t_edges[i]->l); */


  if((!n_data_set) || (!num_tree))
    {
      Print_Banner_Small(fp_out);
    }

  PhyML_Fprintf(fp_out,"\n. Sequence filename: \t\t\t%s", Basename(io->in_align_file));
  PhyML_Fprintf(fp_out,"\n. Data set: \t\t\t\t#%d",n_data_set);

  if(io->mod->s_opt->random_input_tree)
    PhyML_Fprintf(fp_out,"\n. Random init tree: \t\t\t#%d",num_tree+1);
  else if(io->n_trees > 1)
    PhyML_Fprintf(fp_out,"\n. Starting tree number: \t\t\t#%d",num_tree+1);
  
  if(io->mod->s_opt->opt_topo)
    {
      if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Fprintf(fp_out,"\n. Tree topology search : \t\tNNIs");
      else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Fprintf(fp_out,"\n. Tree topology search : \t\tSPRs");
      else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Fprintf(fp_out,"\n. Tree topology search : \t\tBest of NNIs and SPRs");
    }
  else
    {
      PhyML_Fprintf(fp_out,"\n. Tree topology: \t\t\tfixed");
    }


  /* was after Sequence file ; moved here FLT */
  s = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  if(io->in_tree == 2)
    {
      strcat(strcat(strcat(s,"user tree ("),io->in_tree_file),")");
    }
  else
    {
      if(!io->mod->s_opt->random_input_tree)
	{
	  if(io->in_tree == 0)
	    strcat(s,"BioNJ");
	  if(io->in_tree == 1)
	    strcat(s,"parsimony");
	}
      else
	{
	  strcat(s,"random tree");
	}
    }

  PhyML_Fprintf(fp_out,"\n. Initial tree: \t\t\t%s",s);
  Free(s);

  if(tree->io->datatype == NT)
    {
      fprintf(fp_out,"\n. Model of nucleotides substitution: \t%s",io->mod->modelname);
      if(io->mod->whichmodel == CUSTOM)
      fprintf(fp_out," (%s)",io->mod->custom_mod_string);
    }
  else if(tree->io->datatype == AA)
    {
      fprintf(fp_out,"\n. Model of amino acids substitution: \t%s",io->mod->modelname);
    }
  else
    {
      fprintf(fp_out,"\n. Substitution model: \t\t\t%s",io->mod->modelname);
    }


  PhyML_Fprintf(fp_out,"\n. Number of taxa: \t\t\t%d",tree->n_otu);/*added FLT*/

  PhyML_Fprintf(fp_out,"\n. Log-likelihood: \t\t\t%.5f",tree->c_lnL);/*was last ; moved here FLT*/

  Unconstraint_Lk(tree);
  PhyML_Fprintf(fp_out,"\n. Unconstrained likelihood: \t\t%.5f",tree->unconstraint_lk);

  PhyML_Fprintf(fp_out,"\n. Parsimony: \t\t\t\t%d",tree->c_pars);

  PhyML_Fprintf(fp_out,"\n. Tree size: \t\t\t\t%.5f",tree->size);

  if(tree->mod->n_catg > 1 && tree->mod->free_mixt_rates == NO)
    {
      PhyML_Fprintf(fp_out,"\n. Discrete gamma model: \t\t%s","Yes");
      PhyML_Fprintf(fp_out,"\n  - Number of categories: \t\t%d",tree->mod->n_catg);
      PhyML_Fprintf(fp_out,"\n  - Gamma shape parameter: \t\t%.3f",tree->mod->alpha);
    }
  else if(tree->mod->free_mixt_rates == YES)
    {
      PhyML_Fprintf(fp_out,"\n. Discrete gamma model: \t\t%s","No");
      PhyML_Fprintf(fp_out,"\n  - Number of categories: \t\t%d",tree->mod->n_catg);
      For(i,tree->mod->n_catg)
	{
	  PhyML_Fprintf(fp_out,"\n  - Relative rate in class %d: \t\t%.5f [prop=%4f] \t\t",i+1,tree->mod->gamma_rr[i],tree->mod->gamma_r_proba[i]);
	}
    }

  if(tree->mod->invar) PhyML_Fprintf(fp_out,"\n. Proportion of invariant: \t\t%.3f",tree->mod->pinvar);

  /*was before Discrete gamma model ; moved here FLT*/
  if((tree->mod->whichmodel == K80)   ||
     (tree->mod->whichmodel == HKY85) ||
     (tree->mod->whichmodel == F84))
    PhyML_Fprintf(fp_out,"\n. Transition/transversion ratio: \t%.3f",tree->mod->kappa);
  else if(tree->mod->whichmodel == TN93)
    {
      PhyML_Fprintf(fp_out,"\n. Transition/transversion ratio for purines: \t\t\t%.3f",
	      tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
      PhyML_Fprintf(fp_out,"\n. Transition/transversion ratio for pyrimidines: \t\t\t%.3f",
	      tree->mod->kappa*2./(1.+tree->mod->lambda));
    }

  if(tree->io->datatype == NT)
    {
      PhyML_Fprintf(fp_out,"\n. Nucleotides frequencies:");
      PhyML_Fprintf(fp_out,"\n  - f(A)=%8.5f",tree->mod->pi[0]);
      PhyML_Fprintf(fp_out,"\n  - f(C)=%8.5f",tree->mod->pi[1]);
      PhyML_Fprintf(fp_out,"\n  - f(G)=%8.5f",tree->mod->pi[2]);
      PhyML_Fprintf(fp_out,"\n  - f(T)=%8.5f",tree->mod->pi[3]);
    }

  /*****************************************/
  if((tree->mod->whichmodel == GTR) ||
     (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      Update_Qmat_GTR(tree->mod->rr,
		      tree->mod->rr_val,
		      tree->mod->rr_num,
		      tree->mod->pi,
		      tree->mod->qmat);

      PhyML_Fprintf(fp_out,"\n");
      PhyML_Fprintf(fp_out,". GTR relative rate parameters : \n");
      PhyML_Fprintf(fp_out,"  A <-> C   %8.5f\n",  tree->mod->rr[0]);
      PhyML_Fprintf(fp_out,"  A <-> G   %8.5f\n",  tree->mod->rr[1]);
      PhyML_Fprintf(fp_out,"  A <-> T   %8.5f\n",  tree->mod->rr[2]);
      PhyML_Fprintf(fp_out,"  C <-> G   %8.5f\n",  tree->mod->rr[3]);
      PhyML_Fprintf(fp_out,"  C <-> T   %8.5f\n",  tree->mod->rr[4]);
      PhyML_Fprintf(fp_out,"  G <-> T   %8.5f\n",tree->mod->rr[5]);


      PhyML_Fprintf(fp_out,"\n. Instantaneous rate matrix : ");
      PhyML_Fprintf(fp_out,"\n  [A---------C---------G---------T------]\n");
      For(i,4)
	{
	  PhyML_Fprintf(fp_out,"  ");
	  For(j,4)
	    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
	  PhyML_Fprintf(fp_out,"\n");
	}
      PhyML_Fprintf(fp_out,"\n");
    }
  /*****************************************/


  if(io->ratio_test == 1)
    {
      PhyML_Fprintf(fp_out,". aLRT statistics to test branches");
    }
  else if(io->ratio_test == 2)
    {
      PhyML_Fprintf(fp_out,". aLRT branch supports (cubic approximation, mixture of Chi2s distribution)");
    }


  PhyML_Fprintf(fp_out,"\n");
  PhyML_Fprintf(fp_out,"\n. Run ID:\t\t\t\t%s", (io->appebr_run_ID) ? (io->run_id_string): ("none"));
  PhyML_Fprintf(fp_out,"\n. Random seed:\t\t\t\t%d", io->r_seed);
  PhyML_Fprintf(fp_out,"\n. Subtree patterns aliasing:\t\t%s",io->do_alias_subpatt?"yes":"no");
  PhyML_Fprintf(fp_out,"\n. Version:\t\t\t\t%s", VERSION);

  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  PhyML_Fprintf(fp_out,"\n. Time used:\t\t\t\t%dh%dm%ds (%d seconds)", hour.quot,min.quot,(int)(t_end-t_beg)%60,(int)(t_end-t_beg));

  PhyML_Fprintf(fp_out,"\n\n");
  PhyML_Fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp_out," Suggested citations:\n");
  PhyML_Fprintf(fp_out," S. Guindon, JF. Dufayard, V. Lefort, M. Anisimova, W. Hordijk, O. Gascuel\n");
  PhyML_Fprintf(fp_out," \"New algorithms and methods to estimate maximum-likelihood phylogenies: assessing the performance of PhyML 3.0.\"\n");
  PhyML_Fprintf(fp_out," Systematic Biology. 2010. 59(3):307-321.\n");
  PhyML_Fprintf(fp_out,"\n");
  PhyML_Fprintf(fp_out," S. Guindon & O. Gascuel\n");
  PhyML_Fprintf(fp_out," \"A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood\"\n");
  PhyML_Fprintf(fp_out," Systematic Biology. 2003. 52(5):696-704.\n");
  PhyML_Fprintf(fp_out," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/
/*FLT wrote this function*/
void Print_Fp_Out_Lines(FILE *fp_out, time_t t_beg, time_t t_end, t_tree *tree, option *io, int n_data_set)
{
  char *s;
  /*div_t hour,min;*/

  if (n_data_set==1)
      {

	PhyML_Fprintf(fp_out,". Sequence file : [%s]\n\n", Basename(io->in_align_file));

	if((tree->io->datatype == NT) || (tree->io->datatype == AA))
	  {
	    (tree->io->datatype == NT)?
	      (PhyML_Fprintf(fp_out,". Model of nucleotides substitution : %s\n\n",io->mod->modelname)):
	      (PhyML_Fprintf(fp_out,". Model of amino acids substitution : %s\n\n",io->mod->modelname));
	  }

	s = (char *)mCalloc(100,sizeof(char));

	switch(io->in_tree)
	  {
	  case 0: { strcpy(s,"BioNJ");     break; }
	  case 1: { strcpy(s,"parsimony"); break; }
	  case 2: { strcpy(s,"user tree ("); 
	            strcat(s,io->in_tree_file); 
	            strcat(s,")");         break; }
	  }

	PhyML_Fprintf(fp_out,". Initial tree : [%s]\n\n",s);

	Free(s);

	PhyML_Fprintf(fp_out,"\n");

	/*headline 1*/
	PhyML_Fprintf(fp_out, ". Data\t");

	PhyML_Fprintf(fp_out,"Nb of \t");

	PhyML_Fprintf(fp_out,"Likelihood\t");

	PhyML_Fprintf(fp_out, "Discrete   \t");

	if(tree->mod->n_catg > 1)
	  PhyML_Fprintf(fp_out, "Number of \tGamma shape\t");

	PhyML_Fprintf(fp_out,"Proportion of\t");

	if(tree->mod->whichmodel <= 6)
	  PhyML_Fprintf(fp_out,"Transition/ \t");

	PhyML_Fprintf(fp_out,"Nucleotides frequencies               \t");

	if((tree->mod->whichmodel == GTR) ||
	   (tree->mod->whichmodel == CUSTOM))
	  PhyML_Fprintf(fp_out,"Instantaneous rate matrix              \t");

	/*    PhyML_Fprintf(fp_out,"Time\t");*/

	PhyML_Fprintf(fp_out, "\n");


	/*headline 2*/
	PhyML_Fprintf(fp_out, "  set\t");

	PhyML_Fprintf(fp_out,"taxa\t");

	PhyML_Fprintf(fp_out,"loglk     \t");

	PhyML_Fprintf(fp_out, "gamma model\t");

	if(tree->mod->n_catg > 1)
	  PhyML_Fprintf(fp_out, "categories\tparameter  \t");

	PhyML_Fprintf(fp_out,"invariant    \t");

	if(tree->mod->whichmodel <= 6)
	  PhyML_Fprintf(fp_out,"transversion\t");

	PhyML_Fprintf(fp_out,"f(A)      f(C)      f(G)      f(T)    \t");

	if((tree->mod->whichmodel == GTR) ||
	   (tree->mod->whichmodel == CUSTOM))
	  PhyML_Fprintf(fp_out,"[A---------C---------G---------T------]\t");

	/*    PhyML_PhyML_Fprintf(fp_out,"used\t");*/

	PhyML_Fprintf(fp_out, "\n");


	/*headline 3*/
	if(tree->mod->whichmodel == TN93)
	  {
	    PhyML_Fprintf(fp_out,"    \t      \t          \t           \t");
	    if(tree->mod->n_catg > 1) PhyML_Fprintf(fp_out,"         \t         \t");
	    PhyML_Fprintf(fp_out,"             \t");
	    PhyML_Fprintf(fp_out,"purines pyrimid.\t");
	    PhyML_Fprintf(fp_out, "\n");
          }

          PhyML_Fprintf(fp_out, "\n");
      }


  /*line items*/

  PhyML_Fprintf(fp_out,"  #%d\t",n_data_set);

  PhyML_Fprintf(fp_out,"%d   \t",tree->n_otu);

  PhyML_Fprintf(fp_out,"%.5f\t",tree->c_lnL);

  PhyML_Fprintf(fp_out,"%s        \t",
	  (tree->mod->n_catg>1)?("Yes"):("No "));
  if(tree->mod->n_catg > 1)
    {
      PhyML_Fprintf(fp_out,"%d        \t",tree->mod->n_catg);
      PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->alpha);
    }

  /*if(tree->mod->invar)*/
    PhyML_Fprintf(fp_out,"%.3f    \t",tree->mod->pinvar);

  if(tree->mod->whichmodel <= 5)
    {
      PhyML_Fprintf(fp_out,"%.3f     \t",tree->mod->kappa);
    }
  else if(tree->mod->whichmodel == TN93)
    {
      PhyML_Fprintf(fp_out,"%.3f   ",
	      tree->mod->kappa*2.*tree->mod->lambda/(1.+tree->mod->lambda));
      PhyML_Fprintf(fp_out,"%.3f\t",
	      tree->mod->kappa*2./(1.+tree->mod->lambda));
    }


  if(tree->io->datatype == NT)
    {
      PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[0]);
      PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[1]);
      PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->pi[2]);
      PhyML_Fprintf(fp_out,"%8.5f\t",tree->mod->pi[3]);
    }
  /*
  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );

  min.quot -= hour.quot*60;

  PhyML_Fprintf(fp_out,"%dh%dm%ds\t", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  if(t_end-t_beg > 60)
    PhyML_Fprintf(fp_out,". -> %d seconds\t",(int)(t_end-t_beg));
  */

  /*****************************************/
  if((tree->mod->whichmodel == GTR) || (tree->mod->whichmodel == CUSTOM))
    {
      int i,j;

      For(i,4)
	{
	  if (i!=0) {
	    /*format*/
	    PhyML_Fprintf(fp_out,"      \t     \t          \t           \t");
	    if(tree->mod->n_catg > 1) PhyML_Fprintf(fp_out,"          \t           \t");
	    PhyML_Fprintf(fp_out,"             \t                                      \t");
	  }
	  For(j,4)
	    PhyML_Fprintf(fp_out,"%8.5f  ",tree->mod->qmat[i*4+j]);
	  if (i<3) PhyML_Fprintf(fp_out,"\n");
	}
    }
  /*****************************************/

  PhyML_Fprintf(fp_out, "\n\n");
}

/*********************************************************/

matrix *K80_dist(calign *data, phydbl g_shape)
{
  int i,j,k;
  int diff;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;

  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

  For(i,data->c_seq[0]->len)
    {
      For(j,data->n_otu-1)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if(((data->c_seq[j]->state[i] == 'A' || data->c_seq[j]->state[i] == 'G') &&
		  (data->c_seq[k]->state[i] == 'C' || data->c_seq[k]->state[i] == 'T'))||
		 ((data->c_seq[j]->state[i] == 'C' || data->c_seq[j]->state[i] == 'T') &&
		  (data->c_seq[k]->state[i] == 'A' || data->c_seq[k]->state[i] == 'G')))
		{
		  diff++;
		  mat->Q[j][k]+=data->wght[i];
		  len[j][k]+=data->wght[i];
		  len[k][j]=len[j][k];
		}
	      
	      else
		if(((data->c_seq[j]->state[i] == 'A' && data->c_seq[k]->state[i] == 'G') ||
		    (data->c_seq[j]->state[i] == 'G' && data->c_seq[k]->state[i] == 'A'))||
		   ((data->c_seq[j]->state[i] == 'C' && data->c_seq[k]->state[i] == 'T') ||
		    (data->c_seq[j]->state[i] == 'T' && data->c_seq[k]->state[i] == 'C')))
		  {
		    diff++;
		    mat->P[j][k]+=data->wght[i];
		    len[j][k]+=data->wght[i];
		    len[k][j]=len[j][k];
		  }
		else
		  if((data->c_seq[j]->state[i] == 'A' ||
		      data->c_seq[j]->state[i] == 'C' ||
		      data->c_seq[j]->state[i] == 'G' ||
		      data->c_seq[j]->state[i] == 'T')&&
		     (data->c_seq[k]->state[i] == 'A' ||
		      data->c_seq[k]->state[i] == 'C' ||
		      data->c_seq[k]->state[i] == 'G' ||
		      data->c_seq[k]->state[i] == 'T'))
		    {
		      len[j][k]+=data->wght[i];
		      len[k][j]=len[j][k];
		    }
	    }
	}
    }


  For(i,data->n_otu-1)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j] > .0)
	  {
	    mat->P[i][j] /= len[i][j];
	    mat->Q[i][j] /= len[i][j];
	  }
	else
	  {
	    mat->P[i][j] = .5;
	    mat->Q[i][j] = .5;
	  }

	mat->P[j][i] = mat->P[i][j];
	mat->Q[j][i] = mat->Q[i][j];


	if((1-2*mat->P[i][j]-mat->Q[i][j] <= .0) || (1-2*mat->Q[i][j] <= .0))
	  {
	    mat->dist[i][j] = -1.;
	    mat->dist[j][i] = -1.;
	    continue;
	  }

	mat->dist[i][j] = (g_shape/2)*
	  (POW(1-2*mat->P[i][j]-mat->Q[i][j],-1./g_shape) +
	   0.5*POW(1-2*mat->Q[i][j],-1./g_shape) - 1.5);

	if(mat->dist[i][j] > DIST_MAX) mat->dist[i][j] = DIST_MAX;

	mat->dist[j][i] = mat->dist[i][j];
      }

  For(i,data->n_otu) free(len[i]);
  free(len);
  return mat;
}

/*********************************************************/

matrix *JC69_Dist(calign *data, model *mod)
{
  int site,i,j,k;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;
  int datatype;


  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);

  datatype = mod->io->datatype;

  For(site,data->c_seq[0]->len)
    {
      For(j,data->n_otu-1)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if((!Is_Ambigu(data->c_seq[j]->state+site*mod->io->mod->state_len,datatype,mod->io->mod->state_len)) &&
		 (!Is_Ambigu(data->c_seq[k]->state+site*mod->io->mod->state_len,datatype,mod->io->mod->state_len)))
		{
		  len[j][k]+=data->wght[site];
		  len[k][j]=len[j][k];


		  if(strncmp(data->c_seq[j]->state+site*mod->io->mod->state_len,
			     data->c_seq[k]->state+site*mod->io->mod->state_len,mod->io->mod->state_len))
/* 		  if(!Are_Compatible(data->c_seq[j]->state+site*mod->io->mod->state_len, */
/* 				     data->c_seq[k]->state+site*mod->io->mod->state_len, */
/* 				     mod->io->mod->state_len, */
/* 				     mod->io->datatype)) */
		    mat->P[j][k]+=data->wght[site];
		}
	    }
	}
    }
  

  For(i,data->n_otu-1)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j] > .0) mat->P[i][j] /= len[i][j];
	else               mat->P[i][j] = 1.;

	mat->P[j][i] = mat->P[i][j];

	if((1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]) < .0) mat->dist[i][j] = -1.;
	else
	  mat->dist[i][j] = -(mod->ns-1.)/(mod->ns)*(phydbl)LOG(1.-(mod->ns)/(mod->ns-1.)*mat->P[i][j]);

/* 	PhyML_Printf("\n. Incorrect JC distances"); */
/* 	mat->dist[i][j] = len[i][j]; */

	if(mat->dist[i][j] > DIST_MAX) mat->dist[i][j] = DIST_MAX;

	mat->dist[j][i] = mat->dist[i][j];
      }

  For(i,data->n_otu) free(len[i]);
  free(len);

  return mat;
}

/*********************************************************/

matrix *Hamming_Dist(calign *data, model *mod)
{
  int i,j,k;
  phydbl unc_len;
  matrix *mat;
  phydbl **len;
  int datatype;

  len = (phydbl **)mCalloc(data->n_otu,sizeof(phydbl *));
  For(i,data->n_otu)
    len[i] = (phydbl *)mCalloc(data->n_otu,sizeof(phydbl));

  unc_len = .0;

  mat = Make_Mat(data->n_otu);
  Init_Mat(mat,data);
  
  datatype = mod->io->datatype;

  For(i,data->crunch_len)
    {
      For(j,data->n_otu-1)
	{
	  for(k=j+1;k<data->n_otu;k++)
	    {
	      if((!Is_Ambigu(data->c_seq[j]->state+i*mod->io->mod->state_len,datatype,mod->io->mod->state_len)) &&
		 (!Is_Ambigu(data->c_seq[k]->state+i*mod->io->mod->state_len,datatype,mod->io->mod->state_len)))
		{
		  len[j][k]+=data->wght[i];
		  len[k][j]=len[j][k];
/* 		  if(data->c_seq[j]->state[i] != data->c_seq[k]->state[i]) */
		  if(!Are_Compatible(data->c_seq[j]->state+i*mod->io->mod->state_len,
				     data->c_seq[k]->state+i*mod->io->mod->state_len,
				     mod->io->mod->state_len,
				     mod->io->datatype))
		    {
		      mat->P[j][k]+=data->wght[i];
		    }
		}
	    }
	}
    }

  For(i,data->n_otu-1)
    for(j=i+1;j<data->n_otu;j++)
      {
	if(len[i][j] > .0)
	  {
	    mat->P[i][j] /= len[i][j];
	  }
	else
	  {
	    mat->P[i][j] = 1.;
	  }

	mat->P[j][i] = mat->P[i][j];

	mat->dist[i][j] = mat->P[i][j];


	if(mat->dist[i][j] > DIST_MAX)
	  {
	    mat->dist[i][j] = DIST_MAX;
	  }
	mat->dist[j][i] = mat->dist[i][j];
      }

  For(i,data->n_otu) free(len[i]);
  free(len);

  return mat;
}

/*********************************************************/
/* Test if the given site pattern is invariant. Does not handle ambiguities */

int Is_Invar(int patt_num, int stepsize, int datatype, calign *data)
{
  int i, j;

  For(i,data->n_otu)
    {
      For(j,data->n_otu)
	{
	  if(!(Are_Compatible(data->c_seq[i]->state+patt_num,
			      data->c_seq[j]->state+patt_num,
			      stepsize,
			      datatype))) 
	    {
	      break;
	    }
	}
      if(j != data->n_otu) break;
    }
  
  if(i == data->n_otu) return 1;
  else                 return 0;
}


/*********************************************************/

int Is_Ambigu(char *state, int datatype, int stepsize)
{
  int val,i;

  val = -1;
  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  switch(state[i])
	    {
	    case 'A' : case 'C' : case 'G' : case 'T' : case 'U' : { val=0; break; }
	    default : { val=1; break; }
	    }
	  if(val == 1) break;
	}
    }
  else if(datatype == AA)
    {
      switch(state[0])
	{
	case 'X' : case '?' : case '-' : case '.' : {val=1; break; }
	default : { val=0; break; }
	}
    }
  else if(datatype == GENERIC)
    {
      int i;
      For(i,stepsize) if(!isdigit(state[i])) break;
      if(i == stepsize) val = 0;
      else              val = 1;
    }

  return val;
}

/*********************************************************/

void Check_Ambiguities(calign *data, int datatype, int stepsize)
{
  int i,j;

  For(j,data->crunch_len) 
    {
      For(i,data->n_otu)
	{
	  data->ambigu[j]              = 0;
	  data->c_seq[i]->is_ambigu[j] = 0;
	}

      For(i,data->n_otu)
	{
	  if(Is_Ambigu(data->c_seq[i]->state+j*stepsize,
		       datatype,
		       stepsize))
	    {
	      data->ambigu[j]              = 1;
	      data->c_seq[i]->is_ambigu[j] = 1;
	    }
	}
    }
}

/*********************************************************/

int Get_State_From_Ui(int ui, int datatype)
{
  if(datatype == NT)
    {
      switch(ui)
	{
	case 1 : {return 0; break;}
	case 2 : {return 1; break;}
	case 4 : {return 2; break;}
	case 8 : {return 3; break;}
	default : 
	  {
	    PhyML_Printf("\n. ui=%d",ui);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else if(datatype == AA)
    {
      switch(ui)
	{
	case 1 :      {return 0;  break;}
	case 2 :      {return 1;  break;} 
	case 4 :      {return 2;  break;} 
	case 8 :      {return 3;  break;} 
	case 16 :     {return 4;  break;} 
	case 32 :     {return 5;  break;} 
	case 64 :     {return 6;  break;} 
	case 128 :    {return 7;  break;} 
	case 256 :    {return 8;  break;} 
	case 512 :    {return 9;  break;} 
	case 1024 :   {return 10; break;} 
	case 2048 :   {return 11; break;} 
	case 4096 :   {return 12; break;} 
	case 8192 :   {return 13; break;} 
	case 16384 :  {return 14; break;} 
	case 32768 :  {return 15; break;} 
	case 65536 :  {return 16; break;} 
	case 131072 : {return 17; break;} 
	case 262144 : {return 18; break;} 
	case 524288 : {return 19; break;} 
	default : 
	  {
	    PhyML_Printf("\n. ui=%d",ui);
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }
	}
    }
  else
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  return -1;
}

/*********************************************************/

int Assign_State(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;

  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  switch(c[i])
	    {
	    case 'A' : {state[i]=0;  break;}
	    case 'C' : {state[i]=1;  break;}
	    case 'G' : {state[i]=2;  break;}
	    case 'T' : {state[i]=3;  break;}
	    case 'U' : {state[i]=3;  break;}
	    default  : {state[i]=-1; break;}
	    }
	}
      return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
    }
  else if(datatype == AA)
    {
      switch(c[0])
	{
	case 'A' : {state[0]=0 ; break;}
	case 'R' : {state[0]=1 ; break;}
	case 'N' : {state[0]=2 ; break;}
	case 'D' : {state[0]=3 ; break;}
	case 'C' : {state[0]=4 ; break;}
	case 'Q' : {state[0]=5 ; break;}
	case 'E' : {state[0]=6 ; break;}
	case 'G' : {state[0]=7 ; break;}
	case 'H' : {state[0]=8 ; break;}
	case 'I' : {state[0]=9 ; break;}
	case 'L' : {state[0]=10; break;}
	case 'K' : {state[0]=11; break;}
	case 'M' : {state[0]=12; break;}
	case 'F' : {state[0]=13; break;}
	case 'P' : {state[0]=14; break;}
	case 'S' : {state[0]=15; break;}
	case 'T' : {state[0]=16; break;}
	case 'W' : {state[0]=17; break;}
	case 'Y' : {state[0]=18; break;}
	case 'V' : {state[0]=19; break;}

	case 'B' : {state[0] = 2; break;}
	case 'Z' : {state[0] = 5; break;}
	default  : {state[0]=-1;  break;}
	}
      return state[0];
    }
  else if(datatype == GENERIC)
    {
      char format[6];
      int ret;

      sprintf(format,"%%%dd",stepsize);
      ret = sscanf(c,format,state);
      if(!ret) state[0] = -1;      
      return state[0];
    }
  else
    {
      PhyML_Printf("\n. Not implemented yet.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  return -1;
}

/*********************************************************/

char Reciproc_Assign_State(int i_state, int datatype)
{
  
  if(datatype == NT)
    {
      i_state = i_state%4;
      switch(i_state)
	{
	case 0 :   {return 'A';  break;}
	case 1 :   {return 'C';  break;}
	case 2 :   {return 'G';  break;}
	case 3 :   {return 'T';  break;}
	default  : 
	  {
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else if(datatype == AA)
    {
      i_state = i_state%20;
      switch(i_state)
	{
	case 0  : {return 'A' ; break;}
	case 1  : {return 'R' ; break;}
	case 2  : {return 'N' ; break;}
	case 3  : {return 'D' ; break;}
	case 4  : {return 'C' ; break;}
	case 5  : {return 'Q' ; break;}
	case 6  : {return 'E' ; break;}
	case 7  : {return 'G' ; break;}
	case 8  : {return 'H' ; break;}
	case 9  : {return 'I' ; break;}
	case 10 : {return 'L';  break;}
	case 11 : {return 'K';  break;}
	case 12 : {return 'M';  break;}
	case 13 : {return 'F';  break;}
	case 14 : {return 'P';  break;}
	case 15 : {return 'S';  break;}
	case 16 : {return 'T';  break;}
	case 17 : {return 'W';  break;}
	case 18 : {return 'Y';  break;}
	case 19 : {return 'V';  break;}
	default  : 
	  {
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else if(datatype == GENERIC)
    {
      return (char)i_state;
    }
  return -1;
}

/*********************************************************/

int Assign_State_With_Ambiguity(char *c, int datatype, int stepsize)
{
  int state[3];
  int i;

  state[0] = state[1] = state[2] = -1;
  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  switch(c[i])
	    {
	    case 'A' : {state[i]= 0;  break;}
	    case 'C' : {state[i]= 1;  break;}
	    case 'G' : {state[i]= 2;  break;}
	    case 'T' : {state[i]= 3;  break;}
	    case 'U' : {state[i]= 3;  break;}
	    case 'M' : {state[i]= 4;  break;}
	    case 'R' : {state[i]= 5;  break;}
	    case 'W' : {state[i]= 6;  break;}
	    case 'S' : {state[i]= 7;  break;}
	    case 'Y' : {state[i]= 8;  break;}
	    case 'K' : {state[i]= 9;  break;}
	    case 'B' : {state[i]=10;  break;}
	    case 'D' : {state[i]=11;  break;}
	    case 'H' : {state[i]=12;  break;}
	    case 'V' : {state[i]=13;  break;}
	    case 'N' : case 'X' : case '?' : case 'O' : case '-' : {state[i]=T_MAX_ALPHABET-1;  break;}
	    default :
	      {
		PhyML_Printf("\n. Unknown character state : '%c'\n",c[i]);
		Warn_And_Exit("\n. Init failed (check the data type)\n");
		break;
	      }
	    }
	  return (stepsize>1)?(state[0]*16+state[1]*4+state[2]):(state[0]);
	}
    }
  else if(datatype == AA)
    {
      switch(c[0])
	{
	case 'A' : {state[0]= 0; break;}
	case 'R' : {state[0]= 1; break;}
	case 'N' : {state[0]= 2; break;}
	case 'D' : {state[0]= 3; break;}
	case 'C' : {state[0]= 4; break;}
	case 'Q' : {state[0]= 5; break;}
	case 'E' : {state[0]= 6; break;}
	case 'G' : {state[0]= 7; break;}
	case 'H' : {state[0]= 8; break;}
	case 'I' : {state[0]= 9; break;}
	case 'L' : {state[0]=10; break;}
	case 'K' : {state[0]=11; break;}
	case 'M' : {state[0]=12; break;}
	case 'F' : {state[0]=13; break;}
	case 'P' : {state[0]=14; break;}
	case 'S' : {state[0]=15; break;}
	case 'T' : {state[0]=16; break;}
	case 'W' : {state[0]=17; break;}
	case 'Y' : {state[0]=18; break;}
	case 'V' : {state[0]=19; break;}
	case 'B' : {state[0]= 2; break;}
	case 'Z' : {state[0]= 5; break;}
	case 'X' : case '?' : case '-' : {state[0]=T_MAX_ALPHABET-1; break;}
	default  : 
	  {
	    PhyML_Printf("\n. Unknown character state : %c\n",state[0]);
	    Warn_And_Exit("\n. Init failed (check the data type)\n");
	    break;
	  }
	}
      return state[0];
    }
  else if(datatype == GENERIC)
    {
      if(Is_Ambigu(c,GENERIC,stepsize)) state[0] = T_MAX_ALPHABET-1;
      else
	{
	  char format[6];
	  sprintf(format,"%%%dd",stepsize);
	  if(!sscanf(c,format,state))
	    {
	      PhyML_Printf("\n. Error reading character. Was expecting an integer, got '%c' instead.\n",c[0]);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
      return state[0];
    }

  return -1;
}

/*********************************************************/

void Clean_Tree_Connections(t_tree *tree)
{

  int i;
  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->v[0] = NULL;
      tree->noeud[i]->v[1] = NULL;
      tree->noeud[i]->v[2] = NULL;
      tree->noeud[i]->b[0] = NULL;
      tree->noeud[i]->b[1] = NULL;
      tree->noeud[i]->b[2] = NULL;
    }
}

/*********************************************************/

void Bootstrap(t_tree *tree)
{
  int *site_num, n_site;
  int replicate,j,k;
  int position,init_len;
  calign *boot_data;
  t_tree *boot_tree;
  model *boot_mod;
  matrix *boot_mat;
  char *s;
/*   phydbl rf; */

  tree->print_boot_val = 1;
  tree->print_alrt_val = 0;
  boot_tree            = NULL;

  site_num = (int *)mCalloc(tree->data->init_len,sizeof(int));

  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);

  n_site = 0;
  For(j,tree->data->crunch_len) For(k,tree->data->wght[j])
    {
      site_num[n_site] = j;
      n_site++;
    }

  boot_data = Copy_Cseq(tree->data,tree->io);


  PhyML_Printf("\n\n. Non parametric bootstrap analysis \n\n");
  PhyML_Printf("  ["); 
  

  For(replicate,tree->mod->bootstrap)
    {
      For(j,boot_data->crunch_len) boot_data->wght[j] = 0;

      init_len = 0;
      For(j,boot_data->init_len)
	{
	  position = Rand_Int(0,(int)(tree->data->init_len-1.0));
	  boot_data->wght[site_num[position]] += 1;
	  init_len++;
	}
      
      if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");

      init_len = 0;
      For(j,boot_data->crunch_len) init_len += boot_data->wght[j];

      if(init_len != tree->data->init_len) Warn_And_Exit("\n. Pb when copying sequences\n");

      if(tree->io->datatype == NT)      Get_Base_Freqs(boot_data);
      else if(tree->io->datatype == AA) Get_AA_Freqs(boot_data);

      if(tree->io->random_boot_seq_order) Randomize_Sequence_Order(boot_data);

      boot_mod        = Copy_Model(tree->mod);
      boot_mod->s_opt = tree->mod->s_opt; /* WARNING: re-using the same address here instead of creating a copying
					     requires to leave the value of s_opt unchanged during the boostrap. */
      boot_mod->io    = tree->io; /* WARNING: re-using the same address here instead of creating a copying
				     requires to leave the value of io unchanged during the boostrap. */
      Init_Model(boot_data,boot_mod,tree->io);

      if(tree->io->in_tree == 2)
	{
	  rewind(tree->io->fp_in_tree);
	  boot_tree = Read_Tree_File_Phylip(tree->io->fp_in_tree);
	}
      else
	{
	  boot_mat = ML_Dist(boot_data,boot_mod);
	  boot_mat->tree = Make_Tree_From_Scratch(boot_data->n_otu,boot_data);
	  Fill_Missing_Dist(boot_mat);
	  Bionj(boot_mat);
	  boot_tree = boot_mat->tree;
	  boot_tree->mat = boot_mat;
	}

      boot_tree->mod                = boot_mod;
      boot_tree->io                 = tree->io;
      boot_tree->data               = boot_data;
      boot_tree->both_sides         = 1;
      boot_tree->mod->s_opt->print  = 0;
      boot_tree->n_pattern          = boot_tree->data->crunch_len;
      boot_tree->io->print_site_lnl = 0;
      boot_tree->io->print_trace    = 0;

      if((boot_tree->mod->s_opt->random_input_tree) && (boot_tree->mod->s_opt->topo_search == SPR_MOVE)) Random_Tree(boot_tree);
      Order_Tree_CSeq(boot_tree,boot_data);
      Share_Lk_Struct(tree,boot_tree);
      Share_Spr_Struct(tree,boot_tree);
      Share_Pars_Struct(tree,boot_tree);
      Fill_Dir_Table(boot_tree);
      Update_Dirs(boot_tree);

      if(tree->mod->s_opt->greedy) Init_P_Lk_Tips_Double(boot_tree);
      else                         Init_P_Lk_Tips_Int(boot_tree);
      Init_Ui_Tips(boot_tree);
      Init_P_Pars_Tips(boot_tree);
      Br_Len_Not_Involving_Invar(boot_tree);

      
      if(boot_tree->io->do_alias_subpatt)
	{
	  boot_tree->update_alias_subpatt = YES;
	  Lk(boot_tree);
	  boot_tree->update_alias_subpatt = NO;
	}

      if(boot_tree->mod->s_opt->opt_topo)
	{
	  if(boot_tree->mod->s_opt->topo_search == NNI_MOVE) 
	    {
	      Simu_Loop(boot_tree);
	    }
	  else if((boot_tree->mod->s_opt->topo_search == SPR_MOVE) ||
		  (boot_tree->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR))
	    {
	      Speed_Spr_Loop(boot_tree);
	    }
	}
      else
	{
	  if(boot_tree->mod->s_opt->opt_subst_param || boot_tree->mod->s_opt->opt_bl)
	    Round_Optimize(boot_tree,boot_tree->data,ROUND_MAX);
	  else
	    Lk(boot_tree);
	}
      
      Free_Bip(boot_tree);

      Alloc_Bip(boot_tree);

      Match_Tip_Numbers(tree,boot_tree);
      
      Get_Bip(boot_tree->noeud[0],
	      boot_tree->noeud[0]->v[0],
	      boot_tree);


      Compare_Bip(tree,boot_tree);

      Br_Len_Involving_Invar(boot_tree);


      if(tree->io->print_boot_trees)
	{
	  s = Write_Tree(boot_tree);
	  PhyML_Fprintf(tree->io->fp_out_boot_tree,"%s\n",s);
	  Free(s);
          Print_Fp_Out_Lines(tree->io->fp_out_boot_stats,0,0,boot_tree,tree->io,replicate+1);
	}

      /*       rf = .0; */
      /*       For(j,2*tree->n_otu-3)  */
      /* 	rf += tree->t_edges[j]->bip_score; */


      PhyML_Printf("."); 
#ifndef QUIET
fflush(stdout);
#endif
      if(!((replicate+1)%tree->io->boot_prog_every))
	{
	  PhyML_Printf("] %4d/%4d\n  ",replicate+1,tree->mod->bootstrap);
	  if(replicate != tree->mod->bootstrap-1) PhyML_Printf("[");
	}

      if(boot_tree->mat) Free_Mat(boot_tree->mat);
      Free_Tree(boot_tree);      
      Free_Model(boot_mod);
    }

  if(((replicate)%tree->io->boot_prog_every)) PhyML_Printf("] %4d/%4d\n ",replicate,tree->mod->bootstrap);

  tree->lock_topo = 1; /* Topology should not be modified afterwards */

  if(tree->io->print_boot_trees)
    {
      fclose(tree->io->fp_out_boot_tree);
      fclose(tree->io->fp_out_boot_stats);
    }

  Free_Cseq(boot_data);
  Free(site_num);
}

/*********************************************************/

void Br_Len_Involving_Invar(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l *= (1.0-tree->mod->pinvar);
}

/*********************************************************/

void Br_Len_Not_Involving_Invar(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l /= (1.0-tree->mod->pinvar);
}

/*********************************************************/

void Getstring_Stdin(char *s)
{
  if(!fgets(s,T_MAX_LINE,stdin)) Exit("");
  if (strchr(s, '\n') != NULL)
    *strchr(s, '\n') = '\0';
}

/*********************************************************/

void Print_Freq(t_tree *tree)
{

  switch(tree->io->datatype)
    {
    case NT:
      {
	PhyML_Printf("A : %f\n",tree->mod->pi[0]);
	PhyML_Printf("C : %f\n",tree->mod->pi[1]);
	PhyML_Printf("G : %f\n",tree->mod->pi[2]);
	PhyML_Printf("T : %f\n",tree->mod->pi[3]);

	PhyML_Printf("U : %f\n",tree->mod->pi[4]);
	PhyML_Printf("M : %f\n",tree->mod->pi[5]);
	PhyML_Printf("R : %f\n",tree->mod->pi[6]);
	PhyML_Printf("W : %f\n",tree->mod->pi[7]);
	PhyML_Printf("S : %f\n",tree->mod->pi[8]);
	PhyML_Printf("Y : %f\n",tree->mod->pi[9]);
	PhyML_Printf("K : %f\n",tree->mod->pi[10]);
	PhyML_Printf("B : %f\n",tree->mod->pi[11]);
	PhyML_Printf("D : %f\n",tree->mod->pi[12]);
	PhyML_Printf("H : %f\n",tree->mod->pi[13]);
	PhyML_Printf("V : %f\n",tree->mod->pi[14]);
	PhyML_Printf("N : %f\n",tree->mod->pi[15]);
	break;
      }
    case AA:
      {
	PhyML_Printf("A : %f\n",tree->mod->pi[0]);
	PhyML_Printf("R : %f\n",tree->mod->pi[1]);
	PhyML_Printf("N : %f\n",tree->mod->pi[2]);
	PhyML_Printf("D : %f\n",tree->mod->pi[3]);
	PhyML_Printf("C : %f\n",tree->mod->pi[4]);
	PhyML_Printf("Q : %f\n",tree->mod->pi[5]);
	PhyML_Printf("E : %f\n",tree->mod->pi[6]);
	PhyML_Printf("G : %f\n",tree->mod->pi[7]);
	PhyML_Printf("H : %f\n",tree->mod->pi[8]);
	PhyML_Printf("I : %f\n",tree->mod->pi[9]);
	PhyML_Printf("L : %f\n",tree->mod->pi[10]);
	PhyML_Printf("K : %f\n",tree->mod->pi[11]);
	PhyML_Printf("M : %f\n",tree->mod->pi[12]);
	PhyML_Printf("F : %f\n",tree->mod->pi[13]);
	PhyML_Printf("P : %f\n",tree->mod->pi[14]);
	PhyML_Printf("S : %f\n",tree->mod->pi[15]);
	PhyML_Printf("T : %f\n",tree->mod->pi[16]);
	PhyML_Printf("W : %f\n",tree->mod->pi[17]);
	PhyML_Printf("Y : %f\n",tree->mod->pi[18]);
	PhyML_Printf("V : %f\n",tree->mod->pi[19]);

	PhyML_Printf("N : %f\n",tree->mod->pi[20]);
	break;
      }
    default : {break;}
    }
}

/*********************************************************/

phydbl Num_Derivatives_One_Param(phydbl (*func)(t_tree *tree), t_tree *tree,
				 phydbl f0, phydbl *param, phydbl stepsize,
				 phydbl *err, int precise)
{
  int i,j;
  phydbl errt,fac,hh,**a,ans;
  int n_iter;
  a = (phydbl **)mCalloc(11,sizeof(phydbl *));
  For(i,11) a[i] = (phydbl *)mCalloc(11,sizeof(phydbl));


  n_iter = 10; /* */

  ans  = .0;

  if(stepsize < SMALL) Warn_And_Exit("\n. h must be nonzero in Dfridr.");

  hh=stepsize;

  if(!precise)
    {
      *param   = *param+hh;
      a[0][0]  = (*func)(tree);
/*       printf("\n. f0=%f f1=%f hh=%G",f0,a[0][0],hh); */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *param   = *param-hh;

      ans =  a[0][0];

    }
  else
    {
      *param   = *param+hh;
      a[0][0]  = (*func)(tree);
      /*   *param   = *param-2*hh; */
      /*   a[0][0] -= (*func)(tree); */
      /*   a[0][0] /= (2.0*hh); */
      /*   *param   = *param+hh; */
      a[0][0]  -= f0;
      a[0][0]  /= hh;
      *param   = *param-hh;

      *err=1e30;
      for(i=1;i<n_iter;i++)
	{
	  hh /= 1.4;

	  /*       *param   = *param+hh; */
	  /*       a[0][i]  = (*func)(tree); */
	  /*       *param   = *param-2*hh; */
	  /*       a[0][i] -= (*func)(tree); */
	  /*       a[0][i] /= (2.0*hh); */
	  /*       *param   = *param+hh; */


	  *param   = *param+hh;
	  a[0][i]  = (*func)(tree);
	  /*   *param   = *param-2*hh; */
	  /*   a[0][i] -= (*func)(tree); */
	  /*   a[0][i] /= (2.0*hh); */
	  /*   *param   = *param+hh; */
	  a[0][i]  -= f0;
	  a[0][i]  /= hh;
	  *param   = *param-hh;


	  fac=1.4*1.4;
	  for (j=1;j<=i;j++)
	    {
	      a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
	      fac=1.4*1.4*fac;

	      errt=MAX(FABS(a[j][i]-a[j-1][i]),FABS(a[j][i]-a[j-1][i-1]));

	      if (errt <= *err)
		{
		  *err=errt;
		  ans=a[j][i];
		}
	    }

	  if(FABS(a[i][i]-a[i-1][i-1]) >= 2.0*(*err)) break;
	}
    }
  For(i,11) Free(a[i]);
  Free(a);

  return ans;
}

/*********************************************************/

int Num_Derivative_Several_Param(t_tree *tree, phydbl *param, int n_param, phydbl stepsize,
				  phydbl (*func)(t_tree *tree), phydbl *derivatives)
{
  int i;
  phydbl err,f0;

  f0 = (*func)(tree);

  For(i,n_param)
    {
      derivatives[i] = Num_Derivatives_One_Param(func,
						 tree,
						 f0,
						 param+i,
						 stepsize,
						 &err,
						 0
						 );
    }
  return 1;
}

/*********************************************************/

int Compare_Two_States(char *state1, char *state2, int state_size)
{
  /* 1 the two states are identical */
  /* 0 the two states are different */
  int i;

  For(i,state_size) if(state1[i] != state2[i]) break;

  return (i==state_size)?(1):(0);
}

/*********************************************************/

void Copy_One_State(char *from, char *to, int state_size)
{
  int i;
  For(i,state_size) to[i] = from[i];
}

/*********************************************************/

void Make_Custom_Model(model *mod)
{
  mod->rr            = (phydbl *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(phydbl));
  mod->rr_val        = (phydbl *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(phydbl));
  mod->rr_num        = (int *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(int *));
  mod->n_rr_per_cat  = (int *)mCalloc(mod->ns*(mod->ns-1)/2,sizeof(int));
}

/*********************************************************/

model *Make_Model_Basic()
{
  model *mod;

  mod                     = (model *)mCalloc(1,sizeof(model));
  mod->modelname          = (char *)mCalloc(T_MAX_NAME,sizeof(char));
  mod->custom_mod_string  = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
  mod->user_b_freq        = (phydbl *)mCalloc(T_MAX_OPTION,sizeof(phydbl));

  return mod;
}

/*********************************************************/

void Make_Model_Complete(model *mod)
{

  int i;

  if(mod->use_m4mod == YES)
    {
      M4_Make_Complete(mod->m4mod->n_h,mod->m4mod->n_o,mod->m4mod);
      mod->ns = mod->m4mod->n_o * mod->m4mod->n_h;
    }
  
  mod->pi                     = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  mod->pi_unscaled            = (phydbl *)mCalloc(mod->ns,sizeof(phydbl));
  mod->Pij_rr                 = (phydbl *)mCalloc(mod->n_catg*mod->ns*mod->ns,sizeof(phydbl));
  mod->gamma_r_proba          = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->gamma_rr               = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->gamma_rr_unscaled      = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));
  mod->qmat                   = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  mod->qmat_buff              = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  mod->eigen                  = (eigen *)Make_Eigen_Struct(mod->ns);
  mod->gamma_r_proba_unscaled = (phydbl *)mCalloc(mod->n_catg,sizeof(phydbl));

  For(i,mod->n_catg) mod->gamma_rr[i] = 1.0;
  For(i,mod->n_catg) mod->gamma_r_proba_unscaled[i] = 1.0;
  For(i,mod->n_catg) mod->gamma_r_proba[i] = 1.0/(phydbl)mod->n_catg;


  if(mod->n_rr_branch)
    {
      mod->rr_branch   = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
      mod->p_rr_branch = (phydbl *)mCalloc(mod->n_rr_branch,sizeof(phydbl));
    }
}

/*********************************************************/

void Copy_Dist(phydbl **cpy, phydbl **orig, int n)
{
  int i,j;
  For(i,n) For(j,n) cpy[i][j] = orig[i][j];
}

/*********************************************************/

model *Copy_Model(model *ori)
{
  model *cpy;

  cpy         = Make_Model_Basic();
  cpy->ns     = ori->ns;
  cpy->n_catg = ori->n_catg;

  Make_Model_Complete(cpy);
  if((ori->whichmodel == GTR) || (ori->whichmodel == CUSTOM)) 
    Make_Custom_Model(cpy);
  Record_Model(ori,cpy);
  cpy->m4mod = M4_Copy_M4_Model(ori, ori->m4mod);

  return cpy;
}

/*********************************************************/

void Record_Model(model *ori, model *cpy)
{
  int i;

  cpy->ns           = ori->ns;
  cpy->n_catg       = ori->n_catg;
  cpy->alpha_old    = ori->alpha_old;
  cpy->kappa_old    = ori->alpha_old;
  cpy->lambda_old   = ori->lambda_old;
  cpy->pinvar_old   = ori->pinvar_old;
  cpy->whichmodel   = ori->whichmodel;
  cpy->update_eigen = ori->update_eigen;
  cpy->kappa        = ori->kappa;
  cpy->alpha        = ori->alpha;
  cpy->lambda       = ori->lambda;
  cpy->bootstrap    = ori->bootstrap;
  cpy->invar        = ori->invar;
  cpy->pinvar       = ori->pinvar;
  cpy->n_diff_rr    = ori->n_diff_rr;
  cpy->l_min        = ori->l_min;
  cpy->l_max        = ori->l_max;
  cpy->log_l        = ori->log_l;


  cpy->free_mixt_rates = ori->free_mixt_rates;
  cpy->state_len   = ori->state_len;
  cpy->br_len_multiplier = ori->br_len_multiplier;


/*   if(ori->whichmodel == CUSTOM) */
/*     { */
/*       For(i,ori->ns*(ori->ns-1)/2) cpy->rr_num[i] = ori->rr_num[i]; */

/*       For(i,ori->ns*(ori->ns-1)/2) */
/* 	{ */
/* 	  cpy->rr_val[i]  = ori->rr_val[i]; */
/* 	  cpy->rr[i]      = cpy->rr[i]; */
/* 	} */
/*     } */
  
  if((ori->whichmodel == CUSTOM) || (ori->whichmodel == GTR))
    {
      For(i,ori->ns*(ori->ns-1)/2)
	{
	  cpy->rr_num[i]       = ori->rr_num[i];
	  cpy->rr_val[i]       = ori->rr_val[i];
	  cpy->rr[i]           = ori->rr[i];
	  cpy->n_rr_per_cat[i] = ori->n_rr_per_cat[i];
	}
    }
  


  For(i,cpy->ns)
    {
      cpy->pi[i]          = ori->pi[i];
      cpy->pi_unscaled[i] = ori->pi_unscaled[i];
      cpy->user_b_freq[i] = ori->user_b_freq[i];
    }
  
  For(i,cpy->ns*cpy->ns) cpy->qmat[i] = ori->qmat[i];

  For(i,cpy->n_catg)
    {
      cpy->gamma_r_proba[i]          = ori->gamma_r_proba[i];
      cpy->gamma_rr[i]               = ori->gamma_rr[i];
      cpy->gamma_r_proba_unscaled[i] = ori->gamma_r_proba_unscaled[i];
      cpy->gamma_rr_unscaled[i]      = ori->gamma_rr_unscaled[i];
    }
  
  cpy->use_m4mod = ori->use_m4mod;

  cpy->eigen->size = ori->eigen->size;
  For(i,2*ori->ns)       cpy->eigen->space[i]       = ori->eigen->space[i];
  For(i,2*ori->ns)       cpy->eigen->space_int[i]   = ori->eigen->space_int[i];
  For(i,ori->ns)         cpy->eigen->e_val[i]       = ori->eigen->e_val[i];
  For(i,ori->ns)         cpy->eigen->e_val_im[i]    = ori->eigen->e_val_im[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect[i]    = ori->eigen->r_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->r_e_vect_im[i] = ori->eigen->r_e_vect_im[i];
  For(i,ori->ns*ori->ns) cpy->eigen->l_e_vect[i]    = ori->eigen->l_e_vect[i];
  For(i,ori->ns*ori->ns) cpy->eigen->q[i]           = ori->eigen->q[i];
}

/*********************************************************/

option *Make_Input()
{
  int i;
  option* io               = (option *)mCalloc(1,sizeof(option));

  io->in_align_file        = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->in_tree_file         = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_tree_file        = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trees_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_tree_file   = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_boot_stats_file  = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_stats_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_lk_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_ps_file          = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->out_trace_file       = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->nt_or_cd             = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->run_id_string        = (char *)mCalloc(T_MAX_OPTION,sizeof(char));
  io->clade_list_file      = (char *)mCalloc(T_MAX_FILE,sizeof(char));
  io->alphabet             = (char **)mCalloc(T_MAX_ALPHABET,sizeof(char *));
  For(i,T_MAX_ALPHABET) io->alphabet[i] = (char *)mCalloc(T_MAX_STATE,sizeof(char ));
  io->treelist             = (t_treelist *)mCalloc(1,sizeof(t_treelist));
  io->mcmc                 = (t_mcmc *)MCMC_Make_MCMC_Struct();
  io->rates                = (t_rate *)RATES_Make_Rate_Struct(-1);
  return io;
}

/*********************************************************/

void Set_Defaults_Input(option* io)
{
  io->fp_in_align                = NULL;
  io->fp_in_tree                 = NULL;
  io->fp_out_tree                = NULL;
  io->fp_out_trees               = NULL;
  io->fp_out_boot_tree           = NULL;
  io->fp_out_boot_stats          = NULL;
  io->fp_out_stats               = NULL;
  io->long_tax_names             = NULL;
  io->short_tax_names            = NULL;
  io->lon                        = NULL;
  io->lat                        = NULL;
  io->z_scores                   = NULL;

  io->tree                       = NULL;
  io->mod                        = NULL;
  strcpy(io->nt_or_cd,"nucleotides");
  io->n_data_sets                = 1;
  io->interleaved                = 1;
  io->in_tree                    = 0;
  io->out_tree_file_open_mode    = 1;
  io->out_stats_file_open_mode   = 1;
  io->init_len                   = -1;
  io->n_otu                      = -1;
  io->n_data_set_asked           = -1;
  io->print_boot_trees           = 1;
  io->n_part                     = 1;
  io->ratio_test		 = 4;
  io->multigene                  = 0;
  io->config_multigene           = 0;
  io->curr_interface             = 0;
  io->r_seed                     = -1;
  io->collapse_boot              = 0;
  io->random_boot_seq_order      = 1;
  io->print_trace                = 0;
  io->print_site_lnl             = 0;
  io->m4_model                   = NO;
  io->rm_ambigu                  = 0;
  io->appebr_run_ID              = 0;
  io->quiet                      = 0;
  io->datatype                   = NT;
  io->colalias                   = YES;
  io->data_file_format           = PHYLIP;
  io->tree_file_format           = PHYLIP;
  io->boot_prog_every            = 20;
  io->mem_question               = YES;
  io->do_alias_subpatt           = NO;
  io->lk_approx                  = EXACT;

  MCMC_Init_MCMC_Struct(NULL,io->mcmc);
  RATES_Init_Rate_Struct(io->rates,NULL,-1);
  io->rates->model               = THORNE;
}

/*********************************************************/

void Set_Defaults_Model(model *mod)
{
  strcpy(mod->modelname,"HKY85");
  strcpy(mod->custom_mod_string,"000000");
  mod->whichmodel              = HKY85;
  mod->n_catg                  = 4;
  mod->kappa                   = 4.0;
  mod->alpha                   = 1.0;
  mod->lambda                  = 1.0;
  mod->bootstrap               = 0;
  mod->invar                   = 0;
  mod->pinvar                  = 0.0;
  mod->ns                      = 4;
  mod->n_diff_rr               = 0;
  mod->use_m4mod               = NO;
  mod->n_rr_branch             = 0;
  mod->rr_branch_alpha         = 0.1;
  mod->gamma_median            = 0;
  mod->state_len               = 1;
  mod->m4mod                   = NULL;
  mod->rr                      = NULL;
  mod->rr_val                  = NULL;
  mod->n_rr_per_cat            = NULL;
  mod->io                      = NULL;
  mod->log_l                   = NO;
  mod->free_mixt_rates         = NO;
  mod->gamma_mgf_bl            = NO;
  mod->br_len_multiplier       = 1.0;

#ifndef PHYTIME
  mod->l_min = 1.E-8;
  mod->l_max = 100.0;
#else
  mod->l_min = 1.E-8;
  mod->l_max = 20.00;
#endif
}

/*********************************************************/

void Set_Defaults_Optimiz(optimiz *s_opt)
{
  s_opt->print                = 1;
  s_opt->last_opt             = 1;
  s_opt->opt_subst_param      = 1;
  s_opt->opt_alpha            = 1;
  s_opt->opt_kappa            = 1;
  s_opt->opt_bl               = 1;
  s_opt->opt_lambda           = 0;
  s_opt->opt_pinvar           = 0;
  s_opt->opt_cov_delta        = 0;
  s_opt->opt_cov_alpha        = 0;
  s_opt->opt_cov_free_rates   = 0;
  s_opt->opt_rr               = 0;
  s_opt->init_lk              = UNLIKELY;
  s_opt->n_it_max             = 1000;
  s_opt->opt_topo             = 1;
  s_opt->topo_search          = NNI_MOVE;
  s_opt->random_input_tree    = 0;
  s_opt->n_rand_starts        = 5;
  s_opt->brent_it_max         = 500;
  s_opt->steph_spr            = 1;
  s_opt->user_state_freq      = 0;
  s_opt->min_diff_lk_local    = 1.E-04;
  s_opt->min_diff_lk_global   = 1.E-03;
  s_opt->min_diff_lk_move     = 1.E-02;
  s_opt->p_moves_to_examine   = 0.15;
  s_opt->fast_nni             = 0;
  s_opt->greedy               = 0;
  s_opt->general_pars         = 0;
  s_opt->tree_size_mult       = 1;
  s_opt->opt_five_branch      = 1;

  s_opt->pars_thresh          = 5;

  s_opt->hybrid_thresh        = 0;
  s_opt->quickdirty           = 0;
  s_opt->spr_pars             = 1;
  s_opt->spr_lnL              = 0;
  s_opt->min_depth_path       = 0;
  s_opt->max_depth_path       = 20;
  s_opt->deepest_path         = 20;
  s_opt->max_delta_lnL_spr    = 50.;
  s_opt->br_len_in_spr        = 10;
  s_opt->opt_free_mixt_rates  = YES;
  s_opt->constrained_br_len   = NO;
  s_opt->opt_gamma_br_len     = NO;

  s_opt->wim_n_rgrft          = -1;
  s_opt->wim_n_globl          = -1;
  s_opt->wim_max_dist         = -1;
  s_opt->wim_n_optim          = -1;
  s_opt->wim_n_best           = -1;
  s_opt->wim_inside_opt       =  0;
}

/*********************************************************/

void Test_Node_Table_Consistency(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-2)
    {
      if(tree->noeud[i]->num != i)
	{
	  PhyML_Printf("\n. Node table is not consistent with node numbers.");
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
    }
}

/*********************************************************/

void Get_Bip(t_node *a, t_node *d, t_tree *tree)
{
  int i,j;
  t_node *tmp;
  int swapped;

  if(d->tax)
    {
      if(d->common)
	{
	  d->bip_node[0] = (t_node **)mCalloc(1,sizeof(t_node *));
	  d->bip_node[0][0] = d;
	  d->bip_size[0]    = 1;

	  For(i,3)
	    {
	      if(a->v[i] == d)
		{
		  a->bip_size[i] = 0;
		  For(j,tree->n_otu)
		    {
		      if(strcmp(tree->noeud[j]->name,d->name))
			{
			  a->bip_node[i] = (t_node **)realloc(a->bip_node[i],(a->bip_size[i]+1)*sizeof(t_node *));
			  a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
			  a->bip_size[i]++;
			}
		    }
		  
		  /* Sort bipartition */
		  do
		    {
		      swapped = NO;
		      For(j,a->bip_size[i]-1)
			{
			  if(a->bip_node[i][j]->num > a->bip_node[i][j+1]->num)
			    {
			      swapped = YES;
			      tmp                 = a->bip_node[i][j];
			      a->bip_node[i][j]   = a->bip_node[i][j+1];
			      a->bip_node[i][j+1] = tmp;
			    }
			}
		    }while(swapped == YES);
		  
		  break;
		  
		}
	    }
	}
      return;
    }
  else
    {
      int k;
      int d_a;

      d_a = -1;

      For(i,3)
	{
	  if(d->v[i] != a) Get_Bip(d,d->v[i],tree);
	  else if(d->v[i] == a) d_a = i;
	}

      d->bip_size[d_a] = 0;
      For(i,3)
	if(d->v[i] != a)
	  {
	    For(j,3)
	      {
		if(d->v[i]->v[j] == d)
		  {
		    For(k,d->v[i]->bip_size[j])
		      {
			d->bip_node[d_a] = (t_node **)realloc(d->bip_node[d_a],(d->bip_size[d_a]+1)*sizeof(t_node *));
			d->bip_node[d_a][d->bip_size[d_a]] = d->v[i]->bip_node[j][k];
			d->bip_size[d_a]++;
		      }
		    break;
		  }
	      }
	  }
      
      do
	{
	  swapped = NO;
	  For(j,d->bip_size[d_a]-1)
	    {
	      if(d->bip_node[d_a][j]->num > d->bip_node[d_a][j+1]->num)
		{
		  swapped = YES;
		  tmp                   = d->bip_node[d_a][j];
		  d->bip_node[d_a][j]   = d->bip_node[d_a][j+1];
		  d->bip_node[d_a][j+1] = tmp;
		}
	    }
	}while(swapped == YES);
	
      
      For(i,3)
	if(a->v[i] == d)
	  {
	    a->bip_size[i] = 0;
	    For(j,tree->n_otu)
	      {
		For(k,d->bip_size[d_a])
		  {
		    if(d->bip_node[d_a][k] == tree->noeud[j])
		      break;
		  }
		
		if((k == d->bip_size[d_a]) && (tree->noeud[j]->common))
		  {
		    a->bip_node[i] = (t_node **)realloc(a->bip_node[i],(a->bip_size[i]+1)*sizeof(t_node *));
		    a->bip_node[i][a->bip_size[i]] = tree->noeud[j];
		    a->bip_size[i]++;
		  }
	      }

	    do
	      {
		swapped = NO;
		For(j,a->bip_size[i]-1)
		  {
		    if(a->bip_node[i][j]->num > a->bip_node[i][j+1]->num)
		      {
			swapped = YES;
			tmp                 = a->bip_node[i][j];
			a->bip_node[i][j]   = a->bip_node[i][j+1];
			a->bip_node[i][j+1] = tmp;
		      }
		  }
	      }while(swapped == YES);
	    
	    if(a->bip_size[i] != tree->n_otu - d->bip_size[d_a])
	      {
		PhyML_Printf("%d %d \n",a->bip_size[i],tree->n_otu - d->bip_size[d_a]);
		Warn_And_Exit("\n. Problem in counting bipartitions \n");
	      }
	    break;
	  }
    }
}

/*********************************************************/

void Alloc_Bip(t_tree *tree)
{
  int i;

  if(tree->has_bip) return;

  tree->has_bip = YES;

  For(i,2*tree->n_otu-2)
    {
      tree->noeud[i]->bip_size = (int *)mCalloc(3,sizeof(int));
      tree->noeud[i]->bip_node = (t_node ***)mCalloc(3,sizeof(t_node **));

/*       For(j,3) */
/* 	{ */
/* 	  tree->noeud[i]->bip_node[j] = (t_node **)mCalloc(tree->n_otu,sizeof(t_node *)); */
/* 	} */
    }
}

/*********************************************************/

int Sort_Phydbl_Increase(const void *a, const void *b)
{
  if((*(phydbl *)(a)) <= (*(phydbl *)(b))) return -1;
  else return 1;
}

/*********************************************************/

int Sort_String(const void *a, const void *b)
{
  return(strcmp((*(const char **)(a)), (*(const char **)(b))));
}

/*********************************************************/

void Compare_Bip(t_tree *tree1, t_tree *tree2)
{
  int i,j,k;
  t_edge *b1,*b2;
/*   char **bip1,**bip2; */
/*   int *bip1,*bip2; */
  t_node **bip1, **bip2;
  int bip_size1, bip_size2, bip_size;



  For(i,2*tree1->n_otu-3)
    {
      b1 = tree1->t_edges[i];     
      bip_size1 = MIN(b1->left->bip_size[b1->l_r],b1->rght->bip_size[b1->r_l]);
      
      if(bip_size1 > 1)
	{
	  For(j,2*tree2->n_otu-3)
	    {
	      b2 = tree2->t_edges[j];	      
	      bip_size2 = MIN(b2->left->bip_size[b2->l_r],b2->rght->bip_size[b2->r_l]);

	      if(bip_size2 > 1)
		{
		  if(bip_size1 == bip_size2)
		    {
		      bip_size = bip_size1;

		      if(b1->left->bip_size[b1->l_r] == b1->rght->bip_size[b1->r_l])
			{
/* 			  if(b1->left->bip_name[b1->l_r][0][0] < b1->rght->bip_name[b1->r_l][0][0]) */
			  if(b1->left->bip_node[b1->l_r][0]->num < b1->rght->bip_node[b1->r_l][0]->num)
			    {
/* 			      bip1 = b1->left->bip_name[b1->l_r]; */
			      bip1 = b1->left->bip_node[b1->l_r];
			    }
			  else
			    {
/* 			      bip1 = b1->rght->bip_name[b1->r_l]; */
			      bip1 = b1->rght->bip_node[b1->r_l];
			    }
			}
		      else if(b1->left->bip_size[b1->l_r] < b1->rght->bip_size[b1->r_l])
			{
/* 			  bip1 = b1->left->bip_name[b1->l_r]; */
			  bip1 = b1->left->bip_node[b1->l_r];
			}
		      else
			{
/* 			  bip1 = b1->rght->bip_name[b1->r_l]; */
			  bip1 = b1->rght->bip_node[b1->r_l];
			}


		      if(b2->left->bip_size[b2->l_r] == b2->rght->bip_size[b2->r_l])
			{
/* 			  if(b2->left->bip_name[b2->l_r][0][0] < b2->rght->bip_name[b2->r_l][0][0]) */
			  if(b2->left->bip_node[b2->l_r][0]->num < b2->rght->bip_node[b2->r_l][0]->num)
			    {
/* 			      bip2 = b2->left->bip_name[b2->l_r]; */
			      bip2 = b2->left->bip_node[b2->l_r];
			    }
			  else
			    {
/* 			      bip2 = b2->rght->bip_name[b2->r_l]; */
			      bip2 = b2->rght->bip_node[b2->r_l];
			    }
			}
		      else if(b2->left->bip_size[b2->l_r] < b2->rght->bip_size[b2->r_l])
			{
/* 			  bip2 = b2->left->bip_name[b2->l_r]; */
			  bip2 = b2->left->bip_node[b2->l_r];
			}
		      else
			{
/* 			  bip2 = b2->rght->bip_name[b2->r_l]; */
			  bip2 = b2->rght->bip_node[b2->r_l];
			}

		      if(bip_size == 1) Warn_And_Exit("\n. Problem in Compare_Bip\n");

		      For(k,bip_size)
			{
/* 			  if(strcmp(bip1[k],bip2[k])) break; */
			  if(bip1[k]->num != bip2[k]->num) break;
			}

		      if(k == bip_size)
			{
			  b1->bip_score++;
			  b2->bip_score++;
			  break;
			}
		    }
		}
	    }
	}
    }
}

/*********************************************************/

void Match_Tip_Numbers(t_tree *tree1, t_tree *tree2)
{
  int i,j;

  For(i,tree1->n_otu)
    {
      For(j,tree2->n_otu)
	{
	  if(!strcmp(tree1->noeud[i]->name,tree2->noeud[j]->name))
	    {
	      tree2->noeud[j]->num = tree1->noeud[i]->num;
	      break;
	    }
	}
    }

}

/*********************************************************/

void Test_Multiple_Data_Set_Format(option *io)
{
  char *line;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  io->n_trees = 0;

  while(fgets(line,T_MAX_LINE,io->fp_in_tree)) if(strstr(line,";")) io->n_trees++;

  Free(line);

  if((io->mod->bootstrap > 1) && (io->n_trees > 1))
    Warn_And_Exit("\n. Bootstrap option is not allowed with multiple input trees !\n");

  rewind(io->fp_in_tree);

  return;
}

/*********************************************************/

int Are_Compatible(char *statea, char *stateb, int stepsize, int datatype)
{
  int i,j;
  char a,b;

  if(datatype == NT)
    {
      For(i,stepsize)
	{
	  a = statea[i];
	  For(j,stepsize)
	    {
	      b = stateb[j];

	      switch(a)
		{
		case 'A':
		  {
		    switch(b)
		      {
		      case 'A' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'G':
		  {
		    switch(b)
		      {
		      case 'G' :
		      case 'R' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'C':
		  {
		    switch(b)
		      {
		      case 'C' :
		      case 'M' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'T':
		  {
		    switch(b)
		      {
		      case 'T' :
		      case 'W' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'X' :
			{b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'M' :
		  {
		    switch(b)
		      {
		      case 'M' :
		      case 'A' :
		      case 'C' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' :
			{b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'R' :
		  {
		    switch(b)
		      {
		      case 'R' :
		      case 'A' :
		      case 'G' :
		      case 'M' :
		      case 'W' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'W' :
		  {
		    switch(b)
		      {
		      case 'W' :
		      case 'A' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'S' :
		  {
		    switch(b)
		      {
		      case 'S' :
		      case 'C' :
		      case 'G' :
		      case 'M' :
		      case 'R' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'Y' :
		  {
		    switch(b)
		      {
		      case 'Y' :
		      case 'C' :
		      case 'T' :
		      case 'M' :
		      case 'W' :
		      case 'S' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }

		case 'K' :
		  {
		    switch(b)
		      {
		      case 'K' :
		      case 'G' :
		      case 'T' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'B' :
		  {
		    switch(b)
		      {
		      case 'B' :
		      case 'C' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'D' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'D' :
		  {
		    switch(b)
		      {
		      case 'D' :
		      case 'A' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'H' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'H' :
		  {
		    switch(b)
		      {
		      case 'H' :
		      case 'A' :
		      case 'C' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'V' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'V' :
		  {
		    switch(b)
		      {
		      case 'V' :
		      case 'A' :
		      case 'C' :
		      case 'G' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'X' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		case 'X' :
		  {
		    switch(b)
		      {
		      case 'X' :
		      case 'A' :
		      case 'C' :
		      case 'G' :
		      case 'T' :
		      case 'M' :
		      case 'R' :
		      case 'W' :
		      case 'S' :
		      case 'Y' :
		      case 'K' :
		      case 'B' :
		      case 'D' :
		      case 'H' :
		      case 'V' : {b=b; break;}
		      default : return 0;
		      }
		    break;
		  }
		default :
		  {
                      PhyML_Printf("\n. Err. in Are_Compatible\n");
                      PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
                      PhyML_Printf("  correspond to existing nucleotides.\n");
                      Warn_And_Exit("\n");
                      return 0;
		  }
		}
	    }
	}
    }
  else if(datatype == AA)
    {
      a = statea[0]; b = stateb[0];
      switch(a)
	{
	case 'A' :
	  {
	    switch(b)
	      {
	      case 'A' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'R' :
	  {
	    switch(b)
	      {
	      case 'R' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'N' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'B' :
	  {
	    switch(b)
	      {
	      case 'N' :
	      case 'B' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'D' :
	  {
	    switch(b)
	      {
	      case 'D' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'C' :
	  {
	    switch(b)
	      {
	      case 'C' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Q' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Z' :
	  {
	    switch(b)
	      {
	      case 'Q' :
	      case 'Z' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'E' :
	  {
	    switch(b)
	      {
	      case 'E' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'G' :
	  {
	    switch(b)
	      {
	      case 'G' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'H' :
	  {
	    switch(b)
	      {
	      case 'H' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'I' :
	  {
	    switch(b)
	      {
	      case 'I' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'L' :
	  {
	    switch(b)
	      {
	      case 'L' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'K' :
	  {
	    switch(b)
	      {
	      case 'K' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'M' :
	  {
	    switch(b)
	      {
	      case 'M' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'F' :
	  {
	    switch(b)
	      {
	      case 'F' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'P' :
	  {
	    switch(b)
	      {
	      case 'P' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'S' :
	  {
	    switch(b)
	      {
	      case 'S' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'T' :
	  {
	    switch(b)
	      {
	      case 'T' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'W' :
	  {
	    switch(b)
	      {
	      case 'W' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'Y' :
	  {
	    switch(b)
	      {
	      case 'Y' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'V' :
	  {
	    switch(b)
	      {
	      case 'V' :
	      case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	case 'X' :
	  {
	    switch(b)
	      {
	      case 'A':case 'R':case 'N' :case 'B' :case 'D' :
	      case 'C':case 'Q':case 'Z' :case 'E' :case 'G' :
	      case 'H':case 'I':case 'L' :case 'K' :case 'M' :
	      case 'F':case 'P':case 'S' :case 'T' :case 'W' :
	      case 'Y':case 'V': case 'X' : {b=b; break;}
	      default : return 0;
	      }
	    break;
	  }
	default :
	  {
	    PhyML_Printf("\n. Err. in Are_Compatible\n");
            PhyML_Printf("\n. Please check that characters `%c` and `%c`\n",a,b);
            PhyML_Printf("  correspond to existing amino-acids.\n");
            Warn_And_Exit("\n");
	    return 0;
	  }
	}
    }
  else if(datatype == GENERIC)    
    {
      if(Is_Ambigu(statea,GENERIC,stepsize) || Is_Ambigu(stateb,GENERIC,stepsize)) return 1;
      else
	{
	  int a,b;
	  char format[6];      
	  
	  sprintf(format,"%%%dd",stepsize);      

	  if(!sscanf(statea,format,&a))
	    {	    
	      PhyML_Printf("\n. statea = %s",statea);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  if(!sscanf(stateb,format,&b))
	    {	    
	      PhyML_Printf("\n. statea = %s",stateb);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  
/* 	  PhyML_Printf("\n. %s %d a=%d b=%d ",__FILE__,__LINE__,a,b);  */

	  if(a == b) return 1;
	}
      return 0;
    }

  return 1;
}

/*********************************************************/

void Hide_Ambiguities(calign *data)
{
  int i;
  For(i,data->crunch_len) if(data->ambigu[i]) data->wght[i] = 0;
}

/*********************************************************/

void Copy_Tree(t_tree *ori, t_tree *cpy)
{
  int i,j;


  For(i,2*ori->n_otu-2)
    {
      For(j,3)
	{
	  if(ori->noeud[i]->v[j])
	    {
	      cpy->noeud[i]->v[j] = cpy->noeud[ori->noeud[i]->v[j]->num];
	      cpy->noeud[i]->l[j] = ori->noeud[i]->l[j];
	      cpy->noeud[i]->b[j] = cpy->t_edges[ori->noeud[i]->b[j]->num];
	    }
	  else
	    {
	      cpy->noeud[i]->v[j] = NULL;
	      cpy->noeud[i]->b[j] = NULL;
	    }
	}
    }


  For(i,2*ori->n_otu-3) 
    {
      cpy->t_edges[i]->l    = ori->t_edges[i]->l;
      cpy->t_edges[i]->left = cpy->noeud[ori->t_edges[i]->left->num];
      cpy->t_edges[i]->rght = cpy->noeud[ori->t_edges[i]->rght->num];
      cpy->t_edges[i]->l_v1 = ori->t_edges[i]->l_v1;
      cpy->t_edges[i]->l_v2 = ori->t_edges[i]->l_v2;
      cpy->t_edges[i]->r_v1 = ori->t_edges[i]->r_v1;
      cpy->t_edges[i]->r_v2 = ori->t_edges[i]->r_v2;
      cpy->t_edges[i]->l_r  = ori->t_edges[i]->l_r;
      cpy->t_edges[i]->r_l  = ori->t_edges[i]->r_l;
    }


  For(i,ori->n_otu)
    {
      cpy->noeud[i]->tax = 1;
      strcpy(cpy->noeud[i]->name,ori->noeud[i]->name);
    }


  cpy->num_curr_branch_available = 0;
/*   Connect_Edges_To_Nodes_Recur(cpy->noeud[0],cpy->noeud[0]->v[0],cpy); */
/*   Update_Dirs(cpy); */
}

/*********************************************************/

void Prune_Subtree(t_node *a, t_node *d, t_edge **target, t_edge **residual, t_tree *tree)
{
  t_node *v1, *v2;
  t_edge *b1, *b2;
  int dir_v1, dir_v2;
  int i;
/*   phydbl ***buff_p_lk; */
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars, *buff_pars, *buff_p_lk_loc, *buff_patt_id;
  unsigned int *buff_ui;
  short int *buff_p_lk_tip;

  if(a->tax)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  dir_v1 = dir_v2 = -1;
  For(i,3)
    {
      if(a->v[i] != d)
	{
	  if(dir_v1 < 0) dir_v1 = i;
	  else           dir_v2 = i;
	}
    }

  if(a->v[dir_v1]->num < a->v[dir_v2]->num)
    {
      v1 = a->v[dir_v1];
      v2 = a->v[dir_v2];
      b1 = a->b[dir_v1];
      b2 = a->b[dir_v2];
    }
  else
    {
      v1 = a->v[dir_v2];
      v2 = a->v[dir_v1];
      b1 = a->b[dir_v2];
      b2 = a->b[dir_v1];
    }

  if(v1->tax && v2->tax) PhyML_Printf("\n. Pruning is meaningless here.\n");

  a->v[dir_v1] = NULL;
  a->v[dir_v2] = NULL;
  a->b[dir_v1] = NULL;
  a->b[dir_v2] = NULL;

  if(v1 == b1->left)
    {
      b1->rght = v2;

      if(v2 == b2->left)
	{
	  buff_p_lk            = b1->p_lk_rght;
	  b1->p_lk_rght        = b2->p_lk_left;
	  b2->p_lk_left        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_r;
	  b1->p_lk_tip_r       = b2->p_lk_tip_l;
	  b2->p_lk_tip_l       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_rght;
	  b1->sum_scale_rght   = b2->sum_scale_left;
	  b2->sum_scale_left   = buff_scale;

	  buff_scale             = b1->sum_scale_rght_cat;
	  b1->sum_scale_rght_cat = b2->sum_scale_left_cat;
	  b2->sum_scale_left_cat = buff_scale;

	  buff_pars            = b1->pars_r;
	  b1->pars_r           = b2->pars_l;
	  b2->pars_l           = buff_pars;

	  buff_ui              = b1->ui_r;
	  b1->ui_r             = b2->ui_l;
	  b2->ui_l             = buff_ui;

	  buff_p_pars          = b1->p_pars_r;
	  b1->p_pars_r         = b2->p_pars_l;
	  b2->p_pars_l         = buff_p_pars;

	  buff_p_lk_loc        = b1->p_lk_loc_rght;
	  b1->p_lk_loc_rght    = b2->p_lk_loc_left;
	  b2->p_lk_loc_left    = buff_p_lk_loc;

	  buff_patt_id         = b1->patt_id_rght;
	  b1->patt_id_rght     = b2->patt_id_left;
	  b2->patt_id_left     = buff_patt_id;

	}
      else
	{
	  buff_p_lk            = b1->p_lk_rght; /* b1->p_lk_rght = NULL if b1->rght->tax */
	  b1->p_lk_rght        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */ 
	  b2->p_lk_rght        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_r;
	  b1->p_lk_tip_r       = b2->p_lk_tip_r;
	  b2->p_lk_tip_r       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_rght;
	  b1->sum_scale_rght = b2->sum_scale_rght;
	  b2->sum_scale_rght = buff_scale;

	  buff_pars            = b1->pars_r;
	  b1->pars_r           = b2->pars_r;
	  b2->pars_r           = buff_pars;

	  buff_ui              = b1->ui_r;
	  b1->ui_r             = b2->ui_r;
	  b2->ui_r             = buff_ui;

	  buff_p_pars          = b1->p_pars_r;
	  b1->p_pars_r         = b2->p_pars_r;
	  b2->p_pars_r         = buff_p_pars;

	  buff_p_lk_loc        = b1->p_lk_loc_rght;
	  b1->p_lk_loc_rght    = b2->p_lk_loc_rght;
	  b2->p_lk_loc_rght    = buff_p_lk_loc;

	  buff_patt_id         = b1->patt_id_rght;
	  b1->patt_id_rght     = b2->patt_id_rght;
	  b2->patt_id_rght     = buff_patt_id;

	}
    }
  else
    {
      b1->left = v2;

      if(v2 == b2->left)
	{
	  buff_p_lk            = b1->p_lk_left;
	  b1->p_lk_left        = b2->p_lk_left;
	  b2->p_lk_left        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_l;
	  b1->p_lk_tip_l       = b2->p_lk_tip_l;
	  b2->p_lk_tip_l       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_left;
	  b1->sum_scale_left = b2->sum_scale_left;
	  b2->sum_scale_left = buff_scale;

	  buff_scale             = b1->sum_scale_left_cat;
	  b1->sum_scale_left_cat = b2->sum_scale_left_cat;
	  b2->sum_scale_left_cat = buff_scale;

	  buff_pars            = b1->pars_l;
	  b1->pars_l           = b2->pars_l;
	  b2->pars_l           = buff_pars;

	  buff_ui              = b1->ui_l;
	  b1->ui_l             = b2->ui_l;
	  b2->ui_l             = buff_ui;

	  buff_p_pars          = b1->p_pars_l;
	  b1->p_pars_l         = b2->p_pars_l;
	  b2->p_pars_l         = buff_p_pars;

	  buff_p_lk_loc        = b1->p_lk_loc_left;
	  b1->p_lk_loc_left    = b2->p_lk_loc_left;
	  b2->p_lk_loc_left    = buff_p_lk_loc;

	  buff_patt_id         = b1->patt_id_left;
	  b1->patt_id_left     = b2->patt_id_left;
	  b2->patt_id_left     = buff_patt_id;

	}
      else
	{
	  buff_p_lk            = b1->p_lk_left;
	  b1->p_lk_left        = b2->p_lk_rght; /* b2->p_lk_rght = NULL if b2->rght->tax */
	  b2->p_lk_rght        = buff_p_lk;

	  buff_p_lk_tip        = b1->p_lk_tip_l;
	  b1->p_lk_tip_l       = b2->p_lk_tip_r;
	  b2->p_lk_tip_r       = buff_p_lk_tip;

	  buff_scale           = b1->sum_scale_left;
	  b1->sum_scale_left = b2->sum_scale_rght;
	  b2->sum_scale_rght = buff_scale;

	  buff_scale             = b1->sum_scale_left_cat;
	  b1->sum_scale_left_cat = b2->sum_scale_rght_cat;
	  b2->sum_scale_rght_cat = buff_scale;

	  buff_pars            = b1->pars_l;
	  b1->pars_l           = b2->pars_r;
	  b2->pars_r           = buff_pars;

	  buff_ui              = b1->ui_l;
	  b1->ui_l             = b2->ui_r;
	  b2->ui_r             = buff_ui;

	  buff_p_pars          = b1->p_pars_l;
	  b1->p_pars_l         = b2->p_pars_r;
	  b2->p_pars_r         = buff_p_pars;

	  buff_p_lk_loc        = b1->p_lk_loc_left;
	  b1->p_lk_loc_left    = b2->p_lk_loc_rght;
	  b2->p_lk_loc_rght    = buff_p_lk_loc;

	  buff_patt_id         = b1->patt_id_left;
	  b1->patt_id_left     = b2->patt_id_rght;
	  b2->patt_id_rght     = buff_patt_id;
	}
    }

  For(i,3)
    if(v2->v[i] == a)
      {
	v2->v[i] = v1;
	v2->b[i] = b1;
	break;
      }

#ifdef DEBUG
  if(i == 3)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  For(i,3)
    if(v1->v[i] == a)
      {
	v1->v[i] = v2;
	break;
      }

#ifdef DEBUG
  if(i == 3)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif

  b1->l += b2->l;

  (v1 == b1->left)?
    (Make_Edge_Dirs(b1,v1,v2,tree)):
    (Make_Edge_Dirs(b1,v2,v1,tree));

  if(target)   (*target)   = b1;
  if(residual) (*residual) = b2;


#ifdef DEBUG
  if(b1->left->tax)
    {
      PhyML_Printf("\n. b1->left->num = %d",b1->left->num);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
#endif


}

/*********************************************************/

void Graft_Subtree(t_edge *target, t_node *link, t_edge *residual, t_tree *tree)
{
  t_node *v1, *v2;
  int i, dir_v1, dir_v2;
  phydbl *buff_p_lk;
  int *buff_scale;
  int *buff_p_pars, *buff_pars, *buff_p_lk_loc, *buff_patt_id; 
  short int *buff_p_lk_tip;
  unsigned int *buff_ui;
  t_edge *b_up;

  dir_v1 = dir_v2 = -1;
  b_up = NULL;
  For(i,3)
    {
      if(!link->v[i])
	{
	  if(dir_v1 < 0) dir_v1 = i;
	  else           dir_v2 = i;
	}
      else b_up = link->b[i];
    }

  if(target->left->num < target->rght->num)
    {
      v1                           = target->left;
      v2                           = target->rght;

      buff_p_lk                    = residual->p_lk_rght;
      residual->p_lk_rght          = target->p_lk_rght;
      target->p_lk_rght            = buff_p_lk;

      buff_p_lk_tip                = residual->p_lk_tip_r;
      residual->p_lk_tip_r         = target->p_lk_tip_r;
      target->p_lk_tip_r           = buff_p_lk_tip;

      buff_scale                   = residual->sum_scale_rght;
      residual->sum_scale_rght     = target->sum_scale_rght;
      target->sum_scale_rght       = buff_scale;

      buff_pars                    = residual->pars_r;
      residual->pars_r             = target->pars_r;
      target->pars_r               = buff_pars;

      buff_ui                      = residual->ui_r;
      residual->ui_r               = target->ui_r;
      target->ui_r                 = buff_ui;

      buff_p_pars                  = residual->p_pars_r;
      residual->p_pars_r           = target->p_pars_r;
      target->p_pars_r             = buff_p_pars;

      buff_p_lk_loc                = residual->p_lk_loc_rght;
      residual->p_lk_loc_rght      = target->p_lk_loc_rght;
      target->p_lk_loc_rght        = buff_p_lk_loc;

      buff_patt_id                 = residual->patt_id_rght;
      residual->patt_id_rght       = target->patt_id_rght;
      target->patt_id_rght         = buff_patt_id;
    }
  else
    {
      v1                           = target->rght;
      v2                           = target->left;

      buff_p_lk                    = residual->p_lk_rght;
      residual->p_lk_rght          = target->p_lk_left;
      target->p_lk_left            = buff_p_lk;

      buff_p_lk_tip                = residual->p_lk_tip_r;
      residual->p_lk_tip_r         = target->p_lk_tip_l;
      target->p_lk_tip_l           = buff_p_lk_tip;

      buff_scale                   = residual->sum_scale_rght;
      residual->sum_scale_rght     = target->sum_scale_left;
      target->sum_scale_left       = buff_scale;

      buff_scale                   = residual->sum_scale_rght_cat;
      residual->sum_scale_rght_cat = target->sum_scale_left_cat;
      target->sum_scale_left_cat   = buff_scale;

      buff_pars                    = residual->pars_r;
      residual->pars_r             = target->pars_l;
      target->pars_l               = buff_pars;

      buff_ui                      = residual->ui_r;
      residual->ui_r               = target->ui_l;
      target->ui_l                 = buff_ui;

      buff_p_pars                  = residual->p_pars_r;
      residual->p_pars_r           = target->p_pars_l;
      target->p_pars_l             = buff_p_pars;

      buff_p_lk_loc                = residual->p_lk_loc_rght;
      residual->p_lk_loc_rght      = target->p_lk_loc_left;
      target->p_lk_loc_left        = buff_p_lk_loc;

      buff_patt_id                 = residual->patt_id_rght;
      residual->patt_id_rght       = target->patt_id_left;
      target->patt_id_left         = buff_patt_id;
    }

  For(i,3)
    if(v2->b[i] == target)
      {
	v2->v[i] = link;
	v2->b[i] = residual;
	break;
      }

  link->v[dir_v2] = v2;
  link->b[dir_v2] = residual;

  residual->left  = link;
  residual->rght  = v2;

  (v1 == target->left)?(target->rght = link):(target->left = link);

  link->v[dir_v1] = v1;
  link->b[dir_v1] = target;

  For(i,3)
    if(v1->v[i] == v2)
      {
	v1->v[i] = link;
	break;
      }

  target->l /= 2.;
  residual->l = target->l;

  Make_Edge_Dirs(target,target->left,target->rght,tree);
  Make_Edge_Dirs(residual,residual->left,residual->rght,tree);
  Make_Edge_Dirs(b_up,b_up->left,b_up->rght,tree);
}

/*********************************************************/

void Reassign_Node_Nums(t_node *a, t_node *d, int *curr_ext_node, int *curr_int_node, t_tree *tree)
{
  t_node *buff;
  int i;

  if(a->tax)
    {
      buff = tree->noeud[*curr_ext_node];
      tree->noeud[*curr_ext_node] = a;
      tree->noeud[a->num] = buff;
      buff->num = a->num;
      a->num = *curr_ext_node;
      (*curr_ext_node)++;
    }

  if(d->tax)
    {
      buff = tree->noeud[*curr_ext_node];
      tree->noeud[*curr_ext_node] = d;
      tree->noeud[d->num] = buff;
      buff->num = d->num;
      d->num = *curr_ext_node;
      (*curr_ext_node)++;
      return;
    }
  else
    {
      buff = tree->noeud[*curr_int_node];
      tree->noeud[*curr_int_node] = d;
      tree->noeud[d->num] = buff;
      buff->num = d->num;
      d->num = *curr_int_node;
      (*curr_int_node)++;
    }

  For(i,3)
    {
      if(d->v[i] != a)
	Reassign_Node_Nums(d,d->v[i],curr_ext_node,curr_int_node,tree);
    }
}

/*********************************************************/

void Reassign_Edge_Nums(t_node *a, t_node *d, int *curr_br, t_tree *tree)
{
  t_edge *buff;
  int i,j;

  For(i,3)
    if(a->v[i] == d)
      {
	buff = tree->t_edges[*curr_br];
	For(j,2*N_MAX_OTU-3) if(tree->t_edges[j] == a->b[i]) break;
	if(j == 2*N_MAX_OTU-3)
	  {
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	  }
	tree->t_edges[*curr_br] = a->b[i];
	tree->t_edges[j] = buff;
	a->b[i]->num = *curr_br;
	(*curr_br)++;
	break;
      }

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  Reassign_Edge_Nums(d,d->v[i],curr_br,tree);
    }
}

/*********************************************************/

void Find_Mutual_Direction(t_node *n1, t_node *n2, short int *dir_n1_to_n2, short int *dir_n2_to_n1)
{
  int scores[3][3];
  int i,j,k,l;

  if(n1 == n2) return;


  For(i,3)
    {
      For(j,3)
	{
	  scores[i][j] = 0;

	  For(k,n1->bip_size[i])
	    {
	      For(l,n2->bip_size[j])
		{
		  if(n1->bip_node[i][k] == n2->bip_node[j][l])
		    {
		      scores[i][j]++;
		      break;
		    }
		}
	    }
	}
    }

  For(i,3)
    {
      For(j,3)
	{
	  if(!scores[i][j]) 
	    {
	      *dir_n1_to_n2 = i; 
	      *dir_n2_to_n1 = j; 
	      return;
	    } 
	}
    }

  PhyML_Printf("\n. n1=%d n2=%d",n1->num,n2->num);
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  

/*   For(i,3) */
/*     { */
/*       n_zero_line = 0; */
/*       For(j,3) */
/* 	{ */
/* 	  if(!scores[i][j]) n_zero_line++; */
/* 	} */
/*       if(n_zero_line != 2) {*dir_n1_to_n2 = i; break;} */
/*     } */


/*   For(i,3) */
/*     { */
/*       n_zero_col = 0; */
/*       For(j,3) */
/* 	{ */
/* 	  if(!scores[j][i]) n_zero_col++; */
/* 	} */
/*       if(n_zero_col != 2) {*dir_n2_to_n1 = i; break;} */
/*     } */

}

/*********************************************************/

void Update_Dir_To_Tips(t_node *a, t_node *d, t_tree *tree)
{
  int i,j,k;
  short int *inout;
  int d_a;
  int dim;

  dim = 2*tree->n_otu-2;

  inout = (short int *)mCalloc(tree->n_otu,sizeof(short int));

  For(i,3)
    {
      if(a->v[i] == d)
	{
	  For(j,tree->n_otu) inout[j] = 1;
	  For(k,a->bip_size[i]) inout[a->bip_node[i][k]->num] = 0;
	  For(j,tree->n_otu) if(inout[tree->noeud[j]->num]) tree->t_dir[a->num*dim+tree->noeud[j]->num] = i;
	  break;
	}
    }


  if(!d->tax)
    {

      d_a = -1;

      For(i,3)
	{
	  if(d->v[i] != a) Update_Dir_To_Tips(d,d->v[i],tree);
	  else if(d->v[i] == a) d_a = i;
	}

      For(j,tree->n_otu) inout[j] = 1;
      For(k,d->bip_size[d_a]) inout[d->bip_node[d_a][k]->num] = 0;
      For(j,tree->n_otu) if(inout[tree->noeud[j]->num]) tree->t_dir[d->num*dim+tree->noeud[j]->num] = d_a;
    }

  Free(inout);

}

/*********************************************************/

void Fill_Dir_Table(t_tree *tree)
{
  int i,j;
  int dim;
  dim = 2*tree->n_otu-2;
  For(i,dim*dim) tree->t_dir[i] = 0;
  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);
  Update_Dir_To_Tips(tree->noeud[0],tree->noeud[0]->v[0],tree);


  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    for(j=i;j<2*tree->n_otu-2;j++)
      {
	Find_Mutual_Direction(tree->noeud[i],tree->noeud[j],
			      &(tree->t_dir[i*dim+j]),
			      &(tree->t_dir[j*dim+i]));
      }
}

/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/

int Get_Subtree_Size(t_node *a, t_node *d)
{
  int size,i;

  if(d->tax) return 1;
  else
    {
      size = 0;
      For(i,3)
	if(d->v[i] != a)
	  size += Get_Subtree_Size(d,d->v[i]);
    }
  return size;
}

/*********************************************************/

void Fast_Br_Len(t_edge *b, t_tree *tree, int approx)
{
  phydbl sum;
  phydbl *prob, *F;
  int i, j, k, site;
  phydbl v_rght;
  int dim1,dim2,dim3;
  phydbl eps_bl,old_l,new_l;
  int n_iter;
  phydbl scale_rght;

/*   Br_Len_Brent(0.02*b->l,b->l,50.*b->l, */
/* 	       tree->mod->s_opt->min_diff_lk_local, */
/* 	       b,tree, */
/* 	       tree->mod->s_opt->brent_it_max, */
/* 	       tree->mod->s_opt->quickdirty); */
/*   return; */

  n_iter = 0;
  dim1   = tree->mod->ns * tree->mod->n_catg;
  dim2   = tree->mod->ns ;
  dim3   = tree->mod->ns * tree->mod->ns;
  eps_bl = tree->mod->l_min;

  F    = tree->triplet_struct->F_bc;
  prob = tree->triplet_struct->F_cd;

  Update_PMat_At_Given_Edge(b,tree);
  
  For(i,dim1*dim2) F[i] = .0;
  
  For(site,tree->n_pattern)
    {
      /* Joint probabilities of the states at the two ends of the t_edge */
      v_rght = -1.;
      For(i,tree->mod->ns)
	{
	  For(j,tree->mod->ns)
	    {
	      For(k,tree->mod->n_catg)
		{
		  v_rght = (b->rght->tax)?((phydbl)(b->p_lk_tip_r[site*dim2+j])):(b->p_lk_rght[site*dim1+k*dim2+j]);
		  scale_rght = (b->rght->tax)?(0.0):(b->sum_scale_rght[k*tree->n_pattern+site]);

		  prob[dim3*k+dim2*i+j]              =
		    tree->mod->gamma_r_proba[k]      *
		    tree->mod->pi[i]                 *
		    b->Pij_rr[k*dim3+i*dim2+j]       *
		    b->p_lk_left[site*dim1+k*dim2+i] *
		    v_rght *
		    POW(2.,-(b->sum_scale_left[k*tree->n_pattern+site] + scale_rght));
		}
	    }
	}
      
      /* Scaling */
      sum = .0;
      For(k,tree->mod->n_catg) For(i,tree->mod->ns) For(j,tree->mod->ns) sum += prob[dim3*k+dim2*i+j];
      For(k,tree->mod->n_catg) For(i,tree->mod->ns) For(j,tree->mod->ns) prob[dim3*k+dim2*i+j] /= sum;
      
      /* Expected number of each pair of states */
      For(i,tree->mod->ns) For(j,tree->mod->ns) For(k,tree->mod->n_catg)
	F[dim3*k+dim2*i+j] += tree->data->wght[site] * prob[dim3*k+dim2*i+j];
    }
  
  old_l = b->l;
  Opt_Dist_F(&(b->l),F,tree->mod);
  new_l = b->l;
  n_iter++;
  
  if(!approx)
    Br_Len_Brent(0.02*b->l,b->l,50.*b->l,
		 tree->mod->s_opt->min_diff_lk_local,
		 b,tree,
		 tree->mod->s_opt->brent_it_max,
		 tree->mod->s_opt->quickdirty);
  else
    Lk_At_Given_Edge(b,tree);

}

/*********************************************************/

eigen *Make_Eigen_Struct(int ns)
{
  eigen *eig;

  eig              = (eigen *)mCalloc(1,sizeof(eigen));
  eig->size        = ns;
  eig->space       = (phydbl *)mCalloc(2*ns,sizeof(phydbl));
  eig->space_int   = (int *)mCalloc(2*ns,sizeof(int));
  eig->e_val       = (phydbl *)mCalloc(ns,sizeof(phydbl));
  eig->e_val_im    = (phydbl *)mCalloc(ns,sizeof(phydbl));
  eig->r_e_vect    = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->r_e_vect_im = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->l_e_vect    = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));
  eig->q           = (phydbl *)mCalloc(ns*ns,sizeof(phydbl));

  return eig;
}

/*********************************************************/

triplet *Make_Triplet_Struct(model *mod)
{
  int i,j,k;
  triplet *triplet_struct;

  triplet_struct                  = (triplet *)mCalloc(1,sizeof(triplet));
  triplet_struct->size            = mod->ns;
  triplet_struct->pi_bc           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_cd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->pi_bd           = (phydbl *)mCalloc(mod->ns,sizeof(phydbl ));
  triplet_struct->F_bc            = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_catg,sizeof(phydbl));
  triplet_struct->F_cd            = (phydbl *)mCalloc(mod->ns*mod->ns*mod->n_catg,sizeof(phydbl));
  triplet_struct->F_bd            = (phydbl *)mCalloc(mod->ns*mod->ns,sizeof(phydbl));
  triplet_struct->core            = (phydbl ****)mCalloc(mod->n_catg,sizeof(phydbl ***));
  triplet_struct->p_one_site      = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  triplet_struct->sum_p_one_site  = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
  triplet_struct->eigen_struct    = (eigen *)Make_Eigen_Struct(mod->ns);
  triplet_struct->mod             = mod;

  For(k,mod->n_catg)
    {
      triplet_struct->core[k]                = (phydbl ***)mCalloc(mod->ns,sizeof(phydbl **));
      For(i,mod->ns)
	{
	  triplet_struct->core[k][i]         = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
	  For(j,mod->ns)
	    triplet_struct->core[k][i][j]    = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
	}
    }

  For(i,mod->ns)
    {
      triplet_struct->p_one_site[i]          = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
	triplet_struct->p_one_site[i][j]     = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }

  For(i,mod->ns)
    {
      triplet_struct->sum_p_one_site[i]      = (phydbl **)mCalloc(mod->ns,sizeof(phydbl *));
      For(j,mod->ns)
	triplet_struct->sum_p_one_site[i][j] = (phydbl  *)mCalloc(mod->ns,sizeof(phydbl ));
    }
  return triplet_struct;

}

/*********************************************************/

phydbl Triple_Dist(t_node *a, t_tree *tree, int approx)
{
  if(a->tax) return UNLIKELY;
  else
    {
      Update_PMat_At_Given_Edge(a->b[1],tree);
      Update_PMat_At_Given_Edge(a->b[2],tree);

      Update_P_Lk(tree,a->b[0],a);
      Fast_Br_Len(a->b[0],tree,approx);


      Update_P_Lk(tree,a->b[1],a);
      Fast_Br_Len(a->b[1],tree,approx);


      Update_P_Lk(tree,a->b[2],a);
      Fast_Br_Len(a->b[2],tree,approx);

      Update_P_Lk(tree,a->b[1],a);
      Update_P_Lk(tree,a->b[0],a);
    }

  return tree->c_lnL;

}


/*********************************************************/

void Make_Symmetric(phydbl **F, int size)
{
  int i,j;

  For(i,size)
    {
      for(j=i+1;j<size;j++)
	{
	  (*F)[size*i+j] = ((*F)[size*i+j] + (*F)[size*j+i])/2.;
	  (*F)[size*j+i] = (*F)[size*i+j];
	}
    }
}

/*********************************************************/

void Round_Down_Freq_Patt(phydbl **F, t_tree *tree)
{
  int i,j;

  For(i,tree->mod->ns)
    {
      For(j,tree->mod->ns)
	{
	  (*F)[tree->mod->ns*i+j] = RINT((*F)[tree->mod->ns*i+j]);
	}
    }
}

/*********************************************************/

phydbl Get_Sum_Of_Cells(phydbl *F, t_tree *tree)
{
  int i,j;
  phydbl sum = .0;

  For(i,tree->mod->ns)
    For(j,tree->mod->ns)
    sum += F[tree->mod->ns*i+j];

  return sum;
}


/*********************************************************/

void Divide_Cells(phydbl **F, phydbl div, t_tree *tree)
{
  int i,j;

  For(i,tree->mod->ns)
    For(j,tree->mod->ns)
    (*F)[tree->mod->ns*i+j] /= div;
}

/*********************************************************/

void Divide_Mat_By_Vect(phydbl **F, phydbl *vect, int size)
{
  int i,j;
  For(i,size)
    For(j,size)
    (*F)[size*i+j] = (*F)[size*i+j] / vect[j];
}

/*********************************************************/

void Multiply_Mat_By_Vect(phydbl **F, phydbl *vect, int size)
{
  int i,j;
  For(i,size)
    For(j,size)
    (*F)[size*i+j] = (*F)[size*i+j] * vect[j];
}

/*********************************************************/

void Found_In_Subtree(t_node *a, t_node *d, t_node *target, int *match, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      if(d == target) *match =  1;
      For(i,3)
	{
	  if(d->v[i] != a)
	    Found_In_Subtree(d,d->v[i],target,match,tree);
	}
    }
}

/*********************************************************/

void Get_List_Of_Target_Edges(t_node *a, t_node *d, t_edge **list, int *list_size, t_tree *tree)
{
  int i;

  For(i,3)
    {
      if(a->v[i] && a->v[i] == d)
	{
	  list[*list_size] = a->b[i];
	  (*list_size)++;
	}
    }

  if(d->tax) return;
  else
    {
      For(i,3)
	{
	  if(d->v[i] != a)
	    Get_List_Of_Target_Edges(d,d->v[i],list,list_size,tree);
	}
    }
}

/*********************************************************/

void Fix_All(t_tree *tree)
{
  int i;

  tree->mod->pinvar_old = tree->mod->pinvar;
  tree->mod->alpha_old  = tree->mod->alpha;
  tree->mod->kappa_old  = tree->mod->kappa;
  tree->mod->lambda_old = tree->mod->lambda;

  for(i=tree->n_otu;i<2*tree->n_otu-2;i++)
    {
      tree->noeud[i]->b[0]->l_old = tree->noeud[i]->b[0]->l;
      tree->noeud[i]->b[1]->l_old = tree->noeud[i]->b[1]->l;
      tree->noeud[i]->b[2]->l_old = tree->noeud[i]->b[2]->l;
    }
}

/*********************************************************/

void Record_Br_Len(t_tree *tree)
{
  int i;
  
  For(i,2*tree->n_otu-3) 
    {
      tree->t_edges[i]->l_old                = tree->t_edges[i]->l;
      tree->t_edges[i]->gamma_prior_mean_old = tree->t_edges[i]->gamma_prior_mean;
      tree->t_edges[i]->gamma_prior_var_old  = tree->t_edges[i]->gamma_prior_var;
    }
}

/*********************************************************/

void Restore_Br_Len(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-3) 
    {
      tree->t_edges[i]->l                = tree->t_edges[i]->l_old;
      tree->t_edges[i]->gamma_prior_mean = tree->t_edges[i]->gamma_prior_mean_old;
      tree->t_edges[i]->gamma_prior_var  = tree->t_edges[i]->gamma_prior_var_old;
    }
}

/*********************************************************/

void Get_Dist_Btw_Edges(t_node *a, t_node *d, t_tree *tree)
{
  int i;
  t_edge *b_fcus;

  b_fcus = NULL;
  For(i,3) if(a->v[i] == d) {b_fcus = a->b[i]; break;}

  if(d->tax) return;
  else
    {
      For(i,3)
	if(d->v[i] != a)
	  {
	    d->b[i]->topo_dist_btw_edges = b_fcus->topo_dist_btw_edges + 1;
	    d->b[i]->dist_btw_edges      = b_fcus->dist_btw_edges + d->b[i]->l / 2.;
	    Get_Dist_Btw_Edges(d,d->v[i],tree);
	  }
    }


}

/*********************************************************/

void Detect_Polytomies(t_edge *b, phydbl l_thresh, t_tree *tree)
{
  if((b->l < l_thresh) && (!b->left->tax) && (!b->rght->tax))
    {
      b->l               = 0.0;
      b->has_zero_br_len = YES;
    }
  else b->has_zero_br_len = NO;
}

/*********************************************************/

void Get_List_Of_Nodes_In_Polytomy(t_node *a, t_node *d, t_node ***list, int *size_list)
{
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	{
	  if(d->v[i] != a)
	    {
	      if(!d->b[i]->has_zero_br_len)
		{
		  (*list)[*size_list] = d->v[i];
		  (*size_list)++;
		}

	      if(d->b[i]->has_zero_br_len)
		Get_List_Of_Nodes_In_Polytomy(d,d->v[i],list,size_list);
	    }
	}
    }

}


/*********************************************************/

void Check_Path(t_node *a, t_node *d, t_node *target, t_tree *tree)
{
  PhyML_Printf("path---------\n");
  if(d==target) return;
  else Check_Path(d,d->v[tree->t_dir[d->num*(2*tree->n_otu-2)+target->num]],target,tree);
}


/*********************************************************/

void Connect_Two_Nodes(t_node *a, t_node *d)
{
  a->v[0] = d;
  d->v[0] = a;
}

/*********************************************************/

void Get_List_Of_Adjacent_Targets(t_node *a, t_node *d, t_node ***node_list, t_edge ***edge_list, int *list_size)
{
  int i;

  For(i,3)
    if(a->v[i] == d)
      {
	(*node_list)[*list_size] = a;
	(*edge_list)[*list_size] = a->b[i];
	(*list_size)++;
      }
  if(d->tax) return;
  else
    For(i,3)
      if(d->v[i] != a) Get_List_Of_Adjacent_Targets(d,d->v[i],node_list,edge_list,list_size);
}

/*********************************************************/

void Sort_List_Of_Adjacent_Targets(t_edge ***list, int list_size)
{
  t_edge *buff_edge;
  int i,j;

  buff_edge = NULL;

  For(i,list_size-1)
    {
      for(j=i+1;j<list_size;j++)
	if((*list)[j]->topo_dist_btw_edges < (*list)[i]->topo_dist_btw_edges)
	  {
	    buff_edge = (*list)[j];
	    (*list)[j] = (*list)[i];
	    (*list)[i] = buff_edge;
	  }
    }
}

/*********************************************************/


t_node *Common_Nodes_Btw_Two_Edges(t_edge *a, t_edge *b)
{
  if(a->left == b->left)      return b->left;
  else if(a->left == b->rght) return b->rght;
  else if(a->rght == b->left) return b->left;
  else if(a->rght == b->rght) return b->rght;

  PhyML_Printf("\n. First t_edge = %d (%d %d); Second t_edge = %d (%d %d)\n",
	 a->num,a->left->num,a->rght->num,
	 b->num,b->left->num,b->rght->num);
  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
  Warn_And_Exit("");

  return NULL;
}

/*********************************************************/

int KH_Test(phydbl *site_lk_M1, phydbl *site_lk_M2, t_tree *tree)
{
  phydbl *delta,mean,sd,obs_stat,threshold;
  int i;

  
  delta = (phydbl *)mCalloc(tree->data->init_len,sizeof(phydbl));

  threshold = .0;
  mean = .0;
  obs_stat = .0;
  For(i,tree->n_pattern)
    {
      delta[i] = site_lk_M1[i] - site_lk_M2[i];
      mean += ((int)tree->data->wght[i])*delta[i];
    }

  obs_stat = mean;

  mean /= tree->data->init_len;

  For(i,tree->data->init_len) delta[i] -= mean;

  sd = .0;
  For(i,tree->data->init_len) sd += POW(delta[i],2);
  sd /= (phydbl)(tree->data->init_len-1.);

/*   threshold = tree->dnorm_thresh*SQRT(sd*tree->data->init_len); */


/*   PhyML_Printf("\nObs stat = %f Threshold = %f\n",obs_stat,threshold); */
  Free(delta);

  if(obs_stat > threshold) return 1;
  else                     return 0;
}

/*********************************************************/

/*********************************************************/

void Random_Tree(t_tree *tree)
{
  int *is_available,*list_of_nodes;
  int i,node_num,step,n_available;
  phydbl min_edge_len;

  min_edge_len = 1.E-3;

  PhyML_Printf("\n. Randomising the tree...\n");

  is_available  = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  list_of_nodes = (int *)mCalloc(tree->n_otu,    sizeof(int));

  For(i,tree->n_otu) is_available[i]  = 1;
  For(i,tree->n_otu) list_of_nodes[i] = i;

  step = 0;
  do
    {
/*       node_num = (int)RINT(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-1-step)); */
      node_num = Rand_Int(0,tree->n_otu-1-step);
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
      tree->noeud[tree->n_otu+step]->v[1] = tree->noeud[node_num];

/*       node_num = (int)RINT(rand()/(phydbl)(RAND_MAX+1.0)*(tree->n_otu-2-step)); */
      node_num = Rand_Int(0,tree->n_otu-2-step);
      node_num = list_of_nodes[node_num];
      is_available[node_num] = 0;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) {list_of_nodes[n_available++] = i;}

      tree->noeud[node_num]->v[0] = tree->noeud[tree->n_otu+step];
      tree->noeud[tree->n_otu+step]->v[2] = tree->noeud[node_num];

      is_available[tree->n_otu+step] = 1;
      For(i,tree->n_otu) list_of_nodes[i] = -1;
      n_available = 0;
      For(i,2*tree->n_otu-2) if(is_available[i]) list_of_nodes[n_available++] = i;

      step++;
    }while(step < tree->n_otu-2);

  tree->noeud[list_of_nodes[0]]->v[0] = tree->noeud[list_of_nodes[1]];
  tree->noeud[list_of_nodes[1]]->v[0] = tree->noeud[list_of_nodes[0]];

  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
  
  For(i,2*tree->n_otu-3) if(tree->t_edges[i]->l < min_edge_len) tree->t_edges[i]->l = min_edge_len;

  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  
  Free(is_available);
  Free(list_of_nodes);
}

/*********************************************************/


void Random_NNI(int n_moves, t_tree *tree)
{
  int i,j;
  t_edge *b;
  t_node *n1,*n2,*n_target;

  n1 = n2 = NULL;
  b = NULL;
  For(i,n_moves)
    {
      n_target  = tree->noeud[tree->n_otu + (int)((phydbl)rand()/RAND_MAX * (2*tree->n_otu-3-tree->n_otu))];
      For(j,3) if(!n_target->v[j]->tax) {b = n_target->b[j]; break;}


      For(j,3) if(b->left->v[j] != b->rght) {n1 = b->left->v[j]; break;}
      For(j,3) if(b->rght->v[j] != b->left) {n2 = b->rght->v[j]; break;}


      Swap(n1,b->left,b->rght,n2,tree);
    }
}

/*********************************************************/


void Print_Settings(option *io)
{
  int answer;
  char *s;

  s = (char *)mCalloc(100,sizeof(char));
  
  PhyML_Printf("\n\n\n");
  PhyML_Printf("\n\n");

  PhyML_Printf("                                 ..........................                                      \n");
  PhyML_Printf(" ooooooooooooooooooooooooooooo        CURRENT SETTINGS        ooooooooooooooooooooooooooooooooooo\n");
  PhyML_Printf("                                 ..........................                                      \n");

  PhyML_Printf("\n                . Sequence filename:\t\t\t\t %s", Basename(io->in_align_file));

  if(io->datatype == NT) strcpy(s,"dna");
  else if(io->datatype == AA) strcpy(s,"aa");
  else strcpy(s,"generic");

  PhyML_Printf("\n                . Data type:\t\t\t\t\t %s",s);
  PhyML_Printf("\n                . Alphabet size:\t\t\t\t %d",io->mod->ns);

  PhyML_Printf("\n                . Sequence format:\t\t\t\t %s", io->interleaved ? "interleaved": "sequential");
  PhyML_Printf("\n                . Number of data sets:\t\t\t\t %d", io->n_data_sets);

  PhyML_Printf("\n                . Nb of bootstrapped data sets:\t\t\t %d", io->mod->bootstrap);

  if (io->mod->bootstrap > 0)
    PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t no");
  else
    {
      if(io->ratio_test == 1)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (aLRT statistics)");
      else if(io->ratio_test == 2)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (Chi2-based parametric branch supports)");
      else if(io->ratio_test == 3)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (Minimum of SH-like and Chi2-based branch supports)");
      else if(io->ratio_test == 4)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (SH-like branch supports)");
      else if(io->ratio_test == 5)
	PhyML_Printf("\n                . Compute approximate likelihood ratio test:\t yes (aBayes branch supports)");
    }

  PhyML_Printf("\n                . Model name:\t\t\t\t\t %s", io->mod->modelname);

  if (io->datatype == NT)
    {
      if ((io->mod->whichmodel == K80)  ||
	  (io->mod->whichmodel == HKY85)||
	  (io->mod->whichmodel == F84)  ||
	  (io->mod->whichmodel == TN93))
	{
	  if (io->mod->s_opt->opt_kappa)
	    PhyML_Printf("\n                . Ts/tv ratio:\t\t\t\t\t estimated");
	  else
	    PhyML_Printf("\n                . Ts/tv ratio:\t\t\t\t\t %f", io->mod->kappa);
	}
    }

  if (io->mod->s_opt->opt_pinvar)
    PhyML_Printf("\n                . Proportion of invariable sites:\t\t estimated");
  else
    PhyML_Printf("\n                . Proportion of invariable sites:\t\t %f", io->mod->pinvar);


  PhyML_Printf("\n                . Number of subst. rate categs:\t\t\t %d", io->mod->n_catg);
  if(io->mod->n_catg > 1)
    {
      if(io->mod->free_mixt_rates == NO)
	{
	  if(io->mod->s_opt->opt_alpha)
	    PhyML_Printf("\n                . Gamma distribution parameter:\t\t\t estimated");
	  else
	    PhyML_Printf("\n                . Gamma distribution parameter:\t\t\t %f", io->mod->alpha);
	  PhyML_Printf("\n                . 'Middle' of each rate class:\t\t\t %s",(io->mod->gamma_median)?("median"):("mean"));
	}
    }
    
  
  if(io->datatype == AA)
    PhyML_Printf("\n                . Amino acid equilibrium frequencies:\t\t %s", (io->mod->s_opt->opt_state_freq) ? ("empirical"):("model"));
  else if(io->datatype == NT)
    {
      if((io->mod->whichmodel != JC69) &&
	 (io->mod->whichmodel != K80)  &&
	 (io->mod->whichmodel != F81))
	{
	  if(!io->mod->s_opt->user_state_freq)
	    {
	      PhyML_Printf("\n                . Nucleotide equilibrium frequencies:\t\t %s", (io->mod->s_opt->opt_state_freq) ? ("ML"):("empirical"));
	    }
	  else
	    {
	      PhyML_Printf("\n                . Nucleotide equilibrium frequencies:\t\t %s","user-defined");
	    }
	}
    }

  PhyML_Printf("\n                . Optimise tree topology:\t\t\t %s", (io->mod->s_opt->opt_topo) ? "yes": "no");

  switch(io->in_tree)
    {
    case 0: { strcpy(s,"BioNJ");     break; }
    case 1: { strcpy(s,"parsimony"); break; }
    case 2: { strcpy(s,"user tree ("); 
	strcat(s,Basename(io->in_tree_file)); 
	strcat(s,")");         break; }
    }

  if(io->mod->s_opt->opt_topo)
    {
      if(io->mod->s_opt->topo_search == NNI_MOVE) PhyML_Printf("\n                . Tree topology search:\t\t\t\t NNIs");
      else if(io->mod->s_opt->topo_search == SPR_MOVE) PhyML_Printf("\n                . Tree topology search:\t\t\t\t SPRs");
      else if(io->mod->s_opt->topo_search == BEST_OF_NNI_AND_SPR) PhyML_Printf("\n                . Tree topology search:\t\t\t\t Best of NNIs and SPRs");



      PhyML_Printf("\n                . Starting tree:\t\t\t\t %s",s);

      PhyML_Printf("\n                . Add random input tree:\t\t\t %s", (io->mod->s_opt->random_input_tree) ? "yes": "no");
      if(io->mod->s_opt->random_input_tree)
	PhyML_Printf("\n                . Number of random starting trees:\t\t %d", io->mod->s_opt->n_rand_starts);	
    }
  else
    if(!io->mod->s_opt->random_input_tree)
      PhyML_Printf("\n                . Evaluted tree:\t\t\t\t file \"%s\"",s);

  PhyML_Printf("\n                . Optimise branch lengths:\t\t\t %s", (io->mod->s_opt->opt_bl) ? "yes": "no");

  answer = 0;
  if(io->mod->s_opt->opt_alpha  ||
     io->mod->s_opt->opt_kappa  ||
     io->mod->s_opt->opt_lambda ||
     io->mod->s_opt->opt_pinvar ||
     io->mod->s_opt->opt_rr) answer = 1;
  
  PhyML_Printf("\n                . Optimise substitution model parameters:\t %s", (answer) ? "yes": "no");

  PhyML_Printf("\n                . Run ID:\t\t\t\t\t %s", (io->appebr_run_ID) ? (io->run_id_string): ("none"));
  PhyML_Printf("\n                . Random seed:\t\t\t\t\t %d", io->r_seed);
  PhyML_Printf("\n                . Subtree patterns aliasing:\t\t\t %s",io->do_alias_subpatt?"yes":"no");
  PhyML_Printf("\n                . Version:\t\t\t\t\t %s", VERSION);


  PhyML_Printf("\n\n oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");

  PhyML_Printf("\n\n");
  fflush(NULL);
  
  Free(s);
}

/*********************************************************/

void Fill_Missing_Dist(matrix *mat)
{
  int i,j;
  For(i,mat->n_otu)
    {
      for(j=i+1;j<mat->n_otu;j++)
	{
	  if(i != j)
	    {
	      if(mat->dist[i][j] < .0) 
		{
		  Fill_Missing_Dist_XY(i,j,mat);
		  mat->dist[j][i] = mat->dist[i][j];
		}
	    }
	}
    }
}

/*********************************************************/

void Fill_Missing_Dist_XY(int x, int y, matrix *mat)
{

  int i,j;
  phydbl *local_mins,**S1S2;
  int cpt;
  int pos_best_estimate;
  phydbl min_crit, curr_crit;

  local_mins = (phydbl *)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl ));
  S1S2       = (phydbl **)mCalloc(mat->n_otu*mat->n_otu,sizeof(phydbl *));
  For(i,mat->n_otu*mat->n_otu) S1S2[i] = (phydbl *)mCalloc(2,sizeof(phydbl));

  cpt = 0;
  For(i,mat->n_otu)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
	{
	  For(j,mat->n_otu)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
		{
		  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
		    {
		      S1S2[cpt][0] = MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
		      S1S2[cpt][1] = MAX(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]);
		      cpt++;
		    }
		}
	    }
	}
    }

  Qksort_Matrix(S1S2,0,0,cpt-1);

  local_mins[0] = S1S2[0][1];
  for(i=1;i<cpt;i++) local_mins[i] = (i*local_mins[i-1] + S1S2[i][1])/(phydbl)(i+1);
 
  pos_best_estimate = 0;
  min_crit = curr_crit = BIG;
	
  For(i,cpt-1)
    {
      if((local_mins[i] < S1S2[i+1][0]) && (local_mins[i] > S1S2[i][0]))
	{
	  curr_crit = Least_Square_Missing_Dist_XY(x,y,local_mins[i],mat);
	  if(curr_crit < min_crit)
	    {
	      min_crit = curr_crit;
	      pos_best_estimate = i;
	    }
	}
    }

  mat->dist[x][y] = local_mins[pos_best_estimate];
  mat->dist[y][x] = mat->dist[x][y];

  For(i,mat->n_otu*mat->n_otu) Free(S1S2[i]);
  Free(S1S2);
  Free(local_mins);
}

/*********************************************************/

phydbl Least_Square_Missing_Dist_XY(int x, int y, phydbl dxy, matrix *mat)
{
  int i,j;
  phydbl fit;

  fit = .0;
  For(i,mat->n_otu)
    {
      if((mat->dist[i][x] > .0) && (mat->dist[i][y] > .0))
	{
	  For(j,mat->n_otu)
	    {
	      if((mat->dist[j][x] > .0) && (mat->dist[j][y] > .0))
		{
		  if((i != j) && (i != x) && (i != y) && (j != x) && (j != y))
		    {
		      if(dxy < MIN(mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j] , mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]))
			{
			  fit += POW((mat->dist[i][x] + mat->dist[j][y]) - (mat->dist[i][y] + mat->dist[j][x]),2);
			}
		      else if((mat->dist[i][x] + mat->dist[j][y]) < (mat->dist[i][y] + mat->dist[j][x]))
			{
			  fit += POW(dxy - (mat->dist[i][y] + mat->dist[j][x] - mat->dist[i][j]),2);
			}
		      else
			{
			  fit += POW(dxy - (mat->dist[i][x] + mat->dist[j][y] - mat->dist[i][j]),2);
			}
		    }
		}
	    }
	}
    }
  return fit;
}

/*********************************************************/

void Print_Banner(FILE *fp)
{
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                                 ---  PhyML %s  ---                                             \n",VERSION);
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"    A simple, fast, and accurate algorithm to estimate large phylogenies by maximum likelihood    \n");
  PhyML_Fprintf(fp,"                            Stephane Guindon & Olivier Gascuel                                      \n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                           http://www.atgc-montpellier.fr/phyml                                          \n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/

void Print_Banner_Small(FILE *fp)
{
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
  PhyML_Fprintf(fp,"                                  ---  PhyML %s  ---                                             \n",VERSION);
  PhyML_Fprintf(fp,"                            http://www.atgc-montpellier.fr/phyml                                          \n");
  PhyML_Fprintf(fp,"                         Copyright CNRS - Universite Montpellier II                                 \n");
  PhyML_Fprintf(fp," oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/


void Print_Data_Set_Number(option *io, FILE *fp)
{
  PhyML_Fprintf(fp,"\n");
  PhyML_Fprintf(fp,"                                                                                                  \n");
  PhyML_Fprintf(fp,"                                 [ Data set number %3d ]                                           \n",io->curr_gt+1);
  PhyML_Fprintf(fp,"                                                                                                  \n");
}
/*********************************************************/


void Check_Memory_Amount(t_tree *tree)
{
  /* Rough estimate of the amount of memory that has to be used */

  long int nbytes;
  int n_otu;
  model *mod;

  mod    = tree->mod;
  n_otu  = tree->io->n_otu;
  nbytes = 0;

  /* Partial Pars */
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * sizeof(int);
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * sizeof(unsigned int);
  nbytes += (2*n_otu-3) * 2 * tree->data->crunch_len * mod->ns * sizeof(int);

  /* Pmat */
  nbytes += (2*n_otu-3) * mod->n_catg * mod->ns * mod->ns * sizeof(phydbl);
  

  /* Partial Lk */
  nbytes += ((2*n_otu-3) * 2 - tree->n_otu) * tree->data->crunch_len * mod->n_catg * mod->ns * sizeof(phydbl);


  /* Scaling factors */
  nbytes += ((2*n_otu-3) * 2 - tree->n_otu) * tree->data->crunch_len * sizeof(int);
  
      
    
  if(((phydbl)nbytes/(1.E+06)) > 256.)
/*   if(((phydbl)nbytes/(1.E+06)) > 0.) */
    {
      PhyML_Printf("\n. WARNING: this analysis requires at least %.0f MB of memory space.\n",(phydbl)nbytes/(1.E+06));
    }
  else if(((phydbl)nbytes/(1.E+06)) > 100.)
    {
      if(!tree->io->quiet) PhyML_Printf("\n. WARNING: this analysis will use at least %.0f Mo of memory space...\n",(phydbl)nbytes/(1.E+06));
    }
  else if(((phydbl)nbytes/(1.E+06)) > 1.)
    {
      if(!tree->io->quiet) PhyML_Printf("\n. This analysis requires at least %.0f Mo of memory space.\n",(phydbl)nbytes/(1.E+06));
    }
}

/*********************************************************/

int Get_State_From_P_Lk(phydbl *p_lk, int pos, t_tree *tree)
{
  int i;
  For(i,tree->mod->ns) if(p_lk[pos+i] > .0) return i;
  return -1;
}

/*********************************************************/

int Get_State_From_P_Pars(short int *p_pars, int pos, t_tree *tree)
{
  int i;
  For(i,tree->mod->ns) if(p_pars[pos+i] > .0) return i;
  return -1;
}

/*********************************************************/

void Print_Lk(t_tree *tree, char *string)
{
  time(&(tree->t_current));
  PhyML_Printf("\n. (%5d sec) [%15.4f] %s",
	       (int)(tree->t_current-tree->t_beg),tree->c_lnL,
	       string);
#ifndef QUIET 
  fflush(NULL);
#endif
}

/*********************************************************/

void Print_Pars(t_tree *tree)
{
  time(&(tree->t_current));
  PhyML_Printf("\n. (%5d sec) [%5d]",(int)(tree->t_current-tree->t_beg),tree->c_pars);
#ifndef QUIET
  fflush(NULL);
#endif
}

/*********************************************************/

void Print_Lk_And_Pars(t_tree *tree)
{	
  time(&(tree->t_current));

  PhyML_Printf("\n. (%5d sec) [%15.4f] [%5d]",
	 (int)(tree->t_current-tree->t_beg),
	 tree->c_lnL,tree->c_pars);

#ifndef QUIET
  fflush(NULL);
#endif
}

/*********************************************************/

void Check_Dirs(t_tree *tree)
{
  int i;

  For(i,2*tree->n_otu-3)
    {
      if(!tree->t_edges[i]->left->tax)
	{
	  if(tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num <
	     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num)
	    {
	      PhyML_Printf("\n. Edge %d ; v1=%d v2=%d",
		     tree->t_edges[i]->num,
		     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v1]->num,
		     tree->t_edges[i]->left->v[tree->t_edges[i]->l_v2]->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}

      if(!tree->t_edges[i]->rght->tax)
	{
	  if(tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num <
	     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num)
	    {
	      PhyML_Printf("\n. Edge %d ; v3=%d v4=%d",
		     tree->t_edges[i]->num,
		     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v1]->num,
		     tree->t_edges[i]->rght->v[tree->t_edges[i]->r_v2]->num);
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	}
    }
}

/*********************************************************/

void Warn_And_Exit(char *s)
{
  PhyML_Fprintf(stdout,"%s",s);
  fflush(NULL);
#ifndef BATCH
  /* if (! tree->io->quiet) { */
  /* char c; */
  PhyML_Fprintf(stdout,"\n. Type enter to exit.\n");
  /*     if(!fscanf(stdin,"%c",&c))  */
  Exit("");
  /* } */
#endif
  Exit("\n");
}

/*********************************************************/

void Read_Qmat(phydbl *daa, phydbl *pi, FILE *fp)
{
  int i,j;
  phydbl sum;
  double val;

  rewind(fp);

  for(i=1;i<20;i++)
    {
      For(j,19)
	{
/* 	  if(!fscanf(fp,"%lf",&(daa[i*20+j]))) Exit("\n"); */
	  if(!fscanf(fp,"%lf",&val)) 
	    {
	      PhyML_Printf("\n. Rate matrix file does not appear to have a proper format. Please refer to the documentation.");
	      Exit("\n");
	    }
	  daa[i*20+j] = (phydbl)val;
	  daa[j*20+i] = daa[i*20+j];
	  if(j == i-1) break; 
	}
    }


  For(i,20) 
    { 
      if(!fscanf(fp,"%lf",&val)) Exit("\n");
      pi[i] = (phydbl)val;
    }
  sum = .0;
  For(i,20) sum += pi[i];
  if(FABS(sum - 1.) > 1.E-06)
    {
      PhyML_Printf("\n. Sum=%f",sum);
      PhyML_Printf("\n. Scaling amino-acid frequencies...\n");
      For(i,20) pi[i] /= sum;
    }
}

/*********************************************************/

void Print_Qmat_AA(phydbl *daa, phydbl *pi)
{
  int i,j,cpt;

  cpt = 0;
  For(i,20)
    {
      for(j=0;j<i;j++)
	{
	  PhyML_Printf("daa[%2d*20+%2d] = %10f;  ",i,j,daa[i*20+j]);
	  cpt++;
	  if(!(cpt%4)) PhyML_Printf("\n");
	}
    }

  PhyML_Printf("\n\n");
  PhyML_Printf("for (i=0; i<naa; i++)  for (j=0; j<i; j++)  daa[j*naa+i] = daa[i*naa+j];\n\n");
  For(i,20) PhyML_Printf("pi[%d] = %f; ",i,pi[i]);
  PhyML_Printf("\n");
  PhyML_Printf("Ala\tArg\tAsn\tAsp\tCys\tGln\tGlu\tGly\tHis\tIle\tLeu\tLys\tMet\tPhe\tPro\tSer\tThr\tTrp\tTyr\tVal\n");
}


/*********************************************************/

void Randomize_Sequence_Order(calign *cdata)
{
  int i,exchange_with;
  phydbl buff_dbl;
  char *buff_name,*buff_state;
  short int *buff_ambigu;
  
  exchange_with = -1;
  For(i,cdata->n_otu)
    {
      buff_dbl  = rand();
      buff_dbl /= (RAND_MAX+1.);
      buff_dbl *= cdata->n_otu;
      exchange_with = (int)FLOOR(buff_dbl);
      
      buff_name                         = cdata->c_seq[i]->name;
      cdata->c_seq[i]->name             = cdata->c_seq[exchange_with]->name;
      cdata->c_seq[exchange_with]->name = buff_name;

      buff_state                         = cdata->c_seq[i]->state;
      cdata->c_seq[i]->state             = cdata->c_seq[exchange_with]->state;
      cdata->c_seq[exchange_with]->state = buff_state;

      buff_ambigu                            = cdata->c_seq[i]->is_ambigu;
      cdata->c_seq[i]->is_ambigu             = cdata->c_seq[exchange_with]->is_ambigu;
      cdata->c_seq[exchange_with]->is_ambigu = buff_ambigu;
    }
}

/*********************************************************/

void Update_Root_Pos(t_tree *tree)
{
  if(tree->n_root_pos > -1.0)
    {
      tree->n_root->l[0] = tree->e_root->l * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l * (1.-tree->n_root_pos);
    }
  else
    {
/*       tree->n_root->l[0] = tree->e_root->l / 2.; */
/*       tree->n_root->l[1] = tree->e_root->l / 2.; */
    }
}

/*********************************************************/

void Add_Root(t_edge *target, t_tree *tree)
{
  #ifndef PHYML
  PhyML_Printf("\n. Adding root on t_edge %d left = %d right = %d\n.",target->num,target->left->num,target->rght->num); fflush(NULL);
  #endif

  tree->e_root = target;

  /* Create the root t_node if it does not exist yet */
  if((!tree->n_root) || (tree->n_root->num != 2*tree->n_otu-2))
    {      
      tree->n_root = (t_node *)Make_Node_Light(2*tree->n_otu-2);
    }

  tree->noeud[2*tree->n_otu-2] = tree->n_root;

  tree->n_root->tax = 0;

  /* Set the position of the root */
  tree->n_root->v[0] = tree->e_root->left;
  tree->n_root->v[1] = tree->e_root->rght;

  tree->n_root->b[0] = tree->e_root;
  tree->n_root->b[1] = tree->e_root;

  if(tree->n_root_pos > -1.0)
    {
      if(tree->n_root_pos < 1.E-6 &&  tree->n_root_pos > -1.E-6)
	printf("\n. WARNING: you put the root at a weird position...");

/*       tree->n_root->l[0] = tree->e_root->l * (tree->n_root_pos/(1.+tree->n_root_pos)); */
/*       tree->n_root->l[1] = tree->e_root->l - tree->n_root->l[0]; */

      tree->n_root->l[0] = tree->e_root->l * tree->n_root_pos;
      tree->n_root->l[1] = tree->e_root->l * (1. - tree->n_root_pos);
    }
  else
    {
      tree->n_root->l[0] = tree->e_root->l / 2.;
      tree->n_root->l[1] = tree->e_root->l / 2.;
      tree->n_root_pos = 0.5;
    }
  
  Update_Ancestors(tree->n_root,tree->n_root->v[0],tree);
  Update_Ancestors(tree->n_root,tree->n_root->v[1],tree);
  tree->n_root->anc = NULL;
}

/*********************************************************/

void Update_Ancestors(t_node *a, t_node *d, t_tree *tree)
{
  if(a == tree->n_root) a->anc = NULL;
  d->anc = a;
  if(d->tax) return;
  else
    {
      int i;
      For(i,3) 
	if((d->v[i] != d->anc) && (d->b[i] != tree->e_root)) 
	  Update_Ancestors(d,d->v[i],tree);
    }
}

/*********************************************************/
/* Generate a random unrooted tree with 'n_otu' OTUs */
#ifdef PHYTIME
t_tree *Generate_Random_Tree_From_Scratch(int n_otu, int rooted)
{
  t_tree *tree;
  int *connected,*nonconnected,*available_nodes;
  int i,n_connected,n_nonconnected,n1,n2,new_n,n_internal,n_external,n_available;
  t_node *root,*curr_n,**internal_nodes, **external_nodes;
  phydbl *t,*tmp;

  tree = Make_Tree_From_Scratch(n_otu,NULL);

  tree->rates = RATES_Make_Rate_Struct(tree->n_otu);
  RATES_Init_Rate_Struct(tree->rates,tree->io->rates,tree->n_otu);

  For(i,2*tree->n_otu-2) 
    {
      tree->noeud[i]->v[1] = NULL;
      tree->noeud[i]->v[2] = NULL;
    }
  
  root = (t_node *)Make_Node_Light(2*tree->n_otu-2);

  connected       = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  nonconnected    = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  available_nodes = (int *)mCalloc(2*tree->n_otu-2,sizeof(int));
  internal_nodes  = (t_node **)mCalloc(tree->n_otu-2,sizeof(t_node *));
  external_nodes  = (t_node **)mCalloc(tree->n_otu,  sizeof(t_node *));
  t               = (phydbl *)mCalloc(tree->n_otu-1,sizeof(phydbl ));
  tmp             = (phydbl *)mCalloc(2*tree->n_otu-2,sizeof(phydbl ));

  n_nonconnected = 2*n_otu-2;

  For(i,2*tree->n_otu-2) nonconnected[i] = i;

  available_nodes[0] = 2*n_otu-2;

  /* Node times are generated according to a Birth-death process.
     Formulae are as described by Yang and Rannala (1997) */
  phydbl    phi;
  phydbl    rho; /* sampling intensity */
  phydbl     mu; /* birth rate */
  phydbl lambda; /* death rate */
  phydbl      u; /* random U[0,1] */
  phydbl expval;

  /* rho = 1.0 and mu = 0.0 correspond to the Yule process */

  lambda = 6.7;
  mu     = 2.5;
  rho    = 9./150.;

  expval = EXP(mu-lambda);
  phi = (rho*lambda*(expval-1.) + (mu-lambda)*expval)/(expval-1.); /* Equation 16 */

  For(i,tree->n_otu-1)
    {
      u = rand();
      u /= RAND_MAX;

      if(FABS(lambda - mu) > 1.E-4)
	t[i] = (LOG(phi-u*rho*lambda) - LOG(phi-u*rho*lambda + u*(lambda-mu)))/(mu-lambda); /* Equation 15 */
      else
	t[i] = u / (1.+lambda*rho*(1-u)); /* Equation 17 */
    }

  Qksort(t,NULL,0,tree->n_otu-2); /* Node times ordering in ascending order */

  For(i,tree->n_otu-1) tmp[i] =  t[tree->n_otu-2-i];
  For(i,tree->n_otu-1) t[i]   = -tmp[i];
  

  /* Rescale t_node times such that the time at the root t_node is -100 */
  for(i=1;i<tree->n_otu-1;i++) 
    { 
      t[i] /= -t[0]; 
      t[i] *= 1.E+02;
    }
  t[0] = -1.E+02;


  n_available = 1;
  curr_n = root;
  n_connected = 0;
  do
    {
      n1 = Rand_Int(0,n_nonconnected-1);
      n1 = nonconnected[n1];
      connected[n1] = 1;

      n_nonconnected = 0;
      For(i,2*tree->n_otu-2) if(!connected[i]) {nonconnected[n_nonconnected++] = i;}

      n2 = Rand_Int(0,n_nonconnected-1);
      n2 = nonconnected[n2];
      connected[n2] = 1;

      n_nonconnected = 0;
      For(i,2*tree->n_otu-2) if(!connected[i]) {nonconnected[n_nonconnected++] = i;}

      curr_n->v[1] = tree->noeud[n1];
      curr_n->v[2] = tree->noeud[n2];
      tree->noeud[n1]->v[0] = curr_n;
      tree->noeud[n2]->v[0] = curr_n;
      
      tree->rates->nd_t[curr_n->num] = t[n_connected/2];

      available_nodes[n_available] = tree->noeud[n1]->num;  
      For(i,n_available)
	if(available_nodes[i] == curr_n->num) 
	  {
	    available_nodes[i] = tree->noeud[n2]->num;
	    break;
	  }
      n_available++;
      
      new_n = Rand_Int(0,n_available-1);
      curr_n = tree->noeud[available_nodes[new_n]];
      
      n_connected+=2;

    }while(n_connected < 2*tree->n_otu-2);

  For(i,2*tree->n_otu-2) tmp[i] = tree->rates->nd_t[i];

  /* Unroot the tree */
  root->v[1]->v[0] = root->v[2];
  root->v[2]->v[0] = root->v[1];

  n_internal = n_external = 0;
  For(i,2*tree->n_otu-2)
    {
      if(tree->noeud[i]->v[1]) internal_nodes[n_internal++] = tree->noeud[i];
      else                     external_nodes[n_external++] = tree->noeud[i];
    }


  n_internal = n_external = 0;  
  For(i,2*tree->n_otu-2) 
    { 
      if(i < tree->n_otu)
	{
	  tree->noeud[i]      = external_nodes[n_external++];
	  tree->noeud[i]->tax = 1;	  
	}
      else
	{
	  tree->rates->nd_t[i] = tmp[internal_nodes[n_internal]->num];
	  tree->noeud[i]        = internal_nodes[n_internal++];
	  tree->noeud[i]->tax   = 0;
	}

      tree->noeud[i]->num = i;
    }

  For(i,tree->n_otu) tree->rates->nd_t[i] = 0.0;
  
  For(i,tree->n_otu) 
    {
      if(!tree->noeud[i]->name) tree->noeud[i]->name = (char *)mCalloc(T_MAX_NAME,sizeof(char));
      strcpy(tree->noeud[i]->name,"x");
      sprintf(tree->noeud[i]->name+1,"%d",i);
    }


  tree->num_curr_branch_available = 0;
  Connect_Edges_To_Nodes_Recur(tree->noeud[0],tree->noeud[0]->v[0],tree);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);


  /* Add root */
  if(rooted)
    {
      For(i,2*tree->n_otu-3) 
	{
	  if(((tree->t_edges[i]->left == root->v[1]) || (tree->t_edges[i]->rght == root->v[1])) &&
	     ((tree->t_edges[i]->left == root->v[2]) || (tree->t_edges[i]->rght == root->v[2])))
	    {
	      Add_Root(tree->t_edges[i],tree);
	      break;
	    }
	}
    }
  /* Or not... */
  else
    {
      Free_Node(root);
    }

  RATES_Random_Branch_Lengths(tree);
  
  Free(available_nodes);
  Free(connected);
  Free(nonconnected);
  Free(external_nodes);
  Free(internal_nodes);
  Free(t);
  Free(tmp);

  return tree;
}
#endif
/*********************************************************/

void Random_Lineage_Rates(t_node *a, t_node *d, t_edge *b, phydbl stick_prob, phydbl *rates, int curr_rate, int n_rates, t_tree *tree)
{
  phydbl uni;
  int new_rate;
  int i;


  if(b)
    {
      uni  = rand();
      uni /= RAND_MAX;
      
      if(uni > stick_prob) /* Randomly pick a new rate */
	{
	  uni  = rand();
	  uni /= RAND_MAX;
	  uni = (phydbl)(uni * (n_rates-1));	  
	  if(uni-(int)(uni) > 0.5-BIG) new_rate = (int)(uni)+1;
	  else new_rate = (int)(uni);	  
	}
      else
	{
	  new_rate = curr_rate;
	}

      For(i,3) 
	if(a->v[i] == d) 
	  {
	    a->b[i]->l *= rates[new_rate];
	    break;
	  }

      For(i,3)
	if(a->v[i] == d)
	  {
	    if(!(a->b[i]->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(a->b[i]);
	    if(rates[new_rate] > 1.0)      strcpy(a->b[i]->labels[a->b[i]->n_labels],"FAST");
	    else if(rates[new_rate] < 1.0) strcpy(a->b[i]->labels[a->b[i]->n_labels],"SLOW");
	    else                           strcpy(a->b[i]->labels[a->b[i]->n_labels],"MEDIUM");
	    a->b[i]->n_labels++;
	    break;
	  }
      curr_rate = new_rate;
    }
  
  if(d->tax) return;
  else
    {
      For(i,3) 
	if(d->v[i] != a) 
	  Random_Lineage_Rates(d,d->v[i],d->b[i],stick_prob,rates,curr_rate,n_rates,tree);
    }
}

/*********************************************************/

t_edge *Find_Edge_With_Label(char *label, t_tree *tree)
{
  int i,j;

  For(i,2*tree->n_otu-3)
    {
      For(j,tree->t_edges[i]->n_labels)
	{
	  if(!strcmp(tree->t_edges[i]->labels[j],label)) return tree->t_edges[i];
	}
    }
  return NULL;
}

/*********************************************************/

void Print_Square_Matrix_Generic(int n, phydbl *mat)
{
  int i,j;

  PhyML_Printf("\n");
  For(i,n)
    {
      PhyML_Printf("[%3d]",i);
      For(j,n)
	{
	  PhyML_Printf("%7.1G ",mat[i*n+j]);
	}
      PhyML_Printf("\n");
    }
  PhyML_Printf("\n");
}

/*********************************************************/

void Evolve(calign *data, model *mod, t_tree *tree)
{
  int root_state, root_rate_class;
  int site,i;

  data->n_otu = tree->n_otu;

  if(mod->use_m4mod) tree->print_labels = 1;
  
  /* Get the change probability matrices */
  Set_Model_Parameters(mod);
  For(i,2*tree->n_otu-3) Update_PMat_At_Given_Edge(tree->t_edges[i],tree);

  For(site,data->init_len)
    {
      root_state = root_rate_class = -1;

      /* Pick the root nucleotide/aa */
      root_state = Pick_State(mod->ns,mod->pi);
      data->c_seq[0]->state[site] = Reciproc_Assign_State(root_state,tree->io->datatype);

      /* Pick the rate class */
      root_rate_class = Pick_State(mod->n_catg,mod->gamma_r_proba);

      /* tree->noeud[0] is considered as the root t_node */
      Evolve_Recur(tree->noeud[0],
		   tree->noeud[0]->v[0],
		   tree->noeud[0]->b[0],
		   root_state,
		   root_rate_class,
		   site,
		   data,
		   mod,
		   tree);

/*       PhyML_Printf("%s\n",Write_Tree(tree)); */
      
      data->wght[site] = 1;
    }
  data->crunch_len = data->init_len;
}

/*********************************************************/

int Pick_State(int n, phydbl *prob)
{
  int pos;
  phydbl uni;

  do
    {
      pos  = rand();
      pos  = (pos % n);
      uni  = (phydbl)rand();
      uni /= (phydbl)RAND_MAX;
      if(uni < prob[pos]) break;
    }
  while(1);
  
  return (int)pos;
}

/*********************************************************/

void Evolve_Recur(t_node *a, t_node *d, t_edge *b, int a_state, int r_class, int site_num, calign *gen_data, model *mod, t_tree *tree)
{
  int d_state;
  int dim1,dim2;

  dim1 = tree->mod->ns * tree->mod->ns;
  dim2 = tree->mod->ns;

  d_state = Pick_State(mod->ns,b->Pij_rr+r_class*dim1+a_state*dim2);
  
/*   PhyML_Printf("\n>> %c (%d,%d)",Reciproc_Assign_State(d_state,mod->io->datatype),d_state,(int)d_state/mod->m4mod->n_o); */

  if(mod->use_m4mod) 
    {
      phydbl rrate; /* relative rate of substitutions */
      
      rrate = mod->m4mod->multipl[(int)d_state/mod->m4mod->n_o];
      if(!(b->n_labels%BLOCK_LABELS)) Make_New_Edge_Label(b);
      if(rrate > 1.0) strcpy(b->labels[b->n_labels],"FASTER");
      else strcpy(b->labels[b->n_labels],"SLOWER");
      b->n_labels++;
    }

  if(d->tax) 
    {
      gen_data->c_seq[d->num]->state[site_num] = Reciproc_Assign_State(d_state,tree->io->datatype);
      return;
    }
  else
    {
      int i;
      For(i,3)
	if(d->v[i] != a)
	  Evolve_Recur(d,d->v[i],d->b[i],
		       d_state,r_class,site_num,gen_data,
		       mod,tree);
    }
}

/*********************************************************/

void Site_Diversity(t_tree *tree)
{
  int i,j,k,ns;
  int *div,sum;

  ns = tree->mod->ns;

  div = (int *)mCalloc(ns,sizeof(int));

  Site_Diversity_Post(tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0],tree);
  Site_Diversity_Pre (tree->noeud[0],tree->noeud[0]->v[0],tree->noeud[0]->b[0],tree);

  For(i,2*tree->n_otu-3)
    {
      For(j,ns)
	{
	  tree->t_edges[i]->div_post_pred_left[j] = 0;
	  tree->t_edges[i]->div_post_pred_rght[j] = 0;
	}
    }

  For(i,tree->n_pattern)
    {
      For(j,2*tree->n_otu-3)
	{
	  Binary_Decomposition(tree->t_edges[j]->ui_l[i],div,ns);
	  sum = 0;
	  For(k,ns) sum += div[k];
	  tree->t_edges[j]->div_post_pred_left[sum-1] += tree->data->wght[i];

	  Binary_Decomposition(tree->t_edges[j]->ui_r[i],div,ns);
	  sum = 0;
	  For(k,ns) sum += div[k];
	  tree->t_edges[j]->div_post_pred_rght[sum-1] += tree->data->wght[i];
	}
    }

/*   For(j,2*tree->n_otu-3) */
/*     { */
/*       PhyML_Printf("\n. Edge %4d   div_left = %4d %4d %4d %4d -- div_rght = %4d %4d %4d %4d", */
/* 	     j, */
/* 	     tree->t_edges[j]->div_post_pred_left[0], */
/* 	     tree->t_edges[j]->div_post_pred_left[1], */
/* 	     tree->t_edges[j]->div_post_pred_left[2], */
/* 	     tree->t_edges[j]->div_post_pred_left[3], */
/* 	     tree->t_edges[j]->div_post_pred_rght[0], */
/* 	     tree->t_edges[j]->div_post_pred_rght[1], */
/* 	     tree->t_edges[j]->div_post_pred_rght[2], */
/* 	     tree->t_edges[j]->div_post_pred_rght[3]); */
/*     } */
  
  Free(div);
}

/*********************************************************/

void Site_Diversity_Post(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a)
	  Site_Diversity_Post(d,d->v[i],d->b[i],tree);

      Subtree_Union(d,b,tree);
    }
}

/*********************************************************/

void Site_Diversity_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	if(d->v[i] != a)
	  {
	    Subtree_Union(d,d->b[i],tree);
	    Site_Diversity_Pre(d,d->v[i],d->b[i],tree);
	  }
    }
}

/*********************************************************/

void Subtree_Union(t_node *n, t_edge *b_fcus, t_tree *tree)
{
/*  
           |
	   |<- b_cus
	   |
	   n
          / \
       	 /   \
       	/     \
*/

  int site;
  unsigned int *ui, *ui_v1, *ui_v2;
  
  ui = ui_v1 = ui_v2 = NULL;

  if(n == b_fcus->left)
    {	     
      ui = b_fcus->ui_l;

      ui_v1 = 
      (n == n->b[b_fcus->l_v1]->left)?
      (n->b[b_fcus->l_v1]->ui_r):
      (n->b[b_fcus->l_v1]->ui_l);

      ui_v2 = 
      (n == n->b[b_fcus->l_v2]->left)?
      (n->b[b_fcus->l_v2]->ui_r):
      (n->b[b_fcus->l_v2]->ui_l);
    }
  else
    {
      ui = b_fcus->ui_r;
      
      ui_v1 = 
      (n == n->b[b_fcus->r_v1]->left)?
      (n->b[b_fcus->r_v1]->ui_r):
      (n->b[b_fcus->r_v1]->ui_l);

      ui_v2 = 
      (n == n->b[b_fcus->r_v2]->left)?
      (n->b[b_fcus->r_v2]->ui_r):
      (n->b[b_fcus->r_v2]->ui_l);
    }

  For(site,tree->n_pattern) ui[site] = ui_v1[site] | ui_v2[site];

}

/*********************************************************/

void Binary_Decomposition(int value, int *bit_vect, int size)
{
  int i,cumul;

  For(i,size) bit_vect[i] = 0;
  
  cumul = 0;
  for(i=size-1;i>=0;i--)
    {
      if(value - cumul < (int)POW(2,i))
	{
	  bit_vect[i] = 0;
	}
      else
	{
	  bit_vect[i] = 1;
	  cumul += (int)POW(2,i);
	}
    }
}

/*********************************************************/

void Print_Diversity_Header(FILE *fp, t_tree *tree)
{
/*   PhyML_Fprintf(fp,"t_edge side mean\n");  */
  PhyML_Fprintf(fp,"t_edge side diversity count\n"); 
}

/*********************************************************/

void Print_Diversity(FILE *fp, t_tree *tree)
{
  int ns;
  
  if(tree->io->datatype == NT)      ns = 4;
  else if(tree->io->datatype == AA) ns = 20;

  Print_Diversity_Pre(tree->noeud[0],
		      tree->noeud[0]->v[0],
		      tree->noeud[0]->b[0],
		      fp,
		      tree);

/*       mean_div_left = .0; */
/*       For(k,ns)  */
/* 	{ */
/* 	  mean_div_left += (k+1) * tree->t_edges[j]->div_post_pred_left[k]; */
/* 	} */
/*       mean_div_rght = .0; */
/*       For(k,ns) mean_div_rght += (k+1) * tree->t_edges[j]->div_post_pred_rght[k]; */

/*       mean_div_left /= (phydbl)tree->data->init_len; */
/*       mean_div_rght /= (phydbl)tree->data->init_len; */

/*       PhyML_Fprintf(fp,"%4d 0 %f\n",j,mean_div_left); */
/*       PhyML_Fprintf(fp,"%4d 1 %f\n",j,mean_div_rght); */


/*       mean_div_left = .0; */
/*       For(k,ns) mean_div_left += tree->t_edges[j]->div_post_pred_left[k]; */

/*       mean_div_rght = .0; */
/*       For(k,ns)  */
/* 	{ */
/* 	  mean_div_rght += tree->t_edges[j]->div_post_pred_rght[k]; */
/* 	} */

/*       if((mean_div_left != tree->data->init_len) || (mean_div_rght != tree->data->init_len)) */
/* 	{ */
/* 	  PhyML_Printf("\n. mean_div_left = %f mean_div_rght = %f init_len = %d", */
/* 		 mean_div_left,mean_div_rght,tree->data->init_len); */
/* 	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__); */
/* 	  Warn_And_Exit(""); */
/* 	} */
}

/*********************************************************/

void Print_Diversity_Pre(t_node *a, t_node *d, t_edge *b, FILE *fp, t_tree *tree)
{
  int k,ns;

  ns = -1;

  if(d->tax) return;
  else
    {

      if(tree->io->datatype == NT)      ns = 4;
      else if(tree->io->datatype == AA) ns = 20;

      if(d == b->left) For(k,ns) PhyML_Fprintf(fp,"%4d 0 %2d %4d\n",b->num,k,b->div_post_pred_left[k]);
      else             For(k,ns) PhyML_Fprintf(fp,"%4d 1 %2d %4d\n",b->num,k,b->div_post_pred_rght[k]);

      For(k,3) if(d->v[k] != a) Print_Diversity_Pre(d,d->v[k],d->b[k],fp,tree);
    }

}

/*********************************************************/
/* Estimation of density using kernel smoothing. 
- where : point where I want to estimate the density,
- x : data vector, 
- sample_size :  number of data points in x
*/
phydbl Univariate_Kernel_Density_Estimate(phydbl where, phydbl *x, int sample_size)
{
  phydbl sd,h;
  phydbl density,sqrt2pi,cons;
  int i;

  sqrt2pi = 2.506628;

  sd = SQRT(Var(x,sample_size));
  h = 1.06 * sd * POW(sample_size,-1./5.); /* Quick and dirty way to set the bandwidth */
  
  cons = (1./sample_size) * (1./h) * (1./sqrt2pi);

  density = .0;
  For(i,sample_size) density += EXP(-0.5 * POW((x[i] - where)/h,2));
  density *= cons;

  return density;
}

/*********************************************************/

/* Estimation of a multivariate density using kernel smoothing. 

- where : vector where I want to estimate the density,
- x : data matrix, i.e., sample of vectors, 
- sample_size : number of vectors,
- vect_size : vector length. 

See "Multivariate Density Estimation" by David Scott. pp 150.
*/
phydbl Multivariate_Kernel_Density_Estimate(phydbl *where, phydbl **x, int sample_size, int vect_size)
{
  phydbl sd,*h,cons,density,tmp;
  phydbl _2pi;
  int i,j;

  h = (phydbl *)mCalloc(vect_size,sizeof(phydbl));

  _2pi = 6.283185;
  
  For(i,vect_size)
    {
      sd = SQRT(Var(x[i],sample_size));
/*       h[i] = POW(4./(vect_size+2.),1./(vect_size+4)) * sd * POW(sample_size,-1./(vect_size+4)); */
      h[i] = sd * POW(sample_size,-1./(vect_size+4));
/*       PhyML_Printf("\n. sd = %f, h[i] = %f",sd,h[i]); */
    }

  cons = sample_size;
  For(i,vect_size) cons *= h[i];
  cons *= POW(_2pi,vect_size/2.);
  cons = 1./cons;

  density = .0;
  For(i,sample_size)
    {
      tmp = 1.0;
      For(j,vect_size) 
	{
	  tmp *= EXP(-0.5 * POW((x[j][i] - where[j])/h[j],2));
	}
      density += tmp;
    }
  
  density *= cons;

  Free(h);

  return density;
}

/*********************************************************/

phydbl Var(phydbl *x, int n)
{
  phydbl mean, sum2;
  int i;

  mean = Mean(x,n);

  sum2 = .0;
  For(i,n) sum2 += x[i] * x[i];
  
  return (1./n) * (sum2 - n * POW(mean,2));
}

/*********************************************************/

phydbl Mean(phydbl *x, int n)
{
  int i;
  phydbl sum;

  sum = .0;

  For(i,n) sum += x[i];

  return sum / n;
}

/*********************************************************/

void Best_Of_NNI_And_SPR(t_tree *tree)
{
  if(tree->mod->s_opt->random_input_tree) Speed_Spr_Loop(tree); /* Don't do simultaneous NNIs if starting tree is random */
  else
    {
      t_tree *ori_tree,*best_tree;
      model *ori_mod,*best_mod;
      phydbl *ori_bl,*best_bl;
      phydbl best_lnL,ori_lnL,nni_lnL,spr_lnL;
      int i;

      ori_bl = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
      best_bl = (phydbl *)mCalloc(2*tree->n_otu-3,sizeof(phydbl));
      
      ori_mod   = Copy_Model(tree->mod);
      best_mod  = Copy_Model(tree->mod);
      

      ori_tree = Make_Tree_From_Scratch(tree->n_otu,tree->data);

/*       ori_tree = Make_Tree(tree->n_otu); */
/*       Init_Tree(ori_tree,ori_tree->n_otu); */
/*       Make_All_Tree_Nodes(ori_tree); */
/*       Make_All_Tree_Edges(ori_tree); */
      
      best_tree = Make_Tree_From_Scratch(tree->n_otu,tree->data);

/*       best_tree = Make_Tree(tree->n_otu); */
/*       Init_Tree(best_tree,best_tree->n_otu); */
/*       Make_All_Tree_Nodes(best_tree); */
/*       Make_All_Tree_Edges(best_tree); */

      Copy_Tree(tree,ori_tree);
      Record_Br_Len(tree);
      For(i,2*tree->n_otu-3) ori_bl[i] = tree->t_edges[i]->l;

      best_lnL = UNLIKELY;
      Lk(tree);
      ori_lnL = tree->c_lnL; /* Record likelihood of the starting tree */


      Simu_Loop(tree); /* Perform simultaneous NNIs */
      best_lnL = tree->c_lnL; /* Record the likelihood */
      nni_lnL = tree->c_lnL;
      Copy_Tree(tree,best_tree); /* Record the tree topology and branch lengths */
      Record_Br_Len(tree);
      For(i,2*tree->n_otu-3) best_bl[i] = tree->t_edges[i]->l;
      For(i,2*tree->n_otu-3) best_tree->t_edges[i]->l = best_bl[i];
      Record_Model(tree->mod,best_mod);
      
      Copy_Tree(ori_tree,tree); /* Back to the original tree topology */
      For(i,2*tree->n_otu-3) tree->t_edges[i]->l = ori_bl[i]; /* Back to the original branch lengths */
      Record_Model(ori_mod,tree->mod); /* Back to the original model */
      
      /* Make sure the tree is in its original form */
      Lk(tree);

      if(FABS(tree->c_lnL - ori_lnL) > tree->mod->s_opt->min_diff_lk_global)
	{
	  PhyML_Printf("\n. ori_lnL = %f, c_lnL = %f",ori_lnL,tree->c_lnL);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      Speed_Spr_Loop(tree);
      spr_lnL = tree->c_lnL;
      if(tree->c_lnL > best_lnL)
	{
	  best_lnL = tree->c_lnL;
	  Copy_Tree(tree,best_tree); /* Record tree topology, branch lengths and model parameters */
	  Record_Br_Len(tree);
	  For(i,2*tree->n_otu-3) best_bl[i] = tree->t_edges[i]->l;
	  For(i,2*tree->n_otu-3) best_tree->t_edges[i]->l = best_bl[i];
	  Record_Model(tree->mod,best_mod);
	}
      
      Copy_Tree(best_tree,tree);
      For(i,2*tree->n_otu-3) tree->t_edges[i]->l = best_bl[i];
      Record_Model(best_mod,tree->mod);
      
      /* Make sure the current tree has the best topology, branch lengths and model parameters */
      Lk(tree);
      if(FABS(tree->c_lnL - best_lnL) > tree->mod->s_opt->min_diff_lk_global)
	{
	  PhyML_Printf("\n. best_lnL = %f, c_lnL = %f",best_lnL,tree->c_lnL);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      if(tree->mod->s_opt->print)
	{
	  PhyML_Printf("\n\n. Log likelihood obtained after NNI moves : %f",nni_lnL);
	  PhyML_Printf("\n. Log likelihood obtained after SPR moves : %f",spr_lnL);
	}
      
      Free(ori_bl);
      Free(best_bl);
      
      Free_Tree(ori_tree);
      Free_Tree(best_tree);
      
      Free_Model(ori_mod);
      Free_Model(best_mod);
    }
}

/*********************************************************/

/* Polynomial interpolation. Adapted from "Numerical Recipes in C". 
Press, Flannery, Teukolsky, Vetterling, 1988. 
*/
int Polint(phydbl *xa, phydbl *ya, int n, phydbl x, phydbl *y, phydbl *dy)
{
   int i,m,ns=1;
   phydbl den,dif,dift,ho,hp,w;
   phydbl *c,*d;

   dif=FABS(x-xa[1]);

   c = (phydbl *)mCalloc(n,sizeof(phydbl));
   d = (phydbl *)mCalloc(n,sizeof(phydbl));

   for(i=1;i<=n;i++) 
     {
       if((dift=FABS(x-xa[i])) < dif) 
	 {
	   ns=i;
	   dif=dift;
	 }
       c[i]=ya[i];
       d[i]=ya[i];
     }
   
   *y=ya[ns--];
   
   for (m=1;m<n;m++) 
     {
       for (i=1;i<=n-m;i++) 
	 {
	   ho=xa[i]-x;
	   hp=xa[i+m]-x;
	   w=c[i+1]-d[i];
	   if((den=ho-hp) < SMALL && (den=ho-hp) > -SMALL )
	     {
/* 	       Rprintf("\n. Error in routine POLINT.\n"); */
	       Exit("\n. Error in routine POLINT.\n");
	       return(-1);
	     }
	   den=w/den;
	   d[i]=hp*den;
	   c[i]=ho*den;
	 }
       *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
     }

   Free(d);
   Free(c);
   return(0);
}

/*********************************************************/

void JF(t_tree *tree)
{
  //printing loglk for each site, to compute SH-like tests */
  phydbl sum=0.0;
  PhyML_Printf("\n\nSITES LKS:\n");
  int n_patterns = (int)FLOOR(tree->n_pattern);
  int site=0;
  For(site,n_patterns) {
    int wei=0;
    For(wei,tree->data->wght[site]) {
      PhyML_Printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]);
      sum+=tree->c_lnL_sorted[site] / tree->data->wght[site];
    }
  }
  
  PhyML_Printf("\n\nsum=%f\n\n",sum);
  int i=0;
  For(i,2*tree->n_otu-3)
    {
      if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax))
	{
	  PhyML_Printf("%3d %f %f %f\n",
		 tree->t_edges[i]->bip_score,tree->t_edges[i]->alrt_statistic, tree->t_edges[i]->ratio_test,tree->t_edges[i]->l);
	}
    }
  
  
/*   //printing loglk for each site, to compute SH-like tests */
/*   phydbl sum=0.0; */
/*   PhyML_Printf("\n\nSITES LKS:\n"); */
/*   int n_patterns = (int)FLOOR(tree->n_pattern); */
/*   int site=0; */
/*   For(site,n_patterns) { */
/*     int wei=0; */
/*     For(wei,tree->data->wght[site]) { */
/*       PhyML_Printf("%f\n",tree->c_lnL_sorted[site] / tree->data->wght[site]); */
/*       sum+=tree->c_lnL_sorted[site] / tree->data->wght[site]; */
/*     } */
/*   } */
  
/*   PhyML_Printf("\n\nsum=%f\n\n",sum); */
  
/*   int i=0; */
/*   For(i,2*tree->n_otu-3) */
/*     { */
/*       if((!tree->t_edges[i]->left->tax) && (!tree->t_edges[i]->rght->tax)) */
/* 	{ */
/* 	  PhyML_Printf("%3d %f %f %f\n", */
/* 		 tree->t_edges[i]->bip_score,tree->t_edges[i]->alrt_statistic, tree->t_edges[i]->ratio_test,tree->t_edges[i]->l); */
/* 	} */
/*     } */
}

/*********************************************************/

t_tree *Dist_And_BioNJ(calign *cdata, model *mod, option *io)
{
  t_tree *tree;
  matrix *mat;

  if(!io->quiet) PhyML_Printf("\n. Computing pairwise distances...\n");

  mat = ML_Dist(cdata,mod);
  Fill_Missing_Dist(mat);

  if(!io->quiet) PhyML_Printf("\n. Building BioNJ tree...\n");

  mat->tree = Make_Tree_From_Scratch(cdata->n_otu,cdata);
  Bionj(mat);
  tree      = mat->tree;
  tree->mat = mat;

  return tree;
}

/*********************************************************/

void Add_BioNJ_Branch_Lengths(t_tree *tree, calign *cdata, model *mod)
{
  matrix *mat;

  PhyML_Printf("\n. Computing branch length estimates...\n");

  Order_Tree_CSeq(tree,cdata);
  mat = ML_Dist(cdata,mod);
  mat->tree = tree;
  mat->method = 0;
  Bionj_Br_Length(mat);

  Free_Mat(mat);
}

/*********************************************************/

t_tree *Read_User_Tree(calign *cdata, model *mod, option *io)
{
  t_tree *tree;

  
  PhyML_Printf("\n. Reading tree...\n"); fflush(NULL);
  if(io->n_trees == 1) rewind(io->fp_in_tree);
  tree = Read_Tree_File_Phylip(io->fp_in_tree);
  if(!tree) Exit("\n. Input tree not found...\n");
  /* Add branch lengths if necessary */
  if(!tree->has_branch_lengths) Add_BioNJ_Branch_Lengths(tree,cdata,mod);

  return tree;
}

/*********************************************************/

void Print_Time_Info(time_t t_beg, time_t t_end)
{
  div_t hour,min;

  hour = div(t_end-t_beg,3600);
  min  = div(t_end-t_beg,60  );
  min.quot -= hour.quot*60;

  PhyML_Printf("\n. Time used %dh%dm%ds\n", hour.quot,min.quot,(int)(t_end-t_beg)%60);
  PhyML_Printf("\noooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo\n");
}

/*********************************************************/

char *Bootstrap_From_String(char *s_tree, calign *cdata, model *mod, option *io)
{
  t_tree *tree;

  tree = Read_Tree(s_tree);

  if(!tree)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->mod         = mod;
  tree->io          = io;
  tree->data        = cdata;
  tree->both_sides  = 1;
  tree->n_pattern   = tree->data->crunch_len;

  Order_Tree_CSeq(tree,cdata);
  if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
  Make_Tree_4_Lk(tree,cdata,cdata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);  
  Unscale_Br_Len_Multiplier_Tree(tree);
  Br_Len_Not_Involving_Invar(tree);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);

  printf("\n. %s",Write_Tree(tree));
  fflush(NULL);

  tree->both_sides = 1;
  Lk(tree);



#ifdef MPI
  Bootstrap_MPI(tree);
#else
  Bootstrap(tree);
#endif

  Free(s_tree);  
  s_tree = Write_Tree(tree);
  
  Free_Spr_List(tree);
  Free_One_Spr(tree->best_spr);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);

  return s_tree;
}

/*********************************************************/

char *aLRT_From_String(char *s_tree, calign *cdata, model *mod, option *io)
{
  t_tree *tree;

  tree = Read_Tree(s_tree);

  if(!tree)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  tree->mod         = mod;
  tree->io          = io;
  tree->data        = cdata;
  tree->both_sides  = 1;
  tree->n_pattern   = tree->data->crunch_len;

  Order_Tree_CSeq(tree,cdata);
  if(tree->mod->s_opt->random_input_tree) Random_Tree(tree);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,cdata,cdata->init_len);
  Make_Tree_4_Lk(tree,cdata,cdata->init_len);
  tree->triplet_struct = Make_Triplet_Struct(mod);
  Unscale_Br_Len_Multiplier_Tree(tree);
  Br_Len_Not_Involving_Invar(tree);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);

  tree->both_sides = 1;
  Lk(tree);

  aLRT(tree);
  
  Free(s_tree);  
  s_tree = Write_Tree(tree);

  Free_Spr_List(tree);
  Free_One_Spr(tree->best_spr);
  Free_Triplet(tree->triplet_struct);
  Free_Tree_Pars(tree);
  Free_Tree_Lk(tree);
  Free_Tree(tree);
  
  return s_tree;
}

/*********************************************************/

void Prepare_Tree_For_Lk(t_tree *tree)
{
  Order_Tree_CSeq(tree,tree->data);
  Fill_Dir_Table(tree);
  Update_Dirs(tree);
  Make_Tree_4_Pars(tree,tree->data,tree->data->init_len);
  Make_Tree_4_Lk(tree,tree->data,tree->data->init_len);
  tree->triplet_struct = Make_Triplet_Struct(tree->mod);
  Br_Len_Not_Involving_Invar(tree);
  Unscale_Br_Len_Multiplier_Tree(tree);
  Make_Spr_List(tree);
  Make_Best_Spr(tree);
}

/*********************************************************/

void PhyML_Printf(char *format, ...)
{
  va_list ptr;

 #ifdef MPI
  if(Global_myRank == 0)
    {
      va_start (ptr, format);
      vprintf (format, ptr);
      va_end(ptr);
    }
#else
      va_start (ptr, format);
      vprintf (format, ptr);
      va_end(ptr);
#endif
  
  fflush (NULL);
}

/*********************************************************/

void PhyML_Fprintf(FILE *fp, char *format, ...)
{
  va_list ptr;

#ifdef MPI
  if(Global_myRank == 0)
    {
      va_start (ptr, format);
      vfprintf (fp,format, ptr);
      va_end(ptr);
    }
#else
      va_start (ptr, format);
      vfprintf (fp,format, ptr);
      va_end(ptr);
#endif
  
  fflush (NULL);
}

/*********************************************************/

void Find_Common_Tips(t_tree *tree1, t_tree *tree2)
{
  int i,j;

  For(i,tree1->n_otu) tree1->noeud[i]->common = 0;
  For(i,tree2->n_otu) tree2->noeud[i]->common = 0;

  For(i,tree1->n_otu)
    {
      For(j,tree2->n_otu)
	{
	  if(!strcmp(tree1->noeud[i]->name,tree2->noeud[j]->name))
	    {
	      tree1->noeud[i]->common = 1;
	      tree2->noeud[j]->common = 1;
	      break;
	    }
	}
    }
}

/*********************************************************/

phydbl Get_Tree_Size(t_tree *tree)
{
  int i;
  phydbl tree_size;

  tree_size = 0.0;
  For(i,2*tree->n_otu-3) tree_size += tree->t_edges[i]->l;
/*   tree_size = 0.0; */
/*   For(i,2*tree->n_otu-3) tree_size += tree->rates->u_cur_l[i]; */

/*   For(i,2*tree->n_otu-3)  */
/*     tree_size +=  */
/*     FABS(tree->rates->nd_t[tree->t_edges[i]->left->num] -  */
/* 	 tree->rates->nd_t[tree->t_edges[i]->rght->num]); */

  tree->size = tree_size;
  return tree_size;
}

/*********************************************************/

void Dist_To_Root_Pre(t_node *a, t_node *d, t_edge *b, t_tree *tree)
{
  int i;

  if(b) d->dist_to_root = a->dist_to_root + b->l;

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root)) 
	  Dist_To_Root_Pre(d,d->v[i],d->b[i],tree);
    }
}

/*********************************************************/

void Dist_To_Root(t_node *n_root, t_tree *tree)
{  
/*   n_root->v[0]->dist_to_root = tree->rates->cur_l[n_root->v[0]->num]; */
/*   n_root->v[1]->dist_to_root = tree->rates->cur_l[n_root->v[1]->num]; */
  n_root->v[0]->dist_to_root = tree->e_root->l * tree->n_root_pos;
  n_root->v[1]->dist_to_root = tree->e_root->l * (1. - tree->n_root_pos);
  Dist_To_Root_Pre(n_root,n_root->v[0],NULL,tree);
  Dist_To_Root_Pre(n_root,n_root->v[1],NULL,tree);
}

/*********************************************************/
/* 'Borrowed' fromn libgen */
char *Basename(char *path)
{
  char *p;
  
  if( path == NULL || *path == '\0' ) return ".";

  p = path + strlen(path) - 1;
  
  while( *p == '/' ) 
    {
      if( p == path ) return path;
      *p-- = '\0';
    }

  while( p >= path && *p != '/' ) p--;

  return p + 1;
}

/*********************************************************/


/* Find the Last Common Ancestor of n1 and n2 */
t_node *Find_Lca(t_node *n1, t_node *n2, t_tree *tree)
{
  t_node **list1, **list2, *lca;
  int size1, size2;

  if(!tree->n_root)
    {
      PhyML_Printf("\n. The tree must be rooted in this function.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  list1 = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));
  list2 = (t_node **)mCalloc(2*tree->n_otu-1,sizeof(t_node *));

  Get_List_Of_Ancestors(n1,list1,&size1,tree);
  Get_List_Of_Ancestors(n2,list2,&size2,tree);

  while(list1[size1] == list2[size2])
    {
      size1--;
      size2--;

      if(size1 < 0 || size2 < 0) break;
    }
  
  lca = list1[size1+1];

  Free(list1);
  Free(list2);


/*   t_node *lca; */
/*   int score; */
/*   int i,j; */

/*   lca = n1->anc; */
/*   do */
/*     { */
/*       if(lca == tree->n_root) break; */
/*       else */
/* 	{ */
/* 	  score = 0; */
/* 	  For(i,3) if(lca->v[i] == lca->anc) break; */
/* 	  For(j,lca->bip_size[i]) */
/* 	    { */
/* 	      if(lca->bip_node[i][j] == n1 || lca->bip_node[i][j] == n2) */
/* 		{ */
/* 		  score++; */
/* 		  if(score == 2) break; */
/* 		} */
/* 	    } */
/* 	  if(score != 2) lca = lca->anc; */
/* 	} */
/*     }while(score != 2); */
  

  return lca;
}

/*********************************************************/

/* Returns the list of the ancestors of ref_t_node from ref_t_node to the root included */
void Get_List_Of_Ancestors(t_node *ref_node, t_node **list, int *size, t_tree *tree)
{
  t_node *n;

  *size = 0;
  n = ref_node;
  list[0] = n;

  while(n != tree->n_root)
    {
      n = n->anc;
      *size = *size+1;
      list[*size] = n;
    }
}

/*********************************************************/

int Edge_Num_To_Node_Num(int edge_num, t_tree *tree)
{
  int node_num;
  t_edge *b;

  b = tree->t_edges[edge_num];

  node_num = (b->left == b->rght->anc)?(b->rght->num):(b->left->num);
  
  return node_num;
}

/*********************************************************/

void Branch_Lengths_To_Time_Lengths(t_tree *tree)
{
  Branch_Lengths_To_Time_Lengths_Pre(tree->n_root,tree->n_root->v[0],tree);
  Branch_Lengths_To_Time_Lengths_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void Branch_Lengths_To_Time_Lengths_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  tree->rates->cur_l[d->num] = fabs(tree->rates->nd_t[d->num] - tree->rates->nd_t[a->num]);

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  Branch_Lengths_To_Time_Lengths_Pre(d,d->v[i],tree);     
    }
}

/*********************************************************/

void Branch_Lengths_To_Rate_Lengths(t_tree *tree)
{
  Branch_Lengths_To_Rate_Lengths_Pre(tree->n_root,tree->n_root->v[0],tree);
  Branch_Lengths_To_Rate_Lengths_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void Branch_Lengths_To_Rate_Lengths_Pre(t_node *a, t_node *d, t_tree *tree)
{
  int i;

  tree->rates->cur_l[d->num] = 
    tree->rates->br_r[d->num] * 
    tree->rates->clock_r *
    tree->rates->norm_fact;

  if(d->tax) return;
  else
    {
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  Branch_Lengths_To_Rate_Lengths_Pre(d,d->v[i],tree);     
    }
}

/*********************************************************/

int Find_Clade(char **tax_name_list, int list_size, t_tree *tree)
{
  int *tax_num_list;
  int i,j;
  int n_matches;

  tax_num_list = (int *)mCalloc(list_size,sizeof(int));
  
  n_matches = 0;
  For(i,list_size)
    {
      For(j,tree->n_otu)
	{
	  if(!strcmp(tax_name_list[i],tree->noeud[j]->name))
	    {
	      tax_num_list[i] = tree->noeud[j]->num;
	      n_matches++;
	    }
	}
    }

  if(n_matches != list_size)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(list_size == tree->n_otu)
    {
      int i,j;
      int score;
      
      score = 0;
      For(i,list_size)
	{
	  For(j,tree->n_otu)
	    {
/* 	      if(!strcmp(tax_name_list[i],tree->noeud[j]->name)) score++; */
	      if(tax_num_list[i] == tree->noeud[j]->num) score++;
	    }
	}

      Free(tax_num_list);
      if(score == tree->n_otu) return tree->n_root->num;
      else return -1;
    }
  else
    {
      int num;
      num = -1;
      Free_Bip(tree);
      Alloc_Bip(tree);
      Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);
      Find_Clade_Pre(tree->n_root,tree->n_root->v[0],tax_num_list,list_size,&num,tree);
      Find_Clade_Pre(tree->n_root,tree->n_root->v[1],tax_num_list,list_size,&num,tree);
      Free(tax_num_list);
      return num;
    }
  
  Free(tax_num_list);
  return -1;
}

/*********************************************************/

void Find_Clade_Pre(t_node *a, t_node *d, int *tax_num_list, int list_size, int *num, t_tree *tree)
{
  int i,j,k;
  int score;
  
  
  For(i,3)
    if((d->v[i] == a) || (d->b[i] == tree->e_root))
      {
	if(list_size == d->bip_size[i])
	  {
	    score = 0;
	    For(j,d->bip_size[i])
	      {
		For(k,list_size)
		  {
		    if(tax_num_list[k] == d->bip_node[i][j]->num)
		      {
			score++;
			break;
		      }
		  }
	      }
	    if(score == list_size) *num = d->num;
	  }
	break;
      }

  if(d->tax) return;
  else
    For(i,3)
      if((d->v[i] != a) && (d->b[i] != tree->e_root))
	Find_Clade_Pre(d,d->v[i],tax_num_list,list_size,num,tree);
}

/*********************************************************/

void Read_Clade_Priors(char *file_name, t_tree *tree)
{
  FILE *fp;
  char *s,*line;
  int n_clade_priors;
  int clade_size;
  char **clade_list;
  int i,pos;
  phydbl prior_low,prior_up;
  int node_num;

  PhyML_Printf("\n. Reading prior on node ages.\n");

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));
  s    = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  clade_list = (char **)mCalloc(tree->n_otu,sizeof(char *));
  For(i,tree->n_otu) clade_list[i] = (char *)mCalloc(T_MAX_NAME,sizeof(char));

  fp = Openfile(file_name,0);
  
  n_clade_priors = 0;  
  do
    {
      if(!fgets(line,T_MAX_LINE,fp)) break;

      clade_size = 0;
      pos = 0;
      do
	{
	  i = 0;

	  while(line[pos] == ' ') pos++;

	  while((line[pos] != ' ') && (line[pos] != '\n') && line[pos] != '#')
	    {
	      s[i] = line[pos];
	      i++;
	      pos++;
	    }
	  s[i] = '\0';
/* 	  PhyML_Printf("\n. s = %s\n",s); */
	  
	  if(line[pos] == '\n' || line[pos] == '#') break;
	  pos++;

	  if(strcmp(s,"|"))
	    {
	      strcpy(clade_list[clade_size],s);
	      clade_size++;
	    }
	  else break;	    
	}
      while(1);

      if(line[pos] != '#' && line[pos] != '\n')
	{
	  double val1, val2;
/* 	  sscanf(line+pos,"%lf %lf",&prior_up,&prior_low); */
	  sscanf(line+pos,"%lf %lf",&val1,&val2);
	  prior_up = (phydbl)val1;
	  prior_low = (phydbl)val2;
	  node_num = -1;
	  if(!strcmp("@root@",clade_list[0])) node_num = tree->n_root->num;
	  else node_num = Find_Clade(clade_list, clade_size, tree);
	  
	  n_clade_priors++;  

	  if(node_num < 0)
	    {
	      PhyML_Printf("\n");
	      PhyML_Printf("\n");
	      PhyML_Printf("\n. >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>");
	      PhyML_Printf("\n. Could not find any clade in the tree with the following taxa names:");
	      For(i,clade_size) PhyML_Printf("\n. \"%s\"",clade_list[i]);
	      PhyML_Printf("\n. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<");
	      Exit("\n");
	    }
	  else
	    {	      
	      tree->rates->t_has_prior[node_num] = YES;
	      tree->rates->t_prior_min[node_num] = MIN(prior_low,prior_up);
	      tree->rates->t_prior_max[node_num] = MAX(prior_low,prior_up);

	      if(FABS(prior_low - prior_up) < 1.E-6 && tree->noeud[node_num]->tax == YES)
		tree->rates->nd_t[node_num] = prior_low;

/* 	      PhyML_Printf("\n"); */
/* 	      PhyML_Printf("\n. [%3d]>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>",n_clade_priors); */
/* 	      PhyML_Printf("\n. Node %4d matches the clade with the following taxa names:",node_num); */
/* 	      For(i,clade_size) PhyML_Printf("\n. - \"%s\"",clade_list[i]); */
/* 	      PhyML_Printf("\n. Lower bound set to: %15f time units.",MIN(prior_low,prior_up));  */
/* 	      PhyML_Printf("\n. Upper bound set to: %15f time units.",MAX(prior_low,prior_up));  */
/* 	      PhyML_Printf("\n. <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"); */
	    }
	}
    }
  while(1);
  
  
  PhyML_Printf("\n. Read prior information on %d %s.",n_clade_priors,n_clade_priors > 1 ? "clades":"clade");

  if(!n_clade_priors)
    {
      PhyML_Printf("\n. PhyTime could not find any prior on node age.");
      PhyML_Printf("\n. This is likely due to a problem in the calibration ");
      PhyML_Printf("\n. file format. Make sure, for instance, that there is ");
      PhyML_Printf("\n. a blank character between the end of the last name");
      PhyML_Printf("\n. of each clade and the character `|'. Otherwise, ");
      PhyML_Printf("\n. please refer to the example file.\n");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  For(i,tree->n_otu) Free(clade_list[i]);
  Free(clade_list);
  Free(line);
  Free(s);
  fclose(fp);
}

/*********************************************************/

t_edge *Find_Root_Edge(FILE *fp_input_tree, t_tree *tree)
{
  char **subs;
  int degree;
  int i,j;
  t_node *left, *rght;
  int l_r, r_l;
  int score;
  char *line;
  char c;
  t_edge *root_edge;

  line = (char *)mCalloc(T_MAX_LINE,sizeof(char));

  rewind(fp_input_tree);

  do c=fgetc(fp_input_tree);
  while((c != '(') && (c != EOF));
  
  if(c==EOF)
    {
      Free(line);
      return NULL;
    }
  
  i=0;
  for(;;)
    {
      if((c == ' ') || (c == '\n'))
	{
	  c=fgetc(fp_input_tree);
	  if(c==EOF) break;
	  else continue;
	}
      
      line[i]=c;
      i++;
      c=fgetc(fp_input_tree);
      if(c==EOF || c==';') break;
    }
  

  Free_Bip(tree);
  Alloc_Bip(tree);
  Get_Bip(tree->noeud[0],tree->noeud[0]->v[0],tree);

  subs = Sub_Trees(line,&degree);
  Clean_Multifurcation(subs,degree,3);
  if(degree != 2) 
    {
      PhyML_Printf("\n. The tree does not seem to be rooted...");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  left = rght = NULL;
  l_r = r_l = -1;

  For(i,2*tree->n_otu-3)
    {
      left = tree->t_edges[i]->left;
      rght = tree->t_edges[i]->rght;
      l_r  = tree->t_edges[i]->l_r;
      r_l  = tree->t_edges[i]->r_l;
      
      score = 0;
      For(j,left->bip_size[l_r]) if(strstr(subs[1],left->bip_node[l_r][j]->name)) score++;
      if(score == left->bip_size[l_r]) break;

      score = 0;
      For(j,rght->bip_size[r_l]) if(strstr(subs[1],rght->bip_node[r_l][j]->name)) score++;
      if(score == rght->bip_size[r_l]) break;
    }

  root_edge = tree->t_edges[i];

  For(i,NODE_DEG_MAX) Free(subs[i]);
  Free(subs);
  Free(line);

  if(i == 2*tree->n_otu-3)
    {
      PhyML_Printf("\n. Could not find the root edge...");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  
  return root_edge;
}

/*********************************************************/

void Copy_Tree_Topology_With_Labels(t_tree *ori, t_tree *cpy)
{
  int i,j;

  For(i,2*ori->n_otu-2)
    {
      For(j,3)
        {
          if(ori->noeud[i]->v[j])
            {
              cpy->noeud[i]->v[j] = cpy->noeud[ori->noeud[i]->v[j]->num];
              cpy->noeud[i]->l[j] = ori->noeud[i]->l[j];
            }
          else
            cpy->noeud[i]->v[j] = NULL;
        }
      cpy->noeud[i]->num = ori->noeud[i]->num;
      cpy->noeud[i]->tax = 0;
    }

  For(i,2*ori->n_otu-3)
    {
      cpy->t_edges[i]->l = ori->t_edges[i]->l;
    }

  For(i,ori->n_otu)
    {
      cpy->noeud[i]->tax = 1;
      strcpy(cpy->noeud[i]->name,ori->noeud[i]->name);
    }

}

/*********************************************************/

option *Get_Input(int argc, char **argv)
{

  option *io;
  model *mod;
  optimiz *s_opt;
  m4 *m4mod;

  io    = (option *)Make_Input();
  mod   = (model *)Make_Model_Basic();
  s_opt = (optimiz *)Make_Optimiz();
  m4mod = (m4 *)M4_Make_Light();

  Set_Defaults_Input(io);
  Set_Defaults_Model(mod);
  Set_Defaults_Optimiz(s_opt);

  io->mod        = mod;
  io->mod->m4mod = m4mod;

  mod->io        = io;
  mod->s_opt     = s_opt;


#ifdef MPI
  Read_Command_Line(io,argc,argv);
#elif defined (PHYTIME)
  Read_Command_Line(io,argc,argv);
#else
  putchar('\n');

  switch (argc)
    {
    case 1:
      {
	Launch_Interface(io);
	break;
      }
      /*
	case 2:
	Usage();
	break;
      */
    default:
      Read_Command_Line(io,argc,argv);
    }
#endif

  return io;
}

/*********************************************************/

void Set_Model_Name(model *mod)
{
  if(mod->io->datatype == NT)
    {
      switch(mod->whichmodel)
	{
	case JC69:
	  {
	    strcpy(mod->modelname, "JC69");
	    break;
	  }
	case K80:
	  {
	    strcpy(mod->modelname, "K80");
	    break;
	  }
	case F81:
	  {
	    strcpy(mod->modelname, "F81");
	    break;
	  }
	case HKY85:
	  {
	    strcpy(mod->modelname, "HKY85");
	    break;
	  }
	case F84:
	  {
	    strcpy(mod->modelname, "F84");
	    break;
	  }
	case TN93:
	  {
	    strcpy(mod->modelname, "TN93");
	    break;
	  }
	case GTR:
	  {
	    strcpy(mod->modelname, "GTR");
	    break;
	  }
	case CUSTOM:
	  {
	    strcpy(mod->modelname, "Custom");
	    break;
	  }
	default:
	  {
	    PhyML_Printf("\n. Unknown model name.\n");
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else if(mod->io->datatype == AA)
    {
      switch(mod->whichmodel)
	{
	case DAYHOFF:
	  {
	    strcpy(mod->modelname, "Dayhoff");
	    break;
	  }
	case JTT:
	  {
	    strcpy(mod->modelname, "JTT");
	    break;
	  }
	case MTREV:
	  {
	    strcpy(mod->modelname, "MtREV");
	    break;
	  }
	case LG:
	  {
	    strcpy(mod->modelname, "LG");
	    break;
	  }
	case WAG:
	  {
	    strcpy(mod->modelname, "WAG");
	    break;
	  }
	case DCMUT:
	  {
	    strcpy(mod->modelname, "DCMut");
	    break;
	  }
	case RTREV:
	  {
	    strcpy(mod->modelname, "RtREV");
	    break;
	  }
	case CPREV:
	  {
	    strcpy(mod->modelname, "CpREV");
	    break;
	  }
	case VT:
	  {
	    strcpy(mod->modelname, "VT");
	    break;
	  }
	case BLOSUM62:
	  {
	    strcpy(mod->modelname, "Blosum62");
	    break;
	  }
	case MTMAM:
	  {
	    strcpy(mod->modelname, "MtMam");
	    break;
	  }
	case MTART:
	  {
	    strcpy(mod->modelname, "MtArt");
	    break;
	  }
	case HIVW:
	  {
	    strcpy(mod->modelname, "HIVb");
	    break;
	  }
	case HIVB:
	  {
	    strcpy(mod->modelname, "HIVb");
	    break;
	  }
	case CUSTOMAA:
	  {
	    strcpy(mod->modelname, "Custom");
	    break;
	  }
	default:
	  {
	    PhyML_Printf("\n. Unknown model name.\n");
	    PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	    Warn_And_Exit("");
	    break;
	  }
	}
    }
  else if(mod->io->datatype == GENERIC)
    {
      strcpy(mod->modelname, "JC69");
    }
}

/*********************************************************/

void Adjust_Min_Diff_Lk(t_tree *tree)
{
  int exponent;
  
  exponent = (int)FLOOR(log10(FABS(tree->c_lnL)));

  if(sizeof(phydbl) == 4)
    {
      tree->mod->s_opt->min_diff_lk_global = POW(10.,exponent - FLT_DIG + 1);
      tree->mod->s_opt->min_diff_lk_local  = tree->mod->s_opt->min_diff_lk_global;
      tree->mod->s_opt->min_diff_lk_move   = tree->mod->s_opt->min_diff_lk_global;
    }
/*   PhyML_Printf("\n. Exponent = %d Precision = %E DIG = %d",exponent,tree->mod->s_opt->min_diff_lk_global,FLT_DIG); */
}

/*********************************************************/

/*!
  tree->noeud[i]->name is initially a number. It is translated into
  a string of characters using the names provided in the tax_name
  array.
 */
void Translate_Tax_Names(char **tax_names, t_tree *tree)
{
  int i;
  int tax_num;

  For(i,tree->n_otu)
    {
      tax_num = strtol(tree->noeud[i]->name,NULL,10);
      tree->noeud[i]->name = tax_names[tax_num-1];
    }
}

/*********************************************************/

/*!
  Skip coment in NEXUS file.
 */
void Skip_Comment(FILE *fp)
{
  int in_comment;
  char c;

  in_comment = 1;
  do
    {
      c = fgetc(fp);
      if(c == EOF) break;
      if(c == '[')      in_comment++;
      else if(c == ']') in_comment--;
    }
  while(in_comment);
}

/*********************************************************/
/*!
  Determine the most appropriate position of the root if outgroup taxa are specified. 
 */

void Get_Best_Root_Position(t_tree *tree)
{
  int i,j;
  phydbl eps;
  phydbl s, s_max;
  t_edge *best_edge;
  int has_outgrp;
  
  best_edge = NULL;

  if(tree->n_root)
    {
      PhyML_Printf("\n. Tree already has a root.");
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(strstr(tree->noeud[0]->name,"*"))
    {
      /* PhyML_Printf("\n. Found outgroup taxon: %s",tree->noeud[0]->name); */
      tree->noeud[0]->s_ingrp[0]  = 0;
      tree->noeud[0]->s_outgrp[0] = 1;
    }
  else
    {
      tree->noeud[0]->s_ingrp[0]  = 1;
      tree->noeud[0]->s_outgrp[0] = 0;
    }

  has_outgrp = NO;
  Get_Best_Root_Position_Post(tree->noeud[0],tree->noeud[0]->v[0],&has_outgrp,tree);  
  Get_Best_Root_Position_Pre(tree->noeud[0],tree->noeud[0]->v[0],tree);  

  if(has_outgrp == YES)
    {
      eps = 1.E-10;
      s = s_max = 0.0;
      For(i,2*tree->n_otu-2)
	{
	  For(j,3)
	    {
	      s = (tree->noeud[i]->s_outgrp[j]+eps) / (tree->noeud[i]->s_ingrp[j] + eps) ;
	      /* printf("\n. [%d %d] %d %d",i,j,tree->noeud[i]->s_outgrp[j],tree->noeud[i]->s_ingrp[j]); */
	      if(s > s_max) 
		{
		  s_max = s;
		  best_edge = tree->noeud[i]->b[j];
		}
	    }
	}
      
      Add_Root(best_edge,tree);
    }
}

/*********************************************************/

/*!
  Determine the most appropriate position of the root if outgroup taxa are specified. 
  Post-traversal.
 */
void Get_Best_Root_Position_Post(t_node *a, t_node *d, int *has_outgrp, t_tree *tree)
{
  if(d->tax) 
    {
      if(strstr(d->name,"*"))
	{
	  *has_outgrp = YES;
	  /* PhyML_Printf("\n. Found outgroup taxon: %s",d->name); */
	  d->s_ingrp[0]  = 0;
	  d->s_outgrp[0] = 1;
	}
      else
	{
	  d->s_ingrp[0]  = 1;
	  d->s_outgrp[0] = 0;
	}
      return;
    }
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a && (d->b[i] != tree->e_root))
	  Get_Best_Root_Position_Post(d,d->v[i],has_outgrp,tree);

      Get_OutIn_Scores(a,d);

    }
}

/*********************************************************/

/*!
  Determine the most appropriate position of the root if outgroup taxa are specified. 
  Pre-traversal.
 */
void Get_Best_Root_Position_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) 
    {
      return;
    }
  else
    {
      int i;

      For(i,3)
	if(d->v[i] != a && (d->b[i] != tree->e_root))
	  {
	    Get_OutIn_Scores(d->v[i],d);
	    Get_Best_Root_Position_Pre(d,d->v[i],tree);
	  }
    }
}

/*********************************************************/

/*!
  Determine the most appropriate position of the root if outgroup taxa are specified. 
  Core.
 */
void Get_OutIn_Scores(t_node *a, t_node *d)
{
  int i,d_v1,d_v2,v1_d,v2_d,d_a;
  
  d_a = v1_d = v2_d = -1;
  d_v1 = d_v2 = -1;
  For(i,3)
    {
      if(d->v[i] != a)
	{
	  if(d_v1 < 0) d_v1 = i;
	  else         d_v2 = i;
	}
    }
  
  For(i,3) if(d->v[i] == a) { d_a = i; break; }
  For(i,3) if(d->v[d_v1]->v[i] == d) { v1_d = i; break; }
  For(i,3) if(d->v[d_v2]->v[i] == d) { v2_d = i; break; }
  
  d->s_ingrp[d_a] = 
    d->v[d_v1]->s_ingrp[v1_d] +
    d->v[d_v2]->s_ingrp[v2_d] ;
  
  d->s_outgrp[d_a] = 
    d->v[d_v1]->s_outgrp[v1_d] +
    d->v[d_v2]->s_outgrp[v2_d] ;
}

/*********************************************************/

int Check_Sequence_Name(char *s)
{
  if(rindex(s,':'))
    {
      PhyML_Printf("\n. Character ':' is not permitted in sequence name (%s).",s);
      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(rindex(s,','))
    {
      PhyML_Printf("\n. Character ',' is not permitted in sequence name (%s).",s);
      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  if(rindex(s,' '))
    {
      PhyML_Printf("\n. Character ' ' is not permitted in sequence name (%s).",s);
      PhyML_Printf("\n. Err in file %s at line %d",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  return 1;
}

/*********************************************************/

int Scale_Subtree_Height(t_node *a, phydbl K, phydbl floor, int *n_nodes, t_tree *tree)
{
  phydbl new_height;

  *n_nodes = 0;
  
  new_height = .0;


  if(!(tree->rates->nd_t[a->num] > floor))
    new_height = K*(tree->rates->nd_t[a->num]-floor)+floor;

  if(a == tree->n_root)
    {
      tree->rates->nd_t[tree->n_root->num] = new_height;
      *n_nodes = 1;
      Scale_Node_Heights_Post(tree->n_root,tree->n_root->v[0],K,floor,n_nodes,tree);
      Scale_Node_Heights_Post(tree->n_root,tree->n_root->v[1],K,floor,n_nodes,tree);
    }
  else
    {
      int i;
      
      if(new_height < tree->rates->nd_t[a->anc->num]) return 0;
      else 
	{
	  tree->rates->nd_t[a->num] = new_height;
	  *n_nodes = 1;
	}

      For(i,3)
	if(a->v[i] != a->anc && a->b[i] != tree->e_root)
	  {
	    Scale_Node_Heights_Post(a,a->v[i],K,floor,n_nodes,tree);
	  }
    }
  
  return 1;
}

/*********************************************************/

void Scale_Node_Heights_Post(t_node *a, t_node *d, phydbl K, phydbl floor, int *n_nodes, t_tree *tree)
{
  if(d == tree->n_root)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }

  if(d->tax) 
    {      
      if(!(tree->rates->nd_t[d->num] > floor)) 
	{
	  /* Scaling does not change the node height here 
	     but, in theory this tip is among the scaled
	     nodes. Therefore needs to count it in for 
	     working out correct Hastings ratios
	  */
	  /* *n_nodes = *n_nodes+1;  */
	}
      return;
    }
  else
    {
      int i;
      
      /* It is tempting to set floor = tree->rates->t_prior_max[d->num]; but
	 it then becomes possible for nodes with different floor values
	 to have their orders interverted (i.e., ancestor below descendant)
      */
      if(!(tree->rates->nd_t[d->num] > floor))
	{
	  tree->rates->nd_t[d->num] = K*(tree->rates->nd_t[d->num]-floor)+floor;
	  *n_nodes = *n_nodes+1;
	}

      if(tree->rates->nd_t[d->num] < tree->rates->nd_t[a->num])
	{
	  PhyML_Printf("\n. K = %f floor = %f t_prior_max(a) = %f t_prior_max(d) = %f a->t = %f d->t %f",
		       K,floor,tree->rates->t_prior_max[a->num],tree->rates->t_prior_max[d->num],
		       tree->rates->nd_t[a->num],tree->rates->nd_t[d->num]);
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      For(i,3)
	if(d->v[i] != a && d->b[i] != tree->e_root)
	  Scale_Node_Heights_Post(d,d->v[i],K,floor,n_nodes,tree);

    }
}

/*********************************************************/

int Scale_Subtree_Rates(t_node *a, phydbl mult, int *n_nodes, t_tree *tree)
{
  int res;
  int i;

  *n_nodes = 0;
  res      = 1;

  if(a == tree->n_root)
    {
      res = Scale_Subtree_Rates_Post(a,a->v[0],mult,n_nodes,tree);
      if(res) res = Scale_Subtree_Rates_Post(a,a->v[1],mult,n_nodes,tree);
      return res;
    }
  else
    {
      For(i,3) if((a->v[i] != a->anc) && 
		  (a->b[i] != tree->e_root) && 
		  (res == 1)) res = Scale_Subtree_Rates_Post(a,a->v[i],mult,n_nodes,tree);
      return res;
    }
}

/*********************************************************/

void Check_Br_Len_Bounds(t_tree *tree)
{
  int i;
  t_edge *b;

  b = NULL;

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      if(b->l > tree->mod->l_max) b->l = tree->mod->l_max;
      if(b->l < tree->mod->l_min) b->l = tree->mod->l_min;
    }

  if(tree->rates)
    {
      For(i,2*tree->n_otu-3)
	{
	  if(tree->rates->u_cur_l[i] > tree->mod->l_max) tree->rates->u_cur_l[i] = tree->mod->l_max;
	  if(tree->rates->u_cur_l[i] < tree->mod->l_min) tree->rates->u_cur_l[i] = tree->mod->l_min;
	}
    }
}

/*********************************************************/

int Scale_Subtree_Rates_Post(t_node *a, t_node *d, phydbl mult, int *n_nodes, t_tree *tree)
{

  if(tree->rates->model_log_rates == YES)
    {
      tree->rates->br_r[d->num] += LOG(mult);
    }
  else
    {
      tree->rates->br_r[d->num] *= mult;
    }

  *n_nodes = *n_nodes+1;

  if(tree->rates->br_r[d->num] < tree->rates->min_rate) return 0;
  if(tree->rates->br_r[d->num] > tree->rates->max_rate) return 0;

  if(d->tax) return 1;
  else
    {
      int i,res;

      res = 1;
      For(i,3)
	{
	  if((d->v[i] != a) && 
	     (d->b[i] != tree->e_root) && 
	     (res == 1))
	    {
	      res = Scale_Subtree_Rates_Post(d,d->v[i],mult,n_nodes,tree);
	    }
	}
      return res;
    }
}

/*********************************************************/
/*********************************************************/


void Get_Node_Ranks(t_tree *tree)
{
  tree->n_root->rank = 1;
  Get_Node_Ranks_Pre(tree->n_root,tree->n_root->v[0],tree);
  Get_Node_Ranks_Pre(tree->n_root,tree->n_root->v[1],tree);
}

/*********************************************************/

void Get_Node_Ranks_Pre(t_node *a, t_node *d,t_tree *tree)
{
  d->rank = a->rank+1;

  if(d->tax) return;
  else
    {
      int i;
      
      For(i,3)
	{
	  if(d->v[i] != a && d->b[i] != tree->e_root)
	    {
	      Get_Node_Ranks_Pre(d,d->v[i],tree);
	    }
	}
    }
}

/*********************************************************/

void Log_Br_Len(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l = LOG(tree->t_edges[i]->l);
}

/*********************************************************/

void Make_Short_L(t_tree *tree)
{
  if(!tree->short_l)
    tree->short_l = (phydbl *)mCalloc(tree->n_short_l,sizeof(phydbl));
}

/*********************************************************/

phydbl Diff_Lk_Norm_At_Given_Edge(t_edge *b, t_tree *tree)
{
  int i,dim,err;
  phydbl lk_exact,lk_norm,sum;

  Record_Br_Len(tree);

  dim = 2*tree->n_otu-3;
  sum = 0.0;

  For(i,tree->n_short_l)
    {
      b->l = tree->short_l[i];

      lk_exact = Lk_At_Given_Edge(b,tree);
      lk_norm = tree->norm_scale + Log_Dnorm(b->l,tree->rates->mean_l[b->num],
					     tree->rates->cov_l[b->num*dim+b->num],&err);

      if(err)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}

      sum += pow(lk_exact - lk_norm,2);      
    }

  Restore_Br_Len(tree);
  Lk_At_Given_Edge(b,tree);

  return(sum);
}


/*********************************************************/

void Adjust_Variances(t_tree *tree)
{
  int i;
  phydbl new_diff,curr_diff;

  Make_Short_L(tree);
  For(i,tree->n_short_l)
    {
      tree->short_l[i] = tree->mod->l_min + i*(0.1 - tree->mod->l_min)/tree->n_short_l;
    }
  

  For(i,2*tree->n_otu-3)
    {
      if(tree->t_edges[i]->l < 1.1*tree->mod->l_min)
	{
	  tree->rates->mean_l[i]                     = -1.00;
	  tree->rates->cov_l[i*(2*tree->n_otu-3)+i]  =  0.1;
	  tree->norm_scale                           = -100;
	  
	  
	  new_diff = curr_diff = 10.0;
	  do
	    {
	      curr_diff = new_diff;
	      
	      Generic_Brent_Lk(&(tree->norm_scale),
			       -1E+6,
			       0.0,
			       1.E-10,
			       10000,
			       NO,
			       Wrap_Diff_Lk_Norm_At_Given_Edge,tree->t_edges[i],tree,NULL);
	      
	      /* 		      Generic_Brent_Lk(&(tree->rates->mean_l[0]), */
	      /* 				       -100., */
	      /* 				       10*tree->mod->l_min, */
	      /* 				       1.E-3, */
	      /* 				       10000, */
	      /* 				       NO, */
	      /* 				       Wrap_Diff_Lk_Norm_At_Given_Edge,tree->t_edges[0],tree,NULL); */
	      
	      Generic_Brent_Lk(&(tree->rates->cov_l[i*(2*tree->n_otu-3)+i]),
			       0.0,
			       10.0,
			       1.E-10,
			       10000,
			       NO,
			       Wrap_Diff_Lk_Norm_At_Given_Edge,tree->t_edges[i],tree,NULL);
	      	  
	      new_diff = Diff_Lk_Norm_At_Given_Edge(tree->t_edges[i],tree);
	    }while(FABS(new_diff-curr_diff) > 1.E-3);
	}
    }
}

/*********************************************************/

phydbl Effective_Sample_Size(phydbl first_val, phydbl last_val, phydbl sum, phydbl sumsq, phydbl sumcurnext, int n)
{
  phydbl numerator,denom;
  phydbl mean;
  phydbl r;

  mean = sum / n;
  denom = sumsq - n * POW(mean,2);
  numerator = sumcurnext - (n+1.)*POW(mean,2) + (first_val+last_val)*mean;

  r = numerator/denom;
  
  return (phydbl)n * (1.-r)/(1.+r);
}

/*********************************************************/

phydbl Rescale_Free_Rate_Tree(t_tree *tree)
{
  int i;
  phydbl sum;

  sum = 0.0;
  For(i,tree->mod->n_catg) sum += tree->mod->gamma_rr[i]*tree->mod->gamma_r_proba[i];
  For(i,tree->mod->n_catg) tree->mod->gamma_rr[i] /= sum;
  
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l *= sum;
  
  return(-1.);
}

/*********************************************************/

phydbl Rescale_Br_Len_Multiplier_Tree(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l *= tree->mod->br_len_multiplier;
  return(-1.);
}

/*********************************************************/

phydbl Unscale_Free_Rate_Tree(t_tree *tree)
{
  int i;
  phydbl sum;

  sum = 0.0;
  For(i,tree->mod->n_catg) sum += tree->mod->gamma_rr[i]*tree->mod->gamma_r_proba[i];
  For(i,tree->mod->n_catg) tree->mod->gamma_rr[i] *= sum;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l /= sum;
  
  return(-1.);
}

/*********************************************************/

phydbl Unscale_Br_Len_Multiplier_Tree(t_tree *tree)
{
  int i;
  For(i,2*tree->n_otu-3) tree->t_edges[i]->l /= tree->mod->br_len_multiplier;
  return(-1.);
}

/*********************************************************/

phydbl Reflect(phydbl x, phydbl l, phydbl u)
{
  int rounds;
  phydbl tmp;
  int k;

  if(u < l)
    {
      tmp = u;
      u   = l;
      l   = tmp;
    }

  if(x < l) x = x + 2.*(l - x);
  
  if(((x-u) > (u-l)) && (x > u))
    {
      k = (x - (2.*u-l))/(2.*(u-l));
      x = x - 2.*k*(u-l);
    }

  rounds = 0;
  do
    {
      rounds++;
      /* printf("\n. l=%f u=%f x=%f",l,u,x); */
      if(x > u || x < l) 
	{
	  if(x > u) x = x - 2.*(x - u);	  
	  else      x = x + 2.*(l - x);
	}
      else break;
      /* printf(" x'=%f",x); */
    }
  while(rounds < 100);
    
  if(rounds == 100 && (x > u || x < l))
    {
      PhyML_Printf("\n. u=%f l=%f x=%f",u,l,x);
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Exit("\n");
    }

  return x;
}

/*********************************************************/

int Are_Equal(phydbl a, phydbl b, phydbl eps)
{
  if(FABS(a-b) < eps) return TRUE; /* a==b */
  else return FALSE;
}


/*********************************************************/
