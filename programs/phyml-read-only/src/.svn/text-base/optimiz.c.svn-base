/*

PHYML :  a program that  computes maximum likelihood  phyLOGenies from
DNA or AA homoLOGous sequences 

Copyright (C) Stephane Guindon. Oct 2003 onward

All parts of  the source except where indicated  are distributed under
the GNU public licence.  See http://www.opensource.org for details.

*/

#include "optimiz.h"


/*********************************************************/

void Optimize_Single_Param_Generic(t_tree *tree, phydbl *param, phydbl lim_inf, phydbl lim_sup, phydbl tol, int n_max_iter, int quickdirty)
{
  phydbl ax,bx,cx;
  phydbl lk_init;
  
  lk_init = tree->c_lnL;

  ax =  lim_inf;
  bx = (*param);
  cx =  lim_sup;
  
  Generic_Brent(ax,bx,cx,tol,param,tree,n_max_iter,quickdirty);

  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_global) 
    {
      PhyML_Printf("\n. %.10f < %.10f --> diff=%.10f param value = %f initial value = %f\n",
	     tree->c_lnL,lk_init,
	     tree->c_lnL-lk_init,
	     *param,bx);
      Exit("\n. Optimisation failed !\n");
    }
}

/*********************************************************/

int Generic_Brak(phydbl *param,
		 phydbl *ax, phydbl *bx, phydbl *cx, 
		 phydbl *fa, phydbl *fb, phydbl *fc,
		 phydbl lim_inf, phydbl lim_sup,
		 t_tree *tree)
{
   phydbl ulim,u,r,q,fu,dum;

   u = 0.0;
   *param = *ax;

   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fa=-Lk(tree);
   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fb=-Lk(tree);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = FABS(*cx);
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fc=-Lk(tree); 
   while (*fb > *fc) 
     {
        
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;
       if(u   > lim_sup) u   = lim_sup;
       if(u   < lim_inf) u   = lim_inf;

       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > lim_inf) 
	 {
	   *param = FABS(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Lk(tree);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
	       (*ax)=FABS(*ax);
	       (*bx)=FABS(*bx);
	       (*cx)=FABS(*cx);
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;	
	       (*ax)=FABS(*ax);
	       (*bx)=FABS(*bx);
	       (*cx)=FABS(*cx);
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = FABS(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Lk(tree);
	 } 
       else if ((*cx-u)*(u-ulim) > lim_inf) 
	 {
	   *param = FABS(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Lk(tree);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       *param = FABS(u); 
	       SHFT(*fb,*fc,fu,-Lk(tree))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= lim_inf) 
	 {
	   u=ulim;
	   *param = FABS(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Lk(tree);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = FABS(u);
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu=-Lk(tree);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)


     }
   (*ax)=FABS(*ax);
   (*bx)=FABS(*bx);
   (*cx)=FABS(*cx);
   return(0);
}

/*********************************************************/

phydbl Generic_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		     phydbl *xmin, t_tree *tree, int n_iter_max, 
		     int quickdirty)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL;

  
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*xmin) = FABS(bx);
  fw=fv=fx=-Lk(tree);
  init_lnL = -fw;

/*   PhyML_Printf("\n. init_lnL = %f a=%f b=%f c=%f\n",init_lnL,ax,bx,cx); */

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);


      if(FABS(x - xm) <= (tol2 - 0.5 * (b - a)))
	{
	  *xmin = x;
	  Lk(tree);
	  if(tree->c_lnL < init_lnL - tree->mod->s_opt->min_diff_lk_local)
	    {
	      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	      Warn_And_Exit("");
	    }
	  return tree->c_lnL;
	}
      
      if(FABS(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
/* 	      PhyML_Printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
/* 	      PhyML_Printf("Parabolic step [e=%f]\n",e); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
/* 	  PhyML_Printf("Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f]\n",e,tol1,a,b,d); */
	}
      
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*xmin) = FABS(u);
      old_lnL = tree->c_lnL;
      fu = -Lk(tree);
      
/*       PhyML_Printf("\n. iter=%d/%d param=%f LOGlk=%f",iter,BRENT_ITMAX,*xmin,tree->c_lnL); */

/*       if(fu <= fx) */
      if(fu < fx)
	{
/* 	  if(u >= x) a=x; else b=x; */
	  if(u > x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x) */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
/* 	  else if (fu <= fv || v == x || v == w) */
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }

  Exit("\n. Too many iterations in BRENT !");
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

phydbl RRparam_GTR_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
			  phydbl *xmin, t_tree *tree, calign *cdata, phydbl *param, int n_iter_max)
{
   phydbl f1,f2,x0,x1,x2,x3;
   int n_iter;


   x0=ax;
   x3=cx;
   if (FABS(cx-bx) > FABS(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   (*param)=x1;

   Lk(tree);
   f1=-tree->c_lnL;
   (*param)=x2;

   Lk(tree);
   f2=-tree->c_lnL;

   n_iter = 0;
   while (FABS(x3-x0) > tol*(FABS(x1)+FABS(x2))) 
     {

       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   (*param)=x2;
	   Lk(tree);
	   SHFT2(f1,f2,-tree->c_lnL)
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   (*param)=x1;
	   Lk(tree);
	   SHFT2(f2,f1,-tree->c_lnL)
	 }
       
       if(n_iter++ > n_iter_max) break;
       
/*        PhyML_Printf("p=%E %f\n",(*param),tree->c_lnL); */
     }
   if (f1 < f2) 
    {
       *xmin=x1;
       return f1;
     } 
   else 
     {
       *xmin=x2;
       return f2;
     }
}

/*********************************************************/

phydbl Br_Len_Golden(phydbl ax, phydbl bx, phydbl cx, phydbl tol, 
		     phydbl *xmin, t_edge *b_fcus, t_tree *tree)
{
   phydbl f1,f2,x0,x1,x2,x3;

   x0=ax;
   x3=cx;
   if (FABS(cx-bx) > FABS(bx-ax)) 
     {
       x1=bx;
       x2=bx+GOLDEN_C*(cx-bx);
     } 
   else 
     {
       x2=bx;
       x1=bx-GOLDEN_C*(bx-ax);
     }
   
   b_fcus->l=x1;
   f1 = -Lk_At_Given_Edge(b_fcus,tree);
   b_fcus->l=x2;
   f2 = -Lk_At_Given_Edge(b_fcus,tree);
   while (FABS(x3-x0) > tol*(FABS(x1)+FABS(x2))) 
     {
       if (f2 < f1) 
	 {
	   SHFT3(x0,x1,x2,GOLDEN_R*x1+GOLDEN_C*x3)
	   b_fcus->l=x2;
	   SHFT2(f1,f2,-Lk_At_Given_Edge(b_fcus,tree))
	 } 
       else 
	 {
	   SHFT3(x3,x2,x1,GOLDEN_R*x2+GOLDEN_C*x0)
	   b_fcus->l=x1;
	   SHFT2(f2,f1,-Lk_At_Given_Edge(b_fcus,tree))
	 }
     }
   if (f1 < f2) 
     {
       *xmin=FABS(x1);
       return -f1;
     } 
   else 
     {
       *xmin=FABS(x2);
       return -f2;
     }
}

/*********************************************************/

int Br_Len_Brak(phydbl *ax, phydbl *bx, phydbl *cx, 
		phydbl *fa, phydbl *fb, phydbl *fc, 
		t_edge *b_fcus, t_tree *tree)
{
   phydbl ulim,u,r,q,fu,dum;

   b_fcus->l = *ax;
   *fa=-Lk_At_Given_Edge(b_fcus,tree);
   b_fcus->l = *bx;
   *fb=-Lk_At_Given_Edge(b_fcus,tree);
   if (*fb > *fa) {
      SHFT(dum,*ax,*bx,dum)
      SHFT(dum,*fb,*fa,dum)
   }
   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   b_fcus->l = *cx;
   *fc=-Lk_At_Given_Edge(b_fcus,tree);
   while (*fb > *fc + tree->mod->s_opt->min_diff_lk_local) 
     {
       PhyML_Printf("fb=%f fc=%f\n",*fb,*fc);
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;
/* 	       (*ax)=FABS(*ax); */
/* 	       (*bx)=FABS(*bx); */
/* 	       (*cx)=FABS(*cx); */
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;	
/* 	       (*ax)=FABS(*ax); */
/* 	       (*bx)=FABS(*bx); */
/* 	       (*cx)=FABS(*cx); */
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   b_fcus->l = FABS(u);
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       b_fcus->l = u; 
	       SHFT(*fb,*fc,fu,-Lk_At_Given_Edge(b_fcus,tree))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u=ulim;
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   b_fcus->l = u;
	   fu=-Lk_At_Given_Edge(b_fcus,tree);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)
      }
   (*ax)=FABS(*ax);
   (*bx)=FABS(*bx);
   (*cx)=FABS(*cx);
   return(0);
}

/*********************************************************/

phydbl Br_Len_Brent_Default(t_edge *b_fcus, t_tree *tree)
{
  return Br_Len_Brent(10.*b_fcus->l,b_fcus->l,.10*b_fcus->l,tree->mod->s_opt->min_diff_lk_local,b_fcus,tree,1000,0);
}

/*********************************************************/

phydbl Br_Len_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol,
		    t_edge *b_fcus, t_tree *tree, int n_iter_max, int quickdirty)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL, init_lnL;
  
  if(tree->mod->gamma_mgf_bl == YES)
    {
      Generic_Brent_Lk(&(b_fcus->gamma_prior_var),
      		       0.0,1000.,
      		       tree->mod->s_opt->min_diff_lk_local,
      		       tree->mod->s_opt->brent_it_max,
      		       tree->mod->s_opt->quickdirty,
      		       Wrap_Lk_At_Given_Edge,b_fcus,tree,NULL);

    }
      
  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  b_fcus->l = bx;
  fw=fv=fx=fu=-Lk_At_Given_Edge(b_fcus,tree);
  init_lnL = -fw;
  
  /*   PhyML_Printf("\n. INIT BRENT t_edge %3d l=%f lnL=%20f",b_fcus->num,b_fcus->l,fu); */
  
  for(iter=1;iter<=BRENT_ITMAX;iter++)
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*x+BRENT_ZEPS);
      
      if((tree->c_lnL > init_lnL + tol) && (quickdirty))
	{
	  b_fcus->l = x;
	  Lk_At_Given_Edge(b_fcus,tree);
	  /* 	  PhyML_Printf("\n> iter=%3d max=%3d v=%f lnL=%f init_lnL=%f tol=%f",iter,n_iter_max,(*xmin),tree->c_lnL,init_lnL,tol); */
	  /* 	  Exit("\n"); */
	  return tree->c_lnL;	  
	}
      
      /*       if(((FABS(tree->c_lnL-old_lnL) < tol) && (tree->c_lnL > init_lnL - tol)) || (iter > n_iter_max - 1)) */
      if((FABS(tree->c_lnL-old_lnL) < tol) || (iter > n_iter_max - 1))
	{
	  b_fcus->l=x;
	  Lk_At_Given_Edge(b_fcus,tree);
	  /* 	  PhyML_Printf("\n. iter=%3d max=%3d l=%f lnL=%f init_lnL=%f",iter,n_iter_max,b_fcus->l,tree->c_lnL,init_lnL); */
	  /* 	  Exit("\n"); */
	  return tree->c_lnL;
	}
      
      if(FABS(e) > tol1)
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  else{
	    d=p/q;
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  }
	}
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	}
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      /*       if(u<tree->l_min) u = tree->l_min; */
      b_fcus->l=u;
      old_lnL = tree->c_lnL;
      fu=-Lk_At_Given_Edge(b_fcus,tree);
      
      /*       PhyML_Printf("\n. BRENT t_edge %3d l=%f lnL=%20f iter=%3d",b_fcus->num,b_fcus->l,fu,iter); */
      
      /*       if(fu <= fx) */
      if(fu < fx)
	{
	  if(u > x) a=x; else b=x;
	  /* 	  if(u >= x) a=x; else b=x; */
	  SHFT(v,w,x,u)
	    SHFT(fv,fw,fx,fu)
	    }
      else
	{
	  if (u < x) a=u; else b=u;
	  /* 	  if (fu <= fw || w == x) */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    }
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    /* 	  else if (fu <= fv || v == x || v == w)  */
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  if(iter > BRENT_ITMAX) PhyML_Printf("\n. Too many iterations in BRENT (%d) (%f)",iter,b_fcus->l);      
  return(-1);
  /* Not Reached ??  *xmin=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

void Round_Optimize(t_tree *tree, calign *data, int n_round_max)
{
  int n_round,each;
  phydbl lk_old, lk_new, tol;
  t_node *root;

  lk_new = tree->c_lnL;
  lk_old = UNLIKELY;
  n_round = 0;
  each = 0;
  tol = 1.e-2;
  root = tree->noeud[0];
  
  tree->both_sides = 1;
  Lk(tree);

  while(n_round < n_round_max)
    {
      (!((n_round+2)%2))?(root=tree->noeud[0]):(root=tree->noeud[tree->n_otu-1]);
      
      if((tree->mod->s_opt->opt_bl) && 
	 (tree->mod->s_opt->print) && 
	 (!tree->io->quiet)) 
	{
	  Print_Lk(tree,"[Branch lengths     ]");
	  Optimize_Br_Len_Serie(root,root->v[0],root->b[0],tree,data);
	}

      tree->both_sides = 1;
      Lk(tree);

      if(!each)
	{
	  each = 1;
	  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
	  tree->both_sides = 1;
	  Lk(tree);
	}
      
      lk_new = tree->c_lnL;      
      if(lk_new < lk_old - tree->mod->s_opt->min_diff_lk_global) 
	{
	  PhyML_Printf("\n. lk_new = %f lk_old = %f diff = %f",lk_new,lk_old,lk_new-lk_old);
	  Exit("\n. Optimisation failed ! (Round_Optimize)\n");
	}
      if(FABS(lk_new - lk_old) < tree->mod->s_opt->min_diff_lk_global)  break;
      else lk_old  = lk_new;
      n_round++;
      each--;
    }
  
  Optimiz_All_Free_Param(tree,(tree->io->quiet)?(0):(tree->mod->s_opt->print));
}

/*********************************************************/

void Optimize_Br_Len_Serie(t_node *a, t_node *d, t_edge *b_fcus, t_tree *tree, calign *cdata)
{
  int i;
  phydbl l_infa,l_max,l_infb;
  phydbl lk_init;
  
  lk_init = tree->c_lnL;

  if(tree->mod->s_opt->constrained_br_len == YES)
    {
      Generic_Brent_Lk(&(tree->mod->br_len_multiplier),
		       1.E-2,1.E+2,
		       tree->mod->s_opt->min_diff_lk_global,
		       tree->mod->s_opt->brent_it_max,
		       tree->mod->s_opt->quickdirty,
		       Wrap_Lk,NULL,tree,NULL);
      
      
      if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
	{
	  PhyML_Printf("\n. %f %f %f %f",l_infa,l_max,l_infb,b_fcus->l);
	  PhyML_Printf("\n. %f -- %f",lk_init,tree->c_lnL);
	  Warn_And_Exit("\n. Err. in Optimize_Br_Len_Serie\n");
	}

      return;
    }
  
   
  l_max  = b_fcus->l;
  l_infa = tree->mod->l_max;
  l_infb = tree->mod->l_min;
  
  Br_Len_Brent(l_infa,l_max,l_infb,
	       tree->mod->s_opt->min_diff_lk_local,
	       b_fcus,tree,
	       tree->mod->s_opt->brent_it_max,
	       tree->mod->s_opt->quickdirty);

  if(tree->mod->s_opt->opt_gamma_br_len == YES)
    {      
      Generic_Brent_Lk(&(b_fcus->gamma_prior_var),
		       0.0,1000.,
		       tree->mod->s_opt->min_diff_lk_local,
		       tree->mod->s_opt->brent_it_max,
		       tree->mod->s_opt->quickdirty,
		       Wrap_Lk_At_Given_Edge,b_fcus,tree,NULL);
    }

/*   Generic_Brent_Lk(&(b_fcus->l), */
/* 		   l_infa,l_infb, */
/* 		   tree->mod->s_opt->min_diff_lk_local, */
/* 		   tree->mod->s_opt->brent_it_max, */
/* 		   tree->mod->s_opt->quickdirty, */
/* 		   Wrap_Lk_At_Given_Edge,b_fcus,tree,NULL); */
  

  if(tree->c_lnL < lk_init - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. %f %f %f %f",l_infa,l_max,l_infb,b_fcus->l);
      PhyML_Printf("\n. %f -- %f",lk_init,tree->c_lnL);
      Warn_And_Exit("\n. Err. in Optimize_Br_Len_Serie\n");
    }
    
  if(d->tax) return;
  else For(i,3) if(d->v[i] != a)
    {
      Update_P_Lk(tree,d->b[i],d);
      Optimize_Br_Len_Serie(d,d->v[i],d->b[i],tree,cdata);
    }
  For(i,3) if((d->v[i] == a) && !(d->v[i]->tax)) Update_P_Lk(tree,d->b[i],d);
}

/*********************************************************/

void Optimiz_Ext_Br(t_tree *tree)
{
  int i;
  t_edge *b;
  phydbl l_infa,l_max,l_infb;
  phydbl lk, lk_init,l_init;
  
  lk_init = tree->c_lnL;
  

  For(i,2*tree->n_otu-3)
    {
      b = tree->t_edges[i];
      if((b->left->tax) || (b->rght->tax))
	{

	  l_init = b->l;

/* 	  Fast_Br_Len(b,tree); */
/* 	  lk = Lk_At_Given_Edge(tree,b); */

	  l_infa = 10.*b->l;
	  l_max  = b->l;
	  l_infb = tree->mod->l_min;

	  lk = Br_Len_Brent(l_infa,l_max,l_infb,
			    tree->mod->s_opt->min_diff_lk_local,
			    b,tree,
			    tree->mod->s_opt->brent_it_max,
			    tree->mod->s_opt->quickdirty);

	  b->nni->best_l    = b->l;
	  b->nni->l0        = b->l;
	  b->nni->best_conf = 0;
	  b->l              = l_init;

	}
    }
  tree->c_lnL = lk_init; 
}

/*********************************************************/

void Optimiz_All_Free_Param(t_tree *tree, int verbose)
{
  int  init_both_sides;

  init_both_sides  = tree->both_sides;
  tree->both_sides = 0;



  if((tree->mod->whichmodel == GTR) ||
     ((tree->mod->whichmodel == CUSTOM) && 
      (tree->mod->s_opt->opt_rr) && 
      (tree->mod->n_diff_rr > 1)))
    {
      int failed;
      
      failed = 0;
      
      tree->mod->update_eigen = 1;
/*       BFGS(tree,tree->mod->rr_val,tree->mod->n_diff_rr,1.e-5,1.e-5, */
      BFGS(tree,tree->mod->rr_val,tree->mod->n_diff_rr,1.e-5,1.e-3,
	   &Return_Abs_Lk,
	   &Num_Derivative_Several_Param,
	   &Lnsrch_RR_Param,&failed);
      
      if(failed)
	{
	  int i;
	  
	  For(i,tree->mod->n_diff_rr)
	    if(i != 5)
	      {
/* 		Optimize_Single_Param_Generic(tree,&(tree->mod->rr_val[i]), */
/* 					      1.E-2,1.E+2, */
/* 					      tree->mod->s_opt->min_diff_lk_global, */
/* 					      tree->mod->s_opt->brent_it_max, */
/* 					      tree->mod->s_opt->quickdirty); */

		Generic_Brent_Lk(&(tree->mod->br_len_multiplier),
				 1.E-2,1.E+2,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Wrap_Lk,NULL,tree,NULL);
		
	      }

	}
      if(verbose) Print_Lk(tree,"[GTR parameters     ]");
      tree->mod->update_eigen = 0;
    }
  
  if(tree->mod->s_opt->opt_kappa)
    {
      tree->mod->update_eigen = 1;
/*       Optimize_Single_Param_Generic(tree,&(tree->mod->kappa), */
/* 				    .1,100., */
/* 				    tree->mod->s_opt->min_diff_lk_global, */
/* 				    tree->mod->s_opt->brent_it_max, */
/* 				    tree->mod->s_opt->quickdirty); */

      Generic_Brent_Lk(&(tree->mod->kappa),
		       0.1,100.,
		       tree->mod->s_opt->min_diff_lk_global,
		       tree->mod->s_opt->brent_it_max,
		       tree->mod->s_opt->quickdirty,
		       Wrap_Lk,NULL,tree,NULL);
      
      if(verbose) 
	{
	  Print_Lk(tree,"[Ts/ts ratio        ]");
	  PhyML_Printf("[%10f]",tree->mod->kappa);
	}
      tree->mod->update_eigen = 0;
    }

  if(tree->mod->s_opt->opt_lambda) 
    {
/*       Optimize_Single_Param_Generic(tree,&(tree->mod->lambda),.001,100., */
/* 				    tree->mod->s_opt->min_diff_lk_global, */
/* 				    tree->mod->s_opt->brent_it_max, */
/* 				    tree->mod->s_opt->quickdirty); */
      Generic_Brent_Lk(&(tree->mod->lambda),
		       0.001,100.,
		       tree->mod->s_opt->min_diff_lk_global,
		       tree->mod->s_opt->brent_it_max,
		       tree->mod->s_opt->quickdirty,
		       Wrap_Lk,NULL,tree,NULL);

      if(verbose) 
	{
	  Print_Lk(tree,"[Lambda             ]");
	  PhyML_Printf("[%10f]",tree->mod->lambda);
	}
    }
  
  if((tree->mod->s_opt->opt_pinvar) && (tree->mod->s_opt->opt_alpha))
    {
      Optimiz_Alpha_And_Pinv(tree);
      if(verbose) 
	{
	  Print_Lk(tree,"[Alpha              ]"); 
	  PhyML_Printf("[%10f]",tree->mod->alpha);
	  Print_Lk(tree,"[P-inv              ]"); 
	  PhyML_Printf("[%10f]",tree->mod->pinvar);
	}
    }
  else
    {
      if(tree->mod->s_opt->opt_pinvar)
	{
/*  	  Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar), */
/* 					.0001,0.9999, */
/* 					tree->mod->s_opt->min_diff_lk_global, */
/* 					tree->mod->s_opt->brent_it_max, */
/* 					tree->mod->s_opt->quickdirty); */

	  Generic_Brent_Lk(&(tree->mod->pinvar),
			   0.0001,0.9999,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Wrap_Lk,NULL,tree,NULL);

	  if(verbose) 
	    {
	      Print_Lk(tree,"[P-inv              ]");
	      PhyML_Printf("[%10f]",tree->mod->pinvar);
	    }
	}
      
      if(tree->mod->s_opt->opt_alpha && tree->mod->free_mixt_rates == NO)
	{
	  if(tree->mod->n_catg > 1)
/* 	    Optimize_Single_Param_Generic(tree,&(tree->mod->alpha), */
/* /\* 					  .01,100., *\/ */
/* 					  tree->mod->alpha/2.,tree->mod->alpha*2., */
/* 					  tree->mod->s_opt->min_diff_lk_global, */
/* 					  tree->mod->s_opt->brent_it_max, */
/* 					  tree->mod->s_opt->quickdirty); */

	    Generic_Brent_Lk(&(tree->mod->alpha),
			     tree->mod->alpha/2.,tree->mod->alpha*2.,
			     tree->mod->s_opt->min_diff_lk_global,
			     tree->mod->s_opt->brent_it_max,
			     tree->mod->s_opt->quickdirty,
			     Wrap_Lk,NULL,tree,NULL);
	  if(verbose) 
	    {
	      Print_Lk(tree,"[Alpha              ]");
	      PhyML_Printf("[%10f]",tree->mod->alpha);
	    }
	}
    }

  if((tree->mod->s_opt->opt_state_freq) && (tree->io->datatype == NT))
    {
        int failed,i;
        
        failed = 0;
        tree->mod->update_eigen = 1;
        BFGS(tree,tree->mod->pi_unscaled,tree->mod->ns,1.e-5,1.e-5,
	     &Return_Abs_Lk,
	     &Num_Derivative_Several_Param,
	     &Lnsrch_Nucleotide_Frequencies,&failed);

        if(failed)
	  {
	    For(i,tree->mod->ns) 
	      {
/* 		Optimize_Single_Param_Generic(tree,&(tree->mod->pi_unscaled[i]), */
/* 					      -1000.,1000., */
/* 					      tree->mod->s_opt->min_diff_lk_global, */
/* 					      tree->mod->s_opt->brent_it_max, */
/* 					      tree->mod->s_opt->quickdirty); */
 		Generic_Brent_Lk(&(tree->mod->pi_unscaled[i]),
				 -1000.,1000.,
				 tree->mod->s_opt->min_diff_lk_global,
				 tree->mod->s_opt->brent_it_max,
				 tree->mod->s_opt->quickdirty,
				 Wrap_Lk,NULL,tree,NULL);
	      }
	  }
	if(verbose) Print_Lk(tree,"[Nucleotide freqs.  ]");
        tree->mod->update_eigen = 0;
    }



  if((tree->mod->s_opt->opt_free_mixt_rates) && (tree->mod->free_mixt_rates == YES))
    {
      /* int failed; */
      int i;
        
      /* failed = 0; */
      /* tree->mod->update_eigen = 1; */
      /* BFGS(tree,tree->mod->gamma_r_proba_unscaled,tree->mod->n_catg,1.e-5,1.e-5, */
      /* 	   &Return_Abs_Lk, */
      /* 	   &Num_Derivative_Several_Param, */
      /* 	   &Lnsrch_Free_Mixt_Rates,&failed); */


      /* if(failed) */
      /* 	{ */
      	  For(i,tree->mod->n_catg)
      	    {
      	      Generic_Brent_Lk(&(tree->mod->gamma_r_proba_unscaled[i]),
      			       -1000.,1000.,
      			       tree->mod->s_opt->min_diff_lk_global,
      			       tree->mod->s_opt->brent_it_max,
      			       tree->mod->s_opt->quickdirty,
      			       Wrap_Lk,NULL,tree,NULL);
      	    }
      	/* } */
      
      
      if(verbose) Print_Lk(tree,"[Rate class freqs.  ]");

      For(i,tree->mod->n_catg) 
	{
	  Generic_Brent_Lk(&(tree->mod->gamma_rr_unscaled[i]),
			   -1000.,1000.,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Wrap_Lk,NULL,tree,NULL);
	}

      if(verbose) Print_Lk(tree,"[Rate class values  ]");

      tree->mod->update_eigen = 0;
    }



  if(tree->mod->use_m4mod)
    {
      int failed,i;

      if(tree->mod->s_opt->opt_cov_delta) 
	{
	  tree->mod->update_eigen = 1;
/* 	  Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->delta), */
/* 					0.01,10., */
/* 					tree->mod->s_opt->min_diff_lk_global, */
/* 					tree->mod->s_opt->brent_it_max, */
/* 					tree->mod->s_opt->quickdirty); */

	  Generic_Brent_Lk(&(tree->mod->m4mod->delta),
			   0.01,10.,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Wrap_Lk,NULL,tree,NULL);

	  if(verbose) 
	    {
	      Print_Lk(tree,"[Switching param.   ]");
	      PhyML_Printf("[%10f]",tree->mod->m4mod->delta);
	    }
	  tree->mod->update_eigen = 0;
	}
      
      if(tree->mod->s_opt->opt_cov_free_rates) 
	{
	  int rcat;
	  

	  tree->mod->update_eigen = 1;
	  For(rcat,tree->mod->m4mod->n_h)
	    {
/* 	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->multipl_unscaled[rcat]), */
/* 					    .01,10., */
/* 					    tree->mod->s_opt->min_diff_lk_global, */
/* 					    tree->mod->s_opt->brent_it_max, */
/* 					    tree->mod->s_opt->quickdirty); */

	      Generic_Brent_Lk(&(tree->mod->m4mod->multipl_unscaled[rcat]),
			       0.01,10.,
			       tree->mod->s_opt->min_diff_lk_global,
			       tree->mod->s_opt->brent_it_max,
			       tree->mod->s_opt->quickdirty,
			       Wrap_Lk,NULL,tree,NULL);
	      
	      if(verbose) 
		{
		  Print_Lk(tree,"[Rel. subst. rate   ]");
		  PhyML_Printf("[%10f]",tree->mod->m4mod->multipl[rcat]);
		}
	    }
	  
	  For(rcat,tree->mod->m4mod->n_h)
	    {
/*  	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->h_fq_unscaled[rcat]), */
/* 					    .01,100., */
/* 					    tree->mod->s_opt->min_diff_lk_global, */
/* 					    tree->mod->s_opt->brent_it_max, */
/* 					    tree->mod->s_opt->quickdirty); */

	      Generic_Brent_Lk(&(tree->mod->m4mod->h_fq_unscaled[rcat]),
			       0.01,100.,
			       tree->mod->s_opt->min_diff_lk_global,
			       tree->mod->s_opt->brent_it_max,
			       tree->mod->s_opt->quickdirty,
			       Wrap_Lk,NULL,tree,NULL);

	      
	      if(verbose)
		{
		  Print_Lk(tree,"[Subst. class freq  ]");
		  PhyML_Printf("[%10f]",tree->mod->m4mod->h_fq[rcat]);
		}
	    }
	  tree->mod->update_eigen = 0;      
	}
      
      if(tree->mod->s_opt->opt_cov_alpha) 
	{
	  tree->mod->update_eigen = 1;
/* 	  Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->alpha), */
/* 					.01,10., */
/* 					tree->mod->s_opt->min_diff_lk_global, */
/* 					tree->mod->s_opt->brent_it_max, */
/* 					tree->mod->s_opt->quickdirty); */

	  Generic_Brent_Lk(&(tree->mod->m4mod->alpha),
			   0.01,10.,
			   tree->mod->s_opt->min_diff_lk_global,
			   tree->mod->s_opt->brent_it_max,
			   tree->mod->s_opt->quickdirty,
			   Wrap_Lk,NULL,tree,NULL);


	  if(verbose) 
	    {
	      Print_Lk(tree,"[Alpha (covarion)   ]");
	      PhyML_Printf("[%10f]",tree->mod->m4mod->alpha);
	    }
	  tree->mod->update_eigen = 0;      	  
	}
      
	  
      /* Substitutions between nucleotides are considered to follow a 
	 GTR model */

      failed = 0;
      tree->mod->update_eigen = 1;
      BFGS(tree,tree->mod->m4mod->o_rr,5,1.e-5,1.e-5,
	   &Return_Abs_Lk,
	   &Num_Derivative_Several_Param,
	   &Lnsrch_RR_Cov_Param,&failed);
      
      if(failed)
	{
	  For(i,5) 
	    {
/* 	      Optimize_Single_Param_Generic(tree,&(tree->mod->m4mod->o_rr[i]), */
/* 					    1.E-20,1.E+10, */
/* 					    tree->mod->s_opt->min_diff_lk_global, */
/* 					    tree->mod->s_opt->brent_it_max, */
/* 					    tree->mod->s_opt->quickdirty); */

	      Generic_Brent_Lk(&(tree->mod->m4mod->o_rr[i]),
			       1.E-20,1.E+10,
			       tree->mod->s_opt->min_diff_lk_global,
			       tree->mod->s_opt->brent_it_max,
			       tree->mod->s_opt->quickdirty,
			       Wrap_Lk,NULL,tree,NULL);
	      
	    }
	}
      if(verbose) Print_Lk(tree,"[GTR parameters     ]");
      tree->mod->update_eigen = 0;
    }

  tree->both_sides = init_both_sides;

  if(tree->both_sides) Lk(tree); /* Needed to update all partial likelihoods */
}



#define ITMAX 200
#define EPS   3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
static phydbl sqrarg;
#define SQR(a) ((sqrarg=(a)) < SMALL ? 0.0 : sqrarg*sqrarg)

void BFGS(t_tree *tree, 
	  phydbl *p, 
	  int n, 
	  phydbl gtol, 
	  phydbl step_size,
	  phydbl(*func)(t_tree *tree), 
	  int(*dfunc)(t_tree *tree,phydbl *param,int n_param,phydbl stepsize,phydbl(*func)(t_tree *tree),phydbl *derivatives), 
	  int(*lnsrch)(t_tree *tree, int n, phydbl *xold, phydbl fold,phydbl *g, phydbl *p, phydbl *x,phydbl *f, phydbl stpmax, int *check),
	  int *failed)
{

  int check,i,its,j;
  phydbl den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test,fret;
  phydbl *dg,*g,*hdg,**hessin,*pnew,*xi;
  
  hessin = (phydbl **)mCalloc(n,sizeof(phydbl *));
  For(i,n) hessin[i] = (phydbl *)mCalloc(n,sizeof(phydbl));
  dg   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  g    = (phydbl *)mCalloc(n,sizeof(phydbl ));
  pnew = (phydbl *)mCalloc(n,sizeof(phydbl ));
  hdg  = (phydbl *)mCalloc(n,sizeof(phydbl ));
  xi   = (phydbl *)mCalloc(n,sizeof(phydbl ));
  

/*   PhyML_Printf("\n. ENTER BFGS WITH: %f\n",Lk(tree)); */

  fp=(*func)(tree);
  (*dfunc)(tree,p,n,step_size,func,g);

  for (i=0;i<n;i++) 
    {
      for (j=0;j<n;j++) hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }

  stpmax=STPMX*MAX(SQRT(sum),(phydbl)n);

  for(its=1;its<=ITMAX;its++) 
    {
      lnsrch(tree,n,p,fp,g,xi,pnew,&fret,stpmax,&check);

/*       PhyML_Printf("BFGS -> %f\n",tree->c_lnL); */

      fp = fret;
      
      for (i=0;i<n;i++) 
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}

      test=0.0;
      for (i=0;i<n;i++) 
	{
	  temp=FABS(xi[i])/MAX(FABS(p[i]),1.0);
	  if (temp > test) test=temp;
	}
      if (test < TOLX) 
	{
	  (*func)(tree);
	  For(i,n) Free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   

	  if(its == 1) 
	    {
/* 	      PhyML_Printf("\n. WARNING : BFGS failed ! \n"); */
	      *failed = 1;
	    }
	  return;
	}

      for (i=0;i<n;i++) dg[i]=g[i];

      (*dfunc)(tree,p,n,step_size,func,g);

      test=0.0;
      den=MAX(fret,1.0);
      for (i=0;i<n;i++) 
	{
	  temp=FABS(g[i])*MAX(FABS(p[i]),1.0)/den;
	  if (temp > test) test=temp;
	}
      if (test < gtol) 
	{
	  (*func)(tree);
	  For(i,n) Free(hessin[i]);
	  free(hessin);
	  free(xi);
	  free(pnew);
	  free(hdg);
	  free(g);
	  free(dg);   
	  return;
	}

    for (i=0;i<n;i++) dg[i]=g[i]-dg[i];

    for (i=0;i<n;i++) 
      {
	hdg[i]=0.0;
	for (j=0;j<n;j++) hdg[i] += hessin[i][j]*dg[j];
      }

    fac=fae=sumdg=sumxi=0.0;
    for (i=0;i<n;i++) 
      {
	fac += dg[i]*xi[i];
	fae += dg[i]*hdg[i];
	sumdg += SQR(dg[i]);
	sumxi += SQR(xi[i]);
      }
    
    if(fac*fac > EPS*sumdg*sumxi) 
      {
	fac=1.0/fac;
	fad=1.0/fae;
	for (i=0;i<n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
	for (i=0;i<n;i++) 
	  {
	    for (j=0;j<n;j++) 
	      {
		hessin[i][j] += fac*xi[i]*xi[j]
		  -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
	      }
	  }
      }
    for (i=0;i<n;i++) 
      {
	xi[i]=0.0;
	for (j=0;j<n;j++) xi[i] -= hessin[i][j]*g[j];
      }
    }
/*   PhyML_Printf("\n. Too many iterations in BFGS...\n"); */
  *failed = 1;
  For(i,n) Free(hessin[i]);
  free(hessin);
  free(xi);
  free(pnew);
  free(hdg);
  free(g);
  free(dg);   
}

#undef ITMAX
#undef EPS
#undef TOLX
#undef STPMX

/*********************************************************/


#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch_RR_Param(t_tree *tree, int n, phydbl *xold, phydbl fold, 
		    phydbl *g, phydbl *p, phydbl *x,
		    phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];

  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=FABS(p[i])/MAX(FABS(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) 
	{
	  x[i]=FABS(local_xold[i]+alam*p[i]);
	}

      /**/
      for(i=0;i<n;i++)
	{
	  tree->mod->rr_val[i] = FABS(local_xold[i]+alam*p[i]);
	}
      /**/

      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
	}
      else *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) 
	    {
	      x[i]=FABS(local_xold[i]);
	    }
	  /**/      
	  for(i=0;i<n;i++)
	    {
	      tree->mod->rr_val[i] = FABS(local_xold[i]);
	    }
	  /**/

	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return 0;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return 0;
	}
      else 
	{
/* 	  if (alam == 1.0) */
	  if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a < SMALL && a > -SMALL) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) tmplam = 0.5*alam;
		  else if(b <= 0.0) tmplam=(-b+SQRT(disc))/(3.0*a);
		  else tmplam = -slope/(b+SQRT(disc));
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
  return 1;
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch_RR_Cov_Param(t_tree *tree, int n, phydbl *xold, phydbl fold, 
			phydbl *g, phydbl *p, phydbl *x,
			phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=FABS(p[i])/MAX(FABS(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) 
	{
	  x[i]=FABS(local_xold[i]+alam*p[i]);
	}

      /**/
      for(i=0;i<n;i++)
	{
	  tree->mod->m4mod->o_rr[i] = FABS(local_xold[i]+alam*p[i]);
	}
      /**/

      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
/* 	  PhyML_Printf("LOGlk = %f\n",*f); */
	}
      else *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) 
	    {
	      x[i]=FABS(local_xold[i]);
	    }
	  /**/      
	  for(i=0;i<n;i++)
	    {
	      tree->mod->m4mod->o_rr[i] = FABS(local_xold[i]);
	    }
	  /**/

	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return 0;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return 0;
	}
      else 
	{
/* 	  if (alam == 1.0) */
	  if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
/* 	      if (a == 0.0) tmplam = -slope/(2.0*b); */
	      if (a < SMALL && a > -SMALL) tmplam = -slope/(2.0*b);
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) tmplam = 0.5*alam;
		  else if(b <= 0.0) tmplam=(-b+SQRT(disc))/(3.0*a);
		  else tmplam = -slope/(b+SQRT(disc));
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
  return 1;
}

#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch_Nucleotide_Frequencies(t_tree *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,
				   phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=FABS(p[i])/MAX(FABS(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=FABS(local_xold[i]+alam*p[i]);
      /**/      
      for(i=0;i<n;i++) 
	{
	  tree->mod->pi_unscaled[i]=FABS(local_xold[i]+alam*p[i]);
/* 	  if( */
/* 	     (tree->mod->pi[i] < 0.001) || */
/* 	     (tree->mod->pi[i] > 0.999) */
/* 	     ) */
/* 	    break; */
	}
      /**/
      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->pi_unscaled[i]=local_xold[i];
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return 0;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return 0;
	}
      else 
	{
	  if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
/* 	  if (alam == 1.0) */
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a < SMALL && a > -SMALL) tmplam = -slope/(2.0*b);
/* 	      if (a == 0.0) tmplam = -slope/(2.0*b); */
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+SQRT(disc))/(3.0*a);
		      else tmplam = -slope/(b+SQRT(disc));
		    }
		  else tmplam=(-b+SQRT(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
  return 1;
}

/*********************************************************/


#define ALF 1.0e-4
#define TOLX 1.0e-7

int Lnsrch_Free_Mixt_Rates(t_tree *tree, int n, phydbl *xold, phydbl fold, phydbl *g, phydbl *p, phydbl *x,
			   phydbl *f, phydbl stpmax, int *check)
{
  int i;
  phydbl a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,test,tmplam;
  phydbl *local_xold;

  alam = alam2 = f2 = fold2 = tmplam = .0;

  local_xold = (phydbl *)mCalloc(n,sizeof(phydbl));
  For(i,n) local_xold[i] = xold[i];


  *check=0;
  for(sum=0.0,i=0;i<n;i++) sum += p[i]*p[i];
  sum=SQRT(sum);
  if(sum > stpmax)
    for(i=0;i<n;i++) p[i] *= stpmax/sum;
  for(slope=0.0,i=0;i<n;i++)
    slope += g[i]*p[i];
  test=0.0;
  for(i=0;i<n;i++) 
    {
      temp=FABS(p[i])/MAX(FABS(local_xold[i]),1.0);
      if (temp > test) test=temp;
    }
  alamin=TOLX/test;
  alam=1.0;
  for (;;) 
    {
      for(i=0;i<n;i++) x[i]=FABS(local_xold[i]+alam*p[i]);
      /**/      
      for(i=0;i<n;i++) 
	{
	  tree->mod->gamma_r_proba[i]=FABS(local_xold[i]+alam*p[i]);
	}
      /**/
      if(i==n) 
	{
	  *f=Return_Abs_Lk(tree);
	}
      else     *f=1.+fold+ALF*alam*slope;
      if (alam < alamin)
	{
	  for (i=0;i<n;i++) x[i]=local_xold[i];
	  for (i=0;i<n;i++) tree->mod->gamma_r_proba[i]=local_xold[i];
	  *check=1;
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold);
	  return 0;
	} 
      else if (*f <= fold+ALF*alam*slope) 
	{
	  For(i,n) xold[i] = local_xold[i];
	  Free(local_xold); 
	  return 0;
	}
      else 
	{
	  if ((alam < 1.0+SMALL) && (alam > 1.0-SMALL))
/* 	  if (alam == 1.0) */
	    tmplam = -slope/(2.0*(*f-fold-slope));
	  else 
	    {
	      rhs1 = *f-fold-alam*slope;
	      rhs2=f2-fold2-alam2*slope;
	      a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
	      b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
	      if (a < SMALL && a > -SMALL) tmplam = -slope/(2.0*b);
/* 	      if (a == 0.0) tmplam = -slope/(2.0*b); */
	      else 
		{
		  disc=b*b-3.0*a*slope;
		  if (disc<0.0) 
		    {
		      disc=b*b-3.0*a*slope;
		      if (disc<0.0) tmplam = 0.5*alam;
		      else if(b <= 0.0) tmplam=(-b+SQRT(disc))/(3.0*a);
		      else tmplam = -slope/(b+SQRT(disc));
		    }
		  else tmplam=(-b+SQRT(disc))/(3.0*a);
		}
	      if (tmplam>0.5*alam) tmplam=0.5*alam;
	    }
	}
      alam2=alam;
      f2 = *f;
      fold2=fold;
      alam=MAX(tmplam,0.1*alam);
    }
  Free(local_xold);
  return 1;
}


/* void Optimize_Global_Rate(t_tree *tree) */
/* { */
/*     PhyML_Printf("\n. Global rate (%f->)",tree->c_lnL); */
/*     Optimize_Single_Param_Generic(tree,&(tree->tbl),tree->tbl,tree->l_min,1.E+4,100); */
/*     PhyML_Printf("%f)\n",tree->c_lnL); */
/* } */


#undef ALF
#undef TOLX
#undef NRANSI

/*********************************************************/

int Dist_F_Brak(phydbl *ax, phydbl *bx, phydbl *cx, phydbl *F, phydbl *param, model *mod)
{
   phydbl ulim,u,r,q,dum;
   phydbl fa, fb, fc, fu;

   fa = -Lk_Dist(F,FABS(*ax),mod);
   fb = -Lk_Dist(F,FABS(*bx),mod);

   if(fb > fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,fb,fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   fc = -Lk_Dist(F,FABS(*cx),mod);

   while (fb > fc) 
     {
       r=(*bx-*ax)*(fb-fc);
       q=(*bx-*cx)*(fb-fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   fu = -Lk_Dist(F,FABS(u),mod);
	   if (fu < fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       fa=fb;
	       fb=fu;
	       return(0);
	     } 
	   else if (fu > fb) 
	     {
	       *cx=u;
	       fc=fu;
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = -Lk_Dist(F,FABS(u),mod);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   fu = -Lk_Dist(F,FABS(u),mod);
	   if (fu < fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
	       SHFT(fb,fc,fu,-Lk_Dist(F,FABS(u),mod))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u  = ulim;
	   fu = -Lk_Dist(F,FABS(u),mod);
	 } 
       else 
	 {
	   u  =(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = -Lk_Dist(F,FABS(u),mod);
	 }

       SHFT(*ax,*bx,*cx,u)
       SHFT(fa,fb,fc,fu)
      }
   return(0);
}

/*********************************************************/

phydbl Dist_F_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
		    phydbl *param, phydbl *F, model *mod)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL, curr_lnL;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  old_lnL = UNLIKELY;
  fw = fv = fx = -Lk_Dist(F,FABS(bx),mod);
  curr_lnL = init_lnL = -fw;

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);

      tol2=2.0*(tol1=tol*FABS(x)+BRENT_ZEPS);

      if(
	 ((FABS(curr_lnL-old_lnL) < mod->s_opt->min_diff_lk_local) && 
	  (curr_lnL > init_lnL - mod->s_opt->min_diff_lk_local)) ||	 
	  (iter > n_iter_max - 1)
	 )	 
	{
	  *param = x;
	  curr_lnL = Lk_Dist(F,*param,mod);
	  return -curr_lnL;
	}
      
      if(FABS(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      /*                   PhyML_Printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-x);
	      /*                   PhyML_Printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  /*               PhyML_Printf("Golden section step (default)\n"); */
	}
      
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = FABS(u);
      old_lnL = curr_lnL;
      fu = -Lk_Dist(F,FABS(u),mod);
      curr_lnL = -fu;      
/*       PhyML_Printf("param=%f LOGlk=%f\n",*param,fu); */
      
/*       if(fu <= fx)  */
      if(fu < fx) 
	{
	  if(iter > n_iter_max) return -fu;

	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x)  */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
/* 	  else if (fu <= fv || v == x || v == w)  */
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/

void Opt_Dist_F(phydbl *dist, phydbl *F, model *mod)
{
  phydbl ax,bx,cx;

  if(*dist < mod->l_min) *dist = mod->l_min;

  ax = mod->l_min;
  bx =  (*dist);
  cx = mod->l_max;

/*   Dist_F_Brak(&ax,&bx,&cx,F,dist,mod); */
  Dist_F_Brent(ax,bx,cx,1.E-10,1000,dist,F,mod);
}

/*********************************************************/

int Missing_Dist_Brak(phydbl *ax, phydbl *bx, phydbl *cx, int x, int y, matrix *mat)
{
   phydbl ulim,u,r,q,dum;
   phydbl fa, fb, fc, fu;

   fa = Least_Square_Missing_Dist_XY(x,y,FABS(*ax),mat);
   fb = Least_Square_Missing_Dist_XY(x,y,FABS(*bx),mat);

   if(fb > fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,fb,fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*((*bx)-(*ax));
   fc = Least_Square_Missing_Dist_XY(x,y,FABS(*cx),mat);

   while (fb > fc) 
     {
       r=((*bx)-(*ax))*(fb-fc);
       q=((*bx)-(*cx))*(fb-fa);
       u=(*bx)-(((*bx)-(*cx))*q-((*bx)-(*ax))*r)/
               (2.0*SIGN(MAX(FABS(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
	   if (fu < fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       fa=fb;
	       fb=fu;
	       return(0);
	     } 
	   else if (fu > fb) 
	     {
	       *cx=u;
	       fc=fu;
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
	   if (fu < fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))
		 SHFT(fb,fc,fu,Least_Square_Missing_Dist_XY(x,y,FABS(u),mat))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u  = ulim;
	   fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
	 } 
       else 
	 {
	   u  =(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
	 }

       SHFT(*ax,*bx,*cx,u)
       SHFT(fa,fb,fc,fu)
      }
   return(0);
}

/*********************************************************/

phydbl Missing_Dist_Brent(phydbl ax, phydbl bx, phydbl cx, phydbl tol, int n_iter_max, 
			  int x, int y, matrix *mat)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,xx,xm;
  phydbl e=0.0;
  phydbl init_LOGlk, max_LOGlk;
  phydbl bestx;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  xx=w=v=bx;
  fx=Least_Square_Missing_Dist_XY(x,y,FABS(bx),mat);
  fw=fv=-fx;
  init_LOGlk = fw;
  max_LOGlk = UNLIKELY;
  bestx = bx;

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*FABS(xx)+BRENT_ZEPS);

      if(FABS(xx-xm) <= (tol2-0.5*(b-a))) 
	{
	  mat->dist[x][y] = xx;
	  Least_Square_Missing_Dist_XY(x,y,mat->dist[x][y],mat);
	  return -fx;
	}
      
      if(FABS(e) > tol1) 
	{
	  r=(xx-w)*(fx-fv);
	  q=(xx-v)*(fx-fw);
	  p=(xx-v)*q-(xx-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-xx) || p >= q*(b-xx))
	    {
	      d=BRENT_CGOLD*(e=(xx >= xm ? a-xx : b-xx));
	      /*                   PhyML_Printf("Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=xx+d;
	      if (u-a < tol2 || b-u < tol2)
		d=SIGN(tol1,xm-xx);
	      /*                   PhyML_Printf("Parabolic step\n"); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(xx >= xm ? a-xx : b-xx));
	  /*               PhyML_Printf("Golden section step (default)\n"); */
	}
      
      u=(FABS(d) >= tol1 ? xx+d : xx+SIGN(tol1,d));
      fu = Least_Square_Missing_Dist_XY(x,y,FABS(u),mat);
            
/*       PhyML_Printf("param=%f LOGlk=%f\n",u,fu); */
      
/*       if(fu <= fx)  */
      if(fu < fx) 
	{
	  if(iter > n_iter_max) return -fu;

	  if(u >= xx) a=xx; else b=xx;
	  SHFT(v,w,xx,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < xx) a=u; else b=u;
/* 	  if (fu <= fw || w == xx)  */
	  if (fu < fw || FABS(w-xx) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
/* 	  else if (fu <= fv || v == xx || v == w)  */
	  else if (fu < fv || FABS(v-xx) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }
  Exit("\n. Too many iterations in BRENT !");
  return(-1);
}

/*********************************************************/

void Opt_Missing_Dist(int x, int y, matrix *mat)
{
  phydbl ax,bx,cx;

  ax = DIST_MAX;
  bx = DIST_MAX/4.;

  Missing_Dist_Brak(&ax,&bx,&cx,x,y,mat);
  PhyML_Printf("ax=%f bx=%f cx=%f\n",FABS(ax),FABS(bx),FABS(cx));
  Missing_Dist_Brent(FABS(ax),FABS(bx),FABS(cx),1.E-5,100,x,y,mat);
}

/*********************************************************/

int Optimiz_Alpha_And_Pinv(t_tree *tree)
{

  int    iter;
  phydbl best_alpha, best_pinv, best_mult;
  phydbl slope, intercept;
  phydbl lk_b, lk_a;
  phydbl lk_init, lk_final;
  phydbl f0,f1,f2,f3,x0,x1,x2,x3;
  phydbl pinv0, pinv1;
  phydbl a, b, c;
  phydbl fa, fb, fc;
  phydbl K;
  phydbl alpha0, alpha1;
  phydbl best_lnL;


  lk_final = UNLIKELY;
  lk_b     = UNLIKELY;
  lk_a     = UNLIKELY;


/*   PhyML_Printf("\n\n. Init lnL = %f alpha=%f pinv=%f", */
/* 	 tree->c_lnL, */
/* 	 tree->mod->alpha, */
/* 	 tree->mod->pinvar); */
    
  /* Two (full) steps to compute  pinv_alpha_slope & pinv_alpha_intercept */

  tree->both_sides = 1;
  Lk(tree);
  lk_b = tree->c_lnL;

  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],tree,tree->data);
  

  tree->both_sides = 0;
  Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar),.0001,0.9999,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max,
				tree->mod->s_opt->quickdirty);

  Optimize_Single_Param_Generic(tree,&(tree->mod->alpha),.01,100.,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max,
				tree->mod->s_opt->quickdirty);
  
  pinv0  = tree->mod->pinvar;
  alpha0 = tree->mod->alpha;
  f0 = tree->c_lnL;


  tree->both_sides = 1;
  Lk(tree);

  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
			tree->noeud[0]->b[0],tree,tree->data);
    

  tree->both_sides = 0;
  Optimize_Single_Param_Generic(tree,&(tree->mod->pinvar),.0001,0.9999,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max,
				tree->mod->s_opt->quickdirty);
  Optimize_Single_Param_Generic(tree,&(tree->mod->alpha),.01,100.,
				tree->mod->s_opt->min_diff_lk_global,
				tree->mod->s_opt->brent_it_max,
				tree->mod->s_opt->quickdirty);

  lk_a = tree->c_lnL;

  pinv1  = tree->mod->pinvar;
  alpha1 = tree->mod->alpha;
  f1 = tree->c_lnL;
  best_lnL = f1;


  if(lk_a < lk_b - tree->mod->s_opt->min_diff_lk_local)
    {
      PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
      Warn_And_Exit("");
    }
  else if(FABS(lk_a - lk_b) < tree->mod->s_opt->min_diff_lk_local)
    {
      return 1;
    }
    
  Record_Br_Len(tree);
  best_alpha = tree->mod->alpha;
  best_pinv  = tree->mod->pinvar;
  best_mult  = tree->mod->br_len_multiplier;
  lk_init    = tree->c_lnL;

  /* PhyML_Printf("\n\n. Init lnL after std opt = %f [%f] best_alpha=%f best_pinv=%f",tree->c_lnL,Lk(tree),best_alpha,best_pinv); */
  /* PhyML_Printf("\n. Best_lnL = %f",best_lnL); */

  slope     = (pinv1 - pinv0)/(alpha1 - alpha0);
  intercept = pinv1 - slope * alpha1;
  
  if((slope > 0.001) && (slope < 1./0.001))
    {
      /* PhyML_Printf("\n. pinv0 = %f, pinv1 = %f, alpha0 = %f, alpha1 = %f",pinv0,pinv1,alpha0,alpha1); */
      /* PhyML_Printf("\n. slope = %f intercept = %f",slope,intercept); */
      
      K = 0.381966;
      
      if(alpha1 < alpha0) 
	{
	  c  = alpha0;
	  b  = alpha1;
	  fc = f0;
	  fb = f1;
	  
	  a = (0.1 < alpha1)?(0.1):(0.5*alpha1);
	  tree->mod->alpha = a;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);

	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  

	  fa = tree->c_lnL;
	  
	  iter = 0;	

	  /* PhyML_Printf("\n. a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f (alpha=%f pinv=%f)",a,b,c,fa,fb,fc,tree->mod->alpha,tree->mod->pinvar); */

	  while(fa > fb)
	    {
	      a = a/5.;
	      tree->mod->alpha = a;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      fa = tree->c_lnL;
	      /* PhyML_Printf("\n1 a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f",a,b,c,fa,fb,fc); */
	      if(iter++ > 10) return 0;
	    }
	}
      else
	{
	  a  = alpha0;
	  b  = alpha1;
	  fa = f0;
	  fb = f1;
	  
	  c = (alpha1 < 2.)?(2.0):(2.*alpha1);
	  tree->mod->alpha = c;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  fc = tree->c_lnL;

	  /* PhyML_Printf("\n. a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f (alpha=%f pinv=%f)",a,b,c,fa,fb,fc,tree->mod->alpha,tree->mod->pinvar); */

	  iter = 0;
	  while(fc > fb)
	    {
	      c = c*2.;
	      tree->mod->alpha = c;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      fc = tree->c_lnL;
	      /* PhyML_Printf("\n2 a=%f, b=%f, c=%f, fa=%f, fb=%f, fc=%f",a,b,c,fa,fb,fc); */
	      if(iter++ > 10) return 0;
	    }
	}
      
      
      if(FABS(b - c) > FABS(a - b))
	{
	  x0 = a; x1 = b; x3 = c;
	  x2 = b + K * FABS(b - c);
	  
	  f0 = fa;
	  f1 = fb;
	  f3 = fc;
	  tree->mod->alpha = x2;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  f2 = tree->c_lnL;
	}
      else /* |b -c| < |a - b| */
	{
	  x0 = a; x2 = b; x3 = c;
	  x1 = b - K * FABS(b - a);
	  
	  f0 = fa;
	  f2 = fb;
	  f3 = fc;
	  tree->mod->alpha = x1;
	  tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	  if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	  if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	  tree->both_sides = 1;
	  Lk(tree);
	  Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				tree->noeud[0]->b[0],tree,tree->data);
	  f1 = tree->c_lnL;
	}
      
      iter = 0;
      do
	{
	  /* PhyML_Printf("\n. x0=%f, x1=%f, x2=%f, x3=%f, f0=%f, f1=%f, f2=%f, f3=%f", */
	  /* 	 x0,x1,x2,x3,f0,f1,f2,f3); */
	  
	  if(f1 > f2)
	    {
	      x3 = x2;
	      x2 = x1;
	      x1 = x2 - K * FABS(x2 - x0);
	      
	      f3 = f2;
	      f2 = f1;
	      
	      tree->mod->alpha = x1;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      f1 = tree->c_lnL;
	      if(f1 > best_lnL) 
		{
		  Record_Br_Len(tree);
		  best_alpha = tree->mod->alpha;
		  best_pinv = tree->mod->pinvar;
		  best_mult  = tree->mod->br_len_multiplier;
		  /* PhyML_Printf("\n>.< New alpha=%f pinv=%f",best_alpha,best_pinv); */
		}
	      /* PhyML_Printf("\n> f1=%f",f1); */
	    }
	  else /* f1 < f2 */
	    {
	      x0 = x1;
	      x1 = x2;
	      x2 = x2 + K * FABS(x3 - x2);
	      
	      f0 = f1;
	      f1 = f2;
	      
	      tree->mod->alpha = x2;
	      tree->mod->pinvar = slope * tree->mod->alpha + intercept;
	      if(tree->mod->pinvar > 1.0) tree->mod->pinvar = 0.9;
	      if(tree->mod->pinvar < 0.0) tree->mod->pinvar = 0.001;
	      tree->both_sides = 1;
	      Lk(tree);
	      Optimize_Br_Len_Serie(tree->noeud[0],tree->noeud[0]->v[0],
				    tree->noeud[0]->b[0],tree,tree->data);
	      f2 = tree->c_lnL;
	      if(f2 > best_lnL) 
		{
		  Record_Br_Len(tree);
		  best_alpha = tree->mod->alpha;
		  best_pinv = tree->mod->pinvar;
		  best_mult  = tree->mod->br_len_multiplier;
		  /* PhyML_Printf("\n>o< New alpha=%f pinv=%f",best_alpha,best_pinv); */
		}
	      /* PhyML_Printf("\n> f2=%f",f2); */
	    }
	  
	  if(FABS(f1 - f2) < 0.01) break;
	  
	  iter++;
	  
	}while(iter < 100);
    }
  
  tree->mod->alpha = best_alpha;
  tree->mod->pinvar = best_pinv;
  tree->mod->br_len_multiplier = best_mult;
  Restore_Br_Len(tree);      
  tree->both_sides = 1;
  Lk(tree);
  /* PhyML_Printf("\n\n. Init lnL after golden opt = %f [%f] best_alpha=%f best_pinv=%f",tree->c_lnL,Lk(tree),best_alpha,best_pinv); */
  return 1;
}

/*********************************************************/

phydbl Generic_Brent_Lk(phydbl *param, phydbl ax, phydbl cx, phydbl tol, 
			int n_iter_max, int quickdirty,
			phydbl (*obj_func)(t_edge *,t_tree *,supert_tree *), 
			t_edge *branch, t_tree *tree, supert_tree *stree)
{
  int iter;
  phydbl a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  phydbl e=0.0;
  phydbl old_lnL,init_lnL;
  phydbl bx = *param;

  d=0.0;
  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x=w=v=bx;
  old_lnL = UNLIKELY;
  (*param) = bx;
  fw=fv=fx=fu=-(*obj_func)(branch,tree,stree);
  init_lnL = -fw;

  /* PhyML_Printf("\n. init_lnL = %f a=%f b=%f c=%f",init_lnL,ax,bx,cx); */

  for(iter=1;iter<=BRENT_ITMAX;iter++) 
    {
      xm=0.5*(a+b);
      tol2=2.0*(tol1=tol*x+BRENT_ZEPS);

      if((fu > init_lnL + tol) && (quickdirty))
	{
	  (*param) = x;
	  fu = (*obj_func)(branch,tree,stree);
/* 	  Exit("\n"); */
	  return fu;	  
	}

/*       if(((FABS(cur_lnL-old_lnL) < tol) && (cur_lnL > init_lnL - tol)) || (iter > n_iter_max - 1)) */
      if((FABS(fu-old_lnL) < tol) || (iter > n_iter_max - 1))
	{
	  (*param) = x;
	  fu = (*obj_func)(branch,tree,stree);
/* 	  Exit("\n"); */
	  return fu;	  
	}
      
      if(FABS(e) > tol1) 
	{
	  r=(x-w)*(fx-fv);
	  q=(x-v)*(fx-fw);
	  p=(x-v)*q-(x-w)*r;
	  q=2.0*(q-r);
	  if(q > 0.0) p = -p;
	  q=FABS(q);
	  etemp=e;
	  e=d;
	  if(FABS(p) >= FABS(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
	    {
	      d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	      /* PhyML_Printf(" Golden section step\n"); */
	    }
	  else
	    {
	      d=p/q;
	      u=x+d;
	      if (u-a < tol2 || b-u < tol2) d=SIGN(tol1,xm-x);
	      /* PhyML_Printf(" Parabolic step [e=%f]\n",e); */
	    }
        }
      else
	{
	  d=BRENT_CGOLD*(e=(x >= xm ? a-x : b-x));
	  /* PhyML_Printf(" Golden section step (default) [e=%f tol1=%f a=%f b=%f d=%f]\n",e,tol1,a,b,d); */
	}
      
      u=(FABS(d) >= tol1 ? x+d : x+SIGN(tol1,d));
      (*param) = u;
      old_lnL = fu;
      fu = -(*obj_func)(branch,tree,stree);
      
      /* PhyML_Printf("\n. iter=%d/%d param=%f lnL=%f",iter,BRENT_ITMAX,*param,fu); */

      if(fu <= fx)
	{
	  if(u >= x) a=x; else b=x;
	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	} 
      else
	{
	  if (u < x) a=u; else b=u;
/* 	  if (fu <= fw || w == x) */
	  if (fu < fw || FABS(w-x) < SMALL)
	    {
	      v=w;
	      w=u;
	      fv=fw;
	      fw=fu;
	    } 
/* 	  else if (fu <= fv || v == x || v == w) */
	  else if (fu < fv || FABS(v-x) < SMALL || FABS(v-w) < SMALL)
	    {
	      v=u;
	      fv=fu;
	    }
	}
    }

  Exit("\n. Too many iterations in BRENT !");
  return(-1);
  /* Not Reached ??  *param=x;   */
  /* Not Reached ??  return fx; */
}

/*********************************************************/

/* find ML erstimates of node heights given fixed substitution
   rates on branches. Also optimizes the overall substitution
   rate */
void Round_Optimize_Node_Heights(t_tree *tree)
{
  phydbl cur_lnL, new_lnL;
  int n_iter;


  cur_lnL = UNLIKELY;
  new_lnL = Lk(tree);

  
  printf("\n. cur_lnL = %f new_lnL=%f",cur_lnL,new_lnL);

  
  n_iter = 0;
  while(fabs(new_lnL - cur_lnL) > tree->mod->s_opt->min_diff_lk_global)
    {
      cur_lnL = tree->c_lnL;

      Opt_Node_Heights_Recurr(tree);
            
      Generic_Brent_Lk(&(tree->rates->clock_r),
      		       tree->rates->min_clock,
      		       tree->rates->max_clock,
      		       tree->mod->s_opt->min_diff_lk_global,
      		       tree->mod->s_opt->brent_it_max,
      		       tree->mod->s_opt->quickdirty,
      		       Wrap_Lk,NULL,tree,NULL);

      printf("\n. cur_lnL=%f new_lnL=%f clock_r=%G root height=%f",
	     cur_lnL,new_lnL,tree->rates->clock_r,tree->rates->nd_t[tree->n_root->num]);
      new_lnL = tree->c_lnL;
      n_iter++;
      if(n_iter > 100) break;
    }
}

/*********************************************************/

void Opt_Node_Heights_Recurr(t_tree *tree)
{
  Opt_Node_Heights_Recurr_Pre(tree->n_root,tree->n_root->v[0],tree);
  Opt_Node_Heights_Recurr_Pre(tree->n_root,tree->n_root->v[1],tree);

  Generic_Brent_Lk(&(tree->rates->nd_t[tree->n_root->num]),
		   MIN(tree->rates->t_prior_max[tree->n_root->num],
		       MIN(tree->rates->nd_t[tree->n_root->v[0]->num],
			   tree->rates->nd_t[tree->n_root->v[1]->num])),
		   tree->rates->t_prior_min[tree->n_root->num],
		   tree->mod->s_opt->min_diff_lk_global,
		   tree->mod->s_opt->brent_it_max,
		   tree->mod->s_opt->quickdirty,
		   Wrap_Lk,NULL,tree,NULL);
}

/*********************************************************/

void Opt_Node_Heights_Recurr_Pre(t_node *a, t_node *d, t_tree *tree)
{
  if(d->tax) return;
  else
    {
      int i;
      phydbl t0,t1,t2,t3;
      phydbl t_min,t_max;
      t_node *v2,*v3;

      
      v2 = v3 = NULL;
      For(i,3)
	if((d->v[i] != a) && (d->b[i] != tree->e_root))
	  {
	    if(!v2) { v2 = d->v[i]; }
	    else    { v3 = d->v[i]; }
	  }
      
      Opt_Node_Heights_Recurr_Pre(d,v2,tree);
      Opt_Node_Heights_Recurr_Pre(d,v3,tree);

      t0 = tree->rates->nd_t[a->num];
      t1 = tree->rates->nd_t[d->num];
      t2 = tree->rates->nd_t[v2->num];
      t3 = tree->rates->nd_t[v3->num];
      
      t_min = t0;
      t_max = MIN(t2,t3);
      
      t_min = MAX(t_min,tree->rates->t_prior_min[d->num]);
      t_max = MIN(t_max,tree->rates->t_prior_max[d->num]);
      
      t_min += tree->rates->min_dt;
      t_max -= tree->rates->min_dt;
      
      if(t_min > t_max)
	{
	  PhyML_Printf("\n. Err in file %s at line %d\n",__FILE__,__LINE__);
	  Warn_And_Exit("");
	}
      
      Generic_Brent_Lk(&(tree->rates->nd_t[d->num]),
		       t_min,t_max,
		       tree->mod->s_opt->min_diff_lk_global,
		       tree->mod->s_opt->brent_it_max,
		       tree->mod->s_opt->quickdirty,
		       Wrap_Lk,NULL,tree,NULL);
      

      /* printf("\n. t%d = %f [%f;%f] lnL = %f",d->num,tree->rates->nd_t[d->num],t_min,t_max,tree->c_lnL); */

    }
}


/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
/*********************************************************/
