/*############################################################*/
/*                                                            */
/* streambugs 1.2dev                                          */
/*                                                            */
/* creation:                  09.01.2013                      */
/* modifications:             06.09.2013                      */
/* revision for version 1.2:  15.01.2014                      */
/*                                                            */
/*############################################################*/


#include <R.h>

char buffer[256];


/* ========================================================== */
/* definition of global pointers to dynamically allocated     */
/* state, parameter and input values, indices and names       */
/* ========================================================== */

struct input            {
                          int                   n;
                          double*               x;
                          double*               y;
                        };
                      
struct inpparvector     {
                          int                   n_inp;
                          int*                  inpinds;
                          int                   n_par;
                          double*               parvals;
                          char**                parnames;
                        };

struct inpparmatrix     {
                          int                   n_inp;
                          int*                  inpinds;
                          int                   n_par;
                          int                   n_row;
                          double*               parvals;
                          char**                parnames;
                        };
                    
struct inpparvectorlist {
                          int                   n;
                          char**                names;
                          struct inpparvector** vectors;
                        };

struct inpparmatrixlist {
                          int                   n;
                          char**                names;
                          struct inpparmatrix** matrices;
                        };

struct stoichvector     {
                          int                   n;
                          int*                  inds;
                          double*               coeff;
                          char**                names;
                        };
                    
struct process          {
                          char*                 name;
                          struct inpparvector*  par;
                          struct stoichvector*  stoich;
                        };
                    
struct taxprocesses     {
                          int                   n;
                          struct process*       procs;
                        };

struct webprocesses     {
                          int                   n;
                          struct webprocess*    webprocs;
                        };

struct webprocess       {
                          int                   n;
                          char*                 name;
                          struct inpparvector*  par;
                          struct process*       procs;
                        };
                    
struct indfA            {
                          int                   nreachind;
                          int*                  reachind;
                          int                   nfsthabind;
                          int*                  fsthabind;
                        };

static double                   g_parms;              /* dummy var.  for init. */

static int                      g_ninp;               /* number of inputs      */

static struct input*            g_inp;                /* inputs                */

static struct inpparvector*     g_parglobal      = 0; /* global parameters     */

static struct inpparvectorlist* g_parglobalenvtraits
                                                 = 0; /* trait-rel. env. cond. */

static struct inpparmatrix*     g_parenvcondreach= 0; /* reach-dep. env. cond. */

static struct inpparmatrix*     g_parenvcondhabitat 
                                                 = 0; /* hab.-dep. env. cond.  */

static struct inpparmatrixlist* g_parenvcondhabitatgroup
                                                 = 0; /* grouped hab.-dep. e.c.*/

static struct inpparmatrix*     g_parinitcond    = 0; /* initial conditions    */

static struct inpparmatrix*     g_parinput       = 0; /* input                 */

static struct inpparmatrix*     g_partaxapropdirect 
                                                 = 0; /* dir. def. taxa prop.  */

static struct inpparmatrixlist* g_partaxaproptraits
                                                 = 0; /* trait-dep. taxa prop. */

static struct taxprocesses*     g_proctaxon      = 0; /* taxon-based processes */

static struct webprocesses*     g_procweb        = 0; /* consumption processes */

static int                      g_ny             = 0; /* number of state vars. */

static int                      g_nreach         = 0;

static struct indfA*            g_indfA          = 0;

static int                      g_debug          = 0;

/* ========================================================== */
/* functions to allocate and assign data                       */
/* ========================================================== */


struct inpparvector* create_inpparvector(int*    n_inp, 
                                         int*    inpinds,
                                         int*    n_par, 
                                         double* parvals, 
                                         char**  parnames)
{
   struct inpparvector* s; 
                                           
   s = (struct inpparvector*)malloc(sizeof(struct inpparvector));
   
   (*s).n_inp = *n_inp;
   if ( *n_inp > 0 )
   {
      (*s).inpinds = malloc(2*(*n_inp)*sizeof(int));
      for ( int i=0; i<2*(*n_inp); i++ ) (*s).inpinds[i] = inpinds[i];
   }

   (*s).n_par = *n_par;
   if ( *n_par > 0 )
   {
      (*s).parvals  = malloc(*n_par*sizeof(double));
      (*s).parnames = malloc(*n_par*sizeof(char*));
      for ( int i=0; i<*n_par; i++ )
      {
         (*s).parvals[i] = parvals[i]; 
         (*s).parnames[i] = malloc((strlen(parnames[i])+1)*sizeof(char));
         strcpy((*s).parnames[i],parnames[i]);
      }
   }
   
   return(s);
}


void delete_inpparvector(struct inpparvector* s)
{
   if ( (*s).n_inp > 0 )
   {
      free((*s).inpinds);
   }
   
   if ( (*s).n_par > 0 )
   {
      free((*s).parvals);
      for ( int i=0; i<(*s).n_par; i++ ) free((*s).parnames[i]);
      free((*s).parnames);
   }
   
   free(s);
}


struct inpparmatrix* create_inpparmatrix(int*    n_inp, 
                                         int*    inpinds,
                                         int*    n_par, 
                                         int*    n_row, 
                                         double* parvals, 
                                         char**  parnames)
{
   struct inpparmatrix* s; 

   s = (struct inpparmatrix*)malloc(sizeof(struct inpparmatrix));
   
   (*s).n_inp = *n_inp;
   if ( *n_inp > 0 )
   {
      (*s).inpinds = malloc(3*(*n_inp)*sizeof(int));
      for ( int i=0; i<3*(*n_inp); i++ ) (*s).inpinds[i] = inpinds[i];
   }
   
   (*s).n_par = *n_par;
   if ( *n_par > 0 )
   {
      (*s).parvals  = malloc((*n_par)*(*n_row)*sizeof(double));
      (*s).parnames = malloc(*n_par*sizeof(char*));
      for ( int i=0; i<(*n_par)*(*n_row); i++ )
      {
         (*s).parvals[i] = parvals[i]; 
      }
      for ( int i=0; i<*n_par; i++ )
      {
         (*s).parnames[i] = malloc((strlen(parnames[i])+1)*sizeof(char));
         strcpy((*s).parnames[i],parnames[i]);
      }
      (*s).n_row = *n_row;
   }
   
   return(s);
}


void delete_inpparmatrix(struct  inpparmatrix* s)
{
   if ( (*s).n_inp > 0 )
   {
      free((*s).inpinds);
      (*s).n_inp = 0;
   }
   
   if ( (*s).n_par > 0 )
   {
      free((*s).parvals);
      for ( int i=0; i<(*s).n_par; i++ ) free((*s).parnames[i]);
      free((*s).parnames);
      (*s).n_par = 0;
      (*s).n_row = 0;
   }
   
   free(s);
}


void delete_inpparvectorlist(struct inpparvectorlist* s)
{
   if ( (*s).n > 0 )
   {
      for ( int i=0; i<(*s).n; i++ )
      {
         delete_inpparvector((*s).vectors[i]);
         free((*s).names[i]);
      }
      free((*s).names);
      free((*s).vectors);
   }
   free(s);
}


void delete_inpparmatrixlist(struct inpparmatrixlist* s)
{
   if ( (*s).n > 0 )
   {
      for ( int i=0; i<(*s).n; i++ )
      {
         delete_inpparmatrix((*s).matrices[i]);
         free((*s).names[i]);
      }
      free((*s).names);
      free((*s).matrices);
   }
   free(s);
}


struct stoichvector* create_stoichvector(int*    n, 
                                         int*    inds,
                                         double* coeff, 
                                         char**  names)
{
   struct stoichvector* s;
   
   s = (struct stoichvector*)malloc(sizeof(struct stoichvector));
   
   (*s).n = *n;
   
   if ( *n > 0 )
   {
      (*s).inds  = malloc(*n*sizeof(int));
      (*s).coeff = malloc(*n*sizeof(double));
      (*s).names = malloc(*n*sizeof(char*));
      for ( int i=0; i<*n; i++ )
      {
         (*s).inds[i]  = inds[i]; 
         (*s).coeff[i] = coeff[i]; 
         (*s).names[i] = malloc((strlen(names[i])+1)*sizeof(char));
         strcpy((*s).names[i],names[i]);
      }
   }
   
   return(s);
}


void delete_stoichvector(struct  stoichvector* s)
{
   if ( (*s).n > 0 )
   {
      free((*s).inds);
      free((*s).coeff);
      for ( int i=0; i<(*s).n; i++ ) free((*s).names[i]);
      free((*s).names);
   }
   
   free(s);
}


int get_index(char* name, char** names, int n_names)
{
   for ( int i=0; i<n_names; i++ )
   {
      if ( strcmp(name,names[i]) == 0 )  return(i);
   }
   return(-1);
}


int exist_value_vector(char* name, struct inpparvector* s)
{
    int ind = get_index(name,(*s).parnames,(*s).n_par);
    if ( ind < 0 ) return(0);
    return(1);
}


double get_value_vector(char* name, struct inpparvector* s)
{
    int ind = get_index(name,(*s).parnames,(*s).n_par);
    if ( ind < 0 )
    {
       sprintf(buffer,"Parameter %s not found\n",name);
       error(buffer);
    }
    return( (*s).parvals[ind] );
}


int exist_value_matrix(char* name, struct inpparmatrix* s, int i)
{
    int ind = get_index(name,(*s).parnames,(*s).n_par);
    if ( ind < 0 ) return(0);
    if ( i<0 || i>=(*s).n_row ) return(0);
    return(1);
}


double get_value_matrix(char* name, struct inpparmatrix* s, int i)
{
    int ind = get_index(name,(*s).parnames,(*s).n_par);
    if ( ind < 0 )
    {
       sprintf(buffer,"Parameter %s not found\n",name);
       error(buffer);
    }
    return( (*s).parvals[(*s).n_row*ind+i] );
}


void get_stoich(char* name, struct stoichvector* s, int* indy, double* coeff)
{
   int ind = get_index(name,(*s).names,(*s).n);
   if ( ind < 0 )
   {
      sprintf(buffer,"Stoichiometric coefficient for %s not found\n",name);
      error(buffer);
   }
   *indy  = (*s).inds[ind];
   *coeff = (*s).coeff[ind];
   return;
}


double dmin(double a, double b)
{
   if ( a < b ) return(a);
   return(b);
}

double dmax(double a, double b)
{
   if ( a > b ) return(a);
   return(b);
}

int imin(int a, int b)
{
   if ( a < b ) return(a);
   return(b);
}

int imax(int a, int b)
{
   if ( a > b ) return(a);
   return(b);
}


/* ========================================================== */
/* functions to create global variables                       */
/* ========================================================== */


void streambugs_init_debug(int* streambugs_debug)
{
   g_debug = *streambugs_debug;
}


void streambugs_create_input_structure(int* ninp)
{
   if ( *ninp > 0 )
   {  
      g_inp = (struct input*)malloc(*ninp*sizeof(struct input));
      for ( int i=0; i<*ninp; i++ ) g_inp[i].n = 0;
      g_ninp = *ninp;
   }
}


void streambugs_create_input(int*    ind,
                             int*    n, 
                             double* x, 
                             double* y)
{
   if ( *ind < 1 || *ind > g_ninp ) error("illegal input initialization");

   if ( *n > 0 )
   {
      g_inp[*ind-1].x = malloc(*n*sizeof(double));
      g_inp[*ind-1].y = malloc(*n*sizeof(double));
      for ( int i=0; i<*n; i++ ) 
      {
          g_inp[*ind-1].x[i] = x[i];
          g_inp[*ind-1].y[i] = y[i];
      }
      g_inp[*ind-1].n = *n;
   }
}


void streambugs_delete_inputs()
{
   if ( g_ninp > 0 )
   {
      for ( int i=0; i<g_ninp; i++ )
      {
         if ( g_inp[i].n > 0 )
         {
            free(g_inp[i].x);
            free(g_inp[i].y);
         }
      }     
      free(g_inp);
      g_ninp = 0;
   }   
}


void streambugs_create_parglobal(int*    n_inp, 
                                 int*    inpinds,
                                 int*    n_par, 
                                 double* parvals, 
                                 char**  parnames)
{
   g_parglobal = create_inpparvector(n_inp,
                                     inpinds,
                                     n_par,
                                     parvals,
                                     parnames); 
}


void streambugs_delete_parglobal()
{
   delete_inpparvector(g_parglobal);
}


void streambugs_create_parglobalenvtraits_structure(int* n)
{
   g_parglobalenvtraits = (struct inpparvectorlist*)malloc(sizeof(struct inpparvectorlist));
   (*g_parglobalenvtraits).n = 0;
   if ( *n > 0 )
   {
      (*g_parglobalenvtraits).n = *n;
      (*g_parglobalenvtraits).names   = (char**)malloc(*n*sizeof(char*));
      (*g_parglobalenvtraits).vectors = (struct inpparvector**)malloc(*n*sizeof(struct inpparvector*));
   }
}


void streambugs_create_parglobalenvtrait(int*    ind,
                                         char**  name,
                                         int*    n_inp,
                                         int*    inpinds,
                                         int*    n_par,
                                         double* parvals,
                                         char**  parnames)
{
   if ( *ind < 1 || *ind > (*g_parglobalenvtraits).n ) error("illegal parglobalenvtrait initialization");

   (*g_parglobalenvtraits).names[*ind-1] = malloc((strlen(name[0])+1)*sizeof(char));
   strcpy((*g_parglobalenvtraits).names[*ind-1],name[0]);
   (*g_parglobalenvtraits).vectors[*ind-1] = create_inpparvector(n_inp, 
                                                                 inpinds,
                                                                 n_par, 
                                                                 parvals, 
                                                                 parnames);
}


void streambugs_delete_parglobalenvtraits()
{
   delete_inpparvectorlist(g_parglobalenvtraits);
}


void streambugs_create_parenvcondreach(int*    n_inp, 
                                       int*    inpinds,
                                       int*    n_par, 
                                       int*    n_row, 
                                       double* parvals, 
                                       char**  parnames)
{
   g_parenvcondreach = create_inpparmatrix(n_inp,
                                           inpinds,
                                           n_par,
                                           n_row,
                                           parvals,
                                           parnames);                            
}


void streambugs_delete_parenvcondreach()
{
   delete_inpparmatrix(g_parenvcondreach);
}


void streambugs_create_parenvcondhabitat(int*    n_inp, 
                                         int*    inpinds,
                                         int*    n_par, 
                                         int*    n_row, 
                                         double* parvals, 
                                         char**  parnames)
{
   g_parenvcondhabitat = create_inpparmatrix(n_inp,
                                             inpinds,
                                             n_par,
                                             n_row,
                                             parvals,
                                             parnames);                            
}


void streambugs_delete_parenvcondhabitat()
{
   delete_inpparmatrix(g_parenvcondhabitat);
}


void streambugs_create_parenvcondhabitatgroups_structure(int* n)
{
   g_parenvcondhabitatgroup = (struct inpparmatrixlist*)malloc(sizeof(struct inpparmatrixlist));
   (*g_parenvcondhabitatgroup).n = 0;
   if ( *n > 0 )
   {
      (*g_parenvcondhabitatgroup).n = *n;
      (*g_parenvcondhabitatgroup).names    = (char**)malloc(*n*sizeof(char*));
      (*g_parenvcondhabitatgroup).matrices = (struct inpparmatrix**)malloc(*n*sizeof(struct inpparmatrix*));
   }
}


void streambugs_create_parenvcondhabitatgroup(int*    ind,
                                              char**  name,
                                              int*    n_inp,
                                              int*    inpinds,
                                              int*    n_par,
                                              int*    n_row,
                                              double* parvals,
                                              char**  parnames)
{
   if ( *ind < 1 || *ind > (*g_parenvcondhabitatgroup).n ) error("illegal parenvcondhabitatgroup initialization");

   (*g_parenvcondhabitatgroup).names[*ind-1] = malloc((strlen(name[0])+1)*sizeof(char));
   strcpy((*g_parenvcondhabitatgroup).names[*ind-1],name[0]);
   (*g_parenvcondhabitatgroup).matrices[*ind-1] = create_inpparmatrix(n_inp, 
                                                                      inpinds,
                                                                      n_par, 
                                                                      n_row,
                                                                      parvals, 
                                                                      parnames);
}


void streambugs_delete_parenvcondhabitatgroup()
{
   delete_inpparmatrixlist(g_parenvcondhabitatgroup);
}


void streambugs_create_parinitcond(int*    n_inp, 
                                   int*    inpinds,
                                   int*    n_par, 
                                   int*    n_row, 
                                   double* parvals, 
                                   char**  parnames)
{
   g_parinitcond = create_inpparmatrix(n_inp,
                                       inpinds,
                                       n_par,
                                       n_row,
                                       parvals,
                                       parnames);                            
}


void streambugs_delete_parinitcond()
{
   delete_inpparmatrix(g_parinitcond);
}


void streambugs_create_parinput(int*    n_inp, 
                                int*    inpinds,
                                int*    n_par, 
                                int*    n_row, 
                                double* parvals, 
                                char**  parnames)
{
   g_parinput = create_inpparmatrix(n_inp,
                                    inpinds,
                                    n_par,
                                    n_row,
                                    parvals,
                                    parnames);                            
}


void streambugs_delete_parinput()
{
   delete_inpparmatrix(g_parinput);
}


void streambugs_create_partaxapropdirect(int*    n_inp, 
                                         int*    inpinds,
                                         int*    n_par, 
                                         int*    n_row, 
                                         double* parvals, 
                                         char**  parnames)
{
   g_partaxapropdirect = create_inpparmatrix(n_inp,
                                             inpinds,
                                             n_par,
                                             n_row,
                                             parvals,
                                             parnames);                            
}


void streambugs_delete_partaxapropdirect()
{
   delete_inpparmatrix(g_partaxapropdirect);
}


void streambugs_create_partaxaproptraits_structure(int* n)
{
   g_partaxaproptraits = (struct inpparmatrixlist*)malloc(sizeof(struct inpparmatrixlist));
   (*g_partaxaproptraits).n = 0;
   if ( *n > 0 )
   {
      (*g_partaxaproptraits).n = *n;
      (*g_partaxaproptraits).names    = (char**)malloc(*n*sizeof(char*));
      (*g_partaxaproptraits).matrices = (struct inpparmatrix**)malloc(*n*sizeof(struct inpparmatrix*));
   }
}


void streambugs_create_partaxaproptrait(int*    ind,
                                        char**  name,
                                        int*    n_inp,
                                        int*    inpinds,
                                        int*    n_par,
                                        int*    n_row,
                                        double* parvals,
                                        char**  parnames)
{
   if ( *ind < 1 || *ind > (*g_partaxaproptraits).n ) error("illegal partaxaproptraits initialization");

   (*g_partaxaproptraits).names[*ind-1] = malloc((strlen(name[0])+1)*sizeof(char));
   strcpy((*g_partaxaproptraits).names[*ind-1],name[0]);
   (*g_partaxaproptraits).matrices[*ind-1] = create_inpparmatrix(n_inp, 
                                                                 inpinds,
                                                                 n_par, 
                                                                 n_row,
                                                                 parvals, 
                                                                 parnames);
}


void streambugs_delete_partaxaproptraits()
{
   delete_inpparmatrixlist(g_partaxaproptraits);
}


void streambugs_create_processes_structure(int* ny)
{
   g_ny = *ny;
   if ( *ny > 0 )
   {
      g_proctaxon = (struct taxprocesses*)malloc(*ny*sizeof(struct taxprocesses));
      for ( int i=0; i<*ny; i++ )
      {
         g_proctaxon[i].n     = 0;
         g_proctaxon[i].procs = 0;
      }
      g_procweb   = (struct webprocesses*)malloc(*ny*sizeof(struct webprocesses));
      for ( int i=0; i<*ny; i++ )
      {
         g_procweb[i].n        = 0;  
         g_procweb[i].webprocs = 0; 
      } 
   }
}


void streambugs_create_proctaxon(int*    iy,
                                 int*    jproc,
                                 char**  procname,
                                 int*    n_inp,
                                 int*    inpinds,
                                 int*    n_par,
                                 char**  parnames,
                                 double* parvals,
                                 int*    n_stoich,
                                 char**  stoichnames,
                                 int*    stoichinds,
                                 double* stoichvals)
{
   if ( *iy > g_ny ) 
      error("steambugs_init_proctaxon: i too large");
   if ( g_proctaxon[*iy-1].n != *jproc-1 ) 
      error("streambugs_init_proctaxon: illegal value of j");
      
   /* allocate space for processes and copy existing processes */

   struct process* procs = g_proctaxon[*iy-1].procs;
   g_proctaxon[*iy-1].procs = 
      (struct process*)malloc(*jproc*sizeof(struct process));
   if ( g_proctaxon[*iy-1].n > 0 )
   {
      for ( int i=0; i<g_proctaxon[*iy-1].n; i++ )
      {
         g_proctaxon[*iy-1].procs[i] = procs[i];
      }
      free(procs);     
   }
   
   /* add new process */

   g_proctaxon[*iy-1].procs[*jproc-1].name = 
      malloc((strlen(procname[0])+1)*sizeof(char));
   strcpy(g_proctaxon[*iy-1].procs[*jproc-1].name,procname[0]);
   g_proctaxon[*iy-1].procs[*jproc-1].par = 
      create_inpparvector(n_inp,inpinds,n_par,parvals,parnames);
   g_proctaxon[*iy-1].procs[*jproc-1].stoich =
      create_stoichvector(n_stoich,stoichinds,stoichvals,stoichnames);
      
   /* modify number of processes */
   
   g_proctaxon[*iy-1].n = *jproc;
}


void streambugs_create_procweb(int*    iy,
                               int*    jproc,
                               char**  procname,
                               int*    n_inp,
                               int*    inpinds,
                               int*    n_par,
                               char**  parnames,
                               double* parvals)
{
   if ( *iy > g_ny ) 
      error("steambugs_init_procweb: i too large");
   if ( g_procweb[*iy-1].n != *jproc-1 ) 
      error("streambugs_init_procweb: illegal value of j");

   /* allocate space for webprocesses and copy existing webprocesses */
   
   struct webprocess* webprocs = g_procweb[*iy-1].webprocs;
   g_procweb[*iy-1].webprocs = 
      (struct webprocess*)malloc(*jproc*sizeof(struct webprocess));
   if ( g_procweb[*iy-1].n > 0 )
   {
      for ( int i=0; i<g_procweb[*iy-1].n; i++ )
      {
         g_procweb[*iy-1].webprocs[i] = webprocs[i];
      }
      free(webprocs);     
   }
   
   /* add new webprocess */
   
   g_procweb[*iy-1].webprocs[*jproc-1].n = 0;
   g_procweb[*iy-1].webprocs[*jproc-1].name = 
      malloc((strlen(procname[0])+1)*sizeof(char));
   strcpy(g_procweb[*iy-1].webprocs[*jproc-1].name,procname[0]);
   g_procweb[*iy-1].webprocs[*jproc-1].par = 
      create_inpparvector(n_inp,inpinds,n_par,parvals,parnames);
   g_procweb[*iy-1].webprocs[*jproc-1].procs = 0; 
      
   /* modify number of webprocesses */
   
   g_procweb[*iy-1].n = *jproc;                    
}


void streambugs_create_procwebtaxon(int*    iy,
                                    int*    jproc,
                                    int*    kproctax,
                                    char**  procname,
                                    int*    n_inp,
                                    int*    inpinds,
                                    int*    n_par,
                                    char**  parnames,
                                    double* parvals,
                                    int*    n_stoich,
                                    char**  stoichnames,
                                    int*    stoichinds,
                                    double* stoichvals)
{
   if ( *iy > g_ny ) 
      error("steambugs_init_procwebtaxon: i too large");
   if ( g_procweb[*iy-1].n != *jproc ) 
      error("streambugs_init_procwebtaxon: illegal value of j");
   if ( g_procweb[*iy-1].webprocs[*jproc-1].n != *kproctax-1 ) 
      error("streambugs_init_procwebtaxon: illegal value of k");
      
   /* allocate space for processes and copy existing processes */

   struct process* procs = g_procweb[*iy-1].webprocs[*jproc-1].procs;
   g_procweb[*iy-1].webprocs[*jproc-1].procs = 
      (struct process*)malloc(*kproctax*sizeof(struct process));
   if ( g_procweb[*iy-1].webprocs[*jproc-1].n > 0 )
   {
      for ( int i=0; i<g_procweb[*iy-1].webprocs[*jproc-1].n; i++ )
      {
         g_procweb[*iy-1].webprocs[*jproc-1].procs[i] = procs[i];
      }
      free(procs);     
   }
   
   /* add new process */
   
   g_procweb[*iy-1].webprocs[*jproc-1].procs[*kproctax-1].name = 
      malloc((strlen(procname[0])+1)*sizeof(char));
   strcpy(g_procweb[*iy-1].webprocs[*jproc-1].procs[*kproctax-1].name,procname[0]);
   g_procweb[*iy-1].webprocs[*jproc-1].procs[*kproctax-1].par = 
      create_inpparvector(n_inp,inpinds,n_par,parvals,parnames);
   g_procweb[*iy-1].webprocs[*jproc-1].procs[*kproctax-1].stoich =
      create_stoichvector(n_stoich,stoichinds,stoichvals,stoichnames);
      
   /* modify number of processes */
   
   g_procweb[*iy-1].webprocs[*jproc-1].n = *kproctax;
}


void streambugs_delete_processes()
{
   if ( g_ny > 0 )
   {
      /* free taxon-bases processes */
      
      for ( int i=0; i<g_ny; i++ )  /* loop over state variables */
      {
         if ( g_proctaxon[i].n > 0 )
         {
            for ( int j=0; j<g_proctaxon[i].n; j++ ) /* loop over taxon-based processes */
            {
               free(g_proctaxon[i].procs[j].name);
               delete_inpparvector(g_proctaxon[i].procs[j].par);
               delete_stoichvector(g_proctaxon[i].procs[j].stoich);
            }
            free(g_proctaxon[i].procs);
         }
      }
      free(g_proctaxon);
      
      /* free web processes */
      
      for ( int i=0; i<g_ny; i++ )  /* loop over state variables */
      {
         if ( g_procweb[i].n > 0 )
         {
            for ( int j=0; j<g_procweb[i].n; j++ )  /* loop over web processes */
            {
               free(g_procweb[i].webprocs[j].name);
               delete_inpparvector(g_procweb[i].webprocs[j].par);
               if ( g_procweb[i].webprocs[j].n > 0 )
               {
                  for ( int k=0; k<g_procweb[i].webprocs[j].n; k++ ) /* loop over subprocesses */
                  {
                     free(g_procweb[i].webprocs[j].procs[k].name);
                     delete_inpparvector(g_procweb[i].webprocs[j].procs[k].par);
                     delete_stoichvector(g_procweb[i].webprocs[j].procs[k].stoich);
                  }
                  free(g_procweb[i].webprocs[j].procs);
               }
            }
            free(g_procweb[i].webprocs);
         }
      }
      free(g_procweb);
   }
}


void streambugs_create_fA_structure(int* nreach)
{
   g_nreach = *nreach;
   if ( *nreach > 0 )
   {
      g_indfA = (struct indfA*)malloc(*nreach*sizeof(struct indfA));
      for ( int i=0; i<*nreach; i++ )
      {
         g_indfA[i].nreachind  = 0;
         g_indfA[i].nfsthabind = 0;
      }
   }
}


void streambugs_create_fA(int* ireach,
                          int* nreachind,
                          int* reachind,
                          int* nfsthabind,
                          int* fsthabind)
{
   if ( *ireach < 1 || *ireach > g_nreach ) error("steambugs_create_fA: illegal index");
      
   g_indfA[*ireach-1].nreachind = *nreachind;
   if ( *nreachind > 0 )
   {
      g_indfA[*ireach-1].reachind  = malloc(*nreachind*sizeof(int));
      for ( int i=0; i<*nreachind; i++ )
      {
         g_indfA[*ireach-1].reachind[i] = reachind[i];
      }
   }
   g_indfA[*ireach-1].nfsthabind = *nfsthabind;
   if ( *nfsthabind > 0 )
   {
      g_indfA[*ireach-1].fsthabind = malloc(*nfsthabind*sizeof(int));
      for ( int i=0; i<*nfsthabind; i++ )
      {
         g_indfA[*ireach-1].reachind[i] = fsthabind[i];
      }
   } 
}


void streambugs_delete_fA(int* nreach)
{
   if ( g_nreach > 0 )
   {
      for ( int i=0; i<g_nreach; i++ )
      {
         if ( g_indfA[i].nreachind > 0 )  free(g_indfA[i].reachind);
         if ( g_indfA[i].nfsthabind > 0 ) free(g_indfA[i].fsthabind);
      }
      free(g_indfA);
   }
}


/* ========================================================== */
/* function to update parameters (deSolve)                    */
/* ========================================================== */


void streambugs_rhs_init(void (* odeparms)(int *, double *))
{
   int n = 1;
   odeparms(&n, &g_parms);
}


/* ========================================================== */
/* linear interpolation                                       */
/* ========================================================== */

double linint(int* n, double* x, double* y, double* xout)
{
  /* search indices of points bounding interval containing xout */

  int indleft  = -1;
  int indright = -1;
  for ( int i=0; i<*n; i++ )
  {
    if ( x[i] <= *xout )
    {
      if ( indleft < 0 )
      {
        indleft = i;
      }
      else
      {
         if ( x[i] > x[indleft] ) indleft = i;
      }
    }
    if ( x[i] >= *xout )
    {
      if ( indright < 0 )
      {
        indright = i;
      }
      else
      {
         if ( x[i] < x[indright] ) indright = i;
      }
    }
  }

  /* interpolate within interval */

  double yout = y[0];
  if ( indleft>=0 && indright>=0 )
  {
    if ( x[indleft]!=x[indright] )
    {
      yout = ( (*xout-x[indleft])*y[indright] + (x[indright]-*xout)*y[indleft] )
             / (x[indright]-x[indleft]);
    }
    else
    {
      yout = 0.5 * (y[indleft]+y[indright]);
    }
  }
  else
  {
    if ( indleft<0 && indright>=0 )
    {
       yout = y[indright];
    }
    if ( indleft>=0 && indright<0 ) yout = y[indleft];
  }
  return(yout);
}


/* ========================================================== */
/* polynomial interpolation                                   */
/* ========================================================== */

int find_interval(double* x, int n, double xi)
{
  for (int i = 0; i < n - 1; i++) {
    if (x[i] <= xi && xi < x[i + 1]) {
      return i;
    }
  }
  return n - 2; // handle the case where xi == x[n-1]
}

double polyint(int* n, double* x, double* y, double xout)
{
  int len = *n;
  double* h = (double*)malloc((len - 1) * sizeof(double));
  double* delta = (double*)malloc((len - 1) * sizeof(double));
  double* d = (double*)malloc(len * sizeof(double));
  
  // Compute slopes of the secant lines between consecutive points
  for (int i = 0; i < len - 1; i++) {
    h[i] = x[i + 1] - x[i];
    delta[i] = (y[i + 1] - y[i]) / h[i];
  }
  
  // Compute the tangents at each point
  d[0] = delta[0];
  d[len - 1] = delta[len - 2];
  
  for (int i = 1; i < len - 1; i++) {
    if (delta[i - 1] * delta[i] > 0) {
      d[i] = (2 * h[i] + h[i - 1]) / (3 * (h[i] + h[i - 1])) * delta[i - 1] +
        (2 * h[i - 1] + h[i]) / (3 * (h[i] + h[i - 1])) * delta[i];
    } else {
      d[i] = 0;
    }
  }
  
  // Perform piecewise cubic Hermite interpolation
  int j = find_interval(x, len, xout);
  double h_j = x[j + 1] - x[j];
  double a = (3 * delta[j] - 2 * d[j] - d[j + 1]) / h_j;
  double b = (d[j] + d[j + 1] - 2 * delta[j]) / (h_j * h_j);
  
  double diff = xout - x[j];
  double yout = y[j] + d[j] * diff + a * diff * diff + b * diff * diff * diff;
  
  free(h);
  free(delta);
  free(d);
  
  return yout;
}

/* =========================================================== */
/* functions to calculate trait- and environmental conditions- */
/* dependent modification factors                              */
/* =========================================================== */

double calc_fsapro(int iy)
{
   /* check existence of trait structure and data to interpolate */
   
   int i_x = get_index("saprowqclassval",(*g_parglobalenvtraits).names,(*g_parglobalenvtraits).n);
   if ( i_x < 0 ) return(1);
   int i_y = get_index("saprotolval",(*g_partaxaproptraits).names,(*g_partaxaproptraits).n);
   if ( i_y < 0 ) return(1);
   int n_par = (*(*g_parglobalenvtraits).vectors[i_x]).n_par;
   if ( n_par != (*(*g_partaxaproptraits).matrices[i_y]).n_par ) 
      error("calc_fsapro: x and y vectors to interpolate not of same length");
   int n_row = (*(*g_partaxaproptraits).matrices[i_y]).n_row;
   if ( n_par == 0 ) return(1);

   /* check existence of environmental conditions */
   
   if ( exist_value_matrix("saprowqclass",g_parenvcondhabitat,iy) == 0 ) return(1);    
   double xout = get_value_matrix("saprowqclass",g_parenvcondhabitat,iy);
   
   /* collect data and interpolate */
   
   double* x = (*(*g_parglobalenvtraits).vectors[i_x]).parvals;
   double* y = malloc(n_par*sizeof(double));
   for ( int i=0; i<n_par; i++ ) y[i] = (*(*g_partaxaproptraits).matrices[i_y]).parvals[i*n_row+iy];
   double yout = polyint(&n_par,x,y,&xout);
   if ( yout < 0 ) { free(y); return(1); }
   
   /* transform yout to f */
   
   /* if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left  */
   /* if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right */
   double intercept = get_value_vector("fsapro_intercept",g_parglobal);
   double curv      = get_value_vector("fsapro_curv",g_parglobal);
   double f = 1;
   if ( curv == 0 ) f = intercept - (intercept-1)*yout;
   else             f = intercept - (intercept-1)*(1-exp(-curv*yout))/(1-exp(-curv));

if ( g_debug > 1 )
{
Rprintf("calc_fsapro for state variable %i\n",iy+1);   
for ( int i=0; i<n_par; i++ ) Rprintf("x = %f, y = %f\n",x[i],y[i]);
Rprintf("xout = %f\n",xout);
Rprintf("yout = %f\n",yout);
Rprintf("intercept = %f, curv = %f\n",intercept,curv);
Rprintf("f    = %f\n\n",f);
}
   free(y);
   
   /* return f */   

   return(f);
}

double calc_forgmicropoll(int iy)
{
   /* check existence of trait structure and data to interpolate */
   
   int i_x = get_index("orgmicropollTUval",(*g_parglobalenvtraits).names,(*g_parglobalenvtraits).n);
   if ( i_x < 0 ) return(1);
   int i_y = get_index("orgmicropolltolval",(*g_partaxaproptraits).names,(*g_partaxaproptraits).n);
   if ( i_y < 0 ) return(1);
   int n_par = (*(*g_parglobalenvtraits).vectors[i_x]).n_par;
   if ( n_par != (*(*g_partaxaproptraits).matrices[i_y]).n_par ) 
      error("calc_forgmicropoll: x and y vectors to interpolate not of same length");
   int n_row = (*(*g_partaxaproptraits).matrices[i_y]).n_row;
   if ( n_par == 0 ) return(1);

   /* check existence of environmental conditions */
   
   if ( exist_value_matrix("orgmicropollTU",g_parenvcondhabitat,iy) == 0 ) return(1);    
   double xout = get_value_matrix("orgmicropollTU",g_parenvcondhabitat,iy);
   
   /* collect data and interpolate */
   
   double* x = (*(*g_parglobalenvtraits).vectors[i_x]).parvals;
   double* y = malloc(n_par*sizeof(double));
   for ( int i=0; i<n_par; i++ ) y[i] = (*(*g_partaxaproptraits).matrices[i_y]).parvals[i*n_row+iy];
   double yout = linint(&n_par,x,y,&xout);
   if ( yout < 0 ) { free(y); return(1); }
   
   /* transform yout to f */
   
   /* if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left  */
   /* if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right */
   double intercept = get_value_vector("forgmicropoll_intercept",g_parglobal);
   double curv      = get_value_vector("forgmicropoll_curv",g_parglobal);
   double f = 1;
   if ( curv == 0 ) f = intercept - (intercept-1)*yout;
   else             f = intercept - (intercept-1)*(1-exp(-curv*yout))/(1-exp(-curv));

if ( g_debug > 1 )
{
Rprintf("calc_forgmicropoll for state variable %i\n",iy+1);   
for ( int i=0; i<n_par; i++ ) Rprintf("x = %f, y = %f\n",x[i],y[i]);
Rprintf("xout = %f\n",xout);
Rprintf("yout = %f\n",yout);
Rprintf("intercept = %f, curv = %f\n",intercept,curv);
Rprintf("f    = %f\n\n",f);
}
   free(y);
   
   /* return f */   

   return(f);
}

double calc_fcurrent(int iy)
{
   /* check existence of trait structure and data to interpolate */
   
   int i_x = get_index("currentmsval",(*g_parglobalenvtraits).names,(*g_parglobalenvtraits).n);
   if ( i_x < 0 ) return(1);
   int i_y = get_index("currenttolval",(*g_partaxaproptraits).names,(*g_partaxaproptraits).n);
   if ( i_y < 0 ) return(1);
   int n_par = (*(*g_parglobalenvtraits).vectors[i_x]).n_par;
   if ( n_par != (*(*g_partaxaproptraits).matrices[i_y]).n_par ) 
      error("calc_fcurrent: x and y vectors to interpolate not of same length");
   int n_row = (*(*g_partaxaproptraits).matrices[i_y]).n_row;
   if ( n_par == 0 ) return(1);

   /* check existence tof environmental conditions */
   
   if ( exist_value_matrix("currentms",g_parenvcondhabitat,iy) == 0 ) return(1);    
   double xout = get_value_matrix("currentms",g_parenvcondhabitat,iy);
   
   /* collect data and interpolate */
   
   double* x = (*(*g_parglobalenvtraits).vectors[i_x]).parvals;
   double* y = malloc(n_par*sizeof(double));
   for ( int i=0; i<n_par; i++ ) y[i] = (*(*g_partaxaproptraits).matrices[i_y]).parvals[i*n_row+iy];
   double yout = polyint(&n_par,x,y,&xout);
   if ( yout < 0 ) { free(y); return(1); }
   
   /* transform yout to f */
   
   /* if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left  */
   /* if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right */
   double intercept = get_value_vector("fcurrent_intercept",g_parglobal);
   double curv      = get_value_vector("fcurrent_curv",g_parglobal);
   double f = 1;
   if ( curv == 0 ) f = intercept - (intercept-1)*yout;
   else             f = intercept - (intercept-1)*(1-exp(-curv*yout))/(1-exp(-curv));

if ( g_debug > 1 )
{
Rprintf("calc_fcurrent for state variable %i\n",iy+1);   
for ( int i=0; i<n_par; i++ ) Rprintf("x = %f, y = %f\n",x[i],y[i]);
Rprintf("xout = %f\n",xout);
Rprintf("yout = %f\n",yout);
Rprintf("intercept = %f, curv = %f\n",intercept,curv);
Rprintf("f    = %f\n\n",f);
}
   free(y);
   
   /* return f */   

   return(f);
}

double calc_ftempmax(int iy)
{
   /* check existence of trait structure and data to interpolate */
   
   int i_x = get_index("tempmaxKval",(*g_parglobalenvtraits).names,(*g_parglobalenvtraits).n);
   if ( i_x < 0 ) return(1);
   int i_y = get_index("tempmaxtolval",(*g_partaxaproptraits).names,(*g_partaxaproptraits).n);
   if ( i_y < 0 ) return(1);
   int n_par = (*(*g_parglobalenvtraits).vectors[i_x]).n_par;
   if ( n_par != (*(*g_partaxaproptraits).matrices[i_y]).n_par ) 
      error("calc_ftempmax: x and y vectors to interpolate not of same length");
   int n_row = (*(*g_partaxaproptraits).matrices[i_y]).n_row;
   if ( n_par == 0 ) return(1);

   /* check existence of environmental conditions */
   
   if ( exist_value_matrix("tempmaxK",g_parenvcondhabitat,iy) == 0 ) return(1);    
   double xout = get_value_matrix("tempmaxK",g_parenvcondhabitat,iy);
   
   /* collect data and interpolate */
   
   double* x = (*(*g_parglobalenvtraits).vectors[i_x]).parvals;
   double* y = malloc(n_par*sizeof(double));
   for ( int i=0; i<n_par; i++ ) y[i] = (*(*g_partaxaproptraits).matrices[i_y]).parvals[i*n_row+iy];
   double yout = polyint(&n_par,x,y,&xout);
   if ( yout < 0 ) { free(y); return(1); }
   
   /* transform yout to f */
   
   /* if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left  */
   /* if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right */
   double intercept = get_value_vector("ftempmax_intercept",g_parglobal);
   double curv      = get_value_vector("ftempmax_curv",g_parglobal);
   double f = 1;
   if ( curv == 0 ) f = intercept - (intercept-1)*yout;
   else             f = intercept - (intercept-1)*(1-exp(-curv*yout))/(1-exp(-curv));

if ( g_debug > 1 )
{
Rprintf("calc_ftempmax for state variable %i\n",iy+1);   
for ( int i=0; i<n_par; i++ ) Rprintf("x = %f, y = %f\n",x[i],y[i]);
Rprintf("xout = %f\n",xout);
Rprintf("yout = %f\n",yout);
Rprintf("intercept = %f, curv = %f\n",intercept,curv);
Rprintf("f    = %f\n\n",f);
}
   free(y);
   
   /* return f */   

   return(f);
}

double calc_fmicrohab(int iy)
{
   /* check existence of trait structure and area fractions */
   
   int i_af = get_index("microhabaf",(*g_parenvcondhabitatgroup).names,(*g_parenvcondhabitatgroup).n);
   if ( i_af < 0 ) return(1);   
   int i_sc = get_index("microhabtolval",(*g_partaxaproptraits).names,(*g_partaxaproptraits).n);
   if ( i_sc < 0 ) return(1);   
   int n_par = (*(*g_parenvcondhabitatgroup).matrices[i_af]).n_par;   
   if ( n_par != (*(*g_partaxaproptraits).matrices[i_sc]).n_par ) 
      error("calc_fmicrohab: area fractions and scores not of same length");
   if ( n_par == 0 ) return(1);
   int n_row = (*(*g_parenvcondhabitatgroup).matrices[i_af]).n_row;   
   if ( n_row != (*(*g_partaxaproptraits).matrices[i_sc]).n_row ) 
      error("calc_fmicrohab: area fractions and scores matrices have different numbers of rows");
      
   /* calculate weighted average of transformed scores */
   
   double intercept = get_value_vector("fmicrohab_intercept",g_parglobal);
   double curv      = get_value_vector("fmicrohab_curv",g_parglobal);
   double f = 0;
   for ( int i=0; i<n_par; i++ )
   {
      /* get area fraction and score */
      
      double af = (*(*g_parenvcondhabitatgroup).matrices[i_af]).parvals[i*n_row+iy];
      double sc = (*(*g_partaxaproptraits).matrices[i_sc]).parvals[i*n_row+iy];

      /* transform score */
      
      /* if curv > 0 and intercept <1: function is curved to the right, if curv < 0 and intercept <1 function is curved to the left  */
      /* if curv > 0 and intercept >1: function is curved to the left,  if curv < 0 and intercept >1 function is curved to the right */
      if ( curv == 0 ) sc = intercept - (intercept-1)*sc;
      else             sc = intercept - (intercept-1)*(1-exp(-curv*sc))/(1-exp(-curv));

      /* add to f */
      
      f = f + af*sc;
   }

if ( g_debug > 1 )
{
Rprintf("calc_fmicrohab for state variable %i\n",iy+1);   
for ( int i=0; i<n_par; i++ ) Rprintf("af = %f, sc = %f\n",
                                      (*(*g_parenvcondhabitatgroup).matrices[i_af]).parvals[i*n_row+iy],
                                      (*(*g_partaxaproptraits).matrices[i_sc]).parvals[i*n_row+iy]);
Rprintf("intercept = %f, curv = %f\n",intercept,curv);
Rprintf("f  = %f\n\n",f);
}
       
   return(f);
}

/* ========================================================== */
/* function to update time-dependent parameters               */
/* ========================================================== */

void streambugs_updatepar(double* t)
{
  double* inpvals = 0;
  if ( g_ninp > 0 )
  {
     /* interpolate input */
     
     inpvals = malloc(g_ninp*sizeof(double));
     for ( int i=0; i<g_ninp; i++ )
     {
        inpvals[i] = polyint(&g_inp[i].n, 
                            g_inp[i].x,
                            g_inp[i].y,
                            t);
     } 
     
     /* update global parameters */

     int n_inp = (*g_parglobal).n_inp;
     if ( n_inp > 0 ) 
     {
        for ( int k=0; k<n_inp; k++ )
        {
           int ind = (*g_parglobal).inpinds[n_inp+k];
           (*g_parglobal).parvals[ind-1] 
              = inpvals[(*g_parglobal).inpinds[k]-1];
        }
     }
     
     /* update trait-dependent global parameters */

     if ( (*g_parglobalenvtraits).n > 0 )
     {
        for ( int i=0; i<(*g_parglobalenvtraits).n; i++ )
        {
           int n_inp = (*((*g_parglobalenvtraits).vectors[i])).n_inp;
           if ( n_inp > 0 ) 
           {
              for ( int k=0; k<n_inp; k++ )
              {
                 int ind = (*((*g_parglobalenvtraits).vectors[i])).inpinds[n_inp+k];
                 (*((*g_parglobalenvtraits).vectors[i])).parvals[ind-1] 
                    = inpvals[(*((*g_parglobalenvtraits).vectors[i])).inpinds[k]-1];
              }
           }
        }
     }
          
     /* update reach-dependent environmental conditions */
     
     n_inp = (*g_parenvcondreach).n_inp;
     if ( n_inp > 0 ) 
     {
        int n_row = (*g_parenvcondreach).n_row;
        for ( int k=0; k<n_inp; k++ )
        {
           int i = (*g_parenvcondreach).inpinds[n_inp+k];
           int j = (*g_parenvcondreach).inpinds[2*n_inp+k];
           (*g_parenvcondreach).parvals[(j-1)*n_row+i-1] 
              = inpvals[(*g_parenvcondreach).inpinds[k]-1];
        }
     }
     
     /* update habitat-dependent environmental conditions */
     
     n_inp = (*g_parenvcondhabitat).n_inp;
     if ( n_inp > 0 ) 
     {
        int n_row = (*g_parenvcondhabitat).n_row;
        for ( int k=0; k<n_inp; k++ )
        {
           int i = (*g_parenvcondhabitat).inpinds[n_inp+k];
           int j = (*g_parenvcondhabitat).inpinds[2*n_inp+k];
           (*g_parenvcondhabitat).parvals[(j-1)*n_row+i-1] 
              = inpvals[(*g_parenvcondhabitat).inpinds[k]-1];
        }
     }
     /* normalization of fA */
     int indfA = get_index("fA",(*g_parenvcondhabitat).parnames,
                                (*g_parenvcondhabitat).n_par);
     if ( indfA < 0 ) error("streambugs_rhs: fA not found");
     double* fA = malloc(g_ny*sizeof(double));
     for ( int i=0; i<g_ny; i++ ) fA[i] = (*g_parenvcondhabitat).parvals[indfA*g_ny+i];
     for ( int i=0; i<g_nreach; i++ )
     {
        double f = 0;
        for ( int j=0; j<g_indfA[i].nfsthabind; j++ )
        {
           f = f + fA[g_indfA[i].reachind[g_indfA[i].fsthabind[j]-1]-1];
        } 
        f = 1/f;
        for ( int j=0; j<g_indfA[i].nreachind; j++ )
        {
           int ind = g_indfA[i].reachind[j];
           (*g_parenvcondhabitat).parvals[(indfA-1)*g_ny+ind-1]
              = f * fA[ind-1];
        }
     }
     free(fA);
     
     /* update habitat-dependent environmental conditions belonging to a group */
    
     if ( (*g_parenvcondhabitatgroup).n > 0 )
     {
        for ( int l=0; l<(*g_parenvcondhabitatgroup).n; l++ )
        {
           n_inp = (*((*g_parenvcondhabitatgroup).matrices[l])).n_inp;
           if ( n_inp > 0 ) 
           {
              int n_row = (*((*g_parenvcondhabitatgroup).matrices[l])).n_row;
              for ( int k=0; k<n_inp; k++ )
              {
                 int i = (*((*g_parenvcondhabitatgroup).matrices[l])).inpinds[n_inp+k];
                 int j = (*((*g_parenvcondhabitatgroup).matrices[l])).inpinds[2*n_inp+k];
                 (*((*g_parenvcondhabitatgroup).matrices[l])).parvals[(j-1)*n_row+i-1] 
                    = inpvals[(*((*g_parenvcondhabitatgroup).matrices[l])).inpinds[k]-1];
              }
           }
        }
     }

     /* update initial conditions */
     
     n_inp = (*g_parinitcond).n_inp;
     if ( n_inp > 0 ) 
     {
        int n_row = (*g_parinitcond).n_row;
        for ( int k=0; k<n_inp; k++ )
        {
           int i = (*g_parinitcond).inpinds[n_inp+k];
           int j = (*g_parinitcond).inpinds[2*n_inp+k];
           (*g_parinitcond).parvals[(j-1)*n_row+i-1] 
              = inpvals[(*g_parinitcond).inpinds[k]-1];
        }
     }

     /* update inputs */
     
     n_inp = (*g_parinput).n_inp;
     if ( n_inp > 0 ) 
     {
        int n_row = (*g_parinput).n_row;
        for ( int k=0; k<n_inp; k++ )
        {
           int i = (*g_parinput).inpinds[n_inp+k];
           int j = (*g_parinput).inpinds[2*n_inp+k];
           (*g_parinput).parvals[(j-1)*n_row+i-1] 
              = inpvals[(*g_parinput).inpinds[k]-1];
        }
     }

     /* update directly defined taxa properties */
     
     n_inp = (*g_partaxapropdirect).n_inp;
     if ( n_inp > 0 ) 
     {
        int n_row = (*g_partaxapropdirect).n_row;
        for ( int k=0; k<n_inp; k++ )
        {
           int i = (*g_partaxapropdirect).inpinds[n_inp+k];
           int j = (*g_partaxapropdirect).inpinds[2*n_inp+k];
           (*g_partaxapropdirect).parvals[(j-1)*n_row+i-1] 
              = inpvals[(*g_partaxapropdirect).inpinds[k]-1];
        }
     }

     /* update trait-dependent taxa properties */
     
     if ( (*g_partaxaproptraits).n > 0 )
     {
        for ( int l=0; l<(*g_partaxaproptraits).n; l++ )
        {
           n_inp = (*((*g_partaxaproptraits).matrices[l])).n_inp;
           if ( n_inp > 0 ) 
           {
              int n_row = (*((*g_partaxaproptraits).matrices[l])).n_row;
              for ( int k=0; k<n_inp; k++ )
              {
                 int i = (*((*g_partaxaproptraits).matrices[l])).inpinds[n_inp+k];
                 int j = (*((*g_partaxaproptraits).matrices[l])).inpinds[2*n_inp+k];
                 (*((*g_partaxaproptraits).matrices[l])).parvals[(j-1)*n_row+i-1] 
                    = inpvals[(*((*g_partaxaproptraits).matrices[l])).inpinds[k]-1];
              }
           }
        }
     }

     /* update kinetic parameters of taxon-based processes */
     
     for ( int i=0; i<g_ny; i++ )
     {
        if ( g_proctaxon[i].n > 0 )
        {
           for ( int j=0; j<g_proctaxon[i].n; j++ )
           {
              n_inp = (*g_proctaxon[i].procs[j].par).n_inp;
              if ( n_inp > 0 ) 
              {
                 for ( int k=0; k<n_inp; k++ )
                 {
                    int ind = (*g_proctaxon[i].procs[j].par).inpinds[n_inp+k];
                    (*g_proctaxon[i].procs[j].par).parvals[ind-1] 
                       = inpvals[(*g_proctaxon[i].procs[j].par).inpinds[k]-1];
                 }
              }
           }
        }
     }

     /* update kinetic parameters of taxon-based processes */
     
     for ( int i=0; i<g_ny; i++ )
     {
        if ( g_procweb[i].n > 0 )
        {
           for ( int j=0; j<g_procweb[i].n; j++ )
           {
              n_inp = (*g_procweb[i].webprocs[j].par).n_inp;
              if ( n_inp > 0 ) 
              {
                 for ( int k=0; k<n_inp; k++ )
                 {
                    int ind = (*g_procweb[i].webprocs[j].par).inpinds[n_inp+k];
                    (*g_procweb[i].webprocs[j].par).parvals[ind-1] 
                       = inpvals[(*g_procweb[i].webprocs[j].par).inpinds[k]-1];
                 }
              }
              if ( g_procweb[i].webprocs[j].n > 0 )
              {
                 for ( int k=0; k<g_procweb[i].webprocs[j].n; k++ )
                 {
                    n_inp = (*g_procweb[i].webprocs[j].procs[k].par).n_inp;
                    if ( n_inp > 0 ) 
                    {
                       for ( int l=0; l<n_inp; l++ )
                       {
                          int ind = (*g_procweb[i].webprocs[j].procs[k].par).inpinds[n_inp+l];
                          (*g_procweb[i].webprocs[j].procs[k].par).parvals[ind-1] 
                             = inpvals[(*g_procweb[i].webprocs[j].procs[k].par).inpinds[l]-1];
                       }
                    }
                 }
              }
           }
        }
     }
     if ( g_ninp > 0 ) free(inpvals);
  }   
}


/* ========================================================== */
/* function to calculate right-hand side of odes of the       */
/* streambugs model                                           */
/* ========================================================== */

void streambugs_rhs(int* neq, double* t, double* y, double* ydot, 
                    double* yout, int* ip)
{
  if ( ip[0]<1 ) error("nout should be at least 1");
  
  /* copy state variables */
  
  double* y_original = malloc(g_ny*sizeof(double));
  for ( int i=0; i<g_ny; i++ )
  {
     y_original[i] = y[i];
     if ( y[i] < 0 ) y[i] = 0;
  }
  
  /* update parameters */
  
  streambugs_updatepar(t);
    
  /* get global parameter values */

  double kBoltzmann   = get_value_vector("kBoltzmann",g_parglobal);
  double M0           = get_value_vector("M0",g_parglobal);
  
  /* initialize right hand side */

  for ( int i=0; i<g_ny; i++ ) ydot[i] = 0;

  /* loop over state variables to add rates */

  for ( int i=0; i<g_ny; i++ )
  {
     /* get parameter values for state variable i */
     
     double w         = get_value_matrix("w",g_parenvcondreach,i);
     double fA        = get_value_matrix("fA",g_parenvcondhabitat,i);
     double T         = get_value_matrix("T",g_parenvcondhabitat,i);
          
     double i0        = get_value_matrix("i0",g_partaxapropdirect,i);
     double M         = get_value_matrix("M",g_partaxapropdirect,i);
     double b         = get_value_matrix("b",g_partaxapropdirect,i);
     double Ea        = get_value_matrix("Ea",g_partaxapropdirect,i);
     double EC        = get_value_matrix("EC",g_partaxapropdirect,i);
     double fbasaltax = get_value_matrix("fbasaltax",g_partaxapropdirect,i);

     double Input     = get_value_matrix("Input",g_parinput,i);

     /* add Input to rhs */
          
     ydot[i] = ydot[i] + Input * w*fA;
     
     /* calculate basal metabolic rate */
     
     double r_basal = i0 * fbasaltax
                      * pow(M,b-1)/pow(M0,b) * exp( -Ea/(kBoltzmann*T) ) / EC
                      * y[i]/(w*fA);
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, w = %f, fA = %f\n",i+1,w,fA); 
Rprintf("i = %i, r_basal = %f\n",i+1,r_basal); 
}
}
          
     /* loop over taxon-based processes */
     
     if ( g_proctaxon[i].n > 0 )
     {
        for ( int j=0; j<g_proctaxon[i].n; j++ )
        {
           /* process "Miner" */
           
           if ( strcmp("Miner",g_proctaxon[i].procs[j].name) == 0 )
           {
              double kminer = get_value_vector("kminer",g_proctaxon[i].procs[j].par);  

              double r_miner = kminer * y[i]/(w*fA);
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_miner = %f\n",i+1,r_miner); 
}
}
              
              if ( (*g_proctaxon[i].procs[j].stoich).n > 0 )
              {
                 for ( int k=0; k<(*g_proctaxon[i].procs[j].stoich).n; k++ )
                 {
                    int    indy;
                    double coeff;
                    get_stoich((*g_proctaxon[i].procs[j].stoich).names[k],
                               g_proctaxon[i].procs[j].stoich,
                               &indy,&coeff);
                    ydot[indy-1] = ydot[indy-1] + coeff * r_miner * (w*fA);
                 }     
              }
           }
           else
           {
               
           /* process "Drift" */
           
           if ( strcmp("Drift",g_proctaxon[i].procs[j].name) == 0 )
           {
              double cdet     = get_value_vector("cdet",g_proctaxon[i].procs[j].par);  
              
              double tau      = get_value_matrix("tau",g_parenvcondhabitat,i);
              double taucrit  = get_value_matrix("taucrit",g_parenvcondhabitat,i);
              double r_drift  = 0;

              if (tau >= taucrit)
              {
                r_drift = cdet * y[i]/(w*fA) * pow((tau - taucrit),2);
              } 
           
              
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_drift = %f\n",i+1,r_drift); 
}
}

              if ( (*g_proctaxon[i].procs[j].stoich).n > 0 )
              {
                for ( int k=0; k<(*g_proctaxon[i].procs[j].stoich).n; k++ )
                {
                  int    indy;
                  double coeff;
                  get_stoich((*g_proctaxon[i].procs[j].stoich).names[k],
                              g_proctaxon[i].procs[j].stoich,
                              &indy,&coeff);
                   ydot[indy-1] = ydot[indy-1] + coeff * r_drift * (w*fA);
                }     
              }
           }
           else
           {
               
           /* process "Resp" */
           
           if ( strcmp("Resp",g_proctaxon[i].procs[j].name) == 0 )
           {
              double fresp = get_value_vector("fresp",g_proctaxon[i].procs[j].par);  

              double r_resp = fresp * r_basal;
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_resp = %f\n",i+1,r_resp); 
}
}

              if ( (*g_proctaxon[i].procs[j].stoich).n > 0 )
              {
                 for ( int k=0; k<(*g_proctaxon[i].procs[j].stoich).n; k++ )
                 {
                    int    indy;
                    double coeff;
                    get_stoich((*g_proctaxon[i].procs[j].stoich).names[k],
                               g_proctaxon[i].procs[j].stoich,
                               &indy,&coeff);
                    ydot[indy-1] = ydot[indy-1] + coeff * r_resp * (w*fA);
                 }     
              }
           }
           else
           {
               
           /* process "Death" */
           
           if ( strcmp("Death",g_proctaxon[i].procs[j].name) == 0 )
           {
              double fdeath = get_value_vector("fdeath",g_proctaxon[i].procs[j].par); 
              
              double fsapro     = calc_fsapro(i);
              double forgmicropoll = calc_forgmicropoll(i);
              
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, fsapro = %f, forgmicropoll = %f\n",i+1,fsapro,forgmicropoll);
}
}
              double r_death = fdeath * r_basal * fsapro * forgmicropoll;
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_death = %f\n",i+1,r_death); 
}
}
              
              if ( (*g_proctaxon[i].procs[j].stoich).n > 0 )
              {
                 for ( int k=0; k<(*g_proctaxon[i].procs[j].stoich).n; k++ )
                 {
                    int    indy;
                    double coeff;
                    get_stoich((*g_proctaxon[i].procs[j].stoich).names[k],
                               g_proctaxon[i].procs[j].stoich,
                               &indy,&coeff);
                    ydot[indy-1] = ydot[indy-1] + coeff * r_death * (w*fA);
                 }     
              }
           }
           else
           {
               
           /* process "Prod" */
           
           if ( strcmp("Prod",g_proctaxon[i].procs[j].name) == 0 )
           {
              double fprod   = get_value_vector("fprod",g_proctaxon[i].procs[j].par);  
              double fgrotax = get_value_vector("fgrotax",g_proctaxon[i].procs[j].par);  
              double KI      = get_value_vector("KI",g_proctaxon[i].procs[j].par);  
              double KP      = get_value_vector("KP",g_proctaxon[i].procs[j].par);  
              double KN      = get_value_vector("KN",g_proctaxon[i].procs[j].par);  
              double hdens   = get_value_vector("hdens",g_proctaxon[i].procs[j].par);  

              double I0      = get_value_matrix("I0",g_parenvcondhabitat,i);  
              double CP      = get_value_matrix("CP",g_parenvcondhabitat,i);  
              double CN      = get_value_matrix("CN",g_parenvcondhabitat,i);  
              double fshade  = get_value_matrix("fshade",g_parenvcondhabitat,i);  
              
              double r_prod = fprod * fgrotax
                              * I0/(KI+I0)
                              * dmin(CP/(KP+CP),CN/(KN+CN))
                              * hdens/(hdens+y[i]/(w*fA))
                              * (1-fshade)
                              * r_basal;
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_prod = %f\n",i+1,r_prod); 
}
}
              
              if ( (*g_proctaxon[i].procs[j].stoich).n > 0 )
              {
                 for ( int k=0; k<(*g_proctaxon[i].procs[j].stoich).n; k++ )
                 {
                    int    indy;
                    double coeff;
                    get_stoich((*g_proctaxon[i].procs[j].stoich).names[k],
                               g_proctaxon[i].procs[j].stoich,
                               &indy,&coeff);
                    ydot[indy-1] = ydot[indy-1] + coeff * r_prod * (w*fA);
                 }     
              }
           }
           else
           {
              sprintf("Unknown process: %s",g_proctaxon[i].procs[j].name);
              error(buffer);    
           }
           } 
           }              
           }
           }
        }
     }       
     
     /* loop over food web processes */
     
     if ( g_procweb[i].n > 0 )
     {
        for ( int j=0; j<g_procweb[i].n; j++ )
        {
           /* process "Cons" */
           
           if ( strcmp("Cons",g_procweb[i].webprocs[j].name) == 0 )
           {
              /* get parameter values */
              
              double fcons   = get_value_vector("fcons",g_procweb[i].webprocs[j].par);  
              double hdens   = get_value_vector("hdens",g_procweb[i].webprocs[j].par);  
              double Kfood   = get_value_vector("Kfood",g_procweb[i].webprocs[j].par);  
              double q       = get_value_vector("q",g_procweb[i].webprocs[j].par);  
              double fgrotax = get_value_vector("fgrotax",g_procweb[i].webprocs[j].par);
                
              double fcurrent  = calc_fcurrent(i);
              double ftempmax     = calc_ftempmax(i);
              double fmicrohab = calc_fmicrohab(i);
              
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, fcurrent = %f, ftempmax = %f, fmicrohab = %f\n",
        i+1,fcurrent,ftempmax,fmicrohab);
}
}
            
              double Kdens = hdens * fcurrent * ftempmax * fmicrohab;
              
              /* get preferences and food densities */
              
              int n = g_procweb[i].webprocs[j].n;
              double* yfood = malloc(n*sizeof(double));
              double* pref  = malloc(n*sizeof(double));
              for ( int k=0; k<n; k++ )
              {
                 if ( strcmp(g_procweb[i].webprocs[j].procs[k].name,"SusPOM") == 0 )
                 {
                    double DSusPOM = get_value_matrix("DSusPOM",g_parenvcondhabitat,i);
                    yfood[k] = DSusPOM * (w*fA);
                 }
                 else
                 {
                    int ind_food_stoich = get_index(g_procweb[i].webprocs[j].procs[k].name,
                                                    (*g_procweb[i].webprocs[j].procs[k].stoich).names,
                                                    (*g_procweb[i].webprocs[j].procs[k].stoich).n);
                    if ( ind_food_stoich < 0 )
                    {
                       sprintf(buffer,"Cons process of state var %i on %s: food not found in stoich",
                               i+1,g_procweb[i].webprocs[j].procs[k].name);   
                       error(buffer);
                    }
                    int ind_food = (*g_procweb[i].webprocs[j].procs[k].stoich).inds[ind_food_stoich];
                    yfood[k] = y[ind_food-1];
                 }
                 pref[k] = get_value_vector("Pref",g_procweb[i].webprocs[j].procs[k].par);
              }
              double sum_food = 0;
              double sum_food_pref = 0;
              for ( int k=0; k<n; k++ )
              {
                 sum_food = sum_food + yfood[k];
                 sum_food_pref = sum_food_pref + yfood[k]*pref[k];
              }
              
              /*double ffoodlim = pow(sum_food/(w*fA),q) / ( pow(Kfood,q) + pow(sum_food/(w*fA),q) );*/
              
              double ffoodlim;
              
              if(sum_food >= 0)
              {
                ffoodlim = pow(sum_food/(w*fA),q) / ( pow(Kfood,q) + pow(sum_food/(w*fA),q) );
              } 
              else
              {
                ffoodlim = pow(sum_food/(w*fA),q) / ( pow(Kfood,q) );
              }
              
              /* consumption rate without preference factor */
              
              double r_cons_tot = fcons * fgrotax
                                  * Kdens / (Kdens + y[i]/(w*fA))
                                  * ffoodlim
                                  * r_basal;
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_cons_tot = %f, sum_food = %f, sum_food_pref = %f\n",
        i+1,r_cons_tot,sum_food,sum_food_pref); 
}
}
                        
              /* loop over all food types */

              for ( int k=0; k<n; k++ )
              {
                 /* calculate rate */
                 
                 double r_cons = 0;
                 if(sum_food_pref >  0) 
                 {
                   r_cons = r_cons_tot * yfood[k]*pref[k] / sum_food_pref;
                 } 
                 

                 /* calculate contributions to rate of change of all affected taxa */          
                
                 if ( (*g_procweb[i].webprocs[j].procs[k].stoich).n > 0 )
                 {
                    for ( int l=0; l<(*g_procweb[i].webprocs[j].procs[k].stoich).n; l++ )
                    {
                       int    indy;
                       double coeff;
                       get_stoich((*g_procweb[i].webprocs[j].procs[k].stoich).names[l],
                                  g_procweb[i].webprocs[j].procs[k].stoich,
                                  &indy,&coeff);
                       ydot[indy-1] = ydot[indy-1] + coeff * r_cons * (w*fA);
                    }
                 }
              }
              
              free(yfood);
              free(pref);
           }
           else
           {
               
           /* process "FishPred" */
           
           if ( strcmp("FishPred",g_procweb[i].webprocs[j].name) == 0 )
           {
              /* get parameter values */

              double DFish   = get_value_matrix("DFish",g_parenvcondhabitat,i);
              
              double cfish   = get_value_vector("cfish",g_procweb[i].webprocs[j].par);  
              double Kfood   = get_value_vector("Kfood",g_procweb[i].webprocs[j].par);  
              double q       = get_value_vector("q",g_procweb[i].webprocs[j].par);  
              
              /* get preferences and food densities */
              
              int n = g_procweb[i].webprocs[j].n;
              double* yfood = malloc(n*sizeof(double));
              double* pref  = malloc(n*sizeof(double));
              for ( int k=0; k<n; k++ )
              {
                 if ( strcmp(g_procweb[i].webprocs[j].procs[k].name,"SusPOM") == 0 )
                 {
                    double DSusPOM = get_value_matrix("DSusPOM",g_parenvcondhabitat,i);
                    yfood[k] = DSusPOM * (w*fA);
                 }
                 else
                 {
                    int ind_food_stoich = get_index(g_procweb[i].webprocs[j].procs[k].name,
                                                    (*g_procweb[i].webprocs[j].procs[k].stoich).names,
                                                    (*g_procweb[i].webprocs[j].procs[k].stoich).n);
                    if ( ind_food_stoich < 0 )
                    {
                       sprintf(buffer,"Cons process of state var %i on %s: food not found in stoich",
                               i+1,g_procweb[i].webprocs[j].procs[k].name);   
                       error(buffer);
                    }
                    int ind_food = (*g_procweb[i].webprocs[j].procs[k].stoich).inds[ind_food_stoich];
                    yfood[k] = y[ind_food-1];
                 }
                 pref[k] = get_value_vector("Pref",g_procweb[i].webprocs[j].procs[k].par);
              }
              double sum_food = 0;
              double sum_food_pref = 0;
              for ( int k=0; k<n; k++ )
              {
                 sum_food = sum_food + yfood[k];
                 sum_food_pref = sum_food_pref + yfood[k]*pref[k];
              }
              double ffoodlim = pow(sum_food/(w*fA),q) / ( pow(Kfood,q) + pow(sum_food/(w*fA),q) );
              
              /* consumption rate without preference factor */
              
              double r_fishpred_tot = DFish/10000 * cfish*365.25 * ffoodlim;
if ( g_debug > 0 )
{
if ( *t == 0 )
{
Rprintf("i = %i, r_fishpred_tot = %f, sum_food = %f, sum_food_pref = %f\n",
        i+1,r_fishpred_tot,sum_food,sum_food_pref); 
}
}
                        
              /* loop over all food types */

              for ( int k=0; k<n; k++ )
              {
                 /* calculate rate */
                 
                 double r_fishpred = 0;
                 if(sum_food_pref >  0) 
                 {
                   r_fishpred = r_fishpred_tot * yfood[k]*pref[k] / sum_food_pref;
                 } 
                 
                 /* calculate contributions to rate of change of all affected taxa */          
                
                 if ( (*g_procweb[i].webprocs[j].procs[k].stoich).n > 0 )
                 {
                    for ( int l=0; l<(*g_procweb[i].webprocs[j].procs[k].stoich).n; l++ )
                    {
                       int    indy;
                       double coeff;
                       get_stoich((*g_procweb[i].webprocs[j].procs[k].stoich).names[l],
                                  g_procweb[i].webprocs[j].procs[k].stoich,
                                  &indy,&coeff);
                       ydot[indy-1] = ydot[indy-1] + coeff * r_fishpred * (w*fA);
                    }
                 }
              }
              
              free(yfood);
              free(pref);
           }
           else
           {
              sprintf(buffer,"Unknown process: %s",g_procweb[i].webprocs[j].name);
              error(buffer);  
           }
           }
        }
     }
  }
if ( g_debug > 0 )
{
if ( *t == 0 )
{
for ( int i=0; i<g_ny; i++ )
{
Rprintf("i = %i, ydot = %f\n",i+1,ydot[i]); 
}     
}
}
  for ( int i=0; i<g_ny; i++ )
  {
     if ( y_original[i] < 0 ) ydot[i] = - 5 * y_original[i];
     y[i] = y_original[i];
  }
  free(y_original);
}

