/* Find small drifters against a stable background.   3/20/97 */

#include <stdio.h>

#define COUNT   0   /* If true, then periodically print how often */
                    /* various functions are called.              */

/************************************************************************/
/* These control how big things can be.                                 */
/************************************************************************/

#define MAXHT       81      /* Height and width of space */
#define MAXWD       81
int HT = MAXHT;             /* May be reduced by user */
int WD = MAXWD;

#define MAXGEN      500     /* Max generation that can be computed */
#define CHGLISTLTH  100000  /* Max total number of changed cells and    */
                            /* their neighbors in all gens            */
#define HASHTBLSIZE 2000000 /* Max number of distinct life histories */
#define MAXKNOWN    5000    /* Max number of known rotors */
#define MAXFILESIZE 1000000 /* Max length of known rotor descriptors & names */
#define NUMVARS     200     /* # of vars available for program modification */
#define MAXROTORDESCLTH	1000	/* Max length of rotor descriptors */

/************************************************************************/
/* Various definitions                                                  */
/************************************************************************/

typedef char boolean;
#define FALSE   0
#define TRUE    1

#define OK      FALSE   /* Return values for functions */
#define ERR     TRUE

/* Values of elements of bkgd, curr, cell */
#define OFF 0
#define ON  1
#define UNK 16

/* Additional values used in cell during rotor analysis */
#define STATOR  2
#define CONN    3

/* Names for symmetries */
#define NOSYMM      0
#define HORSYMM     1
#define DIAGSYMM    2
#define ROT90SYMM   3
#define ROT180SYMM  4
#define PLUSSYMM    5
#define XSYMM       6
#define FULLSYMM    7
#define VERTSYMM    8
int SYMM = NOSYMM;

#define err(msg)        { printf("\n\n%s\n",msg); exit(1); }
#define err1(msg,arg)   { printf("\n\n"); printf(msg,arg); printf("\n"); \
                          exit(1); }
#define err2(msg,arg0,arg1) { printf("\n\n"); printf(msg,arg0,arg1); \
                              printf("\n"); exit(1); }

#define min(x,y)    ((x)<(y) ? (x) : (y))
#define max(x,y)    ((x)<(y) ? (y) : (x))
#define diff(a,b)   ((a)<(b) ? (b)-(a) : (a)-(b))

#if COUNT       /* Variables used to count function calls */
int TRACE = 0;
long int countsetbkgd, countconsistify, countconsis9, countproceed,
         countbackup, countgo, countlistneighbors, countchangecurr,
         countnxgen, counttrycell;
#endif
long int countcomporbackuplo,     /* Up to 10^6 */
         countcomporbackuphi;     /* Always do this one */

boolean SKIPSTABLE = FALSE; /* If true, don't display stable outcomes */
boolean NOPICS = FALSE;     /* If true, don't show pictures */
boolean SKIPFIZZLE = FALSE; /* If true, don't display fizzle outcomes */
boolean SHOWFIN = FALSE;    /* If true, display finished patterns */
boolean SHOWALL = FALSE;    /* If true, display all gens */

char knownrotorsfilename[100] = "/Users/doris/Drift/knownrotors";

/************************************************************************/
/* These are the arrays that contain information about the rule.        */
/************************************************************************/

boolean rule[2][9] = {{0,0,0,1,0,0,0,0,0},  /* Birth & survival rules */
                      {0,0,1,1,0,0,0,0,0}};

/* -------------------------------------------------------------------- */

char transtable[2][129];

/* If val is the value of a cell in the current generation (OFF or ON)
   and nbhd is the sum of the values of the 8 surrounding cells (each of
   which may be OFF, ON, or UNK), then transtable[val][nbhd] is the value
   of the cell in the next generation (OFF, ON, or UNK).
*/

/* -------------------------------------------------------------------- */

enum {NOINFO = 0, INCONSIS, NAYSOFF, NAYSON, CELLOFF, CELLON};

char consistable[UNK+1][129];
char oldconsistable[UNK+1][129];

/* If val is the value of a background cell (OFF, ON, or UNK) and
   nbhd is the sum of the values of the 8 surrounding cells, then
   consistable[val][nbhd] is:

   INCONSIS if the values are inconsistent;
   NAYSOFF  if all neighboring UNK cells must be OFF;
   NAYSON   if all neighboring UNK cells must be ON;
   CELLOFF  if the cell is UNK and must be OFF;
   CELLON   if the cell is UNK and must be ON;
   NOINFO   if nothing is forced.
   
   Note that for some rules and values of val and nbhd, it's possible
   for both the cell and its neighbors to be forced.  For example, if
   the rule is B78/S8, val=UNK, oncount=7, offcount=0, and unkcount=1,
   then the cell and its unknown neighbor must be ON.  In such cases,
   the value of consistable[val][nbhd] will be CELLOFF or CELLON; after
   consistify(r,c) sets an UNK cell to OFF or ON, it checks again to
   see if any neighbors are forced.
*/

/************************************************************************/
/* Here are the global variables used in the search.                    */
/************************************************************************/

unsigned char bkgd[MAXHT][MAXWD],       /* Each cell is ON, OFF, or UNK */
              curr[MAXHT][MAXWD];

int	naysum[MAXHT][MAXWD];				/* Sum of 8 neighbors in bkgd */

#define DONTCHANGE  1
#define DONTCOUNT   2
unsigned char flag[MAXHT][MAXWD];
/* Each bit is a flag controlling some aspect of the cell.  At the moment
   only 2 are defined:
   
   Bit 0 = DONTCHANGE   If true, don't allow the cell to change.
   Bit 1 = DONTCOUNT    Don't count changes here in computing the size
                        of the changed region.
*/

typedef struct
  { int row;
    int col;
  }  point;

int gen,                /* Current generation number */
    maxgenreached=0;    /* Largest generation computed */
int maxchng=9;          /* Max # of changed cells in any generation */
int maxwidth=3;         /* Max width of set of changed cells in any gen */
int maxheight=3;        /* Max height of set of changed cells in any gen */
int prob = 50;          /* 100*probability, used when freely setting cells */
boolean found;
point *nay,             /* Ptr into list of neighbors of previous gen */
      *chg;             /* Ptr to end of list of changes in current gen */

point chglist[CHGLISTLTH];  /* Concatenation of lists of changes and
                               neighbors in all gens */
point *(chgd[MAXGEN+1]),
      *(nays[MAXGEN+1]);
/* Ptrs into chglist:
        chglist == chgd[0] < nays[0] < chgd[1] < nays[1] < ...
   chgd[g] points to list of cells that differ from bkgd in gen g.
   nays[g] points to list of neighbors of such.
   chgd[g] is subset of nays[g-1].
*/

int chgcount[MAXGEN];	/* Number of changed cells in each gen */
int agesm[MAXGEN];		/* Age sums.  Only computed if var[129] nonzero */
int width[MAXGEN], height[MAXGEN];	/* Width and height of changed region */

int numcc;					/* Required values of chgcount, specified by */
int reqchgcount[MAXGEN];	/* 'cc' command, are in reqchgcount[0], ..., */
							/* reqchgcount[numcc-1].                     */

typedef struct
  { int row;
    int col;
    unsigned char val;
    boolean free;
    int gen;
    point *nay;
    point *chg;
  }  setting;

setting settinglist[MAXHT*MAXWD];       /* List of background cell settings */

setting *nxstng,        /* Pointer to setting whose consequences are
                           being examined */
        *nwstng;        /* Pointer to setting that's being added to list */

long int var[NUMVARS];  /* Variables available for program modifications */
/* Currently defined:
   var[100]   If nonzero (>= 2), only require change count <= maxchng at
			  least once in every var[100] consecutive gens.
   var[101]   If nonzero, bound on # of cells that are changed this gen
			  but not changed last gen.
   var[102]   If nonzero, bound on diagonal width (NE to SW).
   var[103]   If nonzero, bound on diagonal width (NW to SE).
   var[104]   If nonzero, force pattern to move right at >= 2c/3.
   var[105]   If nonzero, force pattern to move down at >= 2c/3.
   var[106]   Only count changes in gens >= var[106].
   var[107]   If nonzero, don't allow cell to be changed in 2 gens between
   & var[108]     var[107] and var[108] gens apart.
   var[109]   If nonzero, force pattern to move down at >= c/2.
   var[110]   If nonzero, allow change count = 2 in at most var[110]
				  consecutive gens.
   var[111]   All gen g changes in rows between var[132] and var[134] and
              columns between var[133] and var[135] must also be present in
			  gen g+var[111].
   var[112]   If nonzero, changed region may be split vertically into
			  2 pieces of height maxheight.
   var[113]   If nonzero, changed region may be split horizontally into
			  2 pieces of width maxwidth.
   var[114]   If nonzero, HORSYMM only applies to columns >= var[114].
   var[115]   If nonzero, VERTSYMM only applies to rows >= var[115].
   var[116]   If nonzero, don't allow change counts for gens 0-18 to be
			  from 5c/9: 1344465331344465331.
   var[117]   If nonzero, bound on # of cells that are changed this gen
   & var[118]    but not changed var[118] gens ago.
   var[119]   If nonzero, don't recognize known fizzlers except at gen 0.
				 (Use this when not starting from known signals.)
   var[120]   If nonzero, lower bound on change count.
   var[121]   If nonzero, exclude change count = var[122] in gen var[121].
   & var[122]
   var[123]   If nonzero, require change count = var[124] up to but not
   & var[124]    including gen var[123].
   var[125]   If nonzero (>=2), only require height <= maxheight at least
				 once in every var[125] consecutive gens.
   var[126]   If nonzero (>=2), only require width <= maxwidth at least
				 once in every var[126] consecutive gens.
   var[127]   If nonzero, read command sets DONTCOUNT for all "01.o" cells.
				 (Use to find variants of known rotors.)
   var[128]   Change count and bounding box size must be same in generations
				 gen and gen-var[128].
   var[129]   If nonzero, bound on age sums. (cf. 137)
   var[130]   If nonzero, use Day/Night symmetry for HORSYMM, VERTSYMM,
				 and ROT180SYMM.
   var[131]   If nonzero, don't display patterns with periods <= 6.
   var[132]   Min row for var[111]
   var[133]   Min col for var[111]
   var[134]   Max row for var[111]
   var[135]   Max col for var[111]
   var[136]   Disallow birth on 2 adjacent neighbors (for Just Friends)
   var[137]   If nonzero, bound on sum, over all changed cells, of number
                 of consecutive gens in which cell was changed. (cf. 129)
   var[138]   If var[137]>0, max number of gens in which cell can be
                 unchanged but still considered "consecutive".  E.g. if
				 var[138]=2 and cell is changed in gens 1,5,8,11, then
				 it counts as 7 consecutive gens, 5 to 11.
   var[139]   If nonzero, bound on number of neighbors of changed cells.
*/

char knownrotorsandnames[MAXFILESIZE];
char *knownrotor[MAXKNOWN], *name[MAXKNOWN];
    /* Rotor descriptors and names of oscillators read from "knownrotors" */

/************************************************************************/
/*  For portability, I'm defining rand and srand as in Kernighan &      */
/*  Ritchie's book "The c Programming Language", instead of using       */
/*  those in 'standard' libraries.                                      */
/************************************************************************/

unsigned long int next=1;

int KRrand(void)    /* fcn */
{   next=next*1103515245 + 12345;
    return (unsigned int)(next/65536) % 32768;
}

void KRsrand(unsigned int seed) /* fcn */
{   next=seed;
}

int myrandom(int n) /* fcn */
/* Return random integer in range [0,n-1] */
{   return KRrand()%n;
}

/************************************************************************/
/* The following functions deal primarily with making the background    */
/* stable.  However, they also change curr[r][c] whenever they change   */
/* bkgd[r][c].                                                          */
/************************************************************************/

#define set(r0,c0) \
  inc = v - bkgd[r0][c0]; \
  bkgd[nwstng->row = (r0)][nwstng->col = (c0)] = curr[r0][c0] = \
											   nwstng->val = v; \
  naysum[(r0)-1][(c0)-1] += inc; \
  naysum[(r0)-1][(c0)  ] += inc; \
  naysum[(r0)-1][(c0)+1] += inc; \
  naysum[(r0)  ][(c0)-1] += inc; \
  naysum[(r0)  ][(c0)+1] += inc; \
  naysum[(r0)+1][(c0)-1] += inc; \
  naysum[(r0)+1][(c0)  ] += inc; \
  naysum[(r0)+1][(c0)+1] += inc;

void setbkgd(int r, int c, unsigned char v, boolean f)  /* fcn */
/* Set bkgd[r,c] to v and store choice at nwstng.
*/
{ static int inc;
#if COUNT
  countsetbkgd++;
#endif
  if (bkgd[r][c] != UNK || r==0 || r==HT-1 || c==0 || c==WD-1)
    { printf("bkgd[%d,%d] = %d\n",r,c,bkgd[r][c]); 
      err("setbkgd error");
    }

  set(r,c)
  (nwstng++)->free = f;

  /* Note that if r=c and SYMM=DIAGSYMM, a duplicate setting will be   */
  /* added to settinglist.  This won't cause any trouble; it will just */
  /* cause backup to set the cell to UNK twice.  The same will happen  */
  /* if r=HT-1-r and SYMM=HORSYMM, etc.                                */
  switch(SYMM)
    { case NOSYMM:
        break;

      case HORSYMM:
		if (var[130])  v=1-v;
		if (c >= var[114])		/* If var[114]==0, this applies to all c */
		  { set(HT-1-r,c)			(nwstng++)->free = 0; }
        break;

      case VERTSYMM:
		if (var[130])  v=1-v;
		if (r >= var[115])
		  { set(r,WD-1-c)			(nwstng++)->free = 0; }
        break;
        
      case DIAGSYMM:
		set(c,r)				(nwstng++)->free = 0;
        break;

      case ROT90SYMM:      	/* Must have HT = WD */
		set(HT-1-c,r)			(nwstng++)->free = 0;
		set(HT-1-r,HT-1-c)		(nwstng++)->free = 0;
		set(c,HT-1-r)			(nwstng++)->free = 0;
        break;

      case ROT180SYMM:
		if (var[130])  v=1-v;
		set(HT-1-r,WD-1-c)		(nwstng++)->free = 0;
        break;

      case PLUSSYMM:
		set(HT-1-r,c)			(nwstng++)->free = 0;
		set(HT-1-r,WD-1-c)		(nwstng++)->free = 0;
		set(r,WD-1-c)			(nwstng++)->free = 0;
        break;

      case XSYMM:           /* Must have HT = WD */
		set(c,r)				(nwstng++)->free = 0;
		set(HT-1-r,HT-1-c)		(nwstng++)->free = 0;
		set(HT-1-c,HT-1-r)		(nwstng++)->free = 0;
        break;

      case FULLSYMM:        /* Must have HT = WD */
		set(r,HT-1-c)			(nwstng++)->free = 0;
		set(HT-1-r,c)			(nwstng++)->free = 0;
		set(HT-1-r,HT-1-c)		(nwstng++)->free = 0;
		set(c,r)				(nwstng++)->free = 0;
		set(c,HT-1-r)			(nwstng++)->free = 0;
		set(HT-1-c,r)			(nwstng++)->free = 0;
		set(HT-1-c,HT-1-r)		(nwstng++)->free = 0;
        break;
        
      default:
        err("Bad symmetry type");
      }

  return;
}

/* -------------------------------------------------------------------- */

char consistify(int r, int c)   /* fcn */
/* Examine neighborhood of (r,c).  If inconsistent, return ERR.  If value
   of bkgd[r,c] or any of its neighbors is forced, set it and add the
   setting at nwstng.
   Assumes that 0<r<HT-1 and 0<c<WD-1.
*/
{ static int tablevalue;

#if COUNT
  countconsistify++;
#endif

START:

  if ((tablevalue = consistable[bkgd[r][c]][naysum[r][c]]) == NOINFO)
	return OK;
        /* This is the most common case, so do it separately for speed */

  switch(tablevalue)
    { case INCONSIS:
        return ERR;

      case NAYSOFF:
        if (bkgd[r-1][c-1] == UNK)  setbkgd(r-1,c-1,OFF,0);
        if (bkgd[r-1][c  ] == UNK)  setbkgd(r-1,c  ,OFF,0);
        if (bkgd[r-1][c+1] == UNK)  setbkgd(r-1,c+1,OFF,0);
        if (bkgd[r  ][c-1] == UNK)  setbkgd(r  ,c-1,OFF,0);
        if (bkgd[r  ][c+1] == UNK)  setbkgd(r  ,c+1,OFF,0);
        if (bkgd[r+1][c-1] == UNK)  setbkgd(r+1,c-1,OFF,0);
        if (bkgd[r+1][c  ] == UNK)  setbkgd(r+1,c  ,OFF,0);
        if (bkgd[r+1][c+1] == UNK)  setbkgd(r+1,c+1,OFF,0);
        return OK;

      case NAYSON:
        if (bkgd[r-1][c-1] == UNK)  setbkgd(r-1,c-1,ON,0);
        if (bkgd[r-1][c  ] == UNK)  setbkgd(r-1,c  ,ON,0);
        if (bkgd[r-1][c+1] == UNK)  setbkgd(r-1,c+1,ON,0);
        if (bkgd[r  ][c-1] == UNK)  setbkgd(r  ,c-1,ON,0);
        if (bkgd[r  ][c+1] == UNK)  setbkgd(r  ,c+1,ON,0);
        if (bkgd[r+1][c-1] == UNK)  setbkgd(r+1,c-1,ON,0);
        if (bkgd[r+1][c  ] == UNK)  setbkgd(r+1,c  ,ON,0);
        if (bkgd[r+1][c+1] == UNK)  setbkgd(r+1,c+1,ON,0);
        return OK;

      case CELLOFF:
        setbkgd(r,c,OFF,0);
        goto START;

      case CELLON:
        setbkgd(r,c,ON,0);
        goto START;

      case NOINFO:
        return OK;
    }
}

/* -------------------------------------------------------------------- */

char consis9(int r, int c)  /* fcn */
/* Call consistify for (r,c) and each of its 8 neighbors.  If inconsistency
   found, return ERR.
*/
{ 

#if COUNT
  countconsis9++;
#endif

  if (consistify(r,c))  return ERR;
  if (consistify(r-1,c  ))  return ERR;
  if (consistify(r  ,c-1))  return ERR;
  if (consistify(r+1,c  ))  return ERR;
  if (consistify(r  ,c+1))  return ERR;
  if (consistify(r-1,c-1))  return ERR;
  if (consistify(r+1,c-1))  return ERR;
  if (consistify(r-1,c+1))  return ERR;
  if (consistify(r+1,c+1))  return ERR;
  return OK;
}

/* -------------------------------------------------------------------- */

char proceed(int r, int c, unsigned char v, boolean f)  /* fcn */
/* Set bkgd[r,c] to v and examine consequences.  Return ERR if
   inconsistency found.
*/
{

#if COUNT
  countproceed++;
  if (TRACE) {printf("  proceed(%d,%d,%d,%d)\n",r,c,(int)v,(int)f);fflush(stdout);}
#endif

  nxstng = nwstng;
  setbkgd(r,c,v,f);
  while (nxstng != nwstng)
    if (consis9(nxstng->row, nxstng->col))  return ERR;
    else  nxstng++;
  return OK;
}

/* -------------------------------------------------------------------- */


char backup(void)   /* fcn */
/* Back up to last free choice.  Return ERR if none left.
   After return, nwstng still points to previous choice.
*/
{ static int r, c, inc;

#if COUNT
  countbackup++;
  if (TRACE)  {printf("  backup()\n");fflush(stdout);}
#endif

  while (nwstng > settinglist)
    { nwstng--;
      inc = UNK - bkgd[r = nwstng->row][c = nwstng->col];
      bkgd[r][c] = UNK;
      curr[r][c] = UNK;
	  naysum[r-1][c-1] += inc;
	  naysum[r-1][c  ] += inc;
	  naysum[r-1][c+1] += inc;
	  naysum[r  ][c-1] += inc;
	  naysum[r  ][c+1] += inc;
	  naysum[r+1][c-1] += inc;
	  naysum[r+1][c  ] += inc;
	  naysum[r+1][c+1] += inc;
      if (nwstng->free)  return OK;
    }

  return ERR;
}

/* -------------------------------------------------------------------- */

char go(int r, int c, unsigned char v, boolean f,   /* fcn */
        setting **wasfree)  /* fcn */
/* Try to set bkgd[r,c] to v, backing up if necessary.  Return ERR if
   try to back up beyond start of setting list.  wasfree becomes ptr
   into settinglist, to last choice that was originally free.
*/
{ 

#if COUNT
  countgo++;
  if (TRACE)  {printf(" go(%d,%d,%d,%d,-)\n",r,c,(int)v,(int)f);fflush(stdout);}
#endif

  *wasfree = nwstng;
  while (proceed(r,c,v,f))
    { if (backup())  return ERR;
      r = nwstng->row;
      c = nwstng->col;
      v = !nwstng->val;
      f = 0;
      *wasfree = nwstng;
    }
  return OK;
}

/************************************************************************/
/* The following functions handle advancing the pattern.                */
/************************************************************************/

#define append(r,c) \
  for (nay=nays[gen], alreadyinlist=0; nay<chgd[gen+1]; nay++) \
    if (nay->row == r && nay->col == c)  { alreadyinlist=1; break; } \
  if (!alreadyinlist) \
    { chgd[gen+1]->row = r;  (chgd[gen+1]++)->col = c; \
    }
#if 0
      chgd[gen+1]->row = -1;}
#endif
/*  zzz      Is this needed?   */



void listneighbors(int gen) /* fcn */
/* Given list of changed cells in generation gen (from chgd[gen] to
   nays[gen]-1), create list of their neighbors (from nays[gen] to
   chgd[gen+1]).
*/
{ static point *chg, *nay, tmp, *p, *q;
  static boolean alreadyinlist;
  static int r,c,tmpcount;

#if COUNT
  countlistneighbors++;
  if (TRACE) {printf("listneighbors(%d)\n",gen);fflush(stdout);}
#endif

  chgd[gen+1] = nays[gen];
  for (chg=chgd[gen]; chg<nays[gen]; chg++)
    { r = chg->row;  c = chg->col;
      append(r,c);
      append(r-1,c-1);  append(r-1,c  );  append(r-1,c+1);
      append(r  ,c-1);                    append(r  ,c+1);
      append(r+1,c-1);  append(r+1,c  );  append(r+1,c+1);
    }

  if (chgd[gen+1]-chglist>=CHGLISTLTH)				/* From Gabriel Nivasch */
     { printf("chgd[%d]-chglist = %d\n", gen, chgd[gen+1]-chglist);
	   err("Overflow of chglist at listneighbors().\n")	/* From GN */
	 }

  /* Bubble sort into increasing order by number of UNK neighbors */
  for (p=nays[gen]+1; p<chgd[gen+1]; p++)
    { tmp = *p;
      r = tmp.row;  c = tmp.col;
      tmpcount = bkgd[r-1][c-1] + bkgd[r-1][c] + bkgd[r-1][c+1] +
                 bkgd[r  ][c-1] + bkgd[r  ][c] + bkgd[r  ][c+1] +
                 bkgd[r+1][c-1] + bkgd[r+1][c] + bkgd[r+1][c+1];
      for (q=p-1; q>=nays[gen]; q--)
        { r = q->row;  c = q->col;
          if (bkgd[r-1][c-1] + bkgd[r-1][c] + bkgd[r-1][c+1] +
              bkgd[r  ][c-1] + bkgd[r  ][c] + bkgd[r  ][c+1] +
              bkgd[r+1][c-1] + bkgd[r+1][c] + bkgd[r+1][c+1] > tmpcount)
            *(q+1) = *q;
          else  break;
        }
      *(q+1) = tmp;
    }
}

/* -------------------------------------------------------------------- */

void changecurr(unsigned char curr[][MAXWD], int gen)   /* fcn */
/* Toggle values in curr of cells pointed to by chgd[gen].  This is used
   to either make curr equal to bkgd or to make it contain the current
   generation.
*/
{ point *p;
  int r,c;

#if COUNT
  countchangecurr++;
#endif

  for (p=chgd[gen]; p<nays[gen]; p++)
    { r = p->row;  c = p->col;
      if (curr[r][c] == UNK)  err2("BUG in changecurr: r=%d c=%d\n",r,c);
      curr[r][c] ^= 1;
    }
}

/* -------------------------------------------------------------------- */

char nxgen(int r, int c)    /* fcn */
/* Given curr[r-1,c-1], ..., curr[r+1,c+1].  Tries to compute next gen
   of cell (r,c).  Returns UNK if can't tell or if curr[r,c] is UNK.
*/
{ unsigned char val, newval;

#if COUNT
  countnxgen++;
#endif

  if ((val = curr[r][c]) == UNK)  return UNK;
  newval = transtable[val]
    [ curr[r-1][c-1] + curr[r-1][c] + curr[r-1][c+1] +
      curr[r  ][c-1] +                curr[r  ][c+1] +
      curr[r+1][c-1] + curr[r+1][c] + curr[r+1][c+1] ];
  if (!var[136])  return newval;

  /* Kludge for Just Friends rule */
  /* If var[136], then disallow birth if the parents are adjacent */
  if (val!=OFF || newval!=ON)  return newval;
  if ((curr[r-1][c-1] + curr[r-1][c  ] == 2) ||
      (curr[r-1][c  ] + curr[r-1][c+1] == 2) ||
      (curr[r-1][c+1] + curr[r  ][c+1] == 2) ||
      (curr[r  ][c+1] + curr[r+1][c+1] == 2) ||
      (curr[r+1][c+1] + curr[r+1][c  ] == 2) ||
      (curr[r+1][c  ] + curr[r+1][c-1] == 2) ||
      (curr[r+1][c-1] + curr[r  ][c-1] == 2) ||
      (curr[r  ][c-1] + curr[r-1][c-1] == 2))  return OFF;
  return ON;
}

/* -------------------------------------------------------------------- */

void findchgcount(int g)	/* fcn */
/* Compute and store in chgcount[g] the number of changed cells in gen g,
   not counting any for which DONTCOUNT is TRUE.  This should only be
   called when gen g has been finished; i.e. nays[g] has been set.
   Also compute and save width and height of changed region.
*/
{ static point *q, *qm;
  static int minr,maxr,minc,maxc, r,c, qr,qc,g0;

  if (g<var[106])  chgcount[g] = width[g] = height[g] = 0;
  else
    { chgcount[g]=0;
	  maxc = maxr = -1;
	  minr = HT;  minc = WD;
	  for (q=chgd[g]; q<nays[g]; q++)
        if (!(flag[r = q->row][c = q->col] & DONTCOUNT))
	      { chgcount[g]++;
			if (r<minr)  minr=r;
			if (r>maxr)  maxr=r;
			if (c<minc)  minc=c;
			if (c>maxc)  maxc=c;
		  }

	  width[g]  = maxc==-1 ? 0 : maxc - minc + 1;
	  height[g] = maxr==-1 ? 0 : maxr - minr + 1;
	}

  /* Record sum of ages of all changes */
  if (var[129] && g==0)
    { for (q=chgd[g], agesm[g]=0; q<nays[g]; q++)
        { qr = q->row;  qc = q->col;
          for (g0=0; g0<g; g0++)
            for (qm=chgd[g0]; qm<nays[g0]; qm++)
              if (qm->row == qr && qm->col == qc)
                goto FOUNDFIRST;

        FOUNDFIRST:
          agesm[g] += g-g0+1;
        }
    }
}

/* -------------------------------------------------------------------- */

char trycell(void)  /* fcn */
/* Given gen>0, nays[gen-1] <= nay < chgd[gen] <= chg.
   Tries to compute next gen of cell specified by nay, possibly setting
   bkgd of its neighbors to ON or OFF.  Returns ERR if problem occurs, in
   which case we must back up.
*/
{ static int r,c,g,qr,qc,minqr,maxqr,minqc,maxqc,ru,cu,changecount,agesum,pgen;
  static unsigned char val;
  static point *q, *qm;
  static boolean oldchange, recentbig,recentwide,recenttall,changeding;

  r = nay->row;  c = nay->col;

#if COUNT
  counttrycell++;
  if (TRACE)  {printf(" trycell()\n");fflush(stdout);}
#endif

  while ((val = nxgen(r,c)) == UNK)
    { if (bkgd[r][c] == UNK)            { ru = r;    cu = c;   }
      else if (bkgd[r-1][c  ] == UNK)   { ru = r-1;  cu = c;   }
      else if (bkgd[r  ][c-1] == UNK)   { ru = r;    cu = c-1; }
      else if (bkgd[r+1][c  ] == UNK)   { ru = r+1;  cu = c;   }
      else if (bkgd[r  ][c+1] == UNK)   { ru = r;    cu = c+1; }
      else if (bkgd[r-1][c-1] == UNK)   { ru = r-1;  cu = c-1; }
      else if (bkgd[r+1][c-1] == UNK)   { ru = r+1;  cu = c-1; }
      else if (bkgd[r-1][c+1] == UNK)   { ru = r-1;  cu = c+1; }
      else if (bkgd[r+1][c+1] == UNK)   { ru = r+1;  cu = c+1; }

      nwstng->gen = gen;
      nwstng->nay = nay;
      nwstng->chg = chg;

      val = myrandom(100)<prob;
      if (proceed(ru,cu,val,1))  return ERR;
    }
  
  if (val == bkgd[r][c])  return OK;

  /* Now we know that the cell is changed in the next generation, */
  /* so we perform some tests to see if that's permitted.         */

  chg->row = r;  (chg++)->col = c;
  if (chg-chglist>=CHGLISTLTH)					/* From Gabriel Nivasch */
    err("Overflow of chglist at trycell().\n")	/* From Gabriel Nivasch */


  /* Test for failure, based on active region getting too big. */

  if (gen < var[106])  return OK;

  if (flag[r][c] & DONTCHANGE)  return ERR;

  if (flag[r][c] & DONTCOUNT)  return OK;

  recentbig = recentwide = recenttall = TRUE;
  if (var[100])
	for (g=gen-1; g>=0 && g>gen-var[100]; g--)
	  if (chgcount[g] <= maxchng)		{ recentbig = FALSE; break; }
  /* recentbig is TRUE if chgcount too large for last var[100]-1 gens */

  if (var[125])
	for (g=gen-1; g>=0 && g>gen-var[125]; g--)
	  if (height[g] <= maxheight)			{ recenttall = FALSE; break; }
  /* recenttall is TRUE if height too large for last var[125]-1 gens */

  if (var[126])
	for (g=gen-1; g>=0 && g>gen-var[126]; g--)
	  if (width[g] <= maxwidth)			{ recentwide = FALSE; break; }
  /* recentwide is TRUE if width too large for last var[126]-1 gens */

  /* Find min & max row & column of changed cells. */
  for (q=chgd[gen], minqr=HT, minqc=WD, maxqr=maxqc=0; q<chg; q++)
    { qr = q->row;  qc = q->col;
      if (flag[qr][qc] & DONTCOUNT)  continue;
	  if (qr < minqr)  minqr=qr;
	  if (qr > maxqr)  maxqr=qr;
	  if (qc < minqc)  minqc=qc;
	  if (qc > maxqc)  maxqc=qc;
	}

  /* Make sure changed region isn't too tall */
  if (var[112])
	{ for (q=chgd[gen]; q<chg-1; q++)
		if ((qr = q->row) >= minqr+maxheight && qr <= maxqr-maxheight)
		  return ERR;
	}
  else
	if (maxqr-minqr >= maxheight && recenttall && gen>=numcc)  return ERR;
	
  /* Make sure changed region isn't too wide */
  if (var[113])
	{ for (q=chgd[gen]; q<chg-1; q++)
		if ((qc = q->col) >= minqc+maxwidth && qc <= maxqc-maxwidth)
		  return ERR;
	}
  else
	if (maxqc-minqc >= maxwidth && recentwide && gen>=numcc)  return ERR;

  for (q=chgd[gen], changecount=1; q<chg-1; q++)
    { qr = q->row;  qc = q->col;
      if (flag[qr][qc] & DONTCOUNT)  continue;
      if (++changecount > maxchng && recentbig && gen>=numcc)  return ERR;

      /* Put bounds on bounding diamond */
      if ((var[102] && diff(r-c,qr-qc) >= var[102]) ||
		  (var[103] && diff(r+c,qr+qc) >= var[103]))
        return ERR;
    }

  if (var[107])			/* Don't allow cell to be changed in 2 gens */
						/* between var[107] and var[108] gens apart */
	for (g=gen-var[107]; g>=gen-var[108] && g>=0; g--)
	  for (q=chgd[g]; q<nays[g]; q++)
	    if (q->row == r && q->col == c)  return ERR;

/* Put bound on number of changed cells that weren't changed last gen */
  if (var[101] && gen > 0)
	for (q=chgd[gen], changecount=0; q<chg; q++)
	  { qr = q->row;  qc = q->col;
		for (qm=chgd[gen-1], oldchange=FALSE; qm<nays[gen-1]; qm++)
		  if (qm->row == qr && qm->col == qc)
			{ oldchange = TRUE;		/* Change isn't new, so don't count it */
			  break;
			}
		if (!oldchange)  changecount++;
		if (changecount > var[101])  return ERR;
	  }

/* Put bound on # of changed cells that weren't changed var[118] gens ago */
  if (var[117] && gen >= var[118])
	for (q=chgd[gen], changecount=0; q<chg; q++)
	  { qr = q->row;  qc = q->col;
		for (qm=chgd[gen-var[118]], oldchange=FALSE;
			 qm<nays[gen-var[118]];
			 qm++)
		  if (qm->row == qr && qm->col == qc)
			{ oldchange = TRUE;		/* Change isn't new, so don't count it */
			  break;
			}
		if (!oldchange)  changecount++;
		if (changecount > var[117])  return ERR;
	  }

/* Bound sum of ages of all changes */
  if (var[129])
    { for (g=0; g<gen; g++)
        for (qm=chgd[g]; qm<nays[g]; qm++)
          if (qm->row == r && qm->col == c)
            goto FOUNDFIRST;

    FOUNDFIRST:
      agesm[gen] += gen-g+1;
      if (agesm[gen] > var[129])  return ERR;
    }

#if 0
/* Bound sum of # of consec gens in which cell has been changed.  6/27/2008 */
  if (var[137])
    { for (q=chgd[gen], agesum=0; q<chg; q++)
	    { qr = q->row;  qc = q->col;
/* printf("gen=%d  qr=%d  qc=%d  agesum=%d  var[137]=%d\n",gen,qr,qc,agesum,var[137]); */
		  if (++agesum > var[137])    /* Count change in current gen */
		    return ERR;
		  for (g=gen-1; g>=0; g--)
		    { 
/*printf("gen=%d g=%d agesum=%d\n",gen,g,agesum); */
			  for (qm=chgd[g]; qm<nays[g]; qm++)
			    if (qm->row == qr && qm->col == qc)
				  goto CHANGEDING;
			  break;

			CHANGEDING:
/*printf("CHANGEDING:  gen=%d  g=%d  qr=%d  qc=%d  agesum was %d\n",gen,g,qr,qc,agesum); */
			  if (++agesum > var[137])  return ERR;
			}
		}
	}
#endif

/* Bound sum of # of consec gens in which cell has been changed.  6/27/2008 */
  if (var[137])
    { for (q=chgd[gen], agesum=0; q<chg; q++)
	    { qr = q->row;  qc = q->col;
		  pgen = gen;  /* Earliest generation in "consecutive" */
		               /* group in which cell changed          */
		  for (g=gen-1; g>=0 && g>=pgen-var[138]-1; g--)
		    for (qm=chgd[g]; qm<nays[g]; qm++)
			  if (qm->row == qr && qm->col == qc)
				{ pgen = g;
				  break;
				}

		  agesum += gen-pgen+1;
		  if (agesum > var[137])  return ERR;
		}
	}

/* Force signal to move right at >= 2c/3 */
  if (var[104] && 3*c < 2*gen + var[104])  return ERR;

/* Force signal to move down at >= 2c/3 */
  if (var[105] && 3*r < 2*gen + var[105])  return ERR;

/* Force signal to move down at >= c/2 */
  if (var[109] && 2*r < gen + var[109])  return ERR;

  return OK;
}

/* -------------------------------------------------------------------- */

char computecellorbackup(void)  /* fcn */
/* Given gen>0, nays[gen-1] <= nay < chgd[gen] <= chg.
   Tries to compute next generation of cell specified by nay, possibly
   setting bkgd of its neighbors to ON or OFF.  If problem occurs, backs up,
   possibly decreasing gen, nay, and chg.  If can't back up, returns ERR,
   in which case no more objects exist.
*/
{ static setting *wasfree;
  static int g,qr,qc;
  static point *q, *qm;

  if (++countcomporbackuplo == 1000000)
    { countcomporbackuplo = 0;
	  countcomporbackuphi++;
	}

#if COUNT
  if (TRACE)  {printf("computecellorbackup()\n");fflush(stdout);}
#endif

  if (!found &&
      !(var[139] && (chgd[gen]-nays[gen-1])>var[139]) &&
	  trycell()==OK)
    { nay++;
      return OK;
    }

  changecurr(curr, gen-1);  /* curr <- bkgd */
  found = FALSE;

  if (backup())  return ERR;

  if (go(nwstng->row, nwstng->col, !nwstng->val, 0, &wasfree))
    return ERR;
  
  gen = wasfree->gen;
  nay = wasfree->nay;
  chg = wasfree->chg;
  changecurr(curr, gen-1);

  /* Recompute value of agesm[gen] */
  if (var[129])
    { for (q=chgd[gen], agesm[gen]=0; q<chg; q++)
        { qr = q->row;  qc = q->col;
          for (g=0; g<gen; g++)
            for (qm=chgd[g]; qm<nays[g]; qm++)
              if (qm->row == qr && qm->col == qc)
                goto FOUNDFIRST;

        FOUNDFIRST:
          agesm[gen] += gen-g+1;
        }
    }

  return OK;
}

/************************************************************************/
/* The following functions analyze finished patterns and print          */
/* information about them.                                              */
/************************************************************************/

int period(void)    /* fcn */
/* Checks to see if changes in generation gen are the same as in some
   previous generation.  If so, returns period.  Otherwise, returns 0.
*/
{ static int g, numchgs;
  static point *q;

  numchgs = nays[gen] - chgd[gen];

  for (g=gen-1; g>=0; g--)
    { if (nays[g] - chgd[g] != numchgs)  continue;
      for (q=chgd[g]; q<nays[g]; q++)
        if (curr[q->row][q->col] == bkgd[q->row][q->col])  break;
      if (q==nays[g])  return gen-g;
    }
  
  return 0;
}

/* -------------------------------------------------------------------- */

void dispchgcts(int g) /* fcn */
{ int i;
  boolean different;

  different=FALSE;
  for (i=0; i<=g; i++)
    if (nays[i]-chgd[i] != chgcount[i]) different = TRUE;

  if (different)
    { printf("Full change counts:");
      for (i=0; i<=g; i++)
	    { if (i && i%5 == 0)  printf(" ");
	      printf(" %d",nays[i]-chgd[i]);
	    }
      printf("\n");
	}

  printf("Change counts:");
  for (i=0; i<=g; i++)
	{ if (i && i%5 == 0)  printf(" ");
	  printf(" %d",chgcount[i]);
	}
  printf("\n");

  printf("Sizes:");
  for (i=0; i<=g; i++)
	{ if (i && i%5 == 0)  printf(" ");
	  printf(" %dx%d",height[i],width[i]);
	}
  printf("\n");

  if (var[129])
    { printf("Age sums:");
      for (i=0; i<=g; i++)
	    { if (i && i%5 == 0)  printf(" ");
	      printf(" %d",agesm[i]);
	    }
      printf("\n");
	}
}

/* -------------------------------------------------------------------- */

void display(int g) /* fcn */
/* Print non-UNK part of bkgd and gen g. */
{ static int r,c,minr,maxr,minc,maxc,lastc;
  static boolean changed;
  static point *p;

  if (NOPICS)
    { printf("\n");
	  return;
	}

  minr = HT;    maxr = -1;
  minc = WD;    maxc = -1;

  /* Find bounding box of union of ON part of bkgd */
  /* and set of changed cells in gen g. */
  for (r=2; r<HT-2; r++)
    for (c=2; c<WD-2; c++)
      if (bkgd[r][c] == ON)
        { if (r<minr)  minr = r;
          if (r>maxr)  maxr = r;
          if (c<minc)  minc = c;
          if (c>maxc)  maxc = c;
        }
  for (p=chgd[g]; p<nays[g]; p++)
    { if (p->row < minr)  minr = p->row;
      if (p->row > maxr)  maxr = p->row;
      if (p->col < minc)  minc = p->col;
      if (p->col > maxc)  maxc = p->col;
    }

  if (maxr<0)  { minr=maxr=HT/2;  minc=maxc=WD/2; }
  minr-=2;  maxr+=2;  minc-=2;  maxc+=2;
  printf("Gen %d.  Rows %d - %d.  Cols %d - %d.\n",g,minr,maxr,minc,maxc);

  for (r=minr; r<=maxr; r++)
    { for (lastc=maxc; lastc>minc; lastc--)
        if (bkgd[r][lastc] != UNK)  break;

      for (c=minc; c<=lastc; c++)
        { if (bkgd[r][c] != UNK)
            for (p=chgd[g], changed=FALSE; p<nays[g]; p++)
              if (p->row == r && p->col == c)
                { changed = TRUE;
                  break;
                }

          printf("%c", bkgd[r][c]==OFF ? (changed ? '1' : '.') :
                       bkgd[r][c]==ON  ? (changed ? '0' : 'o') :
                       ',');
        }
      printf("\n");
    }
  
  fflush(stdout);
}

/* -------------------------------------------------------------------- */

enum {NWR, NER, SWR, SER, NWC, NEC, SWC, SEC};

void getrotordesc(unsigned char cell[][MAXWD], int period,  /* fcn */
  int minr, int maxr, int minc, int maxc, int orientation,  /* fcn */
  char *string) /* fcn */
/* Given description of bounding box of rotor in one generation, in
   cell[minr:maxr][minc:maxc]  (Each element is OFF, ON, or STATOR).
   Orientation has 1 of 8 values:  NWR, NWC, NER, NEC, SWR, SWC, SER, SEC.
   For example, NWR means to list cells by rows starting at the NW corner.
   Writes rotor description into string.  Description has 1 char per cell,
   giving # of live stator neighbors (found from curr[][]) and current
   state of cell: '0' to '8' represent dead cells, '@' to 'H' represent
   live cells.  '.' represents cell not in rotor.
*/
{ int h, w, maxcount, r0, c0, dr0, dc0, dr1, dc1, lth, count, r, c, rotorsize;
  char *p;

  h = maxr-minr+1;
  w = maxc-minc+1;

  switch(orientation)
    { case NWR: r0=minr; c0=minc; dr0= 0; dc0= 1; dr1= 1; dc1= 0; lth=w; break;
      case NWC: r0=minr; c0=minc; dr0= 1; dc0= 0; dr1= 0; dc1= 1; lth=h; break;
      case NER: r0=minr; c0=maxc; dr0= 0; dc0=-1; dr1= 1; dc1= 0; lth=w; break;
      case NEC: r0=minr; c0=maxc; dr0= 1; dc0= 0; dr1= 0; dc1=-1; lth=h; break;
      case SWR: r0=maxr; c0=minc; dr0= 0; dc0= 1; dr1=-1; dc1= 0; lth=w; break;
      case SWC: r0=maxr; c0=minc; dr0=-1; dc0= 0; dr1= 0; dc1= 1; lth=h; break;
      case SER: r0=maxr; c0=maxc; dr0= 0; dc0=-1; dr1=-1; dc1= 0; lth=w; break;
      case SEC: r0=maxr; c0=maxc; dr0=-1; dc0= 0; dr1= 0; dc1=-1; lth=h; break;
    }

  maxcount = h*w + h + w - lth;

  rotorsize = 0;
  for (r=minr; r<=maxr; r++)
    for (c=minc; c<=maxc; c++)
      rotorsize += (cell[r][c] != STATOR);

  sprintf(string,"p%d r%d %dx%d ", period, rotorsize, h+w-lth, lth);
  for (p=string; *p; p++);

  if ((p-string)+maxcount > MAXROTORDESCLTH)
	{ sprintf(p," ROTOR DESCRIPTOR TOO LONG.");
	  return;
	}

  for (r=r0, c=c0, count=1;
       count<maxcount;
       r+=dr0, c+=dc0, count++, p++)
    if (count%(lth+1) == 0)
      { r+=dr1-(lth+1)*dr0;
        c+=dc1-(lth+1)*dc0;
        *p = ' ';
      }
    else
      { if (cell[r][c] == STATOR)  *p = '.';
        else
          { *p = cell[r][c] ? '@' : '0';
            if (cell[r-1][c-1] == STATOR)  *p += curr[r-1][c-1];
            if (cell[r-1][c  ] == STATOR)  *p += curr[r-1][c  ];
            if (cell[r-1][c+1] == STATOR)  *p += curr[r-1][c+1];
            if (cell[r  ][c-1] == STATOR)  *p += curr[r  ][c-1];
            if (cell[r  ][c+1] == STATOR)  *p += curr[r  ][c+1];
            if (cell[r+1][c-1] == STATOR)  *p += curr[r+1][c-1];
            if (cell[r+1][c  ] == STATOR)  *p += curr[r+1][c  ];
            if (cell[r+1][c+1] == STATOR)  *p += curr[r+1][c+1];
          }
      }

  *p = 0;
}

/* -------------------------------------------------------------------- */

unsigned long int newhash(unsigned long int hash,   /* fcn */
                          unsigned long int update) /* fcn */
{ return 31*hash+update;
}

/* -------------------------------------------------------------------- */

unsigned long int hash(void)    /* fcn */
/* Computes hash function based on list of changes in gens 0 to gen */
{ static int g;
  static point *q;
  static unsigned long int h;

  for (g=h=0; g<=gen; g++)
    for (q=chgd[g]; q<nays[g]; q++)
      h = newhash(newhash(h,q->row), q->col);
  return h;
}

/* -------------------------------------------------------------------- */

unsigned long int hashtable[HASHTBLSIZE];
int hashcount;

char hashnew(unsigned long int h)   /* fcn */
/* Look up h in hashtable.  If found, return FALSE.  If not, add to
   table and return TRUE.
*/
{ int i;

#if COUNT
if (TRACE) {printf("hashnew(%lx).  hashcount = %d\n", h,hashcount);fflush(stdout);}
#endif

  for (i=0; i<hashcount; i++)
    if (h == hashtable[i])
	  return FALSE;

  hashtable[hashcount++] = h;
  if (hashcount>=HASHTBLSIZE)  err("Hash table overflow");
  return TRUE;
}

/* -------------------------------------------------------------------- */

void fillcell(unsigned char cell[][MAXWD], int mingen, int maxgen,	/* fcn */
			  int *minr, int *maxr, int *minc, int *maxc)			/* fcn */
/* Make cell describe the rotor from generation mingen to generation
   maxgen.  I.e. a cell that has the same value in all of those generations
   gets set to STATOR; other cells get set to their values in generation
   mingen.  Set minr, maxr, minc, and maxc to the min and max row and column
   of the rotor.
*/
{ static point *q0, *q1;
  static int g0, g1, count, r, c;

  /* Fill cell with STATOR cells */
  for (r=0; r<WD; r++)  for (c=0; c<HT; c++)  cell[r][c] = STATOR;

  *minr = HT;  *maxr = -1;
  *minc = WD;  *maxc = -1;

  /* For each cell that's changed in some gen, count how many gens it's */
  /* changed in.  If it's changed in all gens, it's not in rotor.       */
  for (g0=mingen; g0<=maxgen; g0++)
    for (q0=chgd[g0]; q0<nays[g0]; q0++)
      { r = q0->row;  c = q0->col;
        count = 0;
        for (g1=mingen; g1<=maxgen; g1++)
          for (q1=chgd[g1]; q1<nays[g1]; q1++)
            if (q1->row == r && q1->col == c)
              { if (g1<g0)  goto NOTFIRSTOCCURRENCE;
                count++;
              }
        if (count < maxgen-mingen+1)
          { /* Found another rotor cell */
            cell[r][c] = bkgd[r][c];
			for (q1=chgd[mingen]; q1<nays[mingen]; q1++)
			  if (q1->row == r && q1->col == c)  cell[r][c] ^= 1;

            if (r < *minr)  *minr = r;
            if (r > *maxr)  *maxr = r;
            if (c < *minc)  *minc = c;
            if (c > *maxc)  *maxc = c;
          }
      NOTFIRSTOCCURRENCE:;
      }
}

/* -------------------------------------------------------------------- */

void printoscinfo(int p, char prefix)  /* fcn */
/* If prefix = 'p', print info about pattern, known to have period p.
   If prefix is 'f', 's', or 'u', then print info about the entire Life
   history of the object, starting at gen 0.
*/
{ static point *q0, *q1;
  static int g0, g1, count, r, c, rotorsize, i, minr,maxr,minc,maxc,
    orientation,minorient,maxorient,g, dist, changesome, unconcount, rn, cn;
  static char rotordesc[MAXROTORDESCLTH], minrotordesc[MAXROTORDESCLTH];
  static unsigned char cell[MAXHT][MAXWD];
  static boolean known;

  if (prefix == 'p')		/* Describe rotor of oscillator */
    { fillcell(cell, gen-p, gen-1, &minr,&maxr,&minc,&maxc);
      if (maxc-minc > maxr-minr)		{ minorient = NWR; maxorient = SER; }
      else if (maxc-minc < maxr-minr)	{ minorient = NWC; maxorient = SEC; }
      else								{ minorient = NWR; maxorient = SEC; }

	  strcpy(minrotordesc, "z");  /* Lexicographically larger than any
                               	     rotor descriptor */
	  for (g=gen; g > gen-p; g--)
        { for (orientation = minorient; orientation <= maxorient; orientation++)
            { getrotordesc(cell, p, minr, maxr, minc, maxc, orientation,
																  rotordesc);
              if (strcmp(rotordesc, minrotordesc) < 0)
                strcpy(minrotordesc, rotordesc);
            }

          /* Done with generation g.  Change cell to generation g-1. */
          changecurr(cell, g);  changecurr(cell, g-1);
		}

      /* Check to see if minrotordesc is in list of known rotors */
      for (i=0, known=FALSE; knownrotor[i]; i++)
        if (strcmp(knownrotor[i], minrotordesc) == 0)
		  { known = TRUE;
			break;
		  }
      }
  else		/* Describe 'rotor' of fizzler */
	/* If, at some generation between 0 and gen-1, the fizzler becomes the
	   same as some known fizzler, we'll print the rotor descriptor at that
	   generation.  Otherwise, we'll print the unknown rotor descriptor at
	   gen 0.  (The loop runs to g==gen, but the last time through we
	   really do g==0 again.)
	*/
	for (g=0; g <= (var[119] ? 0 : gen); g++)
	  { fillcell(cell, (g==gen ? 0 : g), gen, &minr,&maxr,&minc,&maxc);
        if (maxc-minc > maxr-minr)		{ minorient = NWR; maxorient = SER; }
        else if (maxc-minc < maxr-minr)	{ minorient = NWC; maxorient = SEC; }
        else							{ minorient = NWR; maxorient = SEC; }

	    strcpy(minrotordesc, "z");
	    for (orientation = minorient; orientation <= maxorient; orientation++)
          { getrotordesc(cell, gen-(g==gen ? 0 : g), minr, maxr, minc, maxc,
						                              orientation, rotordesc);
		    *rotordesc = prefix;  /* Fizzle or eventually periodic */
            if (strcmp(rotordesc, minrotordesc) < 0)
              strcpy(minrotordesc, rotordesc);
          }

        /* Check to see if minrotordesc is in list of known rotors */
        for (i=0, known=FALSE; knownrotor[i]; i++)
          if (strcmp(knownrotor[i], minrotordesc) == 0)
			{ known = TRUE;
			  break;
		    }
		   
		if (known)  break;	/* Becomes known fizzler in gen g */
      }

  /* Check to see if rotor is disconnected.  If so, find distance between */
  /* components.  I'm using an inefficient algorithm here, but it doesn't */
  /* matter, because this doesn't happen very often.                      */

  /* Find a rotor cell (in top row) and change it to CONN */
  for (c=minc; c<=maxc; c++)
    if (cell[minr][c] <= ON)
      { cell[minr][c] = CONN;
        break;
      }

  /* For dist = 1, 2, ... mark all rotor cells as CONN if they're within */
  /* units of a CONN cell.  Continue until all rotor cells are marked.   */
  for (dist=1, unconcount=TRUE; dist<100 && unconcount; dist++)
    do
      { changesome = FALSE;
        unconcount = 0;
        for (r=minr; r<=maxr; r++)  for (c=minc; c<=maxc; c++)
          if (cell[r][c] <= ON)
            { for (rn=max(r-dist,minr); rn<=r+dist && rn<=maxr; rn++)
                for (cn=max(c-dist,minc); cn<=c+dist && cn<=maxc; cn++)
                  if (cell[rn][cn] == CONN)
                    { cell[r][c] = CONN;
                      changesome = TRUE;
                    }
              if (cell[r][c] <= ON)  unconcount++;
            }
      }
    while (changesome);
  dist--;

  /* Print rotor descriptor and name, if known */
  printf("%s\t", minrotordesc);
  if (known)
    if (dist>1)         printf("   %s (gap = %d)\n", name[i],dist-1);
    else                printf("   %s\n", name[i]);
  else
    if (dist>1)         printf("<- unknown: split rotor (gap = %d)\n",dist-1);
    else
	  if (prefix == 'u')    printf("<- unknown\n");
	  else                  printf("<- UNKNOWN\n");
  fflush(stdout);
}

/* -------------------------------------------------------------------- */

#if COUNT
void printcounts(void)  /* fcn */
/* Print how many times various functions have been called. */
{ printf(
    "setbkgd        %d\n"
    "consistify     %d\n"
    "consis9        %d\n"
    "proceed        %d\n"
    "backup         %d\n"
    "go             %d\n"
    "listneighbors  %d\n"
    "changecurr     %d\n"
    "nxgen          %d\n"
    "trycell        %d\n"
    "comporbackup   %d 000000\n",
    countsetbkgd, countconsistify, countconsis9, countproceed,
    countbackup, countgo, countlistneighbors, countchangecurr, countnxgen,
    counttrycell, countcomporbackuphi);
}
#endif

/* -------------------------------------------------------------------- */

boolean semifizzle(void)    /* fcn */
/* Checks to see if all changed cells in the current generation are
   within the DONTCOUNT region.
*/
{  static point *p;

#if COUNT
  if (TRACE) {printf("semifizzle(). chgd[%d]=%d. nays[%d]=%d.\n",
  gen,(int)(chgd[gen]-chglist),gen,(int)(nays[gen]-chglist));fflush(stdout);}
#endif

  for (p=chgd[gen]; p<nays[gen]; p++)
    if (!(flag[p->row][p->col] & DONTCOUNT))  return FALSE;
  return TRUE;
}

/************************************************************************/
/* These are the functions for handling regions of 4 different          */
/* shapes.  A region is either rectangular (specified by 2 opposite     */
/* corners), triangular (specified by its corners), a line segment      */
/* (specified by its endpoints), or a single point.                     */
/************************************************************************/

typedef struct
  { char regiontype;
    int minrow;             /* Bounding box of region */
    int mincol;
    int maxrow;
    int maxcol;
    int r0; int c0;         /* Vertices of region */
    int r1; int c1;
    int r2; int c2;
    int orient;             /* Orientation of triangle */
  }  region;

/* -------------------------------------------------------------------- */

void readregion(char *p, region *reg, char def)    /* fcn */
/* Reads region descriptor from p into reg.  A region descriptor has
   one of these forms:
         Rectangle:     r <r0> <c0> <r1> <c1>
         Triangle:      t <r0> <c0> <r1> <c1> <r2> <c2>
         Line segment:  l <r0> <c0> <r1> <c1>
         Point:         p <r0> <c0>
   The letter (r, t, l, or p) determines the region type and how many
   numbers will be read; if it is missing, then def specifies the type
   of region.
*/
{   while (*p == ' ')  p++;
    reg->regiontype = *p;
    if (isdigit(*p))  reg->regiontype = def;
    else p++;

    switch(reg->regiontype)
      { case 'p':
          if (sscanf(p, "%d %d", &reg->r0, &reg->c0) != 2)
            err("Readregion (point) error");
          reg->minrow = reg->maxrow = reg->r0;
          reg->mincol = reg->maxcol = reg->c0;
          break;

        case 'l':
        case 'r':
          if (sscanf(p, "%d %d %d %d",
                                &reg->r0, &reg->c0, &reg->r1, &reg->c1) != 4)
            err("Readregion (line or rectangle) error");
          reg->minrow = min(reg->r0, reg->r1);
          reg->mincol = min(reg->c0, reg->c1);
          reg->maxrow = max(reg->r0, reg->r1);
          reg->maxcol = max(reg->c0, reg->c1);
          break;

        case 't':
          if (sscanf(p, "%d %d %d %d %d %d",
            &reg->r0, &reg->c0, &reg->r1, &reg->c1, &reg->r2, &reg->c2) != 6)
            err("Readregion (triangle) error");

          reg->orient = reg->r0*(reg->c1-reg->c2) +
                        reg->r1*(reg->c2-reg->c0) +
                        reg->r2*(reg->c0-reg->c1);
          /* Sign of orient is orientation of triangle ABC */

          if (reg->orient==0)
            err("Readregion error: degenerate triangle descriptor.");

          reg->minrow = min(reg->r0, min(reg->r1, reg->r2));
          reg->maxrow = max(reg->r0, max(reg->r1, reg->r2));
          reg->mincol = min(reg->c0, min(reg->c1, reg->c2));
          reg->maxcol = max(reg->c0, max(reg->c1, reg->c2));
          break;

        default:
          err("Bad region type");
          break;
      }

    if (reg->minrow < 0)  reg->minrow = 0;
    if (reg->maxrow > HT-1)  reg->maxrow = HT-1;
    if (reg->mincol < 0)  reg->mincol = 0;
    if (reg->maxcol > WD-1)  reg->maxcol = WD-1;
    printf("Region type %c.  Rows %d to %d.  Cols %d to %d.\n",
      reg->regiontype, reg->minrow, reg->maxrow, reg->mincol, reg->maxcol);
}

/* -------------------------------------------------------------------- */

boolean pointinregion(int row, int col, region *reg)    /* fcn */
/* Is the specified point in the region? */
{   if (row < reg->minrow || row > reg->maxrow ||
        col < reg->mincol || col > reg->maxcol)  return FALSE;

    switch(reg->regiontype)
      { case 'p':
        case 'r':  return TRUE;

        case 'l':  return
          reg->r0 * (reg->c1 - col) +
          reg->r1 * (col - reg->c0) +
          row * (reg->c0 - reg->c1)   == 0;

        case 't':  return
          reg->orient * (reg->r0 * (reg->c1 - col) +
                         reg->r1 * (col - reg->c0) +
                         row  * (reg->c0 - reg->c1)   ) >= 0 &&
          reg->orient * (reg->r0 * (col - reg->c2) +
                         row  * (reg->c2 - reg->c0) +
                         reg->r2 * (reg->c0 - col)    ) >= 0 &&
          reg->orient * (row * (reg->c1 - reg->c2) +
                         reg->r1 * (reg->c2 - col) +
                         reg->r2 * (col - reg->c1)    ) >= 0;
      }
}

/************************************************************************/
/* Initialization functions, including command handlers.                */
/************************************************************************/

void initarrays(void)   /* fcn */
/* Initialize bkgd and curr to UNK with OFF boundaries, and flag
   to FALSE with TRUE boundaries.
*/
{ static int r,c;

  for (r=0; r<HT; r++)
    for (c=0; c<WD; c++)
       if (r <= 1 || r >= HT-2 || c <= 1 || c >= WD-2)
        { bkgd[r][c] = curr[r][c] = OFF;
          flag[r][c] = 0xFF;
        }
      else
        { bkgd[r][c] = curr[r][c] = UNK;
          flag[r][c] = 0;
        }

  for (r=1; r<HT-1; r++)
	for (c=1; c<WD-1; c++)
	  naysum[r][c] = bkgd[r-1][c-1] + bkgd[r-1][c] + bkgd[r-1][c+1] +
                     bkgd[r  ][c-1] +                bkgd[r  ][c+1] +
                     bkgd[r+1][c-1] + bkgd[r+1][c] + bkgd[r+1][c+1];
}

/* -------------------------------------------------------------------- */

void initconsistable(void)  /* fcn */
/* Initialize consistable. */
{ int val, oncount, offcount, unkcount, nbhd, i,
      survct, deathct, birthct, sterilect;

  for (val=OFF; val<=UNK; val=((val==ON) ? UNK : val+1))
  for (oncount=0; oncount<=8; oncount++)
  for (offcount=0; oncount+offcount<=8; offcount++)
    { unkcount = 8 - oncount - offcount;
      nbhd = oncount + unkcount*UNK;

      survct = deathct = birthct = sterilect = 0;

      if (val==ON || val==UNK)
        for (i=oncount; i<=oncount+unkcount; i++)
          if (rule[ON][i])  survct++;
          else              deathct++;

      if (val==OFF || val==UNK)
        for (i=oncount; i<=oncount+unkcount; i++)
          if (rule[OFF][i]) birthct++;
          else              sterilect++;

      if (sterilect==0 && survct==0)        consistable[val][nbhd] = INCONSIS;
      else if (val==UNK && sterilect==0)    consistable[val][nbhd] = CELLON;
      else if (val==UNK && survct==0)       consistable[val][nbhd] = CELLOFF;
      else if (unkcount==0)                 consistable[val][nbhd] = NOINFO;
      else
        if (val==OFF)
          if (sterilect==1 && !rule[OFF][oncount])
                                            consistable[val][nbhd] = NAYSOFF;
          else if (sterilect==1 && !rule[OFF][oncount+unkcount])
                                            consistable[val][nbhd] = NAYSON;
          else                              consistable[val][nbhd] = NOINFO;
        else if (val==ON)
          if (survct==1 && rule[ON][oncount])
                                            consistable[val][nbhd] = NAYSOFF;
          else if (survct==1 && rule[ON][oncount+unkcount])
                                            consistable[val][nbhd] = NAYSON;
          else                              consistable[val][nbhd] = NOINFO;
        else    /* val==UNK */
          if (sterilect==1 && !rule[OFF][oncount] &&
              survct==1 && rule[ON][oncount])
                                            consistable[val][nbhd] = NAYSOFF;
          else if (sterilect==1 && !rule[OFF][oncount+unkcount] &&
                    survct==1 && rule[ON][oncount+unkcount])
                                            consistable[val][nbhd] = NAYSON;
          else                              consistable[val][nbhd] = NOINFO;
          
    }
}

/* -------------------------------------------------------------------- */

void inittranstable(void)   /* fcn */
/* Initialize transtable. */
{ int val, oncount, offcount, unkcount, nbhd, i;
  boolean maybeon, maybeoff;

  for (val=OFF; val<=ON; val++)
  for (oncount=0; oncount<=8; oncount++)
  for (offcount=0; oncount+offcount<=8; offcount++)
    { unkcount = 8 - oncount - offcount;
      nbhd = oncount + unkcount*UNK;

      for (i=oncount, maybeon=maybeoff=FALSE; i<=oncount+unkcount; i++)
        if (rule[val][i])   maybeon = TRUE;
        else                maybeoff = TRUE;

      if (maybeon)
        if (maybeoff)   transtable[val][nbhd] = UNK;
        else            transtable[val][nbhd] = ON;
      else              transtable[val][nbhd] = OFF;
    }
}

/* -------------------------------------------------------------------- */

void setrule(char *p)   /* fcn */
/* Read rule.  E.g. the rule for Life could be given as "B3/S23" or
   "b3/s23" or "b3s23" or "3s23" or "3/23" or "s2b3S3" or ...
*/
{ int i,val;

  for (val=OFF; val<=ON; val++)
    for (i=0; i<=8; i++)  rule[val][i] = 0;

  val=OFF;  /* Default is birth */
  for (; *p; p++)
    { if (*p == 'b' || *p == 'B')       val=OFF;
      else if (*p == 's' || *p == 'S')  val=ON;
      else if (*p == '/')               val = 1-val;
      else if (*p < '0' || *p > '8')    err("Bad rule")
      else
        rule[val][*p - '0'] = 1;
    }
  inittranstable();
  initconsistable();
  printf("Rule set to B");
  for (i=0; i<=8; i++)  if (rule[OFF][i])  printf("%d",i);
  printf("/S");
  for (i=0; i<=8; i++)  if (rule[ON][i])  printf("%d",i);
  printf("\n");
}

/* -------------------------------------------------------------------- */

void readpattern(int r0, int c0)    /* fcn */
/* Read initial values of bkgd and curr and bits 0 of flag */
{ static int r,c,bg,cr;
  static char ch;

  r = r0;  c = c0;

  while (TRUE)
    { ch = getchar();

	  if (ch==':' || ch=='O' || ch=='s')	flag[r][c] |= DONTCHANGE;

	  if (var[127] && (ch=='.' || ch=='o' || ch=='0' || ch=='1' || ch=='?'))
		flag[r][c] |= DONTCOUNT;

      switch (ch)
        { case '.': case ':':	bg=0;  cr=0;  break;
          case 'o': case 'O':	bg=1;  cr=1;  break;
          case '0':				bg=1;  cr=0;  break;
          case '1':				bg=0;  cr=1;  break;
          case ',': case '?': case 's':		c++;	continue;
          case '\n':			r++;   c=c0;  continue;
          case '!':  while (getchar() != '\n');     return;
          default:   err("Bad character or EOF while reading");
        }

      if (bkgd[r][c] == (1^bg) || (bkgd[r][c] == UNK && proceed(r,c,bg,0)==ERR))
        { printf("readpattern error: bkgd[%d, %d] = %d.  bg = %d\n",
			r,c,(int)bkgd[r][c], bg);
          display(0);
          exit(1);
        }

      if (bg != cr)
        { nays[0]->row = r;  (nays[0]++)->col = c;
        }

      c++;
    }
}

/* -------------------------------------------------------------------- */

int readknownrotors(void)   /* fcn */
/* Reads a list of known rotor descriptors and corresponding names from
   the file "knownrotors".  Each rotor is defined by 1 or more lines of the
   file.  The first line begins with the rotor descriptor, which may
   continue on subsequent lines and/or end with a carriage return if it's
   large.  That's followed by 1 or more tabs.  The tabs are followed by the
   name of the oscillator.  The rotor descriptors and names are read into
   the arrays knownrotor and name.  The number of rotors read is returned.
*/
{ FILE *knownrotorsfile;
  int i;
  boolean readingname;
  char ch, *p;

  if ((knownrotorsfile = fopen(knownrotorsfilename, "r")) == NULL)
    err("Can't open known rotors file.\n");

  i=0;
  readingname = FALSE;
  p = knownrotor[i] = knownrotorsandnames;

  while (ch=fgetc(knownrotorsfile), !feof(knownrotorsfile))
    if (ch == '\t')
      { *p++ = 0;                   /* Mark end of rotor */
        readingname = TRUE;         /* Switch to reading name */
        name[i] = p;
      }
    else if (ch == '\n')
      { if (readingname)
          { *p++ = 0;               /* Mark end of name */
            readingname = FALSE;    /* Switch to reading rotor */
            if (++i >= MAXKNOWN)
              err("Too many known rotors; increase MAXKNOWN.");
            knownrotor[i] = p;
          }
      }
    else    *p++ = ch;

  *p = 0;                           /* Mark end of last name */
  knownrotor[i] = NULL;             /* Mark end of list */
  if (p > knownrotorsandnames + MAXFILESIZE - 5)
    err("Known rotors file is too big.  Increase MAXFILESIZE.\n");
  return i;
}

/* -------------------------------------------------------------------- */

void listcommands(void) /* fcn */
/* Print list of commands */
{ printf(   "Commands are:\n\n"
            "h#          Set max height of changed region\n"
            "w#          Set max width of changed region\n"
            "c#          Set max number of changed cells\n"
            "s#          Set random number seed\n"
            "P#          Set probability for free choices\n"
            "H#          Set height of space (<=81, usually 80 or 81)\n"
            "W#          Set width of space (<=81, usually 80 or 81)\n"
            "Rb###/s###  Set rule\n\n"

            "r# #        Read bkgd and curr\n"
            "C<region>   Clear bkgd and curr within region\n\n"

            "d0<region>  Disallow changes within region\n"
            "D0<region>  Allow changes within region\n"
            "d1<region>  Don't count changes within region\n"
            "D1<region>  Do count changes within region\n\n"

            "skipstable  Don't print stable outcomes (except fizzles)\n"
			"nopics      Don't show patterns, just rotor descriptors\n"
            "skipfizzle  Don't print fizzle outcomes\n"
            "showfin     Show finished patterns\n\n"

            "nosymm      No symmetry\n"
            "horsymm     Symmetry across horizontal line\n"
            "vertsymm    Symmetry across vertical line\n"
            "diagsymm    Symmetry across NW-SE dIagonal\n"
            "rot90symm   90 degree rotational symmetry\n"
            "rot180symm  180 degree rotational symmetry\n"
            "plussymm    Symmetry across horizontal and vertical lines\n"
            "xsymm       Symmetry across both diagonal lines\n"
            "fullsymm    Full symmetry\n\n"

            "v# #        Set variable (for program modifications)\n"
            "?           Print this list\n"
            ";<text>     Comment\n"
            "\n"
        );
}

/* -------------------------------------------------------------------- */

void docommand(char *p, boolean cmdline)    /* fcn */
/* Do the command pointed to by p.  The command may have been given as
   a command line argument (if cmdline is TRUE) or read from standard
   input.  In the case of a read command, additional lines will be
   read to complete the command, so this can't be done from the
   command line.
*/
{ int i,r,c,bitposition,varnum;
  long int varval;
  region reg;
  char cmd;

  cmd = *p;
  if (cmd == ';')
	printf("%s\n",p);   /* Comment */

  else if (!strcmp(p, "?"))    listcommands();

  else if (!strcmp(p, "nosymm"))        SYMM = NOSYMM;
  else if (!strcmp(p, "horsymm"))       SYMM = HORSYMM;
  else if (!strcmp(p, "vertsymm"))      SYMM = VERTSYMM;
  else if (!strcmp(p, "diagsymm"))      SYMM = DIAGSYMM;
  else if (!strcmp(p, "rot90symm"))     SYMM = ROT90SYMM;
  else if (!strcmp(p, "rot180symm"))    SYMM = ROT180SYMM;
  else if (!strcmp(p, "plussymm"))      SYMM = PLUSSYMM;
  else if (!strcmp(p, "xsymm"))         SYMM = XSYMM;
  else if (!strcmp(p, "fullsymm"))      SYMM = FULLSYMM;

  else if (!strcmp(p, "skipstable"))    SKIPSTABLE = TRUE;
  else if (!strcmp(p, "nopics"))        NOPICS = TRUE;
  else if (!strcmp(p, "skipfizzle"))    SKIPFIZZLE = TRUE;
  else if (!strcmp(p, "showfin"))       SHOWFIN = TRUE;
  else if (!strcmp(p, "showall"))       SHOWALL = TRUE;

  else if (*p == 'c' && *(p+1) == 'c')					/* "cc" command */
	{ p+=2;
	  while (*p == ' ')  p++;
	  for (numcc=0; *p; numcc++)
		{ reqchgcount[numcc] = atoi(p);
		  while (*p != ' ' && *p)  p++;
		  while (*p == ' ')  p++;
		}
	  printf("Initial change counts must be:");
	  for (i=0; i<numcc; i++)  printf(" %d", reqchgcount[i]);
	  printf("\n");
	}

  else if (cmd == 's')
    { KRsrand(atoi(p+1));
      printf("Random seed = %ld\n", atoi(p+1));
    }

  else if (cmd == 'P')  prob = atoi(p+1);
  else if (cmd == 'h')  maxheight = atoi(p+1);
  else if (cmd == 'w')  maxwidth = atoi(p+1);
  else if (cmd == 'c')  maxchng = atoi(p+1);
  else if (cmd == 'R')  setrule(p+1);

  else if (cmd == 'H')
    { HT = atoi(p+1);
      if (HT < 0 || HT > MAXHT)  err("Bad height");
      initarrays();
      printf("Height changed, so universe cleared.\n");
    }

  else if (cmd == 'W')
    { WD = atoi(p+1);
      if (WD < 0 || WD > MAXWD)  err("Bad width");
      initarrays();
      printf("Width changed, so universe cleared.\n");
    }

  else if (cmd == 'K')
    { strcpy(knownrotorsfilename, p+1);
    }

  else if (cmd == 'v')
    { sscanf(p+1, "%d %d", &varnum, &varval);
      if (varnum<0 || varnum>=NUMVARS)  err("Bad variable number");
      var[varnum] = varval;
	  printf("var[%d] = %ld\n",varnum,varval);
    }

  else if (cmd == 'r' && !cmdline)
    { sscanf(p+1, "%d %d", &r, &c);
      readpattern(r,c);
    }

  else if (cmd == 'C' || cmd == 'd' || cmd == 'D')
    { if (cmd != 'C')
        { bitposition = *(p+1) - '0';
          if (bitposition<0 || bitposition>7)
            err("Bad bit in 'd' or 'D' command.");
          if (cmd == 'd')  printf("Setting flag %d in region\n",bitposition);
          else  printf("Clearing flag %d in region\n",bitposition);
          p++;
        }
      readregion(p+1, &reg, 'r');
      for (r=reg.minrow; r<=reg.maxrow; r++)
      for (c=reg.mincol; c<=reg.maxcol; c++)
      if (pointinregion(r,c,&reg))
        switch(cmd)
          { case 'C':
              if (bkgd[r][c] == UNK && proceed(r,c,0,0)==ERR)
                { printf("Error in Clear command: r=%d c=%d\n",r,c);
                  display(0);
                  exit(1);
                }
              break;
            case 'd':
              flag[r][c] |= 1<<bitposition;
              break;
            case 'D':
              flag[r][c] &= ~(1<<bitposition);
              break;
          }
    }

  else  err1("Unknown command: %s", p);
}

/************************************************************************/
/* The main program.  It first reads the file "knownrotors" and         */
/* initializes things based on commands in the command line and         */
/* standard input.  Then it calls computecellorbackup repeatedly.       */
/* Whenever a generation is finished, it checks for periodicity,        */
/* printing the pattern if appropriate; otherwise it calls              */
/* listneighbors to set things up for computing the next generation.    */
/************************************************************************/

int main(int argc, char *argv[])   /* fcn */
{ int per, lines, i, g;
  unsigned long int h;
  char buff[500];
  boolean semifzl, toomanytwos, pervar111;
  point *p;

  setrule("B3/S23");
  initarrays();
  nwstng = settinglist;
  chgd[0] = nays[0] = chglist;

  /* Read command line arguments */
  if (argc==1)
    printf("Type '?' for list of commands\n\n");

  for (i=1; i<argc; i++)
    docommand(argv[i], TRUE);

  /* Read commands from standard input.  Input is ended by either */
  /* an empty line or end of file.                                */
  while (gets(buff) && buff[0])
    docommand(buff, FALSE);

  /* Print info about search */
  printf("Reading file '%s'\n",knownrotorsfilename);
  fflush(stdout);
  lines = readknownrotors();
  printf("%d known rotors read\n\n",lines);
  fflush(stdout);

  printf("Max height of changed region = %d\n", maxheight);
  printf("Max width of changed region = %d\n", maxwidth);
  printf("Max number of changed cells = %d\n", maxchng);
  printf("Probability = %d\n", prob);

  if (SKIPSTABLE)  printf("Skipping stable outcomes\n");
  if (NOPICS)      printf("Not showing pictures\n");
  if (SKIPFIZZLE)  printf("Skipping fizzle outcomes\n");
  if (SHOWFIN)     printf("Showing finished patterns\n");
  if (SHOWALL)     printf("Showing all gens\n");

  printf("Height = %d\n", HT);
  printf("Width = %d\n", WD);

  if (SYMM==HORSYMM)            printf("Horizontal symmetry\n");
  else if (SYMM==VERTSYMM)      printf("Vertical symmetry\n");
  else if (SYMM==DIAGSYMM)      printf("Diagonal symmetry\n");
  else if (SYMM==ROT90SYMM)     printf("90 degree rotational symmetry\n");
  else if (SYMM==ROT180SYMM)    printf("180 degree rotational symmetry\n");
  else if (SYMM==PLUSSYMM)      printf("Horizontal & vertical symmetry\n");
  else if (SYMM==XSYMM)         printf("Symmetry across both diagonals\n");
  else if (SYMM==FULLSYMM)      printf("Full symmetry\n");

  /* Check squareness for some symmetries */
  if (HT != WD &&
      (SYMM==DIAGSYMM || SYMM==ROT90SYMM || SYMM==XSYMM || SYMM==FULLSYMM))
    err("Symmetry requires height=width");

  nwstng = settinglist;     /* Make initialization un-backup-able */
  changecurr(curr, 0);      /* Change curr to gen 0 */

  found = FALSE;
  listneighbors(0); /* Init list at nays[0], chgd[1] */
  nay = nays[0];
  chg = chgd[gen = 1];
  hashcount=0;
  findchgcount(0);

  display(0);
  printf("Beginning search\n");

  while (computecellorbackup() == OK)
    { if (countcomporbackuplo == 0 &&
	       (countcomporbackuphi<50 || countcomporbackuphi%10 == 0))
        { printf("computecellorbackup calls: %d 000000\n",
		    countcomporbackuphi);
			display(0);		/* Temporary zzz */
		  dispchgcts(gen-1);		/* Added 5/1/2008 */
		  fflush(stdout);
		}

#if COUNT
if (countcomporbackuplo == 0)  printcounts();
if (TRACE)  {printf("back from computecellorbackup()\n");fflush(stdout);}
#endif

      if (nay == chgd[gen])     /* Done with this gen? */
        { if (gen>maxgenreached)
            { printf("maxgenreached = %d\n", maxgenreached=gen);
			  dispchgcts(gen-1);
              display(0);
			  fflush(stdout);
            }
          nays[gen] = chg;
		  findchgcount(gen);
          changecurr(curr, gen-1);
          changecurr(curr, gen);

          if (var[110] && gen>=var[110])
	        { for (g=gen-var[110], toomanytwos=TRUE; g<=gen; g++)
		        if (chgcount[g] != 2) { toomanytwos=FALSE; break; }
			  if (toomanytwos)
				{ found=TRUE;
				  printf("Too many 2s\n");
				  gen++;
				  continue;
				}
	        }

		  if (var[121] && gen==var[121] && chgcount[gen]==var[122])
            { found=TRUE;
              printf("chgcount[%d] can't be %d\n",gen,var[122]);
              gen++;
              continue;
            }

		  if (var[123])
			if (gen<var[123] && chgcount[gen]!=var[124])
              { found=TRUE;
                gen++;
                continue;
              }
			else if (gen==var[123] && chgcount[gen]==var[124])
              { found=TRUE;
                printf("chgcount[%d] can't be %d\n",gen,var[124]);
                gen++;
                continue;
              }

		  if (gen<numcc &&
			  ((reqchgcount[gen] >= 0 &&
									nays[gen]-chgd[gen] != reqchgcount[gen]) ||
			   (reqchgcount[gen] < 0 &&
									nays[gen]-chgd[gen] == -reqchgcount[gen])))
			{ found=TRUE;
			  /*printf("Enforcing chgcount[%d] = %d\n",gen,reqchgcount[gen]);*/
			  gen++;
			  continue;
			}

		  if (var[128] && gen>=var[128] &&
			   !(chgcount[gen] == chgcount[gen-var[128]] &&
				 ((width[gen] == width[gen-var[128]] &&
					 height[gen] == height[gen-var[128]]) ||
				  (width[gen] == height[gen-var[128]] &&
					 height[gen] == width[gen-var[128]]))))
			{ found=TRUE;
			  gen++;
			  continue;
			}

		  if (var[116] && gen>=18 && chgcount[gen]==1 &&
		    chgcount[gen-1]==3 && chgcount[gen-2]==3 && chgcount[gen-3]==5 &&
			chgcount[gen-4]==6 && chgcount[gen-5]==4 && chgcount[gen-6]==4 &&
			chgcount[gen-7]==4 && chgcount[gen-8]==3 && chgcount[gen-9]==1 &&
			chgcount[gen-10]==3 && chgcount[gen-11]==3 &&
			chgcount[gen-12]==5 && chgcount[gen-13]==6 &&
			chgcount[gen-14]==4 && chgcount[gen-15]==4 &&
			chgcount[gen-16]==4 && chgcount[gen-17]==3 &&
			chgcount[gen-18]==1)
            { found=TRUE;
              printf("5c/9 continues too long\n");
              gen++;
              continue;
            }

          if (var[111] && gen >= var[111])
            { for (p=chgd[gen-var[111]], pervar111=TRUE;
                   p<nays[gen-var[111]];   p++)
                if (p->row >= var[132] && p->col >= var[133] &&
				    p->row <= var[134] && p->col <= var[135] &&
					bkgd[p->row][p->col] == curr[p->row][p->col])
                  { pervar111=FALSE; break; }
			  if (!pervar111)
				{ found=TRUE;
				  gen++;
				  continue;
				}
			}

		  if (var[120] && chgcount[gen] && chgcount[gen]<var[120])
            { found=TRUE;
              /* printf("Change count too small\n"); */
              gen++;
              continue;
            }


          h = hash();

          if (nays[gen] == chgd[gen])
            { found = TRUE;
              if (hashnew(h))
                if (!SKIPFIZZLE)
                  { printf("*****  Fizzle at gen %d\n",gen);
					printoscinfo(gen, 'f');
					dispchgcts(gen);
                    display(0);
					if (SHOWALL)
					  for (g=1; g<=gen; g++)  display(g);
                  }
            }

          else if (per = period())
            { found = TRUE;
              semifzl = semifizzle();
              if (hashnew(h))
                if (per>1)
                  { if (var[131] && per==3)
				      hashcount--;			/* Forget about p3 stuff */
					else
				      { printf("*****  Period %d at gen %d%s\n",
                          per, gen-per, semifzl ? " (semifizzle)" : "");
                        printoscinfo(per, 'p');
                        printoscinfo(gen, 'u');
					    dispchgcts(gen);
					  }
					if (!var[131] || per>6)
					  display(0);
                    if (SHOWFIN)  display(gen);
					if (SHOWALL)
					  for (g=1; g<=gen; g++)  display(g);
                  }
                else
                  if (semifzl || !SKIPSTABLE)
                    { printf("*****  Stable at gen %d%s\n",
                        gen-1, semifzl ? " (semifizzle)" : "");
                      printoscinfo(gen, 's');
					  dispchgcts(gen);
                      display(0);
                      if (SHOWFIN)  display(gen);
					  if (SHOWALL)
					    for (g=1; g<=gen; g++)  display(g);
                    }
            }

          else if (gen == MAXGEN)
            { found = TRUE;
              if (hashnew(h))
                { printf("*****  Max gen (%d) reached\n", MAXGEN);
					dispchgcts(gen);
                  display(0);
                  if (SHOWFIN)  display(gen);
				  if (SHOWALL)
					for (g=1; g<=gen; g++)  display(g);
                }
            }

          else
            { listneighbors(gen);
              nay = nays[gen];
              chg = chgd[gen+1];
            }

          gen++;
		  if (var[129])  agesm[gen] = 0;
        }
    }

#if COUNT
  printcounts();
#endif

  printf("computecellorbackup calls: %d %06d\n",
    countcomporbackuphi, countcomporbackuplo);
  printf("No more objects\n");
}
