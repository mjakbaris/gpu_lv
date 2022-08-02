/* Ordinary differential equation solver, Runge-Kutta-England technique.
   Copyright © 1988 Free Software Foundation, Inc.
   François Pinard <pinard@iro.umontreal.ca>, 1988.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 1, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software Foundation,
   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/



#include <stdio.h>
#include <math.h>
#include <stdlib.h>

typedef struct struct_rke_variables
{

    /* The following are saved from rke_init call arguments. */

    int n_equations;		/* Number of simultaneous equations */
    int (*eval_routine) ();	/* Routine to compute derivatives */

    /* These may be changed by the user between two solve calls. */

    double minimum_step;	/* Minimum allowable step size */
    double maximum_step;	/* Maximum allowable step size */
    double current_step;	/* Current integration step size */
    double error_slope;		/* Slope of maximum error per time unit */
    double error_biais;		/* Biais of maximum error per time unit */
    int accepted_steps;		/* Accumulated number of accepted steps */
    int rejected_steps;		/* Accumulated number of rejected steps */
}
rke_variables;


extern double fabs ();

/* Initialize a new system of equations. */

void rke_init (int number, int (*routine) () , rke_variables * var)	/* Newly allocated reentrancy block */
{

  (var)->n_equations = number;
  (var)->eval_routine = routine;
  (var)->minimum_step = 0.000001;
  (var)->maximum_step = 1000000.0;
  (var)->current_step = 1.0;
  (var)->error_slope = 0.0000001;
  (var)->error_biais = 0.00000001;
  (var)->accepted_steps = 0;
  (var)->rejected_steps = 0;

  return;
}

/* Terminate a set of equations. */

void rke_term (rke_variables * var)
{
  free (var);
}


/* Main routine of the module, ODE solver. */

/* Perform a consistent move of time in the system. */

int	 rke_solve (rke_variables * var, double *time, double variables[], double aimed_time, double * scrap_mem)
{
  double whole_step;		/* Signed integration step size */
  double quarter_step;		/* 0.25 * whole_step */
  double half_step;		/* 0.50 * whole_step */
  double three_quarter_step;	/* 0.75 * whole_step */
  double estimated_error;	/* Error as estimated by England method */
  double allowable_error;	/* Maximum error that user permits */
  int within_tolerance;		/* Allowable error has not been passed */
  int all_errors_small;		/* All errors within 2% of tolerances */
  int length_of_array;		/* Length of temporary arrays, is bytes */
  int k;			/* Index in various arrays */

  double *dp, *vt, *v, *d;
  double *a1, *a2, *a3, *a4, *a5, *a6, *a7;

  /* Allocate the work arrays. */

  length_of_array = var->n_equations * sizeof (double);

  dp = &scrap_mem[0];
  vt = &scrap_mem[length_of_array];
  v = &scrap_mem[2*length_of_array];
  d = &scrap_mem[3*length_of_array];
  a1 = &scrap_mem[4*length_of_array];
  a2 = &scrap_mem[5*length_of_array];
  a3 = &scrap_mem[6*length_of_array];
  a4 = &scrap_mem[7*length_of_array];
  a5 = &scrap_mem[8*length_of_array];
  a6 = &scrap_mem[9*length_of_array];
  a7 = &scrap_mem[10*length_of_array];

  /* The integration will continue if a minimum step could bring the
     system closer to the time that is aimed for, even if we have to
     overshoot it a little. */

  while (2 * fabs (aimed_time - *time) > var->minimum_step)
  {

      /* Evaluate initial step size and direction. */

    if ((whole_step = aimed_time - *time) > 0.0)
    {
     if (whole_step > var->current_step)
       whole_step = var->current_step;
   }
   else
   {
     if (whole_step < - var->current_step)
       whole_step = - var->current_step;
   }

      /*  Evaluate initial differentials. */

   if (! (*(var->eval_routine)) (*time, variables, dp))
     return 0;

   do

	/* Loop integrating at this time point until integration error is
	   within tolerances.  In any case, adjust integration step size. */

   {
	  /* Calculate various step sizes. */

     quarter_step = 0.25 * whole_step;
     half_step = quarter_step + quarter_step;
     three_quarter_step = half_step + quarter_step;

	  /* Perform a partial computation for one step of Runge-Kutta
	     4th order integration, as far as necessary to chain it to
	     England method for estimating integration errors. */

     for (k = 0; k < var->n_equations; ++k)
     {
       a1[k] = half_step * dp[k];
       v[k] = variables[k]
       + 0.5*a1[k];
     }

     if (! (*(var->eval_routine)) (*time + quarter_step, v, d))
       return 0;

     for (k = 0; k < var->n_equations; ++k)
     {
       a2[k] = half_step * d[k];
       v[k] = variables[k]
       + 0.25 * (a1[k] + a2[k]);
     }

     if (! (*(var->eval_routine)) (*time + quarter_step, v, d))
       return 0;

     for (k = 0; k < var->n_equations; ++k)
     {
       a3[k] = half_step * d[k];
       v[k] = variables[k]
       + (-a2[k] + a3[k] + a3[k]);
     }

     if (! (*(var->eval_routine)) (*time + half_step, v, d))
       return 0;

     for (k = 0; k < var->n_equations; ++k)
     {
       a4[k] = half_step * d[k];
       vt[k] = variables[k]
       + (a1[k] + 4.0*a3[k] + a4[k]) / 6.0;
     }

     if (! (*(var->eval_routine)) (*time + half_step, vt, d))
       return 0;

     for (k = 0; k < var->n_equations; ++k)
     {
       a5[k] = half_step * d[k];
       v[k] = vt[k]
       + 0.5*a5[k];
     }

     if (! (*(var->eval_routine)) (*time + three_quarter_step, v, d))
       return 0;

     for (k = 0; k < var->n_equations; ++k)
     {
       a6[k] = half_step * d[k];
       v[k] = vt[k]
       + 0.25*(a5[k] + a6[k]);
     }

     if (! (*(var->eval_routine)) (*time + three_quarter_step, v, d))
       return 0;

     for (k = 0; k < var->n_equations; ++k)
     {
       a7[k] = half_step * d[k];
       v[k] = variables[k]
       + (-a1[k] - 96.0*a2[k] + 92.0*a3[k] - 121.0*a4[k]
         + 144.0*a5[k] + 6.0*a6[k] - 12.0*a7[k]) / 6.0;
     }

	  /* Perform England error analysis on partial integration. */

     if (! (*(var->eval_routine)) (*time + whole_step, v, d))
       return 0;

     within_tolerance = 1;
     all_errors_small = 1;

     for (k = 0; k < var->n_equations; ++k)
     {
       estimated_error
       = fabs ((-a1[k] + 4.0*a3[k] + 17.0*a4[k]
         - 23.0*a5[k] + 4.0*a7[k] - half_step*d[k])
       / 90.0);
       allowable_error = fabs (whole_step)
       * (var->error_slope*fabs (vt[k]) + var->error_biais);
       if (estimated_error > allowable_error)
       {
        within_tolerance = 0;
        break;
      }
      else if (estimated_error > 0.02 * allowable_error)
        all_errors_small = 0;
    }
    if (within_tolerance)
    {
	    ++(var->accepted_steps);

	      /* Complete the Runge-Kutta step and return values. */

	    for (k = 0; k < var->n_equations; ++k)
		    v[k] = vt[k] + (-a6[k] + a7[k] + a7[k]);
	    if (! (*(var->eval_routine)) (*time + whole_step, v, d))
		    return 0;
	    *time += whole_step;
	    for (k = 0; k < var->n_equations; ++k)
		    variables[k] = vt[k]
		    + (a5[k] + 4.0*a7[k] + half_step*d[k]) / 6.0;

	      /* Increment step size if desirable. */
	    if (all_errors_small && fabs (whole_step) == var->current_step){
		    if (2 * var->current_step > var->maximum_step){
			    var->current_step = var->maximum_step;
		    }
		    else{
			    var->current_step *= 2;
		    }
	    }
    }
	else
	{
		++var->rejected_steps;
		/* Decrement step size if possible. */
		if (fabs (whole_step) > var->minimum_step)
		{
			if (var->current_step < 2 * var->minimum_step)
				var->current_step = var->minimum_step;
			else
				var->current_step *= 0.5;
			if (aimed_time > *time)
				whole_step = var->current_step;
			else
				whole_step = - var->current_step;
		}
		else
			return 0;	/* Convergence failed */
	}
   }
	  while (!within_tolerance);
	  {
		  return 1;			/* Convergence succeeded */
	  }




/* Check how close we can get back to our initial conditions. */

void print_return (double back,double initial)
{
  printf ("  returning to %12.6lf, got %12.6lf\n", initial, back);
  return;
}



/* Print statistics about number of steps. */

void print_steps(rke_variables * var)
{
  printf ("    using %3d accepted and %3d rejected steps\n",
    var->accepted_steps, var->rejected_steps);
  return;
}


/* Integration under a normal curve. */



static double example_1_const;	/* 1.0 / sqrt (2 * pi) */


static int problem_function_1 (  double t  , double v[1]  , double d[1])
{
  d[0] = example_1_const * exp (-0.5 * t * t);
  return 1;
}


static void example_1 ()
{

  rke_variables* p = malloc(sizeof(struct struct_rke_variables) );
  double * scrap_mem = (double *) malloc(11*sizeof(double));
  
  double t;
  double v[1];

  example_1_const = 1.0 / sqrt (2 * 3.1415926);

  rke_init (1, problem_function_1, p);

  t = -1.0;			/* Start at -1.0 */
  v[0] = 0.0;			/* Surface is 0.0 at this point */

  /* Now, simply move to +1.0, and collect the answer. */

  if (rke_solve (p, &t, v, 1.0 , scrap_mem))
    printf ("\nProbability	= %12.6lf.\n", v[0]);
  else
    printf ("\nProbability not computed, error.\n");
  print_steps (p);

  /* Just undo this, to see if we get back where we started. */

  if (rke_solve (p, &t, v, -1.0 , scrap_mem))
    print_return (v[0], 0.0);
  else
    printf ("  return to start not computed, error.\n");
  print_steps (p);

  rke_term (p);
  free(scrap_mem);
  return;
}



/* Rediscovering cos and sin. */

static int problem_function_2 (double t, double v[2], double d[2])
{
  d[0] = -v[1];
  d[1] = v[0];
  return 1;
}


static void example_2 ()
{
  rke_variables* p = malloc(sizeof(struct struct_rke_variables) );
  double * scrap_mem = (double *) malloc(2*11*sizeof(double));

  double t;
  double v[2];

  rke_init (2, problem_function_2 , p);

  t = 0.0;			/* Start where we know the values */
  v[0] = 1.0;			/* cos 0 = 1.0 */
  v[1] = 0.0;			/* sin 0 = 0.0 */

  /* Now, simply move to 1.5, and collect the answer. */

  if (rke_solve (p, &t, v, 1.5,scrap_mem))
    printf ("\ncos (1.5)	= %12.6lf.\n", v[0]);
  else
    printf ("\ncos (1.5) not computed, error.\n");
  print_steps (p);

  /* Just undo this, to see if we get back where we started. */

  if (rke_solve (p, &t, v, 0.0 , scrap_mem))
  {
    print_return (v[0], 1.0);
    print_return (v[1], 0.0);
  }
  else
    printf ("  return to start not computed, error.\n");
  print_steps (p);

  rke_term (p);
  free(scrap_mem);
  return;
}


/* Box slowing by friction in air. */

static int problem_function_3 (double t, double v[2], double d[2])
{
  d[0] = v[1];
  d[1] = -0.01 * v[1] * v[1];
  return 1;
}


static void example_3 () {
  
  rke_variables* p = malloc(sizeof(struct struct_rke_variables) );
  double * scrap_mem = (double *) malloc(2*11*sizeof(double));

  double t;
  double v[2];

  rke_init (2, problem_function_3, p);

  t = 0.0;			/* Start the clock... */
  v[0] = 0.0;			/* ... with no distance so far */
  v[1] = 100.0;		/* ... but some initial speed */

  /* Now, simply ask the clock to be 5.0, and collect the answer. */

  

  if (rke_solve (p, &t, v, 5.0 , scrap_mem))
    printf ("\nDistance	= %12.6lf.\n", v[0]);
  else
    printf ("\nDistance not computed, error.\n");
  print_steps (p);

  /* Just undo this, to see if we get back where we started. */

  if (rke_solve (p, &t, v, 0.0 , scrap_mem))
  {
    print_return (v[0], 0.0);
    print_return (v[1], 100.0);
  }
  else
    printf ("  return to start not computed, error.\n");
  print_steps (p);

  rke_term (p);
  free(scrap_mem);
  return;
}


/* Main program. */

int main ()
{
  example_1 ();
  example_2 ();
  example_3 ();
  return 0 ;
}
