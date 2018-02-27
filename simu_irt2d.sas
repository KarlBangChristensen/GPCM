/**************************************************************************


                  ************************************
                  *                                  *
                  *     Jakob Bue Bjorner            *
                  *     Institute of Public Health   *
                  *     University of Copenhagen     *
                  *     May    1999                  *
                  *                                  *
                  ************************************


      This macro simulated responses to dichotomous
      and polytomous items according to a two-item parameter item response
      model, the generalized partial credit model (GPCM) or the graded 
      response model.

      The macro simulates responses to two dimensions. Each dimension has a 
      normal distribution with specified mean and standard deviation. Also,
      the correlation between the two dimensions is specified.

      The macro requires one input file, as specified below

      The following statement is used to call the macro

           %simu_irt2d(parfile=work.itempar
                   , outfile= simu
                   , n_respo=1000
                   , items_d1= v1-v3
                   , items_d2= v4-v6
                   , mean_d1= 0
                   , mean_d2= 0
                   , std_d1= 1
                   , std_d2= 1
                   , correlation= .5
                   , seed_th= 0
				   , seedresp= 1
                   , model=GPCM
                   , min_code = 1
                    )


       parfile:  this file contains the item parameters and number of response
                 response choices for each item
                 you can use any directory reference that you would use in
                 a data= statement
                 in the parfile each record (line) is information for one item
                 the parfile should at least contain the following variables
                 with these exact variable names:

                 NAME     The name of the items
                 SLOPE    The item slopes (discrimination parameter)
                 STEP1    First item step parameter (item difficulty
                          parameter)
                 STEP2    Second item step parameter
                 .        (and so on, for an item with 6 response choices
                          there should be 5 item step parameters. In case you
                          have items with different number of response
                          choices, just set the irrelevant item step
                          parameters to missing (e.g. for an item with 3
                          response choices, set STEP3 and higher to missing))
                 CHOICES  Number of item response choices

       outfile    is the name of the output file

       n_respo    number of responses to be simulated. Default is 1000 

       distribu   specifies the way to generate the latent theta distribution.
                  The default option is NORMAL(0), this generates a normal
                  distribution with mean zero and std one (using the SAS 
                  function NORMAL with time of day as seed (see SAS LANGUAGE
                  MANUAL). To specify other means and standard deviations write
                  e.g. (normal(0)*[STD]+[MEAN])

       standard   if YES (default) the simulated data are standardized to 
                  make sure the mean is 0 and the standard deviation is 1

       model      specifies the statistical model. Valid options are GPCM
                  (the default) and GRM

       !! NB !!   The program use the following temporary file names, which
                  you should avoid:

                  full_icc
 

****************************************************************************/


%macro simu_irt2d(parfile, outfile, n_respo, items_d1, items_d2, correlation
              , mean_d1= 0, mean_d2= 0, std_d1= 1, std_d2= 1
              , seed_th= 0, seedresp= 0, model=GPCM, min_code = 1
                    );

options nomprint nomlogic nosymbolgen nonotes nostimer;
*options mprint mlogic symbolgen notes stimer;

**************************************************************************;
*;
*        Check model statement
*;
**************************************************************************;


%if %upcase(&model) ne GPCM and 
            %upcase(&model) ne GRM %then %do;
%PUT;
  %PUT;
  %PUT **********************************************************;
  %PUT ;
  %PUT NOTE!!!! Model should be GPCM or GRM;
  %PUT ;
  %PUT **********************************************************;
  %goto exit;
%end;



/*
         This program part saves each item name in the item list
         as a separate macro variable and calculates the number
         of items
                                                                    */

data _null_;
array _y (*) &items_d1;
length name $8;
do i=1 to dim(_y);
  call vname(_y{i}, name);
  call symput('_itemd1_'||trim(left(put(i,4.))), name);
end;
p=dim(_y);
call symput('_nitemsd1', trim(left(put(p,4.))));
run;

/*     Saving number of item response choices, item slope,
       and item steps as system variables
                                                                           */

data _null_;
set &parfile end=last;
retain match (0);

array steps (5) step1-step5;

%do _j = 1 %to &_nitemsd1;
if TRIM(upcase(name))= "%upcase(&&_itemd1_&_j)" then
   do;
     call symput("itchd1_&_j", left(choices));
     call symput("itsld1_&_j", left(slope));
         do i=1 TO 5;
         call symput ("st&_j"||'d1_'||left(i), left(steps(i)));
       end;
     match = match + 1;
   end;
%end;

/*        Check that the specified item names matches the items in the parfile
          if not, the program aborts                                        
*/

if last then
  do;
    if match < &_nitemsd1 then abort;
  end;

run;



data _null_;
array _y (*) &items_d2;
length name $8;
do i=1 to dim(_y);
  call vname(_y{i}, name);
  call symput('_itemd2_'||trim(left(put(i,4.))), name);
end;
p=dim(_y);
call symput('_nitemsd2', trim(left(put(p,4.))));
run;

/*     Saving number of item response choices, item slope,
       and item steps as system variables
                                                                           */

data _null_;
set &parfile end=last;
retain match (0);

array steps (5) step1-step5;

%do _j = 1 %to &_nitemsd2;
if TRIM(upcase(name))= "%upcase(&&_itemd2_&_j)" then
   do;
     call symput("itchd2_&_j", left(choices));
     call symput("itsld2_&_j", left(slope));
         do i=1 TO 5;
         call symput ("st&_j"||'d2_'||left(i), left(steps(i)));
       end;
     match = match + 1;
   end;
%end;

/*        Check that the specified item names matches the items in the parfile
          if not, the program aborts                                        
*/

if last then
  do;
    if match < &_nitemsd2 then abort;
  end;

run;

data &outfile;
do person= 1 to &n_respo;

**************************************************************************;
*;
*        Generate theta - distribution function specified by user
*;
**************************************************************************;

_x1 = normal(&seed_th);
_x2 = normal(&seed_th);

theta1 = (1-&correlation**2)** .5 * &std_d1 * _x1 + &correlation * &std_d1 * _x2 + &mean_d1;
theta2 = &std_d2 * _x2 + &mean_d2;

output;
end;

run;


data &outfile (keep= person theta1 theta2
%do _i = 1 %to &_nitemsd1;
       &&_itemd1_&_i
%end;
               );
set &outfile;

%do _j=1 %to &_nitemsd1;
%if %upcase(%left(%trim(&model)))= GPCM %then %do;


**************************************************************************;
*;
*        Generalized partial credit model
*;
**************************************************************************;

/*    Calculate ICC
      first we calculate the denominator,
      which is common for all response choices
                                                            */


  divisor=1;
    do i=1 to (symget("itchd1_&_j")-1);
       sum=0;
         do j=1 to i;
           sum=sum+symget("st&_j.d1_"||left(j));
         end;
       divisor=divisor+exp(symget("itsld1_&_j")*(i*theta1-sum));
    end;

/*   Then we calculate the icc for each response choice          */

  icc1=1/divisor;

  array icca&_j (*) icc2-icc&&itchd1_&_j;
      do i=1 to (symget("itchd1_&_j")-1);
       sum2=0;
         do j=1 to i;
           sum2=sum2+symget("st&_j.d1_"||left(j));
         end;
       icca&_j(i)=exp(symget("itsld1_&_j")*(i*theta1-sum2))/divisor;
      end;

cicc1 = icc1;
  array cicca&_j (*) cicc1-cicc&&itchd1_&_j;
       do k = 2 to &&itchd1_&_j;
         cicca&_j(k) = cicca&_j(k-1) + icca&_j(k-1);
       end;

%end;
%else %do;

**************************************************************************;
*;
*        Graded response model
*;
**************************************************************************;

/*   We calculate the icc for each response choice          */

  array icca&_j (*) icc1-icc&&itchd1_&_j;
  array ciccb&_j (*) cicc1-cicc&&itchd1_&_j;

  do i=symget("itchd1_&_j") to 1 by -1;
    j=i-1;
    if i = symget("itchd1_&_j") then 
      do;
        icca&_j(i) = 1/(1+exp(-&&itsld1_&_j*(theta1- symget("st&_j.d1_"||left(j)))));
        ciccb&_j(i) = 1/(1+exp(-&&itsld1_&_j*(theta1- symget("st&_j.d1_"||left(j)))));
      end;
    else if i = 1 then
      icca&_j(i) = 1 - ciccb&_j(i+1);
    else
      do;
        icca&_j(i) = 1/(1+exp(-&&itsld1_&_j*(theta1-symget("st&_j.d1_"||left(j)))))
          - ciccb&_j(i+1) ;
        ciccb&_j(i) = 1/(1+exp(-&&itsld1_&_j*(theta1- symget("st&_j.d1_"||left(j)))));
      end;
  end;
  drop cicc1-cicc&&itchd1_&_j;

cicc1 = icc1;
  array cicca&_j (*) cicc1-cicc&&itchd1_&_j;
       do k = 2 to &&itchd1_&_j;
         cicca&_j(k) = cicca&_j(k-1) + icca&_j(k);
       end;

%end;

_p&_j = ranuni(&seedresp);

&&_itemd1_&_j = &min_code                    %do _k = 1 %to %eval(&&itchd1_&_j-1);
                      + (_p&_j ge cicc&_k)
                                          %end;
                                           ;

%end;

run;


data &outfile (keep= person theta1 theta2
%do _i = 1 %to &_nitemsd1;
       &&_itemd1_&_i
%end;
%do _i = 1 %to &_nitemsd2;
       &&_itemd2_&_i
%end;
               );

set &outfile;

%do _j=1 %to &_nitemsd2;
%if %upcase(%left(%trim(&model)))= GPCM %then %do;


**************************************************************************;
*;
*        Generalized partial credit model
*;
**************************************************************************;

/*    Calculate ICC
      first we calculate the denominator,
      which is common for all response choices
                                                            */


  divisor=1;
    do i=1 to (symget("itchd2_&_j")-1);
       sum=0;
         do j=1 to i;
           sum=sum+symget("st&_j.d2_"||left(j));
         end;
       divisor=divisor+exp(symget("itsld2_&_j")*(i*theta2-sum));
    end;

/*   Then we calculate the icc for each response choice          */

  icc1=1/divisor;

  array icca&_j (*) icc2-icc&&itchd2_&_j;
      do i=1 to (symget("itchd2_&_j")-1);
       sum2=0;
         do j=1 to i;
           sum2=sum2+symget("st&_j.d2_"||left(j));
         end;
       icca&_j(i)=exp(symget("itsld2_&_j")*(i*theta2-sum2))/divisor;
      end;

cicc1 = icc1;
  array cicca&_j (*) cicc1-cicc&&itchd2_&_j;
       do k = 2 to &&itchd2_&_j;
         cicca&_j(k) = cicca&_j(k-1) + icca&_j(k-1);
       end;

%end;
%else %do;

**************************************************************************;
*;
*        Graded response model
*;
**************************************************************************;

/*   We calculate the icc for each response choice          */

  array icca&_j (*) icc1-icc&&itchd2_&_j;
  array ciccb&_j (*) cicc1-cicc&&itchd2_&_j;

  do i=symget("itchd2_&_j") to 1 by -1;
    j=i-1;
    if i = symget("itchd2_&_j") then 
      do;
        icca&_j(i) = 1/(1+exp(-&&itsld2_&_j*(theta2- symget("st&_j.d2_"||left(j)))));
        ciccb&_j(i) = 1/(1+exp(-&&itsld2_&_j*(theta2- symget("st&_j.d2_"||left(j)))));
      end;
    else if i = 1 then
      icca&_j(i) = 1 - ciccb&_j(i+1);
    else
      do;
        icca&_j(i) = 1/(1+exp(-&&itsld2_&_j*(theta2-symget("st&_j.d2_"||left(j)))))
          - ciccb&_j(i+1) ;
        ciccb&_j(i) = 1/(1+exp(-&&itsld2_&_j*(theta2- symget("st&_j.d2_"||left(j)))));
      end;
  end;
  drop cicc1-cicc&&itchd2_&_j;

cicc1 = icc1;
  array cicca&_j (*) cicc1-cicc&&itchd2_&_j;
       do k = 2 to &&itchd2_&_j;
         cicca&_j(k) = cicca&_j(k-1) + icca&_j(k);
       end;

%end;

_p&_j = ranuni(&seedresp*2 );

&&_itemd2_&_j = &min_code                    %do _k = 1 %to %eval(&&itchd2_&_j-1);
                      + (_p&_j ge cicc&_k)
                                          %end;
                                           ;
%end;

run;


%exit:

options notes stimer;

%mend simu_irt2d;
