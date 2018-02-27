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

      The macro requires one input file, as specified below

      The following statement is used to call the macro

           %simu_irt(par_fil=work.itempar
                   , outfile= simu
                   , n_respo=1000
                   , distribu=normal(0)
				   , seedresp= 1
                   , standard=YES
                   , model=GPCM
                   , min_code = 1
                    )


       par_fil:  this file contains the item parameters and number of response
                 response choices for each item
                 you can use any directory reference that you would use in
                 a data= statement
                 in the par_fil each record (line) is information for one item
                 the par_fil should at least contain the following variables
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


%macro simu_irt(par_fil, outfile, seedresp= 0, n_respo=1000, 
                distribu=normal(0), standard= YES,
                model= GPCM, min_code=1);

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
/*
         This program part saves each item name in the par_fil
         as a separate macro variable, saves the number of response choices
         for each item and calculates the number of items
                                                                    */

data _null_;
set &par_fil end=last;

array steps (5) step1-step5;

call symput('item'||trim(left(put(item_no,4.))), name);
call symput('itch'||trim(left(put(item_no,4.))), left(choices));
call symput('itsl'||trim(left(put(item_no,4.))), left(slope));
  do i=1 TO 5;
    call symput ('st'||trim(left(put(item_no,4.)))||'_'||left(i),
                                                left(steps(i)));
  end;

if last then
    call symput('nitems', trim(left(put(item_no,4.))));
run;




data &outfile;
do person= 1 to &n_respo;

**************************************************************************;
*;
*        Generate theta - distribution function specified by user
*;
**************************************************************************;

theta = &distribu;

output;
end;

run;

%if %upcase(%trim(%left(&standard))) = YES %then %do;

proc means data= &outfile noprint;
var theta;
output out= _meanstd mean=mean std=std;
run;

data _NULL_;
set _meanstd;
call symput('_mean', mean);
call symput('_std', std);
run;

data &outfile;
set &outfile;
 
theta = (theta-&_mean)/&_std;

run;

%end;


data &outfile (keep= person theta 
%do _i = 1 %to &nitems;
       &&item&_i
%end; 
               );
set &outfile;

%do _j=1 %to &nitems;
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
    do i=1 to (symget("itch&_j")-1);
       sum=0;
         do j=1 to i;
           sum=sum+symget("st&_j._"||left(j));
         end;
       divisor=divisor+exp(symget("itsl&_j")*(i*theta-sum));
    end;

/*   Then we calculate the icc for each response choice          */

  icc1=1/divisor;

  array icca&_j (*) icc2-icc&&itch&_j;
      do i=1 to (symget("itch&_j")-1);
       sum2=0;
         do j=1 to i;
           sum2=sum2+symget("st&_j._"||left(j));
         end;
       icca&_j(i)=exp(symget("itsl&_j")*(i*theta-sum2))/divisor;
      end;

cicc1 = icc1;
  array cicca&_j (*) cicc1-cicc&&itch&_j;
       do k = 2 to &&itch&_j;
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

  array icca&_j (*) icc1-icc&&itch&_j;
  array ciccb&_j (*) cicc1-cicc&&itch&_j;

  do i=symget("itch&_j") to 1 by -1;
    j=i-1;
    if i = symget("itch&_j") then 
      do;
        icca&_j(i) = 1/(1+exp(-&&itsl&_j*(theta- symget("st&_j._"||left(j)))));
        ciccb&_j(i) = 1/(1+exp(-&&itsl&_j*(theta- symget("st&_j._"||left(j)))));
      end;
    else if i = 1 then
      icca&_j(i) = 1 - ciccb&_j(i+1);
    else
      do;
        icca&_j(i) = 1/(1+exp(-&&itsl&_j*(theta-symget("st&_j._"||left(j)))))
          - ciccb&_j(i+1) ;
        ciccb&_j(i) = 1/(1+exp(-&&itsl&_j*(theta- symget("st&_j._"||left(j)))));
      end;
  end;
  drop cicc1-cicc&&itch&_j;

cicc1 = icc1;
  array cicca&_j (*) cicc1-cicc&&itch&_j;
       do k = 2 to &&itch&_j;
         cicca&_j(k) = cicca&_j(k-1) + icca&_j(k);
       end;

%end;

_p&_j = ranuni(&seedresp);

&&item&_j = &min_code                    %do _k = 1 %to %eval(&&itch&_j-1);
                      + (_p&_j ge cicc&_k)
                                          %end;
                                           ;

%end;


run;


%exit:

options notes stimer;

%mend simu_irt;
