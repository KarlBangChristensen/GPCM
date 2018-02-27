/**************************************************************************
This macro simulates responses to dichotomous and polytomous items from 
two-parameter IRT models, the generalized partial credit model (GPCM) or 
the graded response model. It requires an input file, as specified below

The following statement is used to call the macro

%simu_irt(	par_file=work.itempar, 
		outfile= simu, 
		n_resp=1000,
		dist=normal(0), 
		seedresp=1,
		model=GPCM)

par_file: 

contains the item parameters and number of response options for each item
any directory reference (e.g. SASUSER.itempar) can be used each record 
contains information about one item. Must contain the following variables
with these exact variable names:

 NAME     The name of the items
 SLOPE    The item slopes (discrimination parameter)
 STEP1    First item step parameter (item parameter)
 STEP2    Second item step parameter (missing for CHOICES<3)
 STEP3    Third item step parameter (missing for CHOICES<4)
 STEP4    Fourth item step parameter (missing for CHOICES<5)
 STEP5    Fifth item step parameter (missing for CHOICES<6)
 CHOICES  Number of item response choices

(set irrelevant item step parameters to missing [e.g. for an item with three
response choices, set STEP3 and higher to missing])

outfile:    

is the name of the output data file

n_resp:

number of responses to be simulated. Default is 1000 

dist:   

specifies how the latent variable is distributed. Default option is 
NORMAL(0) generating a normal distribution with mean zero and SD one 
(using the SAS function NORMAL with time of day as seed. To specify 
another mean and/or SD write, e.g. (normal(0)*[STD]+[MEAN])

model:

specifies IRT model. Valid options are GPCM (default) and GRM

****************************************************************************/

%macro simu_irt(par_file, 
		outfile,
		seedresp=0,
		n_resp=1000, 
                dist=normal(0), 
                model= GPCM);

options nomprint nomlogic nosymbolgen nonotes nostimer;
*options mprint mlogic symbolgen notes stimer;

* Check model statement;

%if %upcase(&model) ne GPCM and %upcase(&model) ne GRM %then %do;
%PUT;
  %PUT;
  %PUT **********************************************************;
  %PUT ;
  %PUT Model must be GPCM or GRM;
  %PUT ;
  %PUT **********************************************************;
  %goto exit;
%end;

* save number of items, item names and n umber of response options as macro variables;

data _null_;
	set &par_file end=last;
	array steps (5) step1-step5;
	call symput('item'||trim(left(put(item_no,4.))), name);
	call symput('itch'||trim(left(put(item_no,4.))), left(choices));
	call symput('itsl'||trim(left(put(item_no,4.))), left(slope));
  	do i=1 TO 5;
    		call symput ('st'||trim(left(put(item_no,4.)))||'_'||left(i),left(steps(i)));
  	end;
	if last then call symput('nitems', trim(left(put(item_no,4.))));
run;

* generate theta;

data &outfile;
	do person= 1 to &n_resp;
		theta = &dist.;
		output;
	end;
run;

* simulate item responses;

data &outfile (
	keep=person theta %do _i = 1 %to &nitems; &&item&_i %end;
	);
	set &outfile;
	%do _j=1 %to &nitems;
		* GPCM;
		%if %upcase(%left(%trim(&model)))=GPCM %then %do;
			divisor=1;
			do i=1 to (symget("itch&_j")-1);
				sum=0;
				do j=1 to i;
					sum=sum+symget("st&_j._"||left(j));
				end;
				divisor=divisor+exp(symget("itsl&_j")*(i*theta-sum));
			end;
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
		* GRM;
		%if %upcase(%left(%trim(&model)))=GRM %then %do;
			array icca&_j (*) icc1-icc&&itch&_j;
			array ciccb&_j (*) cicc1-cicc&&itch&_j;
			do i=symget("itch&_j") to 1 by -1;
				j=i-1;
				if i = symget("itch&_j") then do;
					icca&_j(i) = 1/(1+exp(-&&itsl&_j*(theta- symget("st&_j._"||left(j)))));
					ciccb&_j(i) = 1/(1+exp(-&&itsl&_j*(theta- symget("st&_j._"||left(j)))));
				end;
				else if i = 1 then icca&_j(i) = 1 - ciccb&_j(i+1);
				else do;
					icca&_j(i) = 1/(1+exp(-&&itsl&_j*(theta-symget("st&_j._"||left(j)))))-ciccb&_j(i+1) ;
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
		&&item&_j = 0 %do _k = 1 %to %eval(&&itch&_j-1); + (_p&_j ge cicc&_k) %end;;
	%end;
run;

%exit:
options notes stimer;
%mend simu_irt;
