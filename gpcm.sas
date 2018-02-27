/**************************************************************************

macro to compute item parameters and population parameters in a generalized 
partial credit model (GPCM). Conditional maximum likelihood parameters from
a Rasch model are computed and used as starting values

%gpcm(data, item_names, out=GPCM, write=NO, ICC=NO, PLOTCAT=NO, PLOTMEAN=NO, nsimu=30, nquadr=30);

DATA: the data set containing the items (scored 0,1, .. ,'max'), the number of response categories 'max' can differ between items.  

ITEM_NAMES: data set that contains information about the items. This data set should contain the variables 

		item_name: item names
		item_text (optional): item text for plots 
		max: maximum score on the item. 

OUT: the name (default MML) given to output files. 

		Item parameter estimates are put in a file 'out_itempar'.
		Maximum of the log likelihood is put in a file 'out_logl'. 
		Estimated values of the person locations are put in a file 'out_theta'. 

WRITE: If 'write' is YES the macro writes to the log window. 

ICC: If 'ICC' is YES the macro plots item characteristic curves 

PLOTCAT: Indicates whether observed and expected trace lines should be plotted.

PLOTMEAN: If 'PLOTMEAN' is YES observed and item means are plotted. These observed and expected mean 
scores across the range of score values are plotted together with simulated values sets. 
'nsimu' specifies the number of simulated data sets.

NOTE: item names should not be more than eight characters.

**************************************************************************/

%macro gpcm(data, item_names, out=GPCM, write=YES, ICC=NO, PLOTCAT=NO, PLOTMEAN=NO, nsimu=30, nquadr=30);
	options nonotes nostimer;
	options nomprint;
	ods listing close;

	/* save number of respondents as macro variable */
	data _null_; 
		set &data end=final; 
	 	if final then call symput('N',trim(left(_N_))); 
	run;

	/* numbering items - find the maximum length of values - max(max) is the highest number of response categories */
	proc sql noprint;
		select count(distinct(item_name)), sum(max), max(max), max(length(item_name)), 
					length(strip(put(count(distinct(item_name)),5.))), length(strip(put(max(max),5.)))
		into :_nitems, :_maxtotscore, :_max_max, :_l_item_name, :_l_item_no, :_l_max
		from &item_names;
	quit;

	%let _max_max=&_max_max.;
	%let _nitems=&_nitems.;

	proc sql noprint;
		select distinct(item_name), max
		into :_item1-:_item&_nitems., :_max1-:_max&_nitems.
		from &item_names;
	quit;

	%do _i=1 %to &_nitems.; %let _max&_i=&&_max&_i; %end;

	/* data set with item info - make format for a given variable by finding maximum length of the corresponding values */
	data _item_names;
		format _item_no &_l_item_no.. _item_name $&_l_item_name.. _max &_l_max..;
		%do _i=1 %to &_nitems.; _item_no="&_i."*1; _item_name="&&_item&_i"; _max="&&_max&_i"*1; output;	%end;
	run;
	/* make data set including item scores */
	data _item_scores (rename=(_i=_score));
		format _item_score $%eval(&_l_item_name.+&_l_max.+1).;
		set _item_names;
		do _i=0 to _max; _item_score=strip(_item_name)||'|'||strip(_i); output; end;
	run;

	%let maxscore=&_maxtotscore;

	/* check score distribution */
	%do _sc=0 %to &maxscore; %let _n&_sc=0; %end;
	data _NULL_; 
		set &data; 
		t=0 %do _i=1 %to &_nitems; +&&_item&_i %end;; 
		%do _sc=0 %to &maxscore; if t=&_sc then call symput("_n&_sc",1); %end; 
	run;

	/**********************************/
	/* item parameter estimation part */
	/**********************************/
	
	/* check score distribution */
	%do _sc=0 %to &maxscore; %let _n&_sc=0; %end;

	proc sql;
		create table _t0 as select 
		%do _i=1 %to &_nitems; &&_item&_i , %end;
		count(*) as count
		from &data
		group by %do _i=1 %to &_nitems-1; &&_item&_i , %end; &&_item&_nitems
		;
	quit;

	data _t0; set _t0; length _all_ 3; run;

	/* create (potentially very large) data set with all possible response patterns */
	data _patterns;
		%do _i=1 %to &_nitems; length &&_item&_i 3; %end;
		%do _i=1 %to &_nitems; do &&_item&_i=0 to &&_max&_i; %end;
		output;
		%do i=1 %to &_nitems.; end; %end;
	run;

	data _table; merge _t0 _patterns; by %do _i=1 %to &_nitems; &&_item&_i %end;; if n=. then n=0; run;

	data _table; set _table; length t 3; t=&_item1 %do _i=2 %to &_nitems; +&&_item&_i %end;; run;

	%do _it=1 %to &_nitems; 
		data _table; 
			set _table;
			%do _cat=1 %to &&_max&_it; length v&_it&_cat 3; %end; 
			%do _cat=1 %to &&_max&_it; v&_it&_cat=(&&_item&_it=&_cat); %end; 
		run;
	%end;

	data _table;
		set _table;
		%do _sc=0 %to &maxscore; length n&_sc 3; %end;
		%do _sc=0 %to &maxscore; n&_sc=(t=&_sc); %end;
	run;

	/* Starting values: fit Rasch model using CML */
	data _table; set _table; if count=. then count=0.00000000001; run;

	proc genmod data=_table; 
		ods output estimates=_est;
		ods output parameterestimates=_pf;
		model count=
			%do _it=1 %to &_nitems; %do _cat=1 %to &&_max&_it; v&_it&_cat %end; %end;
	 		%do _sc=1 %to %eval(&maxscore-1); %if &&_n&_sc=1 %then %do; n&_sc %end; %end;/d=p noint link=log maxiter=1000
		; 
		/* ESTIMATE statements to estimate item parameters (eta's) */
		%do _it=1 %to &_nitems; 
			%do _c=1 %to %eval(&&_max&_it-1); 
				estimate "eta&_it&_c" v&_it&_c %eval(&_maxtotscore) 
					%do _i0=1 %to &_nitems; v&_i0.&&_max&_i0. -&_c %end;
				;
			%end; 
			estimate "eta&_it.&&_max&_it" v&_it.&&_max&_it %eval(&_maxtotscore-&&_max&_it) 
				%do _i0=1 %to %eval(&_it-1); v&_i0.&&_max&_i0. -&_c %end;
				%do _i0=%eval(&_it+1) %to &_nitems; v&_i0.&&_max&_i0. -&_c %end;
			;
		%end; 
	run;

	/* save item parameter estimates (eta's) as macro variables */
	data _null_; 
		set _est; 
		estimate=LBetaEstimate/&_maxtotscore;
		%do _it=1 %to &_nitems; 
			call symput("eta&_it.0",trim(left(0))); 
			%do _cat=1 %to &&_max&_it; 
				if label="eta&_it&_cat" then call symput("eta&_it.&_cat",trim(left(estimate)));
			%end; 
		%end;
	run;

	/* write to log file */
	%if %upcase(%left(%trim(&write)))=YES %then %do;
	 	%put ---------------------------------------------------------;
		%put computing start values for item parameters;
	 	%put ---------------------------------------------------------;
	 	%put CML estimation - reading data set &data, &N respondents;
	 	%put &_nitems items; 
	 	%do _it=1 %to &_nitems; %put - &&_item&_it (%eval(&&_max&_it+1) response categories); %end; 
	 	%put ---------------------------------------------------------;
	 	%put item parameter estimates (from CML);
	 	%put ---------------------------------------------------------;
	 	%do _it=1 %to &_nitems; %put &&_item&_it; 
	 	%do _cat=1 %to &&_max&_it; %put cat &_cat &&eta&_it.&_cat; %end; %end;
	 	%put ---------------------------------------------------------;
	 	%put used as starting values, calling PROC NLMIXED to fit GPCM;
	 	%put ---------------------------------------------------------;
	%end;

	/* Fit GPCM using NLMIXED, new data set - one line for each response */
	data _new; 
		set &data;
		%do _it=1 %to &_nitems; 
			item="&&_item&_it                 "; 
			value=&&_item&_it; 
			person=_N_; 
			output; 
		%end;
	run;

	data _newnew; set _new; run;

	/* direct output to files */
	ods listing close;
	ods output nlmixed.parameterestimates=_pe;
	ods output nlmixed.additionalestimates=_ae;
	ods output nlmixed.fitstatistics=_lf;

	/* numerical maximization using previous estimates as start values (PROC NLMIXED) */

	proc nlmixed data=_new;
		parms 
		%do _i=1 %to &_nitems; 
			alpha&_i=1, %do _h=1 %to &&_max&_i; eta&_i&_h=&&eta&_i&_h, %end; 
		%end;;
		_theta=epsilon;
		%do _i=1 %to &_nitems; 
			_denom=1 %do _k=1 %to &&_max&_i; +exp(alpha&_i*(&_k*_theta+eta&_i&_k)) %end;;
			if item="&&_item&_i" and value=0 then ll=-log(_denom);
			%do _h=1 %to &&_max&_i; 
				if item="&&_item&_i" and value=&_h then ll=alpha&_i*(&_h*_theta+eta&_i&_h)-log(_denom);
			%end; 
		%end;
		model value~general(ll);
		random epsilon ~ normal(0,1) subject=person;
		/* estimate thresholds (PCM parametrized item parameters) */
		%do _i=1 %to &_nitems; 
			estimate "&&_item&_i.._slope" alpha&_i;
			%do _h=1 %to &&_max&_i; 
				estimate "&&_item&_i|&_h" -eta&_i&_h 
					%if &_h>1 %then %do; +eta&_i%eval(&_h-1) %end;; 
			%end; 
		%end;
	run;

	/* create data set maximum log likelihood value */
	data &out._logl; 
		set _lf; 
		if Descr='-2 Log Likelihood';
	run;

	/* create two data sets with item parameter estimates */
	data &out._itempar_ci; 
		set _ae;
		keep label estimate lower upper;
	run;

	data _ip0; set _ae; 
		if label in ("&_item1._slope" %do _i=2 %to &_nitems; ,"&&_item&_i.._slope" %end;);
		%do _i=1 %to &_nitems; 
			if label="&&_item&_i.._slope" then do; 
				name="&&_item&_i"; 
				choices=%eval(&&_max&_i+1); 
				slope=estimate; 
			end; 
		%end;
		keep name choices slope;
	run;
	%do _r=1 %to &_max_max;
		data _ip&_r; set _ae; 
			if label in ("&_item1.|&_r" %do _i=2 %to &_nitems; ,"&&_item&_i..|&_r" %end;);
			%do _i=1 %to &_nitems; if label="&&_item&_i..|&_r" then do; name="&&_item&_i"; step&_r=estimate; end; %end;
			keep name step&_r;
		run;
	%end;
	data &out._itempar; 
		merge %do _r=0 %to &_max_max; _ip&_r %end;; 
	run;

	/* save item parameters as macro variables */
	data _null_;
		set _pe;
		%do _it=1 %to &_nitems; 
			if trim(parameter)="alpha&_it" then call symput("_a&_it",trim(left(estimate)));		
			%do _h=1 %to &&_max&_it; 
				if trim(parameter)="eta&_it.&_h" then call symput("_e&_it.&_h",trim(left(estimate)));		
			%end; 
		%end;
	run;

	/* print parameter estimates */
	ods listing; 
		title "GPCM (&out): Item parameters";
		proc print data=&out._itempar noobs round; run;
		title "GPCM (&out): Maximum value of loglikelihood";
		proc print data=&out._logl noobs; run;
	ods listing close;


	/* range of thresholds */

	OPTIONS MPRINT;

	%do _i=1 %to &_nitems;
		proc sql noprint; select min(step1-step&&_max&_i), min(step1-step&&_max&_i) as :__min&_i, :__max&_i from &out._itempar; quit;
		%put;
		%put __min&_i = &&__min&_i;
		%put __max&_i = &&__max&_i;
	%end;

	OPTIONS NOMPRINT;

	/*****************************************/
	/* end of item parameter estimation part */
	/*****************************************/


	/******************************************/
	/* person location estimation part        */
	/******************************************/
	
	/*
	data _data; 
		set &data; 
		person=_N_; 
		_score=0 %do _it=1 %to &_nitems; +&&_item&_it %end;;
	run;
			
	data _data; set _data; if (_score>0 and _score<&_maxtotscore); run;

	* save number of respondents as macro variable;
	data _null_; 
		set _data end=final; 
	 	if final then call symput('_nontrivial',trim(left(_N_))); 
	run;

	%put &_nontrivial persons with a non trivial score;

	* MLE;
	proc nlmixed data=_data tech=newrap; 
		parms theta=0;
		%do _it=1 %to &_nitems;			
			_denom&_it=1 %do _x=1 %to &&_max&_it; +exp(&&_a&_it*(&_x*theta+&&_e&_it.&_x)) %end;;
		%end;
		%do _it=1 %to &_nitems;
			if &&_item&_it=0 then _lognum&_it=1; 
			%do _x=1 %to &&_max&_it; 
				if &&_item&_it=&_x then _lognum&_it=&&_a&_it*(&_x*theta+&&_e&_it.&_x); 
			%end;
		%end;
		* log likelihood;
		logl=0 %do _it=1 %to &_nitems; + _lognum&_it - log(_denom&_it) %end;;
		model person ~ general(logl);
		by person;
		%do _ppar=1 %to  &_nontrivial;
			ods output bygroup&_ppar..parameterestimates=_ppar&_ppar;
		%end;
	run;

	*/


	/******************************************/
	/* end of person location estimation part */
	/******************************************/

	/*************************************************************/
	/* ICC part -  plot item characteristic curves (if required) */
	/*************************************************************/
	

	%if %upcase(%left(%trim(&ICC)))=YES %then %do; 
		%put plotting item characteristic curves - ICCs (data set _plot);
	 	%put ---------------------------------------------------------;
		data _plot;
			array dens [101] dens1-dens101; 
			array norm [101] norm1-norm101; 
			sum_dens=0;
			do i=1 to 101; 
				norm[i]=-4+8*(i-1)/100; 
				dens[i]=exp(-(norm[i])**2); 
				sum_dens=sum_dens+dens[i]; 
				theta=norm[i];
				/* loop over items: compute denominators */
				%do _it=1 %to &_nitems; 
					_denom&_it=1 %do _k=1 %to &&_max&_it; +exp(&&_a&_it*(&_k*theta+&&_e&_it.&_k)) %end;;
					_prob&_it.0=1/_denom&_it;
					%do _h=1 %to &&_max&_it; 
						_prob&_it.&_h=exp(&&_a&_it*(&_h*theta+&&_e&_it.&_h))/_denom&_it;
					%end;
				%end;
				output;
			end;
			label theta='Latent variable'; 
			%do _it=1 %to &_nitems;	%do _h=0 %to &&_max&_it; 
				label _prob&_it.&_h='probability'; 
			%end; %end;
			keep %do _it=1 %to &_nitems; 
				_prob&_it.0 
				%do _h=1 %to &&_max&_it; _prob&_it.&_h %end; 
			%end;
			theta;
		run;
		/* symbol statements  */
		%do _h=1 %to %eval(&_max_max+1);	
			symbol&_h v=none c=blue i=join l=1 w=3; 
		%end;
		/* plotting */
		ods listing;
		%do _it=1 %to &_nitems;	
			title "GPCM (&out): Item characteristic curves (ICC) &&_item&_it";
			axis1 order=-4 to 4 by 1 length=15 cm value=(H=2) minor=NONE label=(H=2 'Latent variable');
			axis2 length=10 cm value=(H=2) minor=NONE label=(H=2 A=90 'Probability');
			proc gplot data=_plot;
				plot (%do _h=0 %to &&_max&_it; _prob&_it.&_h %end;)*theta/overlay haxis=axis1 vaxis=axis2;
			run; quit;
		%end;
		ods listing close;
	%end;


	/* FITPLOT part (if required) */
	%if (%upcase(%left(%trim(&PLOTCAT)))=YES or %upcase(%left(%trim(&PLOTMEAN)))=YES) %then %do; 
		%put computing observed and expected numbers - data set _fit;
		%put simulating &nsimu data sets;
		/* empirical score distribution */
		ods listing close; 

		proc sort data=_table; by t; run;
	 	%do _sc=0 %to &maxscore; 
	  		ods output means.bygroup%eval(&_sc+1).summary=_n&_sc; 
	  		proc means data=_table sum; var count; by t; run;
	  		data _null_; 
	   			set _n&_sc;
	   			call symput('N'||trim(left(&_sc)),trim(left(count_sum))); 
	  		run;
	 	%end;
	 	data _allprob;
			set _table;
			%do _it=1 %to &_nitems; 
				v&_it.0=1 %do _k=1 %to &&_max&_it; -v&_it.&_k %end;;
			%end;	
			keep t %do _it=1 %to &_nitems; v&_it.0 v&_it.1-v&_it.&&_max&_it %end;;
		run;
 
	 	data _allprob;
			set _allprob;
			array dens [101] dens1-dens101; 
			array norm [101] norm1-norm101; 
			sum_dens=0;
			do i=1 to 101; 
				norm[i]=-4+8*(i-1)/100; 
				dens[i]=exp(-(norm[i])**2)/sqrt(2*3.141593); 
			end;
			sum=0; 
			/* loop over quadrature points */
			do i=1 to 101;
				_theta=0+1*norm[i];
				/* loop over items: compute denominators */
				%do _it=1 %to &_nitems; 
					_denom&_it=1 %do _k=1 %to &&_max&_it; +exp(&&_a&_it*(&_k*_theta+&&_e&_it.&_k)) %end;;
				%end;
				logp=0
				%do _it=1 %to &_nitems; 
					%do _h=1 %to &&_max&_it; 
						+v&_it.&_h*(&&_a&_it*(&_h*_theta+&&_e&_it.&_h))
					%end; 
					-log(_denom&_it)
				%end;;
				prob=exp(logp);
				probdens=prob*dens[i];
				output; 
			end; 
			drop dens1-dens101 norm1-norm101 _denom1-_denom&_nitems sum_dens i sum logp;
		run;
		/* the density values phi(theta) */
	 	data _density;
			array dens [101] dens1-dens101; 
			array norm [101] norm1-norm101; 
			sum_dens=0;
			do i=1 to 101; 
				norm[i]=-4+8*(i-1)/100; 
				dens[i]=exp(-(norm[i])**2)/sqrt(2*3.141593); 
			end;
			sum=0; 
			/* loop over quadrature points */
			do i=1 to 101;
				_theta=0+1*norm[i];
				_density=dens[i];
				output; 
			end; 
			keep _theta _density;
		run;

		/* probabilities P(T=t|theta) */
		proc means data=_allprob mean; 
			var probdens;
			class t _theta; 
			output out=_scoreprob mean=_scoreprob;
		run;
		data _scoreprob; 
			set _scoreprob; 
			where t ne . and _theta ne .;
			drop _TYPE_ _FREQ_; 
		run;

		/* _prob P(X_i=k,T=t|theta), _condprob P(X_i=k|T=t,theta)=P(X_i=k,T=t|theta)/P(T=t|theta) */
		%do _it=1 %to &_nitems; %do _k=0 %to &&_max&_it;
			proc means data=_allprob mean;
				where v&_it.&_k=1;
				var probdens;
				class t _theta; 
				output out=_prob mean=_prob;
			run;
			data _prob&_it.&_k; 	
				merge _prob(where=(t ne . and _theta ne .)) _scoreprob; 
				_condprob=_prob/_scoreprob;
				drop _TYPE_ _FREQ_; 
			run;			
			proc sql;
				create table __prob&_it.&_k as select 
				a.t, mean(a._condprob*b._density) as _prob
				from _prob&_it.&_k a left join _density b
				on a._theta=b._theta
				group by a.t
				order by a.t;
			quit;
		%end; %end;

	/* end of PLOT part */
	%end;
	ods listing;
	/* clean up */
	proc datasets nolist;
		delete _table /*_ae*/ _pe _lf _new %do _r=0 %to &max; _ip&_r %end;
		%if %upcase(%left(%trim(&fittest)))=YES %then %do; 
			_prob _obs %do _it=1 %to &_nitems; %do _h=0 %to &max; _prob&_it.&_h _n&_it.&_h %end; %end;
		%end;;
	run; quit;
	title ' '; ods listing; options notes stimer;
%mend gpcm;

