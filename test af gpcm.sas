%include 'c:\dropbox\kach\sas\macro\simu_irt.sas';
data itempar;
input item_no name $ slope step1 step2 step3 step4 step5 choices;
cards;
1 item1 1 -2.5 -1.0 1.5 2.0 . 5
2 item2 1 -1.5 -1.0 1.0 1.5 . 5
3 item3 1 -1.5 -0.5 0.0 2.0 . 5
4 item4 1 -1.5 -1.0 0.5 2.0 . 5
;
%simu_irt(par_fil=work.itempar, outfile=simu, n_respo=1000, distribu=normal(0), seedresp=1, 
standard=NO, model=GPCM, min_code=0);

data in;
input item_no item_name $ item_text $ max;
datalines;
1 item1 item1 4
2 item2 item2 4
3 item3 item3 4
4 item4 item4 4
;
run;


%gpcm(data=simu, ITEM_names=in, ICC=yes, plotcat=YES, plotmean=YES);

