/* Import - MSC division data*/
PROC IMPORT OUT=WORK.MSC_div_data DATAFILE="/folders/myfolders/MSCdivision_data.csv"
       DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;

/*Linear mixed model to compare means*/
ods graphics on;

/*random effects -- USE THIS FORMAT*/
proc glimmix data=MSC_div_data pconv=.0002 plots=all ic=q;
	nloptions maxiter=50;
    class TimeCategory;
  
    model cy_difference=TimeCategory/  dist=poisson link=log ddfm=satterth  solution;
		lsmeans TimeCategory / diff adjust=simulate(seed=0) adjdfe=row ilink;
	run;
	
	proc glimmix data=MSC_div_data pconv=.0002 plots=all ic=q;
	nloptions maxiter=50;
    class TimeCategory;
  
    model GAPDH_difference=TimeCategory/  dist=poisson link=log ddfm=satterth  solution;
		lsmeans TimeCategory / diff adjust=simulate(seed=0) adjdfe=row ilink;
	run;
