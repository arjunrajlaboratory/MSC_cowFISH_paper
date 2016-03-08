/* Import - MSC gels*/
PROC IMPORT OUT=WORK.MSC_gel_data DATAFILE="/folders/myfolders/MSCDiffCore.csv" 
        DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;

proc glimmix data=MSC_gel_data pconv=.0002 plots=all ic=q;
	nloptions maxiter=50;
    class Days_CM Donor Media;
  
    model cy_RNA=Media Days_CM Media*Days_CM/  dist=poisson link=log ddfm=satterth  solution;
    random Intercept  Media / subject=Donor solution;
		lsmeans Media*Days_CM / diff adjust=simulate(seed=0) adjdfe=row ilink;
	
		LSMESTIMATE Media*Days_CM 
   		"CM- vs CM+, Day 01" 1 -1 0 0 0 0 0 0 0 0,
		"CM- vs CM+, Day 04" 0 0 1 -1 0 0 0 0 0 0,
		"CM- vs CM+, Day 07" 0 0 0 0 1 -1 0 0 0 0,
		"CM- vs CM+, Day 14" 0 0 0 0 0 0 1 -1 0 0, 
		"CM- vs CM+, Day 21" 0 0 0 0 0 0 0 0 1 -1, 
		"CM+, Day 1 vs Day 4" 1 0 -1 0 0 0 0 0 0 0,
		"CM+, Day 1 vs Day 7" 1 0 0 0 -1 0 0 0 0 0,
		"CM+, Day 1 vs Day 14" 1 0 0 0 0 0 -1 0 0 0,
		"CM+, Day 1 vs Day 21" 1 0 0 0 0 0 0 0 -1 0,
		"CM+, Day 4 vs Day 7" 0 0 1 0 -1 0 0 0 0 0 ,
		"CM+, Day 4 vs Day 14" 0 0 1 0 0 0 -1 0 0 0,
		"CM+, Day 4 vs Day 21" 0 0 1 0 0 0 0 0 -1 0,
		"CM+, Day 7 vs Day 14" 0 0 0 0 1 0 -1 0 0 0,
		"CM+, Day 7 vs Day 21" 0 0 0 0 1 0 0 0 -1 0,
		"CM+, Day 14 vs Day 21" 0 0 0 0 0 0 1 0 -1 0 / adjust=simulate(seed=1) adjdfe=row;
    
proc glimmix data=MSC_gel_data pconv=.0002 plots=all ic=q;
	nloptions maxiter=50;
    class Days_CM Donor Media;
  
    model GAPDH_RNA=Media Days_CM Media*Days_CM/  dist=poisson link=log ddfm=satterth  solution;
    random Intercept  Media / subject=Donor solution;
		lsmeans Media*Days_CM / diff adjust=simulate(seed=0) adjdfe=row ilink;
	
		LSMESTIMATE Media*Days_CM 
   		"CM- vs CM+, Day 01" 1 -1 0 0 0 0 0 0 0 0,
		"CM- vs CM+, Day 04" 0 0 1 -1 0 0 0 0 0 0,
		"CM- vs CM+, Day 07" 0 0 0 0 1 -1 0 0 0 0,
		"CM- vs CM+, Day 14" 0 0 0 0 0 0 1 -1 0 0, 
		"CM- vs CM+, Day 21" 0 0 0 0 0 0 0 0 1 -1, 
		"CM+, Day 1 vs Day 4" 1 0 -1 0 0 0 0 0 0 0,
		"CM+, Day 1 vs Day 7" 1 0 0 0 -1 0 0 0 0 0,
		"CM+, Day 1 vs Day 14" 1 0 0 0 0 0 -1 0 0 0,
		"CM+, Day 1 vs Day 21" 1 0 0 0 0 0 0 0 -1 0,
		"CM+, Day 4 vs Day 7" 0 0 1 0 -1 0 0 0 0 0 ,
		"CM+, Day 4 vs Day 14" 0 0 1 0 0 0 -1 0 0 0,
		"CM+, Day 4 vs Day 21" 0 0 1 0 0 0 0 0 -1 0,
		"CM+, Day 7 vs Day 14" 0 0 0 0 1 0 -1 0 0 0,
		"CM+, Day 7 vs Day 21" 0 0 0 0 1 0 0 0 -1 0,
		"CM+, Day 14 vs Day 21" 0 0 0 0 0 0 1 0 -1 0 / adjust=simulate(seed=1) adjdfe=row;