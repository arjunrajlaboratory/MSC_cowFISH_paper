/* Import - Chondrocyte Gels*/
PROC IMPORT OUT=WORK.CH_gel_data DATAFILE="/folders/myfolders/ChondRediffCore.csv" 
        DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;

/* Aggrecan*/
	
/*One-way ANOVA with homogeneity of variance test*/
/*  	p ~ .01 */
proc glm data=CH_gel_data;
	    class experimentid;
	    model cy_RNA=experimentid;
	    means experimentid / hovtest welch;
	    run;
/*Linear mixed model to compare means*/
ods graphics on;

/* aggrecan */
proc glimmix data=CH_gel_data pconv=.0002 plots=all;
nloptions maxiter=50;
    class Days_CM Donor Passage;
    model cy_RNA=Passage Days_CM Passage*Days_CM/  dist=poisson link=log ddfm=satterth  solution;
    random Intercept Passage/ subject=Donor solution;
	lsmeans Passage*Days_CM / diff adjust=simulate(seed=1)  adjdfe=row ilink;
	
	LSMESTIMATE Passage*Days_CM 
		"Day 1, P0 vs P5" 1 -1 0 0,
		"P0, Day 1 vs 14" 1 0 -1 0,
		"Day 14, P0 vs P5" 0 0 1 -1, 
		"P5, Day 1 vs 14" 0 1 0 -1 / adjust=simulate(seed=1) adjdfe=row;
    /* https://communities.sas.com/message/205151 */
    /* http://www.unt.edu/rss/class/Jon/SAS_SC/LMM/SAS_Module9.2.htm*/
    /*http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_glimmix_sect028.htm*/
run;


/*GAPDH*/
proc glimmix data=CH_gel_data pconv=.0002 plots=all;
nloptions maxiter=50;
    class Days_CM Donor Passage;
    model GAPDH_RNA=Passage Days_CM Passage*Days_CM/  dist=poisson link=log ddfm=satterth  solution;
    random Intercept Passage/ subject=Donor solution;
	lsmeans Passage*Days_CM / diff adjust=simulate(seed=1)  adjdfe=row ilink;
	
	LSMESTIMATE Passage*Days_CM 
		"Day 1, P0 vs P5" 1 -1 0 0,
		"P0, Day 1 vs 14" 1 0 -1 0,
		"Day 14, P0 vs P5" 0 0 1 -1, 
		"P5, Day 1 vs 14" 0 1 0 -1 / adjust=simulate(seed=1) adjdfe=row;
    /* https://communities.sas.com/message/205151 */
    /* http://www.unt.edu/rss/class/Jon/SAS_SC/LMM/SAS_Module9.2.htm*/
    /*http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_glimmix_sect028.htm*/
run;


/*aggrecan/GAPDH per cell*/  
proc glimmix data=CH_gel_data pconv=.0002 plots=all;
nloptions maxiter=50;
    class Days_CM Donor Passage;
    model AGGperGAPDH=Passage Days_CM Passage*Days_CM/  dist=poisson link=log ddfm=satterth  solution;
    random Intercept Passage/ subject=Donor solution;
	lsmeans Passage*Days_CM / diff adjust=simulate(seed=1)  adjdfe=row ilink;
	
	LSMESTIMATE Passage*Days_CM 
		"Day 1, P0 vs P5" 1 -1 0 0,
		"P0, Day 1 vs 14" 1 0 -1 0,
		"Day 14, P0 vs P5" 0 0 1 -1, 
		"P5, Day 1 vs 14" 0 1 0 -1 / adjust=simulate(seed=1) adjdfe=row;
    /* https://communities.sas.com/message/205151 */
    /* http://www.unt.edu/rss/class/Jon/SAS_SC/LMM/SAS_Module9.2.htm*/
    /*http://support.sas.com/documentation/cdl/en/statug/63962/HTML/default/viewer.htm#statug_glimmix_sect028.htm*/
run;