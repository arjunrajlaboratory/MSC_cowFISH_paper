/* Import - MSC division data*/
PROC IMPORT OUT=WORK.CH_area  DATAFILE="/folders/myfolders/chondrocyteSpreadAreas.csv"
       DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
    run;
    
PROC IMPORT OUT=WORK.CH_volume  DATAFILE="/folders/myfolders/ChondrocyteDiameterData.csv"
       DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2;
    run;


/*Linear mixed model to compare means*/
ods graphics on;

	proc glm data=CH_area;
	    class Passage;
	    model AreaMicrons=Passage;
		lsmeans Passage / pdiff=all;
	    run;


	proc glm data=CH_volume;
	    class Passage;
	    model volumeCubicMicrons=Passage;
		lsmeans Passage / pdiff=all;
	    run;


proc glimmix data =CH_volume;
class Passage;
model volumeCubicMicrons=Passage ;
random _residual_/group=Passage; /* This gives the treatments separate variance estimates  */
covtest 'common variance' homogeneity;  /* This does a likelihood ratio test for homogeneity */
lsmeans Passage / adjust=simulate(nsamp=50000 seed=7654321 report); /* Borrowing from 1zmm's code for control of multiple comparisons */
run;

/* https://communities.sas.com/message/152724 */

proc glimmix data =CH_area;
class Passage;
model AreaMicrons=Passage ;
random _residual_/group=Passage; /* This gives the treatments separate variance estimates  */
covtest 'common variance' homogeneity;  /* This does a likelihood ratio test for homogeneity */
lsmeans Passage / adjust=simulate(nsamp=50000 seed=7654321 report); /* Borrowing from 1zmm's code for control of multiple comparisons */
run;
