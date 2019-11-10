1. Contributed by:

        Roy Smith
        Dept. of Electrical & Computer Engineering
        University of California,
	Santa Barbara, CA 93106
	U.S.A.
        roy@ece.ucsb.edu

2. Process/Description:

	The experiment is a simple SISO heating system.
	The input drives a 300 Watt Halogen lamp, suspended
	several inches above a thin steel plate.  The output
	is a thermocouple measurement taken from the back of
	the plate.

3. Sampling interval: 

	2.0 seconds

4. Number of samples

	801

5. Inputs:
        
        u: input drive voltage
        ...
6. Outputs:

        y: temperature (deg. C)
        ...
7. References:

	The use of this experiment and data for robust
	control model validation is described in:

        "Sampled Data Model Validation: an Algorithm and
        Experimental Application," Geir Dullerud & Roy Smith,
        International Journal of Robust and Nonlinear Control,
        Vol. 6, No. 9/10, pp. 1065-1078, 1996.

8. Known properties/peculiarities

	The data (and nominal model) is the above paper have the
	output expressed in 10's deg. C.  This has been rescaled 
	to the original units of deg. C. in the DaISy data set.
	There is also a -1 volt offset in u in the data shown plotted
	in the original paper.  This has been removed in the
	DaISy dataset.

	The data shows evidence of discrepancies.  One of the
	issues studied in the above paper is the size of these
	discrepancies - measured in this case in terms of the norm
	of the smallest perturbation required to account for the
	difference between the nominal model and the data.
	
	The steady state input (prior to the start of the experiment)
	is u = 6.0 Volts.


