/**********************************************************************************/
/***	Scenario C (U influences stroke for blacks and whites, only 		***/
/***	influences mortality for blacks)influences search loop program.		***/													
/***	This program searches for baseline mortality hazards for ages 45+	***/
/***	for exp=0 and exp=1 to line up with US Life Tables. 			***/  
/***	The search loop is needed because stroke kills people: pstrokedeath	***/
/***	die at stroke and g4 for effect of stroke history on death. 		***/
/***	The program also searches for baseline stroke hazards for whites,	***/
/***	since U causes stroke.							***/
/***	The program also searches for the effect of exposure on log hazard 	***/
/***	of death from birth because U influences mortality from birth for	***/
/***	blacks and stroke (which is higher in blacks) kills people.		***/
/**********************************************************************************/
set more off

timer clear 1

timer on 1

/*Step 1: Set parameters*/ 

*specify prevalence of exposure
global pexp = 0.5 	

*parameters for Sij 
//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
/*effect of race on mortality from birth to be identified from search loop
global g1_0to1 =	0.30
global g1_1to5 = 	0.37
global g1_5to10 = 	0.23
global g1_10to15 = 	0.56
global g1_15to20 = 	0.94
global g1_20to25 = 	0.92
global g1_25to30 = 	0.79
global g1_30to35 = 	0.78
global g1_35to40 = 	0.77
global g1_40to45 = 	0.74
global g1_45to50 = 	0.66
global g1_50to55 = 	0.59
global g1_55to60 = 	0.48
global g1_60to65 = 	0.32
global g1_65to70 = 	0.19
global g1_70to75 = 	0.10
global g1_75to80 = 	-0.05
global g1_80to85 = 	-0.10
global g1_85to90 = 	-0.19
global g1_90to95 = 	-0.24
global g1_95to100 =	-0.33*/

*effects of covariates on mortality risk
global g2 = 0		//log(HR) for effect of U on death 
global g3 = ln(1.5)	//log(HR) for interaction effect of exposure & U on death
global g4 = ln(2)	//log(HR) for effect of stroke history on death
global g5 = 0		//delete later

*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
global lambda_0to1 = 	0.0749
global lambda_1to5 = 	0.0082
global lambda_5to10 = 	0.0028
global lambda_10to15 = 	0.0021
global lambda_15to20 = 	0.0034
global lambda_20to25 = 	0.0048
global lambda_25to30 = 	0.0056
global lambda_30to35 = 	0.0062
global lambda_35to40 = 	0.0068
global lambda_40to45 = 	0.0078
/*baseline hazard of mortality for ages 45-100 to be identified from search loop
global lambda_45to50 = 	
global lambda_50to55 = 	
global lambda_55to60 =	
global lambda_60to65 = 	
global lambda_65to70 =  
global lambda_70to75 = 	
global lambda_75to80 = 	
global lambda_80to85 = 	
global lambda_85to90 = 	
global lambda_90to95 = 	
global lambda_95to100 =	*/

/*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011*/
/*baseline hazard of stroke for exp=0 whitesto be identified from search loop
global stk_lambda_exp0_45to50 = 0.005066581
global stk_lambda_exp0_50to55 = 0.009992424
global stk_lambda_exp0_55to60 = 0.01970135
global stk_lambda_exp0_60to65 = 0.032781652
global stk_lambda_exp0_65to70 = 0.046191677
global stk_lambda_exp0_70to75 = 0.065057384
global stk_lambda_exp0_75to80 = 0.09806763
global stk_lambda_exp0_80to85 = 0.117292282
global stk_lambda_exp0_85to90 = 0.128072578
global stk_lambda_exp0_90to95 = 0.142952161*/

/*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011*/
global stk_lambda_delta = 0.002 
global stk_lambda_exp1_45to50 = 	$stk_lambda_exp0_45to50 + $stk_lambda_delta
global stk_lambda_exp1_50to55 = 	$stk_lambda_exp0_50to55 + $stk_lambda_delta
global stk_lambda_exp1_55to60 =		$stk_lambda_exp0_55to60 + $stk_lambda_delta
global stk_lambda_exp1_60to65 = 	$stk_lambda_exp0_60to65 + $stk_lambda_delta
global stk_lambda_exp1_65to70 =  	$stk_lambda_exp0_65to70 + $stk_lambda_delta
global stk_lambda_exp1_70to75 = 	$stk_lambda_exp0_70to75 + $stk_lambda_delta
global stk_lambda_exp1_75to80 = 	$stk_lambda_exp0_75to80 + $stk_lambda_delta
global stk_lambda_exp1_80to85 = 	$stk_lambda_exp0_80to85 + $stk_lambda_delta
global stk_lambda_exp1_85to90 = 	$stk_lambda_exp0_85to90 + $stk_lambda_delta 
global stk_lambda_exp1_90to95 = 	$stk_lambda_exp0_90to95 + $stk_lambda_delta

*parameters for stroke risk
global b1 = ln(1.5)		//log (HR) for U on stroke
global b2 = 0			//delete later
global b3 = 0 			//delete later
global b4 = 0 			//delete later
global b5 = 0			//delete later

*probability of death at stroke
global pstrokedeath = 0.25 

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


/*Challenge: we want the same number of people having strokes after we add U.

This is a do file to run the data generating do file with different values for 
the baseline mortality odds start with a value for the target that is definitely 
too low for the target stroke or cumulative mortality within each age band. 
The search loop will stop at the first increment where the start value results 
in a stroke rate/cumulative mortality >=the target stroke rate/cumulative 
mortality.*/

/*Target and starting values for baseline mortality hazard for whites*/
*Target values
global target_p_death45to50_exp0 = 0.0469
global target_p_death50to55_exp0 = 0.0623
global target_p_death55to60_exp0 = 0.0890
global target_p_death60to65_exp0 = 0.1267
global target_p_death65to70_exp0 = 0.1846
global target_p_death70to75_exp0 = 0.2726
global target_p_death75to80_exp0 = 0.3840
global target_p_death80to85_exp0 = 0.5259 
global target_p_death85to90_exp0 = 0.6705 
global target_p_death90to95_exp0 = 0.7850
*global target_p_death95to100_exp0 = 0.9036

*Starting values
global lambda_45to50l = 0
global lambda_50to55l = 0
global lambda_55to60l = 0
global lambda_60to65l = 0
global lambda_65to70l = 0.02
global lambda_70to75l = 0.01
global lambda_75to80l = 0.04	
global lambda_80to85l = 0.05	
global lambda_85to90l = 0.09
global lambda_90to95l = 0.2
*global lambda_95to100l = 0.22

/*Target and starting values for baseline stroke hazard for whites*/
*Target values
global target_strokerate45to50_exp0 = 4.8
global target_strokerate50to55_exp0 = 11.3
global target_strokerate55to60_exp0 = 20.3
global target_strokerate60to65_exp0 = 32.9
global target_strokerate65to70_exp0 = 48.0
global target_strokerate70to75_exp0 = 68.8
global target_strokerate75to80_exp0 = 101.6
global target_strokerate80to85_exp0 = 127.3
global target_strokerate85to90_exp0 = 142.2
global target_strokerate90to95_exp0 = 166.6

*Starting values
global stk_lambda_exp0_45to50l = 0
global stk_lambda_exp0_50to55l = 0
global stk_lambda_exp0_55to60l = 0.001
global stk_lambda_exp0_60to65l = 0.001
global stk_lambda_exp0_65to70l = 0.002
global stk_lambda_exp0_70to75l = 0.004
global stk_lambda_exp0_75to80l = 0.007
global stk_lambda_exp0_80to85l = 0.003
global stk_lambda_exp0_85to90l = 0.01
global stk_lambda_exp0_90to95l = 0.01


/*Target and starting values for baseline mortality hazard for blacks*/
*Target values
global target_g1_0to1 = 0.30
global target_g1_1to5 = 0.37
global target_g1_5to10 = 0.23
global target_g1_10to15 = 0.56
global target_g1_15to20 = 0.94
global target_g1_20to25 = 0.92
global target_g1_25to30 = 0.79
global target_g1_30to35 = 0.78
global target_g1_35to40 = 0.77
global target_g1_40to45 = 0.74
global target_g1_45to50 = 0.66
global target_g1_50to55 = 0.59
global target_g1_55to60 = 0.48
global target_g1_60to65 = 0.32
global target_g1_65to70 = 0.19
global target_g1_70to75 = 0.10
global target_g1_75to80 = -0.05
global target_g1_80to85 = -0.10
global target_g1_85to90 = -0.19
global target_g1_90to95 = -0.24
*global target_g1_95to100 = -0.33

*Starting values
global g1_0to1l = 0.1
global g1_1to5l = 0.1 
global g1_5to10l = 0 
global g1_10to15l = 0.1 
global g1_15to20l = 0.4
global g1_20to25l = 0.4 
global g1_25to30l = 0.3
global g1_30to35l = 0.3
global g1_35to40l = 0.3 
global g1_40to45l = 0.3
global g1_45to50l = 0.1
global g1_50to55l = 0.1
global g1_55to60l = 0
global g1_60to65l = 0
global g1_65to70l = -0.05
global g1_70to75l = -0.1
global g1_75to80l = -0.2
global g1_80to85l = -0.3
global g1_85to90l = -0.5
global g1_90to95l = -0.6
*global g1_95to100l = -0.75


/***************************************************************************************************************/
/***************************************************************************************************************/
/***		Search loops for baseline mortality hazard 														 ***/
/***		and stroke hazard for whites ages 45-95											 				 ***/
/***************************************************************************************************************/
/***************************************************************************************************************/


/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 45-50	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death45to50_exp0 = $target_p_death45to50_exp0
   *add lower bound guess here
   local lambda_45to50l = $lambda_45to50l	
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_45to50lt =`lambda_45to50l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_45to50 = `lambda_45to50lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			/*local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			/*local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			/*local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50lt'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

		
		 qui sum death45to50 if (exposure==0)
		 local p_death45to50_exp0 = r(mean) + (`pstrokedeath'*0.0012) //Add in approximate number of stroke deaths
         if `p_death45to50_exp0' > `target_p_death45to50_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_45to50==`lambda_45to50lt', p_death45to50_exp0=`p_death45to50_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_45to50==`lambda_45to50lt', p_death45to50_exp0=`p_death45to50_exp0'
         }
      }
   }
   local lambda_45to50`i' = `lambda_45to50lt'
}
clear
set obs 5
gen lambda_45to50=.
forvalues i=1/5 {
   replace lambda_45to50=`lambda_45to50`i'' in `i'
}
sum lambda_45to50
global lambda_45to50 = `r(mean)'
global lambda_45to50_min = `r(min)'
global lambda_45to50_max = `r(max)'




/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 45-50	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate45to50_exp0 = $target_strokerate45to50_exp0
   *add lower bound guess here
   local stk_lambda_exp0_45to50l = $stk_lambda_exp0_45to50l
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_45to50lt =`stk_lambda_exp0_45to50l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_45to50 = `stk_lambda_exp0_45to50lt'"
         clear
		 /*Create blank data set*/
		set obs 20000 //creates blank dataset with XXXXX observations
		gen id = _n

		/*Step 1: Set parameters*/

		*specify prevalence of exposure
		local pexp = $pexp 	

		*parameters for Sij 
		//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
		local g1_0to1 =		0.30
		local g1_1to5 = 	0.37
		local g1_5to10 = 	0.23
		local g1_10to15 = 	0.56
		local g1_15to20 = 	0.94
		local g1_20to25 = 	0.92
		local g1_25to30 = 	0.79
		local g1_30to35 = 	0.78
		local g1_35to40 = 	0.77
		local g1_40to45 = 	0.74
		local g1_45to50 = 	0.66
		local g1_50to55 = 	0.59
		local g1_55to60 = 	0.48
		local g1_60to65 = 	0.32
		local g1_65to70 = 	0.19
		local g1_70to75 = 	0.10
		local g1_75to80 = 	-0.05
		local g1_80to85 = 	-0.10
		local g1_85to90 = 	-0.19
		local g1_90to95 = 	-0.24
		local g1_95to100 =	-0.33

		*effects of covariates on mortality risk
		local g2 = $g2 //effect of U on log hazard of death
		local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
		local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
		local g5 = $g5 //delete later	


		*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
		local lambda_0to1 = 	$lambda_0to1
		local lambda_1to5 = 	$lambda_1to5
		local lambda_5to10 = 	$lambda_5to10
		local lambda_10to15 = 	$lambda_10to15
		local lambda_15to20 = 	$lambda_15to20
		local lambda_20to25 = 	$lambda_20to25
		local lambda_25to30 = 	$lambda_25to30
		local lambda_30to35 = 	$lambda_30to35
		local lambda_35to40 = 	$lambda_35to40
		local lambda_40to45 = 	$lambda_40to45
		local lambda_45to50 = 	$lambda_45to50
		/*local lambda_50to55 = 	$lambda_50to55
		local lambda_55to60 =	$lambda_55to60
		local lambda_60to65 = 	$lambda_60to65
		local lambda_65to70 =  	$lambda_65to70
		local lambda_70to75 = 	$lambda_70to75 
		local lambda_75to80 = 	$lambda_75to80
		local lambda_80to85 = 	$lambda_80to85
		local lambda_85to90 = 	$lambda_85to90
		local lambda_90to95 = 	$lambda_90to95
		local lambda_95to100 =	$lambda_95to100*/

		*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
		/*local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
		local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
		local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
		local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
		local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
		local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
		local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
		local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
		local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
		local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

		*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
		local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
		/*local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
		local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
		local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
		local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
		local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
		local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
		local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
		local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
		local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

		*effects of covariates on stroke risk
		local b1 = $b1 			//log (HR) for U on stroke
		local b2 = $b2 			//delete later
		local b3 = $b3 			//delete later
		local b4 = $b3 			//delete later
		local b5 = $b4 			//delete later

		*probability of death at stroke
		local pstrokedeath = $pstrokedeath 


		/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
		and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
		gen exposure = runiform()<`pexp' 
		/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


		/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
		gen U = rnormal(0,1)


		/*Step 4: Generate survival time for each person and strokes for people alive
		at each interval. 
		a. Each person's underlying time to death is generated for each age interval, 
		conditional on the past provided the person has not died in a previous interval, 
		under an exponential survival distribtion. If the persons generated survival 
		time exceeds the length of the interval between study visits j and j+1, 
		she is considered alive at study visit j+1 and a new survival time is 
		generated for the next interval conditional on history up to the start of the 
		interval, and the process is repeated until the persons survival time falls 
		within a given interval or the end of the study, whichever comes first. Each 
		persons hazard function is defined as:
		h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
		A persons survival time for a given time interval at risk is generated using 
		the inverse cumulative hazard function transformation formula described by 
		Bender et al. (Stat Med 2011)
		b. Stroke code is adapted for survival time code.*/

		*ia. Generate uniform random variable for generating survival time
		gen U_0to1 = runiform()
		gen U_1to5 = runiform()
		gen U_5to10 = runiform()
		gen U_10to15 = runiform()
		gen U_15to20 = runiform()
		gen U_20to25 = runiform()
		gen U_25to30 = runiform()
		gen U_30to35 = runiform()
		gen U_35to40 = runiform()
		gen U_40to45 = runiform()
		gen U_45to50 = runiform()
		gen U_50to55 = runiform()
		gen U_55to60 = runiform()
		gen U_60to65 = runiform()
		gen U_65to70 = runiform()
		gen U_70to75 = runiform()
		gen U_75to80 = runiform()
		gen U_80to85 = runiform()
		gen U_85to90 = runiform()
		gen U_90to95 = runiform()
		gen U_95to100 = runiform()

		*ib. Generate uniform random variable for generating stroke time
		gen U2_45to50 = runiform()
		gen U2_50to55 = runiform()
		gen U2_55to60 = runiform()
		gen U2_60to65 = runiform()
		gen U2_65to70 = runiform()
		gen U2_70to75 = runiform()
		gen U2_75to80 = runiform()
		gen U2_80to85 = runiform()
		gen U2_85to90 = runiform()
		gen U2_90to95 = runiform()


		*ii. Generate survival time and stroke time for each interval

		gen survage = .
		gen strokeage = .

		/***Ages 0-45: no strokes, so only need to generate survival***/
		/*Interval 0-1*/
		*Generate survival time from time 0
		gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0))
		*Generate death indicator for interval 0-1 
		gen death0to1 = 0
		replace death0to1 = 1 if (survtime0to1 < 1) 
		replace survage = survtime0to1 if (death0to1==1)
		*Generate indicator for death before age 1
		gen death1 = death0to1


		/*Interval 1-5*/
		*Generate survival time from time 1
		gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death0to1==0) 
		*Generate death indicator for interval 1-5 
		gen death1to5 = 0 if (death1==0)
		replace death1to5 = 1 if (survtime1to5 < 5)
		replace survage = 1 + survtime1to5 if (death1to5==1)
		*Generate indicator for death before age 5
		gen death5 = 0
		replace death5 = 1 if survage < 5


		/*Interval 5-10*/
		*Generate survival time from time 2
		gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death5==0)
		*Generate death indicator for interval 5-10 
		gen death5to10 = 0 if (death5==0)
		replace death5to10 = 1 if (survtime5to10 < 5) 
		replace survage = 5 + survtime5to10 if (death5to10==1)
		*Generate indicator for death before age 10
		gen death10 = 0
		replace death10 = 1 if survage < 10


		/*Interval 10-15*/
		*Generate survival time from time 3
		gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death10==0)
		*Generate death indicator for interval 10-15
		gen death10to15 = 0 if (death10==0)
		replace death10to15 = 1 if (survtime10to15 < 5) 
		replace survage = 10 + survtime10to15 if (death10to15==1)
		*Generate indicator for death before age 15
		gen death15 = 0
		replace death15 = 1 if survage < 15

		/*Interval 15-20*/
		*Generate survival time from time 4
		gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death15==0)
		*Generate death indicator for interval 15-20 
		gen death15to20 = 0 if (death15==0)
		replace death15to20 = 1 if (survtime15to20 < 5) 
		replace survage = 15 + survtime15to20 if (death15to20==1)
		*Generate indicator for death before age 20
		gen death20 = 0
		replace death20 = 1 if survage < 20


		/*Interval 20-25*/
		*Generate survival time from time 5
		gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death20==0) 
		*Generate death indicator for interval 20-25 
		gen death20to25 = 0 if (death20==0)
		replace death20to25 = 1 if (survtime20to25 < 5) 
		replace survage = 20 + survtime20to25 if (death20to25==1)
		*Generate indicator for death before age 25
		gen death25 = 0
		replace death25 = 1 if survage < 25


		/*Interval 25-30*/
		*Generate survival time from time 6
		gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death25==0)
		*Generate death indicator for interval 25-30 
		gen death25to30 = 0 if (death25==0)
		replace death25to30 = 1 if (survtime25to30 < 5) 
		replace survage = 25 + survtime25to30 if (death25to30==1)
		*Generate indicator for death before age 30
		gen death30 = 0
		replace death30 = 1 if survage < 30


		/*Interval 30-35*/
		*Generate survival time from time 7
		gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death30==0)
		*Generate death indicator for interval 30-35 
		gen death30to35 = 0 if (death30==0)
		replace death30to35 = 1 if (survtime30to35 < 5) 
		replace survage = 30 + survtime30to35 if (death30to35==1)
		*Generate indicator for death before age 35
		gen death35 = 0
		replace death35 = 1 if survage < 35


		/*Interval 35-40*/
		*Generate survival time from time 8
		gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death35==0)
		*Generate death indicator for interval 35-40
		gen death35to40 = 0 if (death35==0)
		replace death35to40 = 1 if (survtime35to40 < 5) 
		replace survage = 35 + survtime35to40 if (death35to40==1)
		*Generate indicator for death before age 40
		gen death40 = 0
		replace death40 = 1 if survage < 40


		/*Interval 40-45*/
		*Generate survival time from time 9
		gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death40==0)
		*Generate death indicator for interval 40-45
		gen death40to45 = 0 if (death40==0)
		replace death40to45 = 1 if (survtime40to45 < 5) 
		replace survage = 40 + survtime40to45 if (death40to45==1)
		*Generate indicator for death before age 45
		gen death45 = 0
		replace death45 = 1 if survage < 45


		/***Starting at age 45--people are at risk of stroke, and prevalent stroke
		increases mortality risk, so we need to start iteratively generating survival 
		times and strokes***/


		/*Interval 45-50*/
		***a. Survival
		*Generate survival time from time 10
		gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
						if (death45==0)
		*Generate death indicator for interval 45-50
		gen death45to50 = 0 if (death45==0)
		replace death45to50 = 1 if (survtime45to50 < 5) 
		replace survage = 45 + survtime45to50 if (death45to50==1)
		*Generate indicator for death before age 50
		gen death50 = 0
		replace death50 = 1 if survage < 50

		***b. Stroke
		*Generate stroke time from time 10 exp==0
		gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50lt'*exp(`b1'*U)) ///
						if (exp==0 & death45==0)
		*Generate stroke time from time 10 exp==1
		replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
						if (exp==1 & death45==0)
		*Generate stroke indicator for interval 45-50
		gen stroke45to50 = 0 if (death45==0)
		replace stroke45to50 = 1 if (stroketime45to50 < 5) 
		replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
		*Generate prevalent stroke variable (stroke up to age 50)
		gen stroke_history = 0
		replace stroke_history = 1 if (stroke45to50==1)
		replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
		*Generate indicator for stroke before age 50
		gen stroke50 = 0
		replace stroke50 = 1 if strokeage < 50

		*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
		replace death45to50 = 1 if (strokedeath45to50==1) 
		replace survage = strokeage if (strokedeath45to50==1) 
		*Generate indicator for stroke death before age 50
		gen strokedeath50 = 0
		replace strokedeath50 = 1 if (strokedeath45to50==1)
		

		*Step added for search loop: generate survage = 49.999 for people who are still alive
		replace survage = 49.999 if (death50==0 & survage==.)
							
		*iv. Generate variables for stroke between ages 45 to 95
		gen stroke = 0
		replace stroke = 1 if (stroke45to50==1/* | stroke50to55==1 | stroke55to60==1 | ///
		stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80==1 | ///
		stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


		*v. Generate strokeage for people who didn't develop a stroke
		replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

		
		*vi. Generate age-stratified stroke variables
		*ages 45 to 50
		gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
		gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
		replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)


		qui stset strokeage, failure(stroke45to50==1) id(id) enter(strokeage45to50_start) exit(strokeage45to50_end)
		qui stptime if (exposure==0), title(person-years) per(10000)
		local strokerate1000pys45to50_exp0 = r(rate)
        if `r(rate)' > `target_strokerate45to50_exp0'+.001 {
           local toolow=0
           noisily di in red "I stopped at stk_lambda_exp0_45to50==`stk_lambda_exp0_45to50lt', strokerate1000pys45to50_exp0=`r(rate)'
        }
        else {
           noisily di in white "I did NOT stop at stk_lambda_exp0_45to50==`stk_lambda_exp0_45to50lt', strokerate1000pys45to50_exp0=`r(rate)'
        }
      }
   }
   local stk_lambda_exp0_45to50`i' = `stk_lambda_exp0_45to50lt'
}
clear
set obs 5
gen stk_lambda_exp0_45to50=.
forvalues i=1/5 {
   replace stk_lambda_exp0_45to50=`stk_lambda_exp0_45to50`i'' in `i'
}
sum stk_lambda_exp0_45to50
global stk_lambda_exp0_45to50 = `r(mean)'
global stk_lambda_exp0_45to50_min = `r(min)'
global stk_lambda_exp0_45to50_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 50-55	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death50to55_exp0 = $target_p_death50to55_exp0
   *add lower bound guess here
   local lambda_50to55l = $lambda_50to55l 
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_50to55lt =`lambda_50to55l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_50to55 = `lambda_50to55lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			/*local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			/*local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			/*local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)
			
			
			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55lt'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

				
		 qui sum death50to55 if (exposure==0)
		 local p_death50to55_exp0 = r(mean) + (`pstrokedeath'*0.0023) //Add in approximate number of stroke deaths 
         if `p_death50to55_exp0' > `target_p_death50to55_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_50to55==`lambda_50to55lt', p_death50to55_exp0=`p_death50to55_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_50to55==`lambda_50to55lt', p_death50to55_exp0=`p_death50to55_exp0'
         }
      }
   }
   local lambda_50to55`i' = `lambda_50to55lt'
}
clear
set obs 5
gen lambda_50to55=.
forvalues i=1/5 {
   replace lambda_50to55=`lambda_50to55`i'' in `i'
}
sum lambda_50to55
global lambda_50to55 = `r(mean)'
global lambda_50to55_min = `r(min)'
global lambda_50to55_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 50-55	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate50to55_exp0 = $target_strokerate50to55_exp0
   *add lower bound guess here
   local stk_lambda_exp0_50to55l = $stk_lambda_exp0_50to55l
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_50to55lt =`stk_lambda_exp0_50to55l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_50to55 = `stk_lambda_exp0_50to55lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			/*local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			/*local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55lt' + 0.002
			/*local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)
			
			
			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55lt'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)

			
			*Step added for search loop: replace survage = 54.999 for people who are still alive 
			replace survage = 54.999 if (death55==0 & survage==.)
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1/* | stroke55to60==1 | ///
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80==1 | ///
			stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)
			
		 qui stset strokeage, failure(stroke50to55==1) id(id) enter(strokeage50to55_start) exit(strokeage50to55_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys50to55_exp0 = r(rate)
         if `r(rate)' > `target_strokerate50to55_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_50to55==`stk_lambda_exp0_50to55lt', strokerate1000pys50to55_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_50to55==`stk_lambda_exp0_50to55lt', strokerate1000pys50to55_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_50to55`i' = `stk_lambda_exp0_50to55lt'
}
clear
set obs 5
gen stk_lambda_exp0_50to55=.
forvalues i=1/5 {
   replace stk_lambda_exp0_50to55=`stk_lambda_exp0_50to55`i'' in `i'
}
sum stk_lambda_exp0_50to55
global stk_lambda_exp0_50to55 = `r(mean)'
global stk_lambda_exp0_50to55_min = `r(min)'
global stk_lambda_exp0_50to55_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 55-60	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death55to60_exp0 = $target_p_death55to60_exp0 
   *add lower bound guess here
   local lambda_55to60l = $lambda_55to60l 
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_55to60lt =`lambda_55to60l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_55to60 = `lambda_55to60lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			/*local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			/*local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			/*local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			**effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)
			
			
			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55
					
			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)

			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60lt'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60
			

		 qui sum death55to60 if (exposure==0)
		 local p_death55to60_exp0 = r(mean) + (`pstrokedeath'*0.0042) //Add in approximate number of stroke deaths 
         if `p_death55to60_exp0' > `target_p_death55to60_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_55to60==`lambda_55to60lt', p_death55to60_exp0=`p_death55to60_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_55to60==`lambda_55to60lt', p_death55to60_exp0=`p_death55to60_exp0'
         }
      }
   }
   local lambda_55to60`i' = `lambda_55to60lt'
}
clear
set obs 5
gen lambda_55to60=.
forvalues i=1/5 {
   replace lambda_55to60=`lambda_55to60`i'' in `i'
}
sum lambda_55to60
global lambda_55to60 = `r(mean)'
global lambda_55to60_min = `r(min)'
global lambda_55to60_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 55-60	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate55to60_exp0 = $target_strokerate55to60_exp0
   *add lower bound guess here
   local stk_lambda_exp0_55to60l = $stk_lambda_exp0_55to60l
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_55to60lt =`stk_lambda_exp0_55to60l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_55to60 = `stk_lambda_exp0_55to60lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			/*local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			/*local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60lt' + 0.002
			/*local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for L on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)
			
			
			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55
						
			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)

			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60lt'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Step added for search loop: replace survage = 59.999 for people who are still alive*/
			replace survage = 59.999 if (death60==0 & survage==.)
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1/* | /// 
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80==1 | ///
			stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
			
			
		 qui stset strokeage, failure(stroke55to60==1) id(id) enter(strokeage55to60_start) exit(strokeage55to60_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys55to60_exp0 = r(rate)
         if `r(rate)' > `target_strokerate55to60_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_55to60==`stk_lambda_exp0_55to60lt', strokerate1000pys55to60_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_55to60==`stk_lambda_exp0_55to60lt', strokerate1000pys55to60_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_55to60`i' = `stk_lambda_exp0_55to60lt'
}
clear
set obs 5
gen stk_lambda_exp0_55to60=.
forvalues i=1/5 {
   replace stk_lambda_exp0_55to60=`stk_lambda_exp0_55to60`i'' in `i'
}
sum stk_lambda_exp0_55to60
global stk_lambda_exp0_55to60 = `r(mean)'
global stk_lambda_exp0_55to60_min = `r(min)'
global stk_lambda_exp0_55to60_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 60-65		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death60to65_exp0 = $target_p_death60to65_exp0
   *add lower bound guess here
   local lambda_60to65l = $lambda_60to65l 
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_60to65lt =`lambda_60to65l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_60to65 = `lambda_60to65lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			/*local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			/*local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			/*local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			**effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65lt'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			
		 qui sum death60to65 if (exposure==0) 
		 local p_death60to65_exp0 = r(mean) + (`pstrokedeath'*0.0062) //Add in approximate number of stroke deaths 
         if `p_death60to65_exp0' > `target_p_death60to65_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_60to65==`lambda_60to65lt', p_death60to65_exp0=`p_death60to65_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_60to65==`lambda_60to65lt', p_death60to65_exp0=`p_death60to65_exp0'
         }
      }
   }
   local lambda_60to65`i' = `lambda_60to65lt'
}
clear
set obs 5
gen lambda_60to65=.
forvalues i=1/5 {
   replace lambda_60to65=`lambda_60to65`i'' in `i'
}
sum lambda_60to65
global lambda_60to65 = `r(mean)'
global lambda_60to65_min = `r(min)'
global lambda_60to65_max = `r(max)'



/******************************************************************************/
/***		Find baseline STROKE hazard for whites age 60-65	    ***/
/******************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate60to65_exp0 = $target_strokerate60to65_exp0
   *add lower bound guess here
   local stk_lambda_exp0_60to65l = $stk_lambda_exp0_60to65l 
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_60to65lt =`stk_lambda_exp0_60to65l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green " stk_lambda_exp0_60to65 = `stk_lambda_exp0_60to65lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			/*local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			/*local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
		
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65lt' + 0.002
			/*local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65lt'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)


			/*Step added for search loop: replace survage = 64.999 for people who are still alive*/ 
			replace survage = 64.999 if (death65==0 & survage==.)
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1/* | stroke65to70==1 | stroke70to75==1 | stroke75to80==1 | ///
			stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)
			
			
		 qui stset strokeage, failure(stroke60to65==1) id(id) enter(strokeage60to65_start) exit(strokeage60to65_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys60to65_exp0 = r(rate)
         if `r(rate)' > `target_strokerate60to65_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_60to65l==`stk_lambda_exp0_60to65lt', strokerate1000pys60to65_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_60to65==`stk_lambda_exp0_60to65lt', strokerate1000pys60to65_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_60to65`i' = `stk_lambda_exp0_60to65lt'
}
clear
set obs 5
gen stk_lambda_exp0_60to65=.
forvalues i=1/5 {
   replace stk_lambda_exp0_60to65=`stk_lambda_exp0_60to65`i'' in `i'
}
sum stk_lambda_exp0_60to65
global stk_lambda_exp0_60to65 = `r(mean)'	
global stk_lambda_exp0_60to65_min = `r(min)'	
global stk_lambda_exp0_60to65_max = `r(max)'								



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 65-70	  	***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death65to70_exp0 = $target_p_death65to70_exp0
   *add lower bound guess here
   local lambda_65to70l = $lambda_65to70l
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_65to70lt =`lambda_65to70l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_65to70 = `lambda_65to70lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			/*local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			/*local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			/*local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70lt'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70


		 qui sum death65to70 if (exposure==0)
		 local p_death65to70_exp0 = r(mean) + (`pstrokedeath'*0.0073) //Add in approximate number of stroke deaths 
         if `p_death65to70_exp0' > `target_p_death65to70_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_65to70==`lambda_65to70lt', p_death65to70_exp0=`p_death65to70_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_65to70==`lambda_65to70lt', p_death65to70_exp0=`p_death65to70_exp0'
         }
      }
   }
   local lambda_65to70`i' = `lambda_65to70lt'
}
clear
set obs 5
gen lambda_65to70=.
forvalues i=1/5 {
   replace lambda_65to70=`lambda_65to70`i'' in `i'
}
sum lambda_65to70
global lambda_65to70 = `r(mean)'
global lambda_65to70_min = `r(min)'
global lambda_65to70_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 65-70	  	***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate65to70_exp0 = $target_strokerate65to70_exp0 
   *add lower bound guess here
   local stk_lambda_exp0_65to70l = $stk_lambda_exp0_65to70l 
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_65to70lt =`stk_lambda_exp0_65to70l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_65to70 = `stk_lambda_exp0_65to70lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33 

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	


			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			/*local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			/*local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70lt' + 0.002
			/*local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

				/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70lt'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)

				
			/*Step added for search loop: replace survage = 69.999 for people who are still alive*/ 
			replace survage = 69.999 if (death70==0 & survage==.)
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1 | stroke65to70==1/* | stroke70to75==1 | stroke75to80==1 | ///
			stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)	
			
			*ages 65 to 70
			gen strokeage65to70_start = 65 if (strokeage !=. & strokeage>=65)
			gen strokeage65to70_end = 69.99999 if (strokeage !=. & strokeage>=70)
			replace strokeage65to70_end = strokeage if (strokeage !=. & strokeage>=65 & strokeage<70)
			
			
		 qui stset strokeage, failure(stroke65to70==1) id(id) enter(strokeage65to70_start) exit(strokeage65to70_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys65to70_exp0 = r(rate)
         if `r(rate)' > `target_strokerate65to70_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_65to70==`stk_lambda_exp0_65to70lt', strokerate1000pys65to70_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_65to70==`stk_lambda_exp0_65to70lt', strokerate1000pys65to70_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_65to70`i' = `stk_lambda_exp0_65to70lt'
}
clear
set obs 5
gen stk_lambda_exp0_65to70=.
forvalues i=1/5 {
   replace stk_lambda_exp0_65to70=`stk_lambda_exp0_65to70`i'' in `i'
}
sum stk_lambda_exp0_65to70
global stk_lambda_exp0_65to70 = `r(mean)'
global stk_lambda_exp0_65to70_min = `r(min)'
global stk_lambda_exp0_65to70_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 70-75	  	***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death70to75_exp0 = $target_p_death70to75_exp0
   *add lower bound guess here
   local lambda_70to75l = $lambda_70to75l
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_70to75lt =`lambda_70to75l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_70to75 = `lambda_70to75lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			**effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & L on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			/*local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			/*local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			/*local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				

			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75lt'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75


		 qui sum death70to75 if (exposure==0)
		 local p_death70to75_exp0 = r(mean) + (`pstrokedeath'*0.0082) //Add in approximate number of stroke deaths 
         if `p_death70to75_exp0' > `target_p_death70to75_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_70to75==`lambda_70to75lt', p_death70to75_exp0=`p_death70to75_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_70to75==`lambda_70to75lt', p_death70to75_exp0=`p_death70to75_exp0'
         }
      }
   }
   local lambda_70to75`i' = `lambda_70to75lt'
}
clear
set obs 5
gen lambda_70to75=.
forvalues i=1/5 {
   replace lambda_70to75=`lambda_70to75`i'' in `i'
}
sum lambda_70to75
global lambda_70to75 = `r(mean)'
global lambda_70to75_min = `r(min)'
global lambda_70to75_max = `r(max)'



/*******************************************************************************/
/***		Find baseline STROKE hazard for whites age 70-75	     ***/ 
/*******************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate70to75_exp0 = $target_strokerate70to75_exp0
   *add lower bound guess here
   local stk_lambda_exp0_70to75l = $stk_lambda_exp0_70to75l
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_70to75lt =`stk_lambda_exp0_70to75l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_70to75 = `stk_lambda_exp0_70to75'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			/*local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			/*local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75lt' + 0.002
			/*local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
			
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75lt'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
	
				
			/*Step added for search loop: replace survage = 74.999 for people who are still alive*/ 
			replace survage = 74.999 if (death75==0 & survage==.)
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1/* | stroke75to80==1 | ///
			stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)	
			
			*ages 65 to 70
			gen strokeage65to70_start = 65 if (strokeage !=. & strokeage>=65)
			gen strokeage65to70_end = 69.99999 if (strokeage !=. & strokeage>=70)
			replace strokeage65to70_end = strokeage if (strokeage !=. & strokeage>=65 & strokeage<70)
			
			*ages 70 to 75
			gen strokeage70to75_start = 70 if (strokeage !=. & strokeage>=70)
			gen strokeage70to75_end = 74.99999 if (strokeage !=. & strokeage>=75)
			replace strokeage70to75_end = strokeage if (strokeage !=. & strokeage>=70 & strokeage<75)
			
			
		 qui stset strokeage, failure(stroke70to75==1) id(id) enter(strokeage70to75_start) exit(strokeage70to75_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys70to75_exp0 = r(rate)
         if `r(rate)' > `target_strokerate70to75_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_70to75==`stk_lambda_exp0_70to75lt', strokerate1000pys70to75_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_70to75==`stk_lambda_exp0_70to75lt', strokerate1000pys70to75_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_70to75`i' = `stk_lambda_exp0_70to75lt'
}
clear
set obs 5
gen stk_lambda_exp0_70to75=.
forvalues i=1/5 {
   replace stk_lambda_exp0_70to75=`stk_lambda_exp0_70to75`i'' in `i'
}
sum stk_lambda_exp0_70to75
global stk_lambda_exp0_70to75 = `r(mean)'
global stk_lambda_exp0_70to75_min = `r(min)'
global stk_lambda_exp0_70to75_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 75-80	  	***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death75to80_exp0 = $target_p_death75to80_exp0
   *add lower bound guess here
   local lambda_75to80l = $lambda_75to80l 
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_75to80lt =`lambda_75to80l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_75to80 = `lambda_75to80lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			/*local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			/*local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			/*local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				

			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75
			
			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				

			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80lt'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80


		 qui sum death75to80 if (exposure==0)																					
		 local p_death75to80_exp0 = r(mean) + (`pstrokedeath'*0.0084) //Add in approximate number of stroke deaths 
         if `r(mean)' > `target_p_death75to80_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_75to80==`lambda_75to80lt', p_death75to80_exp0=`r(mean)'
         }
         else {
            noisily di in white "I did NOT stop at lambda_75to80==`lambda_75to80lt', p_death75to80_exp0=`r(mean)'
         }
      }
   }
   local lambda_75to80`i' = `lambda_75to80lt'
}
clear
set obs 5
gen lambda_75to80=.
forvalues i=1/5 {
   replace lambda_75to80=`lambda_75to80`i'' in `i'
}
sum lambda_75to80
global lambda_75to80 = `r(mean)'
global lambda_75to80_min = `r(min)'
global lambda_75to80_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 75-80		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate75to80_exp0 = $target_strokerate75to80_exp0 
   *add lower bound guess here
   local stk_lambda_exp0_75to80l = $stk_lambda_exp0_75to80l
   quietly forvalues x = 0(.0001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_75to80lt =`stk_lambda_exp0_75to80l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_75to80 = `stk_lambda_exp0_75to80lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			/*local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			/*local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80lt' + 0.002
			/*local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
			
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
	
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80lt'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
				
	
			/*Step added for search loop: replace survage = 79.999 for people who are still alive*/ 
			replace survage = 79.999 if (death80==0 & survage==.) 									
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80==1/* | ///
			stroke80to85==1 | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)	
			
			*ages 65 to 70
			gen strokeage65to70_start = 65 if (strokeage !=. & strokeage>=65)
			gen strokeage65to70_end = 69.99999 if (strokeage !=. & strokeage>=70)
			replace strokeage65to70_end = strokeage if (strokeage !=. & strokeage>=65 & strokeage<70)
			
			*ages 70 to 75
			gen strokeage70to75_start = 70 if (strokeage !=. & strokeage>=70)
			gen strokeage70to75_end = 74.99999 if (strokeage !=. & strokeage>=75)
			replace strokeage70to75_end = strokeage if (strokeage !=. & strokeage>=70 & strokeage<75)
			
			*ages 75 to 80
			gen strokeage75to80_start = 75 if (strokeage !=. & strokeage>=75)
			gen strokeage75to80_end = 79.99999 if (strokeage !=. & strokeage>=80)
			replace strokeage75to80_end = strokeage if (strokeage !=. & strokeage>=75 & strokeage<80)
			
			
		 qui stset strokeage, failure(stroke75to80==1) id(id) enter(strokeage75to80_start) exit(strokeage75to80_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys75to80_exp0 = r(rate)
         if `r(rate)' > `target_strokerate75to80_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_75to80==`stk_lambda_exp0_75to80lt', strokerate1000pys75to80_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_75to80==`stk_lambda_exp0_75to80lt', strokerate1000pys75to80_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_75to80`i' = `stk_lambda_exp0_75to80lt'
}
clear
set obs 5
gen stk_lambda_exp0_75to80=.
forvalues i=1/5 {
   replace stk_lambda_exp0_75to80=`stk_lambda_exp0_75to80`i'' in `i'
}
sum stk_lambda_exp0_75to80
global stk_lambda_exp0_75to80 = `r(mean)'
global stk_lambda_exp0_75to80_min = `r(min)'
global stk_lambda_exp0_75to80_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 80-85	  	***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death80to85_exp0 = $target_p_death80to85_exp0  
   *add lower bound guess here
   local lambda_80to85l = $lambda_80to85l
   quietly forvalues x = 0(.01)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_80to85lt =`lambda_80to85l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_80to85 = `lambda_80to85lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			/*local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			/*local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			/*local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				

			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75
			
			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				

			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80
			
		***b. Stroke
		*Generate stroke time from time 16 exp==0
		gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
						if (exp==0 & death75==0 & stroke75==0)
		*Generate stroke time from time 16 exp==1
		replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
						if (exp==1 & death75==0 & stroke75==0)
		*Generate stroke indicator for interval 75-80
		gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
		replace stroke75to80 = 1 if (stroketime75to80 < 5)
		replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
		*Update prevalent stroke variable (stroke up to age 80)
		replace stroke_history = 1 if (stroke75to80==1)
		replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
		*Generate indicator for stroke before age 80						
		gen stroke80 = 0
		replace stroke80 = 1 if strokeage < 80

		*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
		replace death75to80 = 1 if (strokedeath75to80==1) 
		replace survage = strokeage if (strokedeath75to80==1) 
		*Generate indicator for stroke death before age 80	
		gen strokedeath80 = 0
		replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			

		/*Interval 80-85*/
		***a. Survival
		*Generate survival time from time 17
		gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85lt'*exp(`g1_80to85'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death80==0)
		*Generate death indicator for interval 80-85
		gen death80to85 = 0 if (death80==0)
		replace death80to85 = 1 if (survtime80to85 < 5)
		replace survage = 80 + survtime80to85 if (death80to85==1)
		*Generate indicator for death before age 85
		gen death85 = 0
		replace death85 = 1 if survage < 85


		 qui sum death80to85 if (exposure==0)																					
		 local p_death80to85_exp0 = r(mean) + (`pstrokedeath'*0.0071) //Add in approximate number of stroke deaths
         if `p_death80to85_exp0' > `target_p_death80to85_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_80to85==`lambda_80to85lt', p_death80to85_exp0=`p_death80to85_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_80to85==`lambda_80to85lt', p_death80to85_exp0=`p_death80to85_exp0'
         }
      }
   }
   local lambda_80to85`i' = `lambda_80to85lt'
}
clear
set obs 5
gen lambda_80to85=.
forvalues i=1/5 {
   replace lambda_80to85=`lambda_80to85`i'' in `i'
}
sum lambda_80to85
global lambda_80to85 = `r(mean)'
global lambda_80to85_min = `r(min)'
global lambda_80to85_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 80-85	  	***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate80to85_exp0 = $target_strokerate80to85_exp0 
   *add lower bound guess here
   local stk_lambda_exp0_80to85l = $stk_lambda_exp0_80to85l
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_80to85lt =`stk_lambda_exp0_80to85l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_80to85 = `stk_lambda_exp0_80to85lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			/*local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			/*local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			/*local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
			
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
	
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
				
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85lt'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)	
							
	
			/*Step added for search loop: replace survage = 84.999 for people who are still alive*/ 
			replace survage = 84.999 if (death85==0 & survage==.) 									
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80== 1 | ///
			stroke80to85==1/* | stroke85to90==1 | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)	
			
			*ages 65 to 70
			gen strokeage65to70_start = 65 if (strokeage !=. & strokeage>=65)
			gen strokeage65to70_end = 69.99999 if (strokeage !=. & strokeage>=70)
			replace strokeage65to70_end = strokeage if (strokeage !=. & strokeage>=65 & strokeage<70)
			
			*ages 70 to 75
			gen strokeage70to75_start = 70 if (strokeage !=. & strokeage>=70)
			gen strokeage70to75_end = 74.99999 if (strokeage !=. & strokeage>=75)
			replace strokeage70to75_end = strokeage if (strokeage !=. & strokeage>=70 & strokeage<75)
			
			*ages 75 to 80
			gen strokeage75to80_start = 75 if (strokeage !=. & strokeage>=75)
			gen strokeage75to80_end = 79.99999 if (strokeage !=. & strokeage>=80)
			replace strokeage75to80_end = strokeage if (strokeage !=. & strokeage>=75 & strokeage<80)
			
			*ages 80 to 85
			gen strokeage80to85_start = 80 if (strokeage !=. & strokeage>=80)
			gen strokeage80to85_end = 84.99999 if (strokeage !=. & strokeage>=85)
			replace strokeage80to85_end = strokeage if (strokeage !=. & strokeage>=80 & strokeage<85)
			
		 qui stset strokeage, failure(stroke80to85==1) id(id) enter(strokeage80to85_start) exit(strokeage80to85_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys80to85_exp0 = r(rate)
         if `r(rate)' > `target_strokerate80to85_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_80to85==`stk_lambda_exp0_80to85lt', strokerate1000pys80to85_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_80to85==`stk_lambda_exp0_80to85lt', strokerate1000pys80to85_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_80to85`i' = `stk_lambda_exp0_80to85lt'
}
clear
set obs 5
gen stk_lambda_exp0_80to85=.
forvalues i=1/5 {
   replace stk_lambda_exp0_80to85=`stk_lambda_exp0_80to85`i'' in `i'
}
sum stk_lambda_exp0_80to85
global stk_lambda_exp0_80to85 = `r(mean)'
global stk_lambda_exp0_80to85_min = `r(min)'
global stk_lambda_exp0_80to85_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 85-90	  	***/	
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death85to90_exp0 = $target_p_death85to90_exp0  
   *add lower bound guess here
   local lambda_85to90l = $lambda_85to90l 
   quietly forvalues x = 0(.01)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_85to90lt =`lambda_85to90l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_85to90 = `lambda_85to90lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33
			
			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			/*local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			/*local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			/*local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				

			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75
			
			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				

			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80
			
		***b. Stroke
		*Generate stroke time from time 16 exp==0
		gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
						if (exp==0 & death75==0 & stroke75==0)
		*Generate stroke time from time 16 exp==1
		replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
						if (exp==1 & death75==0 & stroke75==0)
		*Generate stroke indicator for interval 75-80
		gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
		replace stroke75to80 = 1 if (stroketime75to80 < 5)
		replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
		*Update prevalent stroke variable (stroke up to age 80)
		replace stroke_history = 1 if (stroke75to80==1)
		replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
		*Generate indicator for stroke before age 80						
		gen stroke80 = 0
		replace stroke80 = 1 if strokeage < 80

		*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
		replace death75to80 = 1 if (strokedeath75to80==1) 
		replace survage = strokeage if (strokedeath75to80==1) 
		*Generate indicator for stroke death before age 80	
		gen strokedeath80 = 0
		replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			

		/*Interval 80-85*/
		***a. Survival
		*Generate survival time from time 17
		gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death80==0)
		*Generate death indicator for interval 80-85
		gen death80to85 = 0 if (death80==0)
		replace death80to85 = 1 if (survtime80to85 < 5)
		replace survage = 80 + survtime80to85 if (death80to85==1)
		*Generate indicator for death before age 85
		gen death85 = 0
		replace death85 = 1 if survage < 85

		***b. Stroke
		*Generate stroke time from time 17 exp==0
		gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
						if (exp==0 & death80==0 & stroke80==0)
		*Generate stroke time from time 17 exp==1
		replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
						if (exp==1 & death80==0 & stroke80==0)
		*Generate stroke indicator for interval 80-85
		gen stroke80to85 = 0 if (death80==0 & stroke80==0)
		replace stroke80to85 = 1 if (stroketime80to85 < 5)
		replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
		*Update prevalent stroke variable (stroke up to age 85)
		replace stroke_history = 1 if (stroke80to85==1)
		replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
		*Generate indicator for stroke before age 85						
		gen stroke85 = 0
		replace stroke85 = 1 if strokeage < 85

		*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
		replace death80to85 = 1 if (strokedeath80to85==1) 
		replace survage = strokeage if (strokedeath80to85==1) 
		*Generate indicator for stroke death before age 85	
		gen strokedeath85 = 0
		replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
			strokedeath80to85==1)
			

		/*Interval 85-90*/
		***a. Survival
		*Generate survival time from time 18
		gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90lt'*exp(`g1_85to90'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death85==0)
		*Generate death indicator for interval 85-90
		gen death85to90 = 0 if (death85==0)
		replace death85to90 = 1 if (survtime85to90 < 5)
		replace survage = 85 + survtime85to90 if (death85to90==1)
		*Generate indicator for death before age 90
		gen death90 = 0
		replace death90 = 1 if survage < 90


		 qui sum death85to90 if (exposure==0)																					
		 local p_death85to90_exp0 = r(mean) + (`pstrokedeath'*0.005) //Add in approximate number of stroke deaths
         if `p_death85to90_exp0' > `target_p_death85to90_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_85to90==`lambda_85to90lt', p_death85to90_exp0=`p_death85to90_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_85to90==`lambda_85to90lt', p_death85to90_exp0=`p_death85to90_exp0'
         }
      }
   }
   local lambda_85to90`i' = `lambda_85to90lt'
}
clear
set obs 5
gen lambda_85to90=.
forvalues i=1/5 {
   replace lambda_85to90=`lambda_85to90`i'' in `i'
}
sum lambda_85to90
global lambda_85to90 = `r(mean)'
global lambda_85to90_min = `r(min)'
global lambda_85to90_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 85-90	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate85to90_exp0 = $target_strokerate85to90_exp0
   *add lower bound guess here
   local stk_lambda_exp0_85to90l = $stk_lambda_exp0_85to90l
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_85to90lt =`stk_lambda_exp0_85to90l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_85to90 = `stk_lambda_exp0_85to90lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			/*local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			/*local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90lt' + 0.002
			/*local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
			
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
	
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
				
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)	
					
					
			/*Interval 85-90*/
			***a. Survival
			*Generate survival time from time 18
			gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death85==0)
			*Generate death indicator for interval 85-90
			gen death85to90 = 0 if (death85==0)
			replace death85to90 = 1 if (survtime85to90 < 5)
			replace survage = 85 + survtime85to90 if (death85to90==1)
			*Generate indicator for death before age 90
			gen death90 = 0
			replace death90 = 1 if survage < 90

			***b. Stroke
			*Generate stroke time from time 18 exp==0
			gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90lt'*exp(`b1'*U)) ///
							if (exp==0 & death85==0 & stroke85==0)
			*Generate stroke time from time 18 exp==1
			replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
							if (exp==1 & death85==0 & stroke85==0)
			*Generate stroke indicator for interval 85-90
			gen stroke85to90 = 0 if (death85==0 & stroke85==0)
			replace stroke85to90 = 1 if (stroketime85to90 < 5)
			replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
			*Update prevalent stroke variable (stroke up to age 90)
			replace stroke_history = 1 if (stroke85to90==1)
			replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
			*Generate indicator for stroke before age 90						
			gen stroke90 = 0
			replace stroke90 = 1 if strokeage < 90

			*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
			replace death85to90 = 1 if (strokedeath85to90==1) 
			replace survage = strokeage if (strokedeath85to90==1) 
			*Generate indicator for stroke death before age 90	
			gen strokedeath90 = 0
			replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1)
				
	
			/*Step added for search loop: replace survage = 89.999 for people who are still alive*/ 
			replace survage = 89.999 if (death90==0 & survage==.) 									
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80== 1 | ///
			stroke80to85==11 | stroke85to90==1/* | stroke90to95==1*/)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)	
			
			*ages 65 to 70
			gen strokeage65to70_start = 65 if (strokeage !=. & strokeage>=65)
			gen strokeage65to70_end = 69.99999 if (strokeage !=. & strokeage>=70)
			replace strokeage65to70_end = strokeage if (strokeage !=. & strokeage>=65 & strokeage<70)
			
			*ages 70 to 75
			gen strokeage70to75_start = 70 if (strokeage !=. & strokeage>=70)
			gen strokeage70to75_end = 74.99999 if (strokeage !=. & strokeage>=75)
			replace strokeage70to75_end = strokeage if (strokeage !=. & strokeage>=70 & strokeage<75)
			
			*ages 75 to 80
			gen strokeage75to80_start = 75 if (strokeage !=. & strokeage>=75)
			gen strokeage75to80_end = 79.99999 if (strokeage !=. & strokeage>=80)
			replace strokeage75to80_end = strokeage if (strokeage !=. & strokeage>=75 & strokeage<80)
			
			*ages 80 to 85
			gen strokeage80to85_start = 80 if (strokeage !=. & strokeage>=80)
			gen strokeage80to85_end = 84.99999 if (strokeage !=. & strokeage>=85)
			replace strokeage80to85_end = strokeage if (strokeage !=. & strokeage>=80 & strokeage<85)
		
			*ages 85 to 90
			gen strokeage85to90_start = 85 if (strokeage !=. & strokeage>=85)
			gen strokeage85to90_end = 89.99999 if (strokeage !=. & strokeage>=90)
			replace strokeage85to90_end = strokeage if (strokeage !=. & strokeage>=85 & strokeage<90)
			
		 qui stset strokeage, failure(stroke85to90==1) id(id) enter(strokeage85to90_start) exit(strokeage85to90_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys85to90_exp0 = r(rate)
         if `r(rate)' > `target_strokerate85to90_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_85to90==`stk_lambda_exp0_85to90lt', strokerate1000pys85to90_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_85to90==`stk_lambda_exp0_85to90lt', strokerate1000pys85to90_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_85to90`i' = `stk_lambda_exp0_85to90lt'
}
clear
set obs 5
gen stk_lambda_exp0_85to90=.
forvalues i=1/5 {
   replace stk_lambda_exp0_85to90=`stk_lambda_exp0_85to90`i'' in `i'
}
sum stk_lambda_exp0_85to90
global stk_lambda_exp0_85to90 = `r(mean)'
global stk_lambda_exp0_85to90_min = `r(min)'
global stk_lambda_exp0_85to90_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 90-95	  	***/	 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death90to95_exp0 = $target_p_death90to95_exp0 
   *add lower bound guess here
   local lambda_90to95l = $lambda_90to95l
   quietly forvalues x = 0(.01)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_90to95lt =`lambda_90to95l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_90to95 = `lambda_90to95lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			/*local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			/*local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			/*local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002*/

			*effects of covariates on mortality risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				

			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75
			
			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				

			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80
			
		***b. Stroke
		*Generate stroke time from time 16 exp==0
		gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
						if (exp==0 & death75==0 & stroke75==0)
		*Generate stroke time from time 16 exp==1
		replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
						if (exp==1 & death75==0 & stroke75==0)
		*Generate stroke indicator for interval 75-80
		gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
		replace stroke75to80 = 1 if (stroketime75to80 < 5)
		replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
		*Update prevalent stroke variable (stroke up to age 80)
		replace stroke_history = 1 if (stroke75to80==1)
		replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
		*Generate indicator for stroke before age 80						
		gen stroke80 = 0
		replace stroke80 = 1 if strokeage < 80

		*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
		replace death75to80 = 1 if (strokedeath75to80==1) 
		replace survage = strokeage if (strokedeath75to80==1) 
		*Generate indicator for stroke death before age 80	
		gen strokedeath80 = 0
		replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			

		/*Interval 80-85*/
		***a. Survival
		*Generate survival time from time 17
		gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death80==0)
		*Generate death indicator for interval 80-85
		gen death80to85 = 0 if (death80==0)
		replace death80to85 = 1 if (survtime80to85 < 5)
		replace survage = 80 + survtime80to85 if (death80to85==1)
		*Generate indicator for death before age 85
		gen death85 = 0
		replace death85 = 1 if survage < 85

		***b. Stroke
		*Generate stroke time from time 17 exp==0
		gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
						if (exp==0 & death80==0 & stroke80==0)
		*Generate stroke time from time 17 exp==1
		replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
						if (exp==1 & death80==0 & stroke80==0)
		*Generate stroke indicator for interval 80-85
		gen stroke80to85 = 0 if (death80==0 & stroke80==0)
		replace stroke80to85 = 1 if (stroketime80to85 < 5)
		replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
		*Update prevalent stroke variable (stroke up to age 85)
		replace stroke_history = 1 if (stroke80to85==1)
		replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
		*Generate indicator for stroke before age 85						
		gen stroke85 = 0
		replace stroke85 = 1 if strokeage < 85

		*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
		replace death80to85 = 1 if (strokedeath80to85==1) 
		replace survage = strokeage if (strokedeath80to85==1) 
		*Generate indicator for stroke death before age 85	
		gen strokedeath85 = 0
		replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
			strokedeath80to85==1)
			

		/*Interval 85-90*/
		***a. Survival
		*Generate survival time from time 18
		gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death85==0)
		*Generate death indicator for interval 85-90
		gen death85to90 = 0 if (death85==0)
		replace death85to90 = 1 if (survtime85to90 < 5)
		replace survage = 85 + survtime85to90 if (death85to90==1)
		*Generate indicator for death before age 90
		gen death90 = 0
		replace death90 = 1 if survage < 90

		***b. Stroke
		*Generate stroke time from time 18 exp==0
		gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90'*exp(`b1'*U)) ///
						if (exp==0 & death85==0 & stroke85==0)
		*Generate stroke time from time 18 exp==1
		replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
						if (exp==1 & death85==0 & stroke85==0)
		*Generate stroke indicator for interval 85-90
		gen stroke85to90 = 0 if (death85==0 & stroke85==0)
		replace stroke85to90 = 1 if (stroketime85to90 < 5)
		replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
		*Update prevalent stroke variable (stroke up to age 90)
		replace stroke_history = 1 if (stroke85to90==1)
		replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
		*Generate indicator for stroke before age 90						
		gen stroke90 = 0
		replace stroke90 = 1 if strokeage < 90

		*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
		replace death85to90 = 1 if (strokedeath85to90==1) 
		replace survage = strokeage if (strokedeath85to90==1) 
		*Generate indicator for stroke death before age 90	
		gen strokedeath90 = 0
		replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
			strokedeath80to85==1 | strokedeath85to90==1)
			

		/*Interval 90-95*/
		***a. Survival
		*Generate survival time from time 19
		gen survtime90to95 = -ln(U_90to95)/(`lambda_90to95lt'*exp(`g1_90to95'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death90==0)
		*Generate death indicator for interval 90-95
		gen death90to95 = 0 if (death90==0)
		replace death90to95 = 1 if (survtime90to95 < 5)
		replace survage = 90 + survtime90to95 if (death90to95==1)
		*Generate indicator for death before age 95
		gen death95 = 0
		replace death95 = 1 if survage < 95

		
		 qui sum death90to95 if (exposure==0)																					
		 local p_death90to95_exp0 = r(mean) + (`pstrokedeath'*0.004) //Add in approximate number of stroke deaths
         if `p_death90to95_exp0' > `target_p_death90to95_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_90to95==`lambda_90to95lt', p_death90to95_exp0=`p_death90to95_exp0'
         }
         else {
            noisily di in white "I did NOT stop at lambda_90to95==`lambda_90to95lt', p_death90to95_exp0=`p_death90to95_exp0'
         }
      }
   }
   local lambda_90to95`i' = `lambda_90to95lt'
}
clear
set obs 5
gen lambda_90to95=.
forvalues i=1/5 {
   replace lambda_90to95=`lambda_90to95`i'' in `i'
}
sum lambda_90to95
global lambda_90to95 = `r(mean)'
global lambda_90to95_min = `r(min)'
global lambda_90to95_max = `r(max)'



/**********************************************************************************/
/***		Find baseline STROKE hazard for whites age 90-95	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_strokerate90to95_exp0 = $target_strokerate90to95_exp0
   *add lower bound guess here
   local stk_lambda_exp0_90to95l = $stk_lambda_exp0_90to95l
   quietly forvalues x = 0(.001)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local stk_lambda_exp0_90to95lt =`stk_lambda_exp0_90to95l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "stk_lambda_exp0_90to95 = `stk_lambda_exp0_90to95lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			/*local lambda_95to100 =	$lambda_95to100*/

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			/*local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95*/

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95lt' + 0.002

			*effects of covariates on mortality risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
			
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
	
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
				
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)	
					
					
			/*Interval 85-90*/
			***a. Survival
			*Generate survival time from time 18
			gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death85==0)
			*Generate death indicator for interval 85-90
			gen death85to90 = 0 if (death85==0)
			replace death85to90 = 1 if (survtime85to90 < 5)
			replace survage = 85 + survtime85to90 if (death85to90==1)
			*Generate indicator for death before age 90
			gen death90 = 0
			replace death90 = 1 if survage < 90

			***b. Stroke
			*Generate stroke time from time 18 exp==0
			gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90'*exp(`b1'*U)) ///
							if (exp==0 & death85==0 & stroke85==0)
			*Generate stroke time from time 18 exp==1
			replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
							if (exp==1 & death85==0 & stroke85==0)
			*Generate stroke indicator for interval 85-90
			gen stroke85to90 = 0 if (death85==0 & stroke85==0)
			replace stroke85to90 = 1 if (stroketime85to90 < 5)
			replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
			*Update prevalent stroke variable (stroke up to age 90)
			replace stroke_history = 1 if (stroke85to90==1)
			replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
			*Generate indicator for stroke before age 90						
			gen stroke90 = 0
			replace stroke90 = 1 if strokeage < 90

			*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
			replace death85to90 = 1 if (strokedeath85to90==1) 
			replace survage = strokeage if (strokedeath85to90==1) 
			*Generate indicator for stroke death before age 90	
			gen strokedeath90 = 0
			replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1)
			
			
			/*Interval 90-95*/
			***a. Survival
			*Generate survival time from time 19
			gen survtime90to95 = -ln(U_90to95)/(`lambda_90to95'*exp(`g1_90to95'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death90==0)
			*Generate death indicator for interval 90-95
			gen death90to95 = 0 if (death90==0)
			replace death90to95 = 1 if (survtime90to95 < 5)
			replace survage = 90 + survtime90to95 if (death90to95==1)
			*Generate indicator for death before age 95
			gen death95 = 0
			replace death95 = 1 if survage < 95

			***b. Stroke
			/*Interval 90-95*/
			*Generate stroke time from time 19 exp==0
			gen stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp0_90to95lt'*exp(`b1'*U)) ///
							if (exp==0 & death90==0 & stroke90==0)
			*Generate stroke time from time 19 exp==1
			replace stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp1_90to95'*exp(`b1'*U)) ///
							if (exp==1 & death90==0 & stroke90==0)
			*Generate stroke indicator for interval 90-95
			gen stroke90to95 = 0 if (death90==0 & stroke90==0)
			replace stroke90to95 = 1 if (stroketime90to95 < 5)
			replace stroke90to95 = 0 if (death90to95==1 & stroketime90to95 != . & stroketime90to95 > survtime90to95)
			*Update prevalent stroke variable (stroke up to age 95)
			replace stroke_history = 1 if (stroke90to95==1)
			replace strokeage = 90 + stroketime90to95 if (stroke90to95==1)
			*Generate indicator for stroke before age 95						
			gen stroke95 = 0
			replace stroke95 = 1 if strokeage < 95

			*Generate stroke deaths for interval 90-95, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath90to95 = runiform()<`pstrokedeath' if stroke90to95==1
			replace death90to95 = 1 if (strokedeath90to95==1) 
			replace survage = strokeage if (strokedeath90to95==1) 
			*Generate indicator for stroke death before age 95	
			gen strokedeath95 = 0
			replace strokedeath95 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1 | strokedeath90to95==1)
				
	
			/*Step added for search loop: replace survage = 94.999 for people who are still alive*/ 
			replace survage = 94.999 if (death95==0 & survage==.) 									
								
			*iv. Generate variables for stroke between ages 45 to 95
			gen stroke = 0
			replace stroke = 1 if (stroke45to50==1 | stroke50to55==1 | stroke55to60==1 | /// 
			stroke60to65==1 | stroke65to70==1 | stroke70to75==1 | stroke75to80== 1 | ///
			stroke80to85==11 | stroke85to90==1 | stroke90to95==1)


			*v. Generate strokeage for people who didn't develop a stroke
			replace strokeage = survage if (stroke==0 & death45==0) //no strokeage for people who died before age 45

			
			*vi. Generate age-stratified stroke variables
			*ages 45 to 50
			gen strokeage45to50_start = 45 if (strokeage !=. & strokeage>=45)
			gen strokeage45to50_end = 49.99999 if (strokeage !=. & strokeage>=50)
			replace strokeage45to50_end = strokeage if (strokeage !=. & strokeage>=45 & strokeage<50)

			*ages 50 to 55
			gen strokeage50to55_start = 50 if (strokeage !=. & strokeage>=50)
			gen strokeage50to55_end = 54.99999 if (strokeage !=. & strokeage>=55)
			replace strokeage50to55_end = strokeage if (strokeage !=. & strokeage>=50 & strokeage<55)

			*ages 55 to 60
			gen strokeage55to60_start = 55 if (strokeage !=. & strokeage>=55)
			gen strokeage55to60_end = 59.99999 if (strokeage !=. & strokeage>=60)
			replace strokeage55to60_end = strokeage if (strokeage !=. & strokeage>=55 & strokeage<60)
				
			*ages 60 to 65
			gen strokeage60to65_start = 60 if (strokeage !=. & strokeage>=60)
			gen strokeage60to65_end = 64.99999 if (strokeage !=. & strokeage>=65)
			replace strokeage60to65_end = strokeage if (strokeage !=. & strokeage>=60 & strokeage<65)	
			
			*ages 65 to 70
			gen strokeage65to70_start = 65 if (strokeage !=. & strokeage>=65)
			gen strokeage65to70_end = 69.99999 if (strokeage !=. & strokeage>=70)
			replace strokeage65to70_end = strokeage if (strokeage !=. & strokeage>=65 & strokeage<70)
			
			*ages 70 to 75
			gen strokeage70to75_start = 70 if (strokeage !=. & strokeage>=70)
			gen strokeage70to75_end = 74.99999 if (strokeage !=. & strokeage>=75)
			replace strokeage70to75_end = strokeage if (strokeage !=. & strokeage>=70 & strokeage<75)
			
			*ages 75 to 80
			gen strokeage75to80_start = 75 if (strokeage !=. & strokeage>=75)
			gen strokeage75to80_end = 79.99999 if (strokeage !=. & strokeage>=80)
			replace strokeage75to80_end = strokeage if (strokeage !=. & strokeage>=75 & strokeage<80)
			
			*ages 80 to 85
			gen strokeage80to85_start = 80 if (strokeage !=. & strokeage>=80)
			gen strokeage80to85_end = 84.99999 if (strokeage !=. & strokeage>=85)
			replace strokeage80to85_end = strokeage if (strokeage !=. & strokeage>=80 & strokeage<85)
		
			*ages 85 to 90
			gen strokeage85to90_start = 85 if (strokeage !=. & strokeage>=85)
			gen strokeage85to90_end = 89.99999 if (strokeage !=. & strokeage>=90)
			replace strokeage85to90_end = strokeage if (strokeage !=. & strokeage>=85 & strokeage<90)
			
			*ages 90 to 95
			gen strokeage90to95_start = 90 if (strokeage !=. & strokeage>=90)
			gen strokeage90to95_end = 94.99999 if (strokeage !=. & strokeage>=90)
			replace strokeage90to95_end = strokeage if (strokeage !=. & strokeage>=90 & strokeage<95)

			
		 qui stset strokeage, failure(stroke90to95==1) id(id) enter(strokeage90to95_start) exit(strokeage90to95_end)
		 qui stptime if (exposure==0), title(person-years) per(10000)
		 local strokerate1000pys90to95_exp0 = r(rate)
         if `r(rate)' > `target_strokerate90to95_exp0'+.001 {
            local toolow=0
            noisily di in red "I stopped at stk_lambda_exp0_90to95==`stk_lambda_exp0_90to95lt', strokerate1000pys90to95_exp0=`r(rate)'
         }
         else {
            noisily di in white "I did NOT stop at stk_lambda_exp0_90to95==`stk_lambda_exp0_90to95lt', strokerate1000pys90to95_exp0=`r(rate)'
         }
      }
   }
   local stk_lambda_exp0_90to95`i' = `stk_lambda_exp0_90to95lt'
}
clear
set obs 5
gen stk_lambda_exp0_90to95=.
forvalues i=1/5 {
   replace stk_lambda_exp0_90to95=`stk_lambda_exp0_90to95`i'' in `i'
}
sum stk_lambda_exp0_90to95
global stk_lambda_exp0_90to95 = `r(mean)'
global stk_lambda_exp0_90to95_min = `r(min)'
global stk_lambda_exp0_90to95_max = `r(max)'



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for whites age 95-100  		***/
/**********************************************************************************
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_p_death95to100_exp0 = $target_p_death95to100_exp0
   *add lower bound guess here
   local lambda_95to100l = $lambda_95to100l
   quietly forvalues x = 0(.01)30 { //0(.0001)30 {
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local lambda_95to100lt =`lambda_95to100l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "lambda_95to100 = `lambda_95to100lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		0.30
			local g1_1to5 = 	0.37
			local g1_5to10 = 	0.23
			local g1_10to15 = 	0.56
			local g1_15to20 = 	0.94
			local g1_20to25 = 	0.92
			local g1_25to30 = 	0.79
			local g1_30to35 = 	0.78
			local g1_35to40 = 	0.77
			local g1_40to45 = 	0.74
			local g1_45to50 = 	0.66
			local g1_50to55 = 	0.59
			local g1_55to60 = 	0.48
			local g1_60to65 = 	0.32
			local g1_65to70 = 	0.19
			local g1_70to75 = 	0.10
			local g1_75to80 = 	-0.05
			local g1_80to85 = 	-0.10
			local g1_85to90 = 	-0.19
			local g1_90to95 = 	-0.24
			local g1_95to100 =	-0.33

			*effects of covariates on mortality risk
			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			/*local lambda_95to100 =	$lambda_95to100*/ 

			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95
			
			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50lt' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*effects of covariates on mortality risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)


			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)


			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
	
	
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				

			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75
			
			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				

			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80
			
			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
				

			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)
				

			/*Interval 85-90*/
			***a. Survival
			*Generate survival time from time 18
			gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death85==0)
			*Generate death indicator for interval 85-90
			gen death85to90 = 0 if (death85==0)
			replace death85to90 = 1 if (survtime85to90 < 5)
			replace survage = 85 + survtime85to90 if (death85to90==1)
			*Generate indicator for death before age 90
			gen death90 = 0
			replace death90 = 1 if survage < 90

			***b. Stroke
			*Generate stroke time from time 18 exp==0
			gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90'*exp(`b1'*U)) ///
							if (exp==0 & death85==0 & stroke85==0)
			*Generate stroke time from time 18 exp==1
			replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
							if (exp==1 & death85==0 & stroke85==0)
			*Generate stroke indicator for interval 85-90
			gen stroke85to90 = 0 if (death85==0 & stroke85==0)
			replace stroke85to90 = 1 if (stroketime85to90 < 5)
			replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
			*Update prevalent stroke variable (stroke up to age 90)
			replace stroke_history = 1 if (stroke85to90==1)
			replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
			*Generate indicator for stroke before age 90						
			gen stroke90 = 0
			replace stroke90 = 1 if strokeage < 90

			*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
			replace death85to90 = 1 if (strokedeath85to90==1) 
			replace survage = strokeage if (strokedeath85to90==1) 
			*Generate indicator for stroke death before age 90	
			gen strokedeath90 = 0
			replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1)
				

			/*Interval 90-95*/
			***a. Survival
			*Generate survival time from time 19
			gen survtime90to95 = -ln(U_90to95)/(`lambda_90to95'*exp(`g1_90to95'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death90==0)
			*Generate death indicator for interval 90-95
			gen death90to95 = 0 if (death90==0)
			replace death90to95 = 1 if (survtime90to95 < 5)
			replace survage = 90 + survtime90to95 if (death90to95==1)
			*Generate indicator for death before age 95
			gen death95 = 0
			replace death95 = 1 if survage < 95
			
			***b. Stroke
			/*Interval 90-95*/
			*Generate stroke time from time 19 exp==0
			gen stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp0_90to95'*exp(`b1'*U)) ///
							if (exp==0 & death90==0 & stroke90==0)
			*Generate stroke time from time 19 exp==1
			replace stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp1_90to95'*exp(`b1'*U)) ///
							if (exp==1 & death90==0 & stroke90==0)
			*Generate stroke indicator for interval 90-95
			gen stroke90to95 = 0 if (death90==0 & stroke90==0)
			replace stroke90to95 = 1 if (stroketime90to95 < 5)
			replace stroke90to95 = 0 if (death90to95==1 & stroketime90to95 != . & stroketime90to95 > survtime90to95)
			*Update prevalent stroke variable (stroke up to age 95)
			replace stroke_history = 1 if (stroke90to95==1)
			replace strokeage = 90 + stroketime90to95 if (stroke90to95==1)
			*Generate indicator for stroke before age 95						
			gen stroke95 = 0
			replace stroke95 = 1 if strokeage < 95

			*Generate stroke deaths for interval 90-95, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath90to95 = runiform()<`pstrokedeath' if stroke90to95==1
			replace death90to95 = 1 if (strokedeath90to95==1) 
			replace survage = strokeage if (strokedeath90to95==1) 
			*Generate indicator for stroke death before age 95	
			gen strokedeath95 = 0
			replace strokedeath95 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1 | strokedeath90to95==1)
		

			/*Interval 95-100*/
			***a. Survival
			*Generate survival time from time 20
			gen survtime95to100 = -ln(U_95to100)/(`lambda_95to100lt'*exp(`g1_95to100'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death95==0)
			*Generate death indicator for interval 95-100
			gen death95to100 = 0 if (death95==0)
			replace death95to100 = 1 if (survtime95to100 < 5)
			replace survage = 95 + survtime95to100 if (survtime95to100!=.)
			*Generate indicator for death before age 95
			gen death100 = 0
			replace death100 = 1 if survage < 100

		
		 qui sum death95to100 if (exposure==0)																					
		 local p_death95to100_exp0 = r(mean)
         if `r(mean)' > `target_p_death95to100_exp0'+.0001 {
            local toolow=0
            noisily di in red "I stopped at lambda_95to100==`lambda_95to100lt', p_death95to100_exp0=`r(mean)'
         }
         else {
            noisily di in white "I did NOT stop at lambda_95to100==`lambda_95to100lt', p_death95to100_exp0=`r(mean)'
         }
      }
   }
   local lambda_95to100`i' = `lambda_95to100lt'
}
clear
set obs 5
gen lambda_95to100=.
forvalues i=1/5 {
   replace lambda_95to100=`lambda_95to100`i'' in `i'
}
sum lambda_95to100
global lambda_95to100 = `r(mean)'
global lambda_95to100_min = `r(min)'
global lambda_95to100_max = `r(max)'*/



/***************************************************************************************************************/
/***************************************************************************************************************/
/***		Search loops for effect of race on mortality for ages 45-100									 ***/
/***************************************************************************************************************/
/***************************************************************************************************************/



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 0-1	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_0to1 = $target_g1_0to1
   *add lower bound guess here
   local g1_0to1l = $g1_0to1l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_0to1lt =`g1_0to1l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_0to1 = `g1_0to1lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			/*local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1

	
		 stset survtime0to1, failure(death0to1)
		 qui stcox exposure, nohr
		 local lnHR_0to1 = _b[exposure] 
         if `lnHR_0to1' > `target_g1_0to1'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_0to1==`g1_0to1lt', ln(HR death 0-1)=`lnHR_0to1'
         }
         else {
            noisily di in white "I did NOT stop at g1_0to1==`g1_0to1lt', ln(HR death 0-1)=`lnHR_0to1'
         }
      }
   }
   local g1_0to1`i' = `g1_0to1lt'
}
clear
set obs 5
gen g1_0to1=.
forvalues i=1/5 {
   replace g1_0to1=`g1_0to1`i'' in `i'
}
sum g1_0to1
global g1_0to1 = `r(mean)'	
global g1_0to1_min = `r(min)'	
global g1_0to1_max = `r(max)'	





/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 1-5	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_1to5 = $target_g1_1to5
   *add lower bound guess here
   local g1_1to5l = $g1_1to5l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_1to5lt =`g1_1to5l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_1to5 = `g1_1to5lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			/*local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5

	
		 stset survtime1to5, failure(death1to5)
		 qui stcox exposure, nohr
		 local lnHR_1to5 = _b[exposure] 
         if `lnHR_1to5' > `target_g1_1to5'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_1to5==`g1_1to5lt', ln(HR death 1-5)=`lnHR_1to5'
         }
         else {
            noisily di in white "I did NOT stop at g1_1to5==`g1_1to5lt', ln(HR death 1-5)=`lnHR_1to5'
         }
      }
   }
   local g1_1to5`i' = `g1_1to5lt'
}
clear
set obs 5
gen g1_1to5=.
forvalues i=1/5 {
   replace g1_1to5=`g1_1to5`i'' in `i'
}
sum g1_1to5
global g1_1to5 = `r(mean)'	
global g1_1to5_min = `r(min)'	
global g1_1to5_max = `r(max)'	



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 5-10	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_5to10 = $target_g1_5to10
   *add lower bound guess here
   local g1_5to10l = $g1_5to10l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_5to10lt =`g1_5to10l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_5to10 = `g1_5to10lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			/*local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


		 stset survtime5to10, failure(death5to10)
		 qui stcox exposure, nohr
		 local lnHR_5to10 = _b[exposure] 
         if `lnHR_5to10' > `target_g1_5to10'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_5to10==`g1_5to10lt', ln(HR death 5-10)=`lnHR_5to10'
         }
         else {
            noisily di in white "I did NOT stop at g1_5to10==`g1_5to10lt', ln(HR death 5-10)=`lnHR_5to10'
         }
      }
   }
   local g1_5to10`i' = `g1_5to10lt'
}
clear
set obs 5
gen g1_5to10=.
forvalues i=1/5 {
   replace g1_5to10=`g1_5to10`i'' in `i'
}
sum g1_5to10
global g1_5to10 = `r(mean)'	
global g1_5to10_min = `r(min)'	
global g1_5to10_max = `r(max)'	



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 10-15	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_10to15 = $target_g1_10to15
   *add lower bound guess here
   local g1_10to15l = $g1_10to15l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_10to15lt =`g1_10to15l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_10to15 = `g1_10to15lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			/*local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

	
		 stset survtime10to15, failure(death10to15)
		 qui stcox exposure, nohr
		 local lnHR_10to15 = _b[exposure] 
         if `lnHR_10to15' > `target_g1_10to15'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_10to15==`g1_10to15lt', ln(HR death 10-15)=`lnHR_10to15'
         }
         else {
            noisily di in white "I did NOT stop at g1_10to15==`g1_10to15lt', ln(HR death 10-15)=`lnHR_10to15'
         }
      }
   }
   local g1_10to15`i' = `g1_10to15lt'
}
clear
set obs 5
gen g1_10to15=.
forvalues i=1/5 {
   replace g1_10to15=`g1_10to15`i'' in `i'
}
sum g1_10to15
global g1_10to15 = `r(mean)'	
global g1_10to15_min = `r(min)'	
global g1_10to15_max = `r(max)'		

	
	
/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 15-20	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_15to20 = $target_g1_15to20
   *add lower bound guess here
   local g1_15to20l = $g1_15to20l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_15to20lt =`g1_15to20l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_15to20 = `g1_15to20lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			/*local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			*local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15
			
			
			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20

	
		 stset survtime15to20, failure(death15to20)
		 qui stcox exposure, nohr
		 local lnHR_15to20 = _b[exposure] 
         if `lnHR_15to20' > `target_g1_15to20'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_15to20==`g1_15to20lt', ln(HR death 15-20)=`lnHR_15to20'
         }
         else {
            noisily di in white "I did NOT stop at g1_15to20==`g1_15to20lt', ln(HR death 15-20)=`lnHR_15to20'
         }
      }
   }
   local g1_15to20`i' = `g1_15to20lt'
}
clear
set obs 5
gen g1_15to20=.
forvalues i=1/5 {
   replace g1_15to20=`g1_15to20`i'' in `i'
}
sum g1_15to20
global g1_15to20 = `r(mean)'	
global g1_15to20_min = `r(min)'	
global g1_15to20_max = `r(max)'	




/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 20-25	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_20to25 = $target_g1_20to25
   *add lower bound guess here
   local g1_20to25l = $g1_20to25l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_20to25lt =`g1_20to25l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_20to25 = `g1_20to25lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			/*local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15
			
			
			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20
			
			
			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25

	
		 stset survtime20to25, failure(death20to25)
		 qui stcox exposure, nohr
		 local lnHR_20to25 = _b[exposure] 
         if `lnHR_20to25' > `target_g1_20to25'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_20to25==`g1_20to25lt', ln(HR death 20-25)=`lnHR_20to25'
         }
         else {
            noisily di in white "I did NOT stop at g1_20to25==`g1_20to25lt', ln(HR death 20-25)=`lnHR_20to25'
         }
      }
   }
   local g1_20to25`i' = `g1_20to25lt'
}
clear
set obs 5
gen g1_20to25=.
forvalues i=1/5 {
   replace g1_20to25=`g1_20to25`i'' in `i'
}
sum g1_20to25
global g1_20to25 = `r(mean)'		
global g1_20to25_min = `r(min)'
global g1_20to25_max = `r(max)'
	

	

/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 25-30	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_25to30 = $target_g1_25to30
   *add lower bound guess here
   local g1_25to30l = $g1_25to30l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_25to30lt =`g1_25to30l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_25to30 = `g1_25to30lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			/*local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15
			
			
			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20
			
			
			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25
			
			
			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30

	
		 stset survtime25to30, failure(death25to30)
		 qui stcox exposure, nohr
		 local lnHR_25to30 = _b[exposure] 
         if `lnHR_25to30' > `target_g1_25to30'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_20to25==`g1_25to30lt', ln(HR death 25-30)=`lnHR_25to30'
         }
         else {
            noisily di in white "I did NOT stop at g1_25to30==`g1_25to30lt', ln(HR death 25-30)=`lnHR_25to30'
         }
      }
   }
   local g1_25to30`i' = `g1_25to30lt'
}
clear
set obs 5
gen g1_25to30=.
forvalues i=1/5 {
   replace g1_25to30=`g1_25to30`i'' in `i'
}
sum g1_25to30
global g1_25to30 = `r(mean)'	
global g1_25to30_min = `r(min)'	
global g1_25to30_max = `r(max)'		



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 30-35  		***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_30to35 = $target_g1_30to35
   *add lower bound guess here
   local g1_30to35l = $g1_30to35l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_30to35lt =`g1_30to35l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_30to35 = `g1_30to35lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			/*local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15
			
			
			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20
			
			
			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25
			
			
			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30
			
			
			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35

	
		 stset survtime30to35, failure(death30to35)
		 qui stcox exposure, nohr
		 local lnHR_30to35 = _b[exposure] 
         if `lnHR_30to35' > `target_g1_30to35'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_30to35==`g1_30to35lt', ln(HR death 30-35)=`lnHR_30to35'
         }
         else {
            noisily di in white "I did NOT stop at g1_30to35==`g1_30to35lt', ln(HR death 30-35)=`lnHR_30to35'
         }
      }
   }
   local g1_30to35`i' = `g1_30to35lt'
}
clear
set obs 5
gen g1_30to35=.
forvalues i=1/5 {
   replace g1_30to35=`g1_30to35`i'' in `i'
}
sum g1_30to35
global g1_30to35 = `r(mean)'	
global g1_30to35_min = `r(min)'	
global g1_30to35_max = `r(max)'	



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 35-40  		***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_35to40 = $target_g1_35to40
   *add lower bound guess here
   local g1_35to40l = $g1_35to40l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_35to40lt =`g1_35to40l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_35to40 = `g1_35to40lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			/*local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15
			
			
			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20
			
			
			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25
			
			
			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30
			
			
			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35
			
			
			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40

	
		 stset survtime35to40, failure(death35to40)
		 qui stcox exposure, nohr
		 local lnHR_35to40 = _b[exposure] 
         if `lnHR_35to40' > `target_g1_35to40'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_35to40==`g1_35to40lt', ln(HR death 35-40)=`lnHR_35to40'
         }
         else {
            noisily di in white "I did NOT stop at g1_35to40==`g1_35to40lt', ln(HR death 35-40)=`lnHR_35to40'
         }
      }
   }
   local g1_35to40`i' = `g1_35to40lt'
}
clear
set obs 5
gen g1_35to40=.
forvalues i=1/5 {
   replace g1_35to40=`g1_35to40`i'' in `i'
}
sum g1_35to40
global g1_35to40 = `r(mean)'		
global g1_35to40_min = `r(min)'
global g1_35to40_max = `r(max)'
	


/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 40-45  		***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_40to45 = $target_g1_40to45
   *add lower bound guess here
   local g1_40to45l = $g1_40to45l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_40to45lt =`g1_40to45l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_40to45 = `g1_40to45lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			/*local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1
			
			
			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5
			
			
			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10
			
			
			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15
			
			
			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20
			
			
			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25
			
			
			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30
			
			
			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35
			
			
			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40
			
			
			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45

	
		 stset survtime40to45, failure(death40to45)
		 qui stcox exposure, nohr
		 local lnHR_40to45 = _b[exposure] 
         if `lnHR_40to45' > `target_g1_40to45'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_40to45==`g1_40to45lt', ln(HR death 40-45)=`lnHR_40to45'
         }
         else {
            noisily di in white "I did NOT stop at g1_40to45==`g1_40to45lt', ln(HR death 40-45)=`lnHR_40to45'
         }
      }
   }
   local g1_40to45`i' = `g1_40to45lt'
}
clear
set obs 5
gen g1_40to45=.
forvalues i=1/5 {
   replace g1_40to45=`g1_40to45`i'' in `i'
}
sum g1_40to45
global g1_40to45 = `r(mean)'		
global g1_40to45_min = `r(min)'	
global g1_40to45_max = `r(max)'	



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 45-50	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_45to50 = $target_g1_45to50
   *add lower bound guess here
   local g1_45to50l = $g1_45to50l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_45to50lt =`g1_45to50l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_45to50 = `g1_45to50lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			/*local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50
		
			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)

		
		 stset survtime45to50, failure(death45to50)
		 qui stcox exposure, nohr
		 local lnHR_45to50 = _b[exposure] 
         if `lnHR_45to50' > `target_g1_45to50'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_45to50==`g1_45to50lt', ln(HR death 45-50)=`lnHR_45to50'
         }
         else {
            noisily di in white "I did NOT stop at g1_45to50==`g1_45to50lt', ln(HR death 45-50)=`lnHR_45to50'
         }
      }
   }
   local g1_45to50`i' = `g1_45to50lt'
}
clear
set obs 5
gen g1_45to50=.
forvalues i=1/5 {
   replace g1_45to50=`g1_45to50`i'' in `i'
}
sum g1_45to50
global g1_45to50 = `r(mean)'
global g1_45to50_min = `r(min)'	
global g1_45to50_max = `r(max)'																				



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 50-55	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_50to55 = $target_g1_50to55
   *add lower bound guess here
   local g1_50to55l = $g1_50to55l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_50to55lt =`g1_50to55l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_50to55 = `g1_50to55lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			/*local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			*local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)

		
		 stset survtime50to55, failure(death50to55)
		 qui stcox exposure, nohr
		 local lnHR_50to55 = _b[exposure] 
         if `lnHR_50to55' > `target_g1_50to55'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_50to55==`g1_50to55lt', ln(HR death 50-55)=`lnHR_50to55'
         }
         else {
            noisily di in white "I did NOT stop at g1_50to55==`g1_50to55lt', ln(HR death 50-55)=`lnHR_50to55'
         }
      }
   }
   local g1_50to55`i' = `g1_50to55lt'
}
clear
set obs 5
gen g1_50to55=.
forvalues i=1/5 {
   replace g1_50to55=`g1_50to55`i'' in `i'
}
sum g1_50to55
global g1_50to55 = `r(mean)'
global g1_50to55_min = `r(min)'
global g1_50to55_max = `r(max)'		 



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 55-60	  	***/
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_55to60 = $target_g1_55to60
   *add lower bound guess here
   local g1_55to60l = $g1_55to60l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_55to60lt =`g1_55to60l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_55to60 = `g1_55to60lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			/*local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)

		
		 stset survtime55to60, failure(death55to60)
		 qui stcox exposure, nohr
		 local lnHR_55to60 = _b[exposure] 
         if `lnHR_55to60' > `target_g1_55to60'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_55to60==`g1_55to60lt', ln(HR death 55-60)=`lnHR_55to60'
         }
         else {
            noisily di in white "I did NOT stop at g1_55to60==`g1_55to60lt', ln(HR death 55-60)=`lnHR_55to60'
         }
      }
   }
   local g1_55to60`i' = `g1_55to60lt'
}
clear
set obs 5
gen g1_55to60=.
forvalues i=1/5 {
   replace g1_55to60=`g1_55to60`i'' in `i'
}
sum g1_55to60
global g1_55to60 = `r(mean)'	
global g1_55to60_min = `r(min)'	
global g1_55to60_max = `r(max)'		 



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 60-65  		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_60to65 = $target_g1_60to65
   *add lower bound guess here
   local g1_60to65l = $g1_60to65l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_60to65lt =`g1_60to65l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_60to65 = `g1_60to65lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			/*local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)

		
		 stset survtime60to65, failure(death60to65)
		 qui stcox exposure, nohr
		 local lnHR_60to65 = _b[exposure] 
         if `lnHR_60to65' > `target_g1_60to65'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_60to65==`g1_60to65lt', ln(HR death 60-65)=`lnHR_60to65'
         }
         else {
            noisily di in white "I did NOT stop at g1_60to65==`g1_60to65lt', ln(HR death 60-65)=`lnHR_60to65'
         }
      }
   }
   local g1_60to65`i' = `g1_60to65lt'
}
clear
set obs 5
gen g1_60to65=.
forvalues i=1/5 {
   replace g1_60to65=`g1_60to65`i'' in `i'
}
sum g1_60to65
global g1_60to65 = `r(mean)'	
global g1_60to65_min = `r(min)'
global g1_60to65_max = `r(max)'	




/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 65-70 		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_65to70 = $target_g1_65to70
   *add lower bound guess here
   local g1_65to70l = $g1_65to70l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_65to70lt =`g1_65to70l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_65to70 = `g1_65to70lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			/*local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)

		
		 stset survtime65to70, failure(death65to70)
		 qui stcox exposure, nohr
		 local lnHR_65to70 = _b[exposure] 
         if `lnHR_65to70' > `target_g1_65to70'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_65to70==`g1_65to70lt', ln(HR death 65-70)=`lnHR_65to70'
         }
         else {
            noisily di in white "I did NOT stop at g1_65to70==`g1_65to70lt', ln(HR death 65-70)=`lnHR_65to70'
         }
      }
   }
   local g1_65to70`i' = `g1_65to70lt'
}
clear
set obs 5
gen g1_65to70=.
forvalues i=1/5 {
   replace g1_65to70=`g1_65to70`i'' in `i'
}
sum g1_65to70
global g1_65to70 = `r(mean)'	
global g1_65to70_min = `r(min)'	
global g1_65to70_max = `r(max)'		  



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 70-75		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_70to75 = $target_g1_70to75
   *add lower bound guess here
   local g1_70to75l = $g1_70to75l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_70to75lt =`g1_70to75l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_70to75 = `g1_70to75lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			/*local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
				
			/*Interval 70-75*/
		***a. Survival
		*Generate survival time from time 15
		gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75lt'*exposure +`g2'*U ///
					+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
						if (death70==0)
		*Generate death indicator for interval 70-75
		gen death70to75 = 0 if (death70==0)
		replace death70to75 = 1 if (survtime70to75 < 5) 
		replace survage = 70 + survtime70to75 if (death70to75==1)
		*Generate indicator for death before age 75
		gen death75 = 0
		replace death75 = 1 if survage < 75

		***b. Stroke
		*Generate stroke time from time 15 exp==0
		gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
						if (exp==0 & death70==0 & stroke70==0)
		*Generate stroke time from time 15 exp==1
		replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
						if (exp==1 & death70==0 & stroke70==0)
		*Generate stroke indicator for interval 70-75
		gen stroke70to75 = 0 if (death70==0 & stroke70==0)
		replace stroke70to75 = 1 if (stroketime70to75 < 5) 
		replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
		*Update prevalent stroke variable (stroke up to age 75)
		replace stroke_history = 1 if (stroke70to75==1)
		replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
		*Generate indicator for stroke before age 75						
		gen stroke75 = 0
		replace stroke75 = 1 if strokeage < 75

		*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
		gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
		replace death70to75 = 1 if (strokedeath70to75==1) 
		replace survage = strokeage if (strokedeath70to75==1) 
		*Generate indicator for stroke death before age 75	
		gen strokedeath75 = 0
		replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
			strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)

		
		 stset survtime70to75, failure(death70to75)
		 qui stcox exposure, nohr
		 local lnHR_70to75 = _b[exposure] 
         if `lnHR_70to75' > `target_g1_70to75'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_70to75==`g1_70to75lt', ln(HR death 70-75)=`lnHR_70to75'
         }
         else {
            noisily di in white "I did NOT stop at g1_70to75==`g1_70to75lt', ln(HR death 70-75)=`lnHR_70to75'
         }
      }
   }
   local g1_70to75`i' = `g1_70to75lt'
}
clear
set obs 5
gen g1_70to75=.
forvalues i=1/5 {
   replace g1_70to75=`g1_70to75`i'' in `i'
}
sum g1_70to75
global g1_70to75 = `r(mean)'	
global g1_70to75_min = `r(min)'		  
global g1_70to75_max = `r(max)'		  	  



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 75-80		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_75to80 = $target_g1_75to80
   *add lower bound guess here
   local g1_75to80l = $g1_75to80l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_75to80lt =`g1_75to80l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_75to80 = `g1_75to80lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			/*local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
				
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				
				
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
	
		
		 stset survtime75to80, failure(death75to80)
		 qui stcox exposure, nohr
		 local lnHR_75to80 = _b[exposure] 
         if `lnHR_75to80' > `target_g1_75to80'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_75to80==`g1_75to80lt', ln(HR death 75-80)=`lnHR_75to80'
         }
         else {
            noisily di in white "I did NOT stop at g1_75to80==`g1_75to80lt', ln(HR death 75-80)=`lnHR_75to80'
         }
      }
   }
   local g1_75to80`i' = `g1_75to80lt'
}
clear
set obs 5
gen g1_75to80=.
forvalues i=1/5 {
   replace g1_75to80=`g1_75to80`i'' in `i'
}
sum g1_75to80
global g1_75to80 = `r(mean)'	
global g1_75to80_min = `r(min)'	
global g1_75to80_max = `r(max)'		  



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 80-85		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_80to85 = $target_g1_80to85
   *add lower bound guess here
   local g1_80to85l = $g1_80to85l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_80to85lt =`g1_80to85l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_80to85 = `g1_80to85lt'"
         clear
		 /*Create blank data set*/
			set obs 20000 //creates blank dataset with XXXXX observations
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			/*local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
				
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				
				
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			
			
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)
	
	
		 stset survtime80to85, failure(death80to85)
		 qui stcox exposure, nohr
		 local lnHR_80to85 = _b[exposure] 
         if `lnHR_80to85' > `target_g1_80to85'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_80to85==`g1_80to85lt', ln(HR death 80-85)=`lnHR_80to85'
         }
         else {
            noisily di in white "I did NOT stop at g1_80to85==`g1_80to85lt', ln(HR death 80-85)=`lnHR_80to85'
         }
      }
   }
   local g1_80to85`i' = `g1_80to85lt'
}
clear
set obs 5
gen g1_80to85=.
forvalues i=1/5 {
   replace g1_80to85=`g1_80to85`i'' in `i'
}
sum g1_80to85
global g1_80to85 = `r(mean)'	
global g1_80to85_min = `r(min)'
global g1_80to85_max = `r(max)'	



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 85-90		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_85to90 = $target_g1_85to90 
   *add lower bound guess here
   local g1_85to90l = $g1_85to90l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_85to90lt =`g1_85to90l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_85to90 = `g1_85to90lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			/*local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
				
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				
				
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			
			
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)
			
			
			/*Interval 85-90*/
			***a. Survival
			*Generate survival time from time 18
			gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death85==0)
			*Generate death indicator for interval 85-90
			gen death85to90 = 0 if (death85==0)
			replace death85to90 = 1 if (survtime85to90 < 5)
			replace survage = 85 + survtime85to90 if (death85to90==1)
			*Generate indicator for death before age 90
			gen death90 = 0
			replace death90 = 1 if survage < 90

			***b. Stroke
			*Generate stroke time from time 18 exp==0
			gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90'*exp(`b1'*U)) ///
							if (exp==0 & death85==0 & stroke85==0)
			*Generate stroke time from time 18 exp==1
			replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
							if (exp==1 & death85==0 & stroke85==0)
			*Generate stroke indicator for interval 85-90
			gen stroke85to90 = 0 if (death85==0 & stroke85==0)
			replace stroke85to90 = 1 if (stroketime85to90 < 5)
			replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
			*Update prevalent stroke variable (stroke up to age 90)
			replace stroke_history = 1 if (stroke85to90==1)
			replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
			*Generate indicator for stroke before age 90						
			gen stroke90 = 0
			replace stroke90 = 1 if strokeage < 90

			*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
			replace death85to90 = 1 if (strokedeath85to90==1) 
			replace survage = strokeage if (strokedeath85to90==1) 
			*Generate indicator for stroke death before age 90	
			gen strokedeath90 = 0
			replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1)
	
	
		 stset survtime85to90, failure(death85to90)
		 qui stcox exposure, nohr
		 local lnHR_85to90 = _b[exposure] 
         if `lnHR_85to90' > `target_g1_85to90'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_85to90==`g1_85to90lt', ln(HR death 85-90)=`lnHR_85to90'
         }
         else {
            noisily di in white "I did NOT stop at g1_85to90==`g1_85to90lt', ln(HR death 85-90)=`lnHR_85to90'
         }
      }
   }
   local g1_85to90`i' = `g1_85to90lt'
}
clear
set obs 5
gen g1_85to90=.
forvalues i=1/5 {
   replace g1_85to90=`g1_85to90`i'' in `i'
}
sum g1_85to90
global g1_85to90 = `r(mean)'	
global g1_85to90_min = `r(min)'	
global g1_85to90_max = `r(max)'		    



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 90-95		***/ 
/**********************************************************************************/
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_90to95 = $target_g1_90to95
   *add lower bound guess here
   local g1_90to95l = $g1_90to95l 
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_90to95lt =`g1_90to95l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_90to95 = `g1_90to95lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			/*local g1_90to95 = 	$g1_90to95
			local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
				
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				
				
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			
			
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)
			
			
			/*Interval 85-90*/
			***a. Survival
			*Generate survival time from time 18
			gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death85==0)
			*Generate death indicator for interval 85-90
			gen death85to90 = 0 if (death85==0)
			replace death85to90 = 1 if (survtime85to90 < 5)
			replace survage = 85 + survtime85to90 if (death85to90==1)
			*Generate indicator for death before age 90
			gen death90 = 0
			replace death90 = 1 if survage < 90

			***b. Stroke
			*Generate stroke time from time 18 exp==0
			gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90'*exp(`b1'*U)) ///
							if (exp==0 & death85==0 & stroke85==0)
			*Generate stroke time from time 18 exp==1
			replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
							if (exp==1 & death85==0 & stroke85==0)
			*Generate stroke indicator for interval 85-90
			gen stroke85to90 = 0 if (death85==0 & stroke85==0)
			replace stroke85to90 = 1 if (stroketime85to90 < 5)
			replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
			*Update prevalent stroke variable (stroke up to age 90)
			replace stroke_history = 1 if (stroke85to90==1)
			replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
			*Generate indicator for stroke before age 90						
			gen stroke90 = 0
			replace stroke90 = 1 if strokeage < 90

			*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
			replace death85to90 = 1 if (strokedeath85to90==1) 
			replace survage = strokeage if (strokedeath85to90==1) 
			*Generate indicator for stroke death before age 90	
			gen strokedeath90 = 0
			replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1)
			
			
			/*Interval 90-95*/
			***a. Survival
			*Generate survival time from time 19
			gen survtime90to95 = -ln(U_90to95)/(`lambda_90to95'*exp(`g1_90to95lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death90==0)
			*Generate death indicator for interval 90-95
			gen death90to95 = 0 if (death90==0)
			replace death90to95 = 1 if (survtime90to95 < 5)
			replace survage = 90 + survtime90to95 if (death90to95==1)
			*Generate indicator for death before age 95
			gen death95 = 0
			replace death95 = 1 if survage < 95

			***b. Stroke
			/*Interval 90-95*/
			*Generate stroke time from time 19 exp==0
			gen stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp0_90to95'*exp(`b1'*U)) ///
							if (exp==0 & death90==0 & stroke90==0)
			*Generate stroke time from time 19 exp==1
			replace stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp1_90to95'*exp(`b1'*U)) ///
							if (exp==1 & death90==0 & stroke90==0)
			*Generate stroke indicator for interval 90-95
			gen stroke90to95 = 0 if (death90==0 & stroke90==0)
			replace stroke90to95 = 1 if (stroketime90to95 < 5)
			replace stroke90to95 = 0 if (death90to95==1 & stroketime90to95 != . & stroketime90to95 > survtime90to95)
			*Update prevalent stroke variable (stroke up to age 95)
			replace stroke_history = 1 if (stroke90to95==1)
			replace strokeage = 90 + stroketime90to95 if (stroke90to95==1)
			*Generate indicator for stroke before age 95						
			gen stroke95 = 0
			replace stroke95 = 1 if strokeage < 95

			*Generate stroke deaths for interval 90-95, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath90to95 = runiform()<`pstrokedeath' if stroke90to95==1
			replace death90to95 = 1 if (strokedeath90to95==1) 
			replace survage = strokeage if (strokedeath90to95==1) 
			*Generate indicator for stroke death before age 95	
			gen strokedeath95 = 0
			replace strokedeath95 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1 | strokedeath90to95==1)
	
	
		 stset survtime90to95, failure(death90to95)
		 qui stcox exposure, nohr
		 local lnHR_90to95 = _b[exposure] 
         if `lnHR_90to95' > `target_g1_90to95'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_90to95==`g1_90to95lt', ln(HR death 90-95)=`lnHR_90to95'
         }
         else {
            noisily di in white "I did NOT stop at g1_90to95==`g1_90to95lt', ln(HR death 90-95)=`lnHR_90to95'
         }
      }
   }
   local g1_90to95`i' = `g1_90to95lt'
}
clear
set obs 5
gen g1_90to95=.
forvalues i=1/5 {
   replace g1_90to95=`g1_90to95`i'' in `i'
}
sum g1_90to95
global g1_90to95 = `r(mean)'		 
global g1_90to95_min = `r(min)'		    
global g1_90to95_max = `r(max)'		    
   



/**********************************************************************************/
/***		Find baseline MORTALITY hazard for blacks age 95-100		***/ 
/**********************************************************************************
clear

forvalues i=1/5 {
   clear
   local toolow =1
   local target_g1_95to100 = $target_g1_95to100
   *add lower bound guess here
   local g1_95to100l = $g1_95to100l
   quietly forvalues x = 0(.01)30 { 
      if `toolow'==1 {
         local seed = 8675309 + `i'
		 set seed `seed'
         local g1_95to100lt =`g1_95to100l'+`x'
         noisily dis "***next increment***"
         noisily dis "increment = " `x'
         noisily di in green "g1_95to100 = `g1_95to100lt'"
         clear
		 /*Create blank data set*/
			set obs 40000 //creates blank dataset with XXXXX observations ***increased sample size at ages 85+ search loop for better estimate
			gen id = _n


			/*Step 1: Set parameters*/

			*specify prevalence of exposure
			local pexp = $pexp 	

			*parameters for Sij 
			//effect of exposure on log hazard of death, based on US life tables for 1919-1921 birth cohort 
			local g1_0to1 =		$g1_0to1
			local g1_1to5 = 	$g1_1to5
			local g1_5to10 = 	$g1_5to10
			local g1_10to15 = 	$g1_10to15
			local g1_15to20 = 	$g1_15to20
			local g1_20to25 = 	$g1_20to25
			local g1_25to30 = 	$g1_25to30
			local g1_30to35 = 	$g1_30to35
			local g1_35to40 = 	$g1_35to40
			local g1_40to45 = 	$g1_40to45
			local g1_45to50 = 	$g1_45to50
			local g1_50to55 = 	$g1_50to55
			local g1_55to60 = 	$g1_55to60
			local g1_60to65 = 	$g1_60to65
			local g1_65to70 = 	$g1_65to70
			local g1_70to75 = 	$g1_70to75
			local g1_75to80 = 	$g1_75to80
			local g1_80to85 = 	$g1_80to85
			local g1_85to90 = 	$g1_85to90
			local g1_90to95 = 	$g1_90to95
			/*local g1_95to100 =	$g1_95to100*/ 

			local g2 = $g2 //effect of U on log hazard of death
			local g3 = $g3 //interaction effect of exposure & U on log hazard of death	
			local g4 = $g4 //annual increase in mortality risk (probably won't use and will delete later)	
			local g5 = $g5 //delete later	

			*baseline hazard of death (whites), based on US life tables for 1919-1921 birth cohort
			local lambda_0to1 = 	$lambda_0to1
			local lambda_1to5 = 	$lambda_1to5
			local lambda_5to10 = 	$lambda_5to10
			local lambda_10to15 = 	$lambda_10to15
			local lambda_15to20 = 	$lambda_15to20
			local lambda_20to25 = 	$lambda_20to25
			local lambda_25to30 = 	$lambda_25to30
			local lambda_30to35 = 	$lambda_30to35
			local lambda_35to40 = 	$lambda_35to40
			local lambda_40to45 = 	$lambda_40to45
			local lambda_45to50 = 	$lambda_45to50
			local lambda_50to55 = 	$lambda_50to55
			local lambda_55to60 =	$lambda_55to60
			local lambda_60to65 = 	$lambda_60to65
			local lambda_65to70 =  	$lambda_65to70
			local lambda_70to75 = 	$lambda_70to75 
			local lambda_75to80 = 	$lambda_75to80
			local lambda_80to85 = 	$lambda_80to85
			local lambda_85to90 = 	$lambda_85to90
			local lambda_90to95 = 	$lambda_90to95
			*local lambda_95to100 =	$lambda_95to100


			*baseline hazard of stroke (exp=0 whites), based on Howard Ann Neurol 2011
			local stk_lambda_exp0_45to50 = 	$stk_lambda_exp0_45to50
			local stk_lambda_exp0_50to55 = 	$stk_lambda_exp0_50to55 
			local stk_lambda_exp0_55to60 =	$stk_lambda_exp0_55to60
			local stk_lambda_exp0_60to65 = 	$stk_lambda_exp0_60to65
			local stk_lambda_exp0_65to70 =  $stk_lambda_exp0_65to70
			local stk_lambda_exp0_70to75 = 	$stk_lambda_exp0_70to75
			local stk_lambda_exp0_75to80 = 	$stk_lambda_exp0_75to80
			local stk_lambda_exp0_80to85 = 	$stk_lambda_exp0_80to85
			local stk_lambda_exp0_85to90 = 	$stk_lambda_exp0_85to90
			local stk_lambda_exp0_90to95 = 	$stk_lambda_exp0_90to95

			*baseline hazard of stroke (exp=1 blacks), based on Howard Ann Neurol 2011
			local stk_lambda_exp1_45to50 = 	`stk_lambda_exp0_45to50' + 0.002
			local stk_lambda_exp1_50to55 = 	`stk_lambda_exp0_50to55' + 0.002
			local stk_lambda_exp1_55to60 =	`stk_lambda_exp0_55to60' + 0.002
			local stk_lambda_exp1_60to65 = 	`stk_lambda_exp0_60to65' + 0.002
			local stk_lambda_exp1_65to70 =  `stk_lambda_exp0_65to70' + 0.002
			local stk_lambda_exp1_70to75 = 	`stk_lambda_exp0_70to75' + 0.002
			local stk_lambda_exp1_75to80 = 	`stk_lambda_exp0_75to80' + 0.002
			local stk_lambda_exp1_80to85 = 	`stk_lambda_exp0_80to85' + 0.002
			local stk_lambda_exp1_85to90 = 	`stk_lambda_exp0_85to90' + 0.002
			local stk_lambda_exp1_90to95 = 	`stk_lambda_exp0_90to95' + 0.002

			*parameters for stroke risk
			local b1 = $b1 			//log (HR) for U on stroke
			local b2 = $b2 			//delete later
			local b3 = $b3 			//delete later
			local b4 = $b3 			//delete later
			local b5 = $b4 			//delete later

			*probability of death at stroke
			local pstrokedeath = $pstrokedeath 


			/*Step 2: Generate exposure variable by generating a U(0,1) distribution 
			and setting exposure=0 if the random number < pexp, otherwise exposure=1*/
			gen exposure = runiform()<`pexp' 
			/*generates random numbers~U(0,1). exposure=0 if random number< pexp, else exposure=1*/


			/*Step 3: Generate continuous time-constant confounder of death and stroke (U)*/
			gen U = rnormal(0,1)


			/*Step 4: Generate survival time for each person and strokes for people alive
			at each interval. 
			a. Each person's underlying time to death is generated for each age interval, 
			conditional on the past provided the person has not died in a previous interval, 
			under an exponential survival distribtion. If the persons generated survival 
			time exceeds the length of the interval between study visits j and j+1, 
			she is considered alive at study visit j+1 and a new survival time is 
			generated for the next interval conditional on history up to the start of the 
			interval, and the process is repeated until the persons survival time falls 
			within a given interval or the end of the study, whichever comes first. Each 
			persons hazard function is defined as:
			h(tij|x) = lambda*exp(g1*exposurei + g2*Ui + g3*exposurei*Ui + g4*stroke_historyi)
			A persons survival time for a given time interval at risk is generated using 
			the inverse cumulative hazard function transformation formula described by 
			Bender et al. (Stat Med 2011)
			b. Stroke code is adapted for survival time code.*/

			*ia. Generate uniform random variable for generating survival time
			gen U_0to1 = runiform()
			gen U_1to5 = runiform()
			gen U_5to10 = runiform()
			gen U_10to15 = runiform()
			gen U_15to20 = runiform()
			gen U_20to25 = runiform()
			gen U_25to30 = runiform()
			gen U_30to35 = runiform()
			gen U_35to40 = runiform()
			gen U_40to45 = runiform()
			gen U_45to50 = runiform()
			gen U_50to55 = runiform()
			gen U_55to60 = runiform()
			gen U_60to65 = runiform()
			gen U_65to70 = runiform()
			gen U_70to75 = runiform()
			gen U_75to80 = runiform()
			gen U_80to85 = runiform()
			gen U_85to90 = runiform()
			gen U_90to95 = runiform()
			gen U_95to100 = runiform()

			*ib. Generate uniform random variable for generating stroke time
			gen U2_45to50 = runiform()
			gen U2_50to55 = runiform()
			gen U2_55to60 = runiform()
			gen U2_60to65 = runiform()
			gen U2_65to70 = runiform()
			gen U2_70to75 = runiform()
			gen U2_75to80 = runiform()
			gen U2_80to85 = runiform()
			gen U2_85to90 = runiform()
			gen U2_90to95 = runiform()


			*ii. Generate survival time and stroke time for each interval

			gen survage = .
			gen strokeage = .

			/***Ages 0-45: no strokes, so only need to generate survival***/
			/*Interval 0-1*/
			*Generate survival time from time 0
			gen survtime0to1 = -ln(U_0to1)/(`lambda_0to1'*exp(`g1_0to1'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0))
			*Generate death indicator for interval 0-1 
			gen death0to1 = 0
			replace death0to1 = 1 if (survtime0to1 < 1) 
			replace survage = survtime0to1 if (death0to1==1)
			*Generate indicator for death before age 1
			gen death1 = death0to1


			/*Interval 1-5*/
			*Generate survival time from time 1
			gen survtime1to5 = -ln(U_1to5)/(`lambda_1to5'*exp(`g1_1to5'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death0to1==0) 
			*Generate death indicator for interval 1-5 
			gen death1to5 = 0 if (death1==0)
			replace death1to5 = 1 if (survtime1to5 < 5)
			replace survage = 1 + survtime1to5 if (death1to5==1)
			*Generate indicator for death before age 5
			gen death5 = 0
			replace death5 = 1 if survage < 5


			/*Interval 5-10*/
			*Generate survival time from time 2
			gen survtime5to10 = -ln(U_5to10)/(`lambda_5to10'*exp(`g1_5to10'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death5==0)
			*Generate death indicator for interval 5-10 
			gen death5to10 = 0 if (death5==0)
			replace death5to10 = 1 if (survtime5to10 < 5) 
			replace survage = 5 + survtime5to10 if (death5to10==1)
			*Generate indicator for death before age 10
			gen death10 = 0
			replace death10 = 1 if survage < 10


			/*Interval 10-15*/
			*Generate survival time from time 3
			gen survtime10to15 = -ln(U_10to15)/(`lambda_10to15'*exp(`g1_10to15'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death10==0)
			*Generate death indicator for interval 10-15
			gen death10to15 = 0 if (death10==0)
			replace death10to15 = 1 if (survtime10to15 < 5) 
			replace survage = 10 + survtime10to15 if (death10to15==1)
			*Generate indicator for death before age 15
			gen death15 = 0
			replace death15 = 1 if survage < 15

			/*Interval 15-20*/
			*Generate survival time from time 4
			gen survtime15to20 = -ln(U_15to20)/(`lambda_15to20'*exp(`g1_15to20'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death15==0)
			*Generate death indicator for interval 15-20 
			gen death15to20 = 0 if (death15==0)
			replace death15to20 = 1 if (survtime15to20 < 5) 
			replace survage = 15 + survtime15to20 if (death15to20==1)
			*Generate indicator for death before age 20
			gen death20 = 0
			replace death20 = 1 if survage < 20


			/*Interval 20-25*/
			*Generate survival time from time 5
			gen survtime20to25 = -ln(U_20to25)/(`lambda_20to25'*exp(`g1_20to25'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death20==0) 
			*Generate death indicator for interval 20-25 
			gen death20to25 = 0 if (death20==0)
			replace death20to25 = 1 if (survtime20to25 < 5) 
			replace survage = 20 + survtime20to25 if (death20to25==1)
			*Generate indicator for death before age 25
			gen death25 = 0
			replace death25 = 1 if survage < 25


			/*Interval 25-30*/
			*Generate survival time from time 6
			gen survtime25to30 = -ln(U_25to30)/(`lambda_25to30'*exp(`g1_25to30'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death25==0)
			*Generate death indicator for interval 25-30 
			gen death25to30 = 0 if (death25==0)
			replace death25to30 = 1 if (survtime25to30 < 5) 
			replace survage = 25 + survtime25to30 if (death25to30==1)
			*Generate indicator for death before age 30
			gen death30 = 0
			replace death30 = 1 if survage < 30


			/*Interval 30-35*/
			*Generate survival time from time 7
			gen survtime30to35 = -ln(U_30to35)/(`lambda_30to35'*exp(`g1_30to35'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death30==0)
			*Generate death indicator for interval 30-35 
			gen death30to35 = 0 if (death30==0)
			replace death30to35 = 1 if (survtime30to35 < 5) 
			replace survage = 30 + survtime30to35 if (death30to35==1)
			*Generate indicator for death before age 35
			gen death35 = 0
			replace death35 = 1 if survage < 35


			/*Interval 35-40*/
			*Generate survival time from time 8
			gen survtime35to40 = -ln(U_35to40)/(`lambda_35to40'*exp(`g1_35to40'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death35==0)
			*Generate death indicator for interval 35-40
			gen death35to40 = 0 if (death35==0)
			replace death35to40 = 1 if (survtime35to40 < 5) 
			replace survage = 35 + survtime35to40 if (death35to40==1)
			*Generate indicator for death before age 40
			gen death40 = 0
			replace death40 = 1 if survage < 40


			/*Interval 40-45*/
			*Generate survival time from time 9
			gen survtime40to45 = -ln(U_40to45)/(`lambda_40to45'*exp(`g1_40to45'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death40==0)
			*Generate death indicator for interval 40-45
			gen death40to45 = 0 if (death40==0)
			replace death40to45 = 1 if (survtime40to45 < 5) 
			replace survage = 40 + survtime40to45 if (death40to45==1)
			*Generate indicator for death before age 45
			gen death45 = 0
			replace death45 = 1 if survage < 45


			/***Starting at age 45--people are at risk of stroke, and prevalent stroke
			increases mortality risk, so we need to start iteratively generating survival 
			times and strokes***/


			/*Interval 45-50*/
			***a. Survival
			*Generate survival time from time 10
			gen survtime45to50 = -ln(U_45to50)/(`lambda_45to50'*exp(`g1_45to50'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*0 +`g5'*0)) ///
							if (death45==0)
			*Generate death indicator for interval 45-50
			gen death45to50 = 0 if (death45==0)
			replace death45to50 = 1 if (survtime45to50 < 5) 
			replace survage = 45 + survtime45to50 if (death45to50==1)
			*Generate indicator for death before age 50
			gen death50 = 0
			replace death50 = 1 if survage < 50

			***b. Stroke
			*Generate stroke time from time 10 exp==0
			gen stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp0_45to50'*exp(`b1'*U)) ///
							if (exp==0 & death45==0)
			*Generate stroke time from time 10 exp==1
			replace stroketime45to50 = -ln(U2_45to50)/(`stk_lambda_exp1_45to50'*exp(`b1'*U)) ///
							if (exp==1 & death45==0)
			*Generate stroke indicator for interval 45-50
			gen stroke45to50 = 0 if (death45==0)
			replace stroke45to50 = 1 if (stroketime45to50 < 5) 
			replace stroke45to50 = 0 if (death45to50==1 & stroketime45to50 != . & stroketime45to50 > survtime45to50)
			*Generate prevalent stroke variable (stroke up to age 50)
			gen stroke_history = 0
			replace stroke_history = 1 if (stroke45to50==1)
			replace strokeage = 45 + stroketime45to50 if (stroke45to50==1)
			*Generate indicator for stroke before age 50
			gen stroke50 = 0
			replace stroke50 = 1 if strokeage < 50

			*Generate stroke deaths for interval 45-50, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath45to50 = runiform()<`pstrokedeath' if stroke45to50==1
			replace death45to50 = 1 if (strokedeath45to50==1) 
			replace survage = strokeage if (strokedeath45to50==1) 
			*Generate indicator for stroke death before age 50
			gen strokedeath50 = 0
			replace strokedeath50 = 1 if (strokedeath45to50==1)


			/*Interval 50-55*/
			***a. Survival
			*Generate survival time from time 11
			gen survtime50to55 = -ln(U_50to55)/(`lambda_50to55'*exp(`g1_50to55'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death50==0)
			*Generate death indicator for interval 50-55
			gen death50to55 = 0 if (death50==0)
			replace death50to55 = 1 if (survtime50to55 < 5) 
			replace survage = 50 + survtime50to55 if (death50to55==1)
			*Generate indicator for death before age 55
			gen death55 = 0
			replace death55 = 1 if survage < 55

			***b. Stroke
			/*Interval 50-55*/
			*Generate stroke time from time 11 exp==0
			gen stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp0_50to55'*exp(`b1'*U)) ///
							if (exp==0 & death50==0 & stroke50==0)
			*Generate stroke time from time 11 exp==1
			replace stroketime50to55 = -ln(U2_50to55)/(`stk_lambda_exp1_50to55'*exp(`b1'*U)) ///
							if (exp==1 & death50==0 & stroke50==0)
			*Generate stroke indicator for interval 50-55
			gen stroke50to55 = 0 if (death50==0 & stroke50==0)
			replace stroke50to55 = 1 if (stroketime50to55 < 5) 
			replace stroke50to55 = 0 if (death50to55==1 & stroketime50to55 != . & stroketime50to55 > survtime50to55)
			*Update prevalent stroke variable (stroke up to age 55)
			replace stroke_history = 1 if (stroke50to55==1)
			replace strokeage = 50 + stroketime50to55 if (stroke50to55==1)
			*Generate indicator for stroke before age 55						
			gen stroke55 = 0
			replace stroke55 = 1 if strokeage < 55

			*Generate stroke deaths for interval 50-55, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath50to55 = runiform()<`pstrokedeath' if stroke50to55==1
			replace death50to55 = 1 if (strokedeath50to55==1) 
			replace survage = strokeage if (strokedeath50to55==1) 
			*Generate indicator for stroke death before age 55
			gen strokedeath55 = 0
			replace strokedeath55 = 1 if (strokedeath45to50==1 | strokedeath50to55==1)
			
			
			/*Interval 55-60*/
			***a. Survival
			*Generate survival time from time 12
			gen survtime55to60 = -ln(U_55to60)/(`lambda_55to60'*exp(`g1_55to60'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death55==0)
			*Generate death indicator for interval 55-60
			gen death55to60 = 0 if (death55==0)
			replace death55to60 = 1 if (survtime55to60 < 5) 
			replace survage = 55 + survtime55to60 if (death55to60==1)
			*Generate indicator for death before age 60
			gen death60 = 0
			replace death60 = 1 if survage < 60

			***b. Stroke
			*Generate stroke time from time 12 exp==0
			gen stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp0_55to60'*exp(`b1'*U)) ///
							if (exp==0 & death55==0 & stroke55==0)
			*Generate stroke time from time 12 exp==1
			replace stroketime55to60 = -ln(U2_55to60)/(`stk_lambda_exp1_55to60'*exp(`b1'*U)) ///
							if (exp==1 & death55==0 & stroke55==0)				
			*Generate stroke indicator for interval 55-60
			gen stroke55to60 = 0 if (death55==0 & stroke55==0)
			replace stroke55to60 = 1 if (stroketime55to60 < 5) 
			replace stroke55to60 = 0 if (death55to60==1 & stroketime55to60 != . & stroketime55to60 > survtime55to60)
			*Update prevalent stroke variable (stroke up to age 60)
			replace stroke_history = 1 if (stroke55to60==1)
			replace strokeage = 55 + stroketime55to60 if (stroke55to60==1)
			*Generate indicator for stroke before age 60						
			gen stroke60 = 0
			replace stroke60 = 1 if strokeage < 60

			*Generate stroke deaths for interval 55-60, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath55to60 = runiform()<`pstrokedeath' if stroke55to60==1
			replace death55to60 = 1 if (strokedeath55to60==1) 
			replace survage = strokeage if (strokedeath55to60==1) 
			*Generate indicator for stroke death before age 60
			gen strokedeath60 = 0
			replace strokedeath60 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1)
			
			
			/*Interval 60-65*/
			***a. Survival
			*Generate survival time from time 13
			gen survtime60to65 = -ln(U_60to65)/(`lambda_60to65'*exp(`g1_60to65'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death60==0)
			*Generate death indicator for interval 60-65
			gen death60to65 = 0 if (death60==0)
			replace death60to65 = 1 if (survtime60to65 < 5) 
			replace survage = 60 + survtime60to65 if (death60to65==1)
			*Generate indicator for death before age 65
			gen death65 = 0
			replace death65 = 1 if survage < 65

			***b. Stroke
			*Generate stroke time from time 13 exp==0
			gen stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp0_60to65'*exp(`b1'*U)) ///
							if (exp==0 & death60==0 & stroke60==0)
			*Generate stroke time from time 13 exp==1
			replace stroketime60to65 = -ln(U2_60to65)/(`stk_lambda_exp1_60to65'*exp(`b1'*U)) ///
							if (exp==1 & death60==0 & stroke60==0)
			*Generate stroke indicator for interval 60-65
			gen stroke60to65 = 0 if (death60==0 & stroke60==0)
			replace stroke60to65 = 1 if (stroketime60to65 < 5) 
			replace stroke60to65 = 0 if (death60to65==1 & stroketime60to65 != . & stroketime60to65 > survtime60to65)
			*Update prevalent stroke variable (stroke up to age 65)
			replace stroke_history = 1 if (stroke60to65==1)
			replace strokeage = 60 + stroketime60to65 if (stroke60to65==1)
			*Generate indicator for stroke before age 65						
			gen stroke65 = 0
			replace stroke65 = 1 if strokeage < 65

			*Generate stroke deaths for interval 60-65, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath60to65 = runiform()<`pstrokedeath' if stroke60to65==1
			replace death60to65 = 1 if (strokedeath60to65==1) 
			replace survage = strokeage if (strokedeath60to65==1) 
			*Generate indicator for stroke death before age 65
			gen strokedeath65 = 0
			replace strokedeath65 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1)
				
				
			/*Interval 65-70*/
			***a. Survival
			*Generate survival time from time 14
			gen survtime65to70 = -ln(U_65to70)/(`lambda_65to70'*exp(`g1_65to70'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death65==0)
			*Generate death indicator for interval 65-70
			gen death65to70 = 0 if (death65==0)
			replace death65to70 = 1 if (survtime65to70 < 5) 
			replace survage = 65 + survtime65to70 if (death65to70==1)
			*Generate indicator for death before age 70
			gen death70 = 0
			replace death70 = 1 if survage < 70

			***b. Stroke
			*Generate stroke time from time 14 exp==0
			gen stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp0_65to70'*exp(`b1'*U)) ///
							if (exp==0 & death65==0 & stroke65==0)
			*Generate stroke time from time 14 exp==1
			replace stroketime65to70 = -ln(U2_65to70)/(`stk_lambda_exp1_65to70'*exp(`b1'*U)) ///
							if (exp==1 & death65==0 & stroke65==0)
			*Generate stroke indicator for interval 65-70
			gen stroke65to70 = 0 if (death65==0 & stroke65==0)
			replace stroke65to70 = 1 if (stroketime65to70 < 5) 
			replace stroke65to70 = 0 if (death65to70==1 & stroketime65to70 != . & stroketime65to70 > survtime65to70)
			*Update prevalent stroke variable (stroke up to age 70)
			replace stroke_history = 1 if (stroke65to70==1)
			replace strokeage = 65 + stroketime65to70 if (stroke65to70==1)
			*Generate indicator for stroke before age 70						
			gen stroke70 = 0
			replace stroke70 = 1 if strokeage < 70

			*Generate stroke deaths for interval 65-70, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath65to70 = runiform()<`pstrokedeath' if stroke65to70==1
			replace death65to70 = 1 if (strokedeath65to70==1) 
			replace survage = strokeage if (strokedeath65to70==1) 
			*Generate indicator for stroke death before age 70
			gen strokedeath70 = 0
			replace strokedeath70 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1)
				
				
			/*Interval 70-75*/
			***a. Survival
			*Generate survival time from time 15
			gen survtime70to75 = -ln(U_70to75)/(`lambda_70to75'*exp(`g1_70to75'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death70==0)
			*Generate death indicator for interval 70-75
			gen death70to75 = 0 if (death70==0)
			replace death70to75 = 1 if (survtime70to75 < 5) 
			replace survage = 70 + survtime70to75 if (death70to75==1)
			*Generate indicator for death before age 75
			gen death75 = 0
			replace death75 = 1 if survage < 75

			***b. Stroke
			*Generate stroke time from time 15 exp==0
			gen stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp0_70to75'*exp(`b1'*U)) ///
							if (exp==0 & death70==0 & stroke70==0)
			*Generate stroke time from time 15 exp==1
			replace stroketime70to75 = -ln(U2_70to75)/(`stk_lambda_exp1_70to75'*exp(`b1'*U)) ///
							if (exp==1 & death70==0 & stroke70==0)
			*Generate stroke indicator for interval 70-75
			gen stroke70to75 = 0 if (death70==0 & stroke70==0)
			replace stroke70to75 = 1 if (stroketime70to75 < 5) 
			replace stroke70to75 = 0 if (death70to75==1 & stroketime70to75 != . & stroketime70to75 > survtime70to75)
			*Update prevalent stroke variable (stroke up to age 75)
			replace stroke_history = 1 if (stroke70to75==1)
			replace strokeage = 70 + stroketime70to75 if (stroke70to75==1)
			*Generate indicator for stroke before age 75						
			gen stroke75 = 0
			replace stroke75 = 1 if strokeage < 75

			*Generate stroke deaths for interval 70-75, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath70to75 = runiform()<`pstrokedeath' if stroke70to75==1
			replace death70to75 = 1 if (strokedeath70to75==1) 
			replace survage = strokeage if (strokedeath70to75==1) 
			*Generate indicator for stroke death before age 75	
			gen strokedeath75 = 0
			replace strokedeath75 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1)
				
				
			/*Interval 75-80*/
			***a. Survival
			*Generate survival time from time 16
			gen survtime75to80 = -ln(U_75to80)/(`lambda_75to80'*exp(`g1_75to80'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death75==0)
			*Generate death indicator for interval 75-80
			gen death75to80 = 0 if (death75==0)
			replace death75to80 = 1 if (survtime75to80 < 5)
			replace survage = 75 + survtime75to80 if (death75to80==1)
			*Generate indicator for death before age 80
			gen death80 = 0
			replace death80 = 1 if survage < 80

			***b. Stroke
			*Generate stroke time from time 16 exp==0
			gen stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp0_75to80'*exp(`b1'*U)) ///
							if (exp==0 & death75==0 & stroke75==0)
			*Generate stroke time from time 16 exp==1
			replace stroketime75to80 = -ln(U2_75to80)/(`stk_lambda_exp1_75to80'*exp(`b1'*U)) ///
							if (exp==1 & death75==0 & stroke75==0)
			*Generate stroke indicator for interval 75-80
			gen stroke75to80 = 0 if (death75==0 & stroke75==0)	
			replace stroke75to80 = 1 if (stroketime75to80 < 5)
			replace stroke75to80 = 0 if (death75to80==1 & stroketime75to80 != . & stroketime75to80 > survtime75to80)
			*Update prevalent stroke variable (stroke up to age 80)
			replace stroke_history = 1 if (stroke75to80==1)
			replace strokeage = 75 + stroketime75to80 if (stroke75to80==1)
			*Generate indicator for stroke before age 80						
			gen stroke80 = 0
			replace stroke80 = 1 if strokeage < 80

			*Generate stroke deaths for interval 75-80, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath75to80 = runiform()<`pstrokedeath' if stroke75to80==1
			replace death75to80 = 1 if (strokedeath75to80==1) 
			replace survage = strokeage if (strokedeath75to80==1) 
			*Generate indicator for stroke death before age 80	
			gen strokedeath80 = 0
			replace strokedeath80 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1)
			
			
			/*Interval 80-85*/
			***a. Survival
			*Generate survival time from time 17
			gen survtime80to85 = -ln(U_80to85)/(`lambda_80to85'*exp(`g1_80to85'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death80==0)
			*Generate death indicator for interval 80-85
			gen death80to85 = 0 if (death80==0)
			replace death80to85 = 1 if (survtime80to85 < 5)
			replace survage = 80 + survtime80to85 if (death80to85==1)
			*Generate indicator for death before age 85
			gen death85 = 0
			replace death85 = 1 if survage < 85

			***b. Stroke
			*Generate stroke time from time 17 exp==0
			gen stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp0_80to85'*exp(`b1'*U)) ///
							if (exp==0 & death80==0 & stroke80==0)
			*Generate stroke time from time 17 exp==1
			replace stroketime80to85 = -ln(U2_80to85)/(`stk_lambda_exp1_80to85'*exp(`b1'*U)) ///
							if (exp==1 & death80==0 & stroke80==0)
			*Generate stroke indicator for interval 80-85
			gen stroke80to85 = 0 if (death80==0 & stroke80==0)
			replace stroke80to85 = 1 if (stroketime80to85 < 5)
			replace stroke80to85 = 0 if (death80to85==1 & stroketime80to85 != . & stroketime80to85 > survtime80to85)
			*Update prevalent stroke variable (stroke up to age 85)
			replace stroke_history = 1 if (stroke80to85==1)
			replace strokeage = 80 + stroketime80to85 if (stroke80to85==1)
			*Generate indicator for stroke before age 85						
			gen stroke85 = 0
			replace stroke85 = 1 if strokeage < 85

			*Generate stroke deaths for interval 80-85, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath80to85 = runiform()<`pstrokedeath' if stroke80to85==1
			replace death80to85 = 1 if (strokedeath80to85==1) 
			replace survage = strokeage if (strokedeath80to85==1) 
			*Generate indicator for stroke death before age 85	
			gen strokedeath85 = 0
			replace strokedeath85 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1)
			
			
			/*Interval 85-90*/
			***a. Survival
			*Generate survival time from time 18
			gen survtime85to90 = -ln(U_85to90)/(`lambda_85to90'*exp(`g1_85to90'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death85==0)
			*Generate death indicator for interval 85-90
			gen death85to90 = 0 if (death85==0)
			replace death85to90 = 1 if (survtime85to90 < 5)
			replace survage = 85 + survtime85to90 if (death85to90==1)
			*Generate indicator for death before age 90
			gen death90 = 0
			replace death90 = 1 if survage < 90

			***b. Stroke
			*Generate stroke time from time 18 exp==0
			gen stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp0_85to90'*exp(`b1'*U)) ///
							if (exp==0 & death85==0 & stroke85==0)
			*Generate stroke time from time 18 exp==1
			replace stroketime85to90 = -ln(U2_85to90)/(`stk_lambda_exp1_85to90'*exp(`b1'*U)) ///
							if (exp==1 & death85==0 & stroke85==0)
			*Generate stroke indicator for interval 85-90
			gen stroke85to90 = 0 if (death85==0 & stroke85==0)
			replace stroke85to90 = 1 if (stroketime85to90 < 5)
			replace stroke85to90 = 0 if (death85to90==1 & stroketime85to90 != . & stroketime85to90 > survtime85to90)
			*Update prevalent stroke variable (stroke up to age 90)
			replace stroke_history = 1 if (stroke85to90==1)
			replace strokeage = 85 + stroketime85to90 if (stroke85to90==1)
			*Generate indicator for stroke before age 90						
			gen stroke90 = 0
			replace stroke90 = 1 if strokeage < 90

			*Generate stroke deaths for interval 85-90, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath85to90 = runiform()<`pstrokedeath' if stroke85to90==1
			replace death85to90 = 1 if (strokedeath85to90==1) 
			replace survage = strokeage if (strokedeath85to90==1) 
			*Generate indicator for stroke death before age 90	
			gen strokedeath90 = 0
			replace strokedeath90 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1)
			
			
			/*Interval 90-95*/
			***a. Survival
			*Generate survival time from time 19
			gen survtime90to95 = -ln(U_90to95)/(`lambda_90to95'*exp(`g1_90to95'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death90==0)
			*Generate death indicator for interval 90-95
			gen death90to95 = 0 if (death90==0)
			replace death90to95 = 1 if (survtime90to95 < 5)
			replace survage = 90 + survtime90to95 if (death90to95==1)
			*Generate indicator for death before age 95
			gen death95 = 0
			replace death95 = 1 if survage < 95

			***b. Stroke
			/*Interval 90-95*/
			*Generate stroke time from time 19 exp==0
			gen stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp0_90to95'*exp(`b1'*U)) ///
							if (exp==0 & death90==0 & stroke90==0)
			*Generate stroke time from time 19 exp==1
			replace stroketime90to95 = -ln(U2_90to95)/(`stk_lambda_exp1_90to95'*exp(`b1'*U)) ///
							if (exp==1 & death90==0 & stroke90==0)
			*Generate stroke indicator for interval 90-95
			gen stroke90to95 = 0 if (death90==0 & stroke90==0)
			replace stroke90to95 = 1 if (stroketime90to95 < 5)
			replace stroke90to95 = 0 if (death90to95==1 & stroketime90to95 != . & stroketime90to95 > survtime90to95)
			*Update prevalent stroke variable (stroke up to age 95)
			replace stroke_history = 1 if (stroke90to95==1)
			replace strokeage = 90 + stroketime90to95 if (stroke90to95==1)
			*Generate indicator for stroke before age 95						
			gen stroke95 = 0
			replace stroke95 = 1 if strokeage < 95

			*Generate stroke deaths for interval 90-95, and replace death and survage accordingly if stroke death=1																	
			gen strokedeath90to95 = runiform()<`pstrokedeath' if stroke90to95==1
			replace death90to95 = 1 if (strokedeath90to95==1) 
			replace survage = strokeage if (strokedeath90to95==1) 
			*Generate indicator for stroke death before age 95	
			gen strokedeath95 = 0
			replace strokedeath95 = 1 if (strokedeath45to50==1 | strokedeath50to55==1 | strokedeath55to60==1 | ///
				strokedeath60to65==1 | strokedeath65to70==1 | strokedeath70to75==1 | strokedeath75to80==1 | ///
				strokedeath80to85==1 | strokedeath85to90==1 | strokedeath90to95==1)
			
						
			/*Interval 95-100*/
			***a. Survival
			*Generate survival time from time 20
			gen survtime95to100 = -ln(U_95to100)/(`lambda_95to100'*exp(`g1_95to100lt'*exposure +`g2'*U ///
						+`g3'*exposure*U +`g4'*stroke_history +`g5'*0)) ///
							if (death95==0)
			*Generate death indicator for interval 95-100
			gen death95to100 = 0 if (death95==0)
			replace death95to100 = 1 if (survtime95to100 < 5)
			replace survage = 95 + survtime95to100 if (survtime95to100!=.)
			*Generate indicator for death before age 95
			gen death100 = 0
			replace death100 = 1 if survage < 100

			*everyone dies
			gen death = 1

			*top code survage at 100
			*replace survage = 100 if survage > 100

			***b. Stroke: Not generated for ages 95+
				
	
		 stset survtime95to100, failure(death95to100)
		 qui stcox exposure, nohr
		 local lnHR_95to100 = _b[exposure] 
         if `lnHR_95to100' > `target_g1_95to100'+.0001 {
            local toolow=0
            noisily di in red "I stopped at g1_95to100==`g1_95to100lt', ln(HR death 95-100)=`lnHR_95to100'
         }
         else {
            noisily di in white "I did NOT stop at g1_95to100==`g1_95to100lt', ln(HR death 95-100)=`lnHR_95to100'
         }
      }
   }
   local g1_95to100`i' = `g1_95to100lt'
}
clear
set obs 5
gen g1_95to100=.
forvalues i=1/5 {
   replace g1_95to100=`g1_95to100`i'' in `i'
}
sum g1_95to100
global g1_95to100 = `r(mean)'	
global g1_95to100_min = `r(min)'	
global g1_95to100_max = `r(max)'*/		    



/******************************************************************************/
/***		display values lof global variables identified in loops	    ***/
/******************************************************************************/
capture log close
log using global_search_results_ScenarioB2_2016Sept6_corr_per10000PY, replace
*baseline mortality hazard for whites
dis $lambda_45to50 
dis $lambda_50to55 
dis $lambda_55to60 
dis $lambda_60to65 
dis $lambda_65to70 
dis $lambda_70to75 
dis $lambda_75to80 
dis $lambda_80to85 
dis $lambda_85to90 
dis $lambda_90to95 
dis $lambda_95to100

*baseline stroke hazard for whites
dis $stk_lambda_exp0_45to50 
dis $stk_lambda_exp0_50to55 
dis $stk_lambda_exp0_55to60 
dis $stk_lambda_exp0_60to65 
dis $stk_lambda_exp0_65to70 
dis $stk_lambda_exp0_70to75 
dis $stk_lambda_exp0_75to80 
dis $stk_lambda_exp0_80to85 
dis $stk_lambda_exp0_85to90 
dis $stk_lambda_exp0_90to95 

*baseline mortality hazard for blacks
dis $g1_0to1
dis $g1_1to5 
dis $g1_5to10 
dis $g1_10to15 
dis $g1_15to20 
dis $g1_20to25 
dis $g1_25to30 
dis $g1_30to35 
dis $g1_35to40 
dis $g1_40to45
dis $g1_45to50 
dis $g1_50to55 
dis $g1_55to60 
dis $g1_60to65 
dis $g1_65to70 
dis $g1_70to75 
dis $g1_75to80 
dis $g1_80to85 
dis $g1_85to90 
dis $g1_90to95 
dis $g1_95to100


/*mins*/
*baseline mortality hazard for whites
dis $lambda_45to50_min 
dis $lambda_50to55_min  
dis $lambda_55to60_min  
dis $lambda_60to65_min  
dis $lambda_65to70_min  
dis $lambda_70to75_min  
dis $lambda_75to80_min  
dis $lambda_80to85_min  
dis $lambda_85to90_min  
dis $lambda_90to95_min  
dis $lambda_95to100_min 

*baseline stroke hazard for whites
dis $stk_lambda_exp0_45to50_min  
dis $stk_lambda_exp0_50to55_min  
dis $stk_lambda_exp0_55to60_min  
dis $stk_lambda_exp0_60to65_min  
dis $stk_lambda_exp0_65to70_min  
dis $stk_lambda_exp0_70to75_min  
dis $stk_lambda_exp0_75to80_min  
dis $stk_lambda_exp0_80to85_min  
dis $stk_lambda_exp0_85to90_min  
dis $stk_lambda_exp0_90to95_min  

*baseline mortality hazard for blacks
dis $g1_0to1_min
dis $g1_1to5_min 
dis $g1_5to10_min 
dis $g1_10to15_min 
dis $g1_15to20_min 
dis $g1_20to25_min 
dis $g1_25to30_min 
dis $g1_30to35_min 
dis $g1_35to40_min 
dis $g1_40to45_min
dis $g1_45to50_min
dis $g1_50to55_min 
dis $g1_55to60_min 
dis $g1_60to65_min 
dis $g1_65to70_min 
dis $g1_70to75_min 
dis $g1_75to80_min 
dis $g1_80to85_min 
dis $g1_85to90_min 
dis $g1_90to95_min 
dis $g1_95to100_min


/*starting values*/
*baseline mortality hazard for whites
dis $lambda_45to50l 
dis $lambda_50to55l 
dis $lambda_55to60l 
dis $lambda_60to65l 
dis $lambda_65to70l  
dis $lambda_70to75l 
dis $lambda_75to80l 
dis $lambda_80to85l 
dis $lambda_85to90l 
dis $lambda_90to95l 
*dis $lambda_95to100l 

*baseline stroke hazard for whites
dis $stk_lambda_exp0_45to50l  
dis $stk_lambda_exp0_50to55l 
dis $stk_lambda_exp0_55to60l 
dis $stk_lambda_exp0_60to65l 
dis $stk_lambda_exp0_65to70l 
dis $stk_lambda_exp0_70to75l 
dis $stk_lambda_exp0_75to80l 
dis $stk_lambda_exp0_80to85l 
dis $stk_lambda_exp0_85to90l 
dis $stk_lambda_exp0_90to95l 

*baseline mortality hazard for blacks
dis $g1_0to1l 
dis $g1_1to5l 
dis $g1_5to10l
dis $g1_10to15l 
dis $g1_15to20l
dis $g1_20to25l 
dis $g1_25to30l
dis $g1_30to35l 
dis $g1_35to40l
dis $g1_40to45l 
dis $g1_45to50l 
dis $g1_50to55l 
dis $g1_55to60l 
dis $g1_60to65l 
dis $g1_65to70l 
dis $g1_70to75l 
dis $g1_75to80l 
dis $g1_80to85l 
dis $g1_85to90l 
dis $g1_90to95l 
*dis $g1_95to100


capture log close

 timer off 1
 
 timer list 1

