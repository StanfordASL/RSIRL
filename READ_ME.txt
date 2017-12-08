########## RISK-SENSITIVE INVERSE REINFORCEMENT LEARNING VIA SEMI- AND NON-PARAMETRIC METHODS ###########

The directory contains the scripts and data to reproduce the results presented in the paper:
- Learn the risk-envelope of participants from the driving-simulation game, for single-stage or multi-stage decision problems
- Learn the risk-neutral of participants from the driving-simulation game, for single-stage or multi-stage decision problems
- Compute the prediction errors on the test set
- Visualize the results of some participants that present typical behaviors (highly risk-averse or risk-seeking)

##############################################################
########## USAGE
##############################################################
1 - Risk-sensitive model learning: RSIRL.m
	Name of participant: set name_str = 'some_participant' ('p2' by default, risk-averse participant)
	Number of decision stages: set N=1 (single decision stage) or N=2 (two decision stages). By default, N=2.
	Choose parameters of gradient ascent: T (number of time steps), w_initial (initial value of offset parameters) and wc_initial (initial cost weights) 
	Run RSIRL.m. Results are saved in '/Data/Inference/results_N_2_p2.mat'

2 - Risk-neutral model learning: RSIRL_neutral.m
	Name of participant: set name_str = 'some_participant' ('p2' by default, risk-averse participant)
	Number of decision stages: set N=1 (single decision stage) or N=2 (two decision stages). By default, N=2.
	Choose parameters of gradient ascent: T (number of time steps), w_initial (initial value of r) and wc_initial (initial cost weights) 
	Run RSIRL_neutral.m. Results are saved in '/Data/Inference/results_neutral_N_2_p2.mat'

3 - Evaluate performance of risk-sensitive and risk-neutral models on a test set: compute_errors.m
	Name of participant: set name_str = 'some_participant' ('p2' by default, risk-averse participant)
	Number of decision stages: set N=1 (single decision stage) or N=2 (two decision stages). By default, N=2.
	Run compute_errors.m. Prediction errors (for both risk-sensitive and risk-neutral models) are recorded in '/Data/Inference/Performance/perf_N_2_p1.mat'.
	Run analysis.m. Bar plots of errors (as presented in the paper) are recorded in '/Data/Inference/Performance/Plots/Bar_expected/...' and '/Data/Inference/Performance/Plots/Bar_most_likely/...'

4 - Analysis of results of specific participants
	Case-study #1: run /Data/Inference/Performance/rn_participant_analysis.m
	Case-study #2: run /Data/Inference/Performance/rs_participant_analysis.m
	Case-study #3: run /Data/Inference/Performance/rn_participant_analysis.m

##############################################################
### Data
##############################################################
Data recorded from the driving-simulation experiments are located in '/Data'. Data are already split into training and test sets.
For each participant (p1 to p10), the directory contains:
- 'robot_data_10Hz_p1_train.mat': the states and actions of the robot (at 10Hz) in participant p1's experiment.
- 'robot_data_10Hz_p1_test.mat': same as previously but for testing.
- 'human_robot_data_p1_train.mat': states and actions (for training) recorded by the simulator (at 60Hz) for both the participant (p1) and the robot. 
NB: We record twice the data of the robot for the following reasons: for learning purposes, we prefer to use 'robot_data_10Hz_p1_train.mat' since the robot's trajectory generated prior to the experiments (and hence more precise); then we need to align generated robot's data and human's recorded data (by tuning the variable k_delay)
- 'human_robot_data_p1_test.mat': same as previously but for testing.

##############################################################
########## Learning
##############################################################
RSIRL.m (risk-sensitive model) run the semi-parametric algorithm of the paper. 
The gradient-based algorithm is implemented in subgradientdescent.m


##############################################################
########## Dependencies
##############################################################
- Matlab
- Multi-Parametric Toolbox 3.0
- Mosek solver
- Yalmip 











