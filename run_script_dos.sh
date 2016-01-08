#!/bin/bash

: '
Calculates the density of states (DOS) in the nanowire of a superconductor-nanowire-superconductor (SNS) junction. Length of the normal metal part is fixed to L = 1 and is chosen to be along z-axis. There is an exchange field in the xz-plane of the form cos(x) sigma_z + sin(x) sigma_x where x is the angle between the nanowire and the exchange field. The generic spin-orbit term is A = (alpha1 * sigma_x + alpha2 * sigma_z, 0, alpha3 * sigma_x + alpha4 * sigma_z). 

User specified parameters:

delta 			: superconducting order parameter in the units of Thouless energy
phase			: phase difference between the superconductors in the units of pi
alpha_i			: magnitude of the Rashba coupling constant with
				  A = (alpha_1 sigma_x + alpha_2 sigma_z, 0, alpha_3 sigma_x + alpha_4 sigma_z)
h				: magnitude of the exchange field in the units of Thouless energy
h_angle			: angle of the exchange field with h = cos(x) sigma_z + sin(x) sigma_x

phase, alpha_i, h, and h_angle are given as ranges so multiple values can be calculated in a single run. Ranges are of floats and should be given in 2-decimal precision.

User defined directories:

FLOAT_RANGE_DIR	: location of the frange.py function that generates float ranges
SOL_DIR			: directory there the solutions are stored
RES_DIR			: directory where the results are stored

By default the subdirectories for solutions and results are generated in the directory where the script was run and it is assumed that frange.py function is in the same directory as the run_script.

Technical parameters:

max_energy		: the maximum energy for which the DOS is calculated
energy_res		: the energy resolution used in the calculation i.e. the energy spacing
tolerance		: tolerance for the approximate solution
overwrite		: "Yes" if user wants to overwrite existing solutions.
				  "No" if user wants to skip existing solutions.

The solver uses previous solutions as initial guesses. The previous solution for the initial guess is chosen based on the spacing parameters. If previous results do not exist, a simple guess is used. To quarantee convergence, it is best to increase the parameters gradually so that previous solutions provide a good enough initial guess for the next calculation.
'

delta=1000

# Define all the ranges with 2-decimal precision! This makes sure that the solver finds previous solutions and
# that the plot_script will find the results for plotting. Do not use zero as any of the spacings. If you want
# to calculate only one value specify the same start and stop instead.

# The phase difference range (in the units of pi)

phase_start=0.00
phase_stop=0.00
phase_spacing=0.25

# The Rashba coupling strength range (in the units of one over length)

alpha1_start=0.00
alpha1_stop=0.00
alpha1_spacing=0.50

alpha2_start=0.00
alpha2_stop=0.00
alpha2_spacing=0.50

alpha3_start=0.00
alpha3_stop=2.00
alpha3_spacing=0.10

alpha4_start=0.00
alpha4_stop=0.00
alpha4_spacing=0.50

# The exchange field magnitude range (in the units of Thouless energy)

h_start=0.00
h_stop=24.00
h_spacing=8.00

# The exchange field angle range (in degrees)

h_angle_start=0.00
h_angle_stop=0.00
h_angle_spacing=15.00

# The maximum energy used in the calculation.
max_energy=30

# The energy resolution used in the calculation.
energy_res=0.50

# Tolerance for the approximate solution.
tolerance=1.0e-6

# If "Yes", overwrite the solutions that already exists. If "No", the parameters for which the results already exist are skipped.
overwrite="No"


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# User defined directories. You do not need to change #
# these if you're happy with the default directories. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# This points to the function that generates float ranges. You do not need to change this if you do not move the frange.py or the run scripts.
FLOAT_RANGE_DIR="src/frange.py"

# Current directory to be used for the solution and result directories.
CURRENT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# The directory where the solutions are stored. The script generates subdirectories for the results. Give here the root where you want the directories to be generated. By default the subdirectories are generated in the folder where the script was run.
SOL_DIR="$CURRENT_DIR/DOS_Solutions/"

# The directory where the results are stored. The script generates subdirectories for the results. Give here the root where you want the directories to be generated. Make sure to use the same directory in the plot_script. By default the subdirectories are generated in the folder where the script was run.
RES_DIR="$CURRENT_DIR/DOS_Results/"

# Note that a new subdirectory for different delta is not automatically generated. If you change delta, make sure to change the SOL_DIR and RES_DIR if you want to keep the old results. You want to overwrite instead, change overwrite="Yes" and do not change the names of the directories. You can choose for example SOL_DIR="$CURRENT_DIR/Solutions/delta_$delta/" and RES_DIR="$CURRENT_DIR/Results/delta_$delta" if needed. These subdirectories are not automatically generated to avoid an excessive amount of folders.


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is the part that runs the python code and generates the missing  #
# directories for different parameter ranges. Do not change this!       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

for phase in `python $FLOAT_RANGE_DIR $phase_start $phase_stop $phase_spacing`; do		
	echo "phase difference =" $phase
	for alpha_1 in `python $FLOAT_RANGE_DIR $alpha1_start $alpha1_stop $alpha1_spacing`; do
	    echo "alpha_1 =" $alpha_1
    	for alpha_2 in `python $FLOAT_RANGE_DIR $alpha2_start $alpha2_stop $alpha2_spacing`; do
		    echo "alpha_2 =" $alpha_2
			for alpha_3 in `python $FLOAT_RANGE_DIR $alpha3_start $alpha3_stop $alpha3_spacing`; do
				echo "alpha_3 =" $alpha_3
				for alpha_4 in `python $FLOAT_RANGE_DIR $alpha4_start $alpha4_stop $alpha4_spacing`; do
					echo "alpha_4 =" $alpha_4
					for h_angle in `python $FLOAT_RANGE_DIR $h_angle_start $h_angle_stop $h_angle_spacing`; do
						echo "h_angle =" $h_angle
						for h in `python $FLOAT_RANGE_DIR $h_start $h_stop $h_spacing`; do
							echo "h =" $h
							if [ ! -d "$SOL_DIR/Phase_difference_$phase/Rashba_${alpha_1}_${alpha_2}_${alpha_3}_${alpha_4}/h_${h}_angle_${h_angle}/" ]; then
								mkdir -p $SOL_DIR/Phase_difference_$phase/Rashba_${alpha_1}_${alpha_2}_${alpha_3}_${alpha_4}/h_${h}_angle_${h_angle}/
							fi
							if [ ! -d "$RES_DIR/Phase_difference_$phase/Rashba_${alpha_1}_${alpha_2}_${alpha_3}_${alpha_4}/" ]; then
								mkdir -p $RES_DIR/Phase_difference_$phase/Rashba_${alpha_1}_${alpha_2}_${alpha_3}_${alpha_4}/
							fi
							if [ "$overwrite" = "No" ]; then
								if [ `ls -1 "$SOL_DIR/Phase_difference_$phase/Rashba_${alpha_1}_${alpha_2}_${alpha_3}_${alpha_4}/h_${h}_angle_${h_angle}/" |wc -l` != $(bc <<< "$max_energy / $energy_res + 1") ]; then
									python src/DOS.py $phase $alpha_1 $alpha_2 $alpha_3 $alpha_4 $h $h_angle $delta $phase_spacing $alpha1_spacing $alpha2_spacing $alpha3_spacing $alpha4_spacing $h_spacing $h_angle_spacing $SOL_DIR $RES_DIR $max_energy $energy_res $tolerance
								fi
							elif [ "$overwrite" = "Yes" ]; then
								python src/DOS.py $phase $alpha_1 $alpha_2 $alpha_3 $alpha_4 $h $h_angle $delta $phase_spacing $alpha1_spacing $alpha2_spacing $alpha3_spacing $alpha4_spacing $h_spacing $h_angle_spacing $SOL_DIR $RES_DIR $max_energy $energy_res $tolerance
							else
								echo "Please specify the overwrite parameter."
								exit
							fi
						done
		   	        done
				done
			done
		done
	done
done
