#!/bin/bash

: '
Plots supercurrent as a function of the given primary variable. Secondary variable adds multiple curves to a single figure.

Valid arguments for primary_variable and secondary_variable are:

T				: temperature in the units of Thouless energy
phase			: phase difference between the superconductors in the units of pi
alpha_i			: magnitude of the Rashba coupling constant with 
				  A = (alpha_1 sigma_x + alpha_2 sigma_z, 0, alpha_3 sigma_x + alpha_4 sigma_z)
exchange		: magnitude of the exchange field in the units of Thouless energy
exchange_angle	: angle of the exchange field with h = cos(x) sigma_z + sin(x) sigma_x

Note that even if you defined "exchange" as a primary_variable, you must give a value also for the variable "exchange" below. This does not effect the actual figure, but is required for all the parameters for technical reasons. This is true for all primary and secondary variables.

User must also define:

delta		 	: superconducting order parameter
legend_location	: location of the legend in the plot ("none" is also a valid option)
text_location	: location of the text where all the fixed variables are listed ("none" is also a valid option)
suppress_zeros	: do not write down some fixed parameters wiht zero value into to figure
dashed_line		: Draws a dashed line along y=0 if True and does not draw it if "False".

User defined directories:

RES_DIR			: directory where the results are stored (make sure this is the same as in the run script)
FIG_DIR 		: location where the created figure is stored
'

primary_variable="exchange"

primary_range_start=0.00
primary_range_stop=80.00
primary_range_spacing=5.00

secondary_variable="alpha_3"

secondary_range_start=0.00
secondary_range_stop=2.00
secondary_range_spacing=1.00

# Define the rest of the variables required for the plot. Give also arbitrary values for the ones you chose as primary and secondary variables. Give as a float with 2-decimals.

phase=0.50
alpha_1=0.00
alpha_2=0.00
alpha_3=2.00
alpha_4=0.00
exchange=2.00
exchange_angle=0.00
T=0.10
delta=1000


# # # # # # # # # # # # # # # # # # # # # # #
# The following parameters are used to		#
# modify the look of the outcoming figure.  #
# # # # # # # # # # # # # # # # # # # # # # #

# Legend location. Valid locations are: right, center left, upper right, lower right, best, center, lower left, center right, upper left, upper center, lower center.
legend_location="best"

# This adds a text stating the rest of the fixed variables used in the calculation to the plot. Valid locations are: center left, upper right, lower right, lower left, center right, upper left, upper center, lower center, and none for no text.
text_location="center right"

# If "alphas", writes only the non-zero alphas into the text in the figure. If "h_angle", do not write the exchange field angle into the figure if it is zero. If "both", both of the previous apply. If "False", write all parameters regardless the value.
suppress_zeros="both"

# Draws a dashed line along y=0 if "True" and does not draw it if "False".
dashed_line="True"


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# User defined directories. You do not need to change #
# these if you're happy with the default directories. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# The directory where the results are stored. Make sure this is the same one used in the run script.
RES_DIR="Test/Results/"

# The directory where the figures are saved.
FIG_DIR="Test/Figures/"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is the part that runs the python code and generates the missing  #
# directories for different parameter ranges. Do not change this!       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

CURRENT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
PARENT_DIR="$(dirname "$CURRENT_DIR")"

if [ ! -d "$FIG_DIR" ]; then
	mkdir -p $FIG_DIR
fi

python $PARENT_DIR/src/plot_supercurrent.py $phase $alpha_1 $alpha_2 $alpha_3 $alpha_4 $exchange $exchange_angle $T $delta $primary_variable $primary_range_start $primary_range_stop $primary_range_spacing $secondary_variable $secondary_range_start $secondary_range_stop $secondary_range_spacing "$legend_location" "$text_location" $suppress_zeros $dashed_line $RES_DIR $FIG_DIR
