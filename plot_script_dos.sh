#!/bin/bash

: '
Plots the density of states as a function of energy. Secondary variable adds multiple curves to a single figure.

max_energy		: absolute value of the maximum energy displayed in the x-axis in the plot
energy_res		: density of the data points in the plot

Valid arguments for the secondary_variable are:

phase			: phase difference between the superconductors in the units of pi
alpha_i			: magnitude of the Rashba coupling constant with
				  A = (alpha_1 sigma_x + alpha_2 sigma_z, 0, alpha_3 sigma_x + alpha_4 sigma_z)
exchange		: magnitude of the exchange field in the units of Thouless energy
exchange_angle	: angle of the exchange field with h = cos(x) sigma_z + sin(x) sigma_x

Note that even if you defined "exchange" as a primary_variable, you must give a value also for the variable "exchange" below. This does not effect the actual figure, but is required for all the parameters for technical reasons. This is true for all primary and secondary variables.

User must also define:

delta 			: superconducting order parameter
legend_location	: location of the legend in the plot ("none" is also a valid option)
text_location	: location of the text where all the fixed variables are listed ("none" is also a valid option)
suppress_zeros	: do not write down some fixed parameters with zero value into to figure
dashed_line		: Draws a dashed line along y=1 if "True" and does not draw it if "False".

RES_DIR			: directory there the results are stored (make sure this is the same as in the run script dos)
FIG_DIR 		: directory where the created figure is stored
'

# Note that you must have the results to be able to plot them. Remember to define a range for which the results exist.

max_energy=30.00
energy_res=0.50		# Use the same energy_res as in the run_script_dos.

secondary_variable="exchange"

secondary_range_start=0.00
secondary_range_stop=8.00
secondary_range_spacing=8.00

# Define the rest of the variables required for the plot. Give also arbitrary values for the ones you chose as primary and secondary variables. Give as floats with 2-decimals.

phase=0.00
alpha_1=0.00
alpha_2=0.00
alpha_3=0.00
alpha_4=0.00
exchange=8.00
exchange_angle=0.00
delta=1000


# # # # # # # # # # # # # # # # # # # # # # #
# The following parameters are used to		#
# modify the look of the outcoming figure.  #
# # # # # # # # # # # # # # # # # # # # # # #

# Legend location. Valid locations are: right, center left, upper right, lower right, best, center, lower left, center right, upper left, upper center, lower center, and none for no legend.
legend_location="best"

# This adds a text stating the rest of the fixed variables used in the calculation to the plot. Valid locations are: center left, upper right, lower right, lower left, center right, upper left, upper center, lower center, and none for no text.
text_location="lower left"

# If "alphas", writes only the non-zero alphas into the text in the figure. If "h_angle", do not write the exchange field angle into the figure if it is zero. If "both", both of the previous apply. If "False", write all parameters regardless the value.
suppress_zeros="h_angle"

# Draws a dashed line along y=1 if "True" and does not draw it if "False".
dashed_line="False"

# The minimum and maximum of the y-axis in the figure. Give y_min=0.00 and y_max=0.00 to use dynamical scaling on the y-axis.
y_min=0.00
y_max=1.80


# # # # # # # # # # # # # # # # # # # # # # # # # # # #
# User defined directories. You do not need to change #
# these if you're happy with the default directories. #
# # # # # # # # # # # # # # # # # # # # # # # # # # # #

# The directory where the results are stored. Make sure this is the same one used in the run script.
RES_DIR="DOS_results/"

# The directory where the figures are saved
FIG_DIR="DOS_figures/"


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# This is the part that runs the python code and generates the missing  #
# directories for different parameter ranges. Do not change this!       #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if [ ! -d "$FIG_DIR" ]; then
	mkdir -p $FIG_DIR
fi

python src/plot_DOS.py $phase $alpha_1 $alpha_2 $alpha_3 $alpha_4 $exchange $exchange_angle $delta $max_energy $energy_res $secondary_variable $secondary_range_start $secondary_range_stop $secondary_range_spacing "$legend_location" "$text_location" $suppress_zeros $dashed_line $y_min $y_max $RES_DIR $FIG_DIR
