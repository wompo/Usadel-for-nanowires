from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import subprocess

# # # # # # # # # # # # #
# Function definitions. #
# # # # # # # # # # # # #

# Range function for floats.
def frange(start, end=None, inc=None):
    "A range function, that does accept floats"

    if end == None:
        end = start + 0.0
        start = 0.0
    else: start += 0.0 # force it to be a float

    if inc == None:
        inc = 1.0

    count = int((end - start) / inc)
    if start + count * inc != end + inc:
        count += 1

    L = [None,] * count
    for i in xrange(count):
        L[i] = start + i * inc

    return L

# Write text generated into the figure. Used in add_text.
def write_text(loc_x, loc_y, text, ha, va):
	return plt.text(loc_x, loc_y, text,
					ha=ha,
					va=va,
					transform = ax.transAxes,
					fontsize = 18)

# Add all other parameters than the primary and secondary variable into the text.
def add_text(text_location):
	if text_location == "center left":
		loc_x = 0.03
		loc_y = 0.56
		ha='left'
		va='center'
	elif text_location == "upper right":
		loc_x = 0.97
		loc_y = 0.94
		ha='right'
		va='center'
	elif text_location == "lower right":
		loc_x = 0.97
		loc_y = 0.18
		ha = 'right'
		va = 'center'
	elif text_location == "lower left":
		loc_x = 0.03
		loc_y = 0.18
		ha = 'left'
		va = 'center'
	elif text_location == "center right":
		loc_x = 0.97
		loc_y = 0.56
		ha = 'right'
		va = 'center'
	elif text_location == "upper left":
		loc_x = 0.03
		loc_y = 0.94
		ha = 'left'
		va = 'center'
	elif text_location == "upper center":
		loc_x = 0.50
		loc_y = 0.94
		ha = 'center'
		va = 'center'
	elif text_location == "lower center":
		loc_x = 0.50
		loc_y = 0.18
		ha = 'center'
		va = 'center'
	elif text_location != "none":
		sys.exit("Please set the text_location parameter.")
	else:
		return

	parameter_dict = {"phase": r'$\phi = ' + "{0:.2f}".format(phase) + '\pi$', "T": '$T = ' + "{0:.2f}".format(T) + 'E_T$', "delta": r'$\Delta = ' + "{0:.0f}".format(delta) + ' E_T$',  "exchange": '$h = ' + "{0:.2f}".format(exchange) + ' E_T$', "exchange_angle": r'$\theta = ' + "{0:.2f}".format(exchange_angle) + '^\circ$', "alpha_1": "{0:.2f}".format(alpha_1), "alpha_2": "{0:.2f}".format(alpha_2), "alpha_3": "{0:.2f}".format(alpha_3), "alpha_4": "{0:.2f}".format(alpha_4)}

	del parameter_dict[primary_variable]
	del parameter_dict[secondary_variable]

	line1, line2, line3 = ("" for i in range(3))

	for parameter in ["phase", "T", "delta"]:
		if parameter in parameter_dict:
			if len(line1) > 0:
				line1 += '$,$ '
			line1 += parameter_dict[parameter]

	if suppress_zeros == "alphas" or suppress_zeros == "both":
		non_zero = []
		for parameter in ["alpha_1", "alpha_2", "alpha_3", "alpha_4"]:
			if parameter in parameter_dict:
				if float(parameter_dict[parameter]) != 0.00:
					non_zero.append(parameter)
		parameters = non_zero
	elif suppress_zeros != "False" and suppress_zeros != "h_angle":
		sys.exit("Please set suppress_zeros parameter to 'alphas', 'h_angle', 'both' or 'False'.")
	else:
		parameters = ["alpha_1", "alpha_2", "alpha_3", "alpha_4"]

	if len(parameters) == 1:
		line2 += r'$' + '\\' + non_zero[0] + ' = ' + parameter_dict[non_zero[0]] + '/L$'
	else:
		line2 += '$('
		for parameter in parameters:
			if parameter in parameter_dict:
				if len(line2) > 3:
					line2 += ', '
				line2 += '\\' + parameter
		line2 += ') = ('
		k = 0
		for parameter in parameters:
			if parameter in parameter_dict:
				if k > 0:
					line2 += ', '
				line2 += parameter_dict[parameter]
			k += 1
		line2 += ') /L$'

	if suppress_zeros == "h_angle" or suppress_zeros == "both":
		for parameter in ["exchange"]:
			if parameter in parameter_dict:
				line3 += parameter_dict[parameter]
	else:
		for parameter in ["exchange", "exchange_angle"]:
			if parameter in parameter_dict:
				if len(line3) > 0:
					line3 += '$,$ '
				line3 += parameter_dict[parameter]

	if text_location == "lower left" or text_location == "lower center" or text_location == "lower right":
		if len(line2) == 0:
			write_text(loc_x, loc_y - 0.06, line1, ha, va)
			write_text(loc_x, loc_y - 0.12, line3, ha, va)
		elif len(line3) == 0:
			write_text(loc_x, loc_y - 0.06, line1, ha, va)
			write_text(loc_x, loc_y - 0.12, line2, ha, va)
		else:
			write_text(loc_x, loc_y, line1, ha, va)
			write_text(loc_x, loc_y - 0.06, line2, ha, va)
			write_text(loc_x, loc_y - 0.12, line3, ha, va)
	elif text_location == "center left" or text_location == "center right":
		if len(line2) == 0:
			write_text(loc_x, loc_y - 0.03, line1, ha, va)
			write_text(loc_x, loc_y - 0.09, line3, ha, va)
		elif len(line3) == 0:
			write_text(loc_x, loc_y - 0.03, line1, ha, va)
			write_text(loc_x, loc_y - 0.09, line2, ha, va)
		else:
			write_text(loc_x, loc_y, line1, ha, va)
			write_text(loc_x, loc_y - 0.06, line2, ha, va)
			write_text(loc_x, loc_y - 0.12, line3, ha, va)
	else:
		if len(line2) == 0:
			write_text(loc_x, loc_y, line1, ha, va)
			write_text(loc_x, loc_y - 0.06, line3, ha, va)
		elif len(line3) == 0:
			write_text(loc_x, loc_y, line1, ha, va)
			write_text(loc_x, loc_y - 0.06, line2, ha, va)
		else:
			write_text(loc_x, loc_y, line1, ha, va)
			write_text(loc_x, loc_y - 0.06, line2, ha, va)
			write_text(loc_x, loc_y - 0.12, line3, ha, va)

# # # # # # # # # # # # # # # 
# Function definitions end. #
# # # # # # # # # # # # # # #

L = 1  # Lenght of the nanowire.

# Get all the parameters from the run script.

phase	 	   = float(sys.argv[1])
alpha_1 	   = float(sys.argv[2])
alpha_2 	   = float(sys.argv[3])
alpha_3 	   = float(sys.argv[4])
alpha_4 	   = float(sys.argv[5])
exchange 	   = float(sys.argv[6])
exchange_angle = float(sys.argv[7])
T		 	   = float(sys.argv[8])
delta		   = float(sys.argv[9])

primary_variable		= sys.argv[10]
primary_range_start 	= float(sys.argv[11])
primary_range_stop 		= float(sys.argv[12])
primary_range_spacing 	= float(sys.argv[13])

secondary_variable		= sys.argv[14]
secondary_range_start 	= float(sys.argv[15])
secondary_range_stop 	= float(sys.argv[16])
secondary_range_spacing = float(sys.argv[17])

legend_location	= sys.argv[18]
text_location	= sys.argv[19]
suppress_zeros	= sys.argv[20]
dashed_line		= sys.argv[21]

res_dir			= sys.argv[22]
fig_dir			= sys.argv[23]

primary_range	= frange(primary_range_start, primary_range_stop, primary_range_spacing)
secondary_range	= frange(secondary_range_start, secondary_range_stop, secondary_range_spacing)


# Generate data array from the results.

data = [[] for i in range(len(secondary_range))]

l    = 0
for j in secondary_range:
	exec("%s = %.2f" % (secondary_variable, j))
	for i in primary_range:
		exec("%s = %.2f" % (primary_variable, i))
		file = res_dir + 'Phase_difference_' + "{0:.2f}".format(phase) + '/T_' + "{0:.2f}".format(T) + '/Rashba_' + "{0:.2f}".format(alpha_1) + '_' + "{0:.2f}".format(alpha_2) + '_' + "{0:.2f}".format(alpha_3) + '_' + "{0:.2f}".format(alpha_4) + '/h-' + "{0:.2f}".format(exchange) + '-angle-' + "{0:.2f}".format(exchange_angle) + '-supercurrent.txt'
		subprocess.call("sed -n 5p " + file + " | awk '{print $1}' > temp.dat", shell=True)
		data[l].append(np.loadtxt('temp.dat', delimiter=' '))
		subprocess.call("rm temp.dat", shell=True)
	l += 1

# Plot the results with proper labeling.

colors = ['b', 'g', 'r', 'c', 'm', '#FF9900', '#66FF33', 'k', '#FFFF00', '#996633', '#FF99FF']

fig = plt.figure()

for i in range(len(secondary_range)):
	if secondary_variable == "alpha_1":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label=r'$\alpha_1 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "alpha_2":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label=r'$\alpha_2 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "alpha_3":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label=r'$\alpha_3 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "alpha_4":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label=r'$\alpha_4 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "exchange":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label='$h = ' + '{0:.2f}'.format(secondary_range[i]) + ' E_T$')
	elif secondary_variable == "phase":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label='$\phi = ' + '{0:.2f}'.format(secondary_range[i]) + ' \pi$')
	elif secondary_variable == "exchange_angle":
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label=r'$\theta = ' + '{0:.2f}'.format(secondary_range[i]) + '^\circ$')
	else:
		plt.plot(primary_range, data[i], '-', linewidth=1.5, color = colors[i], label='$T = ' + '{0:.2f}'.format(secondary_range[i]) + ' E_T$')

if primary_variable == "alpha_1":
	plt.xlabel(r'$\alpha_1 / L$', fontsize = 19)
elif primary_variable == "alpha_2":
	plt.xlabel(r'$\alpha_2 / L$', fontsize = 19)
elif primary_variable == "alpha_3":
	plt.xlabel(r'$\alpha_3 / L$', fontsize = 19)
elif primary_variable == "alpha_4":
	plt.xlabel(r'$\alpha_4 / L$', fontsize = 19)
elif primary_variable == "exchange":
	plt.xlabel(r'$h E_T$', fontsize = 19)
elif primary_variable == "phase":
	plt.xlabel(r'$\phi \pi$', fontsize = 19)
elif primary_variable == "exchange_angle":
	plt.xlabel(r'$\theta$', fontsize = 19)
else:
	plt.xlabel(r'$T E_T$', fontsize = 19)

plt.ylabel(r'$I_S e R_N / E_T$', fontsize = 19)

# Prints horizontal dashed line to at y = 0

if dashed_line == "True":
	plt.axhline(y=0, xmin=primary_range[0], xmax=primary_range[-1], color='k', linestyle='--')

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

# This adds a text stating the rest of the variables used in the calculation.

ax = fig.add_subplot(111)
add_text(text_location)

# Legend box containing labels.

if legend_location != "none":
	plt.legend(loc=legend_location, handlelength=1.6, handletextpad=0.4, fontsize = 16)

# Finally save the figure into the subdirectory specified in the plot script.

fig.savefig(fig_dir + 'Figure.png')
