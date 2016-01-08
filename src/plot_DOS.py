from __future__ import print_function
import sys
import numpy as np
import matplotlib.pyplot as plt
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

	parameter_dict = {"phase": r'$\phi = ' +  "{0:.2f}".format(phase) + '\pi$', "delta": r'$\Delta = ' + "{0:.0f}".format(delta) + ' E_T$',  "exchange": '$h = ' + "{0:.2f}".format(exchange) + ' E_T$', "exchange_angle": r'$\theta = ' + "{0:.2f}".format(exchange_angle) + '^\circ$', "alpha_1": "{0:.2f}".format(alpha_1), "alpha_2": "{0:.2f}".format(alpha_2), "alpha_3": "{0:.2f}".format(alpha_3), "alpha_4": "{0:.2f}".format(alpha_4)}

	del parameter_dict[secondary_variable]

	line1, line2, line3 = ("" for i in range(3))

	for parameter in ["phase", "delta"]:
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

# Gets data from the result file and append it to the dos array cut to proper size.

def get_data(file, pos, max_energy, energy_res):
	z = pos*1.0/(45 - 1)
	next_z = (pos + 1)*1.0/(45 - 1)

	get_first_line = subprocess.Popen("awk '/z = " + "{0:.2f}".format(z) + "/{ print NR; exit }' " + file, shell=True, stdout=subprocess.PIPE)
	first_line = int(get_first_line.stdout.read()) + 1
	get_last_line = subprocess.Popen("awk '/z = " + "{0:.2f}".format(next_z) + "/{ print NR; exit }' " + file, shell=True, stdout=subprocess.PIPE)
	last_line = int(get_last_line.stdout.read()) - 1

	subprocess.call("sed -n -e " + str(first_line) + "," + str(last_line) + "p " + file + " | awk '{print $1}' > temp.dat", shell=True)
	data = np.genfromtxt('temp.dat', delimiter=' ', usecols = 0, unpack=True)
	data_reversed = data[::-1]
	dos = np.concatenate((data_reversed[:-1] ,data), axis=0)
	
	subprocess.call("rm temp.dat", shell=True)

	cutoff = int(max_energy/energy_res)
	middle_index = len(dos)/2
	start = middle_index - cutoff
	stop = middle_index + cutoff + 1
	dos_cut.append(dos[start:stop])

# # # # # # # # # # # # # # # 
# Function definitions end. #
# # # # # # # # # # # # # # #

L = 1 # Lenght of the nanowire.

# Get all the parameters from the run script.

phase	 	   	= float(sys.argv[1])
alpha_1 	   	= float(sys.argv[2])
alpha_2 	   	= float(sys.argv[3])
alpha_3 	   	= float(sys.argv[4])
alpha_4 	   	= float(sys.argv[5])
exchange 	   	= float(sys.argv[6])
exchange_angle 	= float(sys.argv[7])
delta			= float(sys.argv[8])
max_energy		= float(sys.argv[9])
energy_res 		= float(sys.argv[10])

secondary_variable		= sys.argv[11]
secondary_range_start 	= float(sys.argv[12])
secondary_range_stop 	= float(sys.argv[13])
secondary_range_spacing = float(sys.argv[14])

legend_location	= sys.argv[15]
text_location	= sys.argv[16]
suppress_zeros	= sys.argv[17]
dashed_line		= sys.argv[18]
y_min			= float(sys.argv[19])
y_max			= float(sys.argv[20])

res_dir			= sys.argv[21]
fig_dir			= sys.argv[22]

energy_range	= frange(-max_energy, max_energy, energy_res)
secondary_range	= frange(secondary_range_start, secondary_range_stop, secondary_range_spacing)

pos = 22

# Generate data array from the results.

dos_cut = []

for j in secondary_range:
	exec("%s = %.2f" % (secondary_variable, j))

	file = res_dir + 'Phase_difference_' + "{0:.2f}".format(phase) + '/Rashba_' + "{0:.2f}".format(alpha_1) + '_' + "{0:.2f}".format(alpha_2) + '_' + "{0:.2f}".format(alpha_3) + '_' + "{0:.2f}".format(alpha_4) + '/h-' + "{0:.2f}".format(exchange) + '-angle-' + "{0:.2f}".format(exchange_angle) + '-DOS.txt'

	get_data(file, pos, max_energy, energy_res)

# Plot the results with proper labeling.

colors = ['b', 'g', 'r', 'c', 'm', '#FF9900', '#66FF33', 'k', '#FFFF00', '#996633', '#FF99FF']

fig = plt.figure()

for i in range(len(secondary_range)):
	if secondary_variable == "phase":
		plt.plot(frange(-max_energy, max_energy, energy_res), dos[i], '-', linewidth=1.5, color = colors[i], label = '$\phi = ' + '{0:.2f}'.format(secondary_range[i]) + ' \pi$')
	elif secondary_variable == "alpha_1":
		plt.plot(frange(-max_energy, max_energy, energy_res), dos[i], '-', linewidth=1.5, color = colors[i], label = r'$\alpha_1 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "alpha_2":
		plt.plot(frange(-max_energy, max_energy, energy_res), dos[i], '-', linewidth=1.5, color = colors[i], label = r'$\alpha_2 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "alpha_3":
		plt.plot(frange(-max_energy, max_energy, energy_res), dos[i], '-', linewidth=1.5, color = colors[i], label = r'$\alpha_3 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "alpha_4":
		plt.plot(frange(-max_energy, max_energy, energy_res), dos[i], '-', linewidth=1.5, color = colors[i], label = r'$\alpha_4 = ' + '{0:.2f}'.format(secondary_range[i]) + '/ L$')
	elif secondary_variable == "exchange":
		plt.plot(frange(-max_energy, max_energy, energy_res), dos_cut[i], '-', linewidth=1.5, color = colors[i], label = '$h = ' + '{0:.2f}'.format(secondary_range[i]) + ' E_T$')
	else:
		plt.plot(energy_range, dos_cut[i], '-', linewidth=1.5, color = colors[i], label = r'$\theta = ' + '{0:.2f}'.format(secondary_range[i]) + '^\circ$')

plt.xlabel(r'$E/E_T$', fontsize = 19)
plt.ylabel(r'$N/N_0$', fontsize = 19)

# Prints horizontal dashed line to at y = 0

if dashed_line == "True":
	plt.axhline(y=1, xmin=energy_range[0], xmax=energy_range[-1], color='k', linestyle='--')
elif dashed_line != "False":
	sys.exit("Please set dashed_line parameter to True or False.")

plt.xticks(fontsize = 14)
plt.yticks(fontsize = 14)

plt.xlim(-max_energy, max_energy)
if abs(y_max - y_min) > 0.0001:
	plt.ylim(y_min, y_max)

# This adds a text stating the rest of the variables used in the calculation

ax = fig.add_subplot(111)
add_text(text_location)

# Legend box containing labels.

if legend_location != "none":
	plt.legend(loc=legend_location, handlelength=1.6, handletextpad=0.4, fontsize = 16)

# Finally save the figure into the subdirectory specified in the plot script.

fig.savefig(fig_dir + 'Figure.png')
