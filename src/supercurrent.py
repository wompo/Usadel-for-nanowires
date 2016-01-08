from __future__ import print_function
import scikits.bvp_solver
import numpy as np
from math import pi, sqrt, sin, cos
import fimport
fimport.install(reload_support=True)
import equations
import sys

# # # # # # # # # # # # #
# Function definitions. #
# # # # # # # # # # # # #

# Function to load previous solutions to be used as initial guesses.
def load(sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i):
	return scikits.bvp_solver.Solution.load(sol_dir + 'Phase_difference_%.2f/T_%.2f/Rashba_%.2f_%.2f_%.2f_%.2f/h_%.2f_angle_%.2f/Solution-i-%i.sol' % (phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i))

# Function used to solve the problem.
def solve(runge_kutta, tolerance, initial_guess):
	return scikits.bvp_solver.solve(bvp_problem = problem_definition,
						                    method = runge_kutta,
											tolerance = tolerance,
						                    solution_guess = initial_guess)

# # # # # # # # # # # # # # # 
# Function definitions end. #
# # # # # # # # # # # # # # #

L = 1 # Length of the nanowire.

# Get all the parameters from the run script.

phase			= float(sys.argv[1])

T				= float(sys.argv[2])
alpha_1			= float(sys.argv[3])
alpha_2			= float(sys.argv[4])
alpha_3			= float(sys.argv[5])
alpha_4			= float(sys.argv[6])
h				= float(sys.argv[7])
h_angle			= float(sys.argv[8])
delta			= float(sys.argv[9])
phase_spacing	= float(sys.argv[10])
T_spacing		= float(sys.argv[11])
alpha1_spacing	= float(sys.argv[12])
alpha2_spacing	= float(sys.argv[13])
alpha3_spacing	= float(sys.argv[14])
alpha4_spacing	= float(sys.argv[15])
h_spacing		= float(sys.argv[16])
h_angle_spacing	= float(sys.argv[17])
sol_dir			= sys.argv[18]
res_dir			= sys.argv[19]
matsubara_sum	= int(sys.argv[20])
tolerance		= float(sys.argv[21])

s0 = np.matrix([[1,0],[0,1]])		# Pauli matrices
sx = np.matrix([[0,1],[1,0]])
sy = np.matrix([[0, -1j],[1j, 0]])
sz = np.matrix([[1,0],[0,-1]])

A_1   = alpha_1*sx + alpha_2*sz		# Spin-orbit terms
A_3   = alpha_3*sx + alpha_4*sz

if matsubara_sum == 0:
	i_range=range(int(round(100/(pi*T)-0.5))+1) # Default accuracy omega_max = 200 for all temperatures.
else:
	i_range=range(matsubara_sum)

S, densitysum, supercurrent = ([] for m in range(3))

for n in range(7):
	densitysum.append(np.zeros(45, dtype = complex))
	supercurrent.append(np.zeros(45))

for i in i_range:
	omega = 2 * pi * T * (i + 0.5)
	currentdens = []
	for m in range(7):
		currentdens.append(np.zeros(45, dtype = complex))

	# This function calls the FORTRAN subroutine which contains the Usadel equations and handles
	# all the matrix products. If you need to change something in the Usadel equations, change it
	# in the subroutine.

	def function(X, Y):
		Y = Y.reshape((8, 2, 2))
		res = equations.mysub(Y, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, omega)
		return res.ravel()
    
	def boundary_conditions(Ya, Yb):
		real_left  = delta/(omega+sqrt(delta**2+omega**2))*sin(-0.5*phase*pi)
		imag_left  = delta/(omega+sqrt(delta**2+omega**2))*cos(-0.5*phase*pi)
		real_right = delta/(omega+sqrt(delta**2+omega**2))*sin(0.5*phase*pi)
		imag_right = delta/(omega+sqrt(delta**2+omega**2))*cos(0.5*phase*pi)
         
        # The difference of the current value to the required boundary condition on the left.
		# A clean contact with a bulk superconductor is assumed.
  
		BCa = np.array([Ya[0]               ,
					    Ya[1]  - real_left  ,
					    Ya[2]  + real_left  ,
					    Ya[3]               ,
					    Ya[4]               ,
					    Ya[5]  + imag_left  ,
					    Ya[6]  - imag_left  ,
					    Ya[7]               ,
					    Ya[8]               ,
					    Ya[9]  + real_left,
					    Ya[10] - real_left  ,
					    Ya[11]              ,
					    Ya[12]              ,
					    Ya[13] + imag_left  ,
					    Ya[14] - imag_left  ,
					    Ya[15]              ])

        # The difference of the current value to the required boundary condition on the right
		# A clean contact with a bulk superconductor is assumed.

		BCb = np.array([Yb[0]               ,
					    Yb[1]  - real_right ,
					    Yb[2]  + real_right ,
					    Yb[3]               ,
					    Yb[4]               ,
					    Yb[5]  + imag_right ,
					    Yb[6]  - imag_right ,
					    Yb[7]               ,
					    Yb[8]               ,
					    Yb[9]  + real_right ,
					    Yb[10] - real_right ,
					    Yb[11]              ,
					    Yb[12]              ,
					    Yb[13] + imag_right ,
					    Yb[14] - imag_right ,
					    Yb[15]              ])

		return BCa, BCb

	# The definition of the boundary value problem.

	problem_definition = scikits.bvp_solver.ProblemDefinition(num_ODE = 32,
			                                                  num_parameters = 0,
			                                                  num_left_boundary_conditions = 16,
			                                                  boundary_points = (0, L),
			                                                  function = function,
			                                                  boundary_conditions = boundary_conditions)

	# This is the initial guess that is used if solutions for previous exhance fields/Rashba fields
	# are not available. The current guess is just the half of the boundary values. Terms that are zero
	# are replaced with small numbers to avoid numerical problems. You can change the initial guess to
	# something more appropriate for your problem.

	def guess(X):
		constreal = delta/(omega+sqrt(delta**2+omega**2))*sin(0.5*phase)
		constimag = delta/(omega+sqrt(delta**2+omega**2))*cos(0.5*phase)
		return np.array([ 	0.01 , -constreal/2   ,  constreal/2 ,  0.01 ,
			           		0.01 ,  constimag/2   , -constimag/2 ,  0.01 ,
			           		0.01 ,  constreal/2   , -constreal/2 ,  0.01 ,
			           		0.01 ,  constimag/2   , -constimag/2 ,  0.01 ,
							0.01 ,     0.01       ,     0.01     , -0.01 ,
						   -0.01 ,     0.01       ,     0.01     ,  0.01 ,
						   -0.01 ,     0.01       ,     0.01     ,  0.01 ,
							0.01 ,     0.01       ,     0.01     , -0.01 ])

	# Try to use a range of previous solutions and the above guess as initial guesses.
	# First use Runge-Kutta method of the order 4 and then the order 6 for all initial guesses.

	guess_list = [	[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h - h_spacing, h_angle, i],			# Try with previous exchange field
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i - 1],					# Try with previous Matsubara frequence
					[sol_dir, phase, T, alpha_1 - alpha1_spacing, alpha_2, alpha_3, alpha_4, h, h_angle, i],	# Try with previous alpha_1
					[sol_dir, phase, T, alpha_1, alpha_2 - alpha2_spacing, alpha_3, alpha_4, h, h_angle, i],	# Try with previous alpha_2
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3 - alpha3_spacing, alpha_4, h, h_angle, i],	# Try with previous alpha_3
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4 - alpha4_spacing, h, h_angle, i],	# Try with previous alpha_4
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle - h_angle_spacing, i],	# Try with previous h_angle
					[sol_dir, phase - phase_spacing, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i],		# Try with previous phase
					[sol_dir, phase, T - T_spacing, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i],			# Try with previous T
					guess,																						# Try with the above guess function
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h - 2*h_spacing, h_angle, i],		# Try with even earlier paremeters
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i - 2],					# if everything else fails.
					[sol_dir, phase, T, alpha_1 - 2*alpha1_spacing, alpha_2, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2 - 2*alpha2_spacing, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3 - 2*alpha3_spacing, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4 - 2*alpha4_spacing, h, h_angle, i],
					[sol_dir, phase - 2*phase_spacing, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T - 2*T_spacing, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h - 3*h_spacing, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i - 3],
					[sol_dir, phase, T, alpha_1 - 3*alpha1_spacing, alpha_2, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2 - 3*alpha2_spacing, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3 - 3*alpha3_spacing, alpha_4, h, h_angle, i],
					[sol_dir, phase, T, alpha_1, alpha_2, alpha_3, alpha_4 - 3*alpha4_spacing, h, h_angle, i],
					[sol_dir, phase - 3*phase_spacing, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i],
					[sol_dir, phase, T - 3*T_spacing, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i]		]

	solution = None

	# For zero-energy, try with the above guess function first as using the previous
	# results might give unstable minima.

	for initial_guess in guess_list:
		if solution != None:
			break
		for runge_kutta in [4, 6]:
			try:
				if isinstance(initial_guess, list):
					solution = solve(runge_kutta, tolerance, load(*initial_guess))
				else:
					solution = solve(runge_kutta, tolerance, initial_guess)
				break
			except:
				pass

	# Save the obtained solution into the subdirectory specified in the run script.

	solution.save(sol_dir + 'Phase_difference_%.2f/T_%.2f/Rashba_%.2f_%.2f_%.2f_%.2f/h_%.2f_angle_%.2f/Solution-i-%i.sol' % (phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, i))

	# Set a grid for the solution.

	B = np.linspace(0, L, 45)
	S.append(solution(B))

	# Take the Riccati parameters and their derivatives from the obtained solution to be used in the
	# calculation of the supercurrent and the spin-supercurrents.

	g, gt, dg, dgt, n, nt = ([] for m in range(6))
	a = [[] for m in range(7)]

	for k in range(45):
		g.append(np.matrix([[S[i][0][k]    + 1j*S[i][4][k] , S[i][1][k]  + 1j*S[i][5][k]] , [S[i][2][k]   + 1j*S[i][6][k] , S[i][3][k]  + 1j*S[i][7][k]  ]]))
		gt.append(np.matrix([[S[i][8][k]   + 1j*S[i][12][k], S[i][9][k]  + 1j*S[i][13][k]], [S[i][10][k]  + 1j*S[i][14][k], S[i][11][k] + 1j*S[i][15][k] ]]))
		dg.append(np.matrix([[S[i][16][k]  + 1j*S[i][20][k], S[i][17][k] + 1j*S[i][21][k]], [S[i][18][k]  + 1j*S[i][22][k], S[i][19][k] + 1j*S[i][23][k] ]]))
		dgt.append(np.matrix([[S[i][24][k] + 1j*S[i][28][k], S[i][25][k] + 1j*S[i][29][k]], [S[i][26][k]  + 1j*S[i][30][k], S[i][27][k] + 1j*S[i][31][k] ]]))

		n.append((s0 + g[k] * gt[k]).getI())
		nt.append((s0 + gt[k] * g[k]).getI())

        # Supercurrent

		a[0].append(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k]))

        # Spin-current sigma_1

		a[1].append((s0 + sx)*(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k])))

		a[2].append((s0 - sx)*(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k])))

        # Spin-current sigma_2

		a[3].append((s0 + sy)*(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k])))

		a[4].append((s0 - sy)*(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k])))
        
        # Spin-current sigma_3

		a[5].append((s0 + sz)*(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k])))

		a[6].append((s0 - sz)*(n[k] * (g[k] * dgt[k] - dg[k] * gt[k]) * n[k] - nt[k] * (gt[k] * dg[k] - dgt[k] * g[k]) * nt[k] + 1j * (n[k] * (A_3 * g[k] + g[k] * A_3) * gt[k] * n[k] + n[k] * g[k] * (A_3 * gt[k] + gt[k] * A_3) * n[k] + nt[k] * (A_3 *gt[k] + gt[k] * A_3) * g[k] * nt[k] + nt[k] * gt[k] * (A_3 * g[k] + g[k] * A_3) * nt[k])))

		for m in range(7):
			currentdens[m][k] = a[m][k].trace().item(0).imag*1j
			densitysum[m][k] += currentdens[m][k]

# Write the results into the subdirectory specified in the run script.

for m in range(7):
	for k in range(45):
		supercurrent[m][k] = (2*pi*T*densitysum[m][k]*1j).real/2
	with open(res_dir + 'Phase_difference_%.2f/T_%.2f/Rashba_%.2f_%.2f_%.2f_%.2f/h-%.2f-angle-%.2f-supercurrent.txt' % (phase, T, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle) , 'w') as output_current:
		print(str(h) + '\n' + str(supercurrent[0]) + '\n', file = output_current)
		print('sx_up' + '\n' + str(supercurrent[1]) + '\n', file = output_current)
		print('sx_down' + '\n' + str(supercurrent[2]) + '\n', file = output_current)
		print('sy_up' + '\n' + str(supercurrent[3]) + '\n', file = output_current)
		print('sy_down' + '\n' + str(supercurrent[4]) + '\n', file = output_current)
		print('sz_up' + '\n' + str(supercurrent[5]) + '\n', file = output_current)
		print('sz_down' + '\n' + str(supercurrent[6]) + '\n', file = output_current)
