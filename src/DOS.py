from __future__ import print_function
import scikits.bvp_solver
import numpy as np
from cmath import pi, sqrt, exp
import fimport
fimport.install(reload_support=True)
import equations_real_time, g_bcs
import sys

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

# Function to load previous solutions to be used as initial guesses.
def load(sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps):
	return scikits.bvp_solver.Solution.load(sol_dir + 'Phase_difference_%.2f/Rashba_%.2f_%.2f_%.2f_%.2f/h_%.2f_angle_%.2f/Solution-i-%.2f.sol' % (phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps))

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
alpha_1			= float(sys.argv[2])
alpha_2			= float(sys.argv[3])
alpha_3			= float(sys.argv[4])
alpha_4			= float(sys.argv[5])
h				= float(sys.argv[6])
h_angle			= float(sys.argv[7])
delta			= float(sys.argv[8])
phase_spacing	= float(sys.argv[9])
alpha1_spacing	= float(sys.argv[10])
alpha2_spacing	= float(sys.argv[11])
alpha3_spacing	= float(sys.argv[12])
alpha4_spacing	= float(sys.argv[13])
h_spacing		= float(sys.argv[14])
h_angle_spacing	= float(sys.argv[15])
sol_dir			= sys.argv[16]
res_dir			= sys.argv[17]
max_energy		= float(sys.argv[18])
energy_res		= float(sys.argv[19])
tolerance		= float(sys.argv[20])

s0 = np.matrix([[1,0],[0,1]])		# Pauli matrices
sx = np.matrix([[0,1],[1,0]])
sy = np.matrix([[0, -1j],[1j, 0]])
sz = np.matrix([[1,0],[0,-1]])

A_1   = alpha_1*sx + alpha_2*sz		# Spin-orbit terms
A_3   = alpha_3*sx + alpha_4*sz

i_range = frange(0.00, max_energy, energy_res)

dos = [[] for m in range(len(i_range))]

l = 0
for eps in i_range:
	offset = 0.001	# Offset specifying the location of the poles of the Retarded Green's functions.

	# Get the boundary conditions from FORTRAN subroutines.

	g_bcs_real_left , g_bcs_imag_left , gt_bcs_real_left , gt_bcs_imag_left  = g_bcs.gsub(delta, eps, offset, -0.5*phase*pi)
	g_bcs_real_right, g_bcs_imag_right, gt_bcs_real_right, gt_bcs_imag_right = g_bcs.gsub(delta, eps, offset, 0.5*phase*pi)


	# This function calls the FORTRAN subroutine which contains the Usadel equations and handles
	# all the matrix products. If you need to change something in the Usadel equations, change it
	# in the subroutine.

	def function(X, Y):
		Y = Y.reshape((8, 2, 2))
		res = equations_real_time.mysub(Y, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps)
		return res.ravel()

	def boundary_conditions(Ya, Yb):

        # The difference of the current value to the required boundary condition on the left.
		# A clean contact with a bulk superconductor is assumed.

		BCa = np.array([Ya[0]  - g_bcs_real_left[0,0]   ,
				        Ya[1]  - g_bcs_real_left[0,1]   ,
				        Ya[2]  - g_bcs_real_left[1,0]   ,
				        Ya[3]  - g_bcs_real_left[1,1]   ,
				        Ya[4]  - g_bcs_imag_left[0,0]   ,
				        Ya[5]  - g_bcs_imag_left[0,1]   ,
				        Ya[6]  - g_bcs_imag_left[1,0]   ,
				        Ya[7]  - g_bcs_imag_left[1,1]   ,
				        Ya[8]  - gt_bcs_real_left[0,0]  ,
				        Ya[9]  - gt_bcs_real_left[0,1]  ,
				        Ya[10] - gt_bcs_real_left[1,0]  ,
				        Ya[11] - gt_bcs_real_left[1,1]  ,
				        Ya[12] - gt_bcs_imag_left[0,0]  ,
				        Ya[13] - gt_bcs_imag_left[0,1]  ,
				        Ya[14] - gt_bcs_imag_left[1,0]  ,
				        Ya[15] - gt_bcs_imag_left[1,1]  ])
		  
        # The difference of the current value to the required boundary condition on the right
		# A clean contact with a bulk superconductor is assumed.
		  
		BCb = np.array([Yb[0]  - g_bcs_real_right[0,0]  ,
				        Yb[1]  - g_bcs_real_right[0,1]  ,
				        Yb[2]  - g_bcs_real_right[1,0]  ,
				        Yb[3]  - g_bcs_real_right[1,1]  ,
				        Yb[4]  - g_bcs_imag_right[0,0]  ,
				        Yb[5]  - g_bcs_imag_right[0,1]  ,
				        Yb[6]  - g_bcs_imag_right[1,0]  ,
				        Yb[7]  - g_bcs_imag_right[1,1]  ,
				        Yb[8]  - gt_bcs_real_right[0,0] ,
				        Yb[9]  - gt_bcs_real_right[0,1] ,
				        Yb[10] - gt_bcs_real_right[1,0] ,
				        Yb[11] - gt_bcs_real_right[1,1] ,
				        Yb[12] - gt_bcs_imag_right[0,0] ,
				        Yb[13] - gt_bcs_imag_right[0,1] ,
				        Yb[14] - gt_bcs_imag_right[1,0] ,
				        Yb[15] - gt_bcs_imag_right[1,1] ])
				             
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
		return np.array([ g_bcs_real_left[0,0]  , g_bcs_real_left[0,1]  , g_bcs_real_left[1,0]  , g_bcs_real_left[1,1]  ,
					      g_bcs_imag_left[0,0]  , g_bcs_imag_left[0,1]  , g_bcs_imag_left[1,0]  , g_bcs_imag_left[1,1]  ,
					      gt_bcs_real_left[0,0] , gt_bcs_real_left[0,1] , gt_bcs_real_left[1,0] , gt_bcs_real_left[1,1] ,
					      gt_bcs_imag_left[0,0] , gt_bcs_imag_left[0,1] , gt_bcs_imag_left[1,0] , gt_bcs_imag_left[1,1] ,
					             0.01           ,        0.01           ,        0.01           ,       -0.01           ,
					            -0.01           ,        0.01           ,        0.01           ,        0.01           ,
					            -0.01           ,        0.01           ,        0.01           ,        0.01           ,
					             0.01           ,        0.01           ,        0.01           ,       -0.01           ])

	# Try to use a range of previous solutions and the above guess as initial guesses.
	# First use Runge-Kutta method of the order 4 and then the order 6 for all initial guesses.

	guess_list = [	[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps - energy_res],			# Try with previous energy
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h - h_spacing, h_angle, eps],			# Try with previous exchange field
					[sol_dir, phase, alpha_1 - alpha1_spacing, alpha_2, alpha_3, alpha_4, h, h_angle, eps],		# Try with previous alpha_1
					[sol_dir, phase, alpha_1, alpha_2 - alpha2_spacing, alpha_3, alpha_4, h, h_angle, eps],		# Try with previous alpha_2
					[sol_dir, phase, alpha_1, alpha_2, alpha_3 - alpha3_spacing, alpha_4, h, h_angle, eps],		# Try with previous alpha_3
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4 - alpha4_spacing, h, h_angle, eps],		# Try with previous alpha_4
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle - h_angle_spacing, eps],	# Try with previous h_angle
					[sol_dir, phase - phase_spacing, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps],		# Try with previous phase
					guess,																						# Try with the above guess function
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h - 2*h_spacing, h_angle, eps],		# Try with even earlier paremeters
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps - 2*energy_res],		# if everything else fails.
					[sol_dir, phase, alpha_1 - 2*alpha1_spacing, alpha_2, alpha_3, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2 - 2*alpha2_spacing, alpha_3, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2, alpha_3 - 2*alpha3_spacing, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4 - 2*alpha4_spacing, h, h_angle, eps],
					[sol_dir, phase - 2*phase_spacing, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h - 3*h_spacing, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps - 3*energy_res],
					[sol_dir, phase, alpha_1 - 3*alpha1_spacing, alpha_2, alpha_3, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2 - 3*alpha2_spacing, alpha_3, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2, alpha_3 - 3*alpha3_spacing, alpha_4, h, h_angle, eps],
					[sol_dir, phase, alpha_1, alpha_2, alpha_3, alpha_4 - 3*alpha4_spacing, h, h_angle, eps],
					[sol_dir, phase - 3*phase_spacing, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps] 	]

	solution = None

	# For zero-energy, try with the above guess function first as using the previous
	# results might give unstable minima.

	if eps == 0.00:
		for runge_kutta in [4, 6]:
			try:
				solution = solve(runge_kutta, tolerance, guess)
				break
			except:
				pass

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

	solution.save(sol_dir + 'Phase_difference_%.2f/Rashba_%.2f_%.2f_%.2f_%.2f/h_%.2f_angle_%.2f/Solution-i-%.2f.sol' % (phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle, eps))

	# Set a grid for the solution.

	B = np.linspace(0, L, 45)
	S = solution(B)

	# Take the Riccati parameters and their derivatives from the obtained solution to be used in the
	# calculation of the density of states.

	g, gt, dg, dgt, n, nt = ([] for m in range(6))

	for k in range(len(B)):
		g.append(np.matrix([[S[0][k]    + 1j*S[4][k] , S[1][k]  + 1j*S[5][k]] , [S[2][k]  + 1j*S[6][k] , S[3][k]  + 1j*S[7][k]  ]]))
		gt.append(np.matrix([[S[8][k]   + 1j*S[12][k], S[9][k]  + 1j*S[13][k]], [S[10][k] + 1j*S[14][k], S[11][k] + 1j*S[15][k] ]]))
		dg.append(np.matrix([[S[16][k]  + 1j*S[20][k], S[17][k] + 1j*S[21][k]], [S[18][k] + 1j*S[22][k], S[19][k] + 1j*S[23][k] ]]))
		dgt.append(np.matrix([[S[24][k] + 1j*S[28][k], S[25][k] + 1j*S[29][k]], [S[26][k] + 1j*S[30][k], S[27][k] + 1j*S[31][k] ]]))
		   
		n.append((s0 + g[k] * gt[k]).getI())
		nt.append((s0 + gt[k] * g[k]).getI())

		dos[l].append(0.5*(n[k] * (s0 - g[k] * gt[k]) + nt[k] * (s0 - gt[k] * g[k])).real)

	l += 1

# Write the results into the subdirectory specified in the run script.

with open(res_dir + 'Phase_difference_%.2f/Rashba_%.2f_%.2f_%.2f_%.2f/h-%.2f-angle-%.2f-DOS.txt' % (phase, alpha_1, alpha_2, alpha_3, alpha_4, h, h_angle) , 'w') as output_dos:
	for k in range(len(B)):
		pos = k*1.0/(len(B) - 1)
		print('z = %.2f' % pos, file = output_dos)
		for i in range(len(i_range)):
			print(dos[i][k].item(0), dos[i][k].item(1), dos[i][k].item(2), dos[i][k].item(3), file = output_dos)
