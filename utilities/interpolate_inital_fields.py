"""
Usage: Interpolate initial fields according to grid resolution

ndinterpolator module is part of the Bladex package
https://github.com/mathLab/BladeX

M. Gadalla, M. Tezzele, A. Mola, G. Rozza, (2019).
 BladeX: Python Blade Morphing.
 Journal of Open Source Software, 4(34), 1203, 
 https://doi.org/10.21105/joss.01203
"""
import numpy as np
from ndinterpolator import reconstruct_f

def openfoam_write_internalField_scalar(filename, rbf_points, array):
	outfile = ''
	outfile += str(rbf_points) 
	outfile += '\n(\n'
	for val in array:
		outfile += "{:.7f}".format(val) + '\n'
	outfile += ')'
	f = open(filename, "w")
	f.write(outfile)
	f.close()

rbf_points = 1350
out_states = 'out_states'

grid = np.loadtxt(out_states + '/grid')

for f in ['T', 'CH4', 'O2', 'N2', 'CO2', 'H2O']:
    field = np.loadtxt(out_states + '/' + f)
    xx = np.linspace(grid[0], grid[-1], num=rbf_points)
    yy = np.zeros(rbf_points)
    reconstruct_f(original_input=grid, original_output=field,
				  rbf_input=xx, rbf_output=yy,
				  basis='beckert_wendland_c2_basis', radius=2.0)
    openfoam_write_internalField_scalar(filename=out_states + '/' + f+'.dat',
    									rbf_points=rbf_points, array=yy)
