#!/usr/bin/python3
import numpy as np

from bombanitio.toolbox import constants
from bombanitio.toolbox import cube




def LoadCellData(fn_cube):
    #if fn_cube[-5:] == '.cube':
    data, numbers, coords_au, cell_au, comment1, comment2, origin_au = cube.ReadCubeFile(fn_cube)
    # elif fn_cube[-3:] == '.gz':
    #     data, numbers, coords_au, cell_au, comment1, comment2, origin_au = cube.ReadGZIPCubeFile(fn_cube)
    # else:
    #     raise Exception('Unknown fromat?!')
    tmp = abs(cell_au[0,1])+abs(cell_au[0,2])+abs(cell_au[1,0])+abs(cell_au[1,2])+abs(cell_au[2,0])+abs(cell_au[2,1])
    if tmp > 1E-10:
        raise Exception('Only orthorhombic or simple cubic boxes supported!')
    cell_data = dict()
    cell_data['numbers'  ] = numbers
    cell_data['coords_au'] = coords_au
    cell_data['cell_au'  ] = cell_au
    cell_data['symbols'  ] = [constants.symbols[int(e)-1] for e in cell_data['numbers']]
    cell_data['species'  ] = list(set(cell_data['symbols']))
    cell_data['mesh'     ] = data.shape
    cell_data['d3r_au'   ] = cell_au[0,0]*cell_au[1,1]*cell_au[2,2]
    cell_data['volume_au'] = cell_data['d3r_au']*data.shape[0]*data.shape[1]*data.shape[2]
    cell_data['r_au'     ] = CalcGridPositions(cell_au, data.shape, origin_au=origin_au)
    cell_data['data'     ] = data
    cell_data['comment1' ] = comment1
    cell_data['comment2' ] = comment2
    cell_data['origin_au'] = origin_au
    return cell_data



def CalcGridPositions(cell_au, mesh, origin_au=None):
    """returns r_au of shape (mesh1,mesh2,mesh3,3) in units of cell, i.e. default AU"""
    n_x, n_y, n_z = mesh
    a_x, a_y, a_z = cell_au[0][0], cell_au[1][1], cell_au[2][2]
    if origin_au.all() == None:
        origin_au = np.array([n_x*a_x, n_y*a_y, n_z*a_z])/2
    x_au = np.arange(0, n_x*a_x, a_x) - origin_au[0]
    y_au = np.arange(0, n_y*a_y, a_y) - origin_au[1]
    z_au = np.arange(0, n_z*a_z, a_z) - origin_au[2]
    r_au = np.zeros((3, n_x, n_y, n_z)) #NB! Changed to row-major-order!
    r_au[0,:,:,:] = x_au[:,np.newaxis,np.newaxis]
    r_au[1,:,:,:] = y_au[np.newaxis,:,np.newaxis]
    r_au[2,:,:,:] = z_au[np.newaxis,np.newaxis,:]
    return r_au


def CartesianMoments(r_au, order=1):
    x = r_au[0,:,:,:]
    y = r_au[1,:,:,:]
    z = r_au[2,:,:,:]
    if order == 1:
        sh = np.zeros((3, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = x
        sh[1] = y
        sh[2] = z
    elif order == 2:
        sh = np.zeros((6, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = x**2
        sh[1] = x*y
        sh[2] = x*z
        sh[3] = y**2
        sh[4] = y*z
        sh[5] = z**2
    elif order == 3:
        sh = np.zeros((10, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[0] = x**3
        sh[1] = x**2*y
        sh[2] = x**2*z
        sh[3] = y**2*x
        sh[4] = y**3
        sh[5] = y**2*z
        sh[6] = z**2*x
        sh[7] = z**2*y
        sh[8] = z**3
        sh[9] = x*y*z
    elif order == 4:
        sh = np.zeros((15, r_au.shape[1], r_au.shape[2], r_au.shape[3])) #NB! 16 or 15 components ?! http://onlinelibrary.wiley.com/doi/10.1002/3527602747.app4/pdf
        sh[ 0] = x**4
        sh[ 1] = x**3*y
        sh[ 2] = x**3*z
        sh[ 3] = y**3*x
        sh[ 4] = y**4
        sh[ 5] = y**3*z
        sh[ 6] = z**3*x
        sh[ 7] = z**3*y
        sh[ 8] = z**4
        sh[ 9] = x**2*y**2
        sh[10] = x**2*y*z
        sh[11] = x**2*z**2
        sh[12] = y**2*z**2
        sh[13] = x*y**2*z
        sh[14] = x*y*z**2
    elif order == 5: # NB! not yet tested ...
        sh = np.zeros((21, r_au.shape[1], r_au.shape[2], r_au.shape[3]))
        sh[ 0] = x**5
        sh[ 1] = x**4*y
        sh[ 2] = x**4*z
        sh[ 3] = x**3*y**2
        sh[ 4] = x**3*z**2
        sh[ 5] = x**3*y*z
        sh[ 6] = y**5
        sh[ 7] = y**4*x
        sh[ 8] = y**4*z
        sh[ 9] = y**3*x**2
        sh[10] = y**3*z**2
        sh[11] = y**3*x*z
        sh[12] = z**5
        sh[13] = z**4*x
        sh[14] = z**4*y
        sh[15] = z**3*x**2
        sh[16] = z**3*y**2
        sh[17] = z**3*x*y
        sh[18] = x**2*y**2*z
        sh[19] = x**2*y*z**2
        sh[20] = x*y**2*z**2
    else:
        raise Exception('Order %d not supported!'%order)
    return sh

