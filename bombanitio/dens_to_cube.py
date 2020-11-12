import numpy as np
from bombanitio.toolbox import  cube
from bombanitio.toolbox import  CubeFileTools
from bombanitio.toolbox.constants import  bohr_in_a 
from bombanitio.toolbox import constants
from bombanitio import read_xyz
from bombanitio import create_basis
from bombanitio.toolbox import lib_dme


def basis_func_to_grid(grid, orb_type ,exp, coeff, pos):
   #print(grid.shape)
   #dist1 =  np.linalg.norm(grid - pos, axis=0)
   dist2 = grid.T - pos
   dist2 = dist2.T
   dist1 =  np.linalg.norm(dist2, axis=0)
   #print(grid[0,0:3,0:3,0])
   #print(grid[0,104:108,104:108,0])
   #print()
   #print(grid[0,50:53,50:53,0])
   #print(grid[:,53,53,53])
   #print(dist1.shape)
   #data = np.exp(np.
   #exp1 = np.exp(dist1)
   #print(exp1.shape)
   k = orb_type - 1
   data = np.zeros((grid.shape[1],grid.shape[2],grid.shape[3]))
   for i in range(exp.shape[0]):
           pre_data = np.exp(- exp[i] * dist1**2)
           #print(pre_data.shape, abs(pre_data).sum())
           if orb_type > 0:
               pre_data *= grid[k,:,:,:] - pos[k]
           pre_data *= coeff[i]
           if orb_type == 0:
               pre_data *= (2.0 * exp[i] / np.pi)**0.75
           else:
               pre_data *= (128 *exp[i]**5 / np.pi**3)**0.25
           data += pre_data
           #print(data.shape, abs(data).sum())
   #        print(data[52,52,52])
   #print(data.shape, abs(data).sum())
   #print(data[52,52,52])
   return data

#four *= (128 *exp1**5 / np.pi**3)**0.25
#    four *= (2.0 * exp2 / np.pi)**0.75
#def dens_to_cube(fn1, dens, exp_mat, coeff_mat, orb_typ, orb_pos):
#def dens_to_cube(fn1):
def dens_to_cube():
   #fn1 = "../data/dens_h2o.cube" 
   #fn1 = "/home/dressler/projects/HF/develop_bombanitio/minus9/data/dens_h2o.cube" 
   #fn1 = "/home/dressler/projects/HF/develop_bombanitio/minus9/data/DENSITY.cube" 
   #cube1 = CubeFileTools.LoadCellData(fn1)
   #for i in cube1.keys():
   #    print(i)
   #print(cube1["origin_au"])
   ##print(cube1["cell_au"])
   #print(cube1["data"].shape)
   #print(cube1["coords_au"])
   np.set_printoptions(precision=4, suppress=True) 
   print("Reading coordinates ... ")
   fn_dens = "density_mat0"
   #dens = np.loadtxt("density_mat0")
   dens = np.loadtxt(fn_dens)
   coord, atom = read_xyz.easy_read("coord.xyz", None, False, False)
   coord = coord[0,:,:]
   coord = read_xyz.xyz_to_bohr(coord)
   zoa = read_xyz.atom_to_zoa(atom)
   noa = atom.shape[0]
   noo, orb_coord, orb_type, coef_mat, exp_mat = create_basis.basisset(noa, zoa, coord)
   #read_xyz.plot_initial_orbs(noo, orb_coord, orb_type, coef_mat, exp_mat) 
   
   #n_grid =  108 
   #n_grid =  216
   n_grid =  140
   cube_length_au =  7 / bohr_in_a
   #cube_length_au =  10 / bohr_in_a
   #data = np.zeros((108,108,108))
   data = np.zeros((n_grid,n_grid,n_grid))
   mesh = data.shape
   #origin_au = [-5,-5,-5]/bohr_in_a #punkt in linker unterer ecke, nicht die box mitte!!!
   origin_au = np.array([cube_length_au/2 ,cube_length_au/2,cube_length_au/2]) #punkt in rechter!!!!! unterer ecke, nicht die box mitte!!!
   #origin_au = np.array([cube_length_au ,cube_length_au,cube_length_au]) #punkt in rechter!!!!! unterer ecke, nicht die box mitte!!!
   cell_au = np.eye(3) * cube_length_au / n_grid
#def CalcGridPositions(cell_au, mesh, origin_au=None):
   grid = CubeFileTools.CalcGridPositions(cell_au, mesh, origin_au)

   comment1 =  "bla"
   comment2 =  "blub"
   #numbers = zoa
   
   cell_data = dict()
   cell_data['numbers'  ] = zoa
   #cell_data['coords_au'] = coords_au
   cell_data['coords_au'] = coord
   cell_data['cell_au'  ] = cell_au
   #cell_data['symbols'  ] = [constants.symbols[int(e)-1] for e in cell_data['numbers']]
   cell_data['symbols'  ] = atom
   cell_data['species'  ] = list(set(cell_data['symbols']))
   cell_data['mesh'     ] = data.shape
   cell_data['d3r_au'   ] = cell_au[0,0]*cell_au[1,1]*cell_au[2,2]
   cell_data['volume_au'] = cell_data['d3r_au']*data.shape[0]*data.shape[1]*data.shape[2]
   cell_data['r_au'     ] = grid
   #cell_data['r_au'     ] = CalcGridPositions(cell_au, data.shape, origin_au=origin_au)
   cell_data['data'     ] = data
   cell_data['comment1' ] = comment1
   cell_data['comment2' ] = comment2
   cell_data['origin_au'] = -origin_au
#species  symbols  numbers  comment1  comment2   
   #def WriteCubeFile(filename, comment1, comment2, numbers, coords, cell, data, origin=np.zeros(3)):
   fn_out = "test.cube"
   #cube.WriteCubeFile(fn_out, cell_data['comment1' ], cell_data['comment2' ], cell_data['numbers'  ], cell_data['coords_au'], cell_data['cell_au'  ], cell_data['data'     ], cell_data['origin_au'])
  



   mompol   = dict()
   moments  = dict()


   exp_order = 1
   n_states = 3
   for i_order in range(1, exp_order+1):
        #f cartesian:
        #   mompol['%02d'%i_order] = CubeFileTools.CartesianMoments(r_au, order=i_order)
        #lse:
            #mompol['%02d'%i_order] = CubeFileTools.CartesianMoments(r_au, order=i_order)
            mompol['%02d'%i_order] = CubeFileTools.CartesianMoments(grid, order=i_order)
            #mompol['%02d'%i_order] = SphericalHarmonics.RealRegularSolidHarmonics(r_au, order=i_order)
#moments['%02d'%i_order] = np.zeros((n_states, mompol['%02d'%i_order].shape[0]))
#print(mompol)
#print(moments)
#mom_pol_ar = np.zeros((n_states, 108, 108 , 108))
#mom_pol_ar = np.zeros((n_states + 1 , 108, 108 , 108))
#mom_pol_ar = np.zeros((n_states  , 108, 108 , 108))
   mom_pol_ar = np.zeros((n_states  , n_grid, n_grid, n_grid))
#for i in range(n_states):
#mom_pol_ar[0] = np.ones((108, 108 , 108))

   print(n_states)

   for j in range(n_states):
        #k = j +1
        k = j
        if (j < 3) and (j > -1):
            mom_pol_ar[k] = mompol['01'][j]
        elif (j < 9) and (j > 2):
#            ipdb.set_trace()
            mom_pol_ar[k] = mompol['02'][j-3]
        elif (j < 19) and (j > 8):
            mom_pol_ar[k] = mompol['03'][j-9]
        elif (j < 34) and (j > 18):
            mom_pol_ar[k] = mompol['04'][j-19]
        elif (j < 55) and (j > 33):
            mom_pol_ar[k] = mompol['05'][j-34]
        else:
            print('index error')
#mom_pol_ar.shape



#if 5 > 10:
        #lime.print_states_dict(mom_pol_ar*0.001, cell_data5, n_states, '/shared/dressler/chi/cartesian-storage-fac-0.001/cartesian-functions-%05d', pure = True, pert = True)
#lime.print_states_dict(mom_pol_ar*0.001, cell_data5, n_states + 1, '/shared/dressler/chi/paper_no4_theorem/data/v1/center-chi4-basis-polynoms-shift/cartesian-functions-%05d', pure = True, pert = True)
#lime. print_states_dict(mom_pol_ar*0.001, cell_data, n_states + 1, 'cartesian-functions-%05d', pure = True, pert = True)
   #lime.print_states_dict(mom_pol_ar*0.001, cell_data, n_states, 'cartesian-functions-%05d', pure = True, pert = True)
   for i in range(n_states):
       fnx = "cartesian-functions"+str(i)+".cube"
       cube.WriteCubeFile(fnx, cell_data['comment1' ], cell_data['comment2' ], cell_data['numbers'  ], cell_data['coords_au'], cell_data['cell_au'  ], mom_pol_ar[i], cell_data['origin_au'])








   
   print("start to calc density on grid")
   #print(grid.shape)
   #dist1 =  np.linalg.norm(grid, axis=0)
   #print(dist1.shape)
   ##data = np.exp(np.
   #exp1 = np.exp(dist1)
   #print(exp1.shape)
   basis_at_grid = []
   for i in range(noo):
       #noo, orb_coord, orb_type, coef_mat, exp_mat
       #data = basis_func_to_grid(grid, orb_type ,exp, coeff, pos)
       #if orb_type[i] > 0:
       print("input ",orb_type[i], exp_mat[i], coef_mat[i], orb_coord[i])
       data = basis_func_to_grid(grid, orb_type[i] ,exp_mat[i], coef_mat[i], orb_coord[i])
       #print(data.shape, abs(data).sum())
       print(data.shape, abs(data**2).sum())
       print(data[52,52,52])
       str1 = "basis_function"+str(i)+".cube"
       cube.WriteCubeFile(str1, cell_data['comment1' ], cell_data['comment2' ], cell_data['numbers'  ], cell_data['coords_au'], cell_data['cell_au'  ], data, cell_data['origin_au'])
       basis_at_grid.append(data)
       #fn_dens + 
   #basis_func_to_grid(grid,orb_type ,exp, coeff, pos)
   basis_at_grid = np.array(basis_at_grid)
   basis_overlapp_states = []
   for i in range(noo):
       basis_overlapp_states.append([])
       for j in range(noo):
           basis_overlapp_states[i].append(basis_at_grid[i]*basis_at_grid[j])
   basis_overlapp_states = np.array(basis_overlapp_states)
   print(basis_overlapp_states.shape)
   strange_dens_mat = basis_overlapp_states.sum(axis=(2,3,4))
   print(strange_dens_mat)
   #basis_overlapp_states = basis_overlapp_states*dens
   diff_list = []
   for fn_dens in  [ "density_mat0", "density_mat1", "density_mat2", "density_mat3"]:
       #fn_dens = "density_mat0"
       dens = np.loadtxt(fn_dens)
       basis_overlapp_states1 = basis_overlapp_states*dens[:,:,np.newaxis, np.newaxis, np.newaxis]
       data = basis_overlapp_states1.sum(axis=(0,1))
       str1 = fn_dens+".cube"
       cube.WriteCubeFile(str1, cell_data['comment1' ], cell_data['comment2' ], cell_data['numbers'  ], cell_data['coords_au'], cell_data['cell_au'  ], data, cell_data['origin_au'])
       print("summ of electrons: ", data.sum())
       if fn_dens == "density_mat0":
           ref = np.copy(data)
       else:
        diff = data - ref
        str1 = "diff_"+fn_dens+".cube"
        cube.WriteCubeFile(str1, cell_data['comment1' ], cell_data['comment2' ], cell_data['numbers'  ], cell_data['coords_au'], cell_data['cell_au'  ], diff, cell_data['origin_au'])
        print("summ of diff electrons: ", diff.sum())
        diff_list.append(diff)
   diff_list = np.array(diff_list)
   #mom_mat = lib_dme.create_overlap_mat(diff_list[1:], mom_pol_ar)
   mom_mat = lib_dme.create_overlap_mat(diff_list, mom_pol_ar)
   print("mom mat")
   print(mom_mat)
   print(abs(diff_list[0]).sum())
   print((diff_list[0]*mom_pol_ar[1]).sum())
   special = diff_list[0]*mom_pol_ar[1]
   cube.WriteCubeFile("special.cube", cell_data['comment1' ], cell_data['comment2' ], cell_data['numbers'  ], cell_data['coords_au'], cell_data['cell_au'  ], special, cell_data['origin_au'])
