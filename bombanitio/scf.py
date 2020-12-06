import numpy as np

#def formg(noo, four_mat, p):
#    # ! ! ! calculate g matrix from density matrix and two-electron integrals
#    g = np.zeros((noo, noo))
#    for i in range(noo):
#        for j in range(noo):
#            for k in range(noo):
#                for l in range(noo):
#                    g[i, j] = g[i, j] + p[k, l] * (four_mat[i, j, k, l] - 0.5 * four_mat[i, l, k, j])
#    return g

def formg(nob,two_e_mat, p):
    # ! ! ! calculate g matrix from density matrix and two-electron integrals
    g = np.zeros((nob, nob))
    for i in range(nob):
        for j in range(nob):
            for k in range(nob):
                for l in range(nob):
                    g[i, j] = g[i, j] + p[k, l] * (two_e_mat[i, j, k, l] - 0.5 * two_e_mat[i, l, k, j])
    return g







#def scf(four_mat, noo, noo, x, xt, h, over_mat):
#def scf(four_mat, noe, noo, x, xt, h):
#def scf(four_mat, noe, noo, x, h):
def scf(four_mat, noe, noo, x, h, ext = None):
#def scf(four_mat, n, noo, x, xt, h, s):
    xt = x.T
    # ! ! convergence criterion for density matrix
    #crit = 0.0001
    #deltae 10**-9
    #deltap 10**-5
    crit = 0.00001
    crit_en = 10**-9
    # ! ! maximum number of iterations
    maxit = 2000
    print()
    print('##############################')
    print('start of scf cycle')
    # ! iteration number
    iter = 0
    #oldp = np.zeros((noo, noo))
    en = 0.0
    p = np.zeros((noo, noo))
#    p =np.array([[2.10625, -0.44603, 0.00000, 0.10859, 0.00000, -0.02843, -0.02843],
#     [-0.44603, 1.96730, 0.00000,-0.61771, 0.00000, -0.03406, -0.03406],
#     [0.00000,  0.00000, 0.73554, 0.00000, 0.00000,  0.53974, -0.53974],
#     [0.10859, -0.61771, 0.00000, 1.24023, 0.00000,  0.47272,  0.47272],
#     [0.00000,  0.00000, 0.00000, 0.00000, 2.00000,  0.00000,  0.00000],
#    [-0.02843,-0.03406, 0.53974, 0.47272, 0.00000,  0.60093, -0.19120],
#    [-0.02843,-0.03406,-0.53974, 0.47272, 0.00000, -0.19120,  0.60093]])

    conv1 = False
    for ii in range(maxit):
        if conv1 == False:
            iter = iter + 1
            if (iter >= maxit):
                print('iter_max reached')
            if iter <= maxit:
                print()
                print('##############################')
                print('start of iteration number', iter)
            ##################################################################
            ################IMPORTANT PART OF THE SCF CYCLE STARTS
            #################################################################
            ##! form two-electron part of fock matrix from P
            g = formg(noo, four_mat, p)
            verbose_local = 0
            if verbose_local >= 2:
                print('g', g)
            f = h + g
            if ext is not None:
                f += ext
            # ! calculate electronic energy
            old_en = en
            en = 0.0
            en = en + 0.5 * (p * ( h +f )).sum()
            ######################################################
            #one of the historical simplifications enabled by python in this code compared to fortran:
            #for i in range(noo):
            #    for j in range(noo):
            #        en = en + 0.5 * p[i, j] * (h[i, j] + f[i, j])
            ######################################################
            verbose = 0
            if verbose >= 2:
                print('f', f)
                print()
            print('electronic energy=', en)
            ##! transform Fock matrix using g_tmp for temporal storage
            g_tmp = np.dot(f, x)
            fprime = np.dot(xt, g_tmp)
            e, cprime = np.linalg.eig(fprime)
            idx = e.argsort()[::1]
            # ! transform eigenvectors to get matrix C
            c = np.dot(x, cprime)
            # ! form new density matrix
            ######################################################################
            #### epic python error:
            #### why is this command line wrong:  oldp = p
            ####oldp = p
            ####################################################################### 
            # ! save present density matrix before creating new one
            oldp = np.matrix.copy(p)
            #print("idx", idx)
            for i in range(noo):
                for j in range(noo):
                    ## ! save present density matrix before creating new one
                    ##oldp[i, j] = p[i, j]
                    p[i, j] = 0.0
                    # k laeuft über anzhal der e paare!!!!
                    #real number of electrons devided by 2: rne2
                    #rne2 = 1 #muss verallgemeinert werden, variablen umbenenen
                    #for k in range(rne2):
                    #for k in range(int(noe/2)):
                    #idx[:noe / 2]
                    for k in idx[:int(noe / 2)]:
                        p[i, j] = p[i, j] + 2.0 * c[i, k] * c[j, k]

            ###########################################################################
            ################ important part of scf cycle is finshed##############
            ##########################################################################
            verbose_local = 0
            if verbose_local >= 2:
                print('fprime', fprime)
                print()
                print('cprime', cprime)
                print()
                print('orbital energies', e)
                print()
                print('HF-orbital coefficients', c)
                print()
                print('density matrix p', p)

            delta = 0.0
            for i in range(noo):
                for j in range(noo):
                    delta = delta + (p[i, j] - oldp[i, j]) ** 2
            delta /= noo
            print('delta density', delta)
            
            delta_en = abs(en - old_en)
            print('delta energy', delta_en)


            #if delta <= crit:
            if (delta <= crit) and (delta_en <= crit_en):
                #ent = en + za * zb / r
                print('calculation is converged')
                conv1 = True
                print('electronic energy', en)

            verbose_local = 0
            if verbose_local  == 1:
                print('G', g)
                print('F', f)
                print('E', e)
                print('C', c)
                print('P', p)
            #ps matrix has mulliken populations
                #print('PS', np.dot(oldp, over_mat))

            else:
                if iter >= maxit - 1:
                    print('no convergence in SCF zyklus', maxit)
            
            print('E', e[idx])
            
            #dampong
            #alpha = 0.9
            alpha = 0.5
            p =(1-alpha) *  p + alpha* oldp

    print(f)     
    esort = e[idx]
    csort = c[idx,:]
    return en, p, csort, esort

def calc_tot_en(noa, e_elec, pos_list, charge_list):
    e_tot = e_elec
    for i in range(noa):
        for j in range(noa):
                 if j > i:
        #            print(i,j)
                    e_tot +=  charge_list[i] *  charge_list[j] / np.linalg.norm(pos_list[i]-pos_list[j])
    return e_tot

def ortho_basis(over_mat):
        #diagonalisation of S
        eig, u = np.linalg.eig(over_mat)
        #creation of matrix x via canonical orthogonalisation
        x = u / (eig) ** 0.5
        return x
