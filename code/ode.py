from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# parameters
mu_larvae  = 0.1
mu_adult = 0.1
fecundity = 0.9 # haha can't call this lambda
kk = 9./7. # carrying capacity (to be replaced with spatio-temporal data)
d = 0.1 # ageing rate

# mosquitoes array X[age, sex, genotype, species]
# start with 2 ages (could vary)
# 2 sexes, 3 genotypes (static, I hope!)
# 1 species (could vary)
n_ages = 2
n_species = 1
n_classes = 6 * n_ages * n_species

# inheritance cube: 
# i[x,y,z] = probability of mother genotype x, father genotype y
# producing offspring genotype z
# where index 0, 1, 2 = ww, Gw, GG respectively
# TODO: define as constant
i = np.array([[[1., 0., 0.],# ww x ww
[ .5, .5, 0. ],
[ 0., 1., 0. ]],
[[ .5, .5, 0. ],
[ .25, .5, .25 ], # Gw x Gw
[ 0., .25, .25 ]],
[[ 0., 1., 0. ],
[ 0., .5, .5 ],
[ 0., 0., 1. ]]]) # GG x GG

# calculate fecundity allocated to male and female larvae
# p[z,x,y] = proportion of total fecundity of mother genotype x, 
# father genotype y to producing offspring sex z
# where index of z: 0, 1 = male, female
# note - proportions will not add to 1 for male bias

accuracy = 0.95
fprop = 1./accuracy - 1.
prow = np.array([.5,.5,.5])
p = np.array([[prow, prow, prow],[prow, fprop*prow, fprop*prow]])

# pre-calculate i * p for males and females
ipm = np.zeros((3,3,3))
ipf = np.zeros((3,3,3))

for j in range(3):
    ipm[:,:,j] = i[:,:,j] * p[0,:,:]
    ipf[:,:,j] = i[:,:,j] * p[1,:,:]

# transpose y and z for use in later matrix multiplication
# (useful form to have in case of implicit)
for j in range(3):
    ipm[j,:,:] = np.transpose(ipm[j,:,:])
    ipf[j,:,:] = np.transpose(ipf[j,:,:])

#X = np.ones(n_classes)
X = np.ones(n_classes)*0.01
X[0:2] = 0.1


def f(t, y):
    Y = np.reshape(y, (2,2,3,1)) # example
    n = Y[:-1,:,:,:].sum() # total number of larvae (all but the last age class)
    ratio = Y[-1,0,:,:] # ratio of adult males (last age class) of given genotype and species
    ratio = ratio / ratio.sum()
    mat = np.zeros((n_classes, n_classes))
    # for each father's genotype and species, add the relevant proportion of fecundity for 
    # mothers of each genotype and that species, producing proportions of offspring 
    # of the relevant genotype (and the same species)
    for i in range(3): # fathers' genotypes
        for j in range(n_species): # species
            mat[j:(3*n_species):n_species, (9*(n_ages-1)*n_species + j):n_classes:n_species] += ipm[i,:,:] * ratio[i,j]
            mat[(3*n_species + j):(6*n_species):n_species, (9*(n_ages-1)*n_species + j):n_classes:n_species] += ipf[i,:,:] * ratio[i,j]
    mat *= (1 - n/kk)*fecundity # scaling by fecundity and density dependence
    # mortality
    mat[range(6 * (n_ages-1) * n_species), range(6 * (n_ages-1) * n_species)] = -mu_larvae - d
    mat[range(6 * (n_ages-1) * n_species, n_classes), range(6 * (n_ages-1) * n_species, n_classes)] = -mu_adult
    # ageing
    for i in range(n_ages - 1):
        mat[range(6*(i+1)*n_species, 6*(i+2)*n_species), range(6*i*n_species, 6*(i+1)*n_species)] = d
    dXdt = mat.dot(y)
    return dXdt

# ODE solver
sol = solve_ivp(f, [0, 10], X, t_eval = np.arange(0, 10.00001, 0.1))

# plot output for each class
# order is (currently one species)
# 0: age0 male ww species1
# 1: age0 male Gw species1
# 2: age0 male GG species1
# 3: age0 female ww species1
# 4: age0 female Gw species1
# 5: age0 female GG species1
# 6: age1 male ww species1
# 7: age1 male Gw species1
# 8: age1 male GG species1
# 9: age1 female ww species1
# 10: age1 female Gw species1
# 11: age1 female GG species1
fig, ax = plt.subplots()
for i in range(n_classes):
    ax.plot(sol.t, sol.y[i,:], color = ((i%3)/2.2, ((i//3)%2)*0.9, (i//6)*0.9), label = str(i))
ax.legend()
plt.show()
