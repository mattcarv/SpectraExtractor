import numpy as np

z_12 = [9.21, 9.20, 9.05]

def MetaConv(z_list):
    solar_z_12 = 8.69
    solar_z = 0.0134
    
    z_sol = [(10**(z - solar_z_12)) * solar_z for z in z_list]
    
    return z_sol

z_sol = MetaConv(z_12)
print('Metallicities in solar units: ',z_sol)

def AlphaCO(z_sol):
    
    Sig_mstar = 5.8e5
    Sig_mmol = 2.3e5
    gamma = 0.5
    tot_mass = Sig_mstar + Sig_mmol
    
    alpha_corr = [2.9 * (np.exp((0.4/i) * ((tot_mass/100)**(-gamma)))) for i in z_sol]
    
    return alpha_corr


print('AlphaCO for Z_sol: ',AlphaCO(z_sol))

#but if they use Z` as the normalized Z by the solar value Z/Z_sol with Z_sol = 8.69

def MetaConv_2(z_list):
    solar_z_12 = 8.69
    
    z_sol = [z/solar_z_12 for z in z_list]
    
    return z_sol

z_sol_2 = MetaConv_2(z_12)
print('AlphaCO for Z_norm',AlphaCO(z_sol_2))