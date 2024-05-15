import matplotlib.pyplot as plt

# Data points
x = [5, 15, 22.5]
y_N2_extinction = [9.30, 9.27, 9.28]
y_N2_cardelli = [9.34, 9.30, 9.31]

# Uncertainties
y_err_extinction = [0.17, 0.16, 0.18]
y_err_cardelli = [0.19, 0.17, 0.17]

plt.errorbar(x, y_N2_extinction, yerr=y_err_extinction, fmt='o', capsize=10,
             alpha=0.5, c='k', label='Fitzpatrick (2007) law', capthick=3)
plt.errorbar(x, y_N2_cardelli, yerr=y_err_cardelli, fmt='o', capsize=10,
             alpha=0.5, c='green', label='Cardelli law', capthick=3)

plt.plot(x, y_N2_extinction, marker='o', c='k', ls='--')
plt.plot(x, y_N2_cardelli, marker='o', c='green', ls='--')

plt.xlabel('Radius (kpc)')
plt.ylabel('12 + log(O/H)')

plt.legend()
plt.show()

#%%
plt.rcParams["figure.figsize"] = (12,8)
plt.rcParams.update({'font.size': 18})
# Data points
x = [5, 15, 22.5]
y_NIIOII = [9.30, 9.27, 9.28]
y_OIIINII = [8.73, 8.72, 8.75]
y_R23 = [9.17, 9.12, 9.07]

# Uncertainties
y_err_NIIOII = [0.17, 0.16, 0.18]
y_err_OIIINII = [0.17, 0.18, 0.15]
y_err_R23 = [0.03, 0.08, 0.07]

plt.errorbar(x, y_NIIOII, yerr=y_err_NIIOII, fmt='o', capsize=10,
             alpha=0.5, c='blue', label='NIIOII', capthick=3)
plt.errorbar(x, y_OIIINII, yerr=y_err_OIIINII, fmt='o', capsize=10,
             alpha=0.5, c='orange', label='OIIINII', capthick=3)
plt.errorbar(x, y_R23, yerr=y_err_R23, fmt='o', capsize=10,
             alpha=0.5, c='green', label='R23', capthick=3)

plt.plot(x, y_NIIOII, marker='o', c='blue', ls='--')
plt.plot(x, y_OIIINII, marker='o', c='orange', ls='--')
plt.plot(x, y_R23, marker='o', c='green', ls='--')

plt.xlabel('Radius (kpc)')
plt.ylabel('12 + log(O/H)')

plt.legend(loc='upper center', bbox_to_anchor=(0.85, 1))
plt.show()