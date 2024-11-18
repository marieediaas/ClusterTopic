hist, bin_edges = np.histogram(df_['Fe_H'], bins=30)

# gaussiana dupla
def double_gaussian(x, a1, mu1, sigma1, a2, mu2, sigma2):
    return (a1 * np.exp(-(x - mu1)*2 / (2 * sigma1*2)) / (sigma1 * np.sqrt(2 * np.pi)) +
            a2 * np.exp(-(x - mu2)*2 / (2 * sigma2*2)) / (sigma2 * np.sqrt(2 * np.pi)))

# pega o bin central da coluna desejada F_e
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# paremetros
initial_guess = [1, np.mean(df_['Fe_H']), np.std(df_['Fe_H']), 1, np.mean(df_['Fe_H']), np.std(df_['Fe_H'])]

# ajusta a gaussiana 
params, _ = curve_fit(double_gaussian, bin_centers, hist, p0=initial_guess)

# hist
plt.hist(df_['Fe_H'], bins=30, alpha=0.8)

# 2 * gaussiana
x_range = np.linspace(min(df_['Fe_H']), max(df_['Fe_H']), 100)
plt.plot(x_range, double_gaussian(x_range, *params), 'r-')

plt.text(-0.56,11, f'$\mu_1$={round(params[1],2)}, $\sigma_1={round(params[2],2)}$')

plt.text(-0.56,10.5, f'$\mu_2$={round(params[4],2)}, $\sigma_2={round(params[5],2)}$')

plt.legend()
plt.xlabel('Fe/H')
plt.ylabel('Frequencia')

print("mu1:", params[1])
print("mu2:", params[4])
print("sigma1:", params[2])
print("sigma2:", params[5])
