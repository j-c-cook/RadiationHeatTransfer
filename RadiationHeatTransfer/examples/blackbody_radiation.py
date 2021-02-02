# Jack C. Cook
# Sunday, January 31, 2021

"""
Plancks Law:
Implement Plancks law using Eb,lambda function
"""

import RadiationHeatTransfer as RHT
from scipy.optimize import fminbound
from itertools import count
from itertools import takewhile
import matplotlib.pyplot as plt
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d


def any_range(start, stop, step):
    start = type(start + step)(start)
    return takewhile(lambda n: n < stop, count(start, step))


def main():
    # Recreate Figure 12-9
    # Figure 12-9: The variation of the blackbody emissive power with wavelength for
    # several temperatures
    fig, ax = plt.subplots()

    several_temperatures: list = [100, 300, 500, 1000, 2000, 4000, 5800]  # Kelvin

    wavelengths: list = list(any_range(.01, 1000, .01))

    for _, T in enumerate(several_temperatures):
        # T = several_temperatures[-1]  # the current temperature

        Eblambdas: list = [RHT.blackbody.Eblambda(wavelength, T) for _, wavelength in enumerate(wavelengths)]

        # Find the maximum value over the range of 0.01-1000 micrometers
        # Source: https://stackoverflow.com/a/16781456/11637415
        max_lmbda = fminbound(lambda lmbda: -RHT.blackbody.Eblambda(lmbda, T), 0.01, 1000)
        max_Eblmbda = RHT.blackbody.Eblambda(max_lmbda, T)

        ax.plot(wavelengths, Eblambdas, color='k')
        ax.scatter(max_lmbda, max_Eblmbda, color='r')
        ax.annotate('{0} K ($\lambda={1:.2f}$)'.format(int(T), max_lmbda), xy=(max_lmbda, max_Eblmbda), xytext=(max_lmbda*10, max_Eblmbda),
                    arrowprops=dict(arrowstyle='->'))

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([10**-6, 10**9])

    ax.grid(which='both')
    ax.set_axisbelow(True)

    ax.set_xlabel('Wavelength $\lambda$, $\mu$m')
    ax.set_ylabel('E$_{b \lambda}$, W/m$^2 \cdot \mu$m')

    fig.savefig('blackbody_emissive_power.pdf')
    plt.close(fig)

    # use the Eb function to integrate over the whole wavelength to determine the total black body emissive power
    T = 5800
    integral, _ = RHT.blackbody.Eb(T)
    sigma = integral / T**4
    print('Stephan Boltzman constant, sigma={0:.2E}'.format(sigma))

    # Pick a few points to make sure the 0-lambda is working right
    T = 5000
    lmbda = 0.2
    print('lambda * T = {}'.format(T * lmbda))
    print('f_lambda = {0:.6}'.format(RHT.blackbody.f_lambda(lmbda, T)))
    lmbda = 0.4
    print('lambda * T = {}'.format(T * lmbda))
    print('f_lambda = {0:.6}'.format(RHT.blackbody.f_lambda(lmbda, T)))
    lmbda = 1.
    print('lambda * T = {}'.format(T * lmbda))
    print('f_lambda = {0:.6}'.format(RHT.blackbody.f_lambda(lmbda, T)))

    # Use the effective spectral function to find an effective value
    lmbdas_nm: list = [0.1, 250, 625, 1200, 2200, 2250, 2500, 1000000]  # lambda values in nanometers
    alphas: list = [0.92, 0.92, 0.92, 0.39, 0.39, 0.47, 0.47, 0.47]
    lmbdas_mm: list = [x / 10**3 for _, x in enumerate(lmbdas_nm)]  # 1000 nm = 1micrometer
    # create a linear 1D interpolation function for the absorptivity
    f = interp1d(lmbdas_mm, alphas, kind='linear')

    T_sun: float = 5800
    eff, wavelengths, eff_Eblambdas, Eblambdas = RHT.blackbody.effective_spectral(lmbdas_mm[0], lmbdas_mm[-1], T_sun, f)
    print('The effective absorptivity: {0:.5f}'.format(eff))

    # Plot the absorptivity function
    fig, ax = plt.subplots()

    ax.plot(lmbdas_mm, alphas, color='blue')

    ax.set_ylabel(r'$\alpha$')
    ax.set_xlabel('$\lambda$ ($\mu$m)')

    ax.set_xlim([-0.2, 3])
    ax.set_ylim([0, 1])

    fig.savefig('absorptivity.pdf')
    plt.close(fig)

    # Plot the black body emissive power and the effective black body emissive power
    fig, ax = plt.subplots()

    ax.plot(wavelengths, Eblambdas, label='E$_{b \lambda}$')
    ax.plot(wavelengths, eff_Eblambdas, '--', label=r'$\alpha \cdot$ E$_{b \lambda}$')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([10 ** -6, 10 ** 9])
    ax.set_xlim([10**-2, 10**3])

    fig.legend()

    ax.grid(which='both')
    ax.set_axisbelow(True)

    ax.set_xlabel('Wavelength $\lambda$, $\mu$m')
    ax.set_ylabel('E$_{b \lambda}$, W/m$^2 \cdot \mu$m')

    fig.savefig('alpha_lambda.pdf')
    plt.close(fig)


if __name__ == '__main__':
    main()
