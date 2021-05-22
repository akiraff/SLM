import numpy as np
import matplotlib.pyplot as plt
from zernike import RZern
from scipy.special import eval_jacobi

class Zernike:
    def __init__(self, SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike, percent):
        # This class outputs the Zernike polynomial of given order for wavefront aberration correction
        # SLMResX, SLMResY retrieve the size of the SLM screen
        # pixelpitch = pixel size on the SLM plane
        # aperture_radius defines the aperture size of the imaging system
        # ind_Zernike is the index of the Zernike polynomial. Convention is from https://github.com/jacopoantonello/zernike
        # percent defines the level of the Zernike polynomial. Physically, it corresponds to the rms value of the
        # correspnding aberration
        self.SLMResX = SLMResX
        self.SLMResY = SLMResY
        self.pixelpitch = pixelpitch
        self.aperture_radius = aperture_radius
        self.ind_Zernike = ind_Zernike
        self.percent = percent

    def phase_Zernike(self, Plot = True, Save = False):
        SLMRes = np.min([self.SLMResX, self.SLMResY])
        apertureD = int(round(self.aperture_radius/self.pixelpitch*2))
        Zmatrix_size = np.min([SLMRes, apertureD])
        r0 = self.pixelpitch*Zmatrix_size/2
        r0mm = r0 * 1e3
        cart = RZern(6)
        ddx = np.linspace(-1.0, 1.0, Zmatrix_size)
        ddy = np.linspace(-1.0, 1.0, Zmatrix_size)
        xv, yv = np.meshgrid(ddx, ddy)
        cart.make_cart_grid(xv, yv)
        c = np.zeros(cart.nk)
        c[[self.ind_Zernike]] = 1
        Phi = cart.eval_grid(c, matrix=True)
        # change all nan values in Phi matrix to be zero
        where_are_NaNs = np.isnan(Phi)
        Phi[where_are_NaNs] = 0
        Phi_norm = Phi
        n = cart.ntab[self.ind_Zernike]
        m = cart.mtab[self.ind_Zernike]
        zernike_str = f"Radial: {n}, Angular: {m}"
        print(zernike_str)
        # Launch the Zernike phase pattern to SLM screen
        SLM_screen_aberr = np.zeros((int(self.SLMResY), int(self.SLMResX)))
        col_Phi_norm = np.size(Phi_norm, axis=1)
        row_Phi_norm = np.size(Phi_norm, axis=0)
        startRow_screen = self.SLMResY / 2 - round(row_Phi_norm / 2)
        endRow_screen = self.SLMResY / 2 + round(row_Phi_norm / 2)
        startCol_screen = self.SLMResX / 2 - round(col_Phi_norm / 2)
        endCol_screen = self.SLMResX / 2 + round(col_Phi_norm / 2)
        SLM_screen_aberr[int(startRow_screen):int(endRow_screen), :][:, int(startCol_screen):int(endCol_screen)] = \
            Phi_norm*self.percent
        if Plot:
            im = plt.imshow(Phi_norm, origin='lower', extent=(-r0mm, r0mm, -r0mm, r0mm))
            plt.title(zernike_str)
            plt.colorbar(im)
            plt.show()
        if Save:
            np.savetxt(f"SLM_Zernike_{n}{m}_{percent}.csv", SLM_screen_aberr, delimiter=",")
        return SLM_screen_aberr, m, n

    def radial_polynomial(self, rho, m, n):
        # See https://mathworld.wolfram.com/ZernikePolynomial.html
        if (n - m) % 2 == 1:
            print("ODD")
            return 0

        eff_n = (n - m) // 2
        alpha = m
        beta = 0
        x = 1 - 2 * (rho ** 2.0)
        return (-1) ** (eff_n) * (rho ** m) * eval_jacobi(eff_n, alpha, beta, x)

    def zernike(self, rho, phi, m, n):
        # Zernike polynomials are normalized (according to Wikipedia: OSA/ANSI) such that the integral of the
        # polynomial squared over the unit disk is equal to Pi. Here we enforce this normalization.
        # See opt.indiana.edu/vsg/library/vsia/vsia-2000_taskforce/tops4_2.html
        normalization = np.sqrt(2 * (n + 1))
        if m == 0:
            normalization /= np.sqrt(2)

        if m == 0:
            return normalization * self.radial_polynomial(rho, m, n)
        elif m > 0:
            return normalization * self.radial_polynomial(rho, m, n) * np.cos(m * phi)
        elif m < 0:
            return normalization * self.radial_polynomial(rho, -m, n) * np.sin(m * phi)

    def phase_Zernike_continuous(self, m, n, Plot = True):
        # This function will calculate continuous phase pattern with Zernike polynomials
        # This phase pattern will be continuous over a square area of SLM screen
        # When constraining to unit circle, these are typical Zernike polynomials
        SLMRes = np.min([self.SLMResX, self.SLMResY])
        xs = (np.arange(SLMRes) - SLMRes // 2) / (SLMRes // 2)
        ys = (np.arange(SLMRes) - SLMRes // 2) / (SLMRes // 2)
        xx, yy = np.meshgrid(xs, ys)

        rhos = np.sqrt(xx ** 2.0 + yy ** 2.0)
        phis = np.arctan2(yy, xx)

        phase_profile = self.zernike(rhos, phis, m, n)
       # phase_profile[rhos >= 1.0] = 0.0
        #print(np.max(phase_profile))
        # Launch phase_profile onto SLM screen
        SLM_screen_aberr = np.zeros((int(self.SLMResY), int(self.SLMResX)))
        col_phase_profile = np.size(phase_profile, axis=1)
        row_phase_profile = np.size(phase_profile, axis=0)
        startRow_screen = self.SLMResY / 2 - round(row_phase_profile / 2)
        endRow_screen = self.SLMResY / 2 + round(row_phase_profile / 2)
        startCol_screen = self.SLMResX / 2 - round(col_phase_profile / 2)
        endCol_screen = self.SLMResX / 2 + round(col_phase_profile / 2)
        SLM_screen_aberr[int(startRow_screen):int(endRow_screen), :][:, int(startCol_screen):int(endCol_screen)] = \
            phase_profile * self.percent

        if Plot:
            c = plt.imshow(phase_profile, cmap='Greens', vmin=np.min(phase_profile), vmax=np.max(phase_profile))
            plt.colorbar(c)
            plt.title("n=%d, m=%d: " % (n, m))
            plt.show()
        return SLM_screen_aberr


"""""""""
SLMResX = 1272
SLMResY = 1024
pixelpitch = 12.5e-6
aperture_radius = 6e-3
ind_Zernike = 8
percent = 1
myOberrationCorr = Zernike(SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike, percent)
SLM_screen_aberr, m, n = myOberrationCorr.phase_Zernike(Plot = True, Save = False)
SLM_screen_aberr_continuous = myOberrationCorr.phase_Zernike_continuous(m, n)
print("Phase image shape:", SLM_screen_aberr.shape)
print("Phase image shape:", SLM_screen_aberr_continuous.shape)
print(np.min(SLM_screen_aberr_continuous))
print(np.min(SLM_screen_aberr))

a = np.array([[9, 2, 3],
           [4, 5, 6],
           [7, 0, 5]])
a=a[a[:,1].argsort()]
print(a)
"""""""""