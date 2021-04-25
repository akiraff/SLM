import numpy as np
import matplotlib.pyplot as plt
from zernike import RZern

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
        return SLM_screen_aberr

    def add_Zernke(self, SLM_screen, SLM_screen_aberr):
        SLM_screen_aberrCorr = SLM_screen + SLM_screen_aberr
        return SLM_screen_aberrCorr

"""""
SLMResX = 1272
SLMResY = 1024
pixelpitch = 12.5e-6
aperture_radius = 6.5e-3
ind_Zernike = 1
percent = 0.01
myOberrationCorr = Zernike(SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike, percent)
SLM_screen = myOberrationCorr.phase_Zernike(Plot = True, Save = True)
"""""
"""""
a = np.array([[9, 2, 3],
           [4, 5, 6],
           [7, 0, 5]])
a=a[a[:,1].argsort()]
print(a)
"""""

