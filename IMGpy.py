from typing import Union

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


class IMG:
    def __init__(self,pixelpitch,arraySizeBit,beamwaist,focallength,magnification, wavelength,maskdiameter):
        # pixelpitch = pixel size on the SLM plane
        # arraySizeBit =[size in x, size in y] in bit
        # beamwaist is the input beam waist measured at the SLM plane
        # focallength  = focal length of objective
        # wavelength is the wavelength of tweezer beam
        # maskdiameter is the aperture size in front of the SLM
        self.pixelpitch = pixelpitch
        self.arraySizeBitx = arraySizeBit[0]
        self.arraySizeBity = arraySizeBit[1]
        self.beamwaist = beamwaist
        self.fobj = focallength
        self.magnification = magnification
        self.wavelength = wavelength
        self.maskdiameter = maskdiameter
        self.ImgResX = 2 ** (int(self.arraySizeBitx))
        self.ImgResY = 2 ** (int(self.arraySizeBity))
        self.Xps, self.Yps = np.meshgrid(np.linspace(0, self.ImgResX, self.ImgResX), np.linspace(0, self.ImgResY, self.ImgResY))
        # define unit cell on the focal plane
        self.Focalpitchx = self.wavelength*self.fobj/self.ImgResX/(self.pixelpitch*self.magnification)
        self.Focalpitchy = self.wavelength * self.fobj / self.ImgResY /(self.pixelpitch*self.magnification)

    def initSLMImage(self, mask = True, Plot = True):
        # set the initial Gaussian amplitude confined by a circular aperture
        # set the initial phase to be a random number between [-pi,pi]
        # the normalization is done such that sum(element^2) = 1
        X = self.Xps*self.pixelpitch-self.ImgResX/2*self.pixelpitch
        Y = self.Yps*self.pixelpitch-self.ImgResY/2*self.pixelpitch
        if mask:
            maskAmp = (X**2+Y**2 <= self.maskdiameter)**2*1
            initGaussianAmp = np.multiply(np.sqrt(2/np.pi)/self.beamwaist*np.exp(-(X**2+Y**2)/self.beamwaist**2), maskAmp)
            initGaussianPhase = np.multiply(np.random.rand(self.ImgResX, self.ImgResY)*2*np.pi-np.pi, maskAmp)
        else:
            initGaussianAmp = np.sqrt(2 / np.pi) / self.beamwaist * np.exp(-(X ** 2 + Y ** 2)/self.beamwaist ** 2)
            initGaussianPhase = np.random.rand(self.ImgResX, self.ImgResY)*2*np.pi-np.pi
        initIntensity = np.square(initGaussianAmp)
        initGaussianAmpNorm = initGaussianAmp / np.sqrt(np.sum(initIntensity))
        if Plot:
            self.plotSLMplane(initGaussianAmpNorm)
        return initGaussianAmpNorm, initGaussianPhase

    def plotSLMplane(self, SLM_Field):
        c = plt.pcolormesh(self.Xps, self.Yps, SLM_Field, cmap='Greens', vmin=np.min(SLM_Field),
                           vmax=np.max(SLM_Field))
        plt.colorbar(c)
        plt.title("Field initialization on the SLM plane")
        plt.show()

    def plotFocalplane(self, Focal_Field, location):
        # location =[startRow, endRow, startCol, endCol]
        startRow = location[0]-5
        endRow = location[1]
        startCol = location[2]-5
        endCol = location[3]
        X = self.Xps * self.Focalpitchx * 1e6
        Y = self.Yps * self.Focalpitchy * 1e6
        initIntensityRegion = np.zeros((int(self.ImgResY), int(self.ImgResX)))
        initIntensityRegion[startRow:endRow, :][:, startCol:endCol] = 1
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=True)
        ax1.set_aspect(1)
        ax1.pcolormesh(X[startRow:endRow, :][:, startCol:endCol], Y[startRow:endRow, :][:, startCol:endCol],
                       Focal_Field[startRow:endRow, :][:, startCol:endCol], cmap='Greens', vmin=np.min(Focal_Field),
                       vmax=np.max(Focal_Field))
        ax1.set_title("Focal plane beam amplitude")
        ax1.set_xlabel('x (\u03bcm)')
        ax1.set_ylabel('y (\u03bcm)')

        ax2.set_aspect(1)
        im = ax2.pcolormesh(self.Xps, self.Yps,
                            initIntensityRegion, cmap='Greens', vmin=0, vmax=np.max(Focal_Field))
        ax2.set_title("Tweezer array location")
        fig.colorbar(im, orientation="vertical", pad=-0.6, shrink=0.6)
        fig.tight_layout()
        plt.show()

    def initFocalImage_RecLattice(self, distance, spacing, arraysize, Plot = True):
        # starting by setting the intensity landscape on the focal plane as square lattice offset from the origin
        # This distance is the real physical distance between the nearest point and the origin
        # This spacing is the real physical distance between the lattice spot
        # arraysize =[arraysizex, arraysizey], defines the target array we want on the SLM plane
        # for now, I am using homogeneous trap intensity for all traps. In the end, I need to set them at different
        # intensity to compensate for the possible non-uniform diffraction efficiency
        targetAmp = np.zeros((int(self.ImgResY), int(self.ImgResX)))
        dm = round(distance/np.sqrt(2)/self.Focalpitchx)
        dn = round(distance/np.sqrt(2)/self.Focalpitchy)
        mcenter = self.ImgResX/2
        ncenter = self.ImgResY/2
        m = mcenter-dm
        n = ncenter - dn
        totalsitesnum = arraysize[0]*arraysize[1]
        intensityPerSite = 1/totalsitesnum
        Trap = Tweezer([m,n], intensityPerSite)
        spacingx = round(spacing[0]/self.Focalpitchx)
        spacingy = round(spacing[1]/self.Focalpitchy)
        startRow, endRow, startCol, endCol = Trap.assembleRecLattice(arraysize, [spacingx, spacingy])
        targetAmp[startRow:endRow:spacingy,:][:, startCol:endCol:spacingx]=intensityPerSite**0.5

        if Plot:
            self.plotFocalplane(targetAmp, [startRow, endRow, startCol, endCol])
        location = [startRow, endRow, startCol, endCol]
        return self.Focalpitchx,self.Focalpitchy, targetAmp, location


class Tweezer:
    def __init__(self,location, intensity):
        # Here the location is chosen for the point nearest to the origin (zero order), I choose this point to locate at
        # the bottom left region
        self.m = location[0]
        self.n = location[1]
        self.intensity = intensity

    def assembleRecLattice(self, arraysize, spacing):
        arraysizex = arraysize[0]
        arraysizey = arraysize[1]
        spacingx = spacing[0]
        spacingy = spacing[1]
        startRow = self.n-arraysizey*spacingy
        endRow = self.n
        startCol = self.m-arraysizex*spacingx
        endCol = self.m
        if startRow <0 or startCol < 0:
            raise Exception("Sorry, too big an array, SLM cannot handle, consider shrinking the spacing or the size!")
        return int(startRow), int(endRow), int(startCol), int(endCol)


class WGS:
    def __init__(self, initGaussianAmp, initGaussianPhase, targetAmp):
        self.initGaussianAmp = initGaussianAmp
        self.initGaussianPhase = initGaussianPhase
        self.initGaussian = np.multiply(self.initGaussianAmp, np.exp(1j*self.initGaussianPhase))
        self.targetAmp = targetAmp
        self.targetAmpmask = (self.targetAmp>0)*1
        self.totalsites = np.count_nonzero(self.targetAmp)


    def fftLoop(self, Loop, threshold, Plot=True):
        SLM_Field = self.initGaussian
        count = 0
        g_coeff0 = 1
        non_uniform = np.zeros(Loop)
        while count < Loop:
            fftSLM = sp.fft.fft2(SLM_Field)
            fftSLMShift = sp.fft.fftshift(fftSLM)
            fftSLM_norm = np.sqrt(np.sum(np.square(np.abs(fftSLMShift))))
            fftSLMShift_norm = fftSLMShift/fftSLM_norm
            fftAmp = np.abs(fftSLMShift_norm)
            fftAmp_foci = np.multiply(fftAmp, self.targetAmpmask)
            non_uniform[count] = self.nonUniformity(fftAmp_foci)
            totalsites = np.count_nonzero(self.targetAmp)
            fftAmp_foci_avg = np.multiply(np.sum(fftAmp_foci) / totalsites, self.targetAmpmask)
            g_coeff = np.multiply(np.divide(fftAmp_foci_avg, fftAmp_foci, out=np.zeros_like(fftAmp_foci_avg),
                                          where=fftAmp_foci != 0),g_coeff0)
            Focal_Amp = np.multiply(self.targetAmp, g_coeff)
            if non_uniform[count]/np.max(non_uniform) > threshold:
               Focal_phase0 = np.angle(fftSLMShift_norm)
            else:
               Focal_phase0 = Focal_phase
            Focal_phase = Focal_phase0
            Focal_Field = np.multiply(Focal_Amp, np.exp(1j*Focal_phase))
            SLM_Field = sp.fft.ifft2(sp.fft.ifftshift(Focal_Field))
            SLM_Phase = np.angle(SLM_Field)
            SLM_Amp = self.initGaussianAmp
            SLM_Field = np.multiply(SLM_Amp, np.exp(1j*SLM_Phase))
            g_coeff0 = g_coeff
            count += 1

        SLM_Amp = np.abs(SLM_Field)
        SLM_Phase = np.angle(SLM_Field)
        if Plot:
            new_non_uniform = np.delete(non_uniform,0)
            new_non_uniform_norm = new_non_uniform/np.max(new_non_uniform)
            plt.plot(new_non_uniform_norm)
            plt.grid()
            plt.yscale('log')
            plt.xlabel('Iteration')
            plt.ylabel('Non-uniformity')
            plt.show()
        return SLM_Amp, SLM_Phase, fftAmp, new_non_uniform_norm

    def nonUniformity(self, Amp_foci):
        # This function calculate the uniformity of the tweezer traps at the focal point
        # Amp_foci is the field amplitude at the focal plane. It is a matrix
        Inten_foci = np.square(Amp_foci)/np.sum(np.square(Amp_foci))
        Inten_foci_avg = np.multiply(np.sum(Inten_foci)/self.totalsites, self.targetAmpmask)
        non_Uniform = np.sqrt((np.sum(np.square(Inten_foci-Inten_foci_avg)))/self.totalsites)
        return non_Uniform

    def SLM_IMG(self, SLM_Phase, SLMResX, SLMResY, Plot=True, Save=False):
        # This function is to crop the center pixel area to fit to the SLM screen
        # SLM_Phase is the phase image calculated by WGS
        # SLMResX, SLMResY retrieve the size of the SLM screen
        # The final result will be converted into an 8 bit image
        # X is column
        col = np.size(SLM_Phase, axis=0)
        # Y is row
        row = np.size(SLM_Phase, axis=1)
        centerX = col/2
        centerY = row/2
        startCol = centerX-round(SLMResX/2)
        endCol = centerX+round(SLMResX/2)
        startRow = centerY-round(SLMResY/2)
        endRow = centerY+round(SLMResY/2)
        SLM_IMG = SLM_Phase[int(startRow):int(endRow), :][:, int(startCol):int(endCol)]
       # location = [startRow, endRow, startCol, endCol]
        SLM_bit = self.phaseTobit(SLM_IMG)
        if Plot:
            colIMG = np.size(SLM_IMG, axis=0)
            rowIMG = np.size(SLM_IMG, axis=1)
            # For meshgrid function, X is row, Y is column
            X, Y = np.meshgrid(np.linspace(0, rowIMG, rowIMG), np.linspace(0, colIMG, colIMG))
            c = plt.pcolormesh(X, Y, SLM_bit, cmap='Greens', vmin=np.min(SLM_bit),
                              vmax=np.max(SLM_bit))
            plt.colorbar(c)
            plt.title("Physical SLM plane")
            plt.show()
        if Save:
            np.savetxt("SLM.csv", SLM_bit, delimiter=",")
        return SLM_bit

    def phaseTobit(self, SLM_IMG):
        # This function change a phase number between -pi to pi to an integer between 0 to 255
        SLM_bit = np.around((SLM_IMG+np.pi)/(2*np.pi)*255).astype('uint8')
        return SLM_bit
