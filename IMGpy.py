import numpy as np
import scipy as sp
import matplotlib.pyplot as plt


class IMG:
    def __init__(self,pixelpitch,arraySizeBit,beamwaist,focallength,magnification, wavelength, maskradius, SLMRes):
        # pixelpitch = pixel size on the SLM plane
        # arraySizeBit =[size in x, size in y] in bit
        # beamwaist is the input beam waist measured at the SLM plane
        # focallength  = focal length of objective
        # wavelength is the wavelength of tweezer beam
        # maskradius is the aperture size in front of the SLM
        # SLMres is the min screen resolution (pixel number) of the SLM
        self.pixelpitch = pixelpitch
        self.arraySizeBitx = arraySizeBit[0]
        self.arraySizeBity = arraySizeBit[1]
        self.beamwaist = beamwaist
        self.fobj = focallength
        self.magnification = magnification
        self.wavelength = wavelength
        self.maskradius = maskradius
        self.ImgResX = 2 ** (int(self.arraySizeBitx))
        self.ImgResY = 2 ** (int(self.arraySizeBity))
        self.SLMRes = SLMRes
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
            maskAmp = (X**2+Y**2 <= self.maskradius**2)*1
            initGaussianAmp = np.multiply(np.sqrt(2/np.pi)/self.beamwaist*np.exp(-(X**2+Y**2)/self.beamwaist**2), maskAmp)
            initGaussianPhase = np.multiply(np.random.rand(self.ImgResX, self.ImgResY)*2*np.pi-np.pi, maskAmp)
        else:
            initGaussianAmp = np.sqrt(2/np.pi)/self.beamwaist*np.exp(-(X**2+Y**2)/self.beamwaist**2)
            initGaussianPhase = np.random.rand(self.ImgResX, self.ImgResY)*2*np.pi-np.pi
        initIntensity = np.square(initGaussianAmp)
        initGaussianAmpNorm = initGaussianAmp / np.sqrt(np.sum(initIntensity))
        if Plot:
            self.plotSLMplane(initGaussianAmpNorm)
        return initGaussianAmpNorm, initGaussianPhase

    def plotSLMplane(self, SLM_Field):
        c = plt.pcolormesh(self.Xps, self.Yps, SLM_Field, cmap='Greens', vmin=np.min(SLM_Field),
                           vmax=np.max(SLM_Field), shading='auto')
        plt.colorbar(c)
        plt.title("Field initialization on the SLM plane")
        plt.show()

    def plotFocalplane(self, Focal_Field, location):
        # location =[startRow, endRow, startCol, endCol]
        startRow = location[0]
        endRow = location[1]
        startCol = location[2]
        endCol = location[3]
        # X is column
        colField = np.size(Focal_Field, axis=1)
        focalx = self.Focalpitchx*self.ImgResX/colField
        # Y is row
        rowField = np.size(Focal_Field, axis=0)
        focaly = self.Focalpitchx*self.ImgResY/rowField
        X = self.Xps * focalx * 1e6
        Y = self.Yps * focaly * 1e6
        initIntensityRegion = np.zeros((int(self.ImgResY), int(self.ImgResX)))
        initIntensityRegion[startRow:endRow, :][:, startCol:endCol] = 1
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9, 3.5), constrained_layout=False)
        ax1.set_aspect(1)
        ax1.pcolormesh(X[startRow:endRow, :][:, startCol:endCol], Y[startRow:endRow, :][:, startCol:endCol],
                       Focal_Field[startRow:endRow, :][:, startCol:endCol], cmap='Greens', vmin=np.min(Focal_Field),
                       vmax=np.max(Focal_Field), shading='auto')
        ax1.set_title("Focal plane beam amplitude")
        ax1.set_xlabel('x (\u03bcm)')
        ax1.set_ylabel('y (\u03bcm)')

        ax2.set_aspect(1)
        im = ax2.pcolormesh(self.Xps, self.Yps,
                            initIntensityRegion, cmap='Greens', vmin=0, vmax=np.max(Focal_Field), shading='auto')
        ax2.set_title("Tweezer array location")
        fig.colorbar(im, orientation="vertical", pad=0.1, shrink=0.6)
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
        spacingx = round(spacing[0]/self.Focalpitchx)
        spacingy = round(spacing[1]/self.Focalpitchy)
        print([spacingx,spacingy])
        print(arraysize)
        Trap = Tweezer([m, n],arraysize)
        startRow, endRow, startCol, endCol = Trap.assembleRecLattice([spacingx, spacingy])
        targetAmp[startRow:endRow:spacingy,:][:, startCol:endCol:spacingx]=intensityPerSite**0.5
        startRow_display = int(startRow-spacingy)
        endRow_display = int(self.ImgResY/2)
        startCol_display = int(startCol-spacingx)
        endCol_display = int(self.ImgResX/2)
        location = [startRow_display, endRow_display, startCol_display, endCol_display]
        if Plot:
            self.plotFocalplane(targetAmp, location)
        return self.Focalpitchx,self.Focalpitchy, targetAmp, location

    def initFocalImage_KagomeLattice(self, distance, spacing, arraysize, Triangle = False, Plot = True):
        # This function is used to generate the Kagome geometry in the focal plane
        # Parameters definition are the same as Rec_lattice
        targetAmp = np.zeros((int(self.ImgResY), int(self.ImgResX)))
        targetAmpMask = np.zeros((int(self.ImgResY), int(self.ImgResX)))
        dm = round(distance / np.sqrt(2) / self.Focalpitchx)
        dn = round(distance / np.sqrt(2) / self.Focalpitchy)
        mcenter = self.ImgResX / 2
        ncenter = self.ImgResY / 2
        m = mcenter - dm
        n = ncenter - dn
        totalsitesnum = arraysize[0] * arraysize[1]
        intensityPerSite = 1 / totalsitesnum
        spacingx = round(spacing[0] / self.Focalpitchx)
        spacingy = round(spacing[1] / self.Focalpitchy)
        print([spacingx, spacingy])
        print(arraysize)
        # Build Kagome cell
        Trap = Tweezer([m, n], arraysize)
        p_vector, unitcell, num_cell_h, num_cell_v = Trap.unitcell_Kagome([spacingx, spacingy], Triangle)
        print(p_vector)
        Kagome_cell = Trap.assembleLatticeFromUnitcell([spacingx, spacingy], p_vector, unitcell, num_cell_h, num_cell_v)
       # print(Kagome_cell)
        # Launch Kagome cell onto targetAmp
        for cell_value in Kagome_cell.values():
            site_loc = cell_value[0:2]
            col = site_loc[0]
            row = site_loc[1]
            #targetAmp[row, col] = cell_value[2]
            targetAmpMask[row, col] = 1
        TraplocMask = np.argwhere(targetAmpMask == 1)
        TraplocMask_order = TraplocMask[TraplocMask[:,0].argsort()]
       # print(TraplocMask_order)
        TraplocMask_firstRow = TraplocMask_order[-arraysize[1]:,:]
       # print(TraplocMask_firstRow)
        TraplocMask_firstRow_orderbyCol = TraplocMask_firstRow[TraplocMask_firstRow[:,1].argsort()][:, 1]
       # print(TraplocMask_firstRow_orderbyCol)
        col_min = TraplocMask_firstRow_orderbyCol[0]
        col_spacing = np.abs(TraplocMask_firstRow_orderbyCol[1]-TraplocMask_firstRow_orderbyCol[0])

        for cell_value in Kagome_cell.values():
            site_loc = cell_value[0:2]
            col = site_loc[0]
            row = site_loc[1]
            if col < col_min:
                col = int(col + arraysize[0]*col_spacing)
            targetAmp[row, col] = cell_value[2]

        Traploc = np.argwhere(targetAmp == 1)
        Traprow = Traploc[:, 0]
        Trapcol = Traploc[:, 1]
        startRow_display = int(np.min(Traprow) - spacingy)
        endRow_display = int(self.ImgResY / 2 + spacingy)
        startCol_display = int(np.min(Trapcol) - spacingx)
        endCol_display = int(self.ImgResX / 2 + spacingx)
        location = [startRow_display, endRow_display, startCol_display, endCol_display]
        if Plot:
            self.plotFocalplane(targetAmp, location)
        return self.Focalpitchx,self.Focalpitchy, targetAmp, location

    def diffraction_efficiency(self, location):
        # This function calculates the diffraction efficiency at different location according to the distance from origin
        # Diffraction efficiency = sinc^2(pi/2*(dtrap/dmax)), dmax=fobj/magnification*wavelength/pixelpitch
        # location =[m,n]=[column, row]
        dmax = self.fobj/self.magnification*self.wavelength/self.pixelpitch/2
        mcenter = self.ImgResX / 2
        ncenter = self.ImgResY / 2
        d = np.sqrt(((location[0]-mcenter)*self.Focalpitchx)**2+((location[1]-ncenter)*self.Focalpitchy)**2)
        diffrac_efficiency = (np.sinc(np.pi/2*d/dmax))**2
        return diffrac_efficiency

    def modify_targetAmp(self, targetAmp, location, Plot = True):
        # This function receives an intensity pattern with uniform intensity on different lattice sites, it will
        # output a target intensity pattern taking into account the finite diffraction efficiency.
        # Currently, this function uses the theoretical value, will update by measurements if necessary.
        # X is column
        col = np.size(targetAmp, axis=1)
        # Y is row
        row = np.size(targetAmp, axis=0)
        diffracEff = np.zeros_like(targetAmp)
        for i in range(col):
            for j in range(row):
                diffrac = self.diffraction_efficiency([i,j])
                diffracEff[j][i] = diffrac
        targetAmp_diffrac = np.divide(targetAmp,diffracEff)
        if Plot:
            self.plotFocalplane(targetAmp_diffrac, location)
        return targetAmp_diffrac

class Tweezer:
    def __init__(self,location,arraysize):
        # Here the location is chosen for the point nearest to the origin (zero order), I choose this point to locate at
        # the bottom left region, spacing in this class counts the matrix spacing, not the real spacing in focal plane
        self.m = location[0]
        self.n = location[1]
        self.arraysizex = arraysize[0]
        self.arraysizey = arraysize[1]


    def assembleRecLattice(self, spacing):
        spacingx=spacing[0]
        spacingy=spacing[1]
        startRow = self.n-self.arraysizey*spacingy
        endRow = self.n
        startCol = self.m-self.arraysizex*spacingx
        endCol = self.m
        if startRow < 0 or startCol < 0:
            raise Exception("Sorry, too big an array, SLM cannot handle, consider shrinking the spacing or the size!")
        return int(startRow), int(endRow), int(startCol), int(endCol)

    def assembleLatticeFromUnitcell(self, spacing, p_vector, unitcell, num_cell_h, num_cell_v):
        # This function assembles a lattice from unit cell using primitive vector
        spacingh = spacing[0]
        spacingv = spacing[1]
        xtranslate = int(round(self.arraysizex / num_cell_h - 1))
        ytranslate = int(round(self.arraysizey / num_cell_v - 1))
        cell = {}
        for site_info, site_value in unitcell.items():
            site_loc = site_value[0:2]
            site_inten = site_value[2]
            # Build unit cell
            cell['unit cell: ' + site_info] = site_loc
            cell['unit cell: ' + site_info] = np.append(cell['unit cell: ' + site_info], site_inten)
            for j in range(xtranslate):
                cell['unit cell: ' + site_info + ' cell: ' + str(2 * j + 1)] = np.round(
                    site_loc + p_vector['eh'] * (j + 1) * num_cell_h * spacingh).astype(int)
                cell['unit cell: ' + site_info + ' cell: ' + str(2 * j + 1)] = np.append(
                    cell['unit cell: ' + site_info + ' cell: ' + str(2 * j + 1)], site_inten)
            for j in range(ytranslate):
                cell['unit cell: ' + site_info + ' cell: ' + str(2 * j)] = np.round(
                    site_loc + p_vector['ev'] * (j + 1) * num_cell_v * spacingv).astype(int)
                cell['unit cell: ' + site_info + ' cell: ' + str(2 * j)] = np.append(
                    cell['unit cell: ' + site_info + ' cell: ' + str(2 * j)], site_inten)
                for k in range(xtranslate):
                    cell['unit cell: ' + site_info + ' cell: ' + str(2 * j) + ' Expand Horizontally: ' + str(k)] = \
                        np.round(site_loc + p_vector['ev'] * (j + 1) * num_cell_v * spacingv
                                 + p_vector['eh'] * (k + 1) * num_cell_h * spacingh).astype(int)

                    cell['unit cell: ' + site_info + ' cell: ' + str(2 * j) + ' Expand Horizontally: ' + str(k)] = \
                        np.append(cell['unit cell: ' + site_info + ' cell: ' + str(2 * j) + ' Expand Horizontally: '
                                       + str(k)], site_inten)
        return cell


    def unitcell_Kagome(self, spacing, Triangle = False):
        # This function defines the unit cell for Kagome geometry, it will return the unit cell corrdinate as well as
        # primitive vector.
        # unit cell ['site index'] = [X, Y] = [Col, Row]
        spacingh = spacing[0]
        spacingv = spacing[1]
        p_vector = {}
        p_vector['ev'] = np.array([-1/2, -np.sqrt(3)/2])
        p_vector['eh'] = np.array([-1, 0])

        # Number of sites horizontally and vertically per unit cell
        num_cell_h = 2
        num_cell_v = 2

        unitcell = {}
        unitcell['site 1'] = np.array([self.m, self.n]).astype(int)
        unitcell['site 2'] = (np.round(unitcell['site 1'] + spacingv * p_vector['ev'])).astype(int)
        unitcell['site 3'] = (np.round(unitcell['site 2'] + spacingh * p_vector['eh'])).astype(int)
        unitcell['site 4'] = (np.round(unitcell['site 1'] + spacingh * p_vector['eh'])).astype(int)
        # Append trap intensity as the last value in the array
        unitcell['site 1'] = np.append(unitcell['site 1'], 1)
        unitcell['site 2'] = np.append(unitcell['site 2'], 1)
        unitcell['site 3'] = np.append(unitcell['site 3'], 1)
        if Triangle:
            unitcell['site 4'] = np.append(unitcell['site 4'], 1)
        else:
            unitcell['site 4'] = np.append(unitcell['site 4'], 0)

        return p_vector, unitcell, num_cell_h, num_cell_v



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

    def SLM_IMG(self, SLM_Phase, SLMResX, SLMResY, Plot=True):
        # This function is to crop the center pixel area to fit to the SLM screen
        # SLM_Phase is the phase image calculated by WGS
        # SLMResX, SLMResY retrieve the size of the SLM screen
        # The final result will be converted into an 8 bit image
        # X is column
        col = np.size(SLM_Phase, axis=1)
        # Y is row
        row = np.size(SLM_Phase, axis=0)
        centerX = col/2
        centerY = row/2
        SLMRes = np.min([SLMResX, SLMResY])
        startCol = centerX-round(SLMRes/2)
        endCol = centerX+round(SLMRes/2)
        startRow = centerY-round(SLMRes/2)
        endRow = centerY+round(SLMRes/2)
        SLM_IMG = SLM_Phase[int(startRow):int(endRow), :][:, int(startCol):int(endCol)]
        SLM_gaussianAmp = self.initGaussianAmp[int(startRow):int(endRow), :][:, int(startCol):int(endCol)]
        SLM_Field_final = np.multiply(SLM_gaussianAmp, np.exp(1j*SLM_IMG))
        fftSLM_IMG = sp.fft.fftshift(sp.fft.fft2(SLM_Field_final))
        fftSLM_IMG_normfactor = np.sqrt(np.sum(np.square(np.abs(fftSLM_IMG))))
        fftSLM_IMG_Amp_norm = np.abs(fftSLM_IMG/fftSLM_IMG_normfactor)
       # location = [startRow, endRow, startCol, endCol]
        #SLM_bit = self.phaseTobit(SLM_IMG)
        SLM_bit = SLM_IMG
        SLM_screen = np.zeros((int(SLMResY), int(SLMResX)))
        col_SLM_bit = np.size(SLM_bit, axis=1)
        row_SLM_bit = np.size(SLM_bit, axis=0)
        startRow_screen = SLMResY/2-round(row_SLM_bit/2)
        endRow_screen = SLMResY/2+round(row_SLM_bit/2)
        startCol_screen = SLMResX/2-round(col_SLM_bit/2)
        endCol_screen = SLMResX/2+round(col_SLM_bit/2)
        SLM_screen[int(startRow_screen):int(endRow_screen), :][:, int(startCol_screen):int(endCol_screen)] = SLM_bit
        if Plot:
            colIMG = np.size(SLM_screen, axis=1)
            rowIMG = np.size(SLM_screen, axis=0)
            # For meshgrid function, X is col, Y is row
            X, Y = np.meshgrid(np.linspace(0, colIMG, colIMG), np.linspace(0, rowIMG, rowIMG))
            c = plt.pcolormesh(X, Y, SLM_screen, cmap='Greens', vmin=np.min(SLM_screen),
                              vmax=np.max(SLM_screen), shading='auto')
            plt.colorbar(c)
            plt.title("Physical SLM plane")
            plt.show()
        '''
        if Save:
            np.savetxt("SLM.csv", SLM_screen, delimiter=",")
        '''
        return SLM_bit, fftSLM_IMG_Amp_norm, SLM_screen

    def phaseTobit(self, SLM_IMG):
        # This function change a phase number between -pi to pi to an integer between 0 to 255
        SLM_bit = np.around((SLM_IMG+np.pi)/(2*np.pi)*255).astype('uint8')
        return SLM_bit
