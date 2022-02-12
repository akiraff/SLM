from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication
import numpy as np
from Aberration import Zernike

if __name__ == '__main__':

    import sys
    app = QApplication(sys.argv)
    LS = LoadAndSave()
    print("Now we load SLM config for our previous iteration:")
    SLMconfig = LS.LoadConfigFileDialog()
    print(SLMconfig)
    SLMResX = SLMconfig['SLM resX']
    SLMResY = SLMconfig['SLM resY']
    pixelpitch = SLMconfig['pixel pitch']
    aperture_radius = SLMconfig['aperture radius']
    print("Next, we load the SLM screen phase file:")
    SLM_screen_WGS = LS.LoadFileDialog()
    try:
        sz = np.size(SLM_screen_WGS)
        if sz != SLMResX * SLMResY:
            raise ValueError("Input SLM phase matrix does not match its SLMConfig, you might "
                             "mistakenly load SLM phase file.")
    except ValueError as ve:
        print(ve)
    # Add predetermined aberration correction on top of the pattern, Z(2,-2): 0.3. ind =4, discontinous
    ind_Zernike = 4
    percent = 0.3
    myOberrationCorr = Zernike(SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike, percent)
    SLM_aberr_screen, m, n = myOberrationCorr.phase_Zernike(Plot=True, Save=False)
    # convert total phase between -pi and pi
    SLM_screen_WGS_aberr_add = SLM_screen_WGS + SLM_aberr_screen
    min_add = np.min(SLM_screen_WGS_aberr_add)
    # print(min_add)
    if min_add < 0:
        ind_2pi = np.abs(np.floor(np.divide(min_add, 2 * np.pi)))
        print(ind_2pi)
        SLM_screen_WGS_aberr_add = SLM_screen_WGS_aberr_add + ind_2pi * 2 * np.pi
    SLM_screen_WGS_aberr_mod = np.mod(SLM_screen_WGS_aberr_add, 2 * np.pi)
    SLM_screen_WGS_aberr = np.multiply((SLM_screen_WGS_aberr_mod <= np.pi), SLM_screen_WGS_aberr_mod) \
                           + np.multiply((SLM_screen_WGS_aberr_mod > np.pi), SLM_screen_WGS_aberr_mod - 2 * np.pi)

    # Add Z(4,0): 0.2 ind = 10, discontinous
    SLM_screen_WGS = SLM_screen_WGS_aberr
    ind_Zernike = 10
    percent = 0.2
    myOberrationCorr = Zernike(SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike, percent)
    SLM_aberr_screen, m, n = myOberrationCorr.phase_Zernike(Plot=True, Save=False)
    # convert total phase between -pi and pi
    SLM_screen_WGS_aberr_add = SLM_screen_WGS + SLM_aberr_screen
    min_add = np.min(SLM_screen_WGS_aberr_add)
    # print(min_add)
    if min_add < 0:
        ind_2pi = np.abs(np.floor(np.divide(min_add, 2 * np.pi)))
        print(ind_2pi)
        SLM_screen_WGS_aberr_add = SLM_screen_WGS_aberr_add + ind_2pi * 2 * np.pi
    SLM_screen_WGS_aberr_mod = np.mod(SLM_screen_WGS_aberr_add, 2 * np.pi)
    SLM_screen_WGS_aberr = np.multiply((SLM_screen_WGS_aberr_mod <= np.pi), SLM_screen_WGS_aberr_mod) \
                           + np.multiply((SLM_screen_WGS_aberr_mod > np.pi), SLM_screen_WGS_aberr_mod - 2 * np.pi)


    # Add multiple Zernike correction
    SLM_screen_WGS = SLM_screen_WGS_aberr
    LS = LoadAndSave()
    aberrConfigzind = LS.LoadAberrCorrConfigFileDialog()
    print(aberrConfigzind)
    ind_Zernike_list = list(aberrConfigzind.keys())
    percent_list = list(aberrConfigzind.values())
    SLM_aberr_screen_sum = 0
    for j in range(len(ind_Zernike_list)):
        myOberrationCorr = Zernike(SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike_list[j], percent_list[j])
        SLM_aberr_screen, m, n = myOberrationCorr.phase_Zernike(Plot=False, Save=False)
        SLM_aberr_screen = myOberrationCorr.phase_Zernike_continuous(m, n, Plot=False)
        SLM_aberr_screen_sum = SLM_aberr_screen_sum + SLM_aberr_screen
    SLM_screen_WGS_aberr_add = SLM_screen_WGS + SLM_aberr_screen_sum
    min_add = np.min(SLM_screen_WGS_aberr_add)
    # print(min_add)
    if min_add < 0:
        ind_2pi = np.abs(np.floor(np.divide(min_add, 2 * np.pi)))
        print(ind_2pi)
        SLM_screen_WGS_aberr_add = SLM_screen_WGS_aberr_add + ind_2pi * 2 * np.pi
    SLM_screen_WGS_aberr_mod = np.mod(SLM_screen_WGS_aberr_add, 2 * np.pi)
    SLM_screen_WGS_aberr = np.multiply((SLM_screen_WGS_aberr_mod <= np.pi), SLM_screen_WGS_aberr_mod) \
                           + np.multiply((SLM_screen_WGS_aberr_mod > np.pi),
                                         SLM_screen_WGS_aberr_mod - 2 * np.pi)
    print(np.max(SLM_screen_WGS_aberr))
    print(np.min(SLM_screen_WGS_aberr))
    print("Finally, we add the aberration and save the phase file.")
    LS.SaveFileDialog(SLM_screen_WGS_aberr)
    print("All finished!")