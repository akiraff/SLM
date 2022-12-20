import IMGpy
import numpy as np

# This script is used to perform adapted WGS for ring pattern with user designed reservoir.

# Load SLM parameter
from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication

# Load SLM config file
file_name = 'SLMConfig/SLM_config_ring_130_384433.npy'
SLMconfig = np.load(file_name, allow_pickle=True).item()

pixelpitch = SLMconfig['pixel pitch']
#spacingx = SLMconfig['arr spacing x']
#spacingy = SLMconfig['arr spacing y']
distance = SLMconfig['distance from origin']
target_arr = SLMconfig["lattice coordinate"]
print(target_arr)
location = SLMconfig["ring location"]
locationIMG = SLMconfig["ring location IMG"]
#lattice = SLMconfig['lattice geometry']
#arraysizex = SLMconfig['array size x']
#arraysizey = SLMconfig['array size y']
arraysizeBit = SLMconfig['FFT grid size (bit)']
focallength = SLMconfig['Focal length']
magnification = SLMconfig['Magnification']
wavelength = SLMconfig['wave length']
beamwaist = SLMconfig['beam waist']
maskradius = SLMconfig['aperture radius']
mask = SLMconfig['use aperture?']
resX = SLMconfig['SLM resX']
resY = SLMconfig['SLM resY']
res = np.min([resX, resY])
print("resX =: ", resX)
print("resY=: ", resY)
rot = SLMconfig['rot']
if rot:
    rotAngle = SLMconfig['rotAngle']
print(rotAngle)

# Load measured intensity, SLM phase, targetAmp from last iteration
if __name__ == '__main__':

     import sys
     app = QApplication(sys.argv)
     LS = LoadAndSave()
     print("Load measured intensity array:")
     intenArray = LS.LoadIntensityFileDialog()
     print("intenArray:", intenArray)
     print("Now we load phase pattern from the previous iteration:")
     SLM_Phase = LS.LoadPhaseFileDialog()
     print("Now we load target Amp from last iteration:")
     targetAmp_adapt_lastiter = LS.LoadTargetAmpFileDialog()

# Create targetAmp_foci from intenArray
myRingPattern = IMGpy.arbFocalPattern([arraysizeBit, arraysizeBit])
myRingTarget = myRingPattern.assembleFPfromNonZeroElementAdapt(target_arr, intenArray)
targetAmp_foci = myRingTarget
myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                  wavelength, maskradius, res)
gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)
if rot:
    targetAmp_foci_rot, location_foci_rot = myIMG.rotate_targetAmp(targetAmp_foci, rotAngle, location, Plot=True)
    location = location_foci_rot
    targetAmp_foci = targetAmp_foci_rot

# raise exception if the user loads the wrong file
try:
   sz = np.size(intenArray)
   [r,c] = np.shape(target_arr)
   if sz != r:
        raise ValueError("Input intensity array does not match the array size!")
except ValueError as ve:
 print(ve)

try:
    sz = np.size(SLM_Phase)
    if sz != np.size(targetAmp_foci):
        raise ValueError("Input SLM phase matrix does not match its target amplitude matrix, you might "
                            "mistakenly load SLM screen file.")
except ValueError as ve:
 print(ve)

 try:
     sz = np.size(targetAmp_adapt_lastiter)
     if sz != np.size(targetAmp_foci):
         raise ValueError("target Amp from last iteration does not match current target Amp, you load the "
                          "wrong the file!")
 except ValueError as ve:
     print(ve)

# Calculate adapted WGS
Loop = 15
threshold = 0.02
WGScal = IMGpy.WGS(gaussianAmp, gaussianPhase, targetAmp_foci)
SLM_Amp_adapt, SLM_Phase_adapt, Focal_Amp_adapt, non_uniform_adapt, targetAmp_adapt = WGScal.fftLoop_adapt(
               SLM_Phase, targetAmp_foci, targetAmp_adapt_lastiter, Loop, threshold)
myIMG.plotFocalplane(targetAmp_adapt, location)
myIMG.plotFocalplane(Focal_Amp_adapt, location)
SLM_bit, fftSLM_IMG_Amp_norm, SLM_Screen_WGS = WGScal.SLM_IMG(SLM_Phase_adapt, resX, resY, Plot=True)
myIMG.plotFocalplane(fftSLM_IMG_Amp_norm, locationIMG)
print("location = ", location)
print("locationIMG = ", locationIMG)
print("Phase image shape:", SLM_Screen_WGS.shape)

# Save the SLM phase file and prepare for adapted WGS
save = 1
if save:
   LS = LoadAndSave()
   LS.SaveFileDialog(SLM_Screen_WGS)
   LS.SavePhaseFileDialog(SLM_Phase_adapt)
   LS.SaveTargetAmpFileDialog(targetAmp_adapt)
print("All finished!")