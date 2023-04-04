import IMGpy
import numpy as np

# This script is used to perform adapted WGS for ring pattern with user designed reservoir.
# This script performs spacing feedback on top of intensity feedback
# For simplicity, it deals with R2 intensity feedback

# Load SLM parameter
from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication

# Load SLM config file
file_name = 'SLMConfig/SLM_config_ring_24atom_109_thetaR1.npy'
SLMconfig = np.load(file_name, allow_pickle=True).item()

# Load measured theta
csvfile = 'RingPattern/onAtom_24atomRing109/theta_list_24atom_109_r1.csv'
theta_list = np.genfromtxt(csvfile, delimiter=',')
print(theta_list)

pixelpitch = SLMconfig['pixel pitch']
#spacingx = SLMconfig['arr spacing x']
#spacingy = SLMconfig['arr spacing y']
distance = SLMconfig['distance from origin']
target_arr = SLMconfig["lattice coordinate"]
#print(target_arr)
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
     print("Load measured intensity array R1:")
     intenArray_cur = LS.LoadIntensityFileDialog()
     print("intenArray R1:", intenArray_cur)
     print("Load measured intensity array R0:")
     intenArray_past = LS.LoadIntensityFileDialog()
     print("intenArray R0:", intenArray_past)
     print("Now we load phase pattern  R2:")
     SLM_Phase = LS.LoadPhaseFileDialog()
  #   print("Now we load target Amp from last iteration:")
  #   targetAmp_adapt_last = LS.LoadTargetAmpFileDialog()
print(intenArray_past)
# IntenArray_adapt is the targetAmp for R1 feedback
IntenArray_past_avg = np.mean(intenArray_past)
IntenArray_past_adapt = np.divide(IntenArray_past_avg , intenArray_past)
print(IntenArray_past_adapt)
# Generate target_arr from input theta_list
if rot:
    rotAngle = SLMconfig['rotAngle']
    print(rotAngle)
   #  rotAngle = 38.4433

# Generate ring pattern
# Define ring pattern geometry,unit micron. Here our objective focal length is 500 mm
ring_spacing = 109*1e-6
ring_points = 24
angle = 2*np.pi/ring_points
ring_radius = ring_spacing/2/(np.sin(angle/2))
# First point of the circle determined by the offset
#dm = round(distance*6.22193/12.5 / np.sqrt(2) / pixelpitch)
#dn = round(distance / np.sqrt(2) / pixelpitch)
print(pixelpitch)
ImgResX = 2 ** (int(arraysizeBit))
ImgResY = 2 ** (int(arraysizeBit))
Focalpitch = wavelength*focallength/ImgResX/(pixelpitch*magnification)
print(Focalpitch)
dm = round(distance*6.22193/12.5 / np.sqrt(2) / Focalpitch)
dn = round(distance*6.22193/12.5 / np.sqrt(2) / Focalpitch)
#Focalpitchy = wavelength*focallength/ImgResX/(pixelpitch*magnification)
mcenter = ImgResX / 2
ncenter = ImgResY / 2
m = mcenter - dm
n = ncenter - dn
# Get coordinate for the center of the circle
radius = ring_radius / Focalpitch
c_vecVal = np.sqrt(dm**2+dn**2)
c_vecX = -dm / c_vecVal
c_vecY = -dn / c_vecVal
ring_centerX = m + radius*c_vecX
ring_centerY = n + radius*c_vecY
print("ring center X: = ", ring_centerX)
print("ring_center Y: = ", ring_centerY)
# Initialize and populate numpy array with coordinates of points on a circle
ringcoor_arr = np.zeros([ring_points, 2])
for index in range(ring_points):
    if index == 0:
        ringcoor_arr[index,:] = [m, n]
    else:
       # theta = 2*np.pi / ring_points * index
        theta = theta_list[index - 1] * np.pi / 180
        m_next = round((m-ring_centerX)*np.cos(theta) - (n-ring_centerY)*np.sin(theta) + ring_centerX)
        n_next = round((m-ring_centerX)*np.sin(theta) + (n-ring_centerY)*np.cos(theta) + ring_centerY)
        ringcoor_arr[index, :] = [m_next, n_next]

print(ringcoor_arr)
#print(np.size(ringcoor_arr))
print(target_arr[:24, :])
#print(np.size(target_arr))
print((ringcoor_arr == target_arr[:24,:]).all())

# Replace the first 24 elements corresponding to newly measured theta_list
target_arr[:24, :] = ringcoor_arr
print(target_arr)

# Now generate target Amp according to intensity array
# Create targetAmp_foci from intenArray
myRingPattern = IMGpy.arbFocalPattern([arraysizeBit, arraysizeBit])
myRingTarget = myRingPattern.assembleFPfromNonZeroElementAdapt(target_arr, intenArray_cur)
targetAmp_foci = myRingTarget
myRingTarget_last = myRingPattern.assembleFPfromNonZeroElementAdapt(target_arr, IntenArray_past_adapt)
targetAmp_foci_last = myRingTarget_last
myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                  wavelength, maskradius, res)
gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)

if rot:
    targetAmp_foci_rot_last, location_foci_rot = myIMG.rotate_targetAmp(targetAmp_foci_last, rotAngle, location, Plot=True)
    location = location_foci_rot
    targetAmp_foci_last = targetAmp_foci_rot_last

    targetAmp_foci_rot, location_foci_rot = myIMG.rotate_targetAmp(targetAmp_foci, rotAngle, location, Plot=True)
    location = location_foci_rot
    targetAmp_foci = targetAmp_foci_rot



# raise exception if the user loads the wrong file
try:
   sz = np.size(intenArray_cur)
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
     sz = np.size(intenArray_past)
     [r, c] = np.shape(target_arr)
     if sz != r:
         raise ValueError("Input intensity array does not match the array size!")
 except ValueError as ve:
     print(ve)

# Calculate adapted WGS
Loop = 4
threshold = 0.02
WGScal = IMGpy.WGS(gaussianAmp, gaussianPhase, targetAmp_foci)
SLM_Amp_adapt, SLM_Phase_adapt, Focal_Amp_adapt, non_uniform_adapt, targetAmp_adapt = WGScal.fftLoop_adapt(
               SLM_Phase, targetAmp_foci, targetAmp_foci_last, Loop, threshold)
myIMG.plotFocalplane(targetAmp_adapt, location)
myIMG.plotFocalplane(Focal_Amp_adapt, location)
SLM_bit, fftSLM_IMG_Amp_norm, SLM_Screen_WGS = WGScal.SLM_IMG(SLM_Phase_adapt, resX, resY, Plot=True)
myIMG.plotFocalplane(fftSLM_IMG_Amp_norm, locationIMG)
print("location = ", location)
print("locationIMG = ", locationIMG)
print("Phase image shape:", SLM_Screen_WGS.shape)

# Save the SLM phase file and prepare for adapted WGS
save = 0
if save:
   LS = LoadAndSave()
   LS.SaveFileDialog(SLM_Screen_WGS)
   LS.SavePhaseFileDialog(SLM_Phase_adapt)
   LS.SaveTargetAmpFileDialog(targetAmp_adapt)
print("All finished!")

