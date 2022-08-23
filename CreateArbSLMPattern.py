import IMGpy
from file_dialog_widget import LoadAndSave
import numpy as np

# This script is used to create SLM pattern for an arbitrary focal pattern set by the user
# Currently this function will not be incorporated to our UI interface
# This script will generate an array of non-zero element corresponding to a user defined focal pattern
# This focal pattern will be launched to the SLM focal plane according to predetermined offset from the origin

# Load SLM parameter
from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication

if __name__ == '__main__':

    import sys
    app = QApplication(sys.argv)
    LS = LoadAndSave()
    print("Now we load SLM config:")
    SLMconfig = LS.LoadConfigFileDialog()
    print(SLMconfig)

pixelpitch = SLMconfig['pixel pitch']
#spacingx = SLMconfig['arr spacing x']
#spacingy = SLMconfig['arr spacing y']
distance = SLMconfig['distance from origin']
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

# Generate ring pattern
# Define ring pattern geometry,unit micron. Here our objective focal length is 500 mm
ring_spacing = 120*1e-6
ring_points = 40
angle = 2*np.pi/ring_points
ring_radius = ring_spacing/2/(np.sin(angle/2))
# First point of the circle determined by the offset
dm = round(distance / np.sqrt(2) / pixelpitch)
dn = round(distance / np.sqrt(2) / pixelpitch)
ImgResX = 2 ** (int(arraysizeBit))
ImgResY = 2 ** (int(arraysizeBit))
Focalpitch = wavelength*focallength/ImgResX/(pixelpitch*magnification)
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
# Initialize and populate numpy array with coordinates of points on a circle
ringcoor_arr = np.zeros([ring_points, 2])
for index in range(ring_points):
    if index == 0:
        ringcoor_arr[index,:] = [m, n]
    else:
        theta = 2*np.pi / ring_points * index
        m_next = round((m-ring_centerX)*np.cos(theta) - (n-ring_centerY)*np.sin(theta) + ring_centerX)
        n_next = round((m-ring_centerX)*np.sin(theta) + (n-ring_centerY)*np.cos(theta) + ring_centerY)
        ringcoor_arr[index, :] = [m_next, n_next]
# Create targetAmp with nonzero element list
#print(ringcoor_arr)
myRingPattern = IMGpy.arbFocalPattern([arraysizeBit, arraysizeBit])
myRingTarget = myRingPattern.assembleFPfromNonZeroElement(ringcoor_arr)
# Get ring location
spacing = ring_spacing / pixelpitch
startRow_display = int(n-spacing-2*radius)
endRow_display = int(ImgResY/2)
startCol_display = int(m-spacing-2*radius)
endCol_display = int(ImgResX/2)
location = [startRow_display, endRow_display, startCol_display, endCol_display]
# Create SLM pattern for the ring
myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                  wavelength, maskradius, res)
gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)
# Calculate phase pattern with WGS
targetAmp = myRingTarget
myIMG.plotFocalplane(targetAmp, location)
Loop = 10
threshold = 0.01
WGScal = IMGpy.WGS(gaussianAmp, gaussianPhase, targetAmp)
SLM_Amp, SLM_Phase, Focal_Amp, non_uniform = WGScal.fftLoop(Loop, threshold)
myIMG.plotFocalplane(Focal_Amp, location)
SLM_bit, fftSLM_IMG_Amp_norm, SLM_Screen_WGS = WGScal.SLM_IMG(SLM_Phase, resX, resY, Plot=True)
# X is column
colIMG = np.size(SLM_bit, axis=1)
colFFT = np.size(SLM_Phase, axis=1)
# Y is row
rowIMG = np.size(SLM_bit, axis=0)
rowFFT = np.size(SLM_Phase, axis=0)
endRowIMG = np.max([rowIMG/2, (location[1]-rowFFT/2)*rowIMG/rowFFT + rowIMG/2])
startRowIMG = round(endRowIMG-(location[1]-location[0])*rowIMG/rowFFT)
endColIMG = np.max([colIMG/2, (location[3]-colFFT/2)*colIMG/colFFT + colIMG/2])
startColIMG = round(endColIMG-(location[3]-location[2])*colIMG/colFFT)
locationIMG = [int(startRowIMG), int(endRowIMG), int(startColIMG), int(endColIMG)]
myIMG.plotFocalplane(fftSLM_IMG_Amp_norm, locationIMG)
print("location = ", location)
print("locationIMG = ", locationIMG)
print("Phase image shape:", SLM_Screen_WGS.shape)