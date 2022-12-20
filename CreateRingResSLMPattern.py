import IMGpy
import numpy as np

# This script is used to create SLM pattern for a ring focal pattern with user designed reservoir for rearrangement
# Currently this function will not be incorporated to our UI interface
# This script will generate an array of non-zero element corresponding to a user defined focal pattern
# This focal pattern will be launched to the SLM focal plane according to predetermined offset from the origin

# Load SLM parameter
from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication

#if __name__ == '__main__':

#    import sys
#    app = QApplication(sys.argv)
#    LS = LoadAndSave()
#    print("Now we load SLM config:")
#    SLMconfig = LS.LoadConfigFileDialog()
#    print(SLMconfig)

# Load SLM config file
file_name = 'SLMConfig/SLM_config_ring_130_384433.npy'
SLMconfig = np.load(file_name, allow_pickle=True).item()

pixelpitch = SLMconfig['pixel pitch']
#spacingx = SLMconfig['arr spacing x']
#spacingy = SLMconfig['arr spacing y']
distance = SLMconfig['distance from origin']
#lattice = SLMconfig['lattice geometry']
#arraysizex = SLMconfig['array size x']
#arraysizey = SLMconfig['array size y']
#arraysizeBit = SLMconfig['FFT grid size (bit)']
arraysizeBit = 13
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
    #rotAngle = SLMconfig['rotAngle']
     rotAngle = 38.4433

# Generate ring pattern
# Define ring pattern geometry,unit micron. Here our objective focal length is 500 mm
ring_spacing = 130*1e-6
ring_points = 32
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

#print(ringcoor_arr)
# Add atom reservoir for rearrangement, reservoir will be a rectangular array
# Reservoir start point, rotate from [m,n] by 45 degree
theta_start = -np.pi/4
m_res_start = round((m-ring_centerX)*np.cos(theta_start) - (n-ring_centerY)*np.sin(theta_start) + ring_centerX)
n_res_start = round((m-ring_centerX)*np.sin(theta_start) + (n-ring_centerY)*np.cos(theta_start) + ring_centerY)
theta_end = -np.pi/4 - np.pi
m_res_end = round((m-ring_centerX)*np.cos(theta_end) - (n-ring_centerY)*np.sin(theta_end) + ring_centerX)
n_res_end = round((m-ring_centerX)*np.sin(theta_end) + (n-ring_centerY)*np.cos(theta_end) + ring_centerY)
vecm = m_res_end - m_res_start
vecn = n_res_end - n_res_start
uvecm = vecm/(vecm**2 + vecn**2)**0.5
uvecn = vecn/(vecm**2 + vecn**2)**0.5
# uniform spacing
spacing = ring_spacing / Focalpitch
#print(uvecn)
# Get a vector that is pependicular to [uvecm, uvecn]
# spacing to match ring spacing coordinates
spacingxList = np.zeros(int(round(ring_points/4)))
num_pointx = int(round(ring_points/4))
angle1 = (np.pi- angle)/2
for index in range(num_pointx):
        angle_tmp = angle1 - (np.pi/2 - (np.pi/4 - index*angle))
        spacingxList[index] = spacing * np.cos(angle_tmp)
print(spacingxList)
spacingyList = np.flip(spacingxList[0:int(round(num_pointx/2))])
print(spacingyList)
vertList = np.zeros(np.size(spacingyList))
vert_tmp = 0
for index in range(np.size(vertList)):
    vert_tmp = vert_tmp + spacingyList[index]
    vertList[index] = vert_tmp
print(vertList)
vecnp = 1
vecmp = -uvecn*vecnp/uvecm
uvecmp = vecmp/(vecmp**2 + vecnp**2)**0.5
uvecnp = vecnp/(vecmp**2 + vecnp**2)**0.5
# Get the starting point
m_rec_start = m_res_start + (2*radius - np.sqrt(2)*radius)*0.5*uvecm
n_rec_start = n_res_start + (2*radius - np.sqrt(2)*radius)*0.5*uvecn
reccoor_arr_center = np.zeros([int(round(ring_points/4)+1), 2])

for index in range(int(round(ring_points/4)+1)):
    if index == 0:
        reccoor_arr_center[index,:] = [m_rec_start, n_rec_start]
        spacing_increase = 0
    else:
        spacing_increase = spacing_increase + spacingxList[index-1]
        m_rec_next = m_rec_start + spacing_increase*uvecm
        n_rec_next = n_rec_start + spacing_increase*uvecn
        reccoor_arr_center[index, :] = [m_rec_next, n_rec_next]

#reccoor_arr = np.concatenate((reccoor_arr_top2, reccoor_arr_top1, reccoor_arr_center, reccoor_arr_bottom1, reccoor_arr_bottom2),axis=0)
#reccoor_arr = reccoor_arr_center
#print(reccoor_arr)

# Here write a function to to populate other vertical rows
def CreateResVertical(vert, num_vert, uvec, vertlist, res_center):
     # vert is number + vert means increase the y coordinate by amount of vert from center
     # num_vert: number of atoms in this vertical location
     # uvec = [uvecmp, uvecnp] is the perpendicular vector
     # rec_center is the center coordinate
     # vertlist is the vertical coordinate list
     # find center coordinate
     [r, c] = np.shape(res_center)
     midres = (r-1)/2
     resstart = midres - (num_vert-1)/2
     res = res_center[int(round(resstart)):int(round(resstart))+num_vert,:]
     res_vert = np.zeros(np.shape(res))
     [r,c] = np.shape(res)
     for index in range(r):
         res_row = res[index,:]
         print(res_row)
         print(uvec)
         if vert > 0:
           res_tmp = res_row + uvec*vertlist[int(round(vert-1))]
         else:
           res_tmp = res_row - uvec * vertlist[int(round(-vert - 1))]
         res_vert[index,:] = res_tmp
     return res_vert

# Use the above function to create vertical row
uvec = np.array([uvecmp, uvecnp])
res_vert1 = CreateResVertical(1,9,uvec,vertList,reccoor_arr_center)
res_vert2 = CreateResVertical(2,7,uvec,vertList,reccoor_arr_center)
res_vert3 = CreateResVertical(3,5,uvec,vertList,reccoor_arr_center)
res_vert4 = CreateResVertical(4,3,uvec,vertList,reccoor_arr_center)
res_vertm1 = CreateResVertical(-1,9,uvec,vertList,reccoor_arr_center)
res_vertm2 = CreateResVertical(-2,7,uvec,vertList,reccoor_arr_center)
res_vertm3 = CreateResVertical(-3,5,uvec,vertList,reccoor_arr_center)
res_vertm4 = CreateResVertical(-4,3,uvec,vertList,reccoor_arr_center)
#print(res_vert1)
# Create coordinate lists for reservoir atoms
reccoor_arr = np.concatenate((res_vert4, res_vert3, res_vert2, res_vert1, reccoor_arr_center, res_vertm1, res_vertm2, res_vertm3, res_vertm4),axis=0)
print(reccoor_arr)

# Create targetAmp with nonzero element list
myRingPattern = IMGpy.arbFocalPattern([arraysizeBit, arraysizeBit])
# combine ring target with reservoir
target_arr = np.concatenate((ringcoor_arr, reccoor_arr),axis=0)
myRingTarget = myRingPattern.assembleFPfromNonZeroElement(target_arr)
# Get ring location
spacing = ring_spacing / pixelpitch
startRow_display = int(n-3*spacing-2*radius)
endRow_display = int(ImgResY/2 + 0.5*radius+3*spacing)
startCol_display = int(m-3*spacing-2*radius)
endCol_display = int(ImgResX/2 + 0.5*radius)
location = [startRow_display, endRow_display, startCol_display, endCol_display]
location_ring = location
# Create SLM pattern for the ring
myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                  wavelength, maskradius, res)
gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)
# Calculate phase pattern with WGS
targetAmp = myRingTarget
# rotate the target amp for non-zero rotation angle
if rot:
    targetAmp_rot, location_rot = myIMG.rotate_targetAmp(targetAmp, rotAngle, location, Plot=True)
    location = location_rot
    targetAmp = targetAmp_rot
myIMG.plotFocalplane(targetAmp, location)
Loop = 5
threshold = 0.02
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
# Save the SLM phase file and prepare for adapted WGS
save = 1
saveConfig = 1
if __name__ == '__main__':

     import sys
     app = QApplication(sys.argv)
     if save:
         LS = LoadAndSave()
         LS.SaveFileDialog(SLM_Screen_WGS)
         LS.SavePhaseFileDialog(SLM_Phase)
         LS.SaveTargetAmpFileDialog(targetAmp)
     if saveConfig:
         ConfigFile = {}
         ConfigFile["pixel pitch"] = pixelpitch
         #ConfigFile["arr spacing x"] = spacingx
         #ConfigFile["arr spacing y"] = spacingy
         ConfigFile["distance from origin"] = distance
         ConfigFile["lattice geometry"] = "ring"
         ConfigFile["lattice coordinate"] = target_arr
         ConfigFile["ring location"] = location_ring
         ConfigFile["ring location IMG"] = locationIMG
         #ConfigFile["array size x"] = arraysizex
         #ConfigFile["array size y"] = arraysizey
         ConfigFile["FFT grid size (bit)"] = arraysizeBit
         ConfigFile["WGS threshold"] = threshold
         ConfigFile["Loop"] = Loop
         ConfigFile["Focal length"] = focallength
         ConfigFile["Magnification"] = magnification
         ConfigFile["wave length"] = wavelength
         ConfigFile["beam waist"] = beamwaist
         ConfigFile["aperture radius"] = maskradius
         ConfigFile["use aperture?"] = mask
         ConfigFile["SLM resX"] = resX
         ConfigFile["SLM resY"] = resY
         ConfigFile["rot"] = rot
         if rot:
             ConfigFile["rotAngle"] = rotAngle
         LS = LoadAndSave()
         LS.SaveConfigFileDialog(ConfigFile)
     print("All finished!")
