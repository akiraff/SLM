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
arraysizeBit = 11
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
ring_spacing = 105*1e-6
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
        theta = 2*np.pi / ring_points * index
        m_next = round((m-ring_centerX)*np.cos(theta) - (n-ring_centerY)*np.sin(theta) + ring_centerX)
        n_next = round((m-ring_centerX)*np.sin(theta) + (n-ring_centerY)*np.cos(theta) + ring_centerY)
        ringcoor_arr[index, :] = [m_next, n_next]

print(ringcoor_arr)
# uniform spacing
spacing = ring_spacing / Focalpitch


# generate coordinate value from ring coordinates
xcoor = np.flip(np.sort(np.unique(ringcoor_arr[:,0])))
ycoor = np.flip(np.sort(np.unique(ringcoor_arr[:,1])))
#print('xcoor = :', xcoor)
#print("ycoor = :", ycoor)
# Generate a column of traps outside the ring
def CreateResCol(xlocs, ylocs_arr):
    resCol = np.zeros([np.size(ylocs_arr),2])
    for j in range(np.size(ylocs_arr)):
        resCol[j,:] = [xlocs, ylocs_arr[j]]
    return resCol
# Generate a row of traps outside the ring
def CreateResRow(xlocs_arr, ylocs):
    resRow = np.zeros([np.size(xlocs_arr), 2])
    for j in range(np.size(xlocs_arr)):
        resRow[j, :] = np.array([xlocs_arr[j], ylocs])
    return resRow

# Create reservoir inside the ring
xlocs_arr = xcoor[5:9:2]  # 9 is not included
ylocs = int(ycoor[3])
resRow1 = CreateResRow(xlocs_arr, ylocs)

xlocs_arr = xcoor[4:10:2]
ylocs = int(ycoor[4])
resRow2 = CreateResRow(xlocs_arr, ylocs)

xlocs_arr = xcoor[3:11:2]
ylocs = int(ycoor[5])
resRow3 = CreateResRow(xlocs_arr, ylocs)

xlocs_arr = xcoor[4:10:2]
ylocs = int(ycoor[6])
resRow4 = CreateResRow(xlocs_arr, ylocs)

xlocs_arr = xcoor[3:11:2]
ylocs = int(ycoor[7])
resRow5 = CreateResRow(xlocs_arr, ylocs)

xlocs_arr = xcoor[4:10:2]
ylocs = int(ycoor[8])
resRow6 = CreateResRow(xlocs_arr, ylocs)

xlocs_arr = xcoor[5:9:2]  # 9 is not included
ylocs = int(ycoor[9])
resRow7 = CreateResRow(xlocs_arr, ylocs)

resin_arr = np.concatenate((resRow1, resRow2, resRow3, resRow4, resRow5, resRow6, resRow7), axis=0)

# Create resveroir outside the ring
xlocs = int(xcoor[-1] - round(spacing*1.5))
ylocs_arr = ycoor[4:10:2]
#print("ycoor = :", ycoor)
#print("ycoor site = :", ylocs_arr)
resCol1 = CreateResCol(xlocs, ylocs_arr)
xlocs = int(xcoor[-1] - round(spacing*1.5 + spacing*1.2))
ylocs_arr = ycoor[4:10:2]
resCol2 = CreateResCol(xlocs, ylocs_arr)
#print("resCol1 = :", resCol1)
#print("resCol2 = :", resCol2)
xlocs = int(xcoor[-1] - round(spacing*1.5 + spacing*1.2*2))
ylocs_arr = ycoor[4:10:2]
resCol3 = CreateResCol(xlocs, ylocs_arr)


# Create reservoir row
xlocs_arr = xcoor[3:11:2]
ylocs = int(ycoor[0] + round(spacing*1.5))
resRow1 = CreateResRow(xlocs_arr, ylocs)
xlocs_arr = xcoor[3:11:2]
ylocs = int(ycoor[0] + round(spacing*1.5 + spacing*1.2))
resRow2 = CreateResRow(xlocs_arr, ylocs)
xlocs_arr = xcoor[3:11:2]
ylocs = int(ycoor[0] + round(spacing*1.5 + spacing*1.2*2))
resRow3 = CreateResRow(xlocs_arr, ylocs)
xlocs_arr = xcoor[[3, 9]]
ylocs = int(ycoor[0] + round(spacing*1.5 + spacing*1.2*3))
resRow4 = CreateResRow(xlocs_arr, ylocs)


xlocs_arr = xcoor[4:10:2]
ylocs = int(ycoor[-1] - round(spacing*1.5))
resRow5 = CreateResRow(xlocs_arr, ylocs)
xlocs_arr = xcoor[4:10:2]
ylocs = int(ycoor[-1] - round(spacing*1.5 + spacing*1.2))
resRow6 = CreateResRow(xlocs_arr, ylocs)
xlocs_arr = xcoor[4:10:2]
ylocs = int(ycoor[-1] - round(spacing*1.5 + spacing*1.2*2))
resRow7 = CreateResRow(xlocs_arr, ylocs)

# Create bottom row
xlocs = int(xcoor[0] + round(spacing*1.5))
ylocs_arr = ycoor[[5,7]]
#print("ycoor = :", ycoor)
#print("ycoor site = :", ylocs_arr)
resCol4 = CreateResCol(xlocs, ylocs_arr)
xlocs = int(xcoor[0] + round(spacing*1.5 + spacing*1.2))
ylocs_arr = ycoor[[5,7]]
resCol5 = CreateResCol(xlocs, ylocs_arr)
xlocs = int(xcoor[0] + round(spacing*1.5 + spacing*1.2*2))
ylocs_arr = ycoor[[5,7]]
resCol6 = CreateResCol(xlocs, ylocs_arr)

resout_arr = np.concatenate((resCol1, resCol2, resCol3, resRow1, resRow2, resRow3, resRow4, resRow5, resRow6, resRow7, resCol4, resCol5, resCol6), axis=0)
print(resout_arr)

# Create targetAmp with nonzero element list
myRingPattern = IMGpy.arbFocalPattern([arraysizeBit, arraysizeBit])
# combine ring target with reservoir
target_arr = np.concatenate((ringcoor_arr, resin_arr, resout_arr),axis=0)
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
Loop = 4
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
save = 0
saveConfig = 0
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