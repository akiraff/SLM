import IMGpy
import numpy as np
# This script is used to create SLM pattern for a ring focal pattern with user designed reservoir for rearrangement
# Currently this function will not be incorporated to our UI interface
# This script will generate an array of non-zero element corresponding to a user defined focal pattern
# This focal pattern will be launched to the SLM focal plane according to predetermined offset from the origin

# Load SLM parameter
from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication

# Load SLM config file
file_name = 'SLMConfig/SLM_config_rec.npy'
SLMconfig = np.load(file_name, allow_pickle=True).item()
print(SLMconfig)

pixelpitch = SLMconfig['pixel pitch']
distance = SLMconfig['distance from origin']
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
    # rotAngle = 38.4433
    rotAngle = 0

# Generating rectangular array using arbitray pattern class, in case we want to change the geometry
ImgResX = 2 ** (int(arraysizeBit))
ImgResY = 2 ** (int(arraysizeBit))
Focalpitch = wavelength * focallength / ImgResX / (pixelpitch * magnification)
dm = round(distance * 6.22193 / 12.5 / np.sqrt(2) / Focalpitch)
dn = round(distance * 6.22193 / 12.5 / np.sqrt(2) / Focalpitch)
mcenter = ImgResX / 2
ncenter = ImgResY / 2
m = mcenter - dm
n = ncenter - dn
arraysizex = 17
arraysizey = 10

totalsitesnum = arraysizex * arraysizey
intensityPerSite = 1 / totalsitesnum
rec_spacing = 128*1e-6
spacing = round(rec_spacing/Focalpitch)
Trap = IMGpy.Tweezer([m, n],[arraysizex, arraysizey])
startRow, endRow, startCol, endCol = Trap.assembleRecLattice([spacing, spacing])
print([startRow, endRow, startCol, endCol])
rec_coor = np.zeros([totalsitesnum, 2])
row_coor = np.arange(endRow, startRow, -spacing)
print(row_coor)
col_coor = np.arange(endCol, startCol, -spacing)
print(col_coor)
count = 0
for ind_row in range(np.size(row_coor)):
    row_val = row_coor[ind_row]
    for ind_col in range(np.size(col_coor)):
        col_val = col_coor[ind_col]
        rec_coor[count, :] = [col_val, row_val]
        count = count + 1
#print(rec_coor)

# Create targetAmp with nonzero element list
myRecPattern = IMGpy.arbFocalPattern([arraysizeBit, arraysizeBit])
# combine ring target with reservoir
#target_arr = rec_coor[1:19,:]
target_arr = rec_coor
print(target_arr)
myRecTarget = myRecPattern.assembleFPfromNonZeroElement(target_arr)
# Get rec location
startRow_display = int(startRow-spacing)
endRow_display = int(ImgResY/2)
startCol_display = int(startCol-spacing)
endCol_display = int(ImgResX/2)
location = [startRow_display, endRow_display, startCol_display, endCol_display]
location_rec = location
# Create SLM pattern for the ring
myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                  wavelength, maskradius, res)
gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)
# Calculate phase pattern with WGS
targetAmp = myRecTarget
# rotate the target amp for non-zero rotation angle
if rot:
    targetAmp_rot, location_rot = myIMG.rotate_targetAmp(targetAmp, rotAngle, location, Plot=True)
    location = location_rot
    targetAmp = targetAmp_rot
myIMG.plotFocalplane(targetAmp, location)
Loop = 10
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
         ConfigFile["lattice geometry"] = "rec"
         ConfigFile["lattice coordinate"] = target_arr
         ConfigFile["rec location"] = location_rec
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