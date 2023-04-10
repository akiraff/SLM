import IMGpy
import numpy as np

# This script is used to calculate the intensity array used to generate phase pattern for adapted WGS
# The main usage of this script is for spacing + intensity feedback

# Load SLM parameter
from file_dialog_widget import LoadAndSave
from PyQt5.QtWidgets import QApplication

if __name__ == '__main__':

     import sys
     app = QApplication(sys.argv)
     bfirst_round = 0
     LS = LoadAndSave()
     if bfirst_round:
         print("Load measured intensity array for current round:")
         intenArray_cur = LS.LoadIntensityFileDialog()
         print("intenArray current:", intenArray_cur)
         sz_intenArray = np.size(intenArray_cur)
         print(sz_intenArray)
         IntenArray_adapt = np.ones(sz_intenArray)/sz_intenArray
         print(IntenArray_adapt)
     else:
         # intenArray_past is the intensity array used to create targetAmp in the last feedback iteraction
         print("Load adapted intenArray of previous round:")
         intenArray_past = LS.LoadIntensityFileDialog()
         print(intenArray_past)
         # intenArray_measured is the measured intensity for intensity feedback
         print("Load measured intensity Array of previous round:")
         intenArray_measured = LS.LoadIntensityFileDialog()
         print(intenArray_measured)
         IntenArray_measured_avg = np.mean(intenArray_measured)
         IntenArray_adapt = np.multiply(np.divide(IntenArray_measured_avg, intenArray_measured), intenArray_past)
         print(IntenArray_adapt)
     # IntenArray_adapt is the intensity value used to create targetAmp to generate adapated WGS of current round
     print("Now we save adapted intenArray of this round: ")
     LS.SaveIntensityFileDialog(IntenArray_adapt)