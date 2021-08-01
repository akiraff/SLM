from PyQt5.QtWidgets import QWidget, QFileDialog
import numpy as np
from pathlib import Path
from PIL import Image

class LoadAndSave(QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'SLM files:'
        self.left = 20
        self.top = 20
        self.width = 640
        self.height = 480

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)


    def LoadFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load WGS screen phase files:", "", "SLM screen Files (*.csv)")
        if file_name:
            print(file_name)
        file_path = Path(file_name)
        with open(file_path, newline='') as csvfile:
            SLM_screen_WGS = np.genfromtxt(csvfile, delimiter=',')
        return SLM_screen_WGS

    def LoadPhaseFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load WGS phase files:", "", "SLM WGS phase Files (*.csv)")
        if file_name:
            print(file_name)
        file_path = Path(file_name)
        with open(file_path, newline='') as csvfile:
            SLM_Phase = np.genfromtxt(csvfile, delimiter=',')
        return SLM_Phase

    def LoadIntensityFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load measured focal intensity files:", "", "Focal intensity Files (*.csv)")
        if file_name:
            print(file_name)
        file_path = Path(file_name)
        with open(file_path, newline='') as csvfile:
            intenArray = np.genfromtxt(csvfile, delimiter=',')
        return intenArray

    def LoadTargetAmpFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load target Amp files:", "",
                                                   "Target Amp file from last iteration (*.csv)")
        if file_name:
            print(file_name)
        file_path = Path(file_name)
        with open(file_path, newline='') as csvfile:
            targetAmp_adapt_lastiter = np.genfromtxt(csvfile, delimiter=',')
        return targetAmp_adapt_lastiter

    def LoadConfigFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load SLM config files:", "", "SLM Config Files (*.npy)")
        SLMconfig = np.load(file_name, allow_pickle=True).item()
        return SLMconfig

    def LoadCorrFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load SLM correction files:", "", "SLM correction Files (*.bmp)")
        image = Image.open(file_name)
        corr_data = np.asarray(image)
        return corr_data


    def SaveFileDialog(self, SLM_screen_WGS):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save WGS screen phase file:", "", "SLM screen Files (*.csv)")
        np.savetxt(file_name, SLM_screen_WGS, delimiter=",")

    def SavePhaseFileDialog(self, SLM_phase):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save WGS phase file:", "", "SLM phase WGS Files (*.csv)")
        np.savetxt(file_name, SLM_phase, delimiter=",")

    def SaveConfigFileDialog(self, SLMconfigFile):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save config file:", "", "SLM Config Files (*.npy)")
        np.save(file_name, SLMconfigFile)

    def SaveTargetAmpFileDialog(self, targetAmp):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save target Amp file:", "", "SLM target Amp Files (*.csv)")
        np.savetxt(file_name, targetAmp, delimiter=",")

    def SaveAberrCorrFileDialog(self, aberrConfig):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save config file:", "", "Aberr correction Config Files (*.npy)")
        np.save(file_name, aberrConfig)

# check load config file
#SLMconfig = np.load('SLM_config_rec.npy', allow_pickle=True).item()
#print(SLMconfig)