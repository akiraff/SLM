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
        file_name, _ = QFileDialog.getOpenFileName(self, "Load WGS phase files:", "", "SLM Files (*.csv)")
        if file_name:
            print(file_name)
        file_path = Path(file_name)
        with open(file_path, newline='') as csvfile:
            SLM_screen_WGS = np.genfromtxt(csvfile, delimiter=',')
        return SLM_screen_WGS

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
        file_name, _ = QFileDialog.getSaveFileName(self, "Save WGS phase file:", "", "SLM Files (*.csv)")
        np.savetxt(file_name, SLM_screen_WGS, delimiter=",")

    def SaveConfigFileDialog(self, SLMconfigFile):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save config file:", "", "SLM Config Files (*.npy)")
        np.save(file_name, SLMconfigFile)


# check load config file
#SLMconfig = np.load('SLM_rec_config.npy', allow_pickle=True).item()
#print(SLMconfig['SLM resX'])
#LS = LoadAndSave()
#SLMconfig = LS.LoadConfigFileDialog()
#print(SLMconfig)