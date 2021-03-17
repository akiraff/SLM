from PyQt5.QtWidgets import QWidget, QFileDialog
import numpy as np
from pathlib import Path


class LoadAndSave(QWidget):

    def __init__(self):
        super().__init__()
        self.title = 'SLM files:'
        self.left = 10
        self.top = 10
        self.width = 640
        self.height = 480

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)


    def LoadFileDialog(self):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getOpenFileName(self, "Load WGS files:")
        if file_name:
            print(file_name)
        file_path = Path(file_name)
        with open(file_path, newline='') as csvfile:
            SLM_screen_WGS = np.genfromtxt(csvfile, delimiter=',')
        return SLM_screen_WGS

    def SaveFileDialog(self, SLM_screen_WGS):
        self.initUI()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        file_name, _ = QFileDialog.getSaveFileName(self, "Save file:", "", "SLM Files (*.csv)")
        np.savetxt(file_name, SLM_screen_WGS, delimiter=",")
