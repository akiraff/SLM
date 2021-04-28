# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'UI_SLM.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Initial file generated by PyQt5 UI code, modified by Fang
# Do not re-run pyuic5


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QDoubleValidator, QValidator
import IMGpy
import slmpy
from Aberration import Zernike
from file_dialog_widget import LoadAndSave
import sys
import numpy as np
import time


class Ui_MainWindow(object):
    # UI mode: displayMode = 1, connect to SLM and ready to send data to SLM
    #           displayMode = 0, not connect to SLM, calculate Phase pattern through phase-fixed WGS
    def setupUi(self, MainWindow):
        # Please set displayMode here:
        self.displayMode = 1

        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(743, 531)
        if self.displayMode:
           self.slm = slmpy.SLMdisplay(isImageLock=True)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(0, 30, 55, 16))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.pixelpitch = QtWidgets.QLineEdit(self.centralwidget)
        self.pixelpitch.setGeometry(QtCore.QRect(150, 30, 113, 22))
        self.pixelpitch.setObjectName("pixelpitch")
        self.label_2 = QtWidgets.QLabel(self.centralwidget)
        self.label_2.setGeometry(QtCore.QRect(130, 0, 151, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.label_3 = QtWidgets.QLabel(self.centralwidget)
        self.label_3.setGeometry(QtCore.QRect(330, 30, 55, 16))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_3.setFont(font)
        self.label_3.setObjectName("label_3")
        self.label_4 = QtWidgets.QLabel(self.centralwidget)
        self.label_4.setGeometry(QtCore.QRect(400, 0, 111, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_4.setFont(font)
        self.label_4.setObjectName("label_4")
        self.label_5 = QtWidgets.QLabel(self.centralwidget)
        self.label_5.setGeometry(QtCore.QRect(330, 200, 81, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.label_6 = QtWidgets.QLabel(self.centralwidget)
        self.label_6.setGeometry(QtCore.QRect(430, 110, 141, 20))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.focallength = QtWidgets.QLineEdit(self.centralwidget)
        self.focallength.setGeometry(QtCore.QRect(440, 140, 113, 22))
        self.focallength.setObjectName("focallength")
        self.label_7 = QtWidgets.QLabel(self.centralwidget)
        self.label_7.setGeometry(QtCore.QRect(610, 100, 121, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.magnification = QtWidgets.QLineEdit(self.centralwidget)
        self.magnification.setGeometry(QtCore.QRect(600, 140, 113, 22))
        self.magnification.setObjectName("magnification")
        self.label_8 = QtWidgets.QLabel(self.centralwidget)
        self.label_8.setGeometry(QtCore.QRect(430, 190, 131, 20))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.arraySizeBit = QtWidgets.QComboBox(self.centralwidget)
        self.arraySizeBit.setGeometry(QtCore.QRect(410, 30, 73, 22))
        self.arraySizeBit.setObjectName("arraySizeBit")
        self.arraySizeBit.addItem("")
        self.arraySizeBit.addItem("")
        self.arraySizeBit.addItem("")
        self.arraySizeBit.addItem("")
        self.wavelength = QtWidgets.QLineEdit(self.centralwidget)
        self.wavelength.setGeometry(QtCore.QRect(440, 230, 113, 22))
        self.wavelength.setObjectName("wavelength")
        self.Button_calculate = QtWidgets.QPushButton(self.centralwidget)
        self.Button_calculate.setGeometry(QtCore.QRect(290, 340, 93, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_calculate.setFont(font)
        self.Button_calculate.setObjectName("Button_calculate")
        """""
        self.Button_init = QtWidgets.QPushButton(self.centralwidget)
        self.Button_init.setGeometry(QtCore.QRect(450, 370, 153, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_init.setFont(font)
        self.Button_init.setObjectName("Button_init")
        """""
        self.Button_send = QtWidgets.QPushButton(self.centralwidget)
        self.Button_send.setGeometry(QtCore.QRect(450, 410, 153, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_send.setFont(font)
        self.Button_send.setObjectName("Button_send")

        self.Button_stop = QtWidgets.QPushButton(self.centralwidget)
        self.Button_stop.setGeometry(QtCore.QRect(450, 450, 153, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_stop.setFont(font)
        self.Button_stop.setObjectName("Button_stop")

        self.label_9 = QtWidgets.QLabel(self.centralwidget)
        self.label_9.setGeometry(QtCore.QRect(110, 90, 181, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.label_10 = QtWidgets.QLabel(self.centralwidget)
        self.label_10.setGeometry(QtCore.QRect(0, 200, 91, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.spacingx = QtWidgets.QLineEdit(self.centralwidget)
        self.spacingx.setGeometry(QtCore.QRect(110, 130, 51, 22))
        self.spacingx.setText("")
        self.spacingx.setObjectName("spacingx")
        self.label_11 = QtWidgets.QLabel(self.centralwidget)
        self.label_11.setGeometry(QtCore.QRect(150, 230, 121, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_11.setFont(font)
        self.label_11.setObjectName("label_11")
        self.label_12 = QtWidgets.QLabel(self.centralwidget)
        self.label_12.setGeometry(QtCore.QRect(150, 260, 55, 16))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.arraySizex = QtWidgets.QLineEdit(self.centralwidget)
        self.arraySizex.setGeometry(QtCore.QRect(120, 280, 81, 22))
        self.arraySizex.setText("")
        self.arraySizex.setObjectName("arraySizex")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(260, 260, 55, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.arraySizey = QtWidgets.QLineEdit(self.centralwidget)
        self.arraySizey.setGeometry(QtCore.QRect(230, 280, 81, 22))
        self.arraySizey.setText("")
        self.arraySizey.setObjectName("arraySizey")

        self.label_21 = QtWidgets.QLabel(self.centralwidget)
        self.label_21.setGeometry(QtCore.QRect(15, 260, 100, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_21.setFont(font)
        self.label_21.setObjectName("label_21")
        self.arrayGeometry = QtWidgets.QComboBox(self.centralwidget)
        self.arrayGeometry.setGeometry(QtCore.QRect(10, 280, 81, 22))
        self.arrayGeometry.setObjectName("arrayGeometry")
        self.arrayGeometry.addItem("")
        self.arrayGeometry.addItem("")
        self.arrayGeometry.addItem("")

        self.label_22 = QtWidgets.QLabel(self.centralwidget)
        self.label_22.setGeometry(QtCore.QRect(0, 420, 121, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_22.setFont(font)
        self.label_22.setObjectName("label_22")

        self.label_23 = QtWidgets.QLabel(self.centralwidget)
        self.label_23.setGeometry(QtCore.QRect(130, 400, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_23.setFont(font)
        self.label_23.setObjectName("label_23")
        self.zind = QtWidgets.QLineEdit(self.centralwidget)
        self.zind.setGeometry(QtCore.QRect(150, 430, 51, 21))
        self.zind.setObjectName("zind")

        self.label_24 = QtWidgets.QLabel(self.centralwidget)
        self.label_24.setGeometry(QtCore.QRect(250, 400, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_24.setFont(font)
        self.label_24.setObjectName("label_24")
        self.zindpercent = QtWidgets.QLineEdit(self.centralwidget)
        self.zindpercent.setGeometry(QtCore.QRect(250, 430, 51, 21))
        self.zindpercent.setObjectName("zindpercent")

        self.label_25 = QtWidgets.QLabel(self.centralwidget)
        self.label_25.setGeometry(QtCore.QRect(330, 420, 101, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(12)
        font.setBold(True)
        font.setItalic(True)
        font.setWeight(75)
        self.label_25.setFont(font)
        self.label_25.setObjectName("label_25")

        self.Button_load = QtWidgets.QPushButton(self.centralwidget)
        self.Button_load.setGeometry(QtCore.QRect(140, 470, 93, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_load.setFont(font)
        self.Button_load.setObjectName("Button_load")

        self.label_14 = QtWidgets.QLabel(self.centralwidget)
        self.label_14.setGeometry(QtCore.QRect(430, 270, 151, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.maskradius = QtWidgets.QLineEdit(self.centralwidget)
        self.maskradius.setGeometry(QtCore.QRect(440, 300, 113, 22))
        self.maskradius.setObjectName("maskradius")
        self.label_15 = QtWidgets.QLabel(self.centralwidget)
        self.label_15.setGeometry(QtCore.QRect(520, 10, 81, 16))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.threshold = QtWidgets.QLineEdit(self.centralwidget)
        self.threshold.setGeometry(QtCore.QRect(520, 30, 71, 22))
        self.threshold.setObjectName("threshold")
        self.label_16 = QtWidgets.QLabel(self.centralwidget)
        self.label_16.setGeometry(QtCore.QRect(640, 0, 55, 31))
        font = QtGui.QFont()
        font.setPointSize(10)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.Loop = QtWidgets.QLineEdit(self.centralwidget)
        self.Loop.setGeometry(QtCore.QRect(630, 30, 71, 22))
        self.Loop.setObjectName("Loop")
        self.mask = QtWidgets.QCheckBox(self.centralwidget)
        self.mask.setGeometry(QtCore.QRect(580, 290, 161, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.mask.setFont(font)
        self.mask.setObjectName("mask")
        self.save = QtWidgets.QCheckBox(self.centralwidget)
        self.save.setGeometry(QtCore.QRect(150, 300, 161, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.save.setFont(font)
        self.save.setObjectName("save")

        self.saveConfig = QtWidgets.QCheckBox(self.centralwidget)
        self.saveConfig.setGeometry(QtCore.QRect(150, 330, 161, 41))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.saveConfig.setFont(font)
        self.saveConfig.setObjectName("saveConfig")

        self.label_17 = QtWidgets.QLabel(self.centralwidget)
        self.label_17.setGeometry(QtCore.QRect(90, 170, 231, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_17.setFont(font)
        self.label_17.setObjectName("label_17")
        self.distance = QtWidgets.QLineEdit(self.centralwidget)
        self.distance.setGeometry(QtCore.QRect(150, 200, 71, 22))
        self.distance.setObjectName("distance")
        self.label_18 = QtWidgets.QLabel(self.centralwidget)
        self.label_18.setGeometry(QtCore.QRect(120, 110, 55, 16))
        self.label_18.setObjectName("label_18")
        self.label_19 = QtWidgets.QLabel(self.centralwidget)
        self.label_19.setGeometry(QtCore.QRect(210, 110, 55, 16))
        self.label_19.setObjectName("label_19")
        self.spacingy = QtWidgets.QLineEdit(self.centralwidget)
        self.spacingy.setGeometry(QtCore.QRect(200, 130, 51, 22))
        self.spacingy.setText("")
        self.spacingy.setObjectName("spacingy")
        self.label_20 = QtWidgets.QLabel(self.centralwidget)
        self.label_20.setGeometry(QtCore.QRect(590, 190, 161, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_20.setFont(font)
        self.label_20.setObjectName("label_20")
        self.beamwaist = QtWidgets.QLineEdit(self.centralwidget)
        self.beamwaist.setGeometry(QtCore.QRect(600, 230, 113, 22))
        self.beamwaist.setObjectName("beamwaist")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 743, 26))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

        self.Button_calculate.clicked.connect(self.calculate)
        if self.displayMode:
            self.Button_send.clicked.connect(self.send)
            self.Button_stop.clicked.connect(self.stop)
        self.Button_load.clicked.connect(self.load)
        self.Loop.editingFinished.connect(self.validating_integerLoop)
        self.arraySizex.editingFinished.connect(self.validating_integerSizex)
        self.arraySizey.editingFinished.connect(self.validating_integerSizey)
        self.pixelpitch.editingFinished.connect(self.validating_floatpixelpitch)
        self.spacingx.editingFinished.connect(self.validating_floatspacingx)
        self.spacingy.editingFinished.connect(self.validating_floatspacingy)
        self.distance.editingFinished.connect(self.validating_floatdistance)
        self.threshold.editingFinished.connect(self.validating_floatthreshold)
        self.focallength.editingFinished.connect(self.validating_floatfocallength)
        self.magnification.editingFinished.connect(self.validating_floatmagnification)
        self.wavelength.editingFinished.connect(self.validating_floatwavelength)
        self.beamwaist.editingFinished.connect(self.validating_floatbeamwaist)
        self.maskradius.editingFinished.connect(self.validating_floatmaskradius)
        self.zind.editingFinished.connect(self.validating_integerzind)
        self.zindpercent.editingFinished.connect(self.validating_floatzindpercent)



    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        if self.displayMode:
           MainWindow.setWindowTitle(_translate("MainWindow", "nilab_SLM (Display mode: SLM connected!)"))
        else:
           MainWindow.setWindowTitle(_translate("MainWindow", "nilab_SLM (Calculate mode: SLM disconnected!)"))
        self.label.setText(_translate("MainWindow", "SLM"))
        self.label_2.setText(_translate("MainWindow", " pixel pitch (micron)"))
        self.label_3.setText(_translate("MainWindow", "WGS"))
        self.label_4.setText(_translate("MainWindow", "FFT grid (bit)"))
        self.label_5.setText(_translate("MainWindow", "Imaging"))
        self.label_6.setText(_translate("MainWindow", "Focal length (mm)"))
        self.label_7.setText(_translate("MainWindow", "Magnification"))
        self.label_8.setText(_translate("MainWindow", "Wave length (nm)"))
        self.arraySizeBit.setItemText(0, _translate("MainWindow", "11"))
        self.arraySizeBit.setItemText(1, _translate("MainWindow", "12"))
        self.arraySizeBit.setItemText(2, _translate("MainWindow", "13"))
        self.arraySizeBit.setItemText(3, _translate("MainWindow", "14"))
        self.Button_calculate.setStatusTip(_translate("MainWindow", "Calculate the SLM phase and you can choose to save it in a csv file."))
        self.Button_calculate.setText(_translate("MainWindow", "Calculate"))
        self.Button_send.setStatusTip(_translate("MainWindow", "Load the csv file and send it to SLM."))
        self.Button_send.setText(_translate("MainWindow", "Load and Send"))
        self.Button_stop.setStatusTip(_translate("MainWindow", "Stop SLM."))
        self.Button_stop.setText(_translate("MainWindow", "Stop"))
        self.Button_load.setStatusTip(_translate("MainWindow", "Load the config and phase file, add aberration correction, then save."))
        self.Button_load.setText(_translate("MainWindow", "Load"))
        self.label_9.setText(_translate("MainWindow", "Array spacing (micron)"))
        self.label_10.setText(_translate("MainWindow", "Tweezer "))
        self.label_11.setText(_translate("MainWindow", "Array size"))
        self.label_12.setText(_translate("MainWindow", "x"))
        self.label_13.setText(_translate("MainWindow", "y"))
        self.label_14.setText(_translate("MainWindow", "Aperture size (mm)"))
        self.label_15.setText(_translate("MainWindow", "Threshold"))
        self.label_16.setText(_translate("MainWindow", "Loop"))
        self.mask.setText(_translate("MainWindow", "Use Aperture?"))
        self.save.setText(_translate("MainWindow", "Save?"))
        self.saveConfig.setText(_translate("MainWindow", "Save config?"))
        self.label_17.setText(_translate("MainWindow", "Distance from origin (micron)"))
        self.label_18.setText(_translate("MainWindow", "x"))
        self.label_19.setText(_translate("MainWindow", "y"))
        self.label_20.setText(_translate("MainWindow", "Beam waist (mm)"))
        self.label_21.setText(_translate("MainWindow", "Geometry"))
        self.arrayGeometry.setItemText(0, _translate("MainWindow", "Rec"))
        self.arrayGeometry.setItemText(1, _translate("MainWindow", "Triangle"))
        self.arrayGeometry.setItemText(2, _translate("MainWindow", "Kagome"))
        self.label_22.setText(_translate("MainWindow", "Aberration"))
        self.label_23.setText(_translate("MainWindow", "Zernike index"))
        self.label_24.setText(_translate("MainWindow", "Percent"))
        self.label_25.setText(_translate("MainWindow", "Display"))

    def calculate(self):
        if self.pixelpitch.text() == "":
            print("Please input pixelpitch value of the SLM!")
        elif self.spacingx.text() == "":
            print("Please input array spacing along x!")
        elif self.spacingy.text() == "":
            print("Please input array spacing aling y!")
        elif self.distance.text() == "":
            print("Please specify the distance of the array from the origin!")
        elif self.arraySizex.text() == "":
            print("Please specify the array size along x!")
        elif self.arraySizey.text() == "":
            print("Please specify the array size along y!")
        elif self.threshold.text() == "":
            print("Please specify the threshold for the algorithm!")
        elif self.Loop.text() == "":
            print("Please specify the loop number for the algorithm!")
        elif self.focallength.text()== "":
            print("Please specify the focal length of the objective!")
        elif self.magnification.text()== "":
            print("Please specify the magnification of the imaging system!")
        elif self.wavelength.text() == "":
            print("Please specify the wave length of the light!")
        elif self.beamwaist.text() == "":
            print("Please specify the beam waist of the gaussian beam!")
        elif self.maskradius.text() == "":
            print("Please specify the aperture of the imaging system!")
        else:
            pixelpitch = float(self.pixelpitch.text())*1e-6
            print("pixel pitch: ", pixelpitch)
            spacingx = float(self.spacingx.text())*1e-6
            print("spacingx: ", spacingx)
            spacingy = float(self.spacingy.text()) * 1e-6
            print("spacingy: ", spacingy)
            distance = float(self.distance.text()) * 1e-6
            print("distance: ", distance)
            arraysizex = int(self.arraySizex.text())
            print("array size along x: ", arraysizex)
            arraysizey = int(self.arraySizey.text())
            print("array size along y: ", arraysizey)
            lattice = self.arrayGeometry.currentText()
            print("array geometry: ", lattice)
            arraysizeBit = int(self.arraySizeBit.currentText())
            print("FFT grid (Bit): ", arraysizeBit)
            threshold = float(self.threshold.text())
            print("threshold: ", threshold)
            Loop = int(self.Loop.text())
            print("Loop: ", Loop)
            focallength = float(self.focallength.text())*1e-3
            print("Focal length: ", focallength)
            magnification = float(self.magnification.text())
            print("Magnification: ", magnification)
            wavelength = float(self.wavelength.text())*1e-9
            print("Wave length: ", wavelength)
            beamwaist = float(self.beamwaist.text())*1e-3
            print("Input beam waist: ", beamwaist)
            maskradius = float(self.maskradius.text())*1e-3
            print("Aperture size: ", maskradius)
            if self.mask.isChecked():
                mask = 1
            else:
                mask = 0
            print("Use aperture?", mask)
            if self.save.isChecked():
                save = 1
            else:
                save = 0
            print("Save?", save)
            if self.saveConfig.isChecked():
                saveConfig = 1
            else:
                saveConfig = 0
            print("Save Config?", saveConfig)
            #slm = slmpy.SLMdisplay(isImageLock=False)
           # slm = slmpy.SLMdisplay()
            # Directly specify the resolution corresponding to the Hamamastu SLM dimension
            resX, resY = 1272, 1024
            res = np.min([resX, resY])
           # slm.close()
            print("resX =: ", resX)
            print("resY=: ", resY)
            myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                              wavelength, maskradius, res)
            gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)
            if lattice == "Rec":
                Focalpitchx, Focalpitchy, targetAmp, location = myIMG.initFocalImage_RecLattice(distance, [spacingx, spacingy],
                                                                                            [arraysizex, arraysizey], Plot=False)
            elif lattice == "Kagome":
                Focalpitchx, Focalpitchy, targetAmp, location = myIMG.initFocalImage_KagomeLattice(distance,
                                                                                            [spacingx, spacingy],
                                                                                            [arraysizex, arraysizey],
                                                                                            Plot=True)
            else:
                Focalpitchx, Focalpitchy, targetAmp, location = myIMG.initFocalImage_KagomeLattice(distance,
                                                                                                   [spacingx, spacingy],
                                                                                                   [arraysizex,
                                                                                                    arraysizey], Triangle=True,
                                                                                                   Plot=True)
            print("Focal pitch size: ", Focalpitchx)
         #   targetAmp_diffrac = myIMG.modify_targetAmp(targetAmp, location)
            targetAmp_diffrac = targetAmp
            WGScal = IMGpy.WGS(gaussianAmp, gaussianPhase, targetAmp_diffrac)
            SLM_Amp, SLM_Phase, Focal_Amp, non_uniform = WGScal.fftLoop(Loop, threshold)
            myIMG.plotSLMplane(SLM_Amp)
            myIMG.plotSLMplane(SLM_Phase)
            myIMG.plotFocalplane(Focal_Amp, location)
            SLM_bit, fftSLM_IMG_Amp_norm, SLM_Screen_WGS = WGScal.SLM_IMG(SLM_Phase, resX, resY, Plot=True)
            # X is column
            colIMG = np.size(SLM_bit, axis=1)
            colFFT = np.size(SLM_Phase, axis=1)
            # Y is row
            rowIMG = np.size(SLM_bit, axis=0)
            rowFFT = np.size(SLM_Phase, axis=0)
            endRowIMG = rowIMG/2
            startRowIMG = round(endRowIMG-(location[1]-location[0])*rowIMG/rowFFT)
            endColIMG = colIMG/2
            startColIMG = round(endColIMG-(location[3]-location[2])*colIMG/colFFT)
            locationIMG = [int(startRowIMG), int(endRowIMG), int(startColIMG), int(endColIMG)]
            myIMG.plotFocalplane(fftSLM_IMG_Amp_norm, locationIMG)
            print("location = ", location)
            print("locationIMG = ", locationIMG)
            print("Phase image shape:", SLM_Screen_WGS.shape)
            if save:
                LS = LoadAndSave()
                LS.SaveFileDialog(SLM_Screen_WGS)
            if saveConfig:
                ConfigFile = {}
                ConfigFile["pixel pitch"] = pixelpitch
                ConfigFile["arr spacing x"] = spacingx
                ConfigFile["arr spacing y"] = spacingy
                ConfigFile["distance from origin"] = distance
                ConfigFile["lattice geometry"] = lattice
                ConfigFile["array size x"] = arraysizex
                ConfigFile["array size y"] = arraysizey
                ConfigFile["FFT grid size (bit)"] = arraysizeBit
                ConfigFile["WGS threshold"] = threshold
                ConfigFile["Loop"] = Loop
                ConfigFile["Focal length"] = focallength
                ConfigFile["Magnefication"] = magnification
                ConfigFile["wave length"] = wavelength
                ConfigFile["beam waist"] = beamwaist
                ConfigFile["aperture size"] = maskradius
                ConfigFile["use aperture?"] = mask
                ConfigFile["SLM resX"] = resX
                ConfigFile["SLM resY"] = resY
                LS = LoadAndSave()
                LS.SaveConfigFileDialog(ConfigFile)


    def send(self):
        LS = LoadAndSave()
        SLM_phase_data = LS.LoadFileDialog()
        SLM_bit = np.around((SLM_phase_data + np.pi) / (2 * np.pi) * 255)
        SLM_corrPattern = LS.LoadCorrFileDialog()
        SLM_corrected = SLM_bit + SLM_corrPattern
        SLM_wrappedPattern = np.mod(SLM_corrected, 256)
        # 1038 nm
        value_for2pi = 212
        SLM_displayPattern = np.around(SLM_wrappedPattern*value_for2pi/255).astype('uint8')
        #print(np.max(SLM_displayPattern))
        #print(np.min(SLM_displayPattern))
        self.slm.updateArray(SLM_displayPattern)

    def stop(self):
        self.slm.close()
        time.sleep(1)
        sys.exit()

    def load(self):
        if self.zind.text() == "":
            print("Please input zernike index input!")
        elif self.zindpercent.text() == "":
            print("Please input the percent !")
        ind_Zernike = int(self.zind.text())
        percent = float(self.zindpercent.text())
        LS = LoadAndSave()
        print("We first load SLM config:")
        SLMconfig = LS.LoadConfigFileDialog()
        print(SLMconfig)
        SLMResX = SLMconfig['SLM resX']
        SLMResY = SLMconfig['SLM resY']
        pixelpitch = SLMconfig['pixel pitch']
        aperture_radius = SLMconfig['aperture size']
        print("Next, we load the phase file:")
        SLM_screen_WGS = LS.LoadFileDialog()
        myOberrationCorr = Zernike(SLMResX, SLMResY, pixelpitch, aperture_radius, ind_Zernike, percent)
        SLM_aberr_screen = myOberrationCorr.phase_Zernike(Plot=True, Save=False)
        SLM_screen_WGS_aberr = SLM_screen_WGS + SLM_aberr_screen
        print("Finally, we add the oberration and save the phase file.")
        LS.SaveFileDialog(SLM_screen_WGS_aberr)

    def validating_integerLoop(self):
        validating_rule = QDoubleValidator(0, 100, 0)
     #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Loop.text(), 14)[0] == QValidator.Acceptable:
            self.Loop.setFocus()
        else:
            self.Loop.setText("")

    def validating_integerSizex(self):
        validating_rule = QDoubleValidator(0, 500, 0)
      #  print(validating_rule.validate(self.arraySizex.text(), 14))
        if validating_rule.validate(self.arraySizex.text(), 14)[0] == QValidator.Acceptable:
            self.arraySizex.setFocus()
        else:
            self.arraySizex.setText("")

    def validating_integerSizey(self):
        validating_rule = QDoubleValidator(0, 500, 0)
      #  print(validating_rule.validate(self.arraySizey.text(), 14))
        if validating_rule.validate(self.arraySizey.text(), 14)[0] == QValidator.Acceptable:
            self.arraySizey.setFocus()
        else:
            self.arraySizey.setText("")

    def validating_floatpixelpitch(self):
        validating_rule = QDoubleValidator(0, 50, 4)
      #  print(validating_rule.validate(self.pixelpitch.text(), 14))
        if validating_rule.validate(self.pixelpitch.text(), 14)[0] == QValidator.Acceptable:
            self.pixelpitch.setFocus()
        else:
            self.pixelpitch.setText("")

    def validating_floatspacingx(self):
        validating_rule = QDoubleValidator(0, 100, 4)
       # print(validating_rule.validate(self.spacingx.text(), 14))
        if validating_rule.validate(self.spacingx.text(), 14)[0] == QValidator.Acceptable:
            self.spacingx.setFocus()
        else:
            self.spacingx.setText("")

    def validating_floatspacingy(self):
        validating_rule = QDoubleValidator(0, 100, 4)
      #  print(validating_rule.validate(self.spacingy.text(), 14))
        if validating_rule.validate(self.spacingy.text(), 14)[0] == QValidator.Acceptable:
            self.spacingy.setFocus()
        else:
            self.spacingy.setText("")

    def validating_floatdistance(self):
        validating_rule = QDoubleValidator(0, 500, 4)
       # print(validating_rule.validate(self.distance.text(), 14))
        if validating_rule.validate(self.distance.text(), 14)[0] == QValidator.Acceptable:
            self.distance.setFocus()
        else:
            self.distance.setText("")

    def validating_floatthreshold(self):
        validating_rule = QDoubleValidator(0, 0.2, 5)
     #   print(validating_rule.validate(self.threshold.text(), 0))
        if validating_rule.validate(self.threshold.text(), 0)[0] == QValidator.Acceptable:
            self.threshold.setFocus()
        else:
            self.threshold.setText("")

    def validating_floatfocallength(self):
        validating_rule = QDoubleValidator(0, 1000, 4)
     #   print(validating_rule.validate(self.focallength.text(), 0))
        if validating_rule.validate(self.focallength.text(), 0)[0] == QValidator.Acceptable:
            self.focallength.setFocus()
        else:
            self.focallength.setText("")

    def validating_floatmagnification(self):
        validating_rule = QDoubleValidator(0, 10, 4)
     #   print(validating_rule.validate(self.magnification.text(), 0))
        if validating_rule.validate(self.magnification.text(), 0)[0] == QValidator.Acceptable:
            self.magnification.setFocus()
        else:
            self.magnification.setText("")

    def validating_floatwavelength(self):
        validating_rule = QDoubleValidator(0, 10000, 4)
      #  print(validating_rule.validate(self.wavelength.text(), 0))
        if validating_rule.validate(self.wavelength.text(), 0)[0] == QValidator.Acceptable:
            self.wavelength.setFocus()
        else:
            self.wavelength.setText("")

    def validating_floatbeamwaist(self):
        validating_rule = QDoubleValidator(0, 100, 4)
      #  print(validating_rule.validate(self.beamwaist.text(), 0))
        if validating_rule.validate(self.beamwaist.text(), 0)[0] == QValidator.Acceptable:
            self.beamwaist.setFocus()
        else:
            self.beamwaist.setText("")

    def validating_floatmaskradius(self):
        validating_rule = QDoubleValidator(0, 100, 4)
      #  print(validating_rule.validate(self.maskdiameter.text(), 0))
        if validating_rule.validate(self.maskradius.text(), 0)[0] == QValidator.Acceptable:
            self.maskradius.setFocus()
        else:
            self.maskradius.setText("")

    def validating_integerzind(self):
        validating_rule = QDoubleValidator(1, 27, 0)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.zind.text(), 14)[0] == QValidator.Acceptable:
            self.zind.setFocus()
        else:
            self.zind.setText("")

    def validating_floatzindpercent(self):
        validating_rule = QDoubleValidator(0, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.zindpercent.text(), 14)[0] == QValidator.Acceptable:
            self.zindpercent.setFocus()
        else:
            self.zindpercent.setText("")

if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    SLM = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(SLM)
    SLM.show()
    sys.exit(app.exec_())
