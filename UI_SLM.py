# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'UI_SLM.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtGui import QDoubleValidator, QValidator
import IMGpy
import slmpy
import time
from numpy import genfromtxt
import os.path
import numpy as np


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(743, 431)
       # MainWindow.setWindowTitle("SLM")
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
        self.Button_calculate.setGeometry(QtCore.QRect(140, 340, 93, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_calculate.setFont(font)
        self.Button_calculate.setObjectName("Button_calculate")
        self.Button_send = QtWidgets.QPushButton(self.centralwidget)
        self.Button_send.setGeometry(QtCore.QRect(450, 340, 93, 28))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setBold(True)
        font.setWeight(75)
        self.Button_send.setFont(font)
        self.Button_send.setObjectName("Button_send")
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
        self.label_12.setGeometry(QtCore.QRect(120, 260, 55, 16))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.arraySizex = QtWidgets.QLineEdit(self.centralwidget)
        self.arraySizex.setGeometry(QtCore.QRect(90, 280, 81, 22))
        self.arraySizex.setText("")
        self.arraySizex.setObjectName("arraySizex")
        self.label_13 = QtWidgets.QLabel(self.centralwidget)
        self.label_13.setGeometry(QtCore.QRect(220, 260, 55, 21))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(10)
        self.label_13.setFont(font)
        self.label_13.setObjectName("label_13")
        self.arraySizey = QtWidgets.QLineEdit(self.centralwidget)
        self.arraySizey.setGeometry(QtCore.QRect(200, 280, 81, 22))
        self.arraySizey.setText("")
        self.arraySizey.setObjectName("arraySizey")
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
        self.Button_send.clicked.connect(self.send)
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


    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.label.setText(_translate("MainWindow", "SLM"))
       # self.pixelpitch.setText(_translate("MainWindow", "12.5"))
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
        self.Button_calculate.setStatusTip(_translate("MainWindow", "Calculate the SLM phase and save it in a csv file."))
        self.Button_calculate.setText(_translate("MainWindow", "Calculate"))
        self.Button_send.setStatusTip(_translate("MainWindow", "Load the csv file and send it to SLM."))
        self.Button_send.setText(_translate("MainWindow", "Send"))
        self.label_9.setText(_translate("MainWindow", "Array spacing (micron)"))
        self.label_10.setText(_translate("MainWindow", "Tweezer "))
        self.label_11.setText(_translate("MainWindow", "Array size"))
        self.label_12.setText(_translate("MainWindow", "x"))
        self.label_13.setText(_translate("MainWindow", "y"))
        self.label_14.setText(_translate("MainWindow", "Aperture size (mm)"))
       # self.maskdiameter.setText(_translate("MainWindow", "10"))
        self.label_15.setText(_translate("MainWindow", "Threshold"))
       # self.threshold.setText(_translate("MainWindow", "0.04"))
        self.label_16.setText(_translate("MainWindow", "Loop"))
       # self.Loop.setText(_translate("MainWindow", "25"))
        self.mask.setText(_translate("MainWindow", "Use Aperture?"))
        self.save.setText(_translate("MainWindow", "Save?"))
        self.label_17.setText(_translate("MainWindow", "Distance from origin (micron)"))
        self.label_18.setText(_translate("MainWindow", "x"))
        self.label_19.setText(_translate("MainWindow", "y"))
        self.label_20.setText(_translate("MainWindow", "Beam waist (mm)"))

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
            slm = slmpy.SLMdisplay(isImageLock=False)
            resX, resY = slm.getSize()
            res = np.min([resX, resY])
            slm.close()
            print("resX =: ", resX)
            print("resY=: ", resY)
            myIMG = IMGpy.IMG(pixelpitch, [arraysizeBit, arraysizeBit], beamwaist, focallength, magnification,
                              wavelength, maskradius, res)
            gaussianAmp, gaussianPhase = myIMG.initSLMImage(mask=mask, Plot=False)
            Focalpitchx, Focalpitchy, targetAmp, location = myIMG.initFocalImage_RecLattice(distance, [spacingx, spacingy],
                                                                                            [arraysizex, arraysizey], Plot=False)
            print("Focal pitch size: ", Focalpitchx)
            WGScal = IMGpy.WGS(gaussianAmp, gaussianPhase, targetAmp)
            SLM_Amp, SLM_Phase, Focal_Amp, non_uniform = WGScal.fftLoop(Loop, threshold)
            myIMG.plotSLMplane(SLM_Amp)
            myIMG.plotSLMplane(SLM_Phase)
            myIMG.plotFocalplane(Focal_Amp, location)
            SLM_bit, fftSLM_IMG_Amp_norm = WGScal.SLM_IMG(SLM_Phase, resX, resY, Plot=True, Save=save)
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


    def send(self):
        if os.path.isfile('SLM.csv'):
           slm = slmpy.SLMdisplay(isImageLock=True)
           SLM_data = genfromtxt('SLM.csv', delimiter=',')
           SLM_IMG = SLM_data.astype('uint8')
           slm.updateArray(SLM_IMG)
           time.sleep(5)
           slm.close()
        else:
            print("File not exist! You need to generate a SLM.csv file first!")

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

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    SLM = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(SLM)
    SLM.show()
    sys.exit(app.exec_())
