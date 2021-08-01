from PyQt5.QtWidgets import QApplication, QDialog, QLineEdit, QDialogButtonBox, QGridLayout, QLabel
from PyQt5.QtGui import QDoubleValidator, QValidator

class AberrInputDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Aberration corection")

        layout = QGridLayout(self)
        Z20Label = QLabel('Z(2,0):')
        layout.addWidget(Z20Label, 0, 0)
        self.Z20Input = QLineEdit()
        layout.addWidget(self.Z20Input, 1, 0)
        Z2_2Label = QLabel('Z(2,-2):')
        layout.addWidget(Z2_2Label, 0, 1)
        self.Z2_2Input = QLineEdit()
        layout.addWidget(self.Z2_2Input, 1, 1)
        Z22Label = QLabel('Z(2,2):')
        layout.addWidget(Z22Label, 0, 2)
        self.Z22Input = QLineEdit()
        layout.addWidget(self.Z22Input, 1, 2)
        Z3_1Label = QLabel('Z(3,-1):')
        layout.addWidget(Z3_1Label, 0, 3)
        self.Z3_1Input = QLineEdit()
        layout.addWidget(self.Z3_1Input, 1, 3)

        Z31Label = QLabel('Z(3,1):')
        layout.addWidget(Z31Label, 2, 0)
        self.Z31Input = QLineEdit()
        layout.addWidget(self.Z31Input, 3, 0)
        Z3_3Label = QLabel('Z(3,-3):')
        layout.addWidget(Z3_3Label, 2, 1)
        self.Z3_3Input = QLineEdit()
        layout.addWidget(self.Z3_3Input, 3, 1)
        Z33Label = QLabel('Z(3,3):')
        layout.addWidget(Z33Label, 2, 2)
        self.Z33Input = QLineEdit()
        layout.addWidget(self.Z33Input, 3, 2)

        Z40Label = QLabel('Z(4,0):')
        layout.addWidget(Z40Label, 2, 3)
        self.Z40Input = QLineEdit()
        layout.addWidget(self.Z40Input, 3, 3)

        Z42Label = QLabel('Z(4,2):')
        layout.addWidget(Z42Label, 4, 0)
        self.Z42Input = QLineEdit()
        layout.addWidget(self.Z42Input, 5, 0)
        Z4_2Label = QLabel('Z(4,-2):')
        layout.addWidget(Z4_2Label, 4, 1)
        self.Z4_2Input = QLineEdit()
        layout.addWidget(self.Z4_2Input, 5, 1)
        Z44Label = QLabel('Z(4,4):')
        layout.addWidget(Z44Label, 4, 2)
        self.Z44Input = QLineEdit()
        layout.addWidget(self.Z44Input, 5, 2)
        Z4_4Label = QLabel('Z(4,-4):')
        layout.addWidget(Z4_4Label, 4, 3)
        self.Z4_4Input = QLineEdit()
        layout.addWidget(self.Z4_4Input, 5, 3)

        self.Z20Input.editingFinished.connect(self.validating_floatz20percent)
        self.Z2_2Input.editingFinished.connect(self.validating_floatz2_2percent)
        self.Z22Input.editingFinished.connect(self.validating_floatz22percent)

        self.Z3_1Input.editingFinished.connect(self.validating_floatz3_1percent)
        self.Z31Input.editingFinished.connect(self.validating_floatz31percent)
        self.Z3_3Input.editingFinished.connect(self.validating_floatz3_3percent)
        self.Z33Input.editingFinished.connect(self.validating_floatz33percent)

        self.Z40Input.editingFinished.connect(self.validating_floatz40percent)
        self.Z42Input.editingFinished.connect(self.validating_floatz42percent)
        self.Z4_2Input.editingFinished.connect(self.validating_floatz4_2percent)
        self.Z44Input.editingFinished.connect(self.validating_floatz44percent)
        self.Z4_4Input.editingFinished.connect(self.validating_floatz4_4percent)

        buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel, self)
        layout.addWidget(buttonBox)
        buttonBox.accepted.connect(self.accept)
        buttonBox.rejected.connect(self.reject)

    def getInputs(self):
        # convert no input to zero
        aberrConfig = {}
        aberrConfigzind = {}
        if self.Z20Input.text() == "":
            aberrConfig["Z(2,0)"] = 0
            aberrConfigzind[3] = 0
        else:
            aberrConfig["Z(2,0)"] = float(self.Z20Input.text())
            aberrConfigzind[3] = float(self.Z20Input.text())

        if self.Z2_2Input.text() == "":
            aberrConfig["Z(2,-2)"] = 0
            aberrConfigzind[4] = 0
        else:
            aberrConfig["Z(2,-2)"] = float(self.Z2_2Input.text())
            aberrConfigzind[4] = float(self.Z2_2Input.text())

        if self.Z22Input.text() == "":
            aberrConfig["Z(2,2)"] = 0
            aberrConfigzind[5] = 0
        else:
            aberrConfig["Z(2,2)"] = float(self.Z22Input.text())
            aberrConfigzind[5] = float(self.Z22Input.text())

        if self.Z3_1Input.text() == "":
            aberrConfig["Z(3,-1)"] = 0
            aberrConfigzind[6] = 0
        else:
            aberrConfig["Z(3,-1)"] = float(self.Z3_1Input.text())
            aberrConfigzind[6] = float(self.Z3_1Input.text())

        if self.Z31Input.text() == "":
            aberrConfig["Z(3,1)"] = 0
            aberrConfigzind[7] = 0
        else:
            aberrConfig["Z(3,1)"] = float(self.Z31Input.text())
            aberrConfigzind[7] = float(self.Z31Input.text())

        if self.Z3_3Input.text() == "":
            aberrConfig["Z(3,-3)"] = 0
            aberrConfigzind[8] = 0
        else:
            aberrConfig["Z(3,-3)"] = float(self.Z3_3Input.text())
            aberrConfigzind[8] = float(self.Z3_3Input.text())

        if self.Z33Input.text() == "":
            aberrConfig["Z(3,3)"] = 0
            aberrConfigzind[9] = 0
        else:
            aberrConfig["Z(3,3)"] = float(self.Z33Input.text())
            aberrConfigzind[9] = float(self.Z33Input.text())

        if self.Z40Input.text() == "":
            aberrConfig["Z(4,0)"] = 0
            aberrConfigzind[10] = 0
        else:
            aberrConfig["Z(4,0)"] = float(self.Z40Input.text())
            aberrConfigzind[10] = float(self.Z40Input.text())

        if self.Z42Input.text() == "":
            aberrConfig["Z(4,2)"] = 0
            aberrConfigzind[11] = 0
        else:
            aberrConfig["Z(4,2)"] = float(self.Z42Input.text())
            aberrConfigzind[11] = float(self.Z42Input.text())

        if self.Z4_2Input.text() == "":
            aberrConfig["Z(4,-2)"] = 0
            aberrConfigzind[12] = 0
        else:
            aberrConfig["Z(4,-2)"] = float(self.Z4_2Input.text())
            aberrConfigzind[12] = float(self.Z4_2Input.text())

        if self.Z44Input.text() == "":
            aberrConfig["Z(4,4)"] = 0
            aberrConfigzind[13] = 0
        else:
            aberrConfig["Z(4,4)"] = float(self.Z44Input.text())
            aberrConfigzind[13] = float(self.Z44Input.text())

        if self.Z4_4Input.text() == "":
            aberrConfig["Z(4,-4)"] = 0
            aberrConfigzind[14] = 0
        else:
            aberrConfig["Z(4,-4)"] = float(self.Z4_4Input.text())
            aberrConfigzind[14] = float(self.Z4_4Input.text())
        return aberrConfig, aberrConfigzind

    def validating_floatz20percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z20Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z20Input.setFocus()
        else:
            self.Z20Input.setText("")
            print('Value for Z(2,0) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz2_2percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z2_2Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z2_2Input.setFocus()
        else:
            self.Z2_2Input.setText("")
            print('Value for Z(2,-2) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz22percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z22Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z22Input.setFocus()
        else:
            self.Z22Input.setText("")
            print('Value for Z(2,2) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz3_1percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z3_1Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z3_1Input.setFocus()
        else:
            self.Z3_1Input.setText("")
            print('Value for Z(3,-1) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz31percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z31Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z31Input.setFocus()
        else:
            self.Z31Input.setText("")
            print('Value for Z(3,1) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz3_3percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z3_3Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z3_3Input.setFocus()
        else:
            self.Z3_3Input.setText("")
            print('Value for Z(3,-3) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz33percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z33Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z33Input.setFocus()
        else:
            self.Z33Input.setText("")
            print('Value for Z(3,3) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz40percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z40Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z40Input.setFocus()
        else:
            self.Z40Input.setText("")
            print('Value for Z(4,0) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz42percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z42Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z42Input.setFocus()
        else:
            self.Z42Input.setText("")
            print('Value for Z(4,2) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz4_2percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z4_2Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z4_2Input.setFocus()
        else:
            self.Z4_2Input.setText("")
            print('Value for Z(4,-2) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz44percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z44Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z44Input.setFocus()
        else:
            self.Z44Input.setText("")
            print('Value for Z(4,4) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

    def validating_floatz4_4percent(self):
        validating_rule = QDoubleValidator(-1, 1, 5)
        #   print(validating_rule.validate(self.Loop.text(), 14))
        if validating_rule.validate(self.Z4_4Input.text(), 14)[0] == QValidator.Acceptable:
            self.Z4_4Input.setFocus()
        else:
            self.Z4_4Input.setText("")
            print('Value for Z(4,-4) needs to be in (-1,1)! This input will be set to 0 if you do not'
                  ' change it.')

