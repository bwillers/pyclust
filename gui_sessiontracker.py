# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'pyclust_session_tracker.ui'
#
# Created: Mon May 28 16:45:30 2012
#      by: PyQt4 UI code generator 4.9.1
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    _fromUtf8 = lambda s: s

class Ui_SessionTrackerDialog(object):
    def setupUi(self, SessionTrackerDialog):
        SessionTrackerDialog.setObjectName(_fromUtf8("SessionTrackerDialog"))
        SessionTrackerDialog.setWindowModality(QtCore.Qt.ApplicationModal)
        SessionTrackerDialog.resize(1057, 754)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(SessionTrackerDialog.sizePolicy().hasHeightForWidth())
        SessionTrackerDialog.setSizePolicy(sizePolicy)
        self.horizontalLayout = QtGui.QHBoxLayout(SessionTrackerDialog)
        self.horizontalLayout.setObjectName(_fromUtf8("horizontalLayout"))
        self.mplwidget = MatplotlibWidget(SessionTrackerDialog)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mplwidget.sizePolicy().hasHeightForWidth())
        self.mplwidget.setSizePolicy(sizePolicy)
        self.mplwidget.setAutoFillBackground(True)
        self.mplwidget.setObjectName(_fromUtf8("mplwidget"))
        self.horizontalLayout.addWidget(self.mplwidget)
        self.frame = QtGui.QFrame(SessionTrackerDialog)
        self.frame.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtGui.QFrame.Raised)
        self.frame.setObjectName(_fromUtf8("frame"))
        self.verticalLayout = QtGui.QVBoxLayout(self.frame)
        self.verticalLayout.setObjectName(_fromUtf8("verticalLayout"))
        self.radioButton_nomatch = QtGui.QRadioButton(self.frame)
        self.radioButton_nomatch.setObjectName(_fromUtf8("radioButton_nomatch"))
        self.verticalLayout.addWidget(self.radioButton_nomatch)
        spacerItem = QtGui.QSpacerItem(20, 661, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.buttonBox = QtGui.QDialogButtonBox(self.frame)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Maximum, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.buttonBox.sizePolicy().hasHeightForWidth())
        self.buttonBox.setSizePolicy(sizePolicy)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName(_fromUtf8("buttonBox"))
        self.verticalLayout.addWidget(self.buttonBox)
        self.horizontalLayout.addWidget(self.frame)

        self.retranslateUi(SessionTrackerDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("accepted()")), SessionTrackerDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL(_fromUtf8("rejected()")), SessionTrackerDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(SessionTrackerDialog)

    def retranslateUi(self, SessionTrackerDialog):
        SessionTrackerDialog.setWindowTitle(QtGui.QApplication.translate("SessionTrackerDialog", "Dialog", None, QtGui.QApplication.UnicodeUTF8))
        self.radioButton_nomatch.setText(QtGui.QApplication.translate("SessionTrackerDialog", "No Match", None, QtGui.QApplication.UnicodeUTF8))

from matplotlibwidget import MatplotlibWidget
