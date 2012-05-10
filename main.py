# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 23:43:03 2012

@author: Bernard
"""

#!/usr/bin/python -d


from PyQt4 import QtCore, QtGui

import numpy as np
import scipy.io as sio
import scipy.stats
import matplotlib as mpl
from matplotlibwidget import MatplotlibWidget

import pickle
import os
import sys
import time

from gui import Ui_MainWindow
import spikeset
#import features


class PyClustMainWindow(QtGui.QMainWindow):

    def switch_to_wavecutter(self):
        self.ui.stackedWidget.setCurrentIndex(1)
        self.updateWavecutterPlot()

    def switch_to_maindisplay(self):
        self.ui.stackedWidget.setCurrentIndex(0)
        self.updateFeaturePlot()

    def poopy(self, x=None):
        if x:
                print x, 'is poopy'
        else:
                print 'poopy'

    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.redrawing_proj = True

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # There is no spikeset or feature set on load
        self.spikeset = None
        self.current_feature = None
        self.prof_limits_reference = None
        self.clusters = []
        self.current_filename = None

        self.limit_mode = False
        self.redrawing_proj = False
        self.redrawing_details = False
        self.unsaved = False

        # Connect the handlers
        QtCore.QObject.connect(self.ui.pushButton,
            QtCore.SIGNAL("clicked()"),  self.button_load_click)

        QtCore.QObject.connect(self.ui.pushButton_copy_cluster,
            QtCore.SIGNAL("clicked()"), self.copyCluster)

        QtCore.QObject.connect(self.ui.pushButton_import,
            QtCore.SIGNAL("clicked()"), self.importBounds)

        QtCore.QObject.connect(self.ui.pushButton_save,
            QtCore.SIGNAL("clicked()"), self.save)

        QtCore.QObject.connect(self.ui.comboBox_feature_x_chan,
            QtCore.SIGNAL("currentIndexChanged(int)"),
            self.feature_channel_x_changed)

        QtCore.QObject.connect(self.ui.comboBox_feature_x,
            QtCore.SIGNAL("currentIndexChanged(int)"),
            self.feature_x_changed)

        QtCore.QObject.connect(self.ui.comboBox_feature_y_chan,
            QtCore.SIGNAL("currentIndexChanged(int)"),
            self.feature_channel_y_changed)

        QtCore.QObject.connect(self.ui.comboBox_feature_y,
            QtCore.SIGNAL("currentIndexChanged(int)"),
            self.feature_y_changed)

        QtCore.QObject.connect(self.ui.pushButton_next_projection,
            QtCore.SIGNAL("clicked()"), self.button_next_feature_click)

        QtCore.QObject.connect(self.ui.pushButton_previous_projection,
            QtCore.SIGNAL("clicked()"), self.button_prev_feature_click)

        QtCore.QObject.connect(self.ui.buttonGroup_scatter_plot_type,
            QtCore.SIGNAL("buttonClicked(int)"), self.updateFeaturePlot)

        QtCore.QObject.connect(self.ui.spinBox_markerSize,
            QtCore.SIGNAL("valueChanged(int)"), self.updateFeaturePlot)

        QtCore.QObject.connect(self.ui.pushButton_addLimit,
            QtCore.SIGNAL("clicked()"), self.button_add_limit_click)

        QtCore.QObject.connect(self.ui.pushButton_deleteLimit,
            QtCore.SIGNAL("clicked()"), self.button_delete_limit_click)

        QtCore.QObject.connect(self.ui.checkBox_refractory,
            QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)

        QtCore.QObject.connect(self.ui.pushButton_wavecutter_done,
            QtCore.SIGNAL("clicked()"), self.switch_to_maindisplay)

        QtCore.QObject.connect(self.ui.pushButton_wavecutter_start,
            QtCore.SIGNAL("clicked()"), self.switch_to_wavecutter)

        QtCore.QObject.connect(self.ui.pushButton_wavecutter_redraw,
            QtCore.SIGNAL("clicked()"), self.updateWavecutterPlot)

        QtCore.QObject.connect(self.ui.lineEdit_wavecutter_count,
            QtCore.SIGNAL("editingFinished()"), self.updateWavecutterPlot)

        QtCore.QObject.connect(self.ui.checkBox_wavecutter_refractory,
            QtCore.SIGNAL("stateChanged(int)"), self.updateWavecutterPlot)

        QtCore.QObject.connect(self.ui.spinBox_wavecutter_channel,
            QtCore.SIGNAL("valueChanged(int)"), self.updateWavecutterPlot)

        QtCore.QObject.connect(self.ui.pushButton_wavecutter_add_limit,
            QtCore.SIGNAL("clicked()"), self.wavecutter_add_limit)

        QtCore.QObject.connect(self.ui.pushButton_wavecutter_remove_limit,
            QtCore.SIGNAL("clicked()"), self.wavecutter_remove_limits)

        QtCore.QObject.connect(self.ui.pushButton_hide_all,
            QtCore.SIGNAL("clicked()"),
            lambda: self.hide_show_all_clusters(True))

        QtCore.QObject.connect(self.ui.pushButton_show_all,
            QtCore.SIGNAL("clicked()"),
            lambda: self.hide_show_all_clusters(False))

        QtCore.QObject.connect(self.ui.buttonGroup_trimmer,
            QtCore.SIGNAL("buttonClicked(int)"),
            lambda x: self.ui.stackedWidget_trimmer.setCurrentIndex(-x - 2))

        QtCore.QObject.connect(self.ui.stackedWidget_trimmer,
            QtCore.SIGNAL("currentChanged(int)"),
            lambda x: self.updateWavecutterPlot() if x == 0 \
                else self.updateOutlierPlot())

        # Set up the cluster list area
        self.labels_container = QtGui.QWidget()
        self.ui.scrollArea_cluster_list.setWidget(self.labels_container)

        layout = QtGui.QVBoxLayout(self.labels_container)
        self.buttonGroup_cluster = QtGui.QButtonGroup()
        self.buttonGroup_cluster.setExclusive(True)
        QtCore.QObject.connect(self.buttonGroup_cluster,
            QtCore.SIGNAL("buttonClicked(int)"), self.update_active_cluster)

        layout.addItem(QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum,
            QtGui.QSizePolicy.Expanding))
        QtCore.QObject.connect(self.ui.pushButton_add_cluster,
            QtCore.SIGNAL("clicked()"), self.button_add_cluster_click)
        QtCore.QObject.connect(self.ui.pushButton_delete_cluster,
            QtCore.SIGNAL("clicked()"), self.delete_cluster)

        # Add the unclustered entry
        layout = self.labels_container.layout()

        hlayout = QtGui.QHBoxLayout()

        self.ui.checkBox_show_unclustered = QtGui.QCheckBox(
                self.labels_container)
        self.ui.checkBox_show_unclustered.setChecked(True)
        QtCore.QObject.connect(self.ui.checkBox_show_unclustered,
            QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)
        hlayout.addWidget(self.ui.checkBox_show_unclustered)

        hlayout.addWidget(QtGui.QLabel('Unclustered',
            parent=self.labels_container))

        hlayout.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum))

        self.ui.checkBox_show_unclustered_exclusive = QtGui.QCheckBox(
                self.labels_container)
        self.ui.checkBox_show_unclustered_exclusive.setChecked(True)
        QtCore.QObject.connect(self.ui.checkBox_show_unclustered_exclusive,
            QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)
        hlayout.addWidget(self.ui.checkBox_show_unclustered_exclusive)

        layout.insertLayout(0, hlayout)

        # Create short hands for the plot widgets
        self.mp_wave = self.ui.mplwidget_waveform
        self.mp_isi = self.ui.mplwidget_isi
        self.mp_proj = self.ui.mplwidget_projection
        self.mp_wavecutter = self.ui.mplwidget_wavecutter
        self.mp_outlier = self.ui.mplwidget_outliers
        self.mp_drift = self.ui.mplwidget_drift

        pal = self.palette().window().color()
        bgcolor = (pal.red() / 255.0, pal.blue() / 255.0, pal.green() / 255.0)

        # Set the window background color on plots

        self.mp_wave.figure.clear()
        self.mp_wave.figure.set_facecolor(bgcolor)

        self.mp_isi.figure.clear()
        self.mp_isi.figure.set_facecolor(bgcolor)

        self.mp_proj.figure.clear()
        self.mp_proj.figure.set_facecolor(bgcolor)

        self.mp_wavecutter.figure.clear()
        self.mp_wavecutter.figure.set_facecolor(bgcolor)

        self.mp_outlier.figure.clear()
        self.mp_outlier.figure.set_facecolor(bgcolor)

        self.mp_drift.figure.clear()
        self.mp_drift.figure.set_facecolor(bgcolor)

        self.zoom_active = False
        self.mp_proj.mpl_connect('button_press_event', self.onMousePress)
        self.mp_proj.mpl_connect('button_release_event', self.onMouseRelease)
        self.mp_proj.mpl_connect('motion_notify_event', self.onMouseMove)

        self.wave_limit_mode = False
        self.mp_wavecutter.mpl_connect('button_press_event',
            self.wavecutter_onMousePress)
        self.mp_wavecutter.mpl_connect('button_release_event',
            self.wavecutter_onMouseRelease)
        self.mp_wavecutter.mpl_connect('motion_notify_event',
            self.wavecutter_onMouseMove)

        self._backup_mp_proj_paintEvent = self.mp_proj.paintEvent

        self.mp_proj.paintEvent = lambda e: self.mp_proj_onPaint(e)

        self.ui.comboBox_feature_x.clear()
        self.ui.comboBox_feature_y.clear()
        self.ui.comboBox_feature_x_chan.clear()
        self.ui.comboBox_feature_y_chan.clear()

        self.activeClusterRadioButton = self.buttonGroup_cluster.checkedButton

        self.redrawing_proj = False

        # Set up waveform plot plot axes
        self.mp_wave.figure.clear()
        self.mp_wave.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0, top=1, left=0.0, right=1)
        self.mp_wave.axes = [None, None, None, None]
        for i in range(0, 4):
            if i == 0:
                self.mp_wave.axes[i] = self.mp_wave.figure.add_subplot(2, 2,
                    i + 1)
            else:
                self.mp_wave.axes[i] = self.mp_wave.figure.add_subplot(2, 2,
                    i + 1, sharey=self.mp_wave.axes[0])
            self.mp_wave.axes[i].hold(False)
            self.mp_wave.axes[i].set_xticks([])
            self.mp_wave.axes[i].set_yticks([])

        # Set up ISI plot axes

        self.isi_bins = np.logspace(np.log10(0.1), np.log10(1e5), 100)
        self.isi_bin_centers = (self.isi_bins[0:-1] + self.isi_bins[1:]) / 2

        self.mp_isi.figure.clear()
        self.mp_isi.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0.15, top=1, left=0.0, right=1)
        self.mp_isi.axes = self.mp_isi.figure.add_subplot(1, 1, 1)
        self.mp_isi.axes.hold(False)
        self.mp_isi.axes.set_xscale('log')
        self.mp_isi.axes.set_xlim([0.1, 1e5])
        refractory_line = mpl.lines.Line2D([2, 2], self.mp_isi.axes.get_ylim(),
            color='r', linestyle='--')
        self.mp_isi.axes.add_line(refractory_line)
        burst_line = mpl.lines.Line2D([20, 20], self.mp_isi.axes.get_ylim(),
            color='b', linestyle='--')
        self.mp_isi.axes.add_line(burst_line)
        theta_line = mpl.lines.Line2D([125, 125], self.mp_isi.axes.get_ylim(),
            color='g', linestyle='--')
        self.mp_isi.axes.add_line(theta_line)
        self.mp_isi.axes.set_xticks([1e1, 1e2, 1e3, 1e4])
        self.mp_isi.draw()

        # Set up the drift axes
        self.mp_drift.figure.clear()
        self.mp_drift.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0, top=1, left=0.15, right=1)
        self.mp_drift.axes = self.mp_drift.figure.add_subplot(1, 1, 1)
        self.mp_drift.axes.set_xticks([])
        self.mp_drift.axes.set_yticks([])
        self.mp_drift.draw()

        # Set up the outlier axes
        self.mp_outlier.figure.clear()
        self.mp_outlier.axes = self.mp_outlier.figure.add_subplot(1, 1, 1)
        self.mp_outlier.draw()

        # Set up the feature axes
        self.mp_proj.figure.clear()
        self.mp_proj.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0.0, top=1, left=0.0, right=1)
        self.mp_proj.axes = self.mp_proj.figure.add_subplot(1, 1, 1)

        # Clear the stats labels

        self.ui.label_spike_count.setText('')
        self.ui.label_mean_rate.setText('')
        self.ui.label_burst.setText('')
        self.ui.label_csi.setText('')
        self.ui.label_refr_count.setText('')
        self.ui.label_refr_fp.setText('')
        self.ui.label_refr_frac.setText('')
        self.ui.label_isolation.setText('')

        # Set up the waveform axes
        self.mp_wavecutter.figure.clear()
        self.mp_wavecutter.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0, top=1, left=0.0, right=1)
        self.mp_wavecutter.axes = \
            self.mp_wavecutter.figure.add_subplot(1, 1, 1)
        self.mp_wavecutter.axes.hold(False)
        self.mp_wavecutter.axes.set_xticks([])
        self.mp_wavecutter.axes.set_yticks([])
        self.mp_wavecutter.draw()

        self.ui.stackedWidget.setCurrentIndex(0)

        layout = self.labels_container.layout()

        hlayout = QtGui.QHBoxLayout()

        self.checkBox_junk = QtGui.QCheckBox()
        self.checkBox_junk.setChecked(True)
        QtCore.QObject.connect(self.checkBox_junk,
            QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)

        self.radioButton_junk = QtGui.QRadioButton()
        self.radioButton_junk.setChecked(True)
        spacer = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum)

        hlayout.addWidget(self.checkBox_junk)
        hlayout.addItem(spacer)
        hlayout.addWidget(self.radioButton_junk)
        hlayout.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.buttonGroup_cluster.addButton(self.radioButton_junk)

        label = QtGui.QLabel()
        label.setText('Junk')

        hlayout.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum))

        hlayout.addWidget(label)
        layout.insertLayout(layout.count() - 1, hlayout)

    def hide_show_all_clusters(self, hidden):
        self.redrawing_proj = True
        for cluster in self.clusters:
            if cluster.radio == self.activeClusterRadioButton():
                continue
            cluster.check.setChecked(not hidden)
            if hidden:
                # dont show this on show all
                self.checkBox_junk.setChecked(False)
        self.ui.checkBox_show_unclustered.setChecked(not hidden)
        self.redrawing_proj = False
        self.updateFeaturePlot()

    def update_active_cluster(self):
        self.updateClusterDetailPlots()
        for cluster in self.clusters:
            cluster.check.setEnabled(True)
        if self.activeClusterRadioButton() == self.radioButton_junk:
            return
        check = self.activeClusterRadioButton().cluster_reference.check
        check.setEnabled(False)
        check.setChecked(True)

    def add_cluster(self, color=None):
        new_cluster = spikeset.Cluster(self.spikeset)
        if color:
            new_cluster.color = color

        layout = self.labels_container.layout()

        hlayout = QtGui.QHBoxLayout()

        check = QtGui.QCheckBox()
        check.setChecked(True)
        QtCore.QObject.connect(check, QtCore.SIGNAL("stateChanged(int)"),
            self.updateFeaturePlot)

        radio = QtGui.QRadioButton()
        radio.setChecked(True)
        spacer = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum)

        hlayout.addWidget(check)
        hlayout.addItem(spacer)
        hlayout.addWidget(radio)
        hlayout.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.buttonGroup_cluster.addButton(radio)

        cbut = QtGui.QPushButton()
        cbut.setMaximumSize(20, 20)
        cbut.setText('')
        cbut.setStyleSheet("QPushButton {background-color: rgb(%d, %d, %d)}"
            % new_cluster.color)
        QtCore.QObject.connect(cbut, QtCore.SIGNAL("clicked()"),
            self.button_cluster_color)

        hlayout.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding,
            QtGui.QSizePolicy.Minimum))

        hlayout.addWidget(cbut)
        layout.insertLayout(layout.count() - 1, hlayout)

        # add the gui elements to the cluster reference
        # so we can access them when we need to
        cbut.cluster_reference = new_cluster
        radio.cluster_reference = new_cluster

        new_cluster.radio = radio
        new_cluster.layout = hlayout
        new_cluster.cbut = cbut
        new_cluster.check = check

        self.clusters.append(new_cluster)
        self.update_active_cluster()

        self.unsaved = True

        return new_cluster

    def delete_cluster(self, cluster=None):
        if cluster == None:
            if ((not self.activeClusterRadioButton()) or
                self.activeClusterRadioButton() == self.radioButton_junk):
                return
            cluster = self.activeClusterRadioButton().cluster_reference

        self.buttonGroup_cluster.removeButton(cluster.radio)
        layout = cluster.layout
        self.labels_container.layout().removeItem(layout)
        for i in range(layout.count()):
            if layout.itemAt(i).widget():
                layout.itemAt(i).widget().close()
                layout.itemAt(i).widget().deleteLater()
        # We have some circular references that will drive the GC mad,
        #so take care of those now
        cluster.layout = None
        cluster.radio.cluster_reference = None
        cluster.cbut.cluster_reference = None
        cluster.radio = None
        cluster.cbut = None
        cluster.check = None

        self.clusters.remove(cluster)

        self.updateFeaturePlot()
        self.updateClusterDetailPlots()

        self.unsaved = True

    def button_add_limit_click(self):
        if not self.spikeset:
            return

        if self.activeClusterRadioButton():
            self.limit_mode = True
            self.limit_data = []
            self.unsaved = True

    def button_delete_limit_click(self):
        if not self.spikeset:
            return

        feature_x = self.ui.comboBox_feature_x.currentText()
        feature_y = self.ui.comboBox_feature_y.currentText()
        feature_x_chan = int(self.ui.comboBox_feature_x_chan.currentText()) - 1
        feature_y_chan = int(self.ui.comboBox_feature_y_chan.currentText()) - 1
        if self.activeClusterRadioButton():
            cluster = self.activeClusterRadioButton().cluster_reference

            if cluster != self.junk_cluster:
                cluster.removeBound(feature_x, feature_x_chan,
                    feature_y, feature_y_chan)
            else:
                cluster.removeBound(feature_x, feature_x_chan,
                    feature_y, feature_y_chan, 'add')
            cluster.calculateMembership(self.spikeset)
            self.update_active_cluster()
            self.updateFeaturePlot()
            self.unsaved = True

    def button_cluster_color(self):
        color = QtGui.QColorDialog.getColor(
            QtGui.QColor(*self.sender().cluster_reference.color),
            self, "ColorDialog")

        self.sender().setStyleSheet(
            "QPushButton {background-color: rgb(%d, %d, %d)}"
            % (color.red(), color.green(), color.blue()))

        self.sender().cluster_reference.color = (color.red(),
            color.green(), color.blue())

        self.updateFeaturePlot()

    def button_add_cluster_click(self):
        self.add_cluster()

    def feature_channel_x_changed(self, index):
        if index == -1 or self.setting_up_combos:
            return

        had_proj = self.redrawing_proj
        self.redrawing_proj = True

        # Get the valid y channel options based on px, cx, py
        current_x = str(self.ui.comboBox_feature_x.currentText())
        current_y = str(self.ui.comboBox_feature_y.currentText())

        valid_y = self.spikeset.featureByName(current_x).valid_y_features(
            int(self.ui.comboBox_feature_x_chan.currentText()) - 1)

        valid_y_chans = [chans for (name, chans) in valid_y
            if name == current_y][0]

        # Set the valid projections for y
        current = self.ui.comboBox_feature_y_chan.currentText()
        if current == '':
            current = 0
        else:
            current = int(current)

        # y_chans is None if it should be all channels for a different
        # feature type, since we have no way of knowing how many there
        # will be apriori, we pass none
        if valid_y_chans == None:
            valid_y_chans = range(0,
                self.spikeset.featureByName(current_y).channels)

        self.ui.comboBox_feature_y_chan.clear()
        for i in valid_y_chans:
            self.ui.comboBox_feature_y_chan.addItem(str(i + 1))
            #if i + 1 == current:
            #    self.ui.comboBox_feature_y_chan.setCurrentIndex(i - index-1)
        #self.ui.comboBox_feature_y_chan.setCurrentIndex(0)
        self.prof_limits_reference = None

        if not had_proj:
            self.redrawing_proj = False

        self.updateFeaturePlot()

    def feature_x_changed(self, index):
        if index == -1:
            return

        self.setting_up_combos = True

        # Set the options in the y feature box accordingly
        current_x = str(self.ui.comboBox_feature_x.currentText())
        valid_x = self.spikeset.featureByName(current_x).valid_x_features()

        # Possible x channels
        self.ui.comboBox_feature_x_chan.clear()
        for i in valid_x:
            self.ui.comboBox_feature_x_chan.addItem(str(i + 1))

        # Possible y features
        # is always going to be 0 first
        valid_y = self.spikeset.featureByName(current_x).valid_y_features(0)
        self.ui.comboBox_feature_y.clear()
        for (name, chans) in valid_y:
            self.ui.comboBox_feature_y.addItem(name)

        self.setting_up_combos = False
        self.feature_y_changed(0)

    def feature_y_changed(self, index):
        if index == -1:
            return

        # Set the options in the y channel box accordingly
        current_x = str(self.ui.comboBox_feature_x.currentText())
        current_y = str(self.ui.comboBox_feature_y.currentText())

        # If x and y features are the same, remove the maximal channel number
        chans = self.spikeset.featureByName(current_x).channels
        index = self.ui.comboBox_feature_x_chan.findText(str(chans))
        if current_x == current_y and index >= 0:
            if self.ui.comboBox_feature_x_chan.currentText() == str(chans):
                self.ui.comboBox_feature_x_chan.setCurrentIndex(chans - 2)
            self.ui.comboBox_feature_x_chan.removeItem(index)
        # If they're not the same we might want it back in
        if current_x != current_y and index == -1:
            self.ui.comboBox_feature_x_chan.insertItem(chans - 1, str(chans))

        valid_y = self.spikeset.featureByName(current_x).valid_y_features(
            int(self.ui.comboBox_feature_x_chan.currentText()) - 1)

        valid_y_chans = [chans for (name, chans) in valid_y
            if name == current_y][0]
        # y_chans is None if it should be all channels for a different feature
        # type, since we have no way of knowing how many there will be apriori,
        #we pass none
        if valid_y_chans == None:
            valid_y_chans = range(0,
                self.spikeset.featureByName(current_y).channels)

        self.ui.comboBox_feature_y_chan.clear()
        for i in valid_y_chans:
            self.ui.comboBox_feature_y_chan.addItem(str(i + 1))

    def feature_channel_y_changed(self, index):
        if index == -1:
            return
        self.prof_limits_reference = None
        if not self.redrawing_proj:
            self.updateFeaturePlot()

    def button_next_feature_click(self):
        if (self.ui.comboBox_feature_y_chan.currentIndex()
            == self.ui.comboBox_feature_y_chan.count() - 1):

            if not (self.ui.comboBox_feature_x_chan.currentIndex()
                == self.ui.comboBox_feature_x_chan.count() - 1):

                # step x to the next projection
                self.redrawing_proj = True
                self.ui.comboBox_feature_x_chan.setCurrentIndex(
                    self.ui.comboBox_feature_x_chan.currentIndex() + 1)

                # but this time we do want to loop around the y projection
                self.ui.comboBox_feature_y_chan.setCurrentIndex(0)
                self.redrawing_proj = False
                self.updateFeaturePlot()
        else:
            self.ui.comboBox_feature_y_chan.setCurrentIndex(
                self.ui.comboBox_feature_y_chan.currentIndex() + 1)

    def button_prev_feature_click(self):
        if self.ui.comboBox_feature_y_chan.currentIndex() == 0:
            if not self.ui.comboBox_feature_x_chan.currentIndex() == 0:
                # step x to the previous projection
                self.redrawing_proj = True
                self.ui.comboBox_feature_x_chan.setCurrentIndex(
                    self.ui.comboBox_feature_x_chan.currentIndex() - 1)
                # but this time we do want to loop around the y projection
                self.ui.comboBox_feature_y_chan.setCurrentIndex(
                        self.ui.comboBox_feature_y_chan.count() - 1)
                self.redrawing_proj = False
                self.updateFeaturePlot()
        else:
            self.ui.comboBox_feature_y_chan.setCurrentIndex(
                self.ui.comboBox_feature_y_chan.currentIndex() - 1)

    def load_spikefile(self, fname):
        if self.unsaved:
            reply = QtGui.QMessageBox.question(self, 'Save',
                "Do you want to save before loading?", QtGui.QMessageBox.Yes |
                QtGui.QMessageBox.No | QtGui.QMessageBox.Cancel,
                QtGui.QMessageBox.Yes)

            if reply == QtGui.QMessageBox.Cancel:
                return
            if reply == QtGui.QMessageBox.Yes:
                self.save()
            if reply == QtGui.QMessageBox.No:
                pass

        print 'Clearing current clusters'
        for cluster in self.clusters[:]:
            self.delete_cluster(cluster)

        self.redrawing_proj = True
        self.spikeset = None
        self.current_feature = None
        self.prof_limits_reference = None
        self.clusters = []
        self.current_filename = None
        self.limit_mode = False
        self.redrawing_proj = False
        self.redrawing_details = False
        self.unsaved = False

        print 'Loading ntt file', fname
        t1 = time.clock()
        self.spikeset = spikeset.load(fname)
        t2 = time.clock()
        print 'Loaded', self.spikeset.N, 'spikes in ', (t2 - t1), 'seconds'

        self.t_bins = np.arange(self.spikeset.time[0],
            self.spikeset.time[-1], 60e6)

        self.t_bin_centers = ((self.t_bins[0:-1] + self.t_bins[1:]) / 2
             - self.spikeset.time[0]) / 60e6

        self.redrawing_proj = True

        # Set the combo boxes to the current feature for now

        self.ui.comboBox_feature_x.clear()
        for name in self.spikeset.featureNames():
            self.ui.comboBox_feature_x.addItem(name)

        # Reset the limits
        self.prof_limits_reference = None

        self.junk_cluster = spikeset.Cluster(self.spikeset)
        self.junk_cluster.color = (255, 0, 0)
        self.junk_cluster.check = self.checkBox_junk
        self.checkBox_junk.setChecked(True)
        self.radioButton_junk.cluster_reference = self.junk_cluster

        self.redrawing_proj = False
        self.updateFeaturePlot()
        self.updateClusterDetailPlots()

        self.current_filename = str(fname)

        self.ui.stackedWidget.setCurrentIndex(0)

    def button_load_click(self):
        fname = QtGui.QFileDialog.getOpenFileName(self,
            'Open ntt file', filter='DotSpike (*.spike);;Neuralynx NTT (*.ntt)')

        if fname:
            self.load_spikefile(str(fname))
            (root, ext) = os.path.splitext(str(fname))
            boundfilename = root + os.extsep + 'bounds'
            if os.path.exists(boundfilename):
                print "Found boundary file", boundfilename
                self.importBounds(boundfilename)

    def updateFeaturePlot(self):
        if self.redrawing_proj:
            return
        if not self.spikeset:
            return

        proj_x = int(self.ui.comboBox_feature_x_chan.currentText()) - 1
        proj_y = int(self.ui.comboBox_feature_y_chan.currentText()) - 1

        feature_x = self.spikeset.featureByName(
            self.ui.comboBox_feature_x.currentText())

        feature_y = self.spikeset.featureByName(
            self.ui.comboBox_feature_y.currentText())

        if feature_x == None or feature_y == None:
            return

        self.mp_proj.axes.hold(False)

        if self.prof_limits_reference == None:
            w = np.array([True] * self.spikeset.N)
            if not self.checkBox_junk.isChecked():
                w[self.junk_cluster.member] = False

            temp = ([np.min(feature_x.data[w, proj_x]),
                    np.max(feature_x.data[w, proj_x])],
                [np.min(feature_y.data[w, proj_y]),
                    np.max(feature_y.data[w, proj_y])])

            w_x = (temp[0][1] - temp[0][0]) * 0.05
            w_y = (temp[1][1] - temp[1][0]) * 0.05
            self.prof_limits_reference = ([temp[0][0] - w_x, temp[0][1] + w_x],
                [temp[1][0] - w_y, temp[1][1] + w_y])
            self.prof_limits = self.prof_limits_reference

        # Do a scatter plot
        if self.ui.radioButton_scatter.isChecked():
            #Plot the unclustered spikes
            if self.ui.checkBox_show_unclustered.isChecked():
                w = np.array([True] * self.spikeset.N)

                if self.ui.checkBox_show_unclustered_exclusive.isChecked():
                    for cluster in self.clusters + [self.junk_cluster]:
                            w[cluster.member] = False

                self.ui.mplwidget_projection.axes.plot(
                    feature_x.data[w, proj_x], feature_y.data[w, proj_y],
                    linestyle='None', marker='.',
                    markersize=int(self.ui.spinBox_markerSize.text()),
                    markerfacecolor='k', markeredgecolor='k')

                self.ui.mplwidget_projection.axes.hold(True)

            # Iterate over clusters
            for cluster in self.clusters + [self.junk_cluster]:
                if not cluster.check.isChecked():
                    continue

                col = map(lambda s: s / 255.0, cluster.color)
                # Plot the cluster spikes
                self.mp_proj.axes.plot(feature_x.data[cluster.member, proj_x],
                        feature_y.data[cluster.member, proj_y], marker='.',
                        markersize=int(self.ui.spinBox_markerSize.text()),
                        markerfacecolor=col, markeredgecolor=col,
                        linestyle='None')
                self.ui.mplwidget_projection.axes.hold(True)

                # Plot refractory spikes
                if self.ui.checkBox_refractory.isChecked():
                    self.mp_proj.axes.plot(
                        feature_x.data[cluster.refractory, proj_x],
                        feature_y.data[cluster.refractory, proj_y],
                        marker='o', markersize=5, markerfacecolor='k',
                        markeredgecolor='k', linestyle='None')

        # Do a density plot
        else:
            w = np.array([False] * self.spikeset.N)
            if self.ui.checkBox_show_unclustered.isChecked():
                w = np.array([True] * self.spikeset.N)
                if self.ui.checkBox_show_unclustered_exclusive.isChecked():
                        for cluster in self.clusters + [self.junk_cluster]:
                                w[cluster.member] = False
            for cluster in [
                cluster for cluster in self.clusters + [self.junk_cluster]
                if cluster.check.isChecked()]:
                    w[cluster.member] = True

            if not np.any(w):
                w[0] = True

            bins_x = np.linspace(self.prof_limits[0][0],
                self.prof_limits[0][1], 100)
            bins_y = np.linspace(self.prof_limits[1][0],
                self.prof_limits[1][1], 100)
            count = np.histogram2d(feature_x.data[w, proj_x],
                feature_y.data[w, proj_y], [bins_x, bins_y])[0]

            if self.ui.radioButton_density.isChecked():
                self.mp_proj.axes.pcolor(bins_x, bins_y, np.transpose(count))
            else:
                self.mp_proj.axes.pcolor(bins_x, bins_y,
                    np.transpose(np.log(count + 1)))

            # Iterate over clusters for refractory spikes
            if self.ui.checkBox_refractory.isChecked():
                self.ui.mplwidget_projection.axes.hold(True)
                for cluster in self.clusters + [self.junk_cluster]:
                    if not cluster.check.isChecked():
                        continue

                    # Plot refractory spikes
                    self.mp_proj.axes.plot(
                        feature_x.data[cluster.refractory, proj_x],
                        feature_y.data[cluster.refractory, proj_y],
                        marker='D', markersize=3, markerfacecolor='w',
                        markeredgecolor='w', linestyle='None')

        self.mp_proj.axes.set_xlim(self.prof_limits[0])
        self.mp_proj.axes.set_ylim(self.prof_limits[1])
        self.mp_proj.axes.set_xticks([])
        self.mp_proj.axes.set_yticks([])

        # Now draw the boundaries
        feature_x = self.ui.comboBox_feature_x.currentText()
        feature_y = self.ui.comboBox_feature_y.currentText()
        feature_x_chan = int(self.ui.comboBox_feature_x_chan.currentText()) - 1
        feature_y_chan = int(self.ui.comboBox_feature_y_chan.currentText()) - 1

        for cluster in self.clusters + [self.junk_cluster]:
            if not cluster.check.isChecked():
                continue

            # Limit boundaries with solid line
            bound_polys = cluster.getBoundPolygon(feature_x, feature_x_chan,
                feature_y, feature_y_chan)

            if bound_polys:
                for bound_poly in bound_polys:
                    bound_poly = np.vstack((bound_poly, bound_poly[0, :]))
                    col = map(lambda s: s / 255.0, cluster.color)
                    for i in range(np.size(bound_poly, 0)):
                        line = mpl.lines.Line2D(bound_poly[i:i + 2, 0],
                            bound_poly[i:i + 2, 1], color=col, linestyle='-')
                        self.mp_proj.axes.add_line(line)

            # Addition boundaries with dashed line
            bound_polys = cluster.getBoundPolygon(feature_x, feature_x_chan,
                feature_y, feature_y_chan, 'add')

            if bound_polys:
                for bound_poly in bound_polys:
                    bound_poly = np.vstack((bound_poly, bound_poly[0, :]))
                    col = map(lambda s: s / 255.0, cluster.color)

                    for i in range(np.size(bound_poly, 0)):
                        line = mpl.lines.Line2D(bound_poly[i:i + 2, 0],
                            bound_poly[i:i + 2, 1], color=col, linestyle='--')

                        self.mp_proj.axes.add_line(line)

        self.mp_proj.draw()
        self.redrawing_proj = False

    def updateClusterDetailPlots(self):
        if self.spikeset == None:
            return

        if self.activeClusterRadioButton():
            w = self.activeClusterRadioButton().cluster_reference.member
            cluster = self.activeClusterRadioButton().cluster_reference
        else:
            w = np.array([False] * self.spikeset.N)
        N = np.sum(w)

        # Draw average waveforms
        for i in range(0, 4):
            if N:
                self.mp_wave.axes[i].errorbar(range(0,
                    np.size(self.spikeset.spikes,1)),
                    np.mean(self.spikeset.spikes[w, :, i], axis=0),
                    np.std(self.spikeset.spikes[w, :, i], axis=0), color='k')
            else:
                self.mp_wave.axes[i].cla()

            self.mp_wave.axes[i].set_xticks([])
            self.mp_wave.axes[i].set_yticks([])
        self.mp_wave.draw()

        # Draw ISI distribution
        if N:
            isi_bin_count = np.histogram(cluster.isi, self.isi_bins)[0]
            self.mp_isi.axes.plot(self.isi_bin_centers, isi_bin_count)
        else:
            self.mp_isi.axes.cla()

        self.mp_isi.axes.set_xscale('log')
        self.mp_isi.axes.set_xlim([0.1, 1e5])

        refractory_line = mpl.lines.Line2D([2, 2],
            self.mp_isi.axes.get_ylim(), color='r', linestyle='--')
        self.mp_isi.axes.add_line(refractory_line)

        burst_line = mpl.lines.Line2D([20, 20], self.mp_isi.axes.get_ylim(),
            color='b', linestyle='--')
        self.mp_isi.axes.add_line(burst_line)

        theta_line = mpl.lines.Line2D([125, 125], self.mp_isi.axes.get_ylim(),
            color='g', linestyle='--')
        self.mp_isi.axes.add_line(theta_line)

        self.mp_isi.axes.set_xticks([1e1, 1e2, 1e3, 1e4])
        self.mp_isi.draw()

        # Now update the stats display
        if N:
            self.ui.label_spike_count.setText('%d' %
                cluster.stats['num_spikes'])

            self.ui.label_mean_rate.setText('%3.2f' %
                cluster.stats['mean_rate'])

            self.ui.label_burst.setText('%3.2f' %
                cluster.stats['burst'])

            self.ui.label_csi.setText('%3.0f' %
                cluster.stats['csi'])

            self.ui.label_refr_count.setText('%3.0f' %
                cluster.stats['refr_count'])

            self.ui.label_refr_fp.setText('%3.2f%%' %
                cluster.stats['refr_fp'])

            self.ui.label_refr_frac.setText('%3.3f%%' %
                (100.0 * cluster.stats['refr_frac']))

            self.ui.label_isolation.setText('%3.0f' %
                cluster.stats['isolation'])
        else:
            self.ui.label_spike_count.setText('')
            self.ui.label_mean_rate.setText('')
            self.ui.label_burst.setText('')
            self.ui.label_csi.setText('')
            self.ui.label_refr_count.setText('')
            self.ui.label_refr_fp.setText('')
            self.ui.label_refr_frac.setText('')
            self.ui.label_isolation.setText('')

        # Drift plot
        if N:
            t = self.spikeset.time[cluster.member]
            p = self.spikeset.featureByName('Peak').data
            p = p[cluster.member, :]

            # this is sort of hacky, wish i had histC
            countsPerBin = np.histogram(t, self.t_bins)[0]

            P = np.zeros((np.size(self.t_bin_centers), np.size(p, 1)))
            for i in range(np.size(p, 1)):
                P[:, i] = np.histogram(t, self.t_bins, weights=p[:, i])[0]
                P[countsPerBin == 0, i] = np.NAN
                P[countsPerBin != 0, i] = (P[countsPerBin != 0, i]
                    / countsPerBin[countsPerBin != 0])

            self.mp_drift.axes.hold(False)
            self.mp_drift.axes.plot(self.t_bin_centers, P)
            ylim = self.mp_drift.axes.get_ylim()
            self.mp_drift.axes.set_ylim([0, ylim[1] * 1.25])
            self.mp_drift.axes.set_xticks([])
        else:
            self.mp_drift.axes.cla()
            self.mp_drift.axes.set_yticks([])
            self.mp_drift.axes.set_xticks([])

        self.mp_drift.draw()

    def onMousePress(self, event):
        if self.limit_mode:
            return
        if (event.xdata != None and event.ydata != None) and event.button == 1:
            self.zoom_active = True
            self.zoom_startpos = (event.xdata, event.ydata, event.x, event.y)

    def onMouseRelease(self, event):
        if event.button == 1 and (event.xdata != None and event.ydata != None):
            if self.zoom_active and not self.limit_mode:
                self.zoom_active = False

                temp = (np.sort([event.xdata, self.zoom_startpos[0]]),
                    np.sort([event.ydata, self.zoom_startpos[1]]))

                # do a sanity check, we cant zoom to an area smaller than
                # say 1/20 of the display
                delta_x = ((temp[0][1] - temp[0][0])
                        / (self.prof_limits[0][1] - self.prof_limits[0][0]))
                delta_y = ((temp[1][1] - temp[1][0])
                        / (self.prof_limits[1][1] - self.prof_limits[1][0]))

                if delta_x < 0.025 or delta_y < 0.025:
                    return

                self.prof_limits = temp

                if self.ui.radioButton_scatter.isChecked():
                    self.mp_proj.axes.set_xlim(self.prof_limits[0])
                    self.mp_proj.axes.set_ylim(self.prof_limits[1])
                    self.mp_proj.draw()
                else:  # if its a density plot we should rebin for now
                    self.updateFeaturePlot()
            if self.limit_mode:
                self.limit_data.append((event.xdata, event.ydata,
                    event.x, event.y))

        # reset bounds on right click, or cancel / complete the bounds
        if event.button == 3 and self.limit_mode:
            if event.xdata != None and event.ydata != None:
                self.limit_data.append((event.xdata, event.ydata,
                    event.x, event.y))
                self.limit_mode = False
                #self.mp_proj.repaint()

                # create an array to describe the bound
                temp = np.array([np.array([point[0], point[1]]) for
                    point in self.limit_data])

                feature_x = str(self.ui.comboBox_feature_x.currentText())
                feature_y = str(self.ui.comboBox_feature_y.currentText())
                feature_x_chan = \
                    int(self.ui.comboBox_feature_x_chan.currentText()) - 1
                feature_y_chan = \
                    int(self.ui.comboBox_feature_y_chan.currentText()) - 1

                print "Adding boundary on", feature_x, feature_x_chan,
                print feature_y, feature_y_chan

                bound = (feature_x, feature_x_chan,
                    feature_y, feature_y_chan, temp)

                clust = self.activeClusterRadioButton().cluster_reference

                if clust != self.junk_cluster:
                    clust.addBound(bound)
                else:
                    print "Junk cluster updated, resetting bounds"
                    self.junk_cluster.add_bounds.append(bound)
                    self.prof_limits_reference = None

                clust.calculateMembership(self.spikeset)
                self.update_active_cluster()
                self.updateFeaturePlot()

        elif event.button == 3:
            self.prof_limits = self.prof_limits_reference
            if self.ui.radioButton_scatter.isChecked():
                self.mp_proj.axes.set_xlim(self.prof_limits[0])
                self.mp_proj.axes.set_ylim(self.prof_limits[1])
                self.mp_proj.draw()
            else:
                self.updateFeaturePlot()

    def onMouseMove(self, event):
        if self.limit_mode:
            self._mouse_move_event = event
            self.mp_proj.repaint()
        if (not self.limit_mode) and self.zoom_active:
            # draw rectangle is in the qt coordinates,
            # so we need a little conversion
            height = self.mp_proj.figure.bbox.height
            x0 = self.zoom_startpos[2]
            y0 = height - self.zoom_startpos[3]
            x1 = event.x
            y1 = height - event.y
            rect = [int(val) for val in min(x0, x1),
                min(y0, y1), abs(x1 - x0), abs(y1 - y0)]
            self.mp_proj.drawRectangle(rect)

    def mp_proj_onPaint(self, event):
        MatplotlibWidget.paintEvent(self.mp_proj, event)
        if self.limit_mode:
            qp = QtGui.QPainter()
            height = self.mp_proj.figure.bbox.height

            qp.begin(self.mp_proj)
            qp.setPen(QtGui.QColor(
                *self.activeClusterRadioButton().cluster_reference.color))
            for i in range(len(self.limit_data) - 1):
                qp.drawLine(self.limit_data[i][2],
                    height - self.limit_data[i][3], self.limit_data[i + 1][2],
                    height - self.limit_data[i + 1][3])
            if self.limit_data:
                qp.drawLine(self.limit_data[-1][2],
                    height - self.limit_data[-1][3], self._mouse_move_event.x,
                    height - self._mouse_move_event.y)
            qp.end()

    def save(self):
        if not (self.current_filename and self.clusters):
            return

        (root, ext) = os.path.splitext(self.current_filename)

        # Save the bounds to a format we can easily read back in later
        outfilename = root + os.extsep + 'bounds'
        outfile = open(outfilename, 'wb')

        save_bounds = [(cluster.color, cluster.bounds, cluster.wave_bounds,
            cluster.add_bounds, cluster.del_bounds)
            for cluster in self.clusters]

        pickle.dump(save_bounds, outfile)
        outfile.close()
        print "Saved bounds to", outfilename

        # Save the cluster membership vectors in matlab format
        outfilename = root + os.extsep + 'cluster'

        cluster_member = np.column_stack(tuple([cluster.member
            for cluster in self.clusters]))
        cluster_stats = [cluster.stats for cluster in self.clusters]

        save_data = {'cluster_id': cluster_member,
            'spike_time': self.spikeset.time}

        for key in cluster_stats[0].keys():
            save_data[key] = [stat[key] for stat in cluster_stats]

        sio.savemat(outfilename, save_data, oned_as='column', appendmat=False)
        outfile.close()

        print "Saved cluster membership to", outfilename
        self.unsaved = False

    def importBounds(self, filename=None):
        if filename == None:
            filename = QtGui.QFileDialog.getOpenFileName(self,
                'Open ntt file', filter='*.bounds')

        if not filename:
            return

        if os.path.exists(filename):
            infile = open(filename, 'rb')
            saved_bounds = pickle.load(infile)
            infile.close()
            print "Found", len(saved_bounds),
            print "bounds to import, creating clusters."
            self.redrawing_details = True
            for (col, bound, wave_bound, add_bound, del_bound) in saved_bounds:
                clust = self.add_cluster(col)
                clust.bounds = bound
                clust.wave_bounds = wave_bound
                clust.add_bounds = add_bound
                clust.del_bounds = del_bound
                clust.calculateMembership(self.spikeset)
            self.redrawing_details = False
            self.updateClusterDetailPlots()
            self.updateFeaturePlot()

    def copyCluster(self):
        if self.activeClusterRadioButton():
            self.redrawing_details = True
            backup = self.activeClusterRadioButton().cluster_reference
            clust = self.add_cluster()
            clust.bounds = backup.bounds
            clust.calculateMembership(self.spikeset)
            self.redrawing_details = False
            self.updateClusterDetailPlots()
            self.updateFeaturePlot()

    def closeEvent(self, event):
        if not self.unsaved:
            event.accept()

        reply = QtGui.QMessageBox.question(self, 'Save',
            "Do you want to save before quitting?", QtGui.QMessageBox.Yes |
            QtGui.QMessageBox.No | QtGui.QMessageBox.Cancel,
            QtGui.QMessageBox.Yes)

        if reply == QtGui.QMessageBox.Cancel:
            event.ignore()
        if reply == QtGui.QMessageBox.Yes:
            event.accept()
            self.save()
        if reply == QtGui.QMessageBox.No:
            event.accept()

    def updateWavecutterPlot(self):
        if self.spikeset == None:
            return
        if not self.activeClusterRadioButton():
            return

        self.wave_limit_mode = False
        self.wave_limit_data = []

        w = self.activeClusterRadioButton().cluster_reference.member
        cluster = self.activeClusterRadioButton().cluster_reference

        chan_no = int(self.ui.spinBox_wavecutter_channel.value()) - 1

        #print "Trying to plot channel", chan_no

        wf_perm = np.random.permutation(np.nonzero(w)[0])
        num_wave = np.min([int(self.ui.lineEdit_wavecutter_count.text()),
            np.size(wf_perm)])

        #print "Trying to plot", num_wave, "waveforms"

        wf_perm = wf_perm[:num_wave]

        if not np.size(wf_perm):
            self.mp_wavecutter.axes.clear()
            self.mp_wavecutter.draw()
            return

        self.mp_wavecutter.axes.hold(False)
        self.mp_wavecutter.axes.plot(np.transpose(
            self.spikeset.spikes[wf_perm, :, chan_no]),
            linewidth=0.5, color=(0.7, 0.7, 0.7))

        if self.ui.checkBox_wavecutter_refractory.isChecked():
            # Overlay refractory spikes
            w = np.logical_and(cluster.refractory, cluster.member)
            wf_perm = np.random.permutation(np.nonzero(w))
            num_wave = np.min([int(self.ui.lineEdit_wavecutter_count.text()),
                np.size(wf_perm)])

            #print "Trying to plot", num_wave, "refractory waveforms"
            wf_perm = wf_perm[0, :num_wave]
            if np.size(wf_perm):
                self.mp_wavecutter.axes.hold(True)
                self.mp_wavecutter.axes.plot(np.transpose(
                    self.spikeset.spikes[wf_perm, :, chan_no]),
                    color=(0.3, 0.3, 0.3))

        # Plot boundaries
        if cluster.wave_bounds:
            for (chan, sample, lower, upper) in cluster.wave_bounds:
                if chan != chan_no:
                    continue
                line = mpl.lines.Line2D([sample, sample],
                    [lower, upper], color=(1, 0, 0),
                    linestyle='-', linewidth=4)

                self.mp_wavecutter.axes.add_line(line)

        self.mp_wavecutter.axes.set_xticks([])
        self.mp_wavecutter.axes.set_yticks([])
        ylim = self.mp_wavecutter.axes.get_ylim()
        ylim = np.max(np.abs(ylim))
        self.mp_wavecutter.axes.set_ylim([-ylim, ylim])

        line = mpl.lines.Line2D(self.mp_wavecutter.axes.get_xlim(),
                [0, 0], color=(0, 0, 0), linestyle='-', linewidth=1)
        self.mp_wavecutter.axes.add_line(line)

        self.mp_wavecutter.axes.set_xlim([-0.25, -0.75 +
            np.size(self.spikeset.spikes,1)])
        self.mp_wavecutter.draw()

    def wavecutter_onMousePress(self, event):
        pass

    def wavecutter_onMouseRelease(self, event):
        if not self.spikeset:
            return
        if not self.activeClusterRadioButton():
            return

        if event.button == 1 and (event.xdata != None and event.ydata != None):
            if not self.wave_limit_mode:
                return

            if not self.wave_limit_data:
                self.wave_limit_data = (event.xdata, event.ydata,
                    event.x, event.y)
            else:
                self.wave_limit_mode = False
                chan_no = int(self.ui.spinBox_wavecutter_channel.value()) - 1

                bound = (chan_no, int(np.round(self.wave_limit_data[0])),
                    min(self.wave_limit_data[1], event.ydata),
                    max(self.wave_limit_data[1], event.ydata))

                self.wave_limit_data = []

                cluster = self.activeClusterRadioButton().cluster_reference
                cluster.wave_bounds.append(bound)

                cluster.calculateMembership(self.spikeset)
                self.updateClusterDetailPlots()
                self.updateWavecutterPlot()
                self.unsaved = True

        if event.button == 3:
            pass

    def wavecutter_onMouseMove(self, event):
        if not (self.wave_limit_mode and self.wave_limit_data):
            return

        height = self.mp_wavecutter.figure.bbox.height
        width = self.mp_wavecutter.figure.bbox.width
        #x0 = self.wave_limit_data[2]
        offset = width * 0.5 / (np.size(self.spikeset.spikes,1) - 0.5)
        width = width - offset
        x0 = np.round(self.wave_limit_data[0]) * width / \
                (np.size(self.spikeset.spikes,1) - 1) + offset
        y0 = height - self.wave_limit_data[3]
        x1 = x0
        y1 = height - event.y
        rect = [int(val) for val in min(x0, x1),
            min(y0, y1), abs(x1 - x0), abs(y1 - y0)]
        self.mp_wavecutter.drawRectangle(rect)

    def wavecutter_remove_limits(self):
        if not self.spikeset:
            return
        if not self.activeClusterRadioButton():
            return

        chan_no = int(self.ui.spinBox_wavecutter_channel.value()) - 1
        cluster = self.activeClusterRadioButton().cluster_reference
        cluster.wave_bounds = [bound for bound in
            cluster.wave_bounds if bound[0] != chan_no]

        cluster.calculateMembership(self.spikeset)
        self.updateClusterDetailPlots()
        self.updateWavecutterPlot()
        self.unsaved = True

    def wavecutter_add_limit(self):
        self.wave_limit_mode = True

    def updateOutlierPlot(self):
        if self.spikeset == None:
            return
        if not self.activeClusterRadioButton():
            return

        cluster = self.activeClusterRadioButton().cluster_reference
        if not cluster.mahal_valid:
            cluster.calculateMahal(self.spikeset)

        if cluster.mahal_valid:
            # Mahal histogram
            m1 = np.min(cluster.mahal)
            m2 = np.max(cluster.mahal)
            bins = np.logspace(np.log10(m1), np.log10(m2))
            centers = (bins[0:-1] + bins[1:]) / 2

            count = np.histogram(cluster.mahal, bins)[0]
            count = count.astype(float) / np.sum(count)

            self.mp_outlier.axes.hold(False)
            self.mp_outlier.axes.plot(centers, count)
            self.mp_outlier.axes.hold(True)

            m_ref = cluster.mahal[cluster.refractory[cluster.member]]
            if np.size(m_ref):
                count_ref = np.histogram(
                    cluster.mahal[cluster.refractory[cluster.member]], bins)[0]

                count_ref = (count_ref.astype(float) * np.max(count) /
                        np.max(count_ref))

                self.mp_outlier.axes.plot(centers, count_ref, 'k')

#            chi = spikeset.chi2f(centers, np.size(self.spikeset.spikes,1) *
#                    np.size(self.spikeset.spikes,2) - 1)
            chi = scipy.stats.chi2.pdf(centers,
                    np.size(self.spikeset.spikes,1) *
                    np.size(self.spikeset.spikes,2) - 1)

            chi = chi / np.sum(chi)

            self.mp_outlier.axes.plot(centers, chi, 'r--')

            self.mp_outlier.axes.set_xscale('log')
            self.mp_outlier.axes.set_xlim([m1, m2])

            endpoint = scipy.stats.chi2.ppf(
                ((float(cluster.stats['num_spikes']) - 1.0) /
                cluster.stats['num_spikes']), np.size(self.spikeset.spikes,1) *
                np.size(self.spikeset.spikes,2) - 1)


            endpoint_line = mpl.lines.Line2D([endpoint, endpoint],
                self.mp_outlier.axes.get_ylim(), color='g', linestyle='--')

            self.mp_outlier.axes.add_line(endpoint_line)
        else:
            self.mp_outlier.axes.cla()

        self.mp_outlier.draw()

if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = PyClustMainWindow()
    myapp.show()
    #myapp.load_ntt('Sample2.ntt')
    sys.exit(app.exec_())
