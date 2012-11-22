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
mpl.use('Qt4Agg')

import pickle
import os
import sys
import time
import struct

from gui import Ui_MainWindow
from gui_sessiontracker import Ui_SessionTrackerDialog

import spikeset
import featurewidget
import unique_colors
import mmodel
import boundaries


class PyClustSessionTrackerDialog(QtGui.QDialog):

    def __del__(self):
        print "Session Tracker Dialog being destroy"

    def __init__(self, spikeset_ref, cluster_ref,
            spikeset_comp, cluster_comp, parent=None):
        QtGui.QDialog.__init__(self, parent)

        self.spikeset = spikeset_ref
        self.clusters = cluster_ref

        self.spikeset_comp = spikeset_comp
        self.clusters_comp = cluster_comp

        self.p1 = spikeset_ref.featureByName('Peak').data
        self.p2 = spikeset_comp.featureByName('Peak').data

        self.step_r = 0
        self.matches = [-2] * len(self.clusters)

        self.ui = Ui_SessionTrackerDialog()
        self.ui.setupUi(self)

        self.figure = self.ui.mplwidget.figure
        self.figure.clear()
        self.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0, top=1, left=0.0, right=1)

        self.member_ref = np.any(np.column_stack(tuple([cluster.member
            for cluster in self.clusters])), axis=1)
        self.member_comp = np.any(np.column_stack(tuple([cluster.member
            for cluster in self.clusters_comp])), axis=1)

        self.updatePlot()
        self.group = QtGui.QButtonGroup(self)
        self.group.setExclusive(True)
        self.radio = [self.ui.radioButton_nomatch]
        self.group.addButton(self.radio[0])
        self.group.setId(self.radio[0], -2)

        # Create the radio buttons
        layout = self.ui.frame.layout()
        for i, clust in enumerate(self.clusters_comp):
            radio = QtGui.QRadioButton()
            radio.setText("Cluster %d" % i)
            layout.insertWidget(i + 1, radio)
            self.radio.append(radio)
            self.group.addButton(radio)
            self.group.setId(radio, i)

        self.radio[0].setChecked(True)
        QtCore.QObject.connect(self.group, QtCore.SIGNAL("buttonClicked(int)"),
                self.updatePlot)

    def getMatches(self):
        return ((self.spikeset.subject, self.spikeset_comp.subject),
                (self.spikeset.session, self.spikeset_comp.session),
                self.matches)

    def accept(self):
        self.matches[self.step_r] = self.group.checkedId()
        for radio in self.radio[1:]:
            radio.setEnabled(not self.group.id(radio) in self.matches)

        if self.step_r + 1 == len(self.clusters):
            QtGui.QDialog.accept(self)
        else:
            self.step_r += 1
            self.radio[0].setChecked(True)
            self.updatePlot()

    def updatePlot(self, index=-2):
        gcol = (0.5, 0.5, 0.5)
        dgcol = (0.3, 0.3, 0.3)
        clust1 = self.clusters[self.step_r]
        ccol = map(lambda s: s / 255.0, clust1.color)
        if index != -2:
            clust2 = self.clusters_comp[index]

        counter = 0
        for proj_x in range(0, 4):
            for proj_y in range(proj_x + 1, 4):
                counter = counter + 1
                ax = self.figure.add_subplot(2, 3, counter)
                ax.cla()
                # Plot things that belong to no cluster
                ax.plot(self.p1[np.logical_not(self.member_ref), proj_x],
                        self.p1[np.logical_not(self.member_ref), proj_y],
                                 marker='.', markersize=1,
                                 markerfacecolor=gcol,
                                 markeredgecolor=gcol,
                                 linestyle='None', zorder=0)
                ax.plot(self.p2[np.logical_not(self.member_comp), proj_x],
                        self.p2[np.logical_not(self.member_comp), proj_y],
                                 marker='.', markersize=1,
                                 markerfacecolor=gcol,
                                 markeredgecolor=gcol,
                                 linestyle='None', zorder=0)

                # Plot things belonging to other clusters
                w = np.logical_and(self.member_ref,
                        np.logical_not(clust1.member))
                ax.plot(self.p1[w, proj_x], self.p1[w, proj_y],
                                 marker='.', markersize=1,
                                 markerfacecolor=dgcol,
                                 markeredgecolor=dgcol,
                                 linestyle='None', zorder=1)
                if index != -2:
                    w = np.logical_and(self.member_comp,
                            np.logical_not(clust2.member))
                else:
                    w = self.member_comp
                ax.plot(self.p2[w, proj_x], self.p2[w, proj_y],
                                 marker='.', markersize=1,
                                 markerfacecolor=dgcol,
                                 markeredgecolor=dgcol,
                                 linestyle='None', zorder=1)

                # Plot things belonging to these clusters
                ax.plot(self.p1[clust1.member, proj_x],
                        self.p1[clust1.member, proj_y],
                                 marker='.', markersize=3,
                                 markerfacecolor=ccol,
                                 markeredgecolor=ccol,
                                 linestyle='None', zorder=3)
                if index != -2:
                    ax.plot(self.p2[clust2.member, proj_x],
                            self.p2[clust2.member, proj_y],
                                 marker='.', markersize=3,
                                 markerfacecolor='k',
                                 markeredgecolor='k',
                                 linestyle='None', zorder=3)

                # Draw the boundaries
                for bound in clust1.getBoundaries('Peak', proj_x,
                            'Peak', proj_y):
                    bound.draw(ax, color=ccol, linestyle='--')
                if index != -2:
                    for bound in clust2.getBoundaries('Peak', proj_x,
                            'Peak', proj_y):
                        bound.draw(ax, color='k', linestyle='-.')

        self.ui.mplwidget.draw()


class PyClustMainWindow(QtGui.QMainWindow):

    def switch_to_wavecutter(self):
        self.ui.stackedWidget.setCurrentIndex(1)
        self.updateWavecutterPlot()

    def switch_to_maindisplay(self):
        self.ui.stackedWidget.setCurrentIndex(0)
        self.updateFeaturePlot()

    @QtCore.pyqtSlot('bool')
    def on_actionAutotrim_triggered(self, checked=None):
        active = self.activeClusterRadioButton()
        if active and active.cluster_reference != self.junk_cluster:
            clust = active.cluster_reference
            self.trim_cluster = spikeset.Cluster(self.spikeset)
            self.trim_cluster.color = clust.color
            self.trim_cluster.bounds = clust.bounds
            self.trim_cluster.calculateMembership(self.spikeset)

            self.ui.stackedWidget.setCurrentIndex(2)

            self.ui.pushButton_autotrim_apply.setEnabled(False)
            self.ui.pushButton_autotrim_cancel.setEnabled(False)
            self.ui.mplwidget_autotrim.figure.clear()

            current_x = str(self.ui.comboBox_feature_x.currentText())
            clust.autotrim(
                self.spikeset, fname=current_x,
                canvas=self.ui.mplwidget_autotrim)

            self.updateClusterDetailPlots()

            self.ui.pushButton_autotrim_apply.setEnabled(True)
            self.ui.pushButton_autotrim_cancel.setEnabled(True)

    @QtCore.pyqtSlot('bool')
    def on_actionRefractory_triggered(self, checked=None):
        self.mp_proj.setRefractory(checked)
        self.ui.checkBox_refractory.blockSignals(True)
        self.ui.checkBox_refractory.setChecked(checked)
        self.ui.checkBox_refractory.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionScatter_triggered(self, checked=None):
        if checked:
            self.mp_proj.setPlotType(-2)
            self.ui.radioButton_scatter.blockSignals(True)
            self.ui.radioButton_scatter.setChecked(True)
            self.ui.radioButton_scatter.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionDensity_triggered(self, checked=None):
        if checked:
            self.mp_proj.setPlotType(-3)
            self.ui.radioButton_density.blockSignals(True)
            self.ui.radioButton_density.setChecked(True)
            self.ui.radioButton_density.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionLog_Density_triggered(self, checked=None):
        if checked:
            self.mp_proj.setPlotType(-4)
            self.ui.radioButton_log_density.blockSignals(True)
            self.ui.radioButton_log_density.setChecked(True)
            self.ui.radioButton_log_density.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionMarker_Size1_triggered(self, checked=None):
        if checked:
            self.mp_proj.setMarkerSize(1)
            self.ui.spinBox_markerSize.blockSignals(True)
            self.ui.spinBox_markerSize.setValue(1)
            self.ui.spinBox_markerSize.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionMarker_Size3_triggered(self, checked=None):
        if checked:
            self.mp_proj.setMarkerSize(3)
            self.ui.spinBox_markerSize.blockSignals(True)
            self.ui.spinBox_markerSize.setValue(3)
            self.ui.spinBox_markerSize.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionMarker_Size5_triggered(self, checked=None):
        if checked:
            self.mp_proj.setMarkerSize(5)
            self.ui.spinBox_markerSize.blockSignals(True)
            self.ui.spinBox_markerSize.setValue(5)
            self.ui.spinBox_markerSize.blockSignals(False)

    @QtCore.pyqtSlot('bool')
    def on_actionEliptical_triggered(self, checked=None):
        self.mp_proj.setBoundaryElliptical(checked)
        self.ui.checkBox_ellipse.blockSignals(True)
        self.ui.checkBox_ellipse.setChecked(checked)
        self.ui.checkBox_ellipse.blockSignals(False)

    def autotrim_apply(self):
        self.trim_cluster = None

        self.updateFeaturePlot()
        self.updateClusterDetailPlots()
        self.ui.stackedWidget.setCurrentIndex(0)

    def autotrim_cancel(self):
        clust = self.activeClusterRadioButton().cluster_reference
        clust.bounds = self.trim_cluster.bounds
#        for bound in self.trim_cluster.bounds:
#            clust.addBoundary(bound)
        clust.calculateMembership(self.spikeset)

        self.trim_cluster = None

        self.updateFeaturePlot()
        self.updateClusterDetailPlots()
        self.ui.stackedWidget.setCurrentIndex(0)

    def switch_to_sessiontracker(self):
        #Prompt for a filename for the comparison spikeset
        fname = QtGui.QFileDialog.getOpenFileName(self,
            'Open ntt file',
            filter='DotSpike (*.spike);;Neuralynx NTT (*.ntt)')

        if not fname:
            return

        (root, ext) = os.path.splitext(str(fname))
        boundfilename = root + os.extsep + 'bounds'
        if not os.path.exists(boundfilename):
            # Should probably display a message at some point
            return

        print "Loading comparison spikeset"
        infile = open(boundfilename, 'rb')
        special = pickle.load(infile)
        #compmatches = pickle.load(infile)
        saved_bounds = pickle.load(infile)
        infile.close()

        print "Loading comparison clusters"
        ss_comp = spikeset.load(str(fname))
        ss_comp.calculateFeatures(special)

        clusters_comp = []

        for (col, bound, wave_bound, add_bound, del_bound) in saved_bounds:
            clust = spikeset.Cluster(ss_comp)
            clust.bounds = bound
            clust.wave_bounds = wave_bound
            clust.add_bounds = add_bound
            clust.del_bounds = del_bound
            clust.calculateMembership(ss_comp)
            clusters_comp.append(clust)

        dialog = PyClustSessionTrackerDialog(self.spikeset, self.clusters,
               ss_comp, clusters_comp)
        dialog.showMaximized()
        if dialog.exec_():
            matches = dialog.getMatches()
            if not all(x == -2 for x in matches):
                self.matches = matches
            print "Matches", self.matches

    def __init__(self, parent=None):
        QtGui.QMainWindow.__init__(self, parent)

        self.redrawing_proj = True

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        # There is no spikeset or feature set on load
        self.spikeset = None
        self.current_feature = None
        self.clusters = []
        self.matches = None
        self.current_filename = None
        self.junk_cluster = None

        self.limit_mode = False
        self.redrawing_proj = False
        self.redrawing_details = False
        self.unsaved = False

        # Create action groups for the mnu items
        self.ui.agroup_marker = QtGui.QActionGroup(self, exclusive=True)
        self.ui.agroup_marker.addAction(self.ui.actionMarker_Size1)
        self.ui.agroup_marker.addAction(self.ui.actionMarker_Size3)
        self.ui.agroup_marker.addAction(self.ui.actionMarker_Size5)

        self.ui.agroup_scatter = QtGui.QActionGroup(self, exclusive=True)
        self.ui.agroup_scatter.addAction(self.ui.actionScatter)
        self.ui.agroup_scatter.addAction(self.ui.actionDensity)
        self.ui.agroup_scatter.addAction(self.ui.actionLog_Density)

        # Connect the handlers
        self.ui.pushButton_identify.clicked.connect(
                self.switch_to_sessiontracker)

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
            lambda x: self.updateWavecutterPlot() if x == 0
                else self.updateOutlierPlot())

        self.ui.pushButton_autotrim_apply.clicked.connect(self.autotrim_apply)
        self.ui.pushButton_autotrim_cancel.clicked.connect(
            self.autotrim_cancel)

        self.ui.label_subjectid.setText('')
        self.ui.label_session.setText('')
        self.ui.label_fname.setText('')

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
        hlayout.addWidget(self.ui.checkBox_show_unclustered_exclusive)

        layout.insertLayout(0, hlayout)

        # Create the projection widget
        self.ui.verticalLayout_3.removeWidget(self.ui.mplwidget_projection)
        self.mp_proj = featurewidget.ProjectionWidget()
        self.mp_proj.setAutoFillBackground(True)
        self.mp_proj.setObjectName("mplwidget_projection")
        self.ui.verticalLayout_3.addWidget(self.mp_proj)

        # Connect the relevant signals/slots
        QtCore.QObject.connect(self.mp_proj,
                QtCore.SIGNAL("featureRedrawRequired()"),
                self.updateFeaturePlot)

        # Signals for unclustered show checkbox
        self.ui.checkBox_show_unclustered.stateChanged.connect(
                self.mp_proj.setShowUnclustered)

        # Signals for unclustered exclusive checkbox
        self.ui.checkBox_show_unclustered_exclusive.stateChanged.connect(
                self.mp_proj.setUnclusteredExclusive)

        # Signals for marker size
        self.ui.spinBox_markerSize.valueChanged.connect(
                self.mp_proj.setMarkerSize)

        # Signal for polygon boundary drawn
        QtCore.QObject.connect(self.mp_proj,
                QtCore.SIGNAL("polygonBoundaryDrawn(PyQt_PyObject)"),
                self.addBoundary)
        QtCore.QObject.connect(self.mp_proj,
                QtCore.SIGNAL("ellipseBoundaryDrawn(PyQt_PyObject)"),
                self.addBoundary)

        # Create shorter handles to important widgets
        self.mp_wave = self.ui.mplwidget_waveform
        self.mp_isi = self.ui.mplwidget_isi
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

        self.mp_wavecutter.figure.clear()
        self.mp_wavecutter.figure.set_facecolor(bgcolor)

        self.mp_outlier.figure.clear()
        self.mp_outlier.figure.set_facecolor(bgcolor)

        self.mp_drift.figure.clear()
        self.mp_drift.figure.set_facecolor(bgcolor)

        # Set up autotrim plot
        self.ui.mplwidget_autotrim.figure.clear()
        self.ui.mplwidget_autotrim.figure.set_facecolor(bgcolor)

        self.wave_limit_mode = False
        self.mp_wavecutter.mpl_connect('button_press_event',
            self.wavecutter_onMousePress)
        self.mp_wavecutter.mpl_connect('button_release_event',
            self.wavecutter_onMouseRelease)
        self.mp_wavecutter.mpl_connect('motion_notify_event',
            self.wavecutter_onMouseMove)

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

        # Set up autotrim plot
        self.ui.mplwidget_autotrim.figure.subplots_adjust(bottom=0, top=1,
                left=0, right=1, hspace=0.01, wspace=0.01)

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

        # list containing the cluster check/radio/color buttons
        self.cluster_ui_buttons = []

    def hide_show_all_clusters(self, hidden):
        self.redrawing_proj = True
        radio = self.activeClusterRadioButton()
        for cluster, ui in zip(self.clusters, self.cluster_ui_buttons):
            if ui[0] == radio:
                continue
#            if radio != self.radioButton_junk and \
#                cluster == self.activeClusterRadioButton().cluster_reference:
#                continue
            ui[3].setChecked(not hidden)
#            cluster.check.setChecked(not hidden)
            if hidden:
                # dont show this on show all
                self.checkBox_junk.setChecked(False)
        self.ui.checkBox_show_unclustered.setChecked(not hidden)
        self.redrawing_proj = False
        self.updateFeaturePlot()

    # When we switch clusters, correctly enabled/disable cluster
    # checkboxes, and tell the projection widget to stop drawing
    def update_active_cluster(self):
        self.updateClusterDetailPlots()
        for ui in self.cluster_ui_buttons:
            ui[3].setEnabled(True)  # enable all check boxes
        if self.activeClusterRadioButton() == self.radioButton_junk:
            return
        for ui in self.cluster_ui_buttons:
            if ui[0] == self.activeClusterRadioButton():
                check = ui[3]
                check.setEnabled(False)
                check.setChecked(True)
        self.mp_proj.stopMouseAction()

    # Add a new cluster by generating a color, creating GUI elements, etc.
    def add_cluster(self, color=None):
        new_cluster = spikeset.Cluster(self.spikeset)
        # If we're provided with a cluster color, use it
        if color:
            new_cluster.color = color
        # Otherwise try to intelligently generate one using the existing
        # color list
        else:
            color_list = [map(lambda x: x / 255.0, cluster.color) for
                    cluster in self.clusters]

            # If this is the first color, we can pick whatever we want
            if color_list == []:
                new_cluster.color = (50, 170, 170)
            else:
                c = unique_colors.newcolor(color_list)
                new_cluster.color = tuple(map(lambda x: int(x * 255), c))

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

#        new_cluster.radio = radio
#        new_cluster.layout = hlayout
#        new_cluster.cbut = cbut
#        new_cluster.check = check

        self.cluster_ui_buttons.append((radio, hlayout, cbut, check))
        self.clusters.append(new_cluster)
        self.update_active_cluster()

        self.unsaved = True

        return new_cluster

    @QtCore.pyqtSlot('bool')
    def on_actionDelete_Cluster_triggered(self, checked=None):
        self.delete_cluster()

    def delete_cluster(self, cluster=None):
        if cluster is None:
            if ((not self.activeClusterRadioButton()) or
                    self.activeClusterRadioButton() == self.radioButton_junk):
                return
            cluster = self.activeClusterRadioButton().cluster_reference

        if self.matches is not None:
            retval = QtGui.QMessageBox.question(self, "Warning",
                    "Deleting clusters will invalidate (i.e. delete) "
                    "session tracking data. Continue?",
                    QtGui.QMessageBox.Yes | QtGui.QMessageBox.No,
                    QtGui.QMessageBox.No)
            print retval
            if retval != QtGui.QMessageBox.Yes:
                return
            else:
                print "Clearing cluster matches"
                self.matches = None

        for ui_container in self.cluster_ui_buttons:
            if ui_container[0].cluster_reference == cluster:
                radio = ui_container[0]
                layout = ui_container[1]
                cbut = ui_container[2]

                self.buttonGroup_cluster.removeButton(radio)
                self.labels_container.layout().removeItem(layout)
                for i in range(layout.count()):
                    if layout.itemAt(i).widget():
                        layout.itemAt(i).widget().close()
                        layout.itemAt(i).widget().deleteLater()

                radio.cluster_reference = None
                cbut.cluster_reference = None

                self.cluster_ui_buttons.remove(ui_container)
                break

        self.clusters.remove(cluster)

        self.updateFeaturePlot()
        self.updateClusterDetailPlots()

        self.unsaved = True

    # Tell the projection widget to start drawing a boundary
    @QtCore.pyqtSlot('bool')
    def on_actionAdd_Limit_triggered(self, checked=None):
        if not self.spikeset:
            return

        if self.activeClusterRadioButton():
            self.mp_proj.drawBoundary(
                    self.activeClusterRadioButton().cluster_reference.color)

    # Do something with the boundary that gets returned
    def addBoundary(self, bound):
        if not self.activeClusterRadioButton():
            return

        clust = self.activeClusterRadioButton().cluster_reference
        if clust != self.junk_cluster:
            clust.addBoundary(bound)
        else:
            print "Junk cluster updated, resetting bounds"
            self.junk_cluster.add_bounds.append(bound)
            self.mp_proj.resetLimits()

        clust.calculateMembership(self.spikeset)
        self.update_active_cluster()
        self.updateFeaturePlot()

    # Delete the active cluster boundary on the current projection
    @QtCore.pyqtSlot('bool')
    def on_actionDelete_Limit_triggered(self, checked=None):
        if not self.spikeset:
            return

        feature_x = self.ui.comboBox_feature_x.currentText()
        feature_y = self.ui.comboBox_feature_y.currentText()
        feature_x_chan = int(self.ui.comboBox_feature_x_chan.currentText()) - 1
        feature_y_chan = int(self.ui.comboBox_feature_y_chan.currentText()) - 1
        if self.activeClusterRadioButton():
            cluster = self.activeClusterRadioButton().cluster_reference

            if cluster != self.junk_cluster:
                cluster.removeBound((feature_x, feature_y),
                    (feature_x_chan, feature_y_chan))
            else:
                cluster.removeBound((feature_x, feature_y),
                    (feature_x_chan, feature_y_chan), 'add')
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

    @QtCore.pyqtSlot('bool')
    def on_actionAdd_Cluster_triggered(self, checked=None):
        self.add_cluster()

    def feature_channel_x_changed(self, index):
        if index == -1:
            return
        if self.setting_up_combos:
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
        if valid_y_chans is None:
            valid_y_chans = range(0,
                self.spikeset.featureByName(current_y).channels)

        self.ui.comboBox_feature_y_chan.clear()
        for i in valid_y_chans:
            self.ui.comboBox_feature_y_chan.addItem(str(i + 1))

        if not had_proj:
            self.redrawing_proj = False

        self.mp_proj.setChanX(int(
                self.ui.comboBox_feature_x_chan.currentText()) - 1)

    def feature_x_changed(self, index):
        if index == -1:
            return

        self.redrawing_proj = True
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

        # Update the feature widget
        self.mp_proj.setFeatureX(current_x)

        self.mp_proj.setChanX(int(
                self.ui.comboBox_feature_x_chan.currentText()) - 1)

        self.setting_up_combos = False
        self.redrawing_proj = False
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
        if valid_y_chans is None:
            valid_y_chans = range(0,
                self.spikeset.featureByName(current_y).channels)

        self.ui.comboBox_feature_y_chan.clear()
        for i in valid_y_chans:
            self.ui.comboBox_feature_y_chan.addItem(str(i + 1))

        self.mp_proj.setFeatureY(current_y)

    def feature_channel_y_changed(self, index):
        if index == -1:
            return
        self.mp_proj.setChanY(int(
            self.ui.comboBox_feature_y_chan.currentText()) - 1)

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
                self.on_actionSave_triggered()
            if reply == QtGui.QMessageBox.No:
                pass

        self.redrawing_proj = True

        self.matches = None
        print ''
        print 'Trying to load file: ', fname
        print 'Clearing current clusters'
        for cluster in self.clusters[:]:
            self.delete_cluster(cluster)

        self.redrawing_proj = True
        self.spikeset = None
        self.current_feature = None
        self.mp_proj.resetLimits()
        self.clusters = []
        self.junk_cluster = None
        self.current_filename = None
        self.limit_mode = False
        self.redrawing_details = False
        self.unsaved = False
        self.ui.checkBox_show_unclustered.setChecked(True)
        self.ui.checkBox_show_unclustered_exclusive.setChecked(True)

        print 'Loading file', fname
        t1 = time.clock()
        self.spikeset = spikeset.load(fname)
        t2 = time.clock()
        print 'Loaded', self.spikeset.N, 'spikes in ', (t2 - t1), 'seconds'

        self.ui.label_subjectid.setText(self.spikeset.subject)
        self.ui.label_session.setText(self.spikeset.session)
        self.ui.label_fname.setText(
            os.path.splitext(os.path.split(fname)[1])[0])
        self.curfile = fname

        self.t_bins = np.arange(self.spikeset.time[0],
            self.spikeset.time[-1], 60e6)

        self.t_bin_centers = ((self.t_bins[0:-1] + self.t_bins[1:]) / 2
             - self.spikeset.time[0]) / 60e6

        # Reset the limits
        self.junk_cluster = spikeset.Cluster(self.spikeset)
        self.junk_cluster.color = (255, 0, 0)
        self.junk_cluster.check = self.checkBox_junk
        self.checkBox_junk.setChecked(False)
        self.radioButton_junk.cluster_reference = self.junk_cluster

        # Add some points to the junk cluster since true spikes are rarely >1mv
        d = np.max(self.spikeset.featureByName('Peak').data, axis=0) + 10
        # similarly check for peak < 10
        d2 = np.min(self.spikeset.featureByName('Peak').data, axis=0)
        d2[d2 >= 0] = -0.001
        jthresh = 1e3
        for i in range(d.shape[0]):
            for j in range(i + 1, d.shape[0]):
                jbounds = np.zeros((4, 2))
                # bounds between 1mv and max peak
                jbounds[0, :] = [min([d[i], jthresh]), min([d[j], jthresh])]
                jbounds[1, :] = [max([d[i], jthresh]), min([d[j], jthresh])]
                jbounds[2, :] = [max([d[i], jthresh]), max([d[j], jthresh])]
                jbounds[3, :] = [min([d[i], jthresh]), max([d[j], jthresh])]
                jbound = boundaries.BoundaryPolygon2D(('Peak', 'Peak'),
                        (i, j), jbounds)

                # bounds between negative peaks and zero
                jbounds = np.zeros((6, 2))
                jbounds[0, :] = [d2[i], d2[j]]  # bottom left
                jbounds[1, :] = [d[i], d2[j]]   # bottom right
                jbounds[2, :] = [d[i], 0]       # top right
                jbounds[3, :] = [0, 0]          # origin
                jbounds[4, :] = [0, d[j]]       # top middle
                jbounds[5, :] = [d2[i], d[j]]   # top left
                jbound = boundaries.BoundaryPolygon2D(('Peak', 'Peak'),
                        (i, j), jbounds)

                self.junk_cluster.add_bounds.append(jbound)

        self.junk_cluster.calculateMembership(self.spikeset)
        self.mp_proj.resetLimits()

        # Autoload any existing boundary files for this data set
        imported_bounds = False

        boundfilename = str(fname) + os.extsep + 'bounds'
        if os.path.exists(boundfilename) and (not imported_bounds):
            print "Found boundary file", boundfilename
            self.importBounds(boundfilename)
            imported_bounds = True

        (root, ext) = os.path.splitext(str(fname))
        boundfilename = root + os.extsep + 'bounds'
        if os.path.exists(boundfilename):
            print "Found boundary file", boundfilename
            self.importBounds(boundfilename)
            imported_bounds = True

       # Set the combo boxes to the current feature for now
        self.ui.comboBox_feature_x.clear()
        for name in self.spikeset.featureNames():
            self.ui.comboBox_feature_x.addItem(name)

        self.redrawing_proj = False
        self.updateFeaturePlot()
        self.updateClusterDetailPlots()

        self.current_filename = str(fname)

        self.ui.stackedWidget.setCurrentIndex(0)

        # see if there is a modified klustakwik cluster file to load
        kkwikfilename = fname + os.extsep + 'ackk'
        if (not imported_bounds) and os.path.exists(kkwikfilename):
            print "Found KKwik cluster file", kkwikfilename, ':',
            # Everything is little endian
            # A .ackk record is as folows:
            # Header:
            # uint32 + n x char - subject string
            # uint32 + n x char - session date string
            # uint64 - num spikes
            # uint16 x N - spike cluster IDs
            f = open(kkwikfilename, 'rb')
            subjectstr = spikeset.readStringFromBinary(f)
            datestr = spikeset.readStringFromBinary(f)
            num_samps, = struct.unpack('<Q', f.read(8))

            if num_samps != self.spikeset.N:
                print "Sample counts don't match up, invalid .ackk file"
                print num_samps
                print self.spikeset.N
                f.close()
                return
            elif subjectstr != self.spikeset.subject:
                print "Subject lines don't match up, invalid .ackk file"
                f.close()
                return
            elif datestr != self.spikeset.session:
                print "Session date strings don't match up, invalid .ackk file"
                f.close()
                return

            labels = np.fromfile(f, dtype=np.dtype('<h'), count=num_samps)
            f.close()

            # get a list of cluster numbers
            k = np.unique(labels)
            print np.size(k) - 1, 'components.'
            colors = unique_colors.unique_colors_hsv(np.size(k) - 1)

            # create the clusters for each component in the KK file.
            if np.size(k) > 1:
                print "Importing clusters from .ackk",
                for i in k:
                    if i == 0:
                        continue
                    print i,
                    sys.stdout.flush()
                    cluster = self.add_cluster(
                        color=tuple([int(a * 255) for a in colors[i - 1]]))

                    cluster.membership_model.append(
                            mmodel.PrecalculatedLabelsMembershipModel(
                                labels, [i]))
                    cluster.calculateMembership(self.spikeset)
                print "."

                self.update_active_cluster()
                self.updateFeaturePlot()

    @QtCore.pyqtSlot('bool')
    def on_actionOpen_triggered(self, checked=None):
        fname = QtGui.QFileDialog.getOpenFileName(self,
            'Open ntt file',
            filter='DotSpike (*.spike);;Neuralynx NTT (*.ntt)')

        if fname:
            self.load_spikefile(str(fname))

    def updateFeaturePlot(self):
        if self.redrawing_proj:
            return
        if not self.spikeset or not self.junk_cluster:
            return

        check_boxes = [ui[3] for ui in self.cluster_ui_buttons] + [
                self.junk_cluster.check]
        self.mp_proj.updatePlot(self.spikeset,
                self.clusters, self.junk_cluster, check_boxes)

        self.mp_proj.draw()
        self.redrawing_proj = False

    def updateClusterDetailPlots(self):
        if self.spikeset is None:
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
                    np.size(self.spikeset.spikes, 1)),
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

    @QtCore.pyqtSlot('bool')
    def on_actionSave_triggered(self, checked=None):
        if not (self.current_filename and self.clusters):
            return

        (root, ext) = os.path.splitext(self.current_filename)

        # Save the bounds to a format we can easily read back in later
        #outfilename = root + os.extsep + 'bounds'
        outfilename = self.current_filename + os.extsep + 'bounds'

        outfile = open(outfilename, 'wb')

        versionstr = "0.0.1"
        if versionstr == "0.0.1":
            dumping = dict()
            dumping['feature_special'] = self.spikeset.feature_special
            dumping['matches'] = self.matches

            sb = [{'color': cluster.color, 'bounds': cluster.bounds,
                'wave_bounds': cluster.wave_bounds,
                'add_bounds':  cluster.add_bounds,
                'del_bounds': cluster.del_bounds,
                'mmodel': cluster.membership_model}
                        for cluster in self.clusters]
            dumping['clusters'] = sb

            pickle.dump(versionstr, outfile)
            pickle.dump(dumping, outfile)
        else:
            # save special info about the features, such as PCA coeffs
            pickle.dump(self.spikeset.feature_special, outfile)
            pickle.dump(self.matches, outfile)

            save_bounds = [(cluster.color, cluster.bounds, cluster.wave_bounds,
                cluster.add_bounds, cluster.del_bounds,
                cluster.membership_model)
                for cluster in self.clusters]

            pickle.dump(save_bounds, outfile)
            outfile.close()
            print "Saved bounds to", outfilename

        # Save the cluster membership vectors in matlab format
        #outfilename = root + os.extsep + 'mat'
        outfilename = self.current_filename + os.extsep + 'mat'

        cluster_member = np.column_stack(tuple([cluster.member
            for cluster in self.clusters]))
        cluster_stats = [cluster.stats for cluster in self.clusters]

        save_data = {'cluster_id': cluster_member,
            'spike_time': self.spikeset.time,
            'subject': self.spikeset.subject,
            'session': self.spikeset.session}
        if self.matches is not None:
            save_data['match_subject'] = self.matches[0][1]
            save_data['match_session'] = self.matches[1][1]
            save_data['matches'] = self.matches[2]

        for key in cluster_stats[0].keys():
            save_data[key] = [stat[key] for stat in cluster_stats]

        sio.savemat(outfilename, save_data, oned_as='column', appendmat=False)
        outfile.close()

        print "Saved cluster membership to", outfilename
        self.unsaved = False

    @QtCore.pyqtSlot('bool')
    def on_actionImport_Bound_triggered(self, checked=None):
        filename = QtGui.QFileDialog.getOpenFileName(self,
            'Open ntt file', filter='*.bounds')
        if not filename:
            return
        self.importBounds(filename)

    def importBounds(self, filename=None):
        if os.path.exists(filename):
            infile = open(filename, 'rb')

            versionstr = pickle.load(infile)
            if versionstr == "0.0.1":
                print "Saved bounds version 0.0.1"

                dumped = pickle.load(infile)
                self.spikeset.calculateFeatures(dumped['feature_special'])
                self.matches = dumped['matches']

                self.redrawing_details = True
                print "Founds", len(dumped['clusters']),
                print "bounds to import, creating clusters."
                for cluster in dumped['clusters']:
                    clust = self.add_cluster(cluster['color'])
                    clust.bounds = cluster['bounds']
                    clust.wave_bounds = cluster['wave_bounds']
                    clust.add_bounds = cluster['add_bounds']
                    clust.del_bounds = cluster['del_bounds']
                    clust.membership_model = cluster['mmodel']
                    clust.calculateMembership(self.spikeset)
                self.redrawing_details = False
                self.updateClusterDetailPlots()
                self.updateFeaturePlot()

            else:
                print "Saved bounds version 0.0.0"

                # the old, very inflexible way of doing things
                special = versionstr
                #special = pickle.load(infile)
                self.spikeset.calculateFeatures(special)

                self.matches = pickle.load(infile)

                saved_bounds = pickle.load(infile)
                infile.close()
                print "Found", len(saved_bounds),
                print "bounds to import, creating clusters."
                self.redrawing_details = True
                for (col, bound, wave_bound, add_bound, del_bound) \
                        in saved_bounds:
                    clust = self.add_cluster(col)
                    clust.bounds = bound
                    clust.wave_bounds = wave_bound
                    clust.add_bounds = add_bound
                    clust.del_bounds = del_bound
                    clust.calculateMembership(self.spikeset)
                self.redrawing_details = False
                self.updateClusterDetailPlots()
                self.updateFeaturePlot()

    @QtCore.pyqtSlot('bool')
    def on_actionCopy_Cluster_triggered(self, checked=None):
        if self.activeClusterRadioButton():
            self.redrawing_details = True
            backup = self.activeClusterRadioButton().cluster_reference
            clust = self.add_cluster()
            clust.bounds = backup.bounds
            clust.membership_model = backup.membership_model
            clust.add_bounds = backup.add_bounds
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
            self.on_actionSave_triggered()
        if reply == QtGui.QMessageBox.No:
            event.accept()

    def updateWavecutterPlot(self):
        if self.spikeset is None:
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
            np.size(self.spikeset.spikes, 1)])
        self.mp_wavecutter.draw()

    def wavecutter_onMousePress(self, event):
        pass

    def wavecutter_onMouseRelease(self, event):
        if not self.spikeset:
            return
        if not self.activeClusterRadioButton():
            return

        if event.button == 1 and (event.xdata is not None and
                                  event.ydata is not None):
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
        offset = width * 0.5 / (np.size(self.spikeset.spikes, 1) - 0.5)
        width = width - offset
        x0 = np.round(self.wave_limit_data[0]) * width / \
                (np.size(self.spikeset.spikes, 1) - 1) + offset
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
        if self.spikeset is None:
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
            chi = scipy.stats.chi2.pdf(centers, 4)
#                    np.size(self.spikeset.spikes, 1) *
#                    np.size(self.spikeset.spikes, 2) - 1)

            chi = chi / np.sum(chi)

            self.mp_outlier.axes.plot(centers, chi, 'r--')

            self.mp_outlier.axes.set_xscale('log')
            self.mp_outlier.axes.set_xlim([m1, m2])

            endpoint = scipy.stats.chi2.ppf(
                ((float(cluster.stats['num_spikes']) - 1.0) /
                cluster.stats['num_spikes']), 4)
            #np.size(self.spikeset.spikes, 1)
#                * np.size(self.spikeset.spikes, 2) - 1)

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
