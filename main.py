# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 23:43:03 2012

@author: Bernard
"""

#!/usr/bin/python -d


from PyQt4 import QtCore, QtGui

import numpy as np
import scipy.io as sio
import matplotlib as mpl
from matplotlibwidget import MatplotlibWidget

import pickle, os, sys, time

from gui import Ui_MainWindow
import spikeset
import features

class PyClustMainWindow(QtGui.QMainWindow):    
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
        QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"),  self.button_load_click)
        QtCore.QObject.connect(self.ui.pushButton_copy_cluster, QtCore.SIGNAL("clicked()"), self.copyCluster)
        QtCore.QObject.connect(self.ui.pushButton_import, QtCore.SIGNAL("clicked()"), self.importBounds)
        QtCore.QObject.connect(self.ui.comboBox_feature_x_chan, QtCore.SIGNAL("currentIndexChanged(int)"), self.feature_channel_x_changed)
        QtCore.QObject.connect(self.ui.comboBox_feature_y_chan, QtCore.SIGNAL("currentIndexChanged(int)"), self.feature_channel_y_changed)
        QtCore.QObject.connect(self.ui.pushButton_next_projection, QtCore.SIGNAL("clicked()"), self.button_next_feature_click)
        QtCore.QObject.connect(self.ui.pushButton_previous_projection, QtCore.SIGNAL("clicked()"), self.button_prev_feature_click)
        QtCore.QObject.connect(self.ui.buttonGroup_scatter_plot_type, QtCore.SIGNAL("buttonClicked(int)"), self.updateFeaturePlot)
        QtCore.QObject.connect(self.ui.spinBox_markerSize, QtCore.SIGNAL("valueChanged(int)"), self.updateFeaturePlot)
        QtCore.QObject.connect(self.ui.pushButton_addLimit, QtCore.SIGNAL("clicked()"), self.button_add_limit_click)
        QtCore.QObject.connect(self.ui.pushButton_deleteLimit, QtCore.SIGNAL("clicked()"), self.button_delete_limit_click)
        QtCore.QObject.connect(self.ui.checkBox_refractory, QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)
                
        # Set up the cluster list area
        self.labels_container = QtGui.QWidget()
        self.ui.scrollArea_cluster_list.setWidget(self.labels_container)
        
        layout = QtGui.QVBoxLayout(self.labels_container)
        self.buttonGroup_cluster = QtGui.QButtonGroup()
        self.buttonGroup_cluster.setExclusive(True)        
        QtCore.QObject.connect(self.buttonGroup_cluster, QtCore.SIGNAL("buttonClicked(int)"), self.update_active_cluster)
        
        layout.addItem(QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))        
        QtCore.QObject.connect(self.ui.pushButton_add_cluster, QtCore.SIGNAL("clicked()"), self.button_add_cluster_click)
        QtCore.QObject.connect(self.ui.pushButton_delete_cluster, QtCore.SIGNAL("clicked()"), self.delete_cluster )
        
        # Add the unclustered entry
        layout = self.labels_container.layout()
        
        hlayout = QtGui.QHBoxLayout()                    
        
        self.ui.checkBox_show_unclustered = QtGui.QCheckBox(self.labels_container)
        self.ui.checkBox_show_unclustered.setChecked(True)
        QtCore.QObject.connect(self.ui.checkBox_show_unclustered, QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)
        hlayout.addWidget(self.ui.checkBox_show_unclustered)
        
        hlayout.addWidget(QtGui.QLabel('Unclustered', parent=self.labels_container))
        
        hlayout.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
        
        self.ui.checkBox_show_unclustered_exclusive = QtGui.QCheckBox(self.labels_container)
        self.ui.checkBox_show_unclustered_exclusive.setChecked(True)
        QtCore.QObject.connect(self.ui.checkBox_show_unclustered_exclusive, QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)        
        hlayout.addWidget(self.ui.checkBox_show_unclustered_exclusive)
        
        layout.insertLayout(0, hlayout)
        
        # Create short hands for the plot widgets
        self.mp_wave = self.ui.mplwidget_waveform
        self.mp_isi = self.ui.mplwidget_isi
        self.mp_proj = self.ui.mplwidget_projection
                
        pal = self.palette().window().color()
        bgcolor = (pal.red() / 255.0, pal.blue() / 255.0, pal.green() / 255.0)
        
        self.mp_wave.figure.clear()
        self.mp_wave.figure.set_facecolor(bgcolor)
        
        self.mp_isi.figure.clear()
        self.mp_isi.figure.set_facecolor(bgcolor)
        
        self.mp_proj.figure.clear()
        self.mp_proj.figure.set_facecolor(bgcolor)
        
        self.zoom_active = False
        self.mp_proj.mpl_connect('button_press_event', self.onMousePress)
        self.mp_proj.mpl_connect('button_release_event', self.onMouseRelease)
        self.mp_proj.mpl_connect('motion_notify_event', self.onMouseMove)        
        
        #print self.mp_proj.paintEvent
        self._backup_mp_proj_paintEvent = self.mp_proj.paintEvent
        #print self.mp_proj.paintEvent
        #print self.mp_proj.paintEvent.func_code.co_varnames
        #print self.mp_proj.paintEvent.func_code.co_argcount
        
        self.mp_proj.paintEvent = lambda e: self.mp_proj_onPaint(e)
                
        self.ui.comboBox_feature_x.clear()
        self.ui.comboBox_feature_y.clear()
        self.ui.comboBox_feature_x_chan.clear()
        self.ui.comboBox_feature_y_chan.clear()
        
        self.activeClusterRadioButton = self.buttonGroup_cluster.checkedButton
        
        self.redrawing_proj = False
        
        # Set up waveform plot plot axes
        self.mp_wave.figure.clear()
        self.mp_wave.figure.subplots_adjust(hspace=0.0001, wspace=0.0001, bottom=0, top=1, left=0.0, right=1)
        self.mp_wave.axes = [None, None, None, None]
        for i in range(0,4):
            if i == 0:
                self.mp_wave.axes[i] = self.mp_wave.figure.add_subplot(2,2,i+1)
            else:
                self.mp_wave.axes[i] = self.mp_wave.figure.add_subplot(2,2,i+1, sharey = self.mp_wave.axes[0])
            self.mp_wave.axes[i].hold(False)
            self.mp_wave.axes[i].set_xticks([])
            self.mp_wave.axes[i].set_yticks([])
        
        # Set up ISI plot axes

        self.isi_bins = np.logspace(np.log10(0.1), np.log10(1e5), 100)
        self.isi_bin_centers = (self.isi_bins[0:-1] + self.isi_bins[1:])/2

        self.mp_isi.figure.clear()        
        self.mp_isi.figure.subplots_adjust(hspace=0.0001, wspace=0.0001, bottom=0.15, top=1, left=0.0, right=1)                
        self.mp_isi.axes = self.mp_isi.figure.add_subplot(1,1,1)
        self.mp_isi.axes.hold(False)
        self.mp_isi.axes.set_xscale('log')
        self.mp_isi.axes.set_xlim([0.1, 1e5])
        refractory_line = mpl.lines.Line2D([2, 2], self.mp_isi.axes.get_ylim(), color='r', linestyle='--')
        self.mp_isi.axes.add_line(refractory_line)
        burst_line = mpl.lines.Line2D([20, 20], self.mp_isi.axes.get_ylim(), color='b', linestyle='--')
        self.mp_isi.axes.add_line(burst_line)
        theta_line = mpl.lines.Line2D([125, 125], self.mp_isi.axes.get_ylim(), color='g', linestyle='--')
        self.mp_isi.axes.add_line(theta_line)    
        self.mp_isi.axes.set_xticks([1e1, 1e2, 1e3, 1e4])
        self.mp_isi.draw()
        
        # Set up the feature axes
        self.mp_proj.figure.clear()
        self.mp_proj.figure.subplots_adjust(hspace=0.0001, wspace=0.0001, bottom=0.075, top=1, left=0.075, right=1)
        self.mp_proj.axes = self.mp_proj.figure.add_subplot(1,1,1)        
        
        # Clear the stats labels
        
        self.ui.label_spike_count.setText('')
        self.ui.label_mean_rate.setText('')
        self.ui.label_burst.setText('')
        self.ui.label_csi.setText('')
        self.ui.label_refr_count.setText('')
        self.ui.label_refr_fp.setText('')
        self.ui.label_refr_frac.setText('')
        self.ui.label_isolation.setText('')        
        
        
    def update_active_cluster(self):
        self.ui.pushButton_addLimit.setEnabled(True)
        self.updateClusterDetailPlots()
        for cluster in self.clusters:
            cluster.check.setEnabled(True)
        self.activeClusterRadioButton().cluster_reference.check.setEnabled(False)
        self.activeClusterRadioButton().cluster_reference.check.setChecked(True)
        
    def add_cluster(self, color=None):
        new_cluster = spikeset.Cluster(self.spikeset)
        if color:
            new_cluster.color = color
                
        layout = self.labels_container.layout()
        
        hlayout = QtGui.QHBoxLayout()                    
        
        check = QtGui.QCheckBox()
        check.setChecked(True)
        QtCore.QObject.connect(check, QtCore.SIGNAL("stateChanged(int)"), self.updateFeaturePlot)
        
        radio = QtGui.QRadioButton()
        radio.setChecked(True)
        spacer = QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum)
        
        hlayout.addWidget(check)
        hlayout.addItem(spacer)
        hlayout.addWidget(radio)
        hlayout.setSizeConstraint(QtGui.QLayout.SetMaximumSize)
        self.buttonGroup_cluster.addButton(radio)
        
        cbut = QtGui.QPushButton()
        cbut.setMaximumSize(20,20)
        cbut.setText('')
        cbut.setStyleSheet("QPushButton {background-color: rgb(%d, %d, %d)}" % new_cluster.color)
        QtCore.QObject.connect(cbut, QtCore.SIGNAL("clicked()"), self.button_cluster_color)

        hlayout.addItem(QtGui.QSpacerItem(40, 20, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Minimum))
                
        hlayout.addWidget(cbut)        
        layout.insertLayout(layout.count()-1, hlayout)
        
        # add the gui elements to the cluster reference so we can access them when we need to
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
        
    def delete_cluster(self, cluster = None):
        if cluster == None:
            if self.activeClusterRadioButton() == 0: return
            cluster = self.activeClusterRadioButton().cluster_reference
            
        self.buttonGroup_cluster.removeButton(cluster.radio)
        layout = cluster.layout
        self.labels_container.layout().removeItem(layout)
        for i in range(layout.count()):
            if layout.itemAt(i).widget():
                layout.itemAt(i).widget().close()
                layout.itemAt(i).widget().deleteLater()
        # We have some circular references that will drive the GC mad, so take care of those now
        cluster.layout = None
        cluster.radio.cluster_reference = None
        cluster.cbut.cluster_reference = None
        cluster.radio = None
        cluster.cbut = None
        cluster.check = None
                        
        self.clusters.remove(cluster)
    
        self.updateFeaturePlot()
        self.updateClusterDetailPlots()
        
        self.ui.pushButton_addLimit.setEnabled(False)
        
        self.unsaved = True
        
    def button_add_limit_click(self):
        self.limit_mode = True
        self.ui.pushButton_addLimit.setEnabled(False)
        self.limit_data = []
        self.unsaved = True
        
    def button_delete_limit_click(self):
        feature_x = self.ui.comboBox_feature_x.currentText()
        feature_y = self.ui.comboBox_feature_y.currentText()
        feature_x_chan = int(self.ui.comboBox_feature_x_chan.currentText())-1
        feature_y_chan = int(self.ui.comboBox_feature_y_chan.currentText())-1
        if self.activeClusterRadioButton():
            cluster= self.activeClusterRadioButton().cluster_reference
            cluster.removeBound(feature_x, feature_x_chan, feature_y, feature_y_chan)
            self.activeClusterRadioButton().cluster_reference.calculateMembership(self.spikeset)
            self.update_active_cluster()
            self.updateFeaturePlot()
            self.unsaved = True
        
    def button_cluster_color(self):
        color=QtGui.QColorDialog.getColor(QtGui.QColor(*self.sender().cluster_reference.color), self, "ColorDialog")
        self.sender().setStyleSheet("QPushButton {background-color: rgb(%d, %d, %d)}" % (color.red(),color.green(),color.blue()))
        self.sender().cluster_reference.color = (color.red(),color.green(),color.blue())
        
        self.updateFeaturePlot()
        
    def button_add_cluster_click(self):
        self.add_cluster()
        
    def feature_channel_x_changed(self, index):
        if index == -1: return
        had_proj = self.redrawing_proj
        self.redrawing_proj = True
        # Set the valid projections for y
        current = self.ui.comboBox_feature_y_chan.currentText()
        if current == '': current = 0
        else: current = int(current)        
        self.ui.comboBox_feature_y_chan.clear()
        for i in range(index+2, self.current_feature.channels+1):
            self.ui.comboBox_feature_y_chan.addItem(str(i))
            if i == current:
                self.ui.comboBox_feature_y_chan.setCurrentIndex(i - index-2)
        self.prof_limits_reference = None
        if not had_proj: self.redrawing_proj = False
        self.updateFeaturePlot()
                
    def feature_channel_y_changed(self, index):
        if index == -1: return
        self.prof_limits_reference = None
        if not self.redrawing_proj:
            self.updateFeaturePlot()
    
    def button_next_feature_click(self):
        if self.ui.comboBox_feature_y_chan.currentIndex() == self.ui.comboBox_feature_y_chan.count()-1:
            if not self.ui.comboBox_feature_x_chan.currentIndex() == self.ui.comboBox_feature_x_chan.count()-1:
                # step x to the next projection
                self.redrawing_proj = True
                self.ui.comboBox_feature_x_chan.setCurrentIndex(self.ui.comboBox_feature_x_chan.currentIndex()+1)
                # but this time we do want to loop around the y projection                
                self.ui.comboBox_feature_y_chan.setCurrentIndex(0)
                self.redrawing_proj = False
                self.updateFeaturePlot()
        else:
            self.ui.comboBox_feature_y_chan.setCurrentIndex(self.ui.comboBox_feature_y_chan.currentIndex()+1)
    
    def button_prev_feature_click(self):
        if self.ui.comboBox_feature_y_chan.currentIndex() == 0:
            if not self.ui.comboBox_feature_x_chan.currentIndex() == 0:
                # step x to the previous projection
                self.redrawing_proj = True
                self.ui.comboBox_feature_x_chan.setCurrentIndex(self.ui.comboBox_feature_x_chan.currentIndex()-1)
                # but this time we do want to loop around the y projection                
                self.ui.comboBox_feature_y_chan.setCurrentIndex(self.ui.comboBox_feature_y_chan.count()-1)                
                self.redrawing_proj = False
                self.updateFeaturePlot()
        else:
            self.ui.comboBox_feature_y_chan.setCurrentIndex(self.ui.comboBox_feature_y_chan.currentIndex()-1)

    def load_ntt(self, fname): 
        if self.unsaved:
            reply = QtGui.QMessageBox.question(self, 'Save', "Do you want to save before loading?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No | QtGui.QMessageBox.Cancel, QtGui.QMessageBox.Yes)
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
        self.spikeset = spikeset.loadNtt(fname)
        t2 = time.clock()
        print 'Loaded', self.spikeset.N, 'spikes in ', (t2-t1), 'seconds'
        self.current_feature = features.Feature_Peak(self.spikeset)

        self.redrawing_proj = True

        # Set the combo boxes to the current feature for now
        self.ui.comboBox_feature_x.clear();
        self.ui.comboBox_feature_x.addItem(self.current_feature.name)
        self.ui.comboBox_feature_y.clear();
        self.ui.comboBox_feature_y.addItem(self.current_feature.name)

        # Set the channel combo box options
        self.ui.comboBox_feature_x_chan.clear();
        self.ui.comboBox_feature_y_chan.clear();

        for i in range(1, self.current_feature.channels):
            self.ui.comboBox_feature_x_chan.addItem(str(i))
            
        # Reset the limits
        self.prof_limits_reference = None

        self.redrawing_proj = False
        self.updateFeaturePlot()
        self.updateClusterDetailPlots()
        
        self.current_filename = str(fname)
        
    def button_load_click(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open ntt file', filter='*.ntt')
        if fname:
            self.load_ntt(fname)
            (root, ext) = os.path.splitext(str(fname))
            boundfilename = root + os.extsep + 'bounds'
            self.importBounds(boundfilename)
        
    def updateFeaturePlot(self):
        if self.current_feature == None: return
        if self.redrawing_proj: return

        proj_x = int(self.ui.comboBox_feature_x_chan.currentText())-1
        proj_y = int(self.ui.comboBox_feature_y_chan.currentText())-1

        self.ui.mplwidget_projection.axes.hold(False)
        
        if self.prof_limits_reference == None:
            temp = ([np.min(self.current_feature.data[:,proj_x]), np.max(self.current_feature.data[:,proj_x])], [np.min(self.current_feature.data[:,proj_y]), np.max(self.current_feature.data[:,proj_y])])
            w_x = (temp[0][1] - temp[0][0]) * 0.05
            w_y = (temp[1][1] - temp[1][0]) * 0.05
            self.prof_limits_reference = ([temp[0][0] - w_x, temp[0][1] + w_x], [temp[1][0] - w_y, temp[1][1] + w_y])
            self.prof_limits = self.prof_limits_reference
            
        #in_bounds = np.logical_and(np.logical_and(self.current_feature.data[:,proj_y] >= self.prof_limits[1][0], self.current_feature.data[:,proj_y] <= self.prof_limits[1][1]), np.logical_and(self.current_feature.data[:,proj_x] >= self.prof_limits[0][0], self.current_feature.data[:,proj_x] <= self.prof_limits[0][1]))

        # Do a scatter plot
        if self.ui.radioButton_scatter.isChecked():
            
            #Plot the unclustered spikes
            if self.ui.checkBox_show_unclustered.isChecked():
                w = np.array([True] * self.spikeset.N)
                #w[np.logical_not(in_bounds)] = False
                
                if self.ui.checkBox_show_unclustered_exclusive.isChecked():
                    for cluster in self.clusters:
                        w[cluster.member] = False
                self.ui.mplwidget_projection.axes.plot(self.current_feature.data[w,proj_x], self.current_feature.data[w,proj_y], linestyle='None', marker='.', markersize=int(self.ui.spinBox_markerSize.text()), markerfacecolor='k', markeredgecolor='k')
                self.ui.mplwidget_projection.axes.hold(True)
            
            # Iterate over clusters
            for cluster in self.clusters:
                if not cluster.check.isChecked(): continue
                col = map(lambda s: s / 255.0, cluster.color)
                self.mp_proj.axes.plot(self.current_feature.data[cluster.member,proj_x], self.current_feature.data[cluster.member,proj_y], marker='.', markersize=int(self.ui.spinBox_markerSize.text()), markerfacecolor=col, markeredgecolor=col, linestyle='None')
                self.ui.mplwidget_projection.axes.hold(True)
                
                # Plot refractory spikes
                if self.ui.checkBox_refractory.isChecked():
                    self.mp_proj.axes.plot(self.current_feature.data[cluster.refractory,proj_x], self.current_feature.data[cluster.refractory,proj_y], marker='o', markersize=5, markerfacecolor='k', markeredgecolor='k', linestyle='None')
                
        # Do a density plot
        else:            
            w = np.array([False] * self.spikeset.N)
            if self.ui.checkBox_show_unclustered.isChecked():
                w = np.array([True] * self.spikeset.N)
                if self.ui.checkBox_show_unclustered_exclusive.isChecked():
                    for cluster in self.clusters:
                        w[cluster.member] = False
            for cluster in [cluster for cluster in self.clusters if cluster.check.isChecked()]:
                w[cluster.member] = True
                
            if not np.any(w):
                w[0] = True            
                
            bins_x = np.linspace(self.prof_limits[0][0], self.prof_limits[0][1], 100)
            bins_y = np.linspace(self.prof_limits[1][0], self.prof_limits[1][1], 100)
            count = np.histogram2d(self.current_feature.data[w,proj_x], self.current_feature.data[w,proj_y], [bins_x, bins_y])[0]
            
            if self.ui.radioButton_density.isChecked():
                self.mp_proj.axes.pcolor(bins_x, bins_y, np.transpose(count))
            else:
                self.mp_proj.axes.pcolor(bins_x, bins_y, np.transpose(np.log(count+1)))
                
            # Iterate over clusters for refractory spikes
            if self.ui.checkBox_refractory.isChecked():
                self.ui.mplwidget_projection.axes.hold(True)
                for cluster in self.clusters:
                    if not cluster.check.isChecked(): continue
                    # Plot refractory spikes
                    self.mp_proj.axes.plot(self.current_feature.data[cluster.refractory,proj_x], self.current_feature.data[cluster.refractory,proj_y], marker='D', markersize=3, markerfacecolor='w', markeredgecolor='w', linestyle='None')                
                
        self.mp_proj.axes.set_xlim(self.prof_limits[0])
        self.mp_proj.axes.set_ylim(self.prof_limits[1])
        
        # Now draw the boundaries
        feature_x = self.ui.comboBox_feature_x.currentText()
        feature_y = self.ui.comboBox_feature_y.currentText()
        feature_x_chan = int(self.ui.comboBox_feature_x_chan.currentText())-1
        feature_y_chan = int(self.ui.comboBox_feature_y_chan.currentText())-1        
        for cluster in self.clusters:
            if not cluster.check.isChecked(): continue
            bound_poly = cluster.getBoundPolygon(feature_x, feature_x_chan, feature_y, feature_y_chan)
            if not bound_poly == None:
                col = map(lambda s: s / 255.0, cluster.color)
                for i in range(np.size(bound_poly,0)):
                    line = mpl.lines.Line2D(bound_poly[i:i+2,0], bound_poly[i:i+2,1], color=col, linestyle='-')
                    self.mp_proj.axes.add_line(line)
                line = mpl.lines.Line2D([bound_poly[0,0], bound_poly[-1,0]], [bound_poly[0,1], bound_poly[-1,1]], color=col, linestyle='-')
                self.mp_proj.axes.add_line(line)
                    #self.mp_proj.axes.line()
                    
        self.mp_proj.draw()
        self.redrawing_proj = False
        
    def updateClusterDetailPlots(self):               
        if self.current_feature == None: return
        
        if self.activeClusterRadioButton():
            w = self.activeClusterRadioButton().cluster_reference.member
            cluster = self.activeClusterRadioButton().cluster_reference
        else:
            w = np.array([False] * self.spikeset.N)
        N = np.sum(w)

        # Draw average waveforms
        for i in range(0,4):
            if N:
                self.mp_wave.axes[i].errorbar(range(1,33), np.mean(self.spikeset.spikes[w,:,i], axis=0), np.std(self.spikeset.spikes[w,:,i], axis=0), color='k')
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
        refractory_line = mpl.lines.Line2D([2, 2], self.mp_isi.axes.get_ylim(), color='r', linestyle='--')
        self.mp_isi.axes.add_line(refractory_line)
        burst_line = mpl.lines.Line2D([20, 20], self.mp_isi.axes.get_ylim(), color='b', linestyle='--')
        self.mp_isi.axes.add_line(burst_line)
        theta_line = mpl.lines.Line2D([125, 125], self.mp_isi.axes.get_ylim(), color='g', linestyle='--')
        self.mp_isi.axes.add_line(theta_line)    
        self.mp_isi.axes.set_xticks([1e1, 1e2, 1e3, 1e4])
        self.mp_isi.draw()
        
        # Now update the stats display
        if N:            
            self.ui.label_spike_count.setText('%d' % cluster.stats['num_spikes'])
            self.ui.label_mean_rate.setText('%3.2f' % cluster.stats['mean_rate'])
            self.ui.label_burst.setText('%3.2f' % cluster.stats['burst'])
            self.ui.label_csi.setText('%3.0f' % cluster.stats['csi'])
            self.ui.label_refr_count.setText('%3.0f' % cluster.stats['refr_count'])
            self.ui.label_refr_fp.setText('%3.2f%%' % cluster.stats['refr_fp'])
            self.ui.label_refr_frac.setText('%3.3f%%' % (100.0 * cluster.stats['refr_frac']))
            self.ui.label_isolation.setText('%3.0f' % cluster.stats['isolation'])
        else:
            self.ui.label_spike_count.setText('')
            self.ui.label_mean_rate.setText('')
            self.ui.label_burst.setText('')
            self.ui.label_csi.setText('')
            self.ui.label_refr_count.setText('')
            self.ui.label_refr_fp.setText('')
            self.ui.label_refr_frac.setText('')
            self.ui.label_isolation.setText('')
        
    def onMousePress(self, event):
        if self.limit_mode: return
        if not (event.xdata == None or event.ydata == None) and event.button == 1:
            self.zoom_active = True
            self.zoom_startpos = (event.xdata, event.ydata, event.x, event.y)
        
    def onMouseRelease(self, event):
        if event.button == 1 and not (event.xdata == None or event.ydata == None):
            if self.zoom_active and not self.limit_mode:                
                self.zoom_active = False
                temp = (np.sort([event.xdata, self.zoom_startpos[0]]), np.sort([event.ydata, self.zoom_startpos[1]]))

                # do a sanity check, we cant zoom to an area smaller than say, 1/20 of the display                
                delta_x = (temp[0][1] - temp[0][0]) / (self.prof_limits[0][1] - self.prof_limits[0][0])
                delta_y = (temp[1][1] - temp[1][0]) / (self.prof_limits[1][1] - self.prof_limits[1][0])
                if delta_x < 0.025 or delta_y < 0.025: return
                
                self.prof_limits = temp
                
                if self.ui.radioButton_scatter.isChecked():
                    self.mp_proj.axes.set_xlim(self.prof_limits[0])
                    self.mp_proj.axes.set_ylim(self.prof_limits[1])
                    self.mp_proj.draw()
                else: # if its a density plot we should rebin for now
                    self.updateFeaturePlot()
            if self.limit_mode:
                self.limit_data.append( (event.xdata, event.ydata, event.x, event.y) )
                    
        # reset bounds on right click, or cancel / complete the bounds 
        if event.button == 3:
            if self.limit_mode and not (event.xdata == None or event.ydata == None):
                self.limit_data.append( (event.xdata, event.ydata, event.x, event.y) )
                self.limit_mode = False
                self.ui.pushButton_addLimit.setEnabled(True)
                self.mp_proj.repaint()
                # create an array to describe the bound
                temp = np.array([np.array([point[0], point[1]]) for point in self.limit_data])
                feature_x = str(self.ui.comboBox_feature_x.currentText())
                feature_y = str(self.ui.comboBox_feature_y.currentText())
                feature_x_chan = int(self.ui.comboBox_feature_x_chan.currentText())-1
                feature_y_chan = int(self.ui.comboBox_feature_y_chan.currentText())-1
                print "Adding boundary on", feature_x, feature_x_chan, feature_y, feature_y_chan
                bound = (feature_x, feature_x_chan, feature_y, feature_y_chan, temp)
                self.activeClusterRadioButton().cluster_reference.addBound(bound)
                self.activeClusterRadioButton().cluster_reference.calculateMembership(self.spikeset)
                self.update_active_cluster()
                self.updateFeaturePlot()
            else:                
                self.prof_limits =self.prof_limits_reference
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
            # draw rectangle is in the qt coordinates, so we need a little conversion
            height = self.mp_proj.figure.bbox.height
            x0 = self.zoom_startpos[2]
            y0 = height - self.zoom_startpos[3]
            x1 = event.x
            y1 = height - event.y
            rect = [ int(val) for val in min(x0,x1), min(y0,y1), abs(x1-x0), abs(y1-y0)]
            self.mp_proj.drawRectangle(rect)

        pass
    
    def mp_proj_onPaint(self, event):
        MatplotlibWidget.paintEvent(self.mp_proj, event)
        if self.limit_mode:
            qp = QtGui.QPainter()
            height = self.mp_proj.figure.bbox.height
            
            qp.begin(self.mp_proj)
            qp.setPen(QtGui.QColor(*self.activeClusterRadioButton().cluster_reference.color))
            for i in range(len(self.limit_data)-1):
                qp.drawLine(self.limit_data[i][2], height - self.limit_data[i][3], self.limit_data[i+1][2], height - self.limit_data[i+1][3])
            if self.limit_data:
                qp.drawLine(self.limit_data[-1][2], height - self.limit_data[-1][3], self._mouse_move_event.x, height - self._mouse_move_event.y)
            qp.end()
            
    def save(self):
        if self.current_filename and self.clusters:
            (root, ext) = os.path.splitext(self.current_filename)
            
            # Save the bounds to a format we can easily read back in later
            outfilename = root + os.extsep + 'bounds'
            outfile = open(outfilename, 'wb')            
            save_bounds = [(cluster.color, cluster.bounds) for cluster in self.clusters ]
            pickle.dump(save_bounds, outfile)
            outfile.close()
            print "Saved bounds to", outfilename
            
            # Save the cluster membership vectors in matlab format for people to use
            outfilename = root + os.extsep + 'cluster'
            
            cluster_member = np.column_stack(tuple([cluster.member for cluster in self.clusters]))
            cluster_stats = [cluster.stats for cluster in self.clusters]            
            
            save_data = {'cluster_id':cluster_member, 'spike_time':self.spikeset.time}
            for key in cluster_stats[0].keys():
                save_data[key] = [stat[key] for stat in cluster_stats]
                        
            sio.savemat(outfilename, save_data, oned_as='column', appendmat=False)
            outfile.close()
            print "Saved cluster membership to", outfilename
            self.unsaved = False
            
    def importBounds(self, filename = None):
        if filename == None:
            filename = QtGui.QFileDialog.getOpenFileName(self, 'Open ntt file', filter='*.bounds')
            
        if filename:
            if os.path.exists(filename):
                print "Found bound file", filename, ", importing."
                infile = open(filename, 'rb')
                saved_bounds = pickle.load(infile)
                infile.close()
                print "Found", len(saved_bounds), "bounds to import, creating clusters."
                self.redrawing_details = True
                for (col, bound) in saved_bounds:
                    clust = self.add_cluster(col)
                    clust.bounds = bound
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
        if self.unsaved:
            reply = QtGui.QMessageBox.question(self, 'Save', "Do you want to save before quitting?", QtGui.QMessageBox.Yes | 
            QtGui.QMessageBox.No | QtGui.QMessageBox.Cancel, QtGui.QMessageBox.Yes)
            if reply == QtGui.QMessageBox.Cancel:
                event.ignore()
            if reply == QtGui.QMessageBox.Yes:
                event.accept()
                self.save()
            if reply == QtGui.QMessageBox.No:
                event.accept()
        
    
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = PyClustMainWindow()
    myapp.show()
    #myapp.load_ntt('TT2_neo.ntt')
    sys.exit(app.exec_())    
    
