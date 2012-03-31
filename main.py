# -*- coding: utf-8 -*-
"""
Created on Thu Mar 29 23:43:03 2012

@author: Bernard
"""

#!/usr/bin/python -d
 
 
 
import sys
from PyQt4 import QtCore, QtGui
from gui import Ui_MainWindow
import numpy as np
from spikeset import *
from features import *
import matplotlib as mpl
 
 
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

        self.redrawing_proj = False
        self.redrawing_details = False        
        
        # Connect the handlers        
        QtCore.QObject.connect(self.ui.pushButton, QtCore.SIGNAL("clicked()"),  self.button_load_click)
        #QtCore.QObject.connect(self.ui.lineEdit, QtCore.SIGNAL("returnPressed()"), self.add_entry)
        QtCore.QObject.connect(self.ui.pushButton_quit, QtCore.SIGNAL("clicked()"), self.tryQuit)
        QtCore.QObject.connect(self.ui.comboBox_feature_x_chan, QtCore.SIGNAL("currentIndexChanged(int)"), self.feature_channel_x_changed)
        QtCore.QObject.connect(self.ui.comboBox_feature_y_chan, QtCore.SIGNAL("currentIndexChanged(int)"), self.feature_channel_y_changed)
        QtCore.QObject.connect(self.ui.pushButton_next_projection, QtCore.SIGNAL("clicked()"), self.button_next_feature_click)
        QtCore.QObject.connect(self.ui.pushButton_previous_projection, QtCore.SIGNAL("clicked()"), self.button_prev_feature_click)
        QtCore.QObject.connect(self.ui.buttonGroup_scatter_plot_type, QtCore.SIGNAL("buttonClicked(int)"), self.updateFeaturePlot)
        QtCore.QObject.connect(self.ui.spinBox_markerSize, QtCore.SIGNAL("valueChanged(int)"), self.updateFeaturePlot)
        
        self.updateFeaturePlot()
                
        # Set up the cluster list area
        self.labels_container = QtGui.QWidget()
        self.ui.scrollArea_cluster_list.setWidget(self.labels_container)
        
        layout = QtGui.QVBoxLayout(self.labels_container)
        self.buttonGroup_cluster = QtGui.QButtonGroup()
        self.buttonGroup_cluster.setExclusive(True)        
        QtCore.QObject.connect(self.buttonGroup_cluster, QtCore.SIGNAL("buttonClicked(int)"), self.updateClusterDetailPlots)
        
        layout.addItem(QtGui.QSpacerItem(20, 40, QtGui.QSizePolicy.Minimum, QtGui.QSizePolicy.Expanding))        
        QtCore.QObject.connect(self.ui.pushButton_add_cluster, QtCore.SIGNAL("clicked()"), self.button_add_cluster_click)
        QtCore.QObject.connect(self.ui.pushButton_delete_cluster, QtCore.SIGNAL("clicked()"), self.delete_cluster)
        
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
        self.ui.checkBox_show_unclustered_exclusive.setChecked(False)
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
                
        self.ui.comboBox_feature_x.clear()
        self.ui.comboBox_feature_y.clear()
        self.ui.comboBox_feature_x_chan.clear()
        self.ui.comboBox_feature_y_chan.clear()
        
        self.activeClusterRadioButton = self.buttonGroup_cluster.checkedButton
        
        self.redrawing_proj = False
        
    def add_cluster(self):
        new_cluster = Cluster(self.spikeset)
                
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
        
        self.updateFeaturePlot()
        self.updateClusterDetailPlots()
        
    def delete_cluster(self):
        if self.activeClusterRadioButton():
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
        
        
        
    def button_cluster_color(self):
        color=QtGui.QColorDialog.getColor(QtGui.QColor(*self.sender().cluster_reference.color), self, "ColorDialog")
        self.sender().setStyleSheet("QPushButton {background-color: rgb(%d, %d, %d)}" % (color.red(),color.green(),color.blue()))
        self.sender().cluster_reference.color = (color.red(),color.green(),color.blue())
        
        self.updateFeaturePlot()
        
    def button_add_cluster_click(self):
        self.add_cluster()
        
    def feature_channel_x_changed(self, index):
        if index == -1: return
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
        self.prof_limits_reference == None
        self.redrawing_proj = False
        self.updateFeaturePlot()
                
    def feature_channel_y_changed(self, index):
        if index == -1: return
        if not self.redrawing_proj:
            self.prof_limits_reference == None
            self.updateFeaturePlot()
    
    def button_next_feature_click(self):
        if self.ui.comboBox_feature_y_chan.currentIndex() == self.ui.comboBox_feature_y_chan.count()-1:
            if not self.ui.comboBox_feature_x_chan.currentIndex() == self.ui.comboBox_feature_x_chan.count()-1:
                # step x to the next projection
                self.redrawing_proj = True
                self.ui.comboBox_feature_x_chan.setCurrentIndex(self.ui.comboBox_feature_x_chan.currentIndex()+1)
                # but this time we do want to loop around the y projection
                self.redrawing_proj = False
                self.ui.comboBox_feature_y_chan.setCurrentIndex(0)
        else:
            self.ui.comboBox_feature_y_chan.setCurrentIndex(self.ui.comboBox_feature_y_chan.currentIndex()+1)
    
    def button_prev_feature_click(self):
        if self.ui.comboBox_feature_y_chan.currentIndex() == 0:
            if not self.ui.comboBox_feature_x_chan.currentIndex() == 0:
                # step x to the previous projection
                self.redrawing_proj = True
                self.ui.comboBox_feature_x_chan.setCurrentIndex(self.ui.comboBox_feature_x_chan.currentIndex()-1)
                # but this time we do want to loop around the y projection
                self.redrawing_proj = False
                self.ui.comboBox_feature_y_chan.setCurrentIndex(self.ui.comboBox_feature_y_chan.count()-1)
        else:
            self.ui.comboBox_feature_y_chan.setCurrentIndex(self.ui.comboBox_feature_y_chan.currentIndex()-1)

    def load_ntt(self, fname): 
        print 'Loading ntt file', fname
        self.spikeset = loadNtt(fname)
        print 'Loaded', self.spikeset.N, 'spikes'
        self.current_feature = Feature_Peak(self.spikeset)

        # Set the combo boxes to the current feature for now
        self.ui.comboBox_feature_x.clear();
        self.ui.comboBox_feature_x.addItem(self.current_feature.name)
        self.ui.comboBox_feature_y.clear();
        self.ui.comboBox_feature_y.addItem(self.current_feature.name)

        # Set the channel combo box options
        self.ui.comboBox_feature_x_chan.clear();
        self.ui.comboBox_feature_y_chan.clear();
        print "Defaulting to channels 1,2"
        for i in range(1, self.current_feature.channels):
            self.ui.comboBox_feature_x_chan.addItem(str(i))        

        self.updateFeaturePlot()
        self.updateClusterDetailPlots()
        
    def button_load_click(self):
        fname = QtGui.QFileDialog.getOpenFileName(self, 'Open ntt file', filter='*.ntt')
        self.load_ntt(fname)
        
    def updateFeaturePlot(self):
        if self.current_feature == None: return
        if self.redrawing_proj: return

        proj_x = int(self.ui.comboBox_feature_x_chan.currentText())-1
        proj_y = int(self.ui.comboBox_feature_y_chan.currentText())-1

        print "Plotting features"        
        # Draw projection density map
        self.mp_proj.figure.clear()
        self.mp_proj.figure.subplots_adjust(hspace=0.0001, wspace=0.0001, bottom=0.075, top=1, left=0.075, right=1)
        self.mp_proj.axes = self.mp_proj.figure.add_subplot(1,1,1)
        
        if self.prof_limits_reference == None:
            self.prof_limits_reference = ([np.min(self.current_feature.data[:,proj_x]), np.max(self.current_feature.data[:,proj_x])*1.05], [np.min(self.current_feature.data[:,proj_y]), np.max(self.current_feature.data[:,proj_y])*1.05])
            self.prof_limits = self.prof_limits_reference
                            
        if self.ui.radioButton_scatter.isChecked():
            # Iterate over clusters
            self.ui.mplwidget_projection.axes.hold(False)
            if self.ui.checkBox_show_unclustered.isChecked():
                w = np.array([True] * self.spikeset.N)
                if self.ui.checkBox_show_unclustered_exclusive.isChecked():
                    for cluster in self.clusters:
                        w[cluster.member] = False
                self.ui.mplwidget_projection.axes.scatter(self.current_feature.data[w,proj_x], self.current_feature.data[w,proj_y], marker='.', s=int(self.ui.spinBox_markerSize.text()), color='k')
            self.ui.mplwidget_projection.axes.hold(True)
            for cluster in self.clusters:
                if not cluster.check.isChecked(): continue
                self.mp_proj.axes.scatter(self.current_feature.data[cluster.member,proj_x], self.current_feature.data[cluster.member,proj_y], marker='.', s=int(self.ui.spinBox_markerSize.text()), color=map(lambda s: s / 255.0, cluster.color))
                
        else:
            bins_x = np.linspace(self.prof_limits[0][0], self.prof_limits[0][1], 100)
            print bins_x
            bins_y = np.linspace(self.prof_limits[1][0], self.prof_limits[1][1], 100)
            print bins_y
            count = np.histogram2d(self.current_feature.data[:,proj_x], self.current_feature.data[:,proj_y], [bins_x, bins_y])[0]
            if self.ui.radioButton_density.isChecked():
                self.mp_proj.axes.pcolor(bins_x, bins_y, np.transpose(count))
            else:
                self.mp_proj.axes.pcolor(bins_x, bins_y, np.transpose(np.log(count+1)))
                
        self.mp_proj.axes.set_xlim(self.prof_limits[0])
        self.mp_proj.axes.set_ylim(self.prof_limits[1])
                    
        self.mp_proj.draw()
        self.redrawing_proj = False
        
    def updateClusterDetailPlots(self):               
        if self.current_feature == None: return
        
        if self.activeClusterRadioButton():
            w = self.activeClusterRadioButton().cluster_reference.member
        else:
            w = np.array([True] * self.spikeset.N)
                
        # Draw average waveforms
        self.mp_wave.figure.clear()
        self.mp_wave.figure.subplots_adjust(hspace=0.0001, wspace=0.0001, bottom=0, top=1, left=0.0, right=1)                
        for i in range(0,4):
            ax = self.mp_wave.figure.add_subplot(2,2,i+1)
            ax.errorbar(range(1,33), np.mean(self.spikeset.spikes[w,i,:], axis=0), np.std(self.spikeset.spikes[w,i,:], axis=0), color='k')
            ax.set_xticks([])
            ax.set_yticks([])
        self.mp_wave.draw()

        time = self.spikeset.time[w]
        isi = (time[1:] - time[0:-1])/1e3
        isi_bins = np.logspace(np.log10(1), np.log10(2500), 100)
        isi_bin_centers = (isi_bins[0:-1] + isi_bins[1:])
        isi_bin_count = np.histogram(isi, isi_bins)[0]

        # Draw ISI histogram
        self.mp_isi.figure.clear()
        self.mp_isi.figure.subplots_adjust(hspace=0.0001, wspace=0.0001, bottom=0.1, top=1, left=0.0, right=1)                
        self.mp_isi.axes = self.mp_isi.figure.add_subplot(1,1,1)
        
        self.mp_isi.axes.plot(isi_bin_centers, isi_bin_count)
        self.mp_isi.axes.set_xscale('log')
        self.mp_isi.axes.set_xlim([1, 2500])
        detect_line = mpl.lines.Line2D([2, 2], self.mp_isi.axes.get_ylim(), color='r', linestyle='--')
        self.mp_isi.axes.add_line(detect_line)
        self.mp_isi.axes.set_xticks([2, 10, 100])
        self.mp_isi.axes.set_xticklabels(['2', '10', '100'])
        self.mp_isi.draw()
        
        
        
    def onMousePress(self, event):
        print event.xdata, event.ydata
        if not (event.xdata == None or event.ydata == None) and event.button == 1:
            self.zoom_active = True
            self.zoom_startpos = (event.xdata, event.ydata, event.x, event.y)
        #print 'wee', event
        #print 'zomg', event.xdata, event.ydata
        
    def onMouseRelease(self, event):
        #print event.button, event.xdata, event.ydata
        if event.button == 1 and self.zoom_active and not (event.xdata == None or event.ydata == None):
            self.zoom_active = False
            self.prof_limits = (np.sort([event.xdata, self.zoom_startpos[0]]), np.sort([event.ydata, self.zoom_startpos[1]]))
            #self.mp_proj.axes.set_xlim(np.sort()
            #self.mp_proj.axes.set_ylim(np.sort([event.ydata, self.zoom_startpos[1]]))
            
            if self.ui.radioButton_scatter.isChecked():
                self.mp_proj.axes.set_xlim(self.prof_limits[0])
                self.mp_proj.axes.set_ylim(self.prof_limits[1])
                self.mp_proj.draw()
            else: # if its a density plot we should rebin for now
                self.updateFeaturePlot()
        # reset bounds on right click    
        if event.button == 3:
            self.prof_limits =self.prof_limits_reference
            #self.mp_proj.axes.set_xlim(self.prof_limits_reference[0])
            #self.mp_proj.axes.set_ylim(self.prof_limits_reference[1])
            #self.mp_proj.draw()
            if self.ui.radioButton_scatter.isChecked():
                self.mp_proj.axes.set_xlim(self.prof_limits[0])
                self.mp_proj.axes.set_ylim(self.prof_limits[1])
                self.mp_proj.draw()
            else:
                self.updateFeaturePlot()
        
    def onMouseMove(self, event):
        if self.zoom_active:
            # draw rectangle is in the qt coordinates, so we need a little conversion
            height = self.mp_proj.figure.bbox.height
            x0 = self.zoom_startpos[2]
            y0 = height - self.zoom_startpos[3]
            x1 = event.x
            y1 = height - event.y
            rect = [ int(val) for val in min(x0,x1), min(y0,y1), abs(x1-x0), abs(y1-y0)]
            self.mp_proj.drawRectangle(rect)

        pass
        
    def tryQuit(self):
        app.exit()
     
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    myapp = PyClustMainWindow()
    myapp.show()
    myapp.load_ntt('Sample3.ntt')
    sys.exit(app.exec_())    
    
