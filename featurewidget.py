
__version__ = "1.0.0"

from PyQt4.QtGui import *
from PyQt4.QtCore import *

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as Canvas
from matplotlib.figure import Figure

import numpy as np
import matplotlib as mpl

from matplotlib import rcParams
rcParams['font.size'] = 9


class ProjectionWidget(Canvas):

    __pyqtSignals__ = ("featureRedrawRequired()",
    "polygonBoundaryDrawn(PyObject*)")

    def __init__(self, parent=None):
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.axes.set_title('')
        self.axes.set_xlabel('')
        self.axes.set_ylabel('')
        Canvas.__init__(self, self.figure)
        self.setParent(parent)

        pal = self.palette().window().color()
        bgcolor = (pal.red() / 255.0, pal.blue() / 255.0, pal.green() / 255.0)

        self.figure.clear()
        self.figure.set_facecolor(bgcolor)
        # Set up the feature axes
        self.figure.subplots_adjust(hspace=0.0001, wspace=0.0001,
            bottom=0.0, top=1, left=0.0, right=1)
        self.axes = self.figure.add_subplot(1, 1, 1)

        Canvas.setSizePolicy(self, QSizePolicy.Expanding,
                QSizePolicy.Expanding)
        Canvas.updateGeometry(self)

        # Properties describing how to draw the clusters
        self.unclustered = True
        self.refractory = True
        self.exclusive = True
        self.markersize = 1
        self.ptype = -2

        self.feature_x = 'Peak'
        self.feature_y = 'Peak'
        self.chan_x = 0
        self.chan_y = 1

        # Mouse handlers
        self.zoom_active = False
        self.limit_mode = False
        self.limit_cluster = None
        self.mpl_connect('button_press_event', self.onMousePress)
        self.mpl_connect('button_release_event', self.onMouseRelease)
        self.mpl_connect('motion_notify_event', self.onMouseMove)

    @pyqtSignature("setShowUnclustered(bool)")
    def setShowUnclustered(self, show):
        self.unclustered = show
        self.emit(SIGNAL("featureRedrawRequired()"))

    @pyqtSignature("setUnclusteredExclusive(bool)")
    def setUnclusteredExclusive(self, excl):
        self.exclusive = excl
        self.emit(SIGNAL("featureRedrawRequired()"))

    @pyqtSignature("setRefractory(bool)")
    def setRefractory(self, show):
        self.refractory = show
        self.emit(SIGNAL("featureRedrawRequired()"))

    @pyqtSignature("setMarkerSize(int)")
    def setMarkerSize(self, size):
        self.markersize = size
        self.emit(SIGNAL("featureRedrawRequired()"))

    @pyqtSignature("setPlotType(int)")
    def setPlotType(self, ptype):
        self.ptype = ptype
        self.emit(SIGNAL("featureRedrawRequired()"))

    def resetLimits(self):
        self.prof_limits_reference = None
#        self.emit(SIGNAL("featureRedrawRequired()"))

    # These are not provided as signals since there is some non trivial
    # logic the main form needs to perform first
    def setFeatureX(self, name):
        self.feature_x = name
        self.resetLimits()
        self.stopMouseAction()
        self.emit(SIGNAL("featureRedrawRequired()"))

    def setFeatureY(self, name):
        self.feature_y = name
        self.resetLimits()
        self.stopMouseAction()
        self.emit(SIGNAL("featureRedrawRequired()"))

    def setChanX(self, chan):
        self.chan_x = chan
        self.resetLimits()
        self.stopMouseAction()
        self.emit(SIGNAL("featureRedrawRequired()"))

    def setChanY(self, chan):
        self.chan_y = chan
        self.resetLimits()
        self.stopMouseAction()
        self.emit(SIGNAL("featureRedrawRequired()"))

    # Given a set of spike and cluster data, create the figure
    def updatePlot(self, spikeset, clusters, junk):
        xdata = spikeset.featureByName(self.feature_x).data[:, self.chan_x]
        ydata = spikeset.featureByName(self.feature_y).data[:, self.chan_y]
        self.axes.hold(False)

        #print "Updating feature plot: ", self.feature_x, self.chan_x, \
        #        " vs ", self.feature_y, self.chan_y

        if self.prof_limits_reference == None:
            w = np.array([True] * spikeset.N)
            if not junk.check.isChecked():
                w[junk.member] = False

            temp = ([np.min(xdata[w]), np.max(xdata[w])],
                    [np.min(ydata[w]), np.max(ydata[w])])

            w_x = (temp[0][1] - temp[0][0]) * 0.05
            w_y = (temp[1][1] - temp[1][0]) * 0.05
            self.prof_limits_reference = ([temp[0][0] - w_x, temp[0][1] + w_x],
                [temp[1][0] - w_y, temp[1][1] + w_y])

            self.prof_limits = self.prof_limits_reference

        # Scatter Plot
        if self.ptype == -2:
            #Plot the unclustered spikes
            if self.unclustered:
                w = np.array([True] * spikeset.N)

                if self.exclusive:
                    for cluster in clusters + [junk]:
                        w[cluster.member] = False

                self.axes.plot(xdata[w], ydata[w], linestyle='None',
                    marker='o', markersize=self.markersize,
                    markerfacecolor='k',  markeredgecolor='k')

                self.axes.hold(True)

            # Iterate over clusters
            for cluster in clusters + [junk]:
                if not cluster.check.isChecked():
                    continue

                col = map(lambda s: s / 255.0, cluster.color)
                # Plot the cluster spikes
                self.axes.plot(xdata[cluster.member], ydata[cluster.member],
                        marker='o', markersize=self.markersize,
                        markerfacecolor=col, markeredgecolor=col,
                        linestyle='None')
                self.axes.hold(True)

                # Plot refractory spikes
                if self.refractory:
                    self.axes.plot(xdata[cluster.refractory],
                        ydata[cluster.refractory], marker='o', markersize=5,
                        markerfacecolor='k', markeredgecolor='k',
                        linestyle='None')

        # Do a density plot
        else:
            w = np.array([False] * spikeset.N)
            if self.unclustered:
                w = np.array([True] * spikeset.N)
                if self.exclusive:
                        for cluster in clusters + [junk]:
                                w[cluster.member] = False
            for cluster in \
                [cluster for cluster in \
                    clusters + [junk] if cluster.check.isChecked()]:
                    w[cluster.member] = True

            if not np.any(w):
                w[0] = True  # Histogram routine fails with empty data

            bins_x = np.linspace(self.prof_limits[0][0],
                self.prof_limits[0][1], 100)
            bins_y = np.linspace(self.prof_limits[1][0],
                self.prof_limits[1][1], 100)
            count = np.histogram2d(xdata[w],
                ydata[w], [bins_x, bins_y])[0]

            if self.ptype == -3:
                self.axes.pcolor(bins_x, bins_y, np.transpose(count))
            else:
                self.axes.pcolor(bins_x, bins_y,
                    np.transpose(np.log(count + 1)))

            # Iterate over clusters for refractory spikes
            if self.refractory:
                self.axes.hold(True)
                for cluster in clusters + [junk]:
                    if not cluster.check.isChecked():
                        continue

                    # Plot refractory spikes
                    self.axes.plot(xdata[cluster.refractory],
                        ydata[cluster.refractory], marker='o', markersize=5,
                        markerfacecolor='w', markeredgecolor='w',
                        linestyle='None')

        self.axes.set_xlim(self.prof_limits[0])
        self.axes.set_ylim(self.prof_limits[1])
        self.axes.set_xticks([])
        self.axes.set_yticks([])

        # Now draw the boundaries
        for cluster in clusters + [junk]:
            if not cluster.check.isChecked():
                continue

            # Limit boundaries with solid line
            bound_polys = cluster.getBoundPolygon(self.feature_x, self.chan_x,
                self.feature_y, self.chan_y)

            if bound_polys:
                for bound_poly in bound_polys:
                    bound_poly = np.vstack((bound_poly, bound_poly[0, :]))
                    col = map(lambda s: s / 255.0, cluster.color)
                    for i in range(np.size(bound_poly, 0)):
                        line = mpl.lines.Line2D(bound_poly[i:i + 2, 0],
                            bound_poly[i:i + 2, 1], color=col, linestyle='-')
                        self.axes.add_line(line)

            # Addition boundaries with dashed line
            bound_polys = cluster.getBoundPolygon(self.feature_x, self.chan_x,
                self.feature_y, self.chan_y, 'add')

            if bound_polys:
                for bound_poly in bound_polys:
                    bound_poly = np.vstack((bound_poly, bound_poly[0, :]))
                    col = map(lambda s: s / 255.0, cluster.color)

                    for i in range(np.size(bound_poly, 0)):
                        line = mpl.lines.Line2D(bound_poly[i:i + 2, 0],
                            bound_poly[i:i + 2, 1], color=col, linestyle='--')
                        self.axes.add_line(line)

    def onMousePress(self, event):
        if self.limit_mode:
            return
        if (event.xdata != None and event.ydata != None) and event.button == 1:
            self.zoom_active = True
            self.zoom_startpos = (event.xdata, event.ydata, event.x, event.y)

    def onMouseRelease(self, event):
        # Left mouse button
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

                if self.ptype == -2:
                    self.axes.set_xlim(self.prof_limits[0])
                    self.axes.set_ylim(self.prof_limits[1])
                    self.draw()
                else:  # if its a density plot we should rebin for now
                    self.emit(SIGNAL("featureRedrawRequired()"))
            if self.limit_mode:
                self.limit_data.append((event.xdata, event.ydata,
                    event.x, event.y))

        # reset bounds on right click, or cancel / complete the bounds
        if event.button == 3 and self.limit_mode:
            if event.xdata != None and event.ydata != None:
                self.limit_data.append((event.xdata, event.ydata,
                    event.x, event.y))
                self.limit_mode = False

                # create an array to describe the bound
                temp = np.array([np.array([point[0], point[1]]) for
                    point in self.limit_data])

                print "Adding boundary on", self.feature_x, self.chan_x + 1, \
                    self.feature_y, self.chan_y + 1

                bound = (self.feature_x, self.chan_x,
                    self.feature_y, self.chan_y, temp)

                self.emit(SIGNAL("polygonBoundaryDrawn(PyQt_PyObject)"), bound)
                self.stopMouseAction()

        # Right click resets to original bounds
        elif event.button == 3:
            self.prof_limits = self.prof_limits_reference
            if self.ptype == -2:
                self.axes.set_xlim(self.prof_limits[0])
                self.axes.set_ylim(self.prof_limits[1])
                self.draw()
            else:
                self.emit(SIGNAL("featureRedrawRequired()"))

    def onMouseMove(self, event):
        if self.limit_mode:
            self._mouse_move_event = event
            self.repaint()
        if (not self.limit_mode) and self.zoom_active:
            # draw rectangle is in the qt coordinates,
            # so we need a little conversion
            height = self.figure.bbox.height
            x0 = self.zoom_startpos[2]
            y0 = height - self.zoom_startpos[3]
            x1 = event.x
            y1 = height - event.y
            rect = [int(val) for val in min(x0, x1),
                min(y0, y1), abs(x1 - x0), abs(y1 - y0)]
            self.drawRectangle(rect)

    # Paint using standard Canvas and then draw lines during
    # boundary creation as required
    def paintEvent(self, event):
        Canvas.paintEvent(self, event)
        if self.limit_mode:
            qp = QPainter()
            height = self.figure.bbox.height

            qp.begin(self)
            qp.setPen(QColor(*self.limit_color))
            for i in range(len(self.limit_data) - 1):
                qp.drawLine(self.limit_data[i][2],
                    height - self.limit_data[i][3], self.limit_data[i + 1][2],
                    height - self.limit_data[i + 1][3])
            if self.limit_data:
                qp.drawLine(self.limit_data[-1][2],
                    height - self.limit_data[-1][3], self._mouse_move_event.x,
                    height - self._mouse_move_event.y)
            qp.end()

    # Start drawing a boundary for a cluster
    def drawBoundary(self, color):
        self.limit_color = color
        self.limit_mode = True
        self.limit_data = []

    def stopMouseAction(self):
        self.limit_mode = False
        self.limit_color = None
        self.zoom_active = False
