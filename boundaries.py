
import numpy as np
import matplotlib.nxutils as nx
import matplotlib as mpl

# A simple boundary base class
class Boundary:

    def withinBoundary(self, spikeset):
        return np.array([True] * spikeset.N)

    def draw(self, axes, color, linestyle):
        return

# 2D Polygon boundary
class BoundaryPolygon2D(Boundary):

    # featureNames should be a length 2 tuple (featureX, featureY)
    # bounds is a mx2 set of polygon coordinates
    def __init__(self, featureNames, featureChans, bounds):
        self.features = featureNames
        self.chans = featureChans
        self.bounds = bounds

    def withinBoundary(self, spikeset):
        px = spikeset.featureByName(self.features[0]).data[:, self.chans[0]]
        py = spikeset.featureByName(self.features[1]).data[:, self.chans[1]]

        data = np.column_stack((px, py))

        return nx.points_inside_poly(data, self.bounds)

    def draw(self, axes, color='k', linestyle='-'):
        bound = np.vstack((self.bounds, self.bounds[0, :]))
        for i in range(np.size(bound, 0)):
            line = mpl.lines.Line2D(bound[i:i + 2, 0],
                bound[i:i + 2, 1], color=color, linestyle=linestyle)
            axes.add_line(line)


# 2D Elliptical boundary
class BoundaryEllipse2D(Boundary):

    # size should be (rwidth x rheight) in rotated space
    def __init__(self, featureNames, featureChans, center, angle, size):
        self.features = featureNames
        self.chans = featureChans
        self.center = center
        self.angle = angle
        self.size = size

    def withinBoundary(self, spikeset):
        px = spikeset.featureByName(self.features[0]).data[:, self.chans[0]]
        py = spikeset.featureByName(self.features[1]).data[:, self.chans[1]]

        data = np.column_stack((px, py)) - self.center

        rotmat = np.array([[ np.cos(self.angle), np.sin(self.angle)],
            [-np.sin(self.angle), np.cos(self.angle)]])

        data = np.dot(rotmat, data.T).T

        data = data / self.size
        return np.sum(np.power(data, 2), axis=1) <= 1

    def draw(self, axes, color='k', linestyle='-'):
        if linestyle == '-':
            linestyle = 'solid'
        elif linestyle == '--':
            linestyle = 'dashed'
        else:
            linestyle = 'dashdot'

        ell = mpl.patches.Ellipse(self.center, self.size[0] * 2.0,
                self.size[1] * 2.0, 180 * self.angle / np.pi, color=color,
                linestyle=linestyle,  fill=False)
        axes.add_artist(ell)
