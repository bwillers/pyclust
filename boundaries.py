
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
        self.special_info = None # might need things like PCA coefficients

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

# Helper functions
def ellipseFromCovariance(covar):
    vals, vecs = np.linalg.eigh(covar)
    ind = np.argsort(vals)[::-1]  # descending sort
    vals = vals[ind]
    vecs = vecs[ind,:]
    projv = vecs[0] / np.linalg.norm(vecs[0])
    angle = np.arctan2(projv[1], projv[0])
    radiusx = np.sqrt(vals[0])
    radiusy = np.sqrt(vals[1])
    return(angle, radiusx, radiusy)

def covarianceFromEllipse(angle, radiusx, radiusy):
    # eigen vectors as columns
    u = np.array([[np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]])
    d = np.array([[radiusx ** 2, 0], [0, radiusy ** 2]])
    return np.dot(np.dot(u, d), np.linalg.inv(u))

# 2D Elliptical boundary
class BoundaryEllipse2D(Boundary):

    # size should be (rwidth x rheight) in rotated space
    def __init__(self, featureNames, featureChans, center, angle, size):
        self.features = featureNames
        self.chans = featureChans
        self.center = center
        self.angle = angle
        self.size = size
        self.special_info = None

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

if __name__ == '__main__':
    a = np.pi/5
    rx = 3.0
    ry = 1.0
    print a * 180 / np.pi, rx, ry
    cov = covarianceFromEllipse(a, rx, ry)
    print cov
    angle, rx, ry = ellipseFromCovariance(cov)
    print angle * 180 / np.pi, rx, ry
    print np.mod(np.abs(angle - a), 2 * np.pi) * 180 / np.pi

