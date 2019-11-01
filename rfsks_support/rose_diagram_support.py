
import numpy as np
from matplotlib.projections import PolarAxes
from matplotlib.transforms import Affine2D, Bbox, IdentityTransform

class NorthPolarAxes(PolarAxes):

name = ‘northpolar’

def format_coord(self, angle, radius):
angle = np.pi * 0.5 – angle
if angle < 0:
angle = angle + 2 * np.pi
angle = angle * 180 / np.pi
return u’\u03b8=%f\u00b0, r=%f’ % (angle, radius)

class NorthPolarTransform(PolarAxes.PolarTransform):
def transform(self, tr):
xy = np.zeros(tr.shape, np.float_)
t = tr[:, 0:1]
r = tr[:, 1:2]
x = xy[:, 0:1]
y = xy[:, 1:2]
x[:] = r * np.sin(t)
y[:] = r * np.cos(t)
return xy

transform_non_affine = transform

def inverted(self):
return NorthPolarAxes.InvertedNorthPolarTransform()

class InvertedNorthPolarTransform(PolarAxes.InvertedPolarTransform):
def transform(self, xy):
x = xy[:, 0:1]
y = xy[:, 1:]
r = np.sqrt(x*x + y*y)
theta = np.arctan2(y, x)
return np.concatenate((theta, r), 1)

def inverted(self):
return NorthPolarAxes.NorthPolarTransform()

def _set_lim_and_transforms(self):
PolarAxes._set_lim_and_transforms(self)
self.transProjection = self.NorthPolarTransform()
self.transData = (
self.transScale +
self.transProjection +
(self.transProjectionAffine + self.transAxes))
self._xaxis_transform = (
self.transProjection +
self.PolarAffine(IdentityTransform(), Bbox.unit()) +
self.transAxes)
self._xaxis_text1_transform = (
self._theta_label1_position +
self._xaxis_transform)
self._yaxis_transform = (
Affine2D().scale(np.pi * 2.0, 1.0) +
self.transData)
self._yaxis_text1_transform = (
self._r_label1_position +
Affine2D().scale(1.0 / 360.0, 1.0) +
self._yaxis_transform)