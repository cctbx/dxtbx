import math
import numpy as np

def rotate_and_average(data, angle, deg=False, mask=None):
  """
  Angle is the rotation angle in the coordinate system of x along fast, y along slow, such that angle increases as it rotates around z, where z completes a right-handed system.
  """
  if deg:
    angle *= math.pi/180
  angle = -angle

  ny, nx = np.shape(data)

  xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
  xx_yy = np.row_stack((xx.ravel(), yy.ravel()))
  R = np.array((
      (np.cos(angle), -np.sin(angle)),
      (np.sin(angle), np.cos(angle))
      ))
  xx_yy_rotated = np.matmul(R, xx_yy)
  xx_rotated = xx_yy_rotated[0, :].reshape((ny, nx))
  yy_rotated = xx_yy_rotated[1, :].reshape((ny, nx))
  bin_x = np.arange(-nx, nx)
  bin_x_centers = (bin_x[1:] + bin_x[:-1]) / 2

  if mask is None:
      mask = np.full(data.shape, True)

  counts_rotated, _ = np.histogram(xx_rotated[mask], bins=bin_x)
  summed_rotated, _ = np.histogram(xx_rotated[mask], bins=bin_x, weights=data[mask])
  sel = counts_rotated > 0
  average_rotated = summed_rotated[sel] / counts_rotated[sel]

  return bin_x_centers[sel], average_rotated
