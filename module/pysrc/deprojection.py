import scipy.interpolate as interpolate
import numpy as np

def deproject_image(image, ratio, pa_rad, nx0, ny0):
    nx, ny = image.shape
    # create coordinate system
    x = np.arange(nx)
    y = np.arange(ny)
    xx, yy = np.meshgrid(x, y, indexing="xy")
    # deproject
    dxx1 =  (xx-nx0)*np.cos(pa_rad) + (yy-ny0)*np.sin(pa_rad)
    dyy1 = -(xx-nx0)*np.sin(pa_rad) + (yy-ny0)*np.cos(pa_rad)
    dyy1 = dyy1 / ratio
    xx2 =  dxx1*np.cos(pa_rad) - dyy1*np.sin(pa_rad) + nx0
    yy2 = +dxx1*np.sin(pa_rad) + dyy1*np.cos(pa_rad) + ny0

    # interpolate image
    f = interpolate.RectBivariateSpline(np.arange(nx), np.arange(ny), image.T)
    image_deprojected = f.ev(xx2, yy2)
    return image_deprojected

def deproject_grid(xx, yy, ratio, pa_rad, x0=0, y0=0):
    # deproject
    dxx1 =  (xx-x0)*np.cos(pa_rad) + (yy-y0)*np.sin(pa_rad)
    dyy1 = -(xx-x0)*np.sin(pa_rad) + (yy-y0)*np.cos(pa_rad)
    dyy1 = dyy1 / ratio
    xx2 =  dxx1*np.cos(pa_rad) - dyy1*np.sin(pa_rad) + x0
    yy2 = +dxx1*np.sin(pa_rad) + dyy1*np.cos(pa_rad) + y0
    return xx2, yy2