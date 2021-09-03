import pygmt
from xarray import DataArray
import numpy as np
from Strain_2D.Strain_Tools.strain.pygmt_plots import station_vels_to_arrays
from Strain_2D.Strain_Tools.strain.output_manager import get_grid_eigenvectors

def plot_principal_strains(x, y, e, n, color, linewidth = 0.25, pivot = 'tail', scale = 5E-2):
    import matplotlib.pylab as plt
    plt.quiver(x, y, e*scale, n*scale, scale = 1, color=color, linewidth = linewidth, edgecolor = 'k', pivot = pivot)
    plt.quiver(x, y, -e*scale, -n*scale, scale = 1, color=color, linewidth = linewidth, edgecolor = 'k', pivot = pivot)

def plot_dilatation(lons, lats, exx, exy, eyy, region, value_range = None, width = 2000, do_not_print_value = 2.0, overmax_scale = 2.0, eigs_dec = 1, hillshade = 'data/MammothHillshade.tif', extent = [-119.2, -118.6, 37.47, 37.82], scale = 5E-2, colorbar = True):
    from Strain_2D.Strain_Tools.strain.strain_tensor_toolbox import compute_eigenvectors

    import matplotlib.pylab as plt
    [e1, e2, v00, v01, v10, v11] = compute_eigenvectors(exx, exy, eyy)
    [positive_eigs, negative_eigs] = get_grid_eigenvectors(lons, lats, e1, e2, v00, v01, v10, v11, do_not_print_value=do_not_print_value, overmax_scale=overmax_scale, eigs_dec=eigs_dec)
    dshape = np.shape(exx)
    dilatation = np.zeros(np.shape(exx))
    for j in range(dshape[0]):
        for i in range(dshape[1]):
            dilatation[j][i] = e1[j][i] + e2[j][i]
    if value_range is None:
        value_range = [np.nanmin(dilatation), np.nanmax(dilatation)]
    if hillshade is not None:
        imhs = plt.imread(hillshade)
        plt.imshow(imhs, extent = extent, cmap = 'gist_yarg')
    plt.imshow(np.flipud(dilatation), extent = [np.min(lons), np.max(lons), np.min(lats), np.max(lats)], cmap = 'RdBu', vmin = -overmax_scale, vmax = overmax_scale, alpha = 0.5 if hillshade is not None else 1.0)
    if colorbar:
        plt.colorbar()
    elon, nlat, e, n = station_vels_to_arrays(positive_eigs)
    plot_principal_strains(elon, nlat, e, n, color='b', scale = scale)
    elon, nlat, e, n = station_vels_to_arrays(negative_eigs)
    plot_principal_strains(elon, nlat, e, n, color='r', pivot='tip', scale = scale)
    plt.plot(-119.021559, 37.613324, 'w+')
    plt.text(-119.021559+0.001, 37.613324+0.001, 'HSL', color='w')
    plt.gca().set_aspect('equal', 'box')
    plt.axis([np.min(lons), np.max(lons), np.min(lats), np.max(lats)])

def plot_uplift(lons, lats, uzs, value_range = None, hillshade = 'data/MammothHillshade.tif', extent = [-119.2, -118.6, 37.47, 37.82]):
    import matplotlib.pylab as plt
    from matplotlib.colors import Normalize
    if value_range is None:
        value_range = np.array([np.nanmin(np.array(uzs)), np.nanmax(np.array(uzs))])
    else:
        value_range = np.array(value_range)
    imhs = plt.imread(hillshade)
    plt.imshow(imhs, extent = extent, cmap = 'gist_yarg')
    plt.scatter(lons, lats, c = uzs, norm=Normalize(vmin = value_range[0], vmax = value_range[1]), cmap = 'RdBu_r')
    plt.colorbar()
    plt.plot(-119.021559, 37.613324, 'w+')
    plt.text(-119.021559+0.001, 37.613324+0.001, 'HSL', color='w')
    plt.gca().set_aspect('equal', 'box')
    plt.axis([extent[0], extent[1], extent[2], extent[3]])


