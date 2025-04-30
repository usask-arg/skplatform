from typing import Tuple, Union
from numpy.typing import NDArray

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Rectangle, FancyArrow
import mpl_toolkits.mplot3d.art3d as art3d


def make_plot(rows: int, cols: int) -> Union[Tuple[Figure, NDArray[Axes]], Tuple[Figure, Axes]]:
    """
    Makes a basic grid of axes configured for 3D plotting.
    
    Parameters
    ----------
    rows : int 
        The number of rows in the grid. Must be a value of 1 or greater.
    cols : int
        The number of columns in the grid. Must be a value of 1 or greater.
        
    Returns
    -------
    Figure
        The figure object for the plot created.
    ax : NDArray[Axes]
        An array of axes objects. If a single plot is requested a bare axes object is returned.
    """
    return plt.subplots(rows, cols, figsize=(12.5, 12.5), dpi=150, subplot_kw={'projection': '3d'})


def plot_basic_cube(ax: Axes):
    """
    Adds rectangles on the xy (green), xz (red) and yz (blue) planes. Each rectangle is a 1 by 1
    patch object with an alpha of 0.2, that extends from the origin to 1, 1. 
    """
    xy_plane = Rectangle([0, 0], 1, 1, color='green', alpha=0.2)
    ax.add_patch(xy_plane)
    art3d.pathpatch_2d_to_3d(xy_plane, z=0, zdir="z")

    xz_plane = Rectangle([0, 0], 1, 1, color='blue', alpha=0.2)
    ax.add_patch(xz_plane)
    art3d.pathpatch_2d_to_3d(xz_plane, z=1, zdir="y")

    yz_plane = Rectangle([0, 0], 1, 1, color='red', alpha=0.2)
    ax.add_patch(yz_plane)
    art3d.pathpatch_2d_to_3d(yz_plane, z=0, zdir="x")


def plot_axes_planes(ax: Axes):
    """
    Adds rectangle patches for the xy (green) and yz (red) planes.
    The pathes are centered on the origin and extend from -0.5, -0.5 
    to 0.5, 0.5.
    """
    xy_plane = Rectangle([-1, -1], 2, 2, color='green', alpha=0.2)
    ax.add_patch(xy_plane)
    art3d.pathpatch_2d_to_3d(xy_plane, z=0, zdir="z")

    # xz_plane = Rectangle([0, 0], 1, 1, color='blue', alpha=0.2)
    # ax.add_patch(xz_plane)
    # art3d.pathpatch_2d_to_3d(xz_plane, z=1, zdir="y")

    yz_plane = Rectangle([-1, -1], 2, 2, color='red', alpha=0.2)
    ax.add_patch(yz_plane)
    art3d.pathpatch_2d_to_3d(yz_plane, z=0, zdir="x")


def plot_box_planes(ax: Axes):
    """
    Adds rectangles on the xy (green), xz (red) and yz (blue) planes. Each rectangle is a 1 by 1
    patch object with an alpha of 0.2, that extends from the -1, -1 to 1, 1. 
    """
    xy_plane = Rectangle([-1, -1], 2, 2, color='green', alpha=0.2)
    ax.add_patch(xy_plane)
    art3d.pathpatch_2d_to_3d(xy_plane, z=-1, zdir="z")

    xz_plane = Rectangle([-1, -1], 2, 2, color='blue', alpha=0.2)
    ax.add_patch(xz_plane)
    art3d.pathpatch_2d_to_3d(xz_plane, z=1, zdir="y")

    yz_plane = Rectangle([-1, -1], 2, 2, color='red', alpha=0.2)
    ax.add_patch(yz_plane)
    art3d.pathpatch_2d_to_3d(yz_plane, z=-1, zdir="x")


def plot_axes_arrows(ax: Axes):
    """
    Plots arrows from -1 to 1 along the xaxis (green), y axis (blue) and z axis (red). A dashed line 
    is used for the negative portion the axes. Additionally, rectangle patches are added. The xy plane 
    patch is green and at z=-1. The xz patch is blue and is at y=-1. The yz patch is red and at x=-1.
    """
    x_axis = FancyArrow(0, 0, 1, 0, color='green', width=0.01, head_length=0.02, head_width=0.03)
    y_axis = FancyArrow(0, 0, 0, 1, color='blue', width=0.01, head_length=0.02, head_width=0.03)
    z_axis = FancyArrow(0, 0, 0, 1, color='red', width=0.01, head_length=0.02, head_width=0.03)

    ax.add_patch(x_axis)
    art3d.pathpatch_2d_to_3d(x_axis, z=0, zdir="z")

    ax.add_patch(y_axis)
    art3d.pathpatch_2d_to_3d(y_axis, z=0, zdir="z")

    ax.add_patch(z_axis)
    art3d.pathpatch_2d_to_3d(z_axis, z=0, zdir="x")

    x_axis = FancyArrow(0, 0, -1, 0, color='green', width=0.01, head_length=0.02, head_width=0.03, linestyle='--')
    y_axis = FancyArrow(0, 0, 0, -1, color='blue', width=0.01, head_length=0.02, head_width=0.03, linestyle='--')
    z_axis = FancyArrow(0, 0, 0, -1, color='red', width=0.01, head_length=0.02, head_width=0.03, linestyle='--')

    ax.add_patch(x_axis)
    art3d.pathpatch_2d_to_3d(x_axis, z=0, zdir="z")

    ax.add_patch(y_axis)
    art3d.pathpatch_2d_to_3d(y_axis, z=0, zdir="z")

    ax.add_patch(z_axis)
    art3d.pathpatch_2d_to_3d(z_axis, z=0, zdir="x")


def plot_limits(ax: Axes, lim: float = 1.0):
    """
    Sets all threes plot axes to the same limit. 
    
    Parameters
    ----------
    ax : Axes
        The axes to set the limits of.
    lim : float, optional
        A float indicating the half range to set the plot limits to. All three
        axes will extend from -`lim` to +`lim`. The default is 1.
        
    """
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)


def basic_setup(rows: int = 1, cols: int = 1, use_box: bool = False):
    """
    Creates a basic 3D plot, with arrows along the axes and colored rectange patches 
    created by `plot_box_planes()` or `plot_axes_planes()`. 
    
    Parameters
    ----------
    rows : int, optional
        The number of rows of axes to make. The default is 1.
    cols : int, optional
        The number of columns of axes to make. The default is 1.
    use_box : bool, optional
        If false, then `plot_axes_planes()` is used, otherwise `plot_box_planes()`
        is used. The default is False.
    """
    fig, ax = make_plot(rows, cols)
    if not isinstance(ax, np.ndarray):
        ax = np.array([ax])
    else:
        ax = ax.flatten()
    for a in ax:
        print(a)
        plot_axes_arrows(a)
        if use_box:
            plot_box_planes(a)
        else:
            plot_axes_planes(a)
        plot_limits(a)
    return fig, ax


def plot_vector(ax, x, y, z, ox=0, oy=0, oz=0, fmt='k', trace=True, **kwargs):
    """
    Plots a vector, with optional traces from the xy plane to the tip of the
    vector. The traces 
    
    Parameters
    ----------
    ax : Axes
        The axes to plot to.
    x : float
        The x coordinate of the vector.
    y : float
        The y coordinate of the vector.
    z : float
        The z coordinate of the vector.
    ox : float, optional
        The x coordinate of the vector origin. The default is 0.
    oy : float, optional
        The y coordinate of the vector origin. The default is 0.
    oz : float, optional
        The z coordinate of the vector origin. The default is 0.
    fmt : str, optional
        A format string passed to matplotlib, that formats the 
        line used to plot the vector. The default is 'k', a solid
        black line.
    trace : bool, optional
        If True then thin dotted lines are drawn from the tip of the 
        vector to the xy plane and along the projection of the vector 
        on the xy plane. The default is True.
    kwargs : dict
        A dictionary of arguments to be passed to matplotlib.plot3D().
    """
    ax.plot3D([ox, x], [oy, y], [oz, z], fmt, **kwargs)
    if trace:
        ax.plot3D([ox, x], [oy, y],  [0, 0], ':k', lw=0.3)  # in x,y plane
        ax.plot3D( [x, x],  [y, y], [oz, z], ':k', lw=0.3)  # z component