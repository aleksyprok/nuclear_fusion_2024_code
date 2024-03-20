"""
Module contains the routines for parameterizing the wall.
"""

from multiprocessing import Process, Manager, cpu_count, shared_memory
import numpy as np
from shapely.geometry import Point, Polygon, LineString
from shapely.ops import nearest_points
from scipy.spatial.distance import cdist

def calc_wall_centroid(r_wall, z_wall):
    """
    calculate the centroid of the wall using formula given here:
    https://en.wikipedia.org/wiki/centroid#of_a_polygon
    """

    a = 0.5 * np.sum(r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1])

    cr = np.sum((r_wall[:-1] + r_wall[1:]) * (r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1]))
    cr /= 6 * a
    cz = np.sum((z_wall[:-1] + z_wall[1:]) * (r_wall[:-1] * z_wall[1:] - r_wall[1:] * z_wall[:-1]))
    cz /= 6 * a

    return cr, cz

def get_s_theta_from_rz(r_coord, z_coord, r_wall, z_wall):
    """
    Calculate the s_theta coordinate from a set of 
    r and z coords which lie on the wall

    Args:
        r_coord: np.ndarray
            The r coordinates of the markers
        z_coord: np.ndarray
            The z coordinates of the markers
        r_wall: np.ndarray
            The r coordinates of the wall
        z_wall: np.ndarray
            The z coordinates of the wall
    """
    def calculate_s_theta(n_array):

        for n in n_array:

            # Calculate r_new, z_new which give the closest points on the wall
            point = Point(r_coord[n], z_coord[n])
            p1, _ = nearest_points(poly.boundary, point)
            r_new = p1.x
            z_new = p1.y

            # Determine the wall vertices either side of r_new, z_new
            nearest_node = np.argmin(cdist([[r_new, z_new]], wall_coords))
            if nearest_node == 0:
                coords0 = [r_wall[-2], z_wall[-2]]
            else:
                coords0 = [r_wall[nearest_node - 1], z_wall[nearest_node - 1]]
            coords1 = [r_wall[nearest_node], z_wall[nearest_node]]
            coords2 = [r_wall[nearest_node + 1], z_wall[nearest_node + 1]]

            # Check if the point is on line connecting coords0 with coords1
            if LineString((coords0, coords1)).distance(Point(r_new, z_new)) < 1e-8:
                if nearest_node == 0:
                    s_theta[n] = s_theta_nodes[-1] \
                               + cdist([coords0], [[r_new, z_new]])[0][0]
                else:
                    s_theta[n] = s_theta_nodes[nearest_node - 1] \
                               + cdist([coords0], [[r_new, z_new]])[0][0]

            # Check if the point is on line connecting coords1 with coords2
            elif LineString((coords1, coords2)).distance(Point(r_new, z_new)) < 1e-8:
                s_theta[n] = s_theta_nodes[nearest_node] \
                           + cdist([coords1], [[r_new, z_new]])[0][0]

            # In the unlikely event that the particle is not on the wall connected to the
            # nearest vertices then we check all walls. Note that this is more
            # computationally expensive.
            else:
                index = -1
                for coords1, coords2 in zip(wall_coords[:-1], wall_coords[1:]):
                    index += 1
                    if LineString((coords1, coords2)).distance(Point(r_new, z_new)) < 1e-8:
                        s_theta[n] = s_theta_nodes[index] \
                                   + cdist([coords1], [[r_new, z_new]])[0][0]
                        # print('This point did not lie on walls connected to nearest vertex:', \
                        #       r_coord[n], z_coord[n])
                        break
                    if index == len(r_wall) - 2:
                        print('Error: Unable to calculate s_theta.')

    print('\nCalculating s_theta from r_coords and z_coords.')

    wall_coords = np.vstack([r_wall, z_wall]).T
    poly = Polygon(wall_coords)

    ds = np.sqrt(np.diff(r_wall)**2 + np.diff(z_wall)**2)
    s_theta_nodes = np.zeros(len(r_wall) - 1)
    s_theta_nodes[1:] = np.cumsum(ds[:-1])

    shm = shared_memory.SharedMemory(create = True, \
                                     size = np.zeros(len(r_coord)).nbytes)
    s_theta = np.ndarray(len(r_coord), dtype = float, buffer = shm.buf)
    n_array_list = np.array_split(np.arange(len(r_coord)), cpu_count())
    with Manager():
        processes = []

        for n_array in n_array_list:
            p = Process(target = calculate_s_theta, \
                args = (n_array,))
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

    s_theta = np.array(s_theta)
    shm.close()
    shm.unlink()

    return s_theta
