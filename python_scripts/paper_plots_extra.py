"""
This file contains plots which are not shown in the IAEA FEC 2024 paper,
but are generated at the same time as the paper plots.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

def plot_amise_2d(run, output_dir):
    """
    Plot the Asymptotic Mean Integrated Square error array.
    """
    hx_mesh, hy_mesh = np.meshgrid(run.flux.h_phi_array,
                                   run.flux.h_theta_2d_array)
    fig, ax = plt.subplots()
    im = ax.pcolormesh(hx_mesh, hy_mesh, run.flux.amise_2d,
                       norm=LogNorm())
    ax.plot(run.flux.h_phi, run.flux.h_theta_2d,
            'r.', markersize=10)
    ax.set_title('Asymptotic mean integrated square error')
    ax.set_xscale('log')
    ax.set_yscale('log')
    fig.colorbar(im, ax=ax)
    ax.set_xlabel(r'$h_\phi$ [m]')
    ax.set_ylabel(r'$h_\theta$ [m]')
    fig.savefig(output_dir + '/amise_2d_array.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

    fig, axs = plt.subplots(1, 2)
    fig_size = fig.get_size_inches()
    fig_size[0] *= 2
    fig.set_size_inches(fig_size)

    h_theta_index = np.where(run.flux.h_theta_2d_array == run.flux.h_theta_2d)[0][0]
    h_phi_index = np.where(run.flux.h_phi_array == run.flux.h_phi)[0][0]
    axs[0].plot(run.flux.h_phi_array, run.flux.amise_2d[h_theta_index, :])
    axs[0].plot(run.flux.h_phi, run.flux.amise_2d[h_theta_index, h_phi_index],
                'r.', markersize=10)
    axs[0].set_xlabel(r'$h_\phi$ [m]')
    axs[0].set_ylabel(r'$AMISE$')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_title(r'$h_\theta$ = ' f'{run.flux.h_theta_2d} m')

    axs[1].plot(run.flux.h_theta_2d_array, run.flux.amise_2d[:, h_phi_index])
    axs[1].plot(run.flux.h_theta_2d, run.flux.amise_2d[h_theta_index, h_phi_index],
                'r.', markersize=10)
    axs[1].set_xlabel(r'$h_\theta$ [m]')
    axs[1].set_ylabel(r'$AMISE$')
    axs[1].set_xscale('log')
    axs[1].set_yscale('log')
    axs[1].set_title(r'$h_\phi$ = ' f'{run.flux.h_phi} m')
    fig.savefig(output_dir + '/amise_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def plot_energy_flux_2d(run, output_dir):
    """
    Plot the energy flux array.
    """
    fig, ax = plt.subplots()
    ax.imshow(run.flux.energy_2d.T, origin='lower',
               extent=[run.wall.s_phi_min, run.wall.s_phi_max,
                       run.wall.s_theta_min, run.wall.s_theta_max])
    ax.set_xlabel('s_phi [m]')
    ax.set_ylabel('s_theta [m]')
    ax.set_title('Energy flux [MW/m^2]\n'
                 f'h_phi = {run.flux.h_phi:.2e}, h_theta_2d = {run.flux.h_theta_2d:.2e}')
    fig.savefig(output_dir + '/energy_flux_2d.png',
                bbox_inches='tight', dpi=300)

    fig, ax = plt.subplots()
    ax.plot(run.flux.s_theta_2d, np.max(run.flux.energy_2d, axis=0))
    for s_nod_i in run.wall.special_nodes:
        ax.axvline(x=run.wall.s_nodes[s_nod_i], color='k')
    ax.set_xlabel('s_theta [m]')
    ax.set_title('Max Energy flux [MW/m^2]\n'
                 f'h_phi = {run.flux.h_phi:.2e}, h_theta_2d = {run.flux.h_theta_2d:.2e}')
    fig.savefig(output_dir + '/max_energy_flux_2d_line.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

def save_attributes_to_file(obj, output_dir, indent=0):
    """
    Recursively saves attributes and their values for a given class instance to a file.
    If an attribute is a class instance, it does the same for that class.

    Parameters:
    obj (object): The class instance.
    output_dir (str): Directory to the output file.
    indent (int): The indentation level for pretty printing.
    """
    file_path = os.path.join(output_dir, 'run_i_attributes.txt')
    with open(file_path, "w", encoding="utf-8") as file:
        def write_attributes(obj, indent=0):
            indent_str = '    ' * indent

            if not hasattr(obj, "__dict__"):
                file.write(f"{indent_str}{obj} is not a class instance.\n")
                return

            for attr, value in vars(obj).items():
                if hasattr(value, "__dict__"):
                    file.write(f"{indent_str}{attr}:\n")
                    write_attributes(value, indent + 1)
                else:
                    file.write(f"{indent_str}{attr}: {value}\n")

        write_attributes(obj, indent=indent)
