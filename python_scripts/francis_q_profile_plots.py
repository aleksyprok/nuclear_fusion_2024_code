"""
This file contains functions used to make the plots to compare results from
SPR-045-16 and SPR-046-16b.
"""

import os
import matplotlib.pyplot as plt

from python_scripts import my_gfile_reader

def compare_q_profiles():
    """
    Compare the Q-profiles of SPR-045-16 vs SPR046_16b vs 
    """
    repo_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    plot_dir = os.path.join(repo_dir, 'plots', 'francis_q_profiles_plots')
    os.makedirs(plot_dir, exist_ok=True)
    eqdsk_path = os.path.join(repo_dir, 'input_data', 'SPR-045-16.eqdsk')
    gfile_spr045_16 = my_gfile_reader.getGfile(eqdsk_path)
    eqdsk_path = os.path.join(repo_dir, 'input_data', 'SPR046_16b.eqdsk')
    gfile_spr046_16b = my_gfile_reader.getGfile(eqdsk_path)
    eqdsk_path = os.path.join(repo_dir, 'input_data', 'SPR-045-14.eqdsk')
    gfile_spr045_14 = my_gfile_reader.getGfile(eqdsk_path)

    fig, ax = plt.subplots()
    ax.plot(gfile_spr045_16.psin, gfile_spr045_16.q,
            label='SPR-045-16')
    ax.plot(gfile_spr046_16b.psin, gfile_spr046_16b.q,
            label='SPR-046-16b')
    ax.plot(gfile_spr045_14.psin, gfile_spr045_14.q,
            label='SPR-045-14')
    ax.legend()
    ax.set_xlabel(r'$\psi_n$')
    ax.set_ylabel(r'$q$')
    fig.savefig(plot_dir + '/q_profile_compare.png',
                bbox_inches='tight', dpi=300)
    plt.close('all')

if __name__ == '__main__':
    compare_q_profiles()
    