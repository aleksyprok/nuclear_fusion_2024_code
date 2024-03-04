"""
Test the fstate module.
"""
import os
from python_scripts import fstate

def test_particle_group_init():
    """
    Test the ParticleGroup class initialization.
    """
    test_particle_group = fstate.ParticleGroup()
    assert test_particle_group.r.size < 1
    assert test_particle_group.phi.size < 1
    assert test_particle_group.z.size < 1
    assert test_particle_group.vr.size < 1
    assert test_particle_group.vphi.size < 1
    assert test_particle_group.vz.size < 1
    assert test_particle_group.t.size < 1
    assert test_particle_group.s.size < 1
    assert test_particle_group.particle_id.size < 1
    assert test_particle_group.r0.size < 1
    assert test_particle_group.phi0.size < 1
    assert test_particle_group.z0.size < 1
    assert test_particle_group.vr0.size < 1
    assert test_particle_group.vphi0.size < 1
    assert test_particle_group.vz0.size < 1
    assert test_particle_group.weight.size < 1
    # assert not test_particle_group.phi
    # assert not test_particle_group.z
    # assert not test_particle_group.vr
    # assert not test_particle_group.vphi
    # assert not test_particle_group.vz
    # assert not test_particle_group.t
    # assert not test_particle_group.s
    # assert not test_particle_group.particle_id
    # assert not test_particle_group.r0
    # assert not test_particle_group.phi0
    # assert not test_particle_group.z0
    # assert not test_particle_group.vr0
    # assert not test_particle_group.vphi0
    # assert not test_particle_group.vz0
    # assert not test_particle_group.weight

def test_particle_group_add_particles():
    """
    Test the ParticleGroup class add_particle method.
    """
    def test_particle_group_add_particles():
        """
        Test the ParticleGroup class add_particles method.
        """
        test_particle_group = fstate.ParticleGroup()
        data = [
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0],
            [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
        ]
        test_particle_group.add_particles(data)
        assert test_particle_group.r == [1.0, 2.0]
        assert test_particle_group.phi == [2.0, 3.0]
        assert test_particle_group.z == [3.0, 4.0]
        assert test_particle_group.vr == [4.0, 5.0]
        assert test_particle_group.vphi == [5.0, 6.0]
        assert test_particle_group.vz == [6.0, 7.0]
        assert test_particle_group.t == [7.0, 8.0]
        assert test_particle_group.s == [8.0, 9.0]
        assert test_particle_group.particle_id == [9.0, 10.0]
        assert test_particle_group.r0 == [10.0, 11.0]
        assert test_particle_group.phi0 == [11.0, 12.0]
        assert test_particle_group.z0 == [12.0, 13.0]
        assert test_particle_group.vr0 == [13.0, 14.0]
        assert test_particle_group.vphi0 == [14.0, 15.0]
        assert test_particle_group.vz0 == [15.0, 16.0]
        assert test_particle_group.weight == [16.0, 17.0]

def test_fstate_init():
    """
    Test the Fstate class initialization.
    """
    repo_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    input_dir = os.path.join(repo_path, "input_data", "LOCUST_SPR-045-14_OutputFiles",
                            "axisymmetric", "gpu-q-41")
    fstate_path = os.path.join(input_dir, "FINAL_STATE_13-12-2023_16-51-52.811.dat")
    axisymmetric_fstate = fstate.Fstate(fstate_path)
    assert axisymmetric_fstate.fstate_path == fstate_path
    assert len(axisymmetric_fstate.moving.r) == 18
    assert len(axisymmetric_fstate.thermal.r) == 482343
    assert len(axisymmetric_fstate.stopped.r) == 41927
    assert len(axisymmetric_fstate.unresolved.r) == 0
    assert abs(axisymmetric_fstate.stopped.r[0] - 5.77736) < 1e-5
    assert abs(axisymmetric_fstate.stopped.phi[0] + 1.56278) < 1e-5
    assert abs(axisymmetric_fstate.stopped.z[0] - 8.53939) < 1e-5
    assert abs(axisymmetric_fstate.stopped.vr[0] + 114536.0) < 1e-5
    assert abs(axisymmetric_fstate.stopped.vphi[0] + 577659.0) < 1e-5
    assert abs(axisymmetric_fstate.stopped.vz[0] - 347993.0) < 1e-5
    assert abs(axisymmetric_fstate.stopped.t[0] - 0.833556) < 1e-5
    assert abs(axisymmetric_fstate.stopped.s[0] + 5.0) < 1e-5
    assert int(axisymmetric_fstate.stopped.particle_id[0]) == 64954
    assert abs(axisymmetric_fstate.stopped.r0[0] - 5.12387) < 1e-5
    assert abs(axisymmetric_fstate.stopped.phi0[0] - 2.37106) < 1e-5
    assert abs(axisymmetric_fstate.stopped.z0[0] - 2.5327) < 1e-5
    assert abs(axisymmetric_fstate.stopped.vr0[0] - 8201490.0) < 1e-5
    assert abs(axisymmetric_fstate.stopped.vphi0[0] + 8432840.0) < 1e-5
    assert abs(axisymmetric_fstate.stopped.vz0[0] + 5514400.0) < 1e-5
    assert abs(axisymmetric_fstate.stopped.weight[0] - 23987000000.0) < 1e-5
