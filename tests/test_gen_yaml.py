#! /bin/env python3

from __future__ import division
import unittest

import tempfile
from pathlib import Path

import numpy as np
import pytest

from solsticepy.gen_yaml import Sun, gen_yaml, yamltransform

@pytest.fixture(autouse=True)
def _isolated_workdir(tmp_path, monkeypatch):
	monkeypatch.chdir(tmp_path)


def test_1():
	s = yamltransform([1,2,3],[4,5,6])
	#print(s)
	assert(s=='transform: { translation: [1.000000e+00,2.000000e+00,3.000000e+00], rotation: [4.000000e+00,5.000000e+00,6.000000e+00] }')




def test_quadricxy_uses_azimuth_elevation_tracking():
	sun = Sun(dni=1000, sunshape='pillbox')
	hst_pos = np.array([[0., 100., 8.54], [20., 120., 8.54]])
	hst_aims = np.array([[0., 0., 62.], [0., 0., 62.]])
	hst_foc = np.array([[100., 110.], [120., 130.]])
	rec_param = np.r_[20., 20., 10, 10, 0., 0., 62., 0.]
	with tempfile.TemporaryDirectory() as tmpdir:
		tmp_path = Path(tmpdir)
		outfile_yaml = tmp_path / 'input.yaml'
		outfile_recv = tmp_path / 'input-rcv.yaml'

		gen_yaml(sun, hst_pos, hst_foc, hst_aims, 12.2, 12.2, 1., 1e-3, 'flat',
			rec_param, 1., outfile_yaml=outfile_yaml, outfile_recv=outfile_recv,
			shape='quadricxy')

		text = outfile_yaml.read_text()
		assert text.count('quadricxy:') == 2
		assert 'target_aligned' not in text
		assert 'ax2: 2.500000e-03' in text
		assert 'ay2: 1.923077e-03' in text


def test_canted_quadricxy_uses_azimuth_elevation_tracking():
	sun = Sun(dni=1000, sunshape='pillbox')
	hst_pos = np.array([[0., 100., 8.54]])
	hst_aims = np.array([[0., 0., 62.]])
	hst_foc = np.array([[120., 140.]])
	rec_param = np.r_[20., 20., 10, 10, 0., 0., 62., 0.]
	with tempfile.TemporaryDirectory() as tmpdir:
		tmp_path = Path(tmpdir)
		outfile_yaml = tmp_path / 'input.yaml'
		outfile_recv = tmp_path / 'input-rcv.yaml'

		gen_yaml(sun, hst_pos, hst_foc, hst_aims, 12.2, 12.2, 1., 1e-3, 'flat',
			rec_param, 1., outfile_yaml=outfile_yaml, outfile_recv=outfile_recv,
			cant=True, cant_foc=100., fct_w=6.0, fct_h=6.0, fct_gap=0.1, n_row=2, n_col=2,
			shape='quadricxy')

		text = outfile_yaml.read_text()
		assert text.count('quadricxy:') == 4
		assert 'geometry: *' not in text.split('- entity:', 2)[1]
		assert 'target_aligned' not in text
		assert 'ax2: 2.083333e-03' in text
		assert 'ay2: 1.785714e-03' in text


def test_canted_entity_matches_reference(tmp_path):
    """Match H_0 from the supplied input.yaml, including nested facets."""
    gen_yaml(Sun(), np.array([[-159.2729, 222.1405, 5.596018]]),
             np.array([273.3964]), np.array([[-4.369563, 6.094298, 154.1722]]),
             5.23, 3.4, 0.8096, 0.002, 'cylinder',
             [14.99, 17.56, 50, 20, 0., 0., 150., 0.], 0.95,
             tmp_path / 'input.yaml', tmp_path / 'receiver.yaml',
             cant=True, cant_foc=1e6, 
             fct_w=2.3, fct_h=3.3, fct_gap=0.315, n_row=1, n_col=2)
    text = (tmp_path / 'input.yaml').read_text()
    first_entity = text.split('- entity:', 2)[1]
    expected = (Path(__file__).parent / 'fixtures' / 'inline_heliostat.yaml').read_text()
    # Original entity names, tracking and facet transforms are preserved.
    for line in expected.splitlines():
        if line.strip().startswith(('name:', '- name:', 'transform:', 'target:')):
            assert line in '- entity:' + first_entity
    assert first_entity.count('focal: 2.733964e+02') == 2
    assert 'geometry: *' not in first_entity
    assert '- template:' not in text
    assert 'material_mirror' not in text
    assert '&facet_g_band_0' not in text


@pytest.mark.parametrize('shape', ['flat', 'paraboloid', 'parabolic-cylinder', 'sphere', 'quadricxy'])
@pytest.mark.parametrize('cant', [False, True])
def test_heliostats_have_inline_children(tmp_path, monkeypatch, shape, cant):
    monkeypatch.chdir(tmp_path)
    gen_yaml(Sun(), np.array([[0., 100., 8.], [20., 120., 8.], [30., 140., 8.]]),
             np.array([90., 95., 180.]), np.array([[0., 0., 62.]] * 3),
             12.2, 12.2, 0.95, 0.001, 'flat',
             [20., 20., 10, 10, 0., 0., 62., 0.], 0.9,
             tmp_path / 'input.yaml', tmp_path / 'receiver.yaml',
             cant=cant, shape=shape,
             fct_w=6., fct_h=6., fct_gap=0.1, n_row=2, n_col=2)
    text = (tmp_path / 'input.yaml').read_text()
    assert '- template:' not in text
    assert '*hst_t_' not in text
    assert '*facets_t_' not in text
    for i in range(3):
        assert '    children:\n      - name: hst_%d\n' % i in text
    if cant:
        assert text.count('- name: facets\n') == 3
        assert text.count('- name: facet_0_c0_r0\n') == 1
        assert text.count('- name: facet_1_c0_r0\n') == 1
        assert text.count('- name: facet_2_c0_r0\n') == 1
    entities = text.split('- entity:')[1:4]
    for entity, focal in zip(entities, [90., 95., 180.]):
        assert 'geometry: *' not in entity
        count = 4 if cant else 1
        if shape in ('paraboloid', 'parabolic-cylinder'):
            assert entity.count('focal: %e' % focal) == count
        elif shape == 'sphere':
            assert entity.count('radius: %e' % (2 * focal)) == count
        elif shape == 'quadricxy':
            assert entity.count('ax2: %e' % (1 / (4 * focal))) == count
    assert 'band_' not in text


def test_target_aligned_has_inline_children(tmp_path):
    from solsticepy.gen_yaml import gen_yaml_target_aligned
    gen_yaml_target_aligned(Sun(), np.array([[0., 100., 8.]]),
                            np.array([[0., 0., 62.]]), [100.], [120.],
                            12., 12., 0.95, 0.001, 'flat',
                            [20., 20., 10, 10, 0., 0., 62., 0.], 0.9,
                            tmp_path / 'input.yaml', tmp_path / 'receiver.yaml')
    text = (tmp_path / 'input.yaml').read_text()
    assert '- template:' not in text
    assert '    children:\n      - name: hst_0\n' in text
    assert 'target_aligned: true' in text
    assert 'geometry: *hst_g_0' not in text
    assert 'ax2: 2.500000e-03' in text
    assert 'ay2: 2.083333e-03' in text


@pytest.mark.parametrize('focal', [105., np.array([105., 105.])])
@pytest.mark.parametrize('one_heliostat', [True, False])
def test_actual_focal_without_bands(tmp_path, focal, one_heliostat):
    gen_yaml(Sun(), 
        hst_pos=np.array([[0., 100., 8.], [20., 120., 8.]]), 
        hst_foc=focal,
        hst_aims=np.array([[0., 0., 62.]] * 2), 
        hst_w=12., 
        hst_h=12., 
        rho_refl=.95,
        slope_error= np.array([.001, .002]), 
        receiver='flat', 
        rec_param=[20., 20., 10, 10, 0., 0., 62., 0.], 
        rec_abs=.9,
        outfile_yaml=tmp_path / 'input.yaml', 
        outfile_recv=tmp_path / 'receiver.yaml', 
        one_heliostat=one_heliostat)
    text = (tmp_path / 'input.yaml').read_text()
    assert text.count('focal: 1.050000e+02') == (1 if one_heliostat else 2)
    assert ('slope_error: %15.8e' % .001) in text
    if not one_heliostat:
        assert ('slope_error: %15.8e' % .002) in text
    assert 'band' not in text


def test_first_quadric_uses_first_directional_focal(tmp_path):
    gen_yaml(Sun(), np.array([[0., 100., 8.], [20., 120., 8.]]),
             np.array([[100., 120.], [200., 240.]]), np.array([[0., 0., 62.]] * 2),
             12., 12., .95, .001, 'flat', [20., 20., 10, 10, 0., 0., 62., 0.], .9,
             tmp_path / 'input.yaml', tmp_path / 'receiver.yaml',
             shape='quadricxy', one_heliostat=True)
    text = (tmp_path / 'input.yaml').read_text()
    assert text.count('quadricxy:') == 1
    assert 'ax2: 2.500000e-03' in text
    assert 'ay2: 2.083333e-03' in text


def test_canting_focus_is_independent_of_optical_focus(tmp_path):
    from solsticepy.gen_yaml import heliostat_canted_facets
    gen_yaml(Sun(), np.array([[0., 100., 8.], [20., 120., 8.]]),
             np.array([101., 123.]), np.array([[0., 0., 62.]] * 2),
             12.2, 12.2, .95, np.array([.001, .002]), 'flat',
             [20., 20., 10, 10, 0., 0., 62., 0.], .9,
             tmp_path / 'input.yaml', tmp_path / 'receiver.yaml',
             cant=True, cant_foc=[500., 1000.], fct_w=6., fct_h=6.,
             fct_gap=.1, n_row=2, n_col=2)
    entities = (tmp_path / 'input.yaml').read_text().split('- entity:')[1:3]
    for i, (entity, cant_focus, focal) in enumerate(zip(entities, [500., 1000.], [101., 123.])):
        data = heliostat_canted_facets(12.2, 12.2, 6., 6., .1, 2, 2, cant_focus)
        x, z, rx, ry = data[0]
        assert yamltransform([x, 0, z], [rx, ry, 0]) in entity
        assert entity.count('focal: %e' % focal) == 4
        assert entity.count('slope_error: %15.8e' % [.001, .002][i]) == 5


@pytest.mark.parametrize('mode', ['single', 'canted', 'target_aligned'])
@pytest.mark.parametrize('rho,slope', [(.91, .001), ([.81, .92345678], .001),
                                     (.91, [.001, .003]), ([.81, .92345678], [.001, .003])])
def test_embedded_materials(tmp_path, mode, rho, slope):
    from solsticepy.gen_yaml import gen_yaml_target_aligned
    positions = np.array([[0., 100., 8.], [20., 120., 8.]])
    aims = np.array([[0., 0., 62.]] * 2)
    kwargs = dict(hst_pos=positions, hst_aims=aims, sun=Sun(), hst_w=12.2, hst_h=12.2,
                  rho_refl=rho, slope_error=slope, receiver='flat',
                  rec_param=[20., 20., 10, 10, 0., 0., 62., 0.], rec_abs=.9,
                  outfile_yaml=tmp_path / 'input.yaml', outfile_recv=tmp_path / 'receiver.yaml')
    if mode == 'target_aligned':
        gen_yaml_target_aligned(Fx=[101., 123.], Fy=[111., 133.], **kwargs)
    else:
        gen_yaml(hst_foc=[101., 123.], cant=mode == 'canted', fct_w=6., fct_h=6.,
                 fct_gap=.1, n_row=2, n_col=2, **kwargs)
    text = (tmp_path / 'input.yaml').read_text()
    assert 'material_mirror' not in text
    assert 'mirror:' not in text.split('- entity:')[0]
    for i, entity in enumerate(text.split('- entity:')[1:3]):
        reflectivity = rho if np.ndim(rho) == 0 else rho[i]
        error = slope if np.ndim(slope) == 0 else slope[i]
        expected = 'mirror: {reflectivity: %15.8e, slope_error: %15.8e}' % (reflectivity, error)
        assert entity.count(expected) == (5 if mode == 'canted' else 1)
        assert 'material:\n' in entity


@pytest.mark.parametrize('name', ['rho_refl', 'slope_error'])
def test_invalid_material_array_has_clear_error(tmp_path, name):
    kwargs = dict(rho_refl=.9, slope_error=.001)
    kwargs[name] = [.1, .2, .3]
    with pytest.raises(ValueError, match=name + ' must be a scalar or one value per heliostat'):
        gen_yaml(Sun(), np.array([[0., 100., 8.], [20., 120., 8.]]), [100., 120.],
                 np.array([[0., 0., 62.]] * 2), 12., 12., receiver='flat',
                 rec_param=[20., 20., 10, 10, 0., 0., 62., 0.], rec_abs=.9,
                 outfile_yaml=tmp_path / 'input.yaml', outfile_recv=tmp_path / 'receiver.yaml', **kwargs)
