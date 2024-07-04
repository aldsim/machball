import machball.ballistic.viewfactors as vf

import pytest


def test_par_rec_to_rec():
    xa = (0,1)
    ya = (0,1)
    xb = (0,1)
    yb = (0,1)
    z = 1
    F = vf.par_rec_to_rec(xa, ya, xb, yb, z)
    assert F == pytest.approx(vf.square_to_square(1.0,1.0))


def test_rect_to_side():
    assert vf.rect_to_side(1,1,1) == pytest.approx(0.2000438)


def test_perp_rec_to_rec():
    xa = (1e-7,1)
    ya = (0, 1)
    yb = (0, 1)
    zb = (1e-7,1)
    F = vf.perp_rec_to_rec(xa, ya, yb, zb)
    assert F == pytest.approx(vf.rect_to_side(1,1,1))


def test_cylinder_to_base():
    assert vf.cylinder_base_to_base(1, 1) == pytest.approx(vf.disk_to_disk(1,1))


def test_disk_to_disk():
    assert vf.disk_to_disk(1, 1) == pytest.approx(vf.disk1_to_disk2(1, 1, 1))


def test_strip_to_strip():
    F1 = vf.strip_to_strip(1,1)
    F2 = vf.strip_to_normalstrip(1,1)
    assert F1 + 2*F2 == pytest.approx(1)


def test_square_to_square():
    F1 = vf.square_to_square(1,2)
    F2 = vf.rectangle_to_rectangle(1, 1, 2)
    assert F1 == pytest.approx(F2)
