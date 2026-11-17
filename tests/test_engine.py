"""Offline checks of the Bone Guide engine. No Fusion.

    py -3.10 tests/test_engine.py

The bone is a cylinder of radius 10, so every answer is known in advance: the
bone normal is radial, a base of clearance c and thickness t lies between radii
10 + c and 10 + c + t, its wall must be radial all round, and its volume is the
outline's area in (angle, z) times (R2^2 - R1^2) / 2. The outline wraps 240 deg
round the cylinder - the side tabs a lofted guide cannot keep at thickness.
"""

import math
import os
import sys

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, HERE)

import numpy as np                                            # noqa: E402
import trimesh                                                # noqa: E402

import bg_engine as E                                         # noqa: E402

FAILED = []
R = 10.0


def check(name, cond, detail=''):
    print('  %s %s%s' % ('ok  ' if cond else 'FAIL', name,
                         '   %s' % (detail,) if detail != '' else ''), flush=True)
    if not cond:
        FAILED.append(name)


def bone(h=0.6, length=60.0):
    """An open tube of radius R on a regular grid of ~h mm: every vertex exactly
    on the cylinder. (Subdividing trimesh's cylinder instead splits its long
    sliver sides into 12 million faces.)"""
    nt = int(round(2 * math.pi * R / h))
    nz = int(round(length / h)) + 1
    th = np.linspace(0, 2 * math.pi, nt, endpoint=False)
    z = np.linspace(-length / 2, length / 2, nz)
    T, Z = np.meshgrid(th, z, indexing='ij')
    V = np.column_stack([R * np.cos(T.ravel()), R * np.sin(T.ravel()), Z.ravel()])
    idx = np.arange(nt * nz).reshape(nt, nz)
    a, b = idx[:, :-1], np.roll(idx, -1, axis=0)[:, :-1]
    c, d = idx[:, 1:], np.roll(idx, -1, axis=0)[:, 1:]
    F = np.vstack([np.stack([a.ravel(), b.ravel(), d.ravel()], 1),
                   np.stack([a.ravel(), d.ravel(), c.ravel()], 1)])
    return trimesh.Trimesh(V, F, process=False)      # (a, b, d): normals outward


def rounded_rect(w, h, rc, n=2000):
    """(u, v) outline of a w x h rectangle with corners of radius rc, CCW."""
    pts = []
    cx, cy = w / 2 - rc, h / 2 - rc
    for (sx, sy, a0) in ((1, 1, 0), (-1, 1, 90), (-1, -1, 180), (1, -1, 270)):
        for k in range(n // 4):
            a = math.radians(a0 + 90.0 * k / (n // 4))
            pts.append((sx * cx + rc * math.cos(a), sy * cy + rc * math.sin(a)))
    return np.array(pts)


def to_cyl(uv, radius=R):
    """Wrap a (R theta, z) outline onto the cylinder, densely: a straight side in
    (theta, z) is a helix on the cylinder, not the chord between its ends."""
    uv = E.resample(np.column_stack([uv, np.zeros(len(uv))]), 0.2)[:, :2]
    th = uv[:, 0] / R
    return np.column_stack([radius * np.cos(th), radius * np.sin(th), uv[:, 1]])


def polygon_area(uv):
    x, y = uv[:, 0], uv[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def dist_to_polyline(P, L):
    """Distance from 2D points P to the closed 2D polyline L."""
    A, B = L, np.roll(L, -1, axis=0)
    best = np.full(len(P), np.inf)
    for a, b in zip(A, B):
        ab = b - a
        t = np.clip(((P - a) @ ab) / max(ab @ ab, 1e-18), 0, 1)
        best = np.minimum(best, np.linalg.norm(P - (a + t[:, None] * ab), axis=1))
    return best


def test_attach(m, b):
    print('attach')
    th = np.linspace(0, 2 * math.pi, 50, endpoint=False)
    X = np.column_stack([11.0 * np.cos(th), 11.0 * np.sin(th), 5 * np.sin(3 * th)])
    P, N = E.attach(b, X, 0.1)
    r = np.linalg.norm(P[:, :2], axis=1)
    check('points 1 mm off the cylinder land on it', np.abs(r - R).max() < 0.01,
          '%.4f' % np.abs(r - R).max())
    rad = np.column_stack([P[:, 0], P[:, 1], 0 * P[:, 2]]) / R
    check('and carry the radial normal', np.min(np.einsum('ij,ij->i', N, rad)) > 0.999,
          '%.5f' % np.min(np.einsum('ij,ij->i', N, rad)))


def base_case(b, uv, c, t, r, holes=()):
    Q = to_cyl(uv)
    P, N = E.attach(b, Q, 0.1)
    base = E.Base(b, P, N, c, t, r, 0.0, 0.1, 1.0)
    for o, a, rad in holes:
        base.add_hole(o, a, rad)
    S = base.mesh(0.1)
    tri = S.faces.reshape(-1, 4)[:, 1:]
    tm = trimesh.Trimesh(S.points, tri, process=False)
    return S, tm, base


def test_base(b):
    print('base on a cylinder, outline wrapping 240 deg')
    c, t = 0.2, 3.0
    w = R * math.radians(240.0)
    uv = rounded_rect(w, 20.0, 4.0)
    S, tm, base = base_case(b, uv, c, t, 0.0)
    rr = np.linalg.norm(S.points[:, :2], axis=1)
    check('one closed piece, facing out', tm.is_watertight and tm.volume > 0
          and len(tm.split(only_watertight=False)) == 1,
          (tm.is_watertight, round(float(tm.volume), 1), len(tm.split(only_watertight=False))))
    check('contact face at the clearance', abs(rr.min() - (R + c)) < 0.01, '%.4f' % rr.min())
    check('top at clearance + thickness', abs(rr.max() - (R + c + t)) < 0.01, '%.4f' % rr.max())
    want = polygon_area(uv) / R * ((R + c + t) ** 2 - (R + c) ** 2) / 2.0
    check('volume = outline area x radial shell (walls along the normal)',
          abs(tm.volume / want - 1) < 0.01, '%.1f vs %.1f mm3' % (tm.volume, want))
    # wall vertices: away from both faces; their (angle, z) must lie ON the outline
    mid = (rr > R + c + 0.4) & (rr < R + c + t - 0.4)
    th = np.arctan2(S.points[mid, 1], S.points[mid, 0])
    uvw = np.column_stack([R * th, S.points[mid, 2]])
    d = dist_to_polyline(uvw, uv[::4])
    check('the wall stands along the bone normal (radial) all round',
          d.max() < 0.08, 'max %.3f mm from the outline in (R theta, z)' % d.max())
    tab = np.abs(np.degrees(np.arctan2(S.points[:, 1], S.points[:, 0]))) > 105
    rt = rr[tab]
    check('the wrapped tabs keep the full thickness',
          abs(rt.max() - rt.min() - t) < 0.02, '%.3f mm' % (rt.max() - rt.min()))
    print('rim round')
    S2, tm2, _ = base_case(b, uv, c, t, 1.0)
    rr2 = np.linalg.norm(S2.points[:, :2], axis=1)
    check('the round never touches the contact face', abs(rr2.min() - (R + c)) < 0.01,
          '%.4f' % rr2.min())
    lost = tm.volume - tm2.volume
    per = sum(np.linalg.norm(np.diff(np.vstack([uv, uv[:1]]), axis=0), axis=1))
    # the round removes r^2 (1 - pi/4) per mm of rim, scaled to the top radius
    exp = (1.0 - math.pi / 4.0) * per * (R + c + t) / R
    check('and takes off what a radius-1 round should', abs(lost / exp - 1) < 0.15,
          '%.1f vs %.1f mm3' % (lost, exp))
    print('holes')
    S3, tm3, _ = base_case(b, uv, c, t, 0.0, holes=[((0, 0, 0), (1, 0, 0), 1.2)])
    p = S3.points
    inner = (np.hypot(p[:, 1], p[:, 2]) < 1.15) & (p[:, 0] > 0)
    check('nothing inside the drill radius', not inner.any(), int(inner.sum()))
    want_h = math.pi * 1.2 ** 2 * t
    check('the hole removes about pi r^2 t', abs((tm.volume - tm3.volume) / want_h - 1) < 0.05,
          '%.2f vs %.2f mm3' % (tm.volume - tm3.volume, want_h))


def plate(half=30.0, h=0.3):
    """A flat open "bone" z = 0, normals +z."""
    n = int(round(2 * half / h)) + 1
    x = np.linspace(-half, half, n)
    X, Y = np.meshgrid(x, x, indexing='ij')
    V = np.column_stack([X.ravel(), Y.ravel(), np.zeros(n * n)])
    idx = np.arange(n * n).reshape(n, n)
    a, b, c, d = idx[:-1, :-1], idx[1:, :-1], idx[1:, 1:], idx[:-1, 1:]
    F = np.vstack([np.stack([a.ravel(), b.ravel(), c.ravel()], 1),
                   np.stack([a.ravel(), c.ravel(), d.ravel()], 1)])
    return E.Bone(V, F, step=0.25)


def test_section():
    """Two drill margins of radius 3, 8 mm apart, standing on a flat bone. Apart
    they leave two circles; closed with rho = 5 the gap is bridged by arcs of
    radius 5 tangent to both, centred at (0, +-sqrt(8^2 - 4^2)): the waist at
    x = 0 is |y| = sqrt(48) - 5 = 1.928."""
    print('section: the alpha shape of two drills, on the bone')
    b = plate()
    cyl = [{'kind': 'cyl', 'p0': [sx, 0, -12], 'p1': [sx, 0, 12], 'r': 2.0} for sx in (-4, 4)]
    loops, info = E.section(b, cyl, margin=1.0, closing=0.0, join=False)
    r = [np.hypot(np.abs(L[:, 0]) - 4, L[:, 1]) for L in loops]
    check('no closing: two circles of radius 2 + margin 1',
          len(loops) == 2 and all(np.abs(x - 3).max() < 0.03 for x in r),
          (len(loops), [round(float(np.abs(x - 3).max()), 4) for x in r]))
    want = math.sqrt(48.0) - 5.0
    for planes, tol in ((False, 0.05), (True, 0.15)):
        loops, info = E.section(b, cyl, margin=1.0, closing=5.0, planes=planes)
        L = loops[0] if loops else np.zeros((0, 3))
        outer = np.abs(L[:, 0]) > 4
        ro = np.hypot(np.abs(L[outer, 0]) - 4, L[outer, 1])
        waist = np.abs(L[:, 0]) < 0.3
        how = info['closing']
        check('closing 5 (%s): one loop' % how, len(loops) == 1, info)
        check('(%s) the far sides stay on the margin circles' % how,
              np.abs(ro - 3).max() < tol, '%.3f mm' % np.abs(ro - 3).max())
        check('(%s) the gap is bridged by radius-5 arcs (waist |y| = %.3f)' % (how, want),
              waist.any() and np.abs(np.abs(L[waist, 1]) - want).max() < tol,
              '%.3f' % np.abs(np.abs(L[waist, 1]) - want).max() if waist.any() else 'none')
    print('spline through fit points')
    th = np.linspace(0, 2 * math.pi, 12, endpoint=False)
    fit = np.column_stack([8 * np.cos(th), 8 * np.sin(th), np.full(12, 0.7)])
    out = E.spline_answer(b, fit, 0.25)
    P, F, at = np.asarray(out['pts']), np.asarray(out['fit']), np.asarray(out['fit_at'])
    check('fit points pulled onto the bone and kept exactly by the curve',
          np.abs(F[:, 2]).max() < 1e-9 and np.abs(P[at] - F).max() < 1e-9,
          '%.2e' % np.abs(P[at] - F).max())
    check('the periodic spline through 12 points of a circle stays on it',
          np.abs(np.hypot(P[:, 0], P[:, 1]) - 8).max() < 0.03,
          '%.4f mm' % np.abs(np.hypot(P[:, 0], P[:, 1]) - 8).max())
    th = np.linspace(0, 2 * math.pi, 24, endpoint=False)
    r = 8 + np.random.default_rng(1).normal(0, 0.4, 24)
    noisy = np.column_stack([r * np.cos(th), r * np.sin(th), np.zeros(24)])
    F = np.asarray(E.spline_answer(b, noisy, 0.25, 2.0)['fit'])
    rr = np.hypot(F[:, 0], F[:, 1])
    check('smoothing the fit points takes the noise out without shrinking the loop',
          len(F) == 24 and abs(rr.mean() - r.mean()) < 0.05 and rr.std() < 0.7 * r.std(),
          'mean r %.3f -> %.3f, sd %.3f -> %.3f' % (r.mean(), rr.mean(), r.std(), rr.std()))


def test_tubes():
    print('drill guide tubes on a flat bone')
    rng = np.random.default_rng(3)
    # the pentagon distance against brute force: distance to the five edges,
    # negative inside
    ap = 1.7
    R = ap / math.cos(math.pi / 5)
    ang = -math.pi / 2 + 2 * math.pi * np.arange(5) / 5     # a vertex at -y
    Vp = R * np.column_stack([np.cos(ang), np.sin(ang)])
    X = rng.uniform(-4, 4, (4000, 2))
    d = np.full(len(X), np.inf)
    inside = np.ones(len(X), bool)
    for k in range(5):
        a, b = Vp[k], Vp[(k + 1) % 5]
        ab = b - a
        t = np.clip(((X - a) @ ab) / (ab @ ab), 0, 1)
        d = np.minimum(d, np.linalg.norm(X - (a + t[:, None] * ab), axis=1))
        inside &= (ab[0] * (X[:, 1] - a[1]) - ab[1] * (X[:, 0] - a[0])) > 0
    ref = np.where(inside, -d, d)
    err = np.abs(E.pentagon_sdf(X[:, 0], X[:, 1], ap) - ref).max()
    check('pentagon distance is exact', err < 1e-9, 'max error %.1e' % err)

    b = plate(20.0, 0.3)
    uv = rounded_rect(24.0, 24.0, 3.0)
    Q = np.column_stack([uv, np.zeros(len(uv))])
    P, N = E.attach(b, Q, 0.1)
    base = E.Base(b, P, N, 0.2, 3.5, 1.0, 0.0, 0.1, 1.0)
    tilt = math.radians(20.0)
    # both at the outline's edge, so part of each stands off the base and its
    # bottom meets the bone clearance
    specs = [((-11.5, 0.0, 0.0), (0.0, 0.0, 1.0), 'cyl'),
             ((11.0, 0.0, 0.0), (math.sin(tilt), 0.0, math.cos(tilt)), 'pent')]
    for o, a, kind in specs:
        # the axis given pointing INTO the bone: the tube still stands outward
        tb = E.Tube(b, np.array(o) - 7.0 * np.array(a), -np.array(a), kind, 2.6, 10.0, 0.5,
                    1.0 if kind == 'cyl' else 0.0, 1.0, 0.0, [1.0, 0.0, 0.0], base.kd_in,
                    1.2)
        base.add_tube(tb)
        base.add_hole(o, a, 1.2)
    S = base.mesh(0.1)
    rep = E.tube_report(base, S)
    for r, (o, a, kind) in zip(rep, specs):
        check('%s: stands outward, entry on the bone' % kind,
              np.dot(r['axis_out'], a) > 0.999 and abs(r['entry'][2]) < 1e-6,
              'entry %s' % np.round(r['entry'], 4))
        check('%s: hole diameter 2.4' % kind, abs(2 * r['hole_min_r'] - 2.4) < 0.02,
              '%.4f' % (2 * r['hole_min_r']))
        check('%s: top 10 mm above the bone along the axis' % kind,
              abs(r['top_mm'] - 10.0) < 0.03, '%.3f' % r['top_mm'])
        check('%s: its skin 0.5 mm off the bone' % kind,
              abs(r['skin_to_bone_min_mm'] - 0.5) < 0.05, '%.3f' % r['skin_to_bone_min_mm'])
    # the cross-section half-way up: the cylinder at R, the pentagon's flats at
    # its apothem and its rounded corners at (ap - re) / cos 36 + re
    for tb, (o, a, kind) in zip(base.tubes, specs):
        Y = S.points - tb.hit
        z = Y @ tb.u
        rho = np.linalg.norm(Y - np.outer(z, tb.u), axis=1)
        ring = (np.abs(z - 7.0) < 0.3) & (rho > 1.8) & (rho < 3.2) & tb.near(S.points, 0.2)
        if kind == 'cyl':
            check('cyl: outside radius 2.6', abs(np.median(rho[ring]) - 2.6) < 0.01,
                  'median %.4f' % np.median(rho[ring]))
        else:
            ap = 2.6 * math.cos(math.pi / 5)
            corner = (ap - 1.0) / math.cos(math.pi / 5) + 1.0
            check('pent: flats at the apothem %.3f, rounded corners at %.3f' % (ap, corner),
                  abs(rho[ring].min() - ap) < 0.02 and abs(rho[ring].max() - corner) < 0.02,
                  'min %.4f max %.4f' % (rho[ring].min(), rho[ring].max()))


def test_join():
    """The shapes the curve holds, by name: a cluster of two drills closed into one
    region, a third drill 25 mm out on its own island, a fourth that stops 3 mm
    above the bone, and a body with its own margin."""
    print('section: every shape accounted for')
    b = plate()
    drill = lambda x, y, z0=-12.0, name='', **kw: dict(
        {'kind': 'cyl', 'p0': [x, y, z0], 'p1': [x, y, 12], 'r': 2.0, 'name': name}, **kw)
    shapes = [drill(-4, 0, name='a'), drill(4, 0, name='b'), drill(25, 0, name='far'),
              drill(0, -16, z0=3.0, name='short')]
    for join in (True, False):
        loops, info = E.section(b, shapes, margin=1.0, closing=5.0, planes=False,
                                join=join)
        st = {r['name']: r for r in info['shapes']}
        L = loops[0]
        check('join=%s: the cluster is in the curve' % join,
              st['a']['status'] == 'in' and st['b']['status'] == 'in', info['shapes'])
        check('join=%s: the far drill is %s' % (join, 'joined' if join else 'left out'),
              st['far']['status'] == ('joined' if join else 'out'), st['far'])
        reach = L[:, 0].max()
        check('join=%s: the curve %s round the far drill' % (join, 'goes' if join else 'does not go'),
              (reach > 27.9) if join else (reach < 15.0), '%.2f mm' % reach)
        check('join=%s: the short drill is reported 2 mm short of the bone (3 - margin 1)'
              % join, st['short']['status'] == 'off'
              and abs(st['short'].get('gap_mm', 0) - 2.0) < 0.01, st['short'])
    # per-shape margins: radius 2 drills grown by 1 and by 0.5
    two = [drill(-8, 0, name='m1'), drill(8, 0, name='m05', margin=0.5)]
    loops, info = E.section(b, two, margin=1.0, closing=0.0, join=False)
    r = sorted(float(np.median(np.hypot(np.abs(L[:, 0]) - 8, L[:, 1]))) for L in loops)
    check('each shape grows by its own margin: circles of 2.5 and 3 mm',
          len(r) == 2 and abs(r[0] - 2.5) < 0.03 and abs(r[1] - 3.0) < 0.03, r)


def main():
    test_tubes()
    test_section()
    test_join()
    m = bone()
    b = E.Bone(m.vertices, m.faces, np.array([-30, -30, -30.0]), np.array([30, 30, 30.0]),
               step=0.1)
    test_attach(m, b)
    test_base(b)
    print()
    print('%d failed: %s' % (len(FAILED), FAILED) if FAILED else 'all offline checks passed')
    return 1 if FAILED else 0


if __name__ == '__main__':
    raise SystemExit(main())
