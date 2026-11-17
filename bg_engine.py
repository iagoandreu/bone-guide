"""Bone Guide engine - curves on a bone and the guide base grown from them. py -3.10.

The bone is a triangle mesh. A curve on it is a closed polyline whose points lie
on the mesh ("attached", in 3-matic's words). The base inside a closed curve is
one field expression, contoured once:

    d(x)  signed distance to the bone: nearest dense sample, along its normal
    w(x)  signed distance to the WALL - the ruled surface of bone normals
          standing on the curve - positive outside
    base  = { c < d < c + t }  and  { w < 0 }

so the contact face is the bone offset by c, the top is the bone offset by c + t
(a true offset: a concave spot rounds off instead of folding), and the side wall
stands perpendicular to the bone all along the outline, tabs that wrap round the
bone included. The rim is rounded by a circular blend of the top and wall fields
(Quilez's round intersection): it cannot fail, and it never touches the contact
face. Units are mm throughout.
"""

import json
import math
import sys
import time

import numpy as np
import trimesh
from scipy.ndimage import gaussian_filter1d
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree


class GuideError(RuntimeError):
    pass


def _unit(V):
    n = np.linalg.norm(V, axis=-1, keepdims=True)
    return V / np.where(n > 0, n, 1.0)


# --------------------------------------------------------------------------- #
# the bone
# --------------------------------------------------------------------------- #

class Bone:
    """A bone mesh cropped to a box, with dense samples carrying interpolated
    normals: the KD tree over them is the distance field and the projector."""

    def __init__(self, V, F, lo=None, hi=None, step=0.1):
        m = trimesh.Trimesh(np.asarray(V, float), np.asarray(F, np.int64), process=True)
        if lo is not None:
            c = m.triangles_center
            keep = np.where(np.all((c >= lo) & (c <= hi), axis=1))[0]
            if len(keep) == 0:
                raise GuideError('the curve is not near the bone')
            m = m.submesh([keep], append=True)
        m.merge_vertices()
        self.mesh = m
        self.step = step
        self.P, self.N, self.fid = dense_samples(m, step)
        self.kd = cKDTree(self.P)
        e = m.edges_unique_length
        self.max_edge = float(e.max()) if len(e) else step
        # plain arrays: every read of a trimesh property re-hashes the mesh to
        # validate its cache - 22 ms (triangles) and 340 ms (vertex normals) on
        # a 1.85 M triangle scan, paid on every drag event
        self.T = np.ascontiguousarray(m.triangles)
        self.VN = np.array(m.vertex_normals)
        self.F = np.array(m.faces)
        self.locator = None

    def closest(self, X):
        """Exact closest points on the mesh, their interpolated normals and
        triangles. Candidates: the triangles of the 8 nearest samples."""
        X = np.atleast_2d(np.asarray(X, float))
        _, ii = self.kd.query(X, k=8, workers=-1)
        best_p = np.zeros_like(X)
        best_d = np.full(len(X), np.inf)
        best_f = np.zeros(len(X), np.int64)
        for j in range(ii.shape[1]):
            f = self.fid[ii[:, j]]
            p = trimesh.triangles.closest_point(self.T[f], X)
            d = np.linalg.norm(p - X, axis=1)
            better = d < best_d
            best_p[better], best_d[better], best_f[better] = p[better], d[better], f[better]
        bc = trimesh.triangles.points_to_barycentric(self.T[best_f], best_p)
        N = np.einsum('ik,ikd->id', bc, self.VN[self.F[best_f]])
        return best_p, _unit(N), best_f

    def ray(self, o, d):
        """The first point where the ray o + t d (t > 0) meets the mesh, or None.
        A VTK cell locator, built on the first call: 0.1 ms a ray."""
        import vtk
        if self.locator is None:
            import pyvista as pv
            pd = pv.PolyData(self.mesh.vertices,
                             np.hstack([np.full((len(self.F), 1), 3), self.F]).ravel())
            self.locator = vtk.vtkStaticCellLocator()
            self.locator.SetDataSet(pd)
            self.locator.BuildLocator()
            self._pd = pd
        o = np.asarray(o, float)
        d = np.asarray(d, float) / np.linalg.norm(d)
        far = 2.0 * float(np.linalg.norm(self.mesh.extents)) +             float(np.linalg.norm(o - self.mesh.centroid))
        pts, ids = vtk.vtkPoints(), vtk.vtkIdList()
        self.locator.IntersectWithLine(o, o + far * d, 1e-9, pts, ids)
        best = None
        for i in range(pts.GetNumberOfPoints()):
            q = np.array(pts.GetPoint(i))
            s = float(np.dot(q - o, d))
            if s > 0 and (best is None or s < best[0]):
                best = (s, q)
        return None if best is None else best[1]

    def signed(self, X, bound=None):
        """Signed distance (positive outside the bone) and nearest sample. With
        a bound, points further than it from every sample get +inf and index -1:
        a far query costs the KD tree most of its time."""
        if bound is None:
            _, i = self.kd.query(X, workers=-1)
            return np.einsum('ij,ij->i', X - self.P[i], self.N[i]), i
        dd, i = self.kd.query(X, distance_upper_bound=bound, workers=-1)
        ok = np.isfinite(dd)
        d = np.full(len(X), np.inf)
        d[ok] = np.einsum('ij,ij->i', X[ok] - self.P[i[ok]], self.N[i[ok]])
        return d, np.where(ok, i, -1)


def dense_samples(m, step):
    """Points on every triangle on a barycentric lattice finer than step, with
    the vertex normals interpolated - a smooth normal field, not facets."""
    V, F, Vn = m.vertices, m.faces, m.vertex_normals
    T = V[F]
    e = np.max(np.linalg.norm(T - np.roll(T, 1, axis=1), axis=2), axis=1)
    k_of = np.maximum(1, np.ceil(e / step).astype(int))
    P, N, fid = [], [], []
    for k in np.unique(k_of):
        idx = np.where(k_of == k)[0]
        ij = np.array([(i, j) for i in range(k + 1) for j in range(k + 1 - i)], float) / k
        w = np.column_stack([1.0 - ij.sum(1), ij])
        P.append(np.einsum('bk,fkd->fbd', w, T[idx]).reshape(-1, 3))
        N.append(np.einsum('bk,fkd->fbd', w, Vn[F[idx]]).reshape(-1, 3))
        fid.append(np.repeat(idx, len(w)))
    return np.vstack(P), _unit(np.vstack(N)), np.concatenate(fid)


# --------------------------------------------------------------------------- #
# curves on the bone
# --------------------------------------------------------------------------- #

def arclength(Q, closed=True):
    P = np.vstack([Q, Q[:1]]) if closed else Q
    return np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))])


def resample(Q, step, closed=True):
    s = arclength(Q, closed)
    P = np.vstack([Q, Q[:1]]) if closed else Q
    n = max(8, int(round(s[-1] / step)))
    t = np.linspace(0.0, s[-1], n, endpoint=not closed)
    return np.column_stack([np.interp(t, s, P[:, k]) for k in range(3)])


def attach(bone, Q, step=0.1, smooth_mm=0.0, closed=True, rounds=3):
    """A curve laid on the bone: resampled, optionally smoothed along its
    length, and pulled onto the mesh. Smoothing and pulling alternate, so the
    result is smooth AND on the surface."""
    Q = resample(np.asarray(Q, float), step, closed)
    P, N, _ = bone.closest(Q)
    if smooth_mm > 0:
        sig = smooth_mm / step
        mode = 'wrap' if closed else 'nearest'
        for _ in range(rounds):
            P = gaussian_filter1d(P, sig, axis=0, mode=mode)
            P, N, _ = bone.closest(P)
        P = resample(P, step, closed)
        P, N, _ = bone.closest(P)
    return P, N


def spline(bone, fit, step=0.25):
    """The closed curve through the fit points, laid on the bone: a periodic
    cubic through them (centripetal parameter, so no loops at uneven spacing),
    sampled every step and pulled onto the mesh. The fit points are kept exactly;
    fit_at[k] is the sample that is fit point k."""
    from scipy.interpolate import CubicSpline
    X = np.asarray(fit, float)
    n = len(X)
    if n < 3:
        P, N, _ = bone.closest(X)
        return P, N, np.arange(n)
    C = np.vstack([X, X[:1]])
    chord = np.linalg.norm(np.diff(C, axis=0), axis=1)
    u = np.concatenate([[0.0], np.cumsum(np.sqrt(np.maximum(chord, 1e-9)))])
    cs = CubicSpline(u, C, bc_type='periodic')
    t, fit_at = [], []
    for k in range(n):
        m = max(1, int(math.ceil(chord[k] / step)))
        fit_at.append(len(t))
        t.extend(u[k] + (u[k + 1] - u[k]) * np.arange(m) / m)
    P, N, _ = bone.closest(cs(np.asarray(t)))
    return P, N, np.asarray(fit_at)


def capped_cylinder_sdf(X, p0, p1, r):
    """Exact signed distance to the solid cylinder of radius r from p0 to p1."""
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    a = p1 - p0
    L = np.linalg.norm(a)
    a /= L
    q = X - 0.5 * (p0 + p1)
    t = q @ a
    rad = np.linalg.norm(q - np.outer(t, a), axis=1)
    dx, dy = rad - r, np.abs(t) - 0.5 * L
    out = np.hypot(np.maximum(dx, 0.0), np.maximum(dy, 0.0))
    return out + np.minimum(np.maximum(dx, dy), 0.0)


class Shapes:
    """The picked shapes as one signed distance U(x), negative inside: exact
    for cylinders; for bodies, the distance to the nearest dense sample of their
    mesh signed by its normal, made exact (closest point) near the surface."""

    def __init__(self, shapes, step=0.2):
        self.cyl, self.mesh = [], []
        for sh in shapes:
            if sh['kind'] == 'cyl':
                self.cyl.append((np.asarray(sh['p0'], float), np.asarray(sh['p1'], float),
                                 float(sh['r'])))
            else:
                V = np.fromfile(sh['coords'], '<f8').reshape(-1, 3)
                F = np.fromfile(sh['tris'], '<i4').reshape(-1, 3)
                self.mesh.append(Bone(V, F, step=step))
        self.step = step

    def points(self):
        """Points whose hull holds every shape (for bounds in any frame)."""
        P = [m.mesh.vertices for m in self.mesh]
        for p0, p1, r in self.cyl:
            c = np.array(list(np.ndindex(2, 2, 2)), float) * 2 - 1
            P.append(np.vstack([p0 + r * c, p1 + r * c]))
        return np.vstack(P)

    def sdf(self, X, reach=np.inf):
        """U at X; points further than reach from every shape may be +inf."""
        U = np.full(len(X), np.inf)
        for p0, p1, r in self.cyl:
            lo = np.minimum(p0, p1) - r - reach
            hi = np.maximum(p0, p1) + r + reach
            sel = np.where(np.all((X >= lo) & (X <= hi), axis=1))[0]
            for k0 in range(0, len(sel), 2_000_000):
                s = sel[k0:k0 + 2_000_000]
                U[s] = np.minimum(U[s], capped_cylinder_sdf(X[s], p0, p1, r))
        for m in self.mesh:
            lo, hi = m.mesh.bounds
            sel = np.where(np.all((X >= lo - reach) & (X <= hi + reach), axis=1))[0]
            dd, i = m.kd.query(X[sel], workers=-1)
            d = dd * np.where(np.einsum('ij,ij->i', X[sel] - m.P[i], m.N[i]) >= 0, 1, -1)
            near = np.abs(d) < 3 * self.step
            if near.any():
                Pc, Nc, _ = m.closest(X[sel][near])
                q = X[sel][near] - Pc
                d[near] = np.linalg.norm(q, axis=1) * np.where(
                    np.einsum('ij,ij->i', q, Nc) >= 0, 1, -1)
            U[sel] = np.minimum(U[sel], d)
        return U

    def mean_axis(self):
        """The cylinders' mean axis (signs agreed), or None without cylinders."""
        if not self.cyl:
            return None
        A = np.array([(p1 - p0) / np.linalg.norm(p1 - p0) for p0, p1, _ in self.cyl])
        A *= np.where(A @ A[0] < 0, -1, 1)[:, None]
        n = A.sum(0)
        return n / np.linalg.norm(n)


def frame_across(n):
    """Rows e1, e2, n: an orthonormal frame whose third axis is n."""
    n = np.asarray(n, float) / np.linalg.norm(n)
    ref = np.array([0.0, 0.0, 1.0]) if abs(n[2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    e1 = np.cross(ref, n)
    e1 /= np.linalg.norm(e1)
    return np.array([e1, np.cross(n, e1), n])


def iso_loops(V, F, f):
    """The zero level of a per-vertex field on a triangle mesh, as polylines:
    one point per crossed edge, one segment per crossed triangle, chained.
    Returns (closed loops, open chains)."""
    nV = len(V)
    s = f[F] < 0
    mixed = s.any(1) & ~s.all(1)
    Fm, sm = F[mixed], s[mixed]
    cross = sm[:, [0, 1, 2]] != sm[:, [1, 2, 0]]
    a = np.stack([Fm[:, 0], Fm[:, 1], Fm[:, 2]], 1)
    b = np.stack([Fm[:, 1], Fm[:, 2], Fm[:, 0]], 1)
    key = np.minimum(a, b).astype(np.int64) * nV + np.maximum(a, b)
    pairs = key[cross].reshape(-1, 2)
    uk, inv = np.unique(pairs, return_inverse=True)
    inv = inv.reshape(-1, 2)
    ea, eb = uk // nV, uk % nV
    fa, fb = f[ea], f[eb]
    t = (fa / np.where(fa - fb == 0, 1e-30, fa - fb))[:, None]
    P = V[ea] + np.clip(t, 0, 1) * (V[eb] - V[ea])
    nbr = [[] for _ in range(len(uk))]
    for i, j in inv:
        nbr[i].append(j)
        nbr[j].append(i)
    seen = np.zeros(len(uk), bool)
    loops, chains = [], []
    for start in range(len(uk)):
        if seen[start]:
            continue
        # walk to one end first, so an open chain comes out whole
        cur, prev = start, -1
        while len(nbr[cur]) == 2:
            nxt = nbr[cur][0] if nbr[cur][0] != prev else nbr[cur][1]
            prev, cur = cur, nxt
            if cur == start:
                break
        first = cur
        path, prev = [cur], -1
        seen[cur] = True
        while True:
            nx = [v for v in nbr[cur] if v != prev and not seen[v]]
            if not nx:
                break
            prev, cur = cur, nx[0]
            seen[cur] = True
            path.append(cur)
        closed = len(path) > 2 and first in nbr[path[-1]]
        (loops if closed else chains).append(P[path])
    return loops, chains


def section(bone, shapes, margin=0.0, closing=0.0, planes=True, h=None,
            max_voxels=14e6, log=None):
    """Where the bone surface meets the shapes' alpha shape - cylinders and
    bodies grown by margin, then closed with radius rho - as closed loops on the
    bone, longest first.

    U(x) = min over the shapes of their signed distance, minus margin. With
    rho = 0 the loops are U = 0 on the bone, evaluated at its vertices. With
    rho > 0, on a voxel grid around the bone near the shapes:

    planes=True  the closing is taken section by section in the planes across
                 the cylinders' mean axis - the construction of an alpha shape
                 lofted through closed plane sections (Alpha Shape, 3-matic):
                 side tabs that wrap round the bone edge stay one region.
                 Per-section distance transforms; error about one voxel.
    planes=False a 3D closing by a ball: {x : dist(x, outside {U <= rho}) >=
                 rho}, the erosion read sub-voxel from the outer layer's own
                 exact U.

    The closing at a point depends on nothing further than 2 rho, so the grid
    is the bone's neighbourhood of the shapes plus 2 rho."""
    log = log or (lambda *_: None)
    from scipy import ndimage
    S = Shapes(shapes)
    V, F = bone.mesh.vertices, bone.mesh.faces
    info = {'closing': 'none'}
    if closing <= 0:
        f = np.full(len(V), 1e3)
        U = S.sdf(V, reach=margin + 1.0) - margin
        f[np.isfinite(U)] = U[np.isfinite(U)]
    else:
        n = S.mean_axis() if planes else None
        R = frame_across(n) if n is not None else np.eye(3)
        info['closing'] = 'planes across the drills' if n is not None else '3D'
        pad = margin + 2.0 * closing + 1.0
        Ls, Lb = S.points() @ R.T, V @ R.T
        lo = np.maximum(Ls.min(0) - pad, Lb.min(0) - pad)
        hi = np.minimum(Ls.max(0) + pad, Lb.max(0) + pad)
        if np.any(hi <= lo):
            raise GuideError('the shapes do not reach the bone')
        if h is None:
            h = max(0.12, (np.prod(hi - lo) / max_voxels) ** (1.0 / 3.0))
        dims = np.ceil((hi - lo) / h).astype(int) + 1
        axes = [lo[k] + h * np.arange(dims[k]) for k in range(3)]
        L = np.stack(np.meshgrid(*axes, indexing='ij'), -1).reshape(-1, 3)
        U = (S.sdf(L @ R, reach=margin + closing + 2 * h) - margin).reshape(tuple(dims))
        del L
        log('shape field %d x %d x %d voxels of %.3f mm' % (dims[0], dims[1], dims[2], h))
        g = (Lb - lo) / h
        inside = np.all((g >= 0) & (g <= dims - 1), axis=1)
        if n is not None:
            # sections across n: a huge spacing along n keeps every distance
            # inside its own section; the section's rim counts as outside
            samp = (h, h, 1e7)
            dout = ndimage.distance_transform_edt(~(U <= 0), sampling=samp)
            D = dout <= closing
            Dp = np.pad(D, ((1, 1), (1, 1), (0, 0)))
            ein = ndimage.distance_transform_edt(Dp, sampling=samp)[1:-1, 1:-1, :]
            C = np.where(D, closing - ein + 0.5 * h, closing + h)
        else:
            # erosion of D = {U <= rho} at the voxels the bone samples: the
            # distance from x to the outside of D is min over outer-layer voxels
            # y of |x - y| - (U(y) - rho), sub-voxel since U is exact outside D
            U = np.where(np.isfinite(U), U, 1e3)
            D = U <= closing
            Dp = np.pad(D, 1)
            ring = ndimage.binary_dilation(Dp) & ~Dp
            iy = np.argwhere(ring) - 1
            Up = np.pad(U, 1, constant_values=closing)
            uy = Up[tuple((iy + 1).T)]
            kd = cKDTree(lo + iy * h)
            c0 = np.floor(g[inside]).astype(int)
            corners = np.unique(np.vstack([c0 + np.array(o) for o in np.ndindex(2, 2, 2)]),
                                axis=0)
            corners = corners[np.all((corners >= 0) & (corners < dims), axis=1)]
            cin = corners[D[tuple(corners.T)]]
            dd, ii = kd.query(lo + cin * h, k=8, workers=-1)
            e = np.maximum((dd - (uy[ii] - closing)).min(1), 0.0)
            C = U.copy()
            C[tuple(cin.T)] = closing - e
        f = np.full(len(V), 1e3)
        f[inside] = ndimage.map_coordinates(C, g[inside].T, order=1)
        info.update({'voxel_mm': float(h), 'grid': dims.tolist()})
    loops, chains = iso_loops(V, F, f)
    loops.sort(key=lambda L: -arclength(L)[-1])
    info.update({'loops': len(loops), 'open_chains': len(chains),
                 'loop_mm': [float(arclength(L)[-1]) for L in loops[:8]]})
    return loops, info


def frames(Q, N, step, closed=True, normal_mm=1.5):
    """Tangent, bone normal smoothed along the curve, and binormal T x N (in
    the bone's tangent plane, across the curve)."""
    mode = 'wrap' if closed else 'nearest'
    T = _unit(np.roll(Q, -1, 0) - np.roll(Q, 1, 0)) if closed else \
        _unit(np.gradient(Q, axis=0))
    Ns = _unit(gaussian_filter1d(N, max(normal_mm / step, 0.5), axis=0, mode=mode))
    Ns = _unit(Ns - np.einsum('ij,ij->i', Ns, T)[:, None] * T)
    return T, Ns, np.cross(T, Ns)


def inside_labels(bone, Q, B):
    """Which mesh vertices lie inside the closed curve. Edges the curve crosses
    are cut (endpoints on opposite sides of the local wall), the rest of the
    mesh falls apart into components, and the components on the side away from
    the crop boundary are inside. Returns (vertex labels, outward sign of B)."""
    m = bone.mesh
    V = m.vertices
    delta = min(max(2.0 * bone.max_edge, 1.0), 6.0)
    dq, iq = cKDTree(Q).query(V, distance_upper_bound=delta, workers=-1)
    near = np.isfinite(dq)
    side = np.zeros(len(V), np.int8)
    s = np.einsum('ij,ij->i', V[near] - Q[iq[near]], B[iq[near]])
    side[near] = np.where(s >= 0, 1, -1)
    E = m.edges_unique
    cut = near[E[:, 0]] & near[E[:, 1]] & (side[E[:, 0]] != side[E[:, 1]])
    K = E[~cut]
    g = coo_matrix((np.ones(len(K)), (K[:, 0], K[:, 1])), shape=(len(V), len(V)))
    nc, comp = connected_components(g, directed=False)
    vote = np.bincount(comp, weights=side.astype(float), minlength=nc)
    cside = np.where(vote >= 0, 1, -1)
    on_border = np.zeros(nc, bool)
    be = trimesh.grouping.group_rows(m.edges_sorted, require_count=1)
    if len(be):
        on_border[np.unique(comp[m.edges_sorted[be].ravel()])] = True
    if on_border.any():
        out = 1 if vote[on_border & (cside == 1)].size and \
            np.abs(vote[on_border & (cside == 1)]).sum() >= np.abs(vote[on_border & (cside == -1)]).sum() \
            else -1
    else:
        area = np.bincount(comp[m.faces[:, 0]], weights=m.area_faces, minlength=nc)
        out = 1 if area[cside == 1].sum() >= area[cside == -1].sum() else -1
    inside = cside[comp] != out
    if not inside.any():
        raise GuideError('nothing lies inside the curve')
    if on_border[np.unique(comp[inside])].any():
        raise GuideError('the region inside the curve runs off the cropped bone - '
                         'is the curve closed and on the bone?')
    return inside, out


# --------------------------------------------------------------------------- #
# the base
# --------------------------------------------------------------------------- #

def round_intersection(a, b, r):
    """max(a, b) with the convex edge between the two faces rounded to radius r.
    Exact for distance fields meeting at right angles."""
    if r <= 0:
        return np.maximum(a, b)
    ua, ub = np.maximum(a + r, 0.0), np.maximum(b + r, 0.0)
    return np.sqrt(ua * ua + ub * ub) + np.minimum(np.maximum(a + r, b + r), 0.0) - r


def line_distance(X, origin, axis):
    q = X - origin
    return np.linalg.norm(q - np.outer(q @ axis, axis), axis=1)


# --------------------------------------------------------------------------- #
# drill guides (tubes) standing on the bone
# --------------------------------------------------------------------------- #

_K5 = (0.809016994374947, 0.587785252292473, 0.726542528005361)


def pentagon_sdf(px, py, apothem):
    """Exact 2D signed distance to a regular pentagon with this apothem, a flat
    facing +y and a vertex at -y (Inigo Quilez' construction, vectorised)."""
    kx, ky, kz = _K5
    x = np.abs(px)
    y = np.array(py, float, copy=True)
    for sx in (-kx, kx):
        dot = np.minimum(sx * x + ky * y, 0.0)
        x = x - 2.0 * dot * sx
        y = y - 2.0 * dot * ky
    cx = np.clip(x, -apothem * kz, apothem * kz)
    qx, qy = x - cx, y - apothem
    return np.hypot(qx, qy) * np.sign(qy)


class Tube:
    """A drill guide around one hole axis: a cylinder or a pentagonal prism
    (circumscribed radius R, edges rounded edge_round), standing on the bone at
    the axis' entry point, its top `height` above the bone along the axis and
    rounded top_round, its bottom held `clear` off the bone. The hole itself is
    cut by Base like any other."""

    def __init__(self, bone, origin, axis, kind='cyl', R=2.6, height=10.0, clear=0.5,
                 top_round=1.0, edge_round=1.0, spin_deg=0.0, ref=None, near=None,
                 hole=0.0):
        a = _unit(np.asarray(axis, float))
        o = np.asarray(origin, float)
        span = 2.0 * float(np.linalg.norm(bone.mesh.extents)) + 50.0
        hits = [h for h in (bone.ray(o - span * a, a), bone.ray(o + span * a, -a))
                if h is not None]
        if not hits:
            raise GuideError('a guide axis does not meet the bone')
        if near is not None:
            dn = [float(near.query(h)[0]) for h in hits]
            hit = hits[int(np.argmin(dn))]
        else:
            hit = hits[0]
        _, Nh, _ = bone.closest(hit[None])
        self.u = a if np.dot(a, Nh[0]) > 0 else -a
        tilt = min(math.acos(min(1.0, abs(float(np.dot(a, Nh[0]))))), math.radians(75.0))
        self.hit = hit
        ref = np.asarray(ref if ref is not None else [0.0, 0.0, 1.0], float)
        e1 = ref - np.dot(ref, self.u) * self.u
        if np.linalg.norm(e1) < 1e-6:
            e1 = np.cross(self.u, [1.0, 0.0, 0.0] if abs(self.u[0]) < 0.9 else [0.0, 1.0, 0.0])
        e1 = _unit(e1)
        e2 = np.cross(self.u, e1)
        s = math.radians(spin_deg)
        self.e1, self.e2 = math.cos(s) * e1 + math.sin(s) * e2, -math.sin(s) * e1 + math.cos(s) * e2
        self.kind, self.R, self.height, self.clear = kind, float(R), float(height), float(clear)
        self.rt = min(float(top_round), 0.999 * min(self.R, self.height))
        self.re = min(float(edge_round), 0.999 * self.R * math.cos(math.pi / 5)) \
            if kind == 'pent' else 0.0
        self.bone = bone
        self.hole = float(hole)
        # how far below the entry point the tube reaches before the bone stops
        # it (a tilted axis, plus room for the bone's own curvature); above zcut
        # nothing can be within clear of the bone
        self.below = self.R * math.tan(tilt) + 2.0
        self.zcut = self.below + self.clear + 1.0

    def section_sdf(self, x, y):
        if self.kind == 'pent':
            ap = self.R * math.cos(math.pi / 5)
            # a vertex along e1: the construction has its vertex at -y
            return pentagon_sdf(-y, -x, ap - self.re) - self.re
        return np.hypot(x, y) - self.R

    def sdf(self, X):
        Y = X - self.hit
        z = Y @ self.u
        s2 = self.section_sdf(Y @ self.e1, Y @ self.e2)
        # the section extruded up to the top, its top edge rounded by rt
        wx, wy = s2 + self.rt, z - self.height + self.rt
        f = np.minimum(np.maximum(wx, wy), 0.0) + np.hypot(np.maximum(wx, 0.0),
                                                           np.maximum(wy, 0.0)) - self.rt
        f = np.maximum(f, -(z + self.below))
        low = z < self.zcut
        if low.any():
            d, _ = self.bone.signed(X[low], bound=self.below + self.R + self.clear + 1.0)
            f[low] = np.maximum(f[low], self.clear - d)
        return f

    def segment(self):
        return self.hit - self.below * self.u, self.hit + self.height * self.u

    def near(self, X, pad):
        a, b = self.segment()
        ab = b - a
        t = np.clip((X - a) @ ab / (ab @ ab), 0.0, 1.0)
        return np.linalg.norm(X - (a + t[:, None] * ab), axis=1) <= self.R + pad


class Base:
    """The guide base inside a closed curve on the bone."""

    def __init__(self, bone, Q, N, c=0.2, t=3.5, r=1.0, r_contact=0.0, step=0.1,
                 normal_mm=1.5, log=None):
        self.log = log or (lambda *_: None)
        self.bone, self.c, self.t, self.r, self.rc = bone, c, t, r, r_contact
        T, Ns, B = frames(Q, N, step, True, normal_mm)
        inside_v, out = inside_labels(bone, Q, B)
        self.B = B * out                       # now +B points out of the region
        self.Q, self.Ns = Q, Ns
        m = bone.mesh
        lab = inside_v[m.faces].sum(1) >= 2   # a face is inside with 2 of 3 corners
        self.inside_s = lab[bone.fid]
        # near the curve the face label is too coarse: decide by the local wall
        dq, iq = cKDTree(Q).query(bone.P, distance_upper_bound=2.0 * bone.max_edge,
                                  workers=-1)
        nr = np.isfinite(dq)
        self.inside_s[nr] = np.einsum('ij,ij->i', bone.P[nr] - Q[iq[nr]],
                                      self.B[iq[nr]]) < 0
        # the wall: bone normals standing on the curve, sampled as a ruled band
        lam = np.arange(-1.0, c + t + r + 1.0 + 1e-9, step)
        self.R = (Q[:, None, :] + lam[None, :, None] * Ns[:, None, :]).reshape(-1, 3)
        self.RB = np.repeat(self.B, len(lam), axis=0)
        self.kd_rib = cKDTree(self.R)
        self.kd_in = cKDTree(bone.P[self.inside_s])
        self.area_in = float(m.area_faces[lab].sum())
        self.holes = []
        self.tubes = []
        self.dropped = None

    def add_hole(self, origin, axis, radius):
        self.holes.append((np.asarray(origin, float), _unit(np.asarray(axis, float)),
                           float(radius)))

    def add_tube(self, tube):
        self.tubes.append(tube)

    def field(self, X):
        c, t, r = self.c, self.t, self.r
        F = np.full(len(X), 2.0, float)
        dd, _ = self.kd_in.query(X, distance_upper_bound=c + t + 0.5, workers=-1)
        cand = np.isfinite(dd)
        Xc = X[cand]
        if len(Xc):
            d, i = self.bone.signed(Xc)
            inside = self.inside_s[i]
            dr, ir = self.kd_rib.query(Xc, workers=-1)
            w = np.einsum('ij,ij->i', Xc - self.R[ir], self.RB[ir])
            far = dr > max(2.0 * r, 1.5)
            w[far] = np.where(inside[far], -dr[far], dr[far])
            top = round_intersection(d - (c + t), w, r)
            if self.rc > 0:
                F[cand] = round_intersection(c - d, top, self.rc)
            else:
                F[cand] = np.maximum(c - d, top)
        for tb in self.tubes:
            sel = tb.near(X, 0.5)
            if sel.any():
                F[sel] = np.minimum(F[sel], tb.sdf(X[sel]))
        for o, a, rad in self.holes:
            F = np.maximum(F, rad - line_distance(X, o, a))
        return F

    def mesh(self, h=0.15, block=8):
        """Contour the field on a voxel grid of pitch h. The grid is cut into
        blocks of block^3 voxels and only blocks the base can reach are
        evaluated: most of the box is bone or air."""
        import pyvista as pv
        Pin = self.bone.P[self.inside_s]
        pad = self.c + self.t + 3 * h
        lo = Pin.min(0) - pad
        hi = Pin.max(0) + pad
        for tb in self.tubes:
            ends = np.array(tb.segment())
            lo = np.minimum(lo, ends.min(0) - tb.R - 3 * h)
            hi = np.maximum(hi, ends.max(0) + tb.R + 3 * h)
        nb = np.ceil((hi - lo) / (h * block)).astype(int)
        dims = nb * block + 1
        F = np.full(tuple(dims), 2.0, np.float32)
        bs = h * block
        g = np.stack(np.meshgrid(*[lo[k] + bs * (np.arange(nb[k]) + 0.5) for k in range(3)],
                                 indexing='ij'), -1).reshape(-1, 3)
        rad = 0.5 * math.sqrt(3.0) * bs + h
        dd, _ = self.kd_in.query(g, distance_upper_bound=self.c + self.t + 0.5 + rad,
                                 workers=-1)
        live = np.isfinite(dd)
        db, _ = self.bone.signed(g[live])
        ok = (db > self.c - rad) & (db < self.c + self.t + rad)
        keep = np.zeros(len(g), bool)
        keep[np.where(live)[0][ok]] = True
        for tb in self.tubes:
            keep |= tb.near(g, rad + 0.5)
        idx = np.argwhere(keep.reshape(tuple(nb)))
        o = np.arange(block + 1) * h
        cell = np.stack(np.meshgrid(o, o, o, indexing='ij'), -1).reshape(-1, 3)
        chunk = max(1, 2_000_000 // len(cell))
        for k0 in range(0, len(idx), chunk):
            ib = idx[k0:k0 + chunk]
            X = (lo + ib[:, None, :] * bs + cell[None]).reshape(-1, 3)
            f = self.field(X).reshape(len(ib), block + 1, block + 1, block + 1)
            for j, (a, b, c) in enumerate(ib):
                sl = (slice(a * block, a * block + block + 1),
                      slice(b * block, b * block + block + 1),
                      slice(c * block, c * block + block + 1))
                F[sl] = f[j]
        self.log('field %d x %d x %d voxels, %d of %d blocks evaluated'
                 % (dims[0], dims[1], dims[2], len(idx), int(np.prod(nb))))
        img = pv.ImageData(dimensions=tuple(int(v) for v in dims), spacing=(h, h, h),
                           origin=tuple(lo))
        img.point_data['f'] = F.ravel(order='F')
        surf = img.contour([0.0], scalars='f', method='flying_edges')
        surf = surf.clean().triangulate()
        if surf.n_points == 0:
            raise GuideError('the base came out empty')
        # only the main piece: sub-millimetre specks and voids (a distal radius
        # test case with 12 guides: 24 of them, 0.003 mm3 in all) would each become a body of
        # their own in Convert Mesh
        tri = surf.faces.reshape(-1, 4)[:, 1:]
        parts = trimesh.Trimesh(surf.points, tri, process=False).split(only_watertight=False)
        if len(parts) > 1:
            parts = sorted(parts, key=lambda p: -abs(p.volume))
            self.dropped = {'pieces': len(parts) - 1,
                            'volume_mm3': float(sum(abs(p.volume) for p in parts[1:])),
                            'largest_extent_mm': float(max(p.extents.max() for p in parts[1:]))}
            self.log('dropped %d loose pieces, %.4f mm3' % (self.dropped['pieces'],
                                                            self.dropped['volume_mm3']))
            keep = parts[0]
            surf = pv.PolyData(keep.vertices, np.column_stack(
                [np.full(len(keep.faces), 3), keep.faces]).ravel())
            tri = keep.faces
        # flying edges winds its triangles with the field's gradient, which points
        # out of the solid here (f < 0 inside): the faces come out inside-out
        if trimesh.Trimesh(surf.points, tri, process=False).volume < 0:
            surf.faces = np.column_stack([np.full(len(tri), 3), tri[:, ::-1]]).ravel()
        return surf


def largest_loop(mesh_pv):
    """The longest boundary loop of an open pyvista surface, ordered."""
    fe = mesh_pv.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                       manifold_edges=False, non_manifold_edges=False)
    fe = fe.clean()
    lines = fe.lines.reshape(-1, 3)[:, 1:]
    nbr = {}
    for a, b in lines:
        nbr.setdefault(a, []).append(b)
        nbr.setdefault(b, []).append(a)
    seen, loops = set(), []
    for s in nbr:
        if s in seen:
            continue
        loop, prev, cur = [s], None, s
        seen.add(s)
        while True:
            nx = [v for v in nbr[cur] if v != prev]
            if not nx or nx[0] == s or nx[0] in seen:
                break
            prev, cur = cur, nx[0]
            loop.append(cur)
            seen.add(cur)
        loops.append(fe.points[loop])
    return max(loops, key=lambda L: arclength(L)[-1])


# --------------------------------------------------------------------------- #
# job interface (the Fusion add-in writes a JSON job and runs this file)
# --------------------------------------------------------------------------- #

def read_bin_mesh(job):
    V = np.fromfile(job['bone_coords'], '<f8').reshape(-1, 3)
    F = np.fromfile(job['bone_tris'], '<i4').reshape(-1, 3)
    return V, F


def write_stl(surf, path):
    tri = surf.faces.reshape(-1, 4)[:, 1:]
    trimesh.Trimesh(surf.points, tri, process=False).export(path)


def run(job, log=None):
    log = log or (lambda *_: None)
    t0 = time.time()
    V, F = (read_bin_mesh(job) if 'bone_coords' in job else
            (lambda m: (m.vertices, m.faces))(trimesh.load(job['bone_stl'], force='mesh')))
    Q = np.asarray(job['curve'], float)
    c, t, r = float(job.get('clearance', 0.2)), float(job.get('thickness', 3.5)), \
        float(job.get('round', 1.0))
    pad = c + t + r + 12.0
    bone = Bone(V, F, Q.min(0) - pad, Q.max(0) + pad, step=float(job.get('sample', 0.1)))
    log('bone cropped: %d triangles, %d samples' % (len(bone.mesh.faces), len(bone.P)))
    if job.get('mode') == 'attach':
        P, N = attach(bone, Q, float(job.get('curve_step', 0.25)),
                      float(job.get('smooth', 0.0)))
        return {'curve': P.tolist(), 'normals': N.tolist(), 'secs': time.time() - t0}
    P, N = attach(bone, Q, 0.1, float(job.get('smooth', 0.0)))
    base = Base(bone, P, N, c, t, r, float(job.get('round_contact', 0.0)), 0.1,
                float(job.get('normal_mm', 1.5)), log)
    for hle in job.get('holes', []):
        base.add_hole(hle['origin'], hle['axis'], hle['radius'])
    # the guide's own length direction: the reference a pentagon's vertex points along
    w, v = np.linalg.eigh(np.cov((P - P.mean(0)).T))
    ref = v[:, -1]
    for tj in job.get('tubes', []):
        base.add_tube(Tube(bone, tj['origin'], tj['axis'], tj.get('kind', 'cyl'),
                           float(tj.get('R', 2.6)), float(tj.get('height', 10.0)),
                           float(tj.get('clear', 0.5)), float(tj.get('top_round', 1.0)),
                           float(tj.get('edge_round', 1.0)), float(tj.get('spin', 0.0)),
                           ref, base.kd_in, float(tj['radius'])))
        base.add_hole(tj['origin'], tj['axis'], tj['radius'])
    surf = base.mesh(float(job.get('voxel', 0.15)))
    dec = None
    if float(job.get('decimate_tol', 0.0)) > 0:
        surf, dec = decimate(surf, float(job['decimate_tol']), log)
    write_stl(surf, job['out_stl'])
    d, _ = bone.signed(surf.points)
    contact = d[d < c + 0.05]
    vol = trimesh.Trimesh(surf.points, surf.faces.reshape(-1, 4)[:, 1:], process=False).volume
    near_tube = np.zeros(surf.n_points, bool)
    for tb in base.tubes:
        near_tube |= tb.near(surf.points, 0.2)
    out = {'triangles': int(surf.n_cells), 'volume_mm3': float(vol),
           'footprint_mm2': base.area_in, 'contact_median_mm': float(np.median(contact))
           if len(contact) else None,
           'top_max_mm': float(d[~near_tube].max()) if (~near_tube).any() else float(d.max()),
           'curve_mm': float(arclength(P)[-1]), 'tubes': tube_report(base, surf),
           'decimation': dec, 'dropped': base.dropped, 'secs': time.time() - t0}
    log('done in %.1f s' % out['secs'])
    return out


def tube_report(base, surf):
    """Per tube, measured on the mesh that is written: the narrowest the hole
    gets along the tube (what the drill sees), and on the tube's own skin (the
    vertices where its field is zero) the top above the bone along the axis and
    the closest that skin comes to the bone."""
    X = surf.points
    out = []
    for tb in base.tubes:
        sel = np.where(tb.near(X, 0.2))[0]
        Y = X[sel] - tb.hit
        z = Y @ tb.u
        rho = np.linalg.norm(Y - np.outer(z, tb.u), axis=1)
        bore = (z > 0.3) & (z < tb.height - tb.rt - 0.1) & (rho < 0.5 * tb.R + 0.5)
        # the tube's own outside (the bore wall also sits on its zero level
        # near the bottom, where the bone clearance cuts it)
        skin = (np.abs(tb.sdf(X[sel])) < 0.03) & (rho > tb.hole + 0.3)
        d, _ = base.bone.signed(X[sel][skin]) if skin.any() else (np.array([np.nan]), None)
        out.append({'kind': tb.kind, 'entry': tb.hit.tolist(), 'axis_out': tb.u.tolist(),
                    'hole_min_r': float(rho[bore].min()) if bore.any() else None,
                    'top_mm': float(z[skin].max()) if skin.any() else None,
                    'skin_to_bone_min_mm': float(np.nanmin(d))})
    return out


def decimate(surf, tol, log=None):
    """Fewer triangles for a faceted solid, the surface kept within tol (mm) of
    the contoured one both ways (Hausdorff, measured). Quadric edge collapse in
    steps from the full mesh; the last step within tol is kept. On a distal radius
    test guide with 12 tubes: 432 k -> 65 k triangles at 0.038 mm."""
    import pymeshlab
    import pyvista as pv
    log = log or (lambda *_: None)
    V0 = np.asarray(surf.points, float)
    F0 = surf.faces.reshape(-1, 4)[:, 1:]
    ref = pymeshlab.Mesh(V0, F0)
    best, best_dev, n = (V0, F0), 0.0, len(F0)
    for frac in (0.5, 0.35, 0.25, 0.2, 0.15, 0.12, 0.1, 0.08, 0.06):
        ms = pymeshlab.MeshSet()
        ms.add_mesh(pymeshlab.Mesh(V0, F0))
        ms.meshing_decimation_quadric_edge_collapse(
            targetfacenum=int(n * frac), qualitythr=0.4, preserveboundary=True,
            preservenormal=True, preservetopology=True, optimalplacement=True,
            planarquadric=True)
        m = ms.current_mesh()
        V, F = m.vertex_matrix(), m.face_matrix()
        ms.add_mesh(ref)
        a = ms.get_hausdorff_distance(sampledmesh=0, targetmesh=1, samplevert=True,
                                      sampleface=True, samplenum=max(20000, 3 * len(F)))
        b = ms.get_hausdorff_distance(sampledmesh=1, targetmesh=0, samplevert=True,
                                      sampleface=True, samplenum=max(20000, len(F0)))
        dev = max(a['max'], b['max'])
        tm = trimesh.Trimesh(V, F, process=False)
        log('decimate %.4f: %d triangles, %.4f mm' % (frac, len(F), dev))
        if dev > tol or not tm.is_watertight:
            break
        best, best_dev = (V, F), dev
    V, F = best
    out = pv.PolyData(V, np.column_stack([np.full(len(F), 3), F]).ravel())
    return out, {'triangles_before': int(n), 'triangles': int(len(F)),
                 'hausdorff_mm': float(best_dev), 'tol_mm': float(tol)}


def spline_answer(bone, fit, step, smooth=0.0):
    """Fit points pulled onto the bone, and the spline through them. smooth
    (mm along the curve) first smooths the fit points themselves - their number
    stays, so the spline can still be written back in place."""
    F = np.asarray(fit, float)
    if smooth > 0 and len(F) >= 4:
        # Taubin lambda / mu steps on the closed polygon: noise goes, the loop
        # does not shrink (a Gaussian on a few points pulls it inward)
        gap = float(np.mean(np.linalg.norm(np.roll(F, -1, 0) - F, axis=1)))
        n = int(min(200, max(1, round(2.0 * (smooth / max(gap, 1e-6)) ** 2))))
        for _ in range(n):
            for k in (0.5, -0.53):
                F = F + k * (0.5 * (np.roll(F, 1, 0) + np.roll(F, -1, 0)) - F)
            F, _, _ = bone.closest(F)
    F, _, _ = bone.closest(F)
    P, N, at = spline(bone, F, step)
    return {'fit': F.tolist(), 'pts': P.tolist(), 'nrm': N.tolist(),
            'fit_at': at.tolist(), 'length_mm': float(arclength(P)[-1])}


def serve():
    """Resident mode for the interactive curve editor: one JSON request per line
    on stdin, one JSON answer per line on stdout. The bone is loaded once; an
    attach of a few hundred points then takes milliseconds."""
    bone = None
    for line in sys.stdin:
        try:
            req = json.loads(line)
            op = req.get('op')
            if op == 'quit':
                break
            if op == 'load':
                V = np.fromfile(req['coords'], '<f8').reshape(-1, 3)
                F = np.fromfile(req['tris'], '<i4').reshape(-1, 3)
                bone = Bone(V, F, step=float(req.get('step', 0.25)))
                bone.ray(bone.mesh.centroid, (0.0, 0.0, 1.0))   # builds the locator now
                out = {'triangles': int(len(bone.mesh.faces)), 'samples': int(len(bone.P))}
            elif bone is None:
                out = {'error': 'no bone loaded'}
            elif op == 'ray':
                q = bone.ray(req['o'], req['d'])
                out = {'hit': None if q is None else q.tolist()}
            elif op == 'attach':
                P, N, _ = bone.closest(np.asarray(req['pts'], float))
                out = {'pts': P.tolist(), 'nrm': N.tolist()}
            elif op == 'curve':
                closed = bool(req.get('closed', True))
                P, N = attach(bone, np.asarray(req['pts'], float),
                              float(req.get('step', 0.25)), float(req.get('smooth', 0.0)),
                              closed)
                out = {'pts': P.tolist(), 'nrm': N.tolist(),
                       'length_mm': float(arclength(P, closed)[-1])}
            elif op == 'spline':
                out = spline_answer(bone, req['fit'], float(req.get('step', 0.25)),
                                    float(req.get('smooth', 0.0)))
            elif op == 'fit':
                # a curve given as points: smoothed onto the bone, fit points
                # every `spacing` mm along it, and the spline through them
                P, _ = attach(bone, np.asarray(req['pts'], float), 0.1,
                              float(req.get('smooth', 0.0)))
                fit = resample(P, float(req.get('spacing', 3.0)))
                out = spline_answer(bone, fit, float(req.get('step', 0.25)))
            elif op == 'section':
                t0 = time.time()
                loops, info = section(bone, req['shapes'], float(req.get('margin', 0.0)),
                                      float(req.get('closing', 0.0)),
                                      bool(req.get('planes', True)), max_voxels=8e6)
                if not loops:
                    raise GuideError('the shapes do not cut the bone surface')
                P, _ = attach(bone, loops[0], 0.1, float(req.get('smooth', 0.5)))
                fit = resample(P, float(req.get('spacing', 3.0)))
                out = spline_answer(bone, fit, float(req.get('step', 0.25)))
                info['secs'] = time.time() - t0
                out['info'] = info
            else:
                out = {'error': 'unknown op %r' % op}
        except GuideError as exc:
            out = {'error': str(exc), 'user_error': True}
        except Exception as exc:
            out = {'error': '%s: %s' % (type(exc).__name__, exc)}
        sys.stdout.write(json.dumps(out) + '\n')
        sys.stdout.flush()


def main(argv):
    if argv[1] == '--serve':
        return serve()
    job = json.load(open(argv[1], encoding='utf-8'))
    logp = argv[2] + '.log'

    def log(msg):
        with open(logp, 'a', encoding='utf-8') as fh:
            fh.write(msg + '\n')
    try:
        out = {'ok': True, 'result': run(job, log)}
    except GuideError as exc:
        out = {'ok': False, 'user_error': True, 'error': str(exc)}
    except Exception:
        import traceback
        out = {'ok': False, 'error': traceback.format_exc()}
    with open(argv[2], 'w', encoding='utf-8') as fh:
        json.dump(out, fh)


if __name__ == '__main__':
    main(sys.argv)
