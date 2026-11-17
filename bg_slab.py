"""The guide base as a NURBS solid, fitted in numpy - no mesh in the result.

Fusion splitting, offsetting and thickening a scanned bone of 20 000 faces is the
slow part of any B-Rep guide. Here Fusion never touches the bone: its triangles
are sampled once, and the base comes back as six NURBS faces that Fusion only
assembles (BRepBodyDefinition) and rounds.

    bottom   the bone inside the outline, offset by the clearance along its normals
    top      the same, offset by clearance + thickness
    wall     four ruled faces between the bottom and top edges - the wall stands
             along the bone normal all round, as the outline is projected

The patch inside the outline is a disk; it is mapped onto the unit square with a
Tutte embedding (uniform weights: a convex boundary guarantees no fold), the
outline split into four arcs at quarter lengths. Every edge is fitted ONCE and
shared by the two faces on it, so the shell is watertight by construction, and
the faces are laid out so u x v points outward and every loop is counter-
clockwise - no face needs isParamReversed.

    8 vertices  b0..b3 bottom, t0..t3 top - corners (0,0) (1,0) (1,1) (0,1)
    12 edges    eb0..eb3, et0..et3 along u or v (increasing), ev0..ev3 b -> t
    6 faces     top (u, v), bottom (v, u), side0..side3 ruled
"""

import math

import numpy as np
import trimesh
from scipy.interpolate import BSpline
from scipy.sparse import coo_matrix, csr_matrix, diags, vstack as svstack
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import lsqr, spsolve
from scipy.spatial import cKDTree

import bg_engine as E


# how far the cutting band reaches below the contact face: room for the round there
LO_BELOW = 2.0


class SlabError(ValueError):
    pass


def clamped_knots(n, deg=3):
    inner = n - deg - 1
    return np.r_[np.zeros(deg + 1), np.arange(1, inner + 1) / (inner + 1.0), np.ones(deg + 1)]


def basis(x, knots, deg=3):
    x = np.clip(np.asarray(x, float), 0.0, 1.0 - 1e-12)
    return BSpline.design_matrix(x, knots, deg).tocsr()


def second_diff(n):
    if n < 3:
        return csr_matrix((0, n))
    return diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2],
                 shape=(n - 2, n)).tocsr()


# --------------------------------------------------------------------------- #
# the patch and its parametrisation
# --------------------------------------------------------------------------- #

def region(bone, P, N):
    """Faces of the bone inside the outline P (on the bone, CCW about N),
    the largest piece, and its boundary loop. Returns (V, F, loop vertex ids)."""
    T, Ns, B = E.frames(P, N, 0.25)
    inside_v, out = E.inside_labels(bone, P, B)
    m = bone.mesh
    keep = np.where(inside_v[m.faces].all(1))[0]
    if len(keep) == 0:
        raise SlabError('nothing of the bone lies inside the outline')
    sub = m.submesh([keep], append=True)
    sub.merge_vertices()
    # the largest connected piece
    g = coo_matrix((np.ones(len(sub.edges_unique)), (sub.edges_unique[:, 0],
                                                      sub.edges_unique[:, 1])),
                   shape=(len(sub.vertices),) * 2)
    _n, lab = connected_components(g, directed=False)
    big = np.argmax(np.bincount(lab))
    fk = np.where(lab[sub.faces[:, 0]] == big)[0]
    sub = sub.submesh([fk], append=True)
    sub.merge_vertices()
    be = sub.edges[trimesh.grouping.group_rows(sub.edges_sorted, require_count=1)]
    loops = trimesh.graph.connected_components(be)
    if not loops:
        raise SlabError('the patch inside the outline has no boundary')
    # the outer loop is the longest; any other is a hole in the patch
    def loop_order(vs):
        vs = set(vs)
        nb = {}
        for a, b in be:
            if a in vs:
                nb.setdefault(a, []).append(b)
                nb.setdefault(b, []).append(a)
        start = next(iter(vs))
        seq, prev = [start], None
        while True:
            nx = [w for w in nb[seq[-1]] if w != prev]
            if not nx or nx[0] == start:
                break
            prev = seq[-1]
            seq.append(nx[0])
            if len(seq) > len(vs):
                break
        return seq
    loops = sorted((loop_order(L) for L in loops), key=len, reverse=True)
    return np.asarray(sub.vertices), np.asarray(sub.faces), loops[0], len(loops) - 1


def corners_of(P):
    """Four outline indices at quarter lengths, the first at the outline's far
    end along its long axis - the square's corners (0,0) (1,0) (1,1) (0,1)."""
    s = E.arclength(P, True)
    L = s[-1]
    w, v = np.linalg.eigh(np.cov((P - P.mean(0)).T))
    k0 = int(np.argmax((P - P.mean(0)) @ v[:, -1]))
    ks = [k0]
    for q in (0.25, 0.5, 0.75):
        target = (s[k0] + q * L) % L
        ks.append(int(np.argmin(np.abs(((s - target + L / 2) % L) - L / 2))))
    return ks, s, L


def square_of(t_rel):
    """A position along the outline from corner 0, as a fraction of each side
    (4 sides, t_rel in [0, 4)), onto the unit square's boundary, CCW."""
    k = np.floor(t_rel).astype(int) % 4
    f = t_rel - np.floor(t_rel)
    uv = np.zeros((len(t_rel), 2))
    uv[k == 0] = np.c_[f[k == 0], np.zeros((k == 0).sum())]
    uv[k == 1] = np.c_[np.ones((k == 1).sum()), f[k == 1]]
    uv[k == 2] = np.c_[1.0 - f[k == 2], np.ones((k == 2).sum())]
    uv[k == 3] = np.c_[np.zeros((k == 3).sum()), 1.0 - f[k == 3]]
    return uv


def outline_rel(s, ks, L):
    """Each outline sample's position as side index + fraction (0..4)."""
    s0 = [s[k] for k in ks]
    rel = np.zeros(len(s))
    for i, si in enumerate(s):
        d = [((si - s0[j]) % L) for j in range(4)]
        side_len = [((s0[(j + 1) % 4] - s0[j]) % L) or L for j in range(4)]
        j = int(np.argmin([dj if dj < side_len[jj] else np.inf for jj, dj in enumerate(d)]))
        rel[i] = j + d[j] / side_len[j]
    return rel


def mean_value(V, F, loop, loop_uv):
    """Floater's mean-value embedding: interior vertex i at sum w_ij x_j / sum
    w_ij, w_ij = (tan(a/2) + tan(b/2)) / |x_i - x_j| with a, b the angles at i
    beside edge ij. The weights are positive, so with a convex (square) boundary
    the map cannot fold, and unlike uniform weights they follow the geometry -
    uniform weights on a scan tessellation collapsed whole strips to zero area."""
    n = len(V)
    rows, cols, vals = [], [], []
    for a_, b_, c_ in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        i, j, k = F[:, a_], F[:, b_], F[:, c_]
        eij, eik = V[j] - V[i], V[k] - V[i]
        lij, lik = np.linalg.norm(eij, axis=1), np.linalg.norm(eik, axis=1)
        cosang = np.einsum('ij,ij->i', eij, eik) / np.maximum(lij * lik, 1e-30)
        th = np.arccos(np.clip(cosang, -1.0, 1.0))
        t = np.tan(0.5 * th)
        rows += [i, i]
        cols += [j, k]
        vals += [t / np.maximum(lij, 1e-12), t / np.maximum(lik, 1e-12)]
    W = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(n, n)).tocsr()
    fixed = np.zeros(n, bool)
    fixed[loop] = True
    free = np.where(~fixed)[0]
    UV = np.zeros((n, 2))
    UV[loop] = loop_uv
    deg = np.asarray(W.sum(1)).ravel()
    L = (diags(deg) - W).tocsr()
    Lff = L[free][:, free].tocsc()
    rhs = -L[free][:, fixed] @ UV[fixed]
    UV[free, 0] = spsolve(Lff, rhs[:, 0])
    UV[free, 1] = spsolve(Lff, rhs[:, 1])
    return UV


def loop_square(V, loop, P, ks):
    """The patch's boundary loop onto the square by its OWN arc length (monotone
    by construction), corners at the loop vertices nearest the outline's four."""
    X = V[loop]
    seg = np.linalg.norm(np.roll(X, -1, 0) - X, axis=1)
    sl = np.r_[0.0, np.cumsum(seg)[:-1]]
    Ltot = seg.sum()
    kd = cKDTree(X)
    ci = [int(kd.query(P[k])[1]) for k in ks]
    rel = np.zeros(len(loop))
    for j in range(4):
        a, b = ci[j], ci[(j + 1) % 4]
        span = (sl[b] - sl[a]) % Ltot or Ltot
        idx = np.arange(a, b) if b > a else np.r_[np.arange(a, len(loop)), np.arange(0, b)]
        rel[idx] = j + ((sl[idx] - sl[a]) % Ltot) / span
    return square_of(rel)


def tutte(V, F, loop, loop_uv):
    """Interior vertices at the average of their neighbours, boundary fixed."""
    n = len(V)
    E2 = np.vstack([F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]])
    E2 = np.unique(np.sort(E2, 1), axis=0)
    A = coo_matrix((np.ones(2 * len(E2)), (np.r_[E2[:, 0], E2[:, 1]], np.r_[E2[:, 1], E2[:, 0]])),
                   shape=(n, n)).tocsr()
    deg = np.asarray(A.sum(1)).ravel()
    fixed = np.zeros(n, bool)
    fixed[loop] = True
    free = np.where(~fixed)[0]
    UV = np.zeros((n, 2))
    UV[loop] = loop_uv
    L = (diags(deg) - A).tocsr()
    Lff = L[free][:, free]
    rhs = -L[free][:, fixed] @ UV[fixed]
    UV[free, 0] = spsolve(Lff.tocsc(), rhs[:, 0])
    UV[free, 1] = spsolve(Lff.tocsc(), rhs[:, 1])
    return UV


# --------------------------------------------------------------------------- #
# fitting
# --------------------------------------------------------------------------- #

def fit_curve(t, X, n, first, last, fair):
    """Cubic B-spline of n clamped control points through samples X at params
    t, ends pinned to first/last, with a second-difference fairness term."""
    k = clamped_knots(n)
    A = basis(t, k).toarray()
    free = slice(1, n - 1)
    rhs = X - np.outer(A[:, 0], first) - np.outer(A[:, -1], last)
    D = second_diff(n).toarray()
    lam = fair * len(t) / n
    Dfree = D[:, free]
    rhsD = -(np.outer(D[:, 0], first) + np.outer(D[:, -1], last))
    M = np.vstack([A[:, free], math.sqrt(lam) * Dfree])
    R = np.vstack([rhs, math.sqrt(lam) * rhsD])
    C = np.linalg.lstsq(M, R, rcond=None)[0]
    return np.vstack([first, C, last]), k


def fit_surface(uv, X, nu, nv, rows, fair):
    """Bicubic net (nu x nv, [i][j] i along u) through X at uv, the boundary
    rows given: rows = {'v0': (nu,3), 'u1': (nv,3), 'v1': (nu,3), 'u0': (nv,3)}."""
    ku, kv = clamped_knots(nu), clamped_knots(nv)
    Bu, Bv = basis(uv[:, 0], ku), basis(uv[:, 1], kv)
    # tensor rows: for each point the 16 products
    # every row of a cubic design matrix holds 4 entries: the tensor row is
    # their 4 x 4 outer product
    m = len(uv)
    iu, au = Bu.indices.reshape(m, 4), Bu.data.reshape(m, 4)
    iv, av = Bv.indices.reshape(m, 4), Bv.data.reshape(m, 4)
    cc = (iu[:, :, None] * nv + iv[:, None, :]).reshape(m, 16)
    vv = (au[:, :, None] * av[:, None, :]).reshape(m, 16)
    rr = np.repeat(np.arange(m), 16)
    A = coo_matrix((vv.ravel(), (rr, cc.ravel())), shape=(m, nu * nv)).tocsr()
    net = np.full((nu, nv, 3), np.nan)
    net[:, 0] = rows['v0']
    net[:, -1] = rows['v1']
    net[-1, :] = rows['u1']
    net[0, :] = rows['u0']
    fixed = ~np.isnan(net[:, :, 0]).ravel()
    free = np.where(~fixed)[0]
    flat = net.reshape(-1, 3)
    # fairness: second differences along u and along v of the whole net
    Du = second_diff(nu)
    Dv = second_diff(nv)
    from scipy.sparse import kron, identity
    D = svstack([kron(Du, identity(nv)), kron(identity(nu), Dv)]).tocsr()
    lam = fair * len(uv) / (nu * nv)
    fixed_idx = np.where(fixed)[0]
    out = flat.copy()
    for d in range(3):
        rhs = X[:, d] - A[:, fixed_idx] @ flat[fixed_idx, d]
        rhsD = -(D[:, fixed_idx] @ flat[fixed_idx, d])
        M = svstack([A[:, free], math.sqrt(lam) * D[:, free]]).tocsr()
        R = np.r_[rhs, math.sqrt(lam) * rhsD]
        out[free, d] = lsqr(M, R, atol=1e-12, btol=1e-12, iter_lim=20000)[0]
    return out.reshape(nu, nv, 3), ku, kv, A


def eval_surface(net, ku, kv, uv):
    nu, nv = net.shape[:2]
    Bu, Bv = basis(uv[:, 0], ku).toarray(), basis(uv[:, 1], kv).toarray()
    return np.einsum('pi,pj,ijd->pd', Bu, Bv, net)


# --------------------------------------------------------------------------- #
# the slab
# --------------------------------------------------------------------------- #

def slab(bone, Q, clearance, thickness, spacing=3.0, fair=0.02, log=None):
    """The base as 6 NURBS faces. bone: bg_engine.Bone (cropped round the outline
    is enough); Q: the outline (mm), on or near the bone. Returns a dict for
    Fusion's BRepBodyDefinition plus how far the faces are from what they model."""
    log = log or (lambda *_: None)
    c, t = float(clearance), float(thickness)
    P, N = E.attach(bone, np.asarray(Q, float), 0.25, 0.0)
    T, Ns, B = E.frames(P, N, 0.25)
    _inside, out = E.inside_labels(bone, P, B)
    if out < 0:
        # make the outline counter-clockwise about the bone normal: then the
        # square's (u, v) gives u x v along the normal
        P, N, Ns = P[::-1].copy(), N[::-1].copy(), Ns[::-1].copy()
    V, F, loop, holes = region(bone, P, N)
    log('patch: %d vertices, %d faces, %d hole(s)' % (len(V), len(F), holes))
    ks, s, L = corners_of(P)
    rel = outline_rel(s, ks, L)
    # the loop must run the outline's way round (its nearest outline samples
    # advance on the whole)
    _d, near = cKDTree(P).query(V[loop])
    d = np.diff(rel[near])
    if np.sum((d + 2.0) % 4.0 - 2.0) < 0:
        loop = loop[::-1]
    UV = mean_value(V, F, loop, loop_square(V, loop, P, ks))
    # normals at the patch vertices (smooth: interpolated from the bone)
    _p, Nv, _f = bone.closest(V)
    side = [(s[ks[(j + 1) % 4]] - s[ks[j]]) % L for j in range(4)]
    Lu, Lv = 0.5 * (side[0] + side[2]), 0.5 * (side[1] + side[3])
    nu = int(np.clip(round(Lu / spacing) + 3, 6, 80))
    nv = int(np.clip(round(Lv / spacing) + 3, 6, 80))
    fair_c = fair

    out_faces = {}
    rep = {'nu': nu, 'nv': nv, 'patch_vertices': int(len(V)), 'holes_in_patch': holes,
           'outline_mm': float(L)}
    edges = {}
    for name, off in (('b', c), ('t', c + t)):
        X = V + off * Nv
        Ep = P + off * Ns                      # the edge: outline, offset along the wall
        cor = [Ep[k] for k in ks]
        rows = {}
        # side j runs from corner j to corner j+1 along the outline; as curves the
        # four go in increasing u or v: v0 c0->c1, u1 c1->c2, v1 c3->c2, u0 c0->c3
        for j, (key, n, rev) in enumerate((('v0', nu, False), ('u1', nv, False),
                                           ('v1', nu, True), ('u0', nv, True))):
            a, b = ks[j], ks[(j + 1) % 4]
            idx = np.arange(a, b + 1) if b > a else np.r_[np.arange(a, len(P)), np.arange(0, b + 1)]
            f = rel[idx] - j
            f[0], f[-1] = 0.0, 1.0
            f = np.clip(f, 0.0, 1.0)
            pts = Ep[idx]
            if rev:
                f, pts = 1.0 - f[::-1], pts[::-1]
            first, last = (cor[(j + 1) % 4], cor[j]) if rev else (cor[j], cor[(j + 1) % 4])
            ctrl, kn = fit_curve(f, pts, n, first, last, fair_c)
            rows[key] = ctrl
            edges['e%s%d' % (name, j)] = {'ctrl': ctrl, 'knots': kn,
                                         'dev_mm': float(np.abs(eval_curve(ctrl, kn, f) - pts).max())}
        sel = subsample(len(X), 40000)
        net, ku, kv, _A = fit_surface(UV[sel], X[sel], nu, nv, rows, fair)
        chk = subsample(len(X), 200000, 2)
        dev = np.linalg.norm(eval_surface(net, ku, kv, UV[chk]) - X[chk], axis=1)
        rep['%s_dev_max_mm' % name] = float(dev.max())
        rep['%s_dev_p99_mm' % name] = float(np.percentile(dev, 99))
        out_faces[name] = (net, ku, kv)
        rep['%s_corners' % name] = [x.tolist() for x in cor]
    rep['edge_dev_max_mm'] = max(e['dev_mm'] for e in edges.values())
    return assemble(out_faces, edges), rep


def eval_curve(ctrl, knots, t):
    return basis(t, knots).toarray() @ ctrl


def assemble(faces, edges):
    """Topology for BRepBodyDefinition: vertices, edges, 6 faces with loops."""
    nb, kub, kvb = faces['b']
    nt, kut, kvt = faces['t']
    nu, nv = nb.shape[:2]
    cb = [nb[0, 0], nb[-1, 0], nb[-1, -1], nb[0, -1]]
    ct = [nt[0, 0], nt[-1, 0], nt[-1, -1], nt[0, -1]]
    V = {'b%d' % i: cb[i].tolist() for i in range(4)}
    V.update({'t%d' % i: ct[i].tolist() for i in range(4)})
    # side j of the square: which corners, as curves (increasing u / v)
    ends = [(0, 1), (1, 2), (3, 2), (0, 3)]
    rows_b = [nb[:, 0], nb[-1, :], nb[:, -1], nb[0, :]]
    rows_t = [nt[:, 0], nt[-1, :], nt[:, -1], nt[0, :]]
    kn = [kub, kvb, kub, kvb]
    Ed = {}
    for j in range(4):
        a, b = ends[j]
        Ed['eb%d' % j] = {'ctrl': rows_b[j].tolist(), 'knots': kn[j].tolist(), 'deg': 3,
                          'start': 'b%d' % a, 'end': 'b%d' % b}
        Ed['et%d' % j] = {'ctrl': rows_t[j].tolist(), 'knots': kn[j].tolist(), 'deg': 3,
                          'start': 't%d' % a, 'end': 't%d' % b}
        Ed['ev%d' % j] = {'ctrl': [cb[j].tolist(), ct[j].tolist()], 'knots': [0, 0, 1, 1],
                          'deg': 1, 'start': 'b%d' % j, 'end': 't%d' % j}

    def surf(net, ku, kv, du=3, dv=3):
        return {'du': du, 'dv': dv, 'nu': int(net.shape[0]), 'nv': int(net.shape[1]),
                'ku': list(map(float, ku)), 'kv': list(map(float, kv)),
                'net': net.reshape(-1, 3).tolist()}           # u-major

    lin = np.array([0.0, 0.0, 1.0, 1.0])
    F = []
    # top (u, v): CCW = v0 fwd, u1 fwd, v1 rev, u0 rev
    F.append({'name': 'top', 'surface': surf(nt, kut, kvt),
              'loop': [('et0', False), ('et1', False), ('et2', True), ('et3', True)]})
    # bottom, transposed (v, u) so its u x v points into the bone
    F.append({'name': 'bottom', 'surface': surf(np.transpose(nb, (1, 0, 2)), kvb, kub),
              'loop': [('eb3', False), ('eb2', False), ('eb1', True), ('eb0', True)]})
    # walls: along the edge, then bottom -> top (sides 0, 1) or top -> bottom
    # (sides 2, 3, whose curves run against the loop) - outward either way
    for j in range(4):
        rb, rt = rows_b[j], rows_t[j]
        if j < 2:
            net = np.stack([rb, rt], 1)
            loop = [('eb%d' % j, False), ('ev%d' % ends[j][1], False),
                    ('et%d' % j, True), ('ev%d' % ends[j][0], True)]
        else:
            net = np.stack([rt, rb], 1)
            loop = [('et%d' % j, False), ('ev%d' % ends[j][1], True),
                    ('eb%d' % j, True), ('ev%d' % ends[j][0], False)]
        F.append({'name': 'wall%d' % j, 'surface': surf(net, kn[j], lin, 3, 1), 'loop': loop})
    return {'vertices': V, 'edges': Ed, 'faces': F}


def crop(V, F, lo, hi):
    """The triangles whose centroid lies in [lo, hi], reindexed - in numpy, so a
    scan of millions of triangles is never handed whole to trimesh."""
    c = V[F].mean(1)
    keep = np.all((c >= lo) & (c <= hi), axis=1)
    Fk = F[keep]
    used, inv = np.unique(Fk, return_inverse=True)
    return V[used], inv.reshape(-1, 3)


def subsample(n, cap, seed=1):
    if n <= cap:
        return np.arange(n)
    return np.sort(np.random.default_rng(seed).choice(n, cap, replace=False))


# --------------------------------------------------------------------------- #
# the slab over a conformal map, cut later by the outline's curtain
# --------------------------------------------------------------------------- #

def lscm(V, F):
    """Least-squares conformal map (Levy et al.) with a free boundary: angles are
    kept, nothing is forced onto a square, so tabs and necks of the outline are
    not squeezed. Two boundary-far vertices pinned; scaled to millimetres."""
    p0, p1, p2 = V[F[:, 0]], V[F[:, 1]], V[F[:, 2]]
    e1 = p1 - p0
    n = np.cross(e1, p2 - p0)
    a2 = np.linalg.norm(n, axis=1)
    ok = a2 > 1e-10
    F, p0, p1, p2, e1, n, a2 = F[ok], p0[ok], p1[ok], p2[ok], e1[ok], n[ok], a2[ok]
    ex = e1 / np.linalg.norm(e1, axis=1)[:, None]
    ey = np.cross(n / a2[:, None], ex)
    X = np.stack([np.zeros(len(F)), np.linalg.norm(e1, axis=1), np.einsum('ij,ij->i', p2 - p0, ex)], 1)
    Y = np.stack([np.zeros(len(F)), np.zeros(len(F)), np.einsum('ij,ij->i', p2 - p0, ey)], 1)
    s = 1.0 / np.sqrt(a2)
    nv, nf = len(V), len(F)
    rows, cols, vals = [], [], []
    for j in range(3):
        wr = (X[:, (j + 2) % 3] - X[:, (j + 1) % 3]) * s
        wi = (Y[:, (j + 2) % 3] - Y[:, (j + 1) % 3]) * s
        f = np.arange(nf)
        # Re: wr u - wi v ; Im: wi u + wr v   (unknowns [u; v])
        rows += [f, f, nf + f, nf + f]
        cols += [F[:, j], nv + F[:, j], F[:, j], nv + F[:, j]]
        vals += [wr, -wi, wi, wr]
    A = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(2 * nf, 2 * nv)).tocsr()
    # pins: the two vertices farthest apart (approximately)
    i0 = int(np.argmax(np.linalg.norm(V - V.mean(0), axis=1)))
    i1 = int(np.argmax(np.linalg.norm(V - V[i0], axis=1)))
    d = float(np.linalg.norm(V[i1] - V[i0]))
    pin_idx = np.array([i0, i1, nv + i0, nv + i1])
    pin_val = np.array([0.0, d, 0.0, 0.0])
    free = np.setdiff1d(np.arange(2 * nv), pin_idx)
    rhs = -A[:, pin_idx] @ pin_val
    Af = A[:, free]
    sol = spsolve((Af.T @ Af).tocsc(), Af.T @ rhs)
    x = np.zeros(2 * nv)
    x[free] = sol
    x[pin_idx] = pin_val
    UV = np.c_[x[:nv], x[nv:]]
    # millimetres: total map area = total surface area
    uva = 0.5 * np.abs(np.cross(UV[F[:, 1]] - UV[F[:, 0]], UV[F[:, 2]] - UV[F[:, 0]]))
    UV *= math.sqrt(0.5 * a2.sum() / max(uva.sum(), 1e-30))
    return UV


def slab_cut(bone, Q, clearance, thickness, spacing=3.0, margin=4.0, fair=0.02, log=None,
             rail_spacing=1.0, rim=0.0):
    """The base as a 6-face NURBS slab over the bone around the outline, and the
    curtain that cuts it to the outline in Fusion. Q: the outline (mm); rim: the
    round wanted on the top rim (mm) - convex bends of the outline tighter than
    1.25 x it are eased, since no round of that radius fits inside them."""
    log = log or (lambda *_: None)
    c, t = float(clearance), float(thickness)
    P, N = E.attach(bone, np.asarray(Q, float), 0.25, 0.0)
    T, Ns, B = E.frames(P, N, 0.25)
    inside_v, out = E.inside_labels(bone, P, B)
    rep0 = {}
    if rim > 0:
        P0 = P.copy()
        P, N, eased = ease_convex(bone, P, N, out, 1.25 * float(rim),
                                  heights=(0.0,))
        if eased:
            T, Ns, B = E.frames(P, N, 0.25)
            inside_v, out = E.inside_labels(bone, P, B)
            rep0 = {'eased_mm': eased * 0.25,
                    'eased_moved_max_mm': float(np.linalg.norm(P - P0, axis=1).max())}
    m = bone.mesh
    Vb = np.asarray(m.vertices)
    # the patch and a margin round it: inside, or within `margin` of the outline
    dq, _i = cKDTree(P).query(Vb, distance_upper_bound=margin, workers=-1)
    keepv = inside_v | np.isfinite(dq)
    keep = np.where(keepv[m.faces].all(1))[0]
    sub = m.submesh([keep], append=True)
    sub.merge_vertices()
    g = coo_matrix((np.ones(len(sub.edges_unique)), (sub.edges_unique[:, 0],
                                                      sub.edges_unique[:, 1])),
                   shape=(len(sub.vertices),) * 2)
    _n, lab = connected_components(g, directed=False)
    big = np.argmax(np.bincount(lab))
    sub = sub.submesh([np.where(lab[sub.faces[:, 0]] == big)[0]], append=True)
    sub.merge_vertices()
    V, F = np.asarray(sub.vertices), np.asarray(sub.faces)
    log('patch + %.0f mm margin: %d vertices' % (margin, len(V)))
    UV = lscm(V, F)
    _p, Nv, _f = bone.closest(V)
    # orient the map so u x v follows the bone normal
    tri = F[np.argmax(np.linalg.norm(np.cross(V[F[:, 1]] - V[F[:, 0]], V[F[:, 2]] - V[F[:, 0]]), axis=1))]
    a3 = np.cross(V[tri[1]] - V[tri[0]], V[tri[2]] - V[tri[0]])
    a2d = np.cross(UV[tri[1]] - UV[tri[0]], UV[tri[2]] - UV[tri[0]])
    if (a3 @ Nv[tri].mean(0)) * a2d < 0:
        UV[:, 1] = -UV[:, 1]
    lo, hi = UV.min(0) - 1.0, UV.max(0) + 1.0
    Lu, Lv = hi - lo
    W = (UV - lo) / (hi - lo)
    nu = int(np.clip(round(Lu / spacing) + 3, 6, 90))
    nv = int(np.clip(round(Lv / spacing) + 3, 6, 90))
    # where the fit matters: inside the outline
    _pp, ins = cKDTree(Vb).query(V)
    inside = inside_v[ins]
    faces, rep = {}, {'nu': nu, 'nv': nv, 'map_mm': [float(Lu), float(Lv)], **rep0,
                      'patch_vertices': int(inside.sum()), 'outline_mm': float(E.arclength(P)[-1])}
    # the contact face: the bone offset by the clearance, fitted
    X = V + c * Nv
    sel = subsample(len(X), 40000)
    net, ku, kv = fit_free(W[sel], X[sel], nu, nv, fair)
    chk = np.where(inside)[0]
    chk = chk[subsample(len(chk), 200000, 2)]
    dev = np.linalg.norm(eval_surface(net, ku, kv, W[chk]) - X[chk], axis=1)
    rep['b_dev_max_mm'] = float(dev.max())
    rep['b_dev_p99_mm'] = float(np.percentile(dev, 99))
    faces['b'] = (net, ku, kv)
    # the top: NOT the bone offset by c + t - in a concave spot tighter than the
    # thickness that offset folds over itself and no surface can follow it - but
    # the contact face's control points moved by t along its own normal
    gu = np.array([ku[i + 1:i + 4].mean() for i in range(nu)])
    gv = np.array([kv[j + 1:j + 4].mean() for j in range(nv)])
    G = np.stack(np.meshgrid(gu, gv, indexing='ij'), -1).reshape(-1, 2)
    h = 1e-4
    du = eval_surface(net, ku, kv, np.c_[np.clip(G[:, 0] + h, 0, 1), G[:, 1]]) -         eval_surface(net, ku, kv, np.c_[np.clip(G[:, 0] - h, 0, 1), G[:, 1]])
    dv = eval_surface(net, ku, kv, np.c_[G[:, 0], np.clip(G[:, 1] + h, 0, 1)]) -         eval_surface(net, ku, kv, np.c_[G[:, 0], np.clip(G[:, 1] - h, 0, 1)])
    nrm = np.cross(du, dv)
    nrm /= np.maximum(np.linalg.norm(nrm, axis=1), 1e-30)[:, None]
    top = net + t * nrm.reshape(nu, nv, 3)
    faces['t'] = (top, ku, kv)
    # the thickness as built, inside the outline: top points to the contact face
    Pb = eval_surface(net, ku, kv, W[chk])
    Pt = eval_surface(top, ku, kv, W[chk])
    dense = np.stack(np.meshgrid(np.linspace(0, 1, 600), np.linspace(0, 1, 400),
                                 indexing='ij'), -1).reshape(-1, 2)
    thick, _ = cKDTree(eval_surface(net, ku, kv, dense)).query(Pt, workers=-1)
    rep['thickness_min_mm'] = float(np.percentile(thick, 0.5))
    rep['thickness_max_mm'] = float(np.percentile(thick, 99.5))
    topo = assemble_free(faces)
    # the curtain: rulings along the smoothed bone normal through the projected
    # outline, from just under the contact face to just over the top AS BUILT.
    # Not a fixed height: where the outline crosses a hollow, lines along the
    # normals meet a few mm up and a taller band folds over itself (Fusion then
    # refuses it) - and the top is thinner there anyway.
    dense_t = cKDTree(eval_surface(top, ku, kv, dense))
    _d, kk = dense_t.query(P + (c + t) * Ns, workers=-1)
    lam_top = np.einsum('ij,ij->i', eval_surface(top, ku, kv, dense[kk]) - P, Ns)
    # the curtain: rulings along the smoothed bone normal through the projected
    # outline. Where the outline dips into a groove tighter than the base is tall,
    # rulings standing on it meet at h_f = -1 / (dN/ds . T) - below the top, so no
    # wall along the normals exists there and no fillet can be built on it. Only
    # there, the cut line bridges the groove: smoothed in 3D (lifted off the bone,
    # the band lowered by as much) until h_f clears the top. Elsewhere it IS the
    # projected outline.
    from scipy.ndimage import gaussian_filter1d, maximum_filter1d
    m_ = len(P)
    need = lam_top + 0.4

    def fold_height(C, Nn):
        sC = E.arclength(C)[:m_]
        ds = np.roll(sC, -1) - np.roll(sC, 1)
        ds[0], ds[-1] = ds[1], ds[-2]
        Tc = _unit3(np.roll(C, -1, 0) - np.roll(C, 1, 0))
        dN = (np.roll(Nn, -1, 0) - np.roll(Nn, 1, 0)) / np.maximum(np.abs(ds), 1e-6)[:, None]
        k = np.einsum('ij,ij->i', dN, Tc)
        return np.where(k < -1e-6, -1.0 / np.minimum(k, -1e-6), np.inf), Tc
    hf, Tc = fold_height(P, Ns)
    bad = hf < need
    C = P.copy()
    rep['cut_bridged_mm'] = 0.0
    if bad.any():
        win = maximum_filter1d(bad.astype(float), int(8.0 / 0.25), mode='wrap') > 0
        w = gaussian_filter1d(win.astype(float), 2.0 / 0.25, mode='wrap')
        w = np.clip(w / max(w.max(), 1e-9), 0.0, 1.0)
        Nw = Ns
        for sig in (1.0, 2.0, 3.0, 5.0, 8.0, 12.0):
            Csm = gaussian_filter1d(P, sig / 0.25, axis=0, mode='wrap')
            C = (1 - w[:, None]) * P + w[:, None] * Csm
            Nsm = gaussian_filter1d(Ns, sig / 0.25, axis=0, mode='wrap')
            Nw = _unit3((1 - w[:, None]) * Ns + w[:, None] * Nsm)
            hf, Tc = fold_height(C, Nw)
            if not (hf < need).any():
                break
        Ns = Nw
        rep['cut_bridged_mm'] = float(win.sum() * 0.25)
        rep['cut_bridge_sigma_mm'] = sig
        rep['cut_moved_max_mm'] = float(np.linalg.norm(C - P, axis=1).max())
    rep['cut_fold_clear'] = bool(not (hf < need).any())
    # where each ruling C + lam Ns actually crosses the top and the contact face
    # (signed distance to each surface along lam, first change of sign): the
    # band has to clear both by a margin everywhere, or Fusion's split refuses
    # it as not intersecting
    lam_t = crossing(top, ku, kv, dense, C, Ns, np.arange(1.0, c + t + 3.0, 0.02))
    lam_b = crossing(net, ku, kv, dense, C, Ns, np.arange(-3.0, c + 1.5, 0.02))
    hi = lam_t + 1.0
    lo = lam_b - 1.0
    # rulings meet above the cut where the bone is concave along it, below it
    # where it is convex (at 1 / |dN/ds . T|): stop short of either, keeping at
    # least 0.3 mm past the faces
    sC = E.arclength(C)[:m_]
    ds_ = np.roll(sC, -1) - np.roll(sC, 1)
    ds_[0], ds_[-1] = ds_[1], ds_[-2]
    kk = np.einsum('ij,ij->i', (np.roll(Ns, -1, 0) - np.roll(Ns, 1, 0)) /
                   np.maximum(np.abs(ds_), 1e-6)[:, None], Tc)
    from scipy.ndimage import minimum_filter1d
    up_ok = minimum_filter1d(np.where(kk < -1e-6, 0.85 / np.maximum(-kk, 1e-6), np.inf), 9, mode='wrap')
    down_ok = minimum_filter1d(np.where(kk > 1e-6, 0.85 / np.maximum(kk, 1e-6), np.inf), 9, mode='wrap')
    hi = np.maximum(np.minimum(hi, up_ok), lam_t + 0.3)
    lo = np.minimum(np.maximum(lo, -down_ok), lam_b - 0.3)
    rep['band_over_top_mm'] = float((hi - lam_t).min())
    rep['band_under_contact_mm'] = float((lam_b - lo).min())
    curtain, crep = curtain_rails(C, Ns, Tc, lo, hi, spacing=rail_spacing)
    # the same, as samples every ~1 mm: for a band made in Fusion from one fitted
    # spline through the cut line, its control points offset to lo / hi
    every = 4
    curtain['samples'] = {'centre': C[::every].tolist(), 'normals': Ns[::every].tolist(),
                          'lo': lo[::every].tolist(), 'hi': hi[::every].tolist()}
    rep.update(crep)

    # a point surely inside the base, to keep the right piece after the cut
    far = np.where(inside_v)[0]
    dd, _ = cKDTree(P).query(Vb[far])
    k = far[int(np.argmax(dd))]
    _pk, nk, _fk = bone.closest(Vb[k][None])
    keep_pt = {'inside': (Vb[k] + (c + 0.5 * t) * nk[0]).tolist(),
               'top': (Vb[k] + (c + t) * nk[0]).tolist(),
               'bottom': (Vb[k] + c * nk[0]).tolist()}
    rep['keep_depth_mm'] = float(dd.max())
    return topo, curtain, keep_pt, rep


def fit_free(W, X, nu, nv, fair):
    """Bicubic net through X at W in [0,1]^2, no boundary fixed, faired."""
    ku, kv = clamped_knots(nu), clamped_knots(nv)
    Bu, Bv = basis(W[:, 0], ku), basis(W[:, 1], kv)
    m = len(W)
    iu, au = Bu.indices.reshape(m, 4), Bu.data.reshape(m, 4)
    iv, av = Bv.indices.reshape(m, 4), Bv.data.reshape(m, 4)
    cc = (iu[:, :, None] * nv + iv[:, None, :]).reshape(m, 16)
    vv = (au[:, :, None] * av[:, None, :]).reshape(m, 16)
    A = coo_matrix((vv.ravel(), (np.repeat(np.arange(m), 16), cc.ravel())),
                   shape=(m, nu * nv)).tocsr()
    from scipy.sparse import kron, identity
    D = svstack([kron(second_diff(nu), identity(nv)), kron(identity(nu), second_diff(nv))]).tocsr()
    lam = fair * m / (nu * nv)
    M = svstack([A, math.sqrt(lam) * D]).tocsr()
    out = np.zeros((nu * nv, 3))
    for d in range(3):
        R = np.r_[X[:, d], np.zeros(D.shape[0])]
        out[:, d] = lsqr(M, R, atol=1e-12, btol=1e-12, iter_lim=50000)[0]
    return out.reshape(nu, nv, 3), ku, kv


def assemble_free(faces):
    """The 6-face slab between two free nets of the same size and knots: the
    edges are the nets' own boundary rows (shared exactly by bottom/top and the
    ruled walls)."""
    nb, ku, kv = faces['b']
    nt, _ku, _kv = faces['t']
    return assemble({'b': (nb, ku, kv), 't': (nt, ku, kv)}, {})


# --------------------------------------------------------------------------- #
# the curtain's two rails: closed cubics sharing their knots
# --------------------------------------------------------------------------- #

def periodic_basis(x, n):
    """Uniform periodic cubic B-spline basis at x in [0, n): (index, weight) x 4."""
    i = np.floor(x).astype(int) % n
    u = x - np.floor(x)
    w = np.stack([(1 - u) ** 3 / 6.0, (3 * u ** 3 - 6 * u ** 2 + 4) / 6.0,
                  (-3 * u ** 3 + 3 * u ** 2 + 3 * u + 1) / 6.0, u ** 3 / 6.0], 1)
    idx = np.stack([(i - 1) % n, i % n, (i + 1) % n, (i + 2) % n], 1)
    return idx, w


def periodic_fit(x, X, n, fair=1e-3):
    """Control points c_0..c_n-1 of the periodic uniform cubic through X at x."""
    idx, w = periodic_basis(x, n)
    m = len(x)
    A = coo_matrix((w.ravel(), (np.repeat(np.arange(m), 4), idx.ravel())), shape=(m, n)).tocsr()
    D = coo_matrix((np.tile([1.0, -2.0, 1.0], n),
                    (np.repeat(np.arange(n), 3),
                     np.stack([np.arange(n) - 1, np.arange(n), np.arange(n) + 1], 1).ravel() % n)),
                   shape=(n, n)).tocsr()
    lam = fair * m / n
    M = svstack([A, math.sqrt(lam) * D]).tocsr()
    C = np.zeros((n, 3))
    for d in range(3):
        C[:, d] = lsqr(M, np.r_[X[:, d], np.zeros(n)], atol=1e-12, btol=1e-12, iter_lim=20000)[0]
    return C


def periodic_eval(C, x):
    idx, w = periodic_basis(x, len(C))
    return np.einsum('pk,pkd->pd', w, C[idx])


def clamped_closed(C):
    """The periodic uniform cubic with controls C as an ordinary CLAMPED cubic
    (first control = last), knots normalised to [0, 1]: the open form
    [c_n-1, c_0..c_n-1, c_0, c_1] on knots -3..n+3, clamped at 0 and n by
    inserting each three times (Boehm), then cut to [0, n]."""
    from scipy.interpolate import insert
    n = len(C)
    Cw = np.vstack([C[-1:], C, C[:2]])
    t = np.arange(-3, n + 4, dtype=float)
    out = []
    for d in range(3):
        tck = (t, np.r_[Cw[:, d], np.zeros(4)], 3)
        for x in (0.0, float(n)):
            for _ in range(3):
                tck = insert(x, tck)
        out.append(tck)
    tk = out[0][0]
    i0 = int(np.searchsorted(tk, 0.0, side='left'))
    i1 = int(np.searchsorted(tk, float(n), side='right'))
    knots = tk[i0:i1]
    ncp = len(knots) - 4
    ctrl = np.stack([out[d][1][i0:i0 + ncp] for d in range(3)], 1)
    return ctrl, knots / float(n)


def curtain_rails(P, Ns, T, lo, hi, step=0.25, spacing=1.0):
    """Bottom rail P + lo Ns exactly; top rail P + hi Ns, smoothed along the
    outline ONLY where it would fold back on itself (an outline crossing a hollow
    tighter than the band is tall), until it does not. Both as clamped closed
    cubics on the same knots - a ruling joins equal parameters."""
    m = len(P)
    s = E.arclength(P)[:m]
    L = E.arclength(P)[-1]
    n = max(16, int(round(L / spacing)))
    x = s / L * n
    lo_pts = P + lo[:, None] * Ns
    hi_pts = P + hi[:, None] * Ns

    def folds(Y):
        d = np.roll(Y, -1, 0) - Y
        return np.einsum('ij,ij->i', d, T) < 0
    bad = folds(hi_pts)
    smoothed = 0
    sig = 2.0
    from scipy.ndimage import gaussian_filter1d
    while bad.any() and sig <= 40.0:
        # widen the trouble spots, blend toward a smoothed rail only there
        win = np.convolve(np.r_[bad[-40:], bad, bad[:40]].astype(float),
                          np.ones(81), 'same')[40:-40] > 0
        wgt = gaussian_filter1d(win.astype(float), 8.0, mode='wrap')
        wgt = np.clip(wgt / max(wgt.max(), 1e-9), 0, 1)
        sm = gaussian_filter1d(hi_pts, sig / step, axis=0, mode='wrap')
        hi_pts = (1 - wgt[:, None]) * hi_pts + wgt[:, None] * sm
        smoothed += int(win.sum())
        bad = folds(hi_pts)
        sig *= 1.5
    Clo = periodic_fit(x, lo_pts, n)
    Chi = periodic_fit(x, hi_pts, n)
    ylo, yhi = periodic_eval(Clo, x), periodic_eval(Chi, x)
    rep = {'rail_lo_dev_mm': float(np.linalg.norm(ylo - lo_pts, axis=1).max()),
           'rail_hi_dev_mm': float(np.linalg.norm(yhi - hi_pts, axis=1).max()),
           'rail_hi_smoothed_samples': smoothed, 'rail_hi_folds_left': int(folds(yhi).sum()),
           'rail_lo_folds': int(folds(ylo).sum())}
    a, ka = clamped_closed(Clo)
    b, kb = clamped_closed(Chi)
    return {'lo': {'ctrl': a.tolist(), 'knots': ka.tolist()},
            'hi': {'ctrl': b.tolist(), 'knots': kb.tolist()}}, rep


def _unit3(V):
    return V / np.maximum(np.linalg.norm(V, axis=1), 1e-30)[:, None]


def crossing(net, ku, kv, dense, C, D, lams):
    """For each ruling C_i + lam D_i, the first lam (from lams, ascending) where
    it passes from one side of the surface to the other; linear in between."""
    S = eval_surface(net, ku, kv, dense)
    h = 1e-4
    du = eval_surface(net, ku, kv, np.c_[np.clip(dense[:, 0] + h, 0, 1), dense[:, 1]]) - S
    dv = eval_surface(net, ku, kv, np.c_[dense[:, 0], np.clip(dense[:, 1] + h, 0, 1)]) - S
    Nn = np.cross(du, dv)
    Nn /= np.maximum(np.linalg.norm(Nn, axis=1), 1e-30)[:, None]
    kd = cKDTree(S)
    m = len(C)
    X = C[None, :, :] + lams[:, None, None] * D[None, :, :]
    _d, k = kd.query(X.reshape(-1, 3), workers=-1)
    sd = np.einsum('ij,ij->i', X.reshape(-1, 3) - S[k], Nn[k]).reshape(len(lams), m)
    sgn = np.sign(sd)
    out = np.full(m, np.nan)
    for i in range(m):
        ch = np.nonzero(sgn[:-1, i] != sgn[1:, i])[0]
        if len(ch):
            a = ch[0]
            f = sd[a, i] / (sd[a, i] - sd[a + 1, i])
            out[i] = lams[a] + f * (lams[a + 1] - lams[a])
    bad = np.isnan(out)
    if bad.any():
        out[bad] = np.interp(np.nonzero(bad)[0], np.nonzero(~bad)[0], out[~bad])
    return out


def ease_convex(bone, P, N, out, r_min, step=0.25, heights=(0.0,)):
    """Convex bends of the outline (turning towards the inside of the base, in
    the bone's tangent plane) tighter than r_min, opened up: the outline is
    smoothed only around them and laid back on the bone, a little more each
    round, until none is left. Returns (P, N, samples eased)."""
    from scipy.ndimage import gaussian_filter1d, maximum_filter1d
    eased = 0
    for _ in range(40):
        T, Ns, B = E.frames(P, N, step)
        bad = np.zeros(len(P), bool)
        # the bend as it is at each height of the wall: over a concave bone the
        # rim bends tighter than the outline on the bone
        for h in heights:
            Y = P + h * Ns
            sY = np.linalg.norm(np.roll(Y, -1, 0) - np.roll(Y, 1, 0), axis=1)
            TY = (np.roll(Y, -1, 0) - np.roll(Y, 1, 0)) / np.maximum(sY, 1e-9)[:, None]
            dT = (np.roll(TY, -1, 0) - np.roll(TY, 1, 0)) / np.maximum(sY, 1e-9)[:, None]
            kg = gaussian_filter1d(-np.einsum('ij,ij->i', dT, B * out), 2, mode='wrap')
            bad |= kg > 1.0 / r_min
        if not bad.any():
            break
        win = maximum_filter1d(bad.astype(float), int(4.0 / step), mode='wrap') > 0
        eased = max(eased, int(win.sum()))
        w = gaussian_filter1d(win.astype(float), 1.0 / step, mode='wrap')
        w = np.clip(w / max(w.max(), 1e-9), 0.0, 1.0)
        Psm = gaussian_filter1d(P, 1.0 / step, axis=0, mode='wrap')
        P = (1 - w[:, None]) * P + w[:, None] * Psm
        P, N, _f = bone.closest(P)
    return P, N, eased
