"""Bone Guide - all the working code. RE-IMPORTED ON EVERY COMMAND EVENT.

Because of that reload, nothing may keep state in a module global between events;
per-invocation state lives in the _bg_state module (state()).

TWO COMMANDS, 3-matic's guide workflow in Fusion
    Bone Curve   a closed curve ON the bone: made from drill cylinders and other
                 bodies (their alpha shape - margin + closing - cut with the
                 bone), from a sketch curve, or by clicks. It is a spline through
                 fit points on the bone, edited like a sketch spline: drag a point
                 over the bone (neighbours within the influence distance follow -
                 3-matic's Edit Curve), click the curve to add a point,
                 Shift+click to remove one. OK writes it back into its spline
                 (fit points moved in place) or into a new 3D sketch.
    Guide Base   the base inside a closed curve (3-matic's Create Base): contact
                 face = bone offset by the clearance, top = bone offset by
                 clearance + thickness, side wall along the bone normals at the
                 outline, rim rounded, holes cut along the selected cylinders.
                 Built as one field by bg_engine.py and imported as a mesh body.

THE MATHS LIVES IN bg_engine.py (py -3.10, numpy/scipy/trimesh/pyvista):
    as a subprocess for Guide Base, and as a resident server (one JSON line per
    request) for the interactive editor, so a drag costs milliseconds.

UNITS
    Fusion internal is centimetres; the engine, the dialog and the reports are mm.
"""

import adsk.core
import adsk.fusion

import array
import hashlib
import importlib
import json
import math
import os
import subprocess
import sys
import time
import traceback
import types

ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))
ENGINE = os.path.join(ADDIN_DIR, 'bg_engine.py')
WORK = os.path.join(ADDIN_DIR, 'work')
LOG_PATH = os.path.join(ADDIN_DIR, 'BoneGuide.log')
SETTINGS_PATH = os.path.join(ADDIN_DIR, 'last_settings.json')
RES_DIR = os.path.join(ADDIN_DIR, 'resources')
VERSION = '0.3.0'
ATTR = 'BoneGuide'

CM = 10.0          # internal cm -> mm
MM = 0.1           # mm -> internal cm
NO_WINDOW = 0x08000000 if os.name == 'nt' else 0

CMD_CURVE = 'BoneGuide_Curve'
CMD_BASE = 'BoneGuide_Base'
CMD_NAMES = {CMD_CURVE: 'Bone Curve', CMD_BASE: 'Guide Base'}
PANELS = ['SolidCreatePanel', 'SurfaceCreatePanel']

DEFAULTS = {'influence': 0.0, 'smooth_mm': 1.0, 'clearance': 0.2, 'thickness': 3.5,
            'round': 2.0, 'round_contact': 0.2, 'voxel': 1, 'curve_smooth': 0.5,
            'fit_spacing': 3.0, 'margin': 2.4, 'margin2': 1.0, 'closing': 3.5,
            'closing_mode': 0, 'join': True,
            'guide_d': 7.0, 'guide_h': 12.1, 'guide_clear': 0.2, 'guide_top_round': 1.5,
            'guide_edge_round': 0.0, 'guide_spin': 0.0, 'hole_d': 3.5, 'output': 1,
            'pent_d': 6.0, 'pent_hole_d': 2.0, 'pent_h': 10.5, 'pent_top_round': 1.0,
            'facet_tol': 0.05}
# the drill guides of a distal radius test case, measured: hole = the reference
# cylinder, outside 5.2 mm (a pentagon inscribed in it), top 10 mm above the bone along the
# axis, bottom on the bone offset by 0.5 mm, 1 mm rounds (cylinder: top rim;
# pentagon: its five edges, top left sharp)
OUTPUTS = [('mesh', 'Mesh body'), ('solid', 'Solid body (faceted)')]
VOXELS = [('0.10', 'Fine - 0.10 mm voxels'), ('0.15', 'Medium - 0.15 mm voxels'),
          ('0.25', 'Draft - 0.25 mm voxels')]
CLOSINGS = [('planes', 'In planes across the drills (alpha shape)'),
            ('ball', '3D')]
STEP = 0.25                # curve sampling of the editor, mm
CURVE_RGB = (0, 110, 220)  # a sketch curve's blue
CURVE_WEIGHT = 3.0
POINT_IMG = os.path.join(RES_DIR, 'points', 'fit.png')
HOT_IMG = os.path.join(RES_DIR, 'points', 'hot.png')
PICK_PX = 9.0              # a fit point or the curve is under the cursor within this
LIFT = 0.15                # the overlay is drawn this far off the bone, mm


def log(msg):
    try:
        with open(LOG_PATH, 'a', encoding='utf-8') as fh:
            fh.write('%s  %s\n' % (time.strftime('%H:%M:%S'), msg))
    except Exception:
        pass


def short_error(tb):
    for line in reversed(tb.strip().splitlines()):
        line = line.strip()
        if line and not line.startswith(('File "', 'Traceback', '^', '~')):
            return line
    return tb.strip().splitlines()[-1]


def _esc(t):
    return str(t).replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')


def state():
    st = sys.modules.get('_bg_state')
    if st is None:
        st = types.ModuleType('_bg_state')
        st.session = None
        st.command = None
        st.handlers = []
        st.server = None
        st.bones = {}
        st.python = None
        sys.modules['_bg_state'] = st
    return st


def load_settings():
    o = dict(DEFAULTS)
    try:
        with open(SETTINGS_PATH, encoding='utf-8') as fh:
            o.update({k: v for k, v in json.load(fh).items() if k in DEFAULTS})
    except Exception:
        pass
    return o


def save_settings(o):
    cur = load_settings()
    cur.update({k: v for k, v in o.items() if k in DEFAULTS})
    try:
        with open(SETTINGS_PATH, 'w', encoding='utf-8') as fh:
            json.dump(cur, fh, indent=1)
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# small vector helpers (mm tuples; no numpy inside Fusion)
# --------------------------------------------------------------------------- #

def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _add(a, b, k=1.0):
    return (a[0] + k * b[0], a[1] + k * b[1], a[2] + k * b[2])


def _dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def _mm(p):
    return (p.x * CM, p.y * CM, p.z * CM)


def _pt(q):
    return adsk.core.Point3D.create(q[0] * MM, q[1] * MM, q[2] * MM)


def arc_positions(Q, closed=True):
    s = [0.0]
    for i in range(1, len(Q)):
        s.append(s[-1] + _dist(Q[i - 1], Q[i]))
    total = s[-1] + (_dist(Q[-1], Q[0]) if closed else 0.0)
    return s, total


def resample(Q, step, closed=True):
    P = list(Q) + ([Q[0]] if closed else [])
    s, _ = arc_positions(P, False)
    L = s[-1]
    n = max(8, int(round(L / step)))
    out, j = [], 0
    for k in range(n if closed else n + 1):
        t = L * k / n
        while j < len(s) - 2 and s[j + 1] < t:
            j += 1
        seg = s[j + 1] - s[j]
        u = 0.0 if seg <= 0 else (t - s[j]) / seg
        out.append(_add(P[j], _sub(P[j + 1], P[j]), u))
    return out


# --------------------------------------------------------------------------- #
# the engine: a subprocess for builds, a resident server for the editor
# --------------------------------------------------------------------------- #

def find_python():
    cands = []
    if os.environ.get('BONEGUIDE_PYTHON'):
        cands.append([os.environ['BONEGUIDE_PYTHON']])
    txt = os.path.join(ADDIN_DIR, 'python_path.txt')
    if os.path.exists(txt):
        try:
            with open(txt, encoding='utf-8') as fh:
                line = fh.read().strip()
            if line:
                cands.append([line])
        except Exception:
            pass
    local = os.environ.get('LOCALAPPDATA', '')
    for v in ('310', '311', '312', '313', '39'):
        cands.append([os.path.join(local, 'Programs', 'Python', 'Python%s' % v,
                                   'python.exe')])
    cands.append(['py', '-3.10'])
    for c in cands:
        try:
            r = subprocess.run(c + ['-c', 'import numpy, scipy, trimesh, pyvista'],
                               capture_output=True, timeout=90, creationflags=NO_WINDOW)
            if r.returncode == 0:
                return c
        except Exception:
            continue
    return None


def python_cmd():
    st = state()
    if not st.python:
        st.python = find_python() or []
    if not st.python:
        raise RuntimeError('No Python with numpy, scipy, trimesh and pyvista was found. '
                           'Set BONEGUIDE_PYTHON, or put the interpreter path in '
                           'python_path.txt next to the add-in.')
    return list(st.python)


def _engine_mtime():
    """The newest of the engine's files: the process keeps what it imported."""
    try:
        return max(os.path.getmtime(os.path.join(ADDIN_DIR, f))
                   for f in ('bg_engine.py', 'bg_slab.py') if
                   os.path.exists(os.path.join(ADDIN_DIR, f)))
    except (OSError, ValueError):
        return None


class Server:
    """bg_engine.py --serve: one JSON request per line, one answer per line."""

    def __init__(self):
        self.proc = subprocess.Popen(python_cmd() + [ENGINE, '--serve'],
                                     stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                     stderr=subprocess.DEVNULL, creationflags=NO_WINDOW,
                                     text=True, bufsize=1)
        self.bone_key = None
        self.engine_mtime = _engine_mtime()

    def alive(self):
        return self.proc.poll() is None

    def ask(self, req):
        self.proc.stdin.write(json.dumps(req) + '\n')
        self.proc.stdin.flush()
        line = self.proc.stdout.readline()
        if not line:
            raise RuntimeError('the Bone Guide engine stopped')
        out = json.loads(line)
        if 'error' in out:
            raise (ValueError if out.get('user_error') else RuntimeError)(out['error'])
        return out

    def close(self):
        try:
            self.proc.stdin.write('{"op": "quit"}\n')
            self.proc.stdin.flush()
            self.proc.wait(timeout=3)
        except Exception:
            try:
                self.proc.kill()
            except Exception:
                pass


def engine():
    """The resident engine, with no bone loaded on its account."""
    st = state()
    srv = st.server
    stale = srv is not None and getattr(srv, 'engine_mtime', None) != _engine_mtime()
    if srv is None or stale or not srv.alive():
        if srv is not None:
            srv.close()
        srv = st.server = Server()
    return srv


def server_for(body):
    """The resident engine with this bone loaded (loaded once per bone)."""
    st = state()
    srv = st.server
    # an edited bg_engine.py: the running process holds the old code. Checked
    # here, not in a Server method - a server made by an older bg_core keeps its
    # old class
    stale = srv is not None and getattr(srv, 'engine_mtime', None) != _engine_mtime()
    if srv is None or stale or not srv.alive():
        if srv is not None:
            srv.close()
        srv = st.server = Server()
    key = body_key(body)
    if srv.bone_key != key:
        paths = write_bone(body)
        srv.ask({'op': 'load', 'coords': paths['coords'], 'tris': paths['tris'],
                 'step': 0.25})
        srv.bone_key = key
    return srv


def shutdown():
    st = state()
    if st.server is not None:
        st.server.close()
        st.server = None


def run_engine(job, ui=None, title='Guide Base', timeout=3600):
    os.makedirs(WORK, exist_ok=True)
    jp = os.path.join(WORK, 'job.json')
    op = os.path.join(WORK, 'out.json')
    for f in (op, op + '.log'):
        if os.path.exists(f):
            os.remove(f)
    with open(jp, 'w', encoding='utf-8') as fh:
        json.dump(job, fh)
    proc = subprocess.Popen(python_cmd() + [ENGINE, jp, op], stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, creationflags=NO_WINDOW)
    pd = None
    if ui is not None:
        pd = ui.createProgressDialog()
        pd.isCancelButtonShown = True
        pd.show(title, 'Starting...', 0, 100, 0)
    t0 = time.time()
    try:
        while proc.poll() is None:
            if pd is not None:
                adsk.doEvents()
                if pd.wasCancelled:
                    proc.kill()
                    raise RuntimeError('Cancelled.')
                msg = _last_line(op + '.log')
                pd.message = '%s  (%.0f s)' % (msg or 'Working', time.time() - t0)
            if time.time() - t0 > timeout:
                proc.kill()
                raise RuntimeError('The engine ran past %d s.' % timeout)
            time.sleep(0.15)
    finally:
        if pd is not None:
            pd.hide()
    out, err = proc.communicate()
    if not os.path.exists(op):
        raise RuntimeError('The engine produced nothing\n%s'
                           % (out[-600:] + err[-600:]).decode('utf-8', 'replace'))
    with open(op, encoding='utf-8') as fh:
        res = json.load(fh)
    if not res.get('ok'):
        if res.get('user_error'):
            raise ValueError(res['error'])
        raise RuntimeError('The engine failed:\n%s' % res.get('error', '')[-1500:])
    return res['result']


def _last_line(path):
    try:
        with open(path, encoding='utf-8') as fh:
            lines = fh.read().strip().splitlines()
        return lines[-1].strip() if lines else ''
    except Exception:
        return ''


# --------------------------------------------------------------------------- #
# geometry out of Fusion
# --------------------------------------------------------------------------- #

def body_key(b):
    try:
        return b.entityToken
    except Exception:
        return '%s|%s' % (b.name, id(b))


def _calc_mesh(body):
    """Triangles of a BRep body, tolerance tied to its size (a tight tolerance
    with no side or aspect limit has taken Fusion down)."""
    bb = body.boundingBox
    diag = bb.minPoint.distanceTo(bb.maxPoint) or 1.0
    calc = body.meshManager.createMeshCalculator()
    calc.surfaceTolerance = max(diag / 6000.0, 1.0e-4)
    calc.maxSideLength = diag / 40.0
    calc.maxAspectRatio = 10.0
    calc.maxNormalDeviation = math.radians(8.0)
    tm = calc.calculate()
    return tm.nodeCoordinatesAsDouble, tm.nodeIndices


def body_triangles(body):
    mb = adsk.fusion.MeshBody.cast(body)
    if mb:
        pm = mb.mesh
        if pm and pm.triangleCount > 0 and pm.quadCount == 0 and pm.polygonCount == 0:
            return pm.nodeCoordinatesAsDouble, pm.triangleNodeIndices
        dm = mb.displayMesh
        return dm.nodeCoordinatesAsDouble, dm.nodeIndices
    return _calc_mesh(adsk.fusion.BRepBody.cast(body))


def write_mesh(body):
    """A body's triangles as raw little-endian arrays in mm, cached per body
    (bone or shape) in files named after its token."""
    st = state()
    key = body_key(body)
    hit = st.bones.get(key)
    if hit and all(os.path.exists(p) for p in hit.values()):
        return hit
    c, t = body_triangles(body)
    os.makedirs(WORK, exist_ok=True)
    tag = 'm%s' % hashlib.md5(key.encode('utf-8')).hexdigest()[:12]
    paths = {'coords': os.path.join(WORK, tag + '_v.bin'),
             'tris': os.path.join(WORK, tag + '_f.bin')}
    with open(paths['coords'], 'wb') as fh:
        array.array('d', [v * CM for v in c]).tofile(fh)
    with open(paths['tris'], 'wb') as fh:
        array.array('i', t).tofile(fh)
    st.bones[key] = paths
    return paths


write_bone = write_mesh


def curve_polyline(curve, tol_mm=0.005):
    """A sketch curve as mm points along its true (world) geometry."""
    ev = curve.worldGeometry.evaluator
    ok, t0, t1 = ev.getParameterExtents()
    ok, pts = ev.getStrokes(t0, t1, tol_mm * MM)
    return [_mm(p) for p in pts]


def chain(polys, tol=0.05):
    """Several sketch curves joined end to end into one polyline; closed if its
    ends meet."""
    polys = [list(p) for p in polys if len(p) > 1]
    if not polys:
        return [], False
    out = polys.pop(0)
    while polys:
        best = None
        for i, p in enumerate(polys):
            for rev in (False, True):
                d = _dist(out[-1], p[-1] if rev else p[0])
                if best is None or d < best[0]:
                    best = (d, i, rev)
        d, i, rev = best
        p = polys.pop(i)
        out.extend((p[::-1] if rev else p)[1:])
    closed = _dist(out[0], out[-1]) < max(tol, 1e-3)
    if closed:
        out = out[:-1]
    return out, closed


def cylinders(entities):
    """(origin mm, unit axis, radius mm) of every cylindrical face picked, or of
    every cylindrical face of a picked body; one per axis and radius."""
    faces = []
    for e in entities:
        f = adsk.fusion.BRepFace.cast(e)
        if f:
            faces.append(f)
            continue
        b = adsk.fusion.BRepBody.cast(e)
        if b:
            faces.extend(b.faces.item(i) for i in range(b.faces.count))
    out = []
    for f in faces:
        g = adsk.core.Cylinder.cast(f.geometry)
        if not g:
            continue
        o, a, r = _mm(g.origin), g.axis, g.radius * CM
        n = math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
        a = (a.x / n, a.y / n, a.z / n)
        dup = False
        for o2, a2, r2 in out:
            q = _sub(o, o2)
            cr = (q[1] * a2[2] - q[2] * a2[1], q[2] * a2[0] - q[0] * a2[2],
                  q[0] * a2[1] - q[1] * a2[0])
            if abs(abs(a[0] * a2[0] + a[1] * a2[1] + a[2] * a2[2]) - 1) < 1e-6 and \
                    math.sqrt(sum(c * c for c in cr)) < 1e-3 and abs(r - r2) < 1e-3:
                dup = True
                break
        if not dup:
            out.append((o, a, r))
    return out


def as_cylinder(body):
    """(p0, p1, r) mm when a BRep body is one cylinder (coaxial cylindrical
    faces of one radius, the rest planar), else None."""
    axis = None
    for i in range(body.faces.count):
        f = body.faces.item(i)
        g = adsk.core.Cylinder.cast(f.geometry)
        if g:
            o, a, r = _mm(g.origin), g.axis, g.radius * CM
            n = math.sqrt(a.x * a.x + a.y * a.y + a.z * a.z)
            a = (a.x / n, a.y / n, a.z / n)
            if axis is None:
                axis = (o, a, r)
                continue
            q = _sub(o, axis[0])
            along = sum(q[k] * axis[1][k] for k in range(3))
            off = _dist(q, tuple(along * axis[1][k] for k in range(3)))
            if abs(abs(sum(a[k] * axis[1][k] for k in range(3))) - 1) > 1e-6 or \
                    off > 1e-3 or abs(r - axis[2]) > 1e-3:
                return None
        elif not adsk.core.Plane.cast(f.geometry):
            return None
    if axis is None:
        return None
    o, a, r = axis
    c, _t = _calc_mesh(body)
    ts = [sum((c[i + k] * CM - o[k]) * a[k] for k in range(3)) for i in range(0, len(c), 3)]
    return _add(o, a, min(ts)), _add(o, a, max(ts)), r


def shape_specs(entities, margin):
    """The picked shapes for the engine: a cylinder body as its exact cylinder,
    anything else as its triangles; each with its margin and its name."""
    out = []
    for e in entities:
        b = adsk.fusion.BRepBody.cast(e)
        cyl = as_cylinder(b) if b else None
        if cyl:
            spec = {'kind': 'cyl', 'p0': cyl[0], 'p1': cyl[1], 'r': cyl[2]}
        else:
            paths = write_mesh(e)
            spec = {'kind': 'mesh', 'coords': paths['coords'], 'tris': paths['tris']}
        spec.update({'margin': float(margin), 'name': getattr(e, 'name', '')})
        out.append(spec)
    return out


# --------------------------------------------------------------------------- #
# picking on the bone
# --------------------------------------------------------------------------- #

def ray_from_view(vp, pt2d):
    """Ray through a viewport point (origin, direction), model cm. Both the
    perspective and the orthographic candidates are built and the one that
    reprojects onto the cursor is kept; the orthographic origin is pulled back in
    front of the eye, and the perspective direction forced into the screen (see
    the Fusion painting notes - both traps are silent)."""
    p = vp.viewToModelSpace(pt2d)
    cam = vp.camera
    eye = cam.eye
    vd = adsk.core.Vector3D.create(cam.target.x - eye.x, cam.target.y - eye.y,
                                   cam.target.z - eye.z)
    vd.normalize()
    fwd = (vd.x, vd.y, vd.z)
    focal = max(eye.distanceTo(cam.target), 1.0)
    along = (p.x - eye.x) * fwd[0] + (p.y - eye.y) * fwd[1] + (p.z - eye.z) * fwd[2]
    back = along + focal
    ortho = adsk.core.Point3D.create(p.x - vd.x * back, p.y - vd.y * back,
                                     p.z - vd.z * back)
    cands = []
    to_p = adsk.core.Vector3D.create(p.x - eye.x, p.y - eye.y, p.z - eye.z)
    if to_p.length > 1e-9:
        if to_p.x * fwd[0] + to_p.y * fwd[1] + to_p.z * fwd[2] < 0:
            to_p = adsk.core.Vector3D.create(-to_p.x, -to_p.y, -to_p.z)
        to_p.normalize()
        cands.append((eye, to_p))
    cands.append((ortho, vd))
    scale = max(1.0, eye.distanceTo(cam.target)) * 0.1
    best, err_best = None, None
    for o, d in cands:
        probe = adsk.core.Point3D.create(o.x + d.x * scale, o.y + d.y * scale,
                                         o.z + d.z * scale)
        b = vp.modelToViewSpace(probe)
        if b is None:
            continue
        err = math.hypot(b.x - pt2d.x, b.y - pt2d.y)
        if err_best is None or err < err_best:
            best, err_best = (o, d), err
    return best if best else (eye, vd)


def hit_bone(app, body, origin, direction):
    """Nearest point where the ray meets the bone, mm, or None. Asked of the
    engine's cell locator (0.1 ms): Fusion's own ray test against a scanned
    bone of 1.85 M triangles ran on every drag event and froze the editor."""
    try:
        r = server_for(body).ask({'op': 'ray', 'o': _mm(origin),
                                  'd': (direction.x, direction.y, direction.z)})
    except Exception:
        log(traceback.format_exc())
        return None
    return tuple(r['hit']) if r['hit'] else None


def screen_xy(vp, q):
    p = vp.modelToViewSpace(_pt(q))
    return (p.x, p.y) if p is not None else None


def _view_key(vp):
    cam = vp.camera
    e, t, u = cam.eye, cam.target, cam.upVector
    return (round(e.x, 6), round(e.y, 6), round(e.z, 6), round(t.x, 6), round(t.y, 6),
            round(t.z, 6), round(u.x, 6), round(u.y, 6), round(u.z, 6),
            round(cam.viewExtents, 6), vp.width, vp.height)


def fit_screen(vp, s):
    """Screen positions of the fit points, kept until the view or the points
    change: hovering then costs one camera read, not a projection per point."""
    key = (_view_key(vp), tuple(s['fit']))
    c = s.get('screen')
    if c is None or c[0] != key:
        s['screen'] = c = (key, [screen_xy(vp, f) for f in s['fit']])
    return c[1]


def px_per_mm(vp, q):
    """Screen pixels per mm of model at q, across the view."""
    up = vp.camera.upVector
    n = math.sqrt(up.x * up.x + up.y * up.y + up.z * up.z) or 1.0
    a = screen_xy(vp, q)
    b = screen_xy(vp, _add(q, (up.x / n, up.y / n, up.z / n)))
    if a is None or b is None:
        return 1.0
    return max(math.hypot(b[0] - a[0], b[1] - a[1]), 1e-6)


def fit_under_cursor(vp, s, cursor):
    """Index of the fit point drawn under the cursor (within PICK_PX), or None."""
    best, bi = None, None
    for k, xy in enumerate(fit_screen(vp, s)):
        if xy is None:
            continue
        d = math.hypot(xy[0] - cursor.x, xy[1] - cursor.y)
        if d <= PICK_PX and (best is None or d < best):
            best, bi = d, k
    return bi


# --------------------------------------------------------------------------- #
# the curve editor (Bone Curve)
#
# The curve is a closed spline through FIT POINTS that lie on the bone, as a
# Fusion fitted spline is; the engine samples it and lays every sample on the
# bone. Dragging a fit point slides it over the bone (the ray hit under the
# cursor) and the spline follows; with an influence distance its neighbours
# along the curve follow too, with a cosine fall-off (3-matic's Edit Curve).
# --------------------------------------------------------------------------- #

def new_session(body, curves=None, shapes=None, p=None, shapes2=None):
    """The bone in the engine, and the curve to edit: the shapes' alpha shape
    cut with the bone, or the picked sketch curve laid on the bone, or nothing
    (clicks then draw one). shapes2 take the second margin; a body picked in
    both boxes counts once, with the second margin."""
    p = p or load_settings()
    shapes2 = list(shapes2 or [])
    own = set(body_key(e) for e in shapes2)
    shapes = [e for e in (shapes or []) if body_key(e) not in own]
    s = {'bone': body, 'fit': [], 'Q': [], 'N': [], 'at': [], 'spline': None,
         'curves': list(curves or []), 'shapes': shapes, 'shapes2': shapes2,
         'drag': None, 'hot': None, 'undo': [], 'stats': {}}
    srv = server_for(body)
    if shapes or shapes2:
        r = srv.ask({'op': 'section',
                     'shapes': shape_specs(shapes, p['margin'])
                     + shape_specs(shapes2, p.get('margin2', p['margin'])),
                     'margin': p['margin'], 'closing': p['closing'],
                     'planes': int(p['closing_mode']) == 0, 'spacing': p['fit_spacing'],
                     'join': bool(p.get('join', True)),
                     'smooth': p['curve_smooth'], 'step': STEP})
        apply_answer(s, r)
        s['stats']['section'] = r['info']
    elif curves:
        spl = adsk.fusion.SketchFittedSpline.cast(curves[0]) if len(curves) == 1 else None
        if spl:
            s['spline'] = spl
            fit = [_mm(spl.fitPoints.item(i).worldGeometry)
                   for i in range(spl.fitPoints.count)]
            # a closed spline may list its first point again at the end
            s['src_idx'] = [i for i in range(len(fit)) if not (
                i == len(fit) - 1 and i > 0 and _dist(fit[i], fit[0]) < 1e-4)]
            fit = [fit[i] for i in s['src_idx']]
            apply_answer(s, srv.ask({'op': 'spline', 'fit': fit, 'step': STEP}))
            s['stats']['attached_max_mm'] = max(_dist(a, b) for a, b in zip(fit, s['fit']))
        else:
            P, closed = chain([curve_polyline(c) for c in curves])
            if not closed:
                raise ValueError('The curve must be closed (its ends %.2f mm apart).'
                                 % _dist(P[0], P[-1]))
            apply_answer(s, srv.ask({'op': 'fit', 'pts': P, 'spacing': p['fit_spacing'],
                                     'step': STEP}))
    return s


def nearest_index(Q, p):
    best, bi = None, -1
    for i, q in enumerate(Q):
        d = (q[0] - p[0]) ** 2 + (q[1] - p[1]) ** 2 + (q[2] - p[2]) ** 2
        if best is None or d < best:
            best, bi = d, i
    return bi, math.sqrt(best) if best is not None else None


def apply_answer(s, r):
    s['fit'] = [tuple(q) for q in r['fit']]
    s['Q'] = [tuple(q) for q in r['pts']]
    s['N'] = [tuple(n) for n in r['nrm']]
    s['at'] = list(r['fit_at'])


def respline(s, smooth_mm=0.0):
    if len(s['fit']) < 3:
        r = server_for(s['bone']).ask({'op': 'attach', 'pts': s['fit']}) if s['fit'] \
            else {'pts': [], 'nrm': []}
        s['fit'] = [tuple(q) for q in r['pts']]
        s['Q'], s['N'], s['at'] = list(s['fit']), [tuple(n) for n in r['nrm']], \
            list(range(len(s['fit'])))
        return
    apply_answer(s, server_for(s['bone']).ask({'op': 'spline', 'fit': s['fit'],
                                               'step': STEP, 'smooth': smooth_mm}))


def push_undo(s):
    s['undo'].append((list(s['fit']), list(s['Q']), list(s['N']), list(s['at'])))
    del s['undo'][:-30]


def undo(s):
    if s['undo']:
        s['fit'], s['Q'], s['N'], s['at'] = s['undo'].pop()


def fit_arc(s):
    """Arc position of every fit point along the curve, and the length."""
    sarr, total = arc_positions(s['Q'])
    return [sarr[j] for j in s['at']], total


def begin_drag(s, k):
    push_undo(s)
    arc, total = fit_arc(s) if len(s['fit']) >= 3 else ([0.0] * len(s['fit']), 1.0)
    s['drag'] = {'k': k, 'fit0': list(s['fit']), 'arc': arc, 'total': total}
    s['hot'] = k


def deform(s, target, influence):
    """Fit point k to target (on the bone); with an influence distance, the fit
    points within it along the curve follow with a cosine fall-off."""
    dr = s['drag']
    k, F0 = dr['k'], dr['fit0']
    delta = _sub(target, F0[k])
    fit = list(F0)
    fit[k] = tuple(target)
    if influence > 0:
        for j in range(len(F0)):
            if j == k:
                continue
            d = abs(dr['arc'][j] - dr['arc'][k])
            d = min(d, dr['total'] - d)
            if d < influence:
                fit[j] = _add(F0[j], delta, 0.5 * (1.0 + math.cos(math.pi * d / influence)))
    s['fit'] = fit
    respline(s)


def insert_fit(s, j):
    """A new fit point on the curve at sample j, between the fit points around
    it; returns its index."""
    at = s['at']
    k = max(i for i in range(len(at)) if at[i] <= j) if any(a <= j for a in at) \
        else len(at) - 1
    s['fit'].insert(k + 1, s['Q'][j])
    respline(s)
    return k + 1


def add_ctrl(fit, p):
    """A clicked point joins the polygon where it lengthens it least."""
    if len(fit) < 3:
        fit.append(p)
        return len(fit) - 1
    best, at = None, len(fit)
    for i in range(len(fit)):
        a, b = fit[i], fit[(i + 1) % len(fit)]
        extra = _dist(a, p) + _dist(p, b) - _dist(a, b)
        if best is None or extra < best:
            best, at = extra, i + 1
    fit.insert(at, p)
    return at


def write_back(app, s):
    """Into the picked spline when the number of fit points is unchanged (moved
    in place, so whatever uses the spline follows); otherwise a new 3D sketch
    with a closed spline through exactly these fit points."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    spl = s.get('spline')
    if spl and spl.isValid and len(s['fit']) == len(s.get('src_idx', [])):
        sk = spl.parentSketch
        sk.isComputeDeferred = True
        moved = 0
        try:
            for i, src in enumerate(s['src_idx']):
                fp = spl.fitPoints.item(src)
                cur = sk.modelToSketchSpace(fp.worldGeometry)
                tgt = sk.modelToSketchSpace(_pt(s['fit'][i]))
                v = adsk.core.Vector3D.create(tgt.x - cur.x, tgt.y - cur.y, tgt.z - cur.z)
                if v.length > 1e-7:
                    if fp.isFixed:
                        fp.isFixed = False
                    fp.move(v)
                    moved += 1
        finally:
            sk.isComputeDeferred = False
        return {'spline': sk.name, 'fit_points': spl.fitPoints.count, 'moved': moved,
                'new': False}
    root = des.rootComponent
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'Bone curve'
    coll = adsk.core.ObjectCollection.create()
    sk.isComputeDeferred = True
    try:
        for q in s['fit']:
            coll.add(sk.modelToSketchSpace(_pt(q)))
        spl = sk.sketchCurves.sketchFittedSplines.add(coll)
        spl.isClosed = True
    finally:
        sk.isComputeDeferred = False
    spl.attributes.add(ATTR, 'role', 'bone_curve')
    return {'spline': sk.name, 'fit_points': spl.fitPoints.count, 'moved': 0, 'new': True}


# --------------------------------------------------------------------------- #
# overlay - drawn like a sketch spline: a thin line, round fit points
# --------------------------------------------------------------------------- #

def clear_graphics(app):
    try:
        des = adsk.fusion.Design.cast(app.activeProduct)
        groups = des.rootComponent.customGraphicsGroups
        for i in range(groups.count - 1, -1, -1):
            g = groups.item(i)
            try:
                if g.id == ATTR:
                    g.deleteMe()
            except Exception:
                pass
    except Exception:
        pass


def _solid(entity, rgb, weight, priority=3):
    entity.color = adsk.fusion.CustomGraphicsSolidColorEffect.create(
        adsk.core.Color.create(rgb[0], rgb[1], rgb[2], 255))
    entity.weight = weight
    entity.depthPriority = priority
    return entity


def _flat(points):
    c = []
    for q in points:
        c.extend((q[0] * MM, q[1] * MM, q[2] * MM))
    return adsk.fusion.CustomGraphicsCoordinates.create(c)


def draw(app, s):
    """The overlay, then a viewport refresh - without it the change shows only
    when the view next moves."""
    clear_graphics(app)
    if s is not None:
        des = adsk.fusion.Design.cast(app.activeProduct)
        g = des.rootComponent.customGraphicsGroups.add()
        g.id = ATTR
        Q, N = s['Q'], s['N']
        if len(N) != len(Q):
            N = [(0.0, 0.0, 0.0)] * len(Q)
        # a hair off the bone along its normal, or the bone hides half of it
        lifted = [_add(q, n, LIFT) for q, n in zip(Q, N)]
        if len(lifted) >= 2:
            pts = lifted + ([lifted[0]] if len(lifted) > 2 else [])
            idx = []
            for i in range(len(pts) - 1):
                idx.extend((i, i + 1))
            _solid(g.addLines(_flat(pts), idx, False), CURVE_RGB, CURVE_WEIGHT)
        fit = [lifted[j] if j < len(lifted) else s['fit'][k]
               for k, j in enumerate(s['at'][:len(s['fit'])])] or s['fit']
        hot = s['drag']['k'] if s.get('drag') else s.get('hot')
        plain = [q for k, q in enumerate(fit) if k != hot]
        user = adsk.fusion.CustomGraphicsPointTypes.UserDefinedCustomGraphicsPointType
        if plain:
            ps = g.addPointSet(_flat(plain), [], user, POINT_IMG)
            ps.depthPriority = 4
        if hot is not None and hot < len(fit):
            ps = g.addPointSet(_flat([fit[hot]]), [], user, HOT_IMG)
            ps.depthPriority = 5
    try:
        app.activeViewport.refresh()
    except Exception:
        pass


# --------------------------------------------------------------------------- #
# the base (Guide Base)
# --------------------------------------------------------------------------- #

def base_job(bone, curves, holes, p, tubes=()):
    polys = [curve_polyline(c) for c in curves]
    P, closed = chain(polys)
    if not closed:
        raise ValueError('The outline must be closed (its ends are %.2f mm apart).'
                         % _dist(P[0], P[-1]))
    paths = write_bone(bone)
    os.makedirs(WORK, exist_ok=True)
    job = {'mode': 'base', 'bone_coords': paths['coords'], 'bone_tris': paths['tris'],
           'curve': resample(P, 0.1), 'clearance': p['clearance'],
           'thickness': p['thickness'], 'round': p['round'],
           'round_contact': p['round_contact'], 'voxel': p['voxel_mm'],
           'smooth': p['curve_smooth'], 'out_stl': os.path.join(WORK, 'base.stl'),
           'holes': [{'origin': o, 'axis': a, 'radius': r} for o, a, r in holes],
           'tubes': [{'origin': o, 'axis': a, 'kind': kind,
                      'radius': p['hole_d'] / 2.0 if p['hole_d'] > 0 else r,
                      'R': p['guide_d'] / 2.0, 'height': p['guide_h'],
                      'clear': p['guide_clear'],
                      'top_round': p['guide_top_round'] if kind == 'cyl' else 0.0,
                      'edge_round': p['guide_edge_round'], 'spin': p['guide_spin']}
                     for (o, a, r), kind in tubes],
           'decimate_tol': p['facet_tol'] if OUTPUTS[int(p['output'])][0] == 'solid'
           else 0.0}
    return job


def import_mesh(des, comp, stl, name, summary):
    parametric = des.designType == adsk.fusion.DesignTypes.ParametricDesignType
    base = None
    if parametric:
        base = comp.features.baseFeatures.add()
        base.name = name
        base.startEdit()
    try:
        made = comp.meshBodies.add(stl, adsk.fusion.MeshUnits.MillimeterMeshUnit, base) \
            if base else comp.meshBodies.add(stl, adsk.fusion.MeshUnits.MillimeterMeshUnit)
        body = made.item(0)
        body.name = name
        body.attributes.add(ATTR, 'role', 'guide_base')
        body.attributes.add(ATTR, 'params', json.dumps(summary))
    finally:
        if base:
            base.finishEdit()
    if base:
        body = comp.meshBodies.item(comp.meshBodies.count - 1)
    return body


def to_solid(des, comp, mesh_body, name):
    """The mesh body as a faceted solid (the only conversion this licence
    allows), the mesh body hidden. Returns the solid."""
    ci = comp.features.meshConvertFeatures.createInput([mesh_body])
    ci.meshConvertMethodType = adsk.fusion.MeshConvertMethodTypes.FacetedMeshConvertMethodType
    ci.meshConvertOperationType = adsk.fusion.MeshConvertOperationTypes \
        .ParametricFeatureMeshConvertOperationType
    feat = comp.features.meshConvertFeatures.add(ci)
    made = [feat.bodies.item(i) for i in range(feat.bodies.count)] if feat else []
    solid = max(made, key=lambda b: b.volume) if made else None
    if solid is None or not solid.isSolid:
        raise RuntimeError('Convert Mesh did not give a solid.')
    if len(made) > 1:
        log('Convert Mesh gave %d bodies; the largest is kept' % len(made))
    solid.name = name
    solid.attributes.add(ATTR, 'role', 'guide_base')
    mesh_body.isLightBulbOn = False
    return solid


def build_base(app, ui, bone, curves, holes, p, tubes=()):
    des = adsk.fusion.Design.cast(app.activeProduct)
    res = run_engine(base_job(bone, curves, holes, p, tubes), ui)
    name = 'Guide base t%.2f c%.2f r%.2f' % (p['thickness'], p['clearance'], p['round'])
    body = import_mesh(des, des.rootComponent, os.path.join(WORK, 'base.stl'), name,
                       {'version': VERSION, 'params': p, 'result': res})
    if OUTPUTS[int(p['output'])][0] == 'solid':
        t0 = time.time()
        body = to_solid(des, des.rootComponent, body, name)
        res['solid'] = {'faces': body.faces.count, 'volume_mm3': body.volume * CM ** 3,
                        'secs': time.time() - t0}
    return body, res


def guide_from(app, bone, curves, tubes, holes, p):
    """Guide Base, B-Rep: the engine fits the slab from the bone's triangles
    (written once per session), Fusion assembles and cuts it. The full bone is
    loaded into the engine only when there are tubes to stand on it."""
    # the base needs the bone's triangles (written once per session); tubes on a
    # solid bone need nothing but Fusion, on a mesh bone the engine's ray
    paths = write_mesh(bone) if curves else None
    mesh_bone = adsk.fusion.MeshBody.cast(bone) is not None
    if tubes and mesh_bone:
        ask = server_for(bone).ask
    else:
        ask = engine().ask
    body, rep = _brep().build_guide(app, bone, curves, tubes, holes, p, ask, paths, log)
    log('guide base: %s' % json.dumps({k: v for k, v in rep.items() if k != 'fit'}))
    log('guide base fit: %s' % json.dumps(rep.get('fit')))
    return body, rep


def brep_report(p, rep):
    f = rep.get('fit') or {}
    if rep.get('tubes_only'):
        txt = ('No outline picked: %s only, %s, %.0f mm3, in %.0f s.'
               % (rep.get('body', 'guide tubes'), 'all solid' if rep.get('solid') else
                  'NOT all solid', rep.get('volume_mm3', 0.0), rep.get('secs', 0.0)))
        for r in rep.get('tubes') or []:
            if r.get('error'):
                txt += '\n  %s: %s' % (r['name'], r['error'])
            else:
                txt += '\n  %s: %s, axis %.0f deg off the bone normal%s' % (
                    r['name'], 'pentagon' if r['kind'] == 'pent' else 'cylinder',
                    r['tilt_deg'], '' if r['round'] == 'ok' else ', round NOT built (%s)'
                    % r['round']) + ('' if str(r.get('gap', 'ok')).startswith('ok') else
                                     ', gap NOT applied (%s)' % r.get('gap'))
        if rep.get('tubes_note'):
            txt += '\n  ' + rep['tubes_note']
        if rep.get('holes'):
            txt += '\n%d holes cut.' % rep['holes']
        return txt
    gap = rep.get('bone_gap') or {}
    txt = ('Guide base: %s, %d faces, %.0f mm3, in %.0f s.\n'
           'Built from one surface fitted to the bone (within %.3f mm of it), cut by the '
           'outline, thickened %.2f mm minus %.2f mm.' % (
               'solid' if rep.get('solid') else 'NOT solid', rep.get('faces', 0),
               rep.get('volume_mm3', 0.0), rep.get('secs', 0.0), f.get('dev_max_mm', 0.0),
               p['clearance'] + p['thickness'] + rep.get('clearance_added_mm', 0.0),
               p['clearance'] + rep.get('clearance_added_mm', 0.0)))
    if gap:
        txt += ('\nMeasured: base %.3f mm from the bone at the closest (clearance %.2f), '
                'overlap %.3f mm3.' % (gap['min_mm'], p['clearance'], gap['overlap_mm3']))
    for key, label in (('rim', 'Rim round'), ('contact', 'Contact edge round')):
        rr = rep.get(key)
        if rr:
            txt += ('\n%s %.2f mm: %s' % (label, rr.get('radius_mm', 0.0), 'done' if rr.get('ok')
                                           else 'NOT built - %s' % rr.get('msg', '').strip()))
    if rep.get('attempts', 1) > 1:
        txt += ('\nThe rim round needed %d tries (the outline restarted on another straight '
                'stretch each time).' % rep['attempts'])
    if f.get('open_radius_mm'):
        txt += ('\nConvex bends of the outline tighter than %.2f mm (tightest %.2f) were opened '
                'for the rim round: the outline moved %.2f mm at most there.'
                % (1.25 * p['round'], f.get('rim_tightest_before_mm', 0.0),
                   f.get('eased_moved_max_mm', 0.0)))
    rows = rep.get('tubes') or []
    if rows:
        txt += '\n\n%d guide tubes:' % len(rows)
        for r in rows:
            if r.get('error'):
                txt += '\n  %s: %s' % (r['name'], r['error'])
            else:
                txt += '\n  %s: %s, axis %.0f deg off the bone normal%s' % (
                    r['name'], 'pentagon' if r['kind'] == 'pent' else 'cylinder',
                    r['tilt_deg'], '' if r['round'] == 'ok' else ', round NOT built (%s)'
                    % r['round']) + (
                    (', ' + str(r['gap']).split(', ', 2)[-1]) if str(r.get('gap', 'ok')).startswith('ok')
                    else ', gap NOT applied (%s)' % r.get('gap'))
        if rep.get('tubes_note'):
            txt += '\n  ' + rep['tubes_note']
    if rep.get('holes'):
        txt += '\n%d holes cut.' % rep['holes']
    txt += ('\n\nTime: fit %.0f s, split %.0f s, thicken %.0f s, rounds %.0f s, tubes %.0f s.'
            % (rep.get('engine_s', 0), rep.get('split_s', 0), rep.get('thicken_s', 0),
               rep.get('rim_s', 0), rep.get('tubes_s', 0)))
    return txt


def base_report(p, res):
    txt = ('Guide base: %d triangles, %.0f mm3, footprint %.0f mm2.\n'
           'Contact face %.3f mm from the bone (asked %.2f), base top up to %.2f mm.\n'
           'Outline %.1f mm, %d holes. %.0f s.'
           % (res['triangles'], res['volume_mm3'], res['footprint_mm2'],
              res['contact_median_mm'] or 0.0, p['clearance'], res['top_max_mm'],
              res['curve_mm'], res.get('holes', 0), res['secs']))
    tubes = res.get('tubes') or []
    if tubes:
        hr = [t['hole_min_r'] for t in tubes if t['hole_min_r'] is not None]
        tops = [t['top_mm'] for t in tubes if t['top_mm'] is not None]
        gaps = [t['skin_to_bone_min_mm'] for t in tubes]
        txt += ('\n\n%d guides (%d cylindrical, %d pentagonal), measured on the result:\n'
                'narrowest hole diameter %.3f mm, top %.2f-%.2f mm above the bone, '
                'tube skin at least %.2f mm off the bone (asked %.2f).'
                % (len(tubes), sum(t['kind'] == 'cyl' for t in tubes),
                   sum(t['kind'] == 'pent' for t in tubes), 2 * min(hr) if hr else 0.0,
                   min(tops) if tops else 0.0, max(tops) if tops else 0.0, min(gaps),
                   p['guide_clear']))
    if res.get('dropped'):
        dr = res['dropped']
        txt += ('\n\n%d loose specks left out (%.3f mm3 in all, none over %.2f mm).'
                % (dr['pieces'], dr['volume_mm3'], dr['largest_extent_mm']))
    dec = res.get('decimation')
    if dec:
        txt += ('\n\nFaceted solid from %d -> %d triangles, within %.3f mm of the '
                'contoured surface (asked %.3f).' % (dec['triangles_before'], dec['triangles'],
                                                     dec['hausdorff_mm'], dec['tol_mm']))
    if res.get('solid'):
        txt += ' %d faces, converted in %.0f s.' % (res['solid']['faces'], res['solid']['secs'])
    return txt


# --------------------------------------------------------------------------- #
# dialogs
# --------------------------------------------------------------------------- #

def _sel(inputs, iid):
    it = inputs.itemById(iid)
    return [it.selection(i).entity for i in range(it.selectionCount)] if it else []


def _cmd_id(cmd):
    try:
        return cmd.parentCommandDefinition.id
    except Exception:
        return '?'


def on_created(app, cmd):
    cid = _cmd_id(cmd)
    o = load_settings()
    inputs = cmd.commandInputs
    rl = adsk.core.ValueInput.createByReal
    cmd.okButtonText = 'OK'
    if cid == CMD_CURVE:
        cmd.setDialogInitialSize(380, 560)
        b = inputs.addSelectionInput('bone', 'Bone', 'The bone the curve lies on')
        b.addSelectionFilter('SolidBodies')
        b.addSelectionFilter('MeshBodies')
        # 0, not 1: editing empties the boxes (the session holds the picks), and
        # a box short of its minimum greys OK whatever validate says
        b.setSelectionLimits(0, 1)
        sh = inputs.addSelectionInput('shapes', 'Shapes',
                                      'Drill cylinders and other bodies, or an alpha '
                                      'shape: the curve is where they meet the bone.')
        sh.addSelectionFilter('SolidBodies')
        sh.addSelectionFilter('MeshBodies')
        sh.setSelectionLimits(0, 0)
        vi = inputs.addValueInput('margin', 'Margin', 'mm', rl(float(o['margin']) * MM))
        vi.tooltip = ('Grows every shape in Shapes by this much before the closing: the '
                      'cannula and the wall round it on a drill cylinder. 0 for a picked '
                      'alpha shape.')
        sh2 = inputs.addSelectionInput('shapes2', 'Shapes, own margin',
                                       'Bodies that take the margin below instead - a '
                                       'cutting guide, say. A body picked in both boxes '
                                       'takes this one.')
        sh2.addSelectionFilter('SolidBodies')
        sh2.addSelectionFilter('MeshBodies')
        sh2.setSelectionLimits(0, 0)
        vi = inputs.addValueInput('margin2', 'Their margin', 'mm',
                                  rl(float(o['margin2']) * MM))
        vi.tooltip = 'The margin of the bodies in Shapes, own margin.'
        vi = inputs.addValueInput('closing', 'Closing radius', 'mm',
                                  rl(float(o['closing']) * MM))
        vi.tooltip = ('The alpha radius: gaps narrower than twice this between the '
                      'shapes are filled. 0 cuts the shapes with the bone as they are.')
        cm = inputs.addDropDownCommandInput('closing_mode', 'Closing',
                                            adsk.core.DropDownStyles.TextListDropDownStyle)
        for i, (_k, lab) in enumerate(CLOSINGS):
            cm.listItems.add(lab, i == int(o['closing_mode']))
        cm.tooltip = ('In planes across the drills: section by section, as an alpha '
                      'shape lofted through plane sections is built - side tabs that '
                      'wrap round the bone edge stay with the rest. 3D: a closing by a '
                      'ball, the same in every direction.')
        jn = inputs.addBoolValueInput('join', 'Join every shape', True, '',
                                      bool(o['join']))
        jn.tooltip = ('A shape further from the others than the closing bridges ends '
                      'up on its own island, which the curve would leave out. On: each '
                      'such island is joined to the rest by a strip along the shortest '
                      'path over the bone, as wide as the closing. A drill that also '
                      'cuts the far side of the bone is not bridged there.')
        c = inputs.addSelectionInput('curve', 'Or a curve',
                                     'A closed sketch curve to edit instead. With no '
                                     'shapes and no curve, clicks on the bone draw one.')
        c.addSelectionFilter('SketchCurves')
        c.setSelectionLimits(0, 0)
        inputs.addBoolValueInput('edit', 'Edit curve', True, '', False)
        inputs.itemById('edit').tooltip = (
            'On: the curve is made (from the shapes, the curve, or nothing) and the left '
            'mouse button edits it like a sketch spline that stays on the bone - drag a '
            'point, click the curve to add one, Shift+click a point to remove it. The '
            'selection boxes are emptied while editing; Off puts the picks back.')
        inputs.addValueInput('influence', 'Influence distance', 'mm',
                             rl(float(o['influence']) * MM))
        inputs.itemById('influence').tooltip = (
            '0: only the dragged point moves, as in a sketch spline. Above 0: the fit '
            'points within this distance along the curve follow it, less and less.')
        inputs.addValueInput('fit_spacing', 'Fit point spacing', 'mm',
                             rl(float(o['fit_spacing']) * MM))
        inputs.itemById('fit_spacing').tooltip = (
            'Spacing of the fit points on a curve made from the shapes or from sketch '
            'curves that are not one fitted spline.')
        inputs.addValueInput('smooth_mm', 'Smoothing', 'mm', rl(float(o['smooth_mm']) * MM))
        inputs.addBoolValueInput('smooth', 'Smooth the curve', False, '', False)
        inputs.addBoolValueInput('undo', 'Undo last change', False, '', False)
        info = inputs.addTextBoxCommandInput('info', '', '', 7, True)
        info.isFullWidth = True
        st = state()
        st.session = None
        st.command = cmd
        _wire_mouse(cmd)
        refresh_curve(app, cmd)
    else:
        cmd.setDialogInitialSize(380, 640)
        b = inputs.addSelectionInput('bone', 'Bone', 'The bone the guide sits on: a solid '
                                     'or Form body (a mesh: convert it first, Mesh to Part)')
        b.addSelectionFilter('SolidBodies')
        b.setSelectionLimits(1, 1)
        c = inputs.addSelectionInput('curve', 'Outline', 'The closed curve - a 3D sketch '
                                     'on or near the bone; it is projected onto the bone '
                                     'along the bone\'s normals')
        c.addSelectionFilter('SketchCurves')
        c.setSelectionLimits(0, 0)
        h = inputs.addSelectionInput('holes', 'Holes', 'Cylindrical faces or cylinder '
                                     'bodies: a hole is cut along each, at its radius')
        h.addSelectionFilter('CylindricalFaces')
        h.addSelectionFilter('SolidBodies')
        h.setSelectionLimits(0, 0)
        for iid, lab, tip in (
                ('tubes_cyl', 'Cylindrical guides', 'Drill or wire cylinders (faces or '
                 'bodies): a round guide tube stands on the bone along each, with the '
                 'hole through it.'),
                ('tubes_pent', 'Pentagonal guides', 'Drill or wire cylinders: a '
                 'pentagonal guide tube along each, with the hole through it.')):
            sel = inputs.addSelectionInput(iid, lab, tip)
            sel.addSelectionFilter('CylindricalFaces')
            sel.addSelectionFilter('SolidBodies')
            sel.setSelectionLimits(0, 0)
        g = inputs.addGroupCommandInput('guides', 'Guide tubes')
        g.isExpanded = False
        gi = g.children
        # defaults measured on the hand-made GUIA CORTE of a distal femur case: round
        # tubes Ø7.0 with a Ø3.5 hole and a 1.5 mm top round; pentagons inscribed
        # in Ø6.0 with a Ø2.0 hole, their five top edges rounded 1.0, long edges
        # sharp; tops a median 12.1 / 10.5 mm above the bone along the axis
        for iid, lab, tip in (
                ('guide_d', 'Round: outside diameter', 'Outside diameter of a round tube.'),
                ('hole_d', 'Round: hole diameter', '0: the picked cylinder\'s own diameter.'),
                ('guide_h', 'Round: height above the bone', 'From where the axis enters the '
                 'bone, along the axis, to the top.'),
                ('guide_top_round', 'Round: top round', 'Round on the top rim.'),
                ('pent_d', 'Pentagon: outside diameter', 'The circle the pentagon is inscribed '
                 'in (6.0 mm: flats 4.85 mm apart).'),
                ('pent_hole_d', 'Pentagon: hole diameter', '0: the picked cylinder\'s own '
                 'diameter.'),
                ('pent_h', 'Pentagon: height above the bone', 'From where the axis enters the '
                 'bone, along the axis, to the top.'),
                ('pent_top_round', 'Pentagon: top edge round', 'Round on the five top edges; '
                 'the long edges stay sharp.'),
                ('guide_clear', 'Gap to the bone', 'The tube\'s bottom stays this far off '
                 'the bone.')):
            vi = gi.addValueInput(iid, lab, 'mm', rl(float(o[iid]) * MM))
            vi.tooltip = tip
        sp = gi.addValueInput('guide_spin', 'Pentagon rotation', 'deg',
                              rl(math.radians(float(o['guide_spin']))))
        sp.tooltip = ('Turn about the axis. 0: a corner points along the guide\'s length '
                      '(the outline\'s long direction).')
        for iid, lab, tip in (
                ('clearance', 'Clearance', 'Gap between the bone and the base: the surface '
                 'is offset this far off the bone before it is thickened.'),
                ('thickness', 'Thickness', 'Of the base, from the offset surface outward.'),
                ('round', 'Rim round', 'Round on the outer (top) rim; 0 for none. Convex '
                 'bends of the outline tighter than 1.25 x this are opened first.'),
                ('round_contact', 'Contact edge round', 'Round on the edge at the bone; 0 '
                 'leaves it sharp.')):
            vi = inputs.addValueInput(iid, lab, 'mm', rl(float(o[iid]) * MM))
            vi.tooltip = tip
        info = inputs.addTextBoxCommandInput('info', '', 'The bone round the outline as one '
                                             'surface, cut by the outline, thickened by '
                                             'clearance + thickness minus the clearance, '
                                             'then rounded; tubes and holes as solid '
                                             'bodies. No mesh.', 4, True)
        info.isFullWidth = True


def base_params(inputs):
    p = {k: inputs.itemById(k).value * CM for k in ('clearance', 'thickness', 'round',
                                                     'round_contact')}
    p.update({k: inputs.itemById(k).value * CM for k in
              ('guide_d', 'guide_h', 'guide_clear', 'hole_d', 'guide_top_round',
               'pent_d', 'pent_hole_d', 'pent_h', 'pent_top_round')})
    p['guide_spin'] = math.degrees(inputs.itemById('guide_spin').value)
    return p


def picked_tubes(inputs):
    """(cylinder, kind, name, near) for every guide picked - near: the centre
    of the picked cylindrical face(s), mm, which says which end of the axis the
    guide is on when there is no outline. A cylinder picked both ways stays the
    round one."""
    tubes = []
    for iid, kind in (('tubes_cyl', 'cyl'), ('tubes_pent', 'pent')):
        for e in _sel(inputs, iid):
            name = getattr(e, 'name', '') or getattr(getattr(e, 'body', None), 'name', '')
            for c in cylinders([e]):
                if not any(_same_axis(c, t[0]) for t in tubes):
                    tubes.append((c, kind, name, cylinder_centre(e)))
    return tubes


def cylinder_centre(e):
    """Area centre of the cylindrical face(s) of a picked face or body, mm."""
    f = adsk.fusion.BRepFace.cast(e)
    faces = [f] if f else [x for x in adsk.fusion.BRepBody.cast(e).faces
                           if adsk.core.Cylinder.cast(x.geometry)]
    if not faces:
        return None
    tot = sum(x.area for x in faces)
    return tuple(sum(_mm(x.centroid)[k] * x.area for x in faces) / tot for k in range(3))


def _same_axis(c1, c2):
    (o1, a1, r1), (o2, a2, r2) = c1, c2
    q = _sub(o1, o2)
    cr = (q[1] * a2[2] - q[2] * a2[1], q[2] * a2[0] - q[0] * a2[2], q[0] * a2[1] - q[1] * a2[0])
    return abs(abs(sum(a1[k] * a2[k] for k in range(3))) - 1) < 1e-6 and \
        math.sqrt(sum(x * x for x in cr)) < 1e-3


def curve_params(inputs):
    p = {k: inputs.itemById(k).value * CM for k in
         ('margin', 'margin2', 'closing', 'influence', 'fit_spacing', 'smooth_mm')}
    p['join'] = bool(inputs.itemById('join').value)
    it = inputs.itemById('closing_mode').selectedItem
    p['closing_mode'] = it.index if it else 0
    p['curve_smooth'] = load_settings()['curve_smooth']
    return p


def refresh_curve(app, cmd):
    """Dialog state of Bone Curve, then the overlay - drawn LAST, because writing
    any command input wipes every custom graphics group."""
    inputs = cmd.commandInputs
    st = state()
    s = st.session
    edit = inputs.itemById('edit').value
    for iid in ('bone', 'shapes', 'shapes2', 'curve', 'margin', 'margin2', 'closing',
                'closing_mode', 'join', 'fit_spacing'):
        inputs.itemById(iid).isEnabled = not edit
    for iid in ('smooth', 'undo', 'smooth_mm'):
        inputs.itemById(iid).isEnabled = bool(edit and s)
    if s is None:
        txt = ('Pick the bone, then the shapes (drill cylinders, other bodies, or an '
               'alpha shape) or a closed sketch curve - or neither, to draw one. Then '
               'tick <b>Edit curve</b>.')
    elif len(s['fit']) < 3 and not s['curves'] and not (s['shapes'] or s.get('shapes2')):
        txt = ('<b>Draw:</b> %d points. Click on the bone to add one (it goes where it '
               'lengthens the outline least); drag a point to move it.' % len(s['fit']))
    else:
        _, total = arc_positions(s['Q']) if len(s['Q']) > 1 else (None, 0.0)
        txt = ('<b>Edit:</b> %d fit points, curve %.1f mm, every point on the bone.<br/>'
               'Drag a point. Click the curve to add a point, Shift+click a point to '
               'remove it.' % (len(s['fit']), total))
    if s:
        sec = s['stats'].get('section')
        if sec:
            others = sec['loops'] - 1
            txt += ('<br/>From the shapes (%s closing): %.0f mm loop%s, in %.0f s.'
                    % (sec['closing'], sec['loop_mm'][0],
                       '; %d other loop%s on the bone left out (far side or islands)'
                       % (others, '' if others == 1 else 's') if others else '',
                       sec.get('secs', 0)))
            txt += shapes_note(sec.get('shapes') or [])
        if s['stats'].get('last'):
            txt += '<br/>%s' % _esc(s['stats']['last'])
    inputs.itemById('info').formattedText = txt
    draw(app, s)


def shapes_note(rows):
    """The shapes the curve does not simply hold, by name."""
    lines = []
    for r in rows:
        nm = _esc(r.get('name') or '?')
        if r['status'] == 'joined':
            lines.append('%s: on its own island, joined by a strip %.1f mm long.'
                         % (nm, r.get('bridge_mm', 0.0)))
        elif r['status'] == 'out':
            lines.append('<b>%s: on its own island, left out</b> (tick Join every shape, '
                         'or raise the closing radius).' % nm)
        elif r['status'] == 'off':
            gap = r.get('gap_mm')
            lines.append('<b>%s: does not reach the bone</b> (%s even with its margin).'
                         % (nm, '%.1f mm short' % gap if gap is not None
                            else 'more than 5 mm away'))
    return ''.join('<br/>' + ln for ln in lines)


def start_session(app, cmd):
    inputs = cmd.commandInputs
    st = state()
    bone = _sel(inputs, 'bone')
    if not bone:
        inputs.itemById('edit').value = False
        inputs.itemById('info').formattedText = '<b>Pick the bone first.</b>'
        return
    p = curve_params(inputs)
    save_settings(p)
    try:
        st.session = new_session(bone[0], _sel(inputs, 'curve'), _sel(inputs, 'shapes'), p,
                                 _sel(inputs, 'shapes2'))
    except Exception as exc:
        log(traceback.format_exc())
        inputs.itemById('edit').value = False
        st.session = None
        inputs.itemById('info').formattedText = '<b>%s</b>' % _esc(
            exc if isinstance(exc, ValueError) else short_error(traceback.format_exc()))
        return
    # the picks stay in their boxes (locked while editing): emptying them read as
    # the tool dropping the bodies


def on_input_changed(app, cmd, iid):
    cid = _cmd_id(cmd)
    inputs = cmd.commandInputs
    if cid != CMD_CURVE:
        return
    st = state()
    if iid in ('bone', 'shapes', 'shapes2', 'curve'):
        # picking: once the bone is in, every further click goes to Shapes, and
        # Shapes keeps taking bodies until another box is clicked. Nothing else
        # is rewritten mid-pick.
        if iid == 'bone' and _sel(inputs, 'bone'):
            inputs.itemById('shapes').hasFocus = True
        return
    if iid == 'edit':
        if inputs.itemById('edit').value:
            start_session(app, cmd)
        else:
            # back to picking: put the picks back in their boxes
            s = st.session
            if s and not _sel(inputs, 'bone'):
                inputs.itemById('bone').addSelection(s['bone'])
                for e in s['shapes']:
                    inputs.itemById('shapes').addSelection(e)
                for e in s.get('shapes2', []):
                    inputs.itemById('shapes2').addSelection(e)
                for c in s['curves']:
                    inputs.itemById('curve').addSelection(c)
    s = st.session
    if iid == 'smooth' and inputs.itemById('smooth').value:
        inputs.itemById('smooth').value = False
        if s and len(s['fit']) >= 4:
            push_undo(s)
            respline(s, inputs.itemById('smooth_mm').value * CM)
            s['stats']['last'] = 'Smoothed; Undo takes it back.'
    elif iid == 'undo' and inputs.itemById('undo').value:
        inputs.itemById('undo').value = False
        if s:
            undo(s)
    refresh_curve(app, cmd)


def on_validate(cmd):
    inputs = cmd.commandInputs
    if _cmd_id(cmd) == CMD_CURVE:
        s = state().session
        return bool(s and len(s['fit']) >= 3)
    # a base needs the outline; tubes alone need only the cylinders
    return bool(_sel(inputs, 'bone') and (_sel(inputs, 'curve') or _sel(inputs, 'tubes_cyl')
                                          or _sel(inputs, 'tubes_pent')))


def on_execute(app, ui, cmd):
    inputs = cmd.commandInputs
    if _cmd_id(cmd) == CMD_CURVE:
        s = state().session
        save_settings(curve_params(inputs))
        if s and len(s['fit']) >= 3:
            out = write_back(app, s)
            log('curve written: %s' % out)
        return
    p = base_params(inputs)
    save_settings(p)
    bone = _sel(inputs, 'bone')[0]
    holes = cylinders(_sel(inputs, 'holes'))
    tubes = picked_tubes(inputs)
    try:
        body, res = guide_from(app, bone, _sel(inputs, 'curve'), tubes, holes, p)
    except ValueError as exc:
        ui.messageBox(str(exc), CMD_NAMES[CMD_BASE])
        return
    ui.messageBox(brep_report(p, res), CMD_NAMES[CMD_BASE])


def on_destroy(app, cmd):
    st = state()
    if _cmd_id(cmd) == CMD_CURVE and (st.command is None or
                                      _cmd_id(st.command) == CMD_CURVE):
        clear_graphics(app)
        try:
            app.activeViewport.refresh()
        except Exception:
            pass
        st.session = None
        st.command = None
        st.handlers = []


# --------------------------------------------------------------------------- #
# the mouse (Bone Curve)
# --------------------------------------------------------------------------- #

def _wire_mouse(cmd):
    st = state()
    st.handlers = []
    for h, ev in ((_MouseDown(), cmd.mouseDown), (_MouseDrag(), cmd.mouseDrag),
                  (_MouseUp(), cmd.mouseUp), (_MouseMove(), cmd.mouseMove)):
        ev.add(h)
        st.handlers.append(h)


class _MouseDown(adsk.core.MouseEventHandler):
    def notify(self, args):
        on_mouse(args, 'down')


class _MouseDrag(adsk.core.MouseEventHandler):
    def notify(self, args):
        on_mouse(args, 'drag')


class _MouseUp(adsk.core.MouseEventHandler):
    def notify(self, args):
        on_mouse(args, 'up')


class _MouseMove(adsk.core.MouseEventHandler):
    def notify(self, args):
        on_mouse(args, 'move')


def press(app, s, vp, cursor, shift, influence, hit):
    """What a left click at cursor does; hit is the bone point under it (mm) or
    None. Returns a note for the dialog, or None if nothing happened."""
    k = fit_under_cursor(vp, s, cursor)
    if k is not None and hit is not None and \
            _dist(hit, s['fit'][k]) > max(2.0, 3.0 * PICK_PX / px_per_mm(vp, hit)):
        k = None                                  # that point is behind the bone
    if k is not None:
        if shift:
            if len(s['fit']) <= 3:
                return 'A closed curve keeps at least 3 points.'
            push_undo(s)
            del s['fit'][k]
            respline(s)
            s['hot'] = None
            return 'Point removed.'
        begin_drag(s, k)
        return None
    if hit is None:
        return None
    if len(s['Q']) > 2 and len(s['fit']) >= 3:
        j, dj = nearest_index(s['Q'], hit)
        if dj is not None and dj * px_per_mm(vp, hit) <= PICK_PX:
            push_undo(s)
            k = insert_fit(s, j)
            s['drag'] = None
            begin_drag(s, k)
            s['undo'].pop()                       # one undo for insert + drag
            return 'Point added.'
        if s['curves'] or s['shapes'] or s.get('shapes2'):
            return None
    push_undo(s)
    k = add_ctrl(s['fit'], hit)
    respline(s)
    begin_drag(s, k)
    s['undo'].pop()
    return None


def on_mouse(args, phase):
    """Never raises into Fusion. The button is latched on mouseDown: mouseDrag
    does not reliably carry it."""
    try:
        st = state()
        s, cmd = st.session, st.command
        if s is None or cmd is None or not cmd.commandInputs.itemById('edit').value:
            return
        app = adsk.core.Application.get()
        inputs = cmd.commandInputs
        vp = args.viewport
        cursor = args.viewportPosition
        influence = inputs.itemById('influence').value * CM
        if phase == 'down':
            if not (args.button & adsk.core.MouseButtons.LeftMouseButton):
                return
            o, d = ray_from_view(vp, cursor)
            hit = hit_bone(app, s['bone'], o, d)
            shift = bool(args.keyboardModifiers &
                         adsk.core.KeyboardModifiers.ShiftKeyboardModifier)
            note = press(app, s, vp, cursor, shift, influence, hit)
            if note:
                s['stats']['last'] = note
                refresh_curve(app, cmd)
            else:
                draw(app, s)
        elif phase in ('drag', 'move'):
            if s.get('drag'):
                o, d = ray_from_view(vp, cursor)
                hit = hit_bone(app, s['bone'], o, d)
                if hit is not None:
                    deform(s, hit, influence)
                    draw(app, s)
            elif phase == 'move':
                k = fit_under_cursor(vp, s, cursor)
                if k != s.get('hot'):
                    s['hot'] = k
                    draw(app, s)
        elif phase == 'up':
            if s.get('drag'):
                s['drag'] = None
                s['stats']['last'] = 'Moved; Undo takes it back.'
                refresh_curve(app, cmd)
    except Exception:
        log(traceback.format_exc())


# --------------------------------------------------------------------------- #
# installation
# --------------------------------------------------------------------------- #

_TOOLTIPS = {
    CMD_CURVE: ('A closed curve on the bone, edited by dragging it over the bone.\n\n'
                'Drag any point of the curve: it follows the cursor on the bone, and '
                'its neighbours within the influence distance follow, still on the '
                'bone. With no curve picked, clicks on the bone draw a new one.'),
    CMD_BASE: ('The guide base inside a closed curve on the bone.\n\n'
               'Contact face at the clearance, top at clearance + thickness, both '
               'measured along the bone normal; the wall stands along the bone normal '
               'all round, so side tabs keep their thickness; rim rounded by a field '
               'blend that cannot fail; holes along the picked cylinders.'),
}


def install(app, ui, created_handler_factory):
    placed = []
    for cmd_id, name in CMD_NAMES.items():
        old = ui.commandDefinitions.itemById(cmd_id)
        if old:
            old.deleteMe()
        folder = os.path.join(RES_DIR, 'curve' if cmd_id == CMD_CURVE else 'base')
        cd = ui.commandDefinitions.addButtonDefinition(
            cmd_id, name, _TOOLTIPS[cmd_id], folder if os.path.isdir(folder) else '')
        cd.commandCreated.add(created_handler_factory())
        for pid in PANELS:
            panel = ui.allToolbarPanels.itemById(pid)
            if panel and not panel.controls.itemById(cmd_id):
                panel.controls.addCommand(cd)
                placed.append('%s -> %s' % (name, pid))
    return '; '.join(placed)


def uninstall(ui):
    for cmd_id in CMD_NAMES:
        for pid in PANELS:
            panel = ui.allToolbarPanels.itemById(pid)
            c = panel.controls.itemById(cmd_id) if panel else None
            if c:
                c.deleteMe()
        cd = ui.commandDefinitions.itemById(cmd_id)
        if cd:
            cd.deleteMe()


# --------------------------------------------------------------------------- #
# dev bridge - the same code paths, driven without a mouse
# --------------------------------------------------------------------------- #

def _design(app):
    return adsk.fusion.Design.cast(app.activeProduct)


def _all_bodies(app):
    des = _design(app)
    out = []
    for comp in des.allComponents:
        out.extend(comp.bRepBodies.item(i) for i in range(comp.bRepBodies.count))
        out.extend(comp.meshBodies.item(i) for i in range(comp.meshBodies.count))
    return out


def _body(app, name):
    for b in _all_bodies(app):
        if b.name == name:
            return b
    raise ValueError('no body %r' % name)


def _sketch_curves(app, sketch_name):
    root = _design(app).rootComponent
    sk = root.sketches.itemByName(sketch_name)
    if sk is None:
        raise ValueError('no sketch %r' % sketch_name)
    cs = sk.sketchCurves
    return [cs.item(i) for i in range(cs.count)]


def dev_ping(app, job):
    des = _design(app)
    return {'doc': app.activeDocument.name if app.activeDocument else None,
            'version': VERSION, 'parametric': des.designType ==
            adsk.fusion.DesignTypes.ParametricDesignType if des else None}


def dev_new_test_doc(app, job):
    """A NEW document with a bone STL as a mesh body and a curve as a closed
    fitted spline - the only kind of document the dev commands modify."""
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    des = _design(app)
    root = des.rootComponent
    base = root.features.baseFeatures.add()
    base.startEdit()
    try:
        mb = root.meshBodies.add(job['bone_stl'], adsk.fusion.MeshUnits.MillimeterMeshUnit,
                                 base).item(0)
        mb.name = 'Bone'
    finally:
        base.finishEdit()
    sk = root.sketches.add(root.xYConstructionPlane)
    sk.name = 'Outline'
    coll = adsk.core.ObjectCollection.create()
    sk.isComputeDeferred = True
    for p in job['curve']:
        coll.add(sk.modelToSketchSpace(_pt(p)))
    spl = sk.sketchCurves.sketchFittedSplines.add(coll)
    spl.isClosed = True
    sk.isComputeDeferred = False
    for cyl in job.get('cylinders', []):
        o, a, r, L = cyl['origin'], cyl['axis'], cyl['radius'], cyl.get('length', 40.0)
        tbm = adsk.fusion.TemporaryBRepManager.get()
        p0 = _pt(_add(o, a, -L / 2))
        p1 = _pt(_add(o, a, L / 2))
        tmp = tbm.createCylinderOrCone(p0, r * MM, p1, r * MM)
        bf = root.features.baseFeatures.add()
        bf.startEdit()
        try:
            b = root.bRepBodies.add(tmp, bf)
            b.name = cyl.get('name', 'Drill')
        finally:
            bf.finishEdit()
    return {'doc': doc.name, 'mesh_triangles': mb.mesh.triangleCount,
            'fit_points': spl.fitPoints.count}


def dev_close_doc(app, job):
    doc = app.activeDocument
    name = doc.name
    if not name.startswith('Untitled'):
        return {'refused': 'only an unsaved test document is closed: %s' % name}
    doc.close(False)
    return {'closed': name}


def dev_session(app, job):
    """Open an editing session exactly as the dialog does: from shapes (body
    names), a sketch, or nothing."""
    st = state()
    curves = _sketch_curves(app, job['sketch']) if job.get('sketch') else None
    shapes = [_body(app, n) for n in job.get('shapes', [])]
    shapes2 = [_body(app, n) for n in job.get('shapes2', [])]
    p = dict(load_settings())
    p.update(job.get('params', {}))
    t0 = time.time()
    st.session = new_session(_body(app, job['bone']), curves, shapes, p, shapes2)
    s = st.session
    _, total = arc_positions(s['Q']) if len(s['Q']) > 1 else (None, 0)
    return {'fit': len(s['fit']), 'samples': len(s['Q']), 'length_mm': total,
            'section': s['stats'].get('section'),
            'attached_max_mm': s['stats'].get('attached_max_mm'),
            'secs': round(time.time() - t0, 2)}


def dev_drag(app, job):
    """A drag without a mouse: grab the fit point nearest `grab` (mm), pull it
    towards `to` in `steps` events, each target projected onto the bone first
    (what the ray hit does), then release."""
    s = state().session
    srv = server_for(s['bone'])
    k, dgrab = nearest_index(s['fit'], tuple(job['grab']))
    fit0 = list(s['fit'])
    begin_drag(s, k)
    t0 = time.time()
    n = int(job.get('steps', 10))
    a, b = s['fit'][k], tuple(job['to'])
    for i in range(1, n + 1):
        tgt = srv.ask({'op': 'attach', 'pts': [_add(a, _sub(b, a), i / n)]})['pts'][0]
        deform(s, tuple(tgt), float(job.get('influence', 0.0)))
    per_event = (time.time() - t0) / n
    s['drag'] = None
    chk = srv.ask({'op': 'attach', 'pts': s['Q']})['pts']
    off = max(_dist(p, q) for p, q in zip(chk, s['Q']))
    moved = [j for j in range(len(fit0)) if _dist(fit0[j], s['fit'][j]) > 1e-6]
    return {'grabbed': k, 'grab_dist_mm': dgrab, 'moved_fit_points': moved,
            'dragged_to_target_mm': _dist(s['fit'][k], tgt),
            'ms_per_event': round(1000 * per_event, 1), 'max_off_bone_mm': off,
            'fit': len(s['fit'])}


def dev_click(app, job):
    """press() at the screen position of a model point (mm), as a click there
    does: the dev path to adding or removing a point."""
    s = state().session
    vp = app.activeViewport
    p = tuple(job['at'])
    xy = vp.modelToViewSpace(_pt(p))
    hit = server_for(s['bone']).ask({'op': 'attach', 'pts': [p]})['pts'][0]
    n0 = len(s['fit'])
    note = press(app, s, vp, xy, bool(job.get('shift')), 0.0, tuple(hit))
    s['drag'] = None
    draw(app, s)
    return {'note': note, 'fit_before': n0, 'fit_after': len(s['fit']),
            'px_per_mm': px_per_mm(vp, p)}


def dev_selftest(app, job):
    """The whole Bone Curve path in ONE job, in a scratch document closed at
    the end (the previous document comes back): bone + drill margins -> curve
    from the shapes -> drag, click-add, shift-click-remove -> overlay (captured)
    -> write-back read again. Refused unless Fusion is idle, so a command the
    user has open is never ended."""
    ui = app.userInterface
    if ui.activeCommand != 'SelectCommand':
        return {'refused': 'Fusion is busy (%s)' % ui.activeCommand}
    prev = app.activeDocument
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    st = state()
    out = {}
    try:
        des = _design(app)
        root = des.rootComponent
        base = root.features.baseFeatures.add()
        base.startEdit()
        try:
            bone = root.meshBodies.add(job['bone_stl'],
                                       adsk.fusion.MeshUnits.MillimeterMeshUnit, base).item(0)
        finally:
            base.finishEdit()
        bone = root.meshBodies.item(root.meshBodies.count - 1)
        tbm = adsk.fusion.TemporaryBRepManager.get()
        shapes = []
        bf = root.features.baseFeatures.add()
        bf.startEdit()
        try:
            for c in job['cylinders']:
                tmp = tbm.createCylinderOrCone(_pt(c['p0']), c['r'] * MM, _pt(c['p1']),
                                               c['r'] * MM)
                shapes.append(root.bRepBodies.add(tmp, bf))
        finally:
            bf.finishEdit()
        shapes = [root.bRepBodies.item(i) for i in range(root.bRepBodies.count)]
        specs = [as_cylinder(b) for b in shapes]
        out['cylinders_recognised'] = sum(1 for x in specs if x)
        out['cylinder_axis_err_mm'] = max(
            max(_dist(x[0], tuple(c['p0'])), _dist(x[1], tuple(c['p1'])), abs(x[2] - c['r']))
            for x, c in zip(specs, job['cylinders']) if x)
        p = dict(load_settings())
        p.update(job.get('params', {}))
        t0 = time.time()
        st.session = s = new_session(bone, None, shapes, p)
        out['session'] = {'secs': round(time.time() - t0, 1), 'fit': len(s['fit']),
                          'section': s['stats']['section']}
        srv = server_for(bone)
        # the camera on the curve, so screen picks are meaningful
        vp = app.activeViewport
        cam = vp.camera
        c0 = [sum(q[i] for q in s['Q']) / len(s['Q']) for i in range(3)]
        nn = [sum(nv[i] for nv in s['N']) / len(s['N']) for i in range(3)]
        cam.isSmoothTransition = False
        cam.cameraType = adsk.core.CameraTypes.OrthographicCameraType
        cam.target = _pt(c0)
        cam.eye = _pt(_add(c0, nn, 150.0 / (math.sqrt(sum(x * x for x in nn)) or 1)))
        cam.upVector = adsk.core.Vector3D.create(0, 0, 1)
        cam.viewExtents = float(job.get('extent', 70.0)) * MM
        vp.camera = cam
        adsk.doEvents()
        # a drag of one fit point 3 mm across the curve, by the mouse's own path:
        # screen point -> ray -> bone hit -> deform -> overlay redrawn
        k = len(s['fit']) // 3
        Q = s['fit']
        tang = _sub(Q[(k + 1) % len(Q)], Q[k - 1])
        n = s['N'][s['at'][k]]
        side = (tang[1] * n[2] - tang[2] * n[1], tang[2] * n[0] - tang[0] * n[2],
                tang[0] * n[1] - tang[1] * n[0])
        ln = math.sqrt(sum(x * x for x in side)) or 1.0
        to = _add(Q[k], side, 3.0 / ln)
        fit0 = list(s['fit'])
        begin_drag(s, k)
        t0 = time.time()
        misses = 0
        for i in range(1, 11):
            xy = vp.modelToViewSpace(_pt(_add(fit0[k], _sub(to, fit0[k]), i / 10)))
            o, d = ray_from_view(vp, xy)
            hit = hit_bone(app, bone, o, d)
            if hit is None:
                misses += 1
                continue
            deform(s, hit, 0.0)
            draw(app, s)
        per_event = (time.time() - t0) / 10
        s['drag'] = None
        chk = srv.ask({'op': 'attach', 'pts': s['Q']})['pts']
        out['drag'] = {'ms_per_event': round(1000 * per_event, 1), 'ray_misses': misses,
                       'moved': [j for j in range(len(fit0)) if _dist(fit0[j], s['fit'][j]) > 1e-6],
                       'moved_mm': round(_dist(fit0[k], s['fit'][k]), 3),
                       'curve_off_bone_max_mm': max(_dist(a, b) for a, b in zip(chk, s['Q']))}
        # hovering: the fit point under the cursor, at 40 cursor positions
        t0 = time.time()
        for i in range(40):
            xy = vp.modelToViewSpace(_pt(s['fit'][i % len(s['fit'])]))
            fit_under_cursor(vp, s, xy)
        out['hover_ms_per_event'] = round(25 * (time.time() - t0), 2)
        # click on the curve halfway between two fit points: a point is added
        j = (s['at'][5] + s['at'][6]) // 2
        r1 = dev_click(app, {'at': s['Q'][j]})
        # shift+click on a fit point: removed
        r2 = dev_click(app, {'at': s['fit'][10], 'shift': True})
        out['click_add'], out['shift_click_remove'] = r1, r2
        s['hot'] = 3
        draw(app, s)
        groups = root.customGraphicsGroups
        g = [groups.item(i) for i in range(groups.count) if groups.item(i).id == ATTR]
        out['overlay'] = [type(g[0].item(i)).__name__ for i in range(g[0].count)] if g else None
        vp.refresh()
        out['capture'] = bool(vp.saveAsImageFile(job['png'], 1100, 950))
        out['write'] = dev_write(app, {})
    finally:
        clear_graphics(app)
        st.session = None
        doc.close(False)
        if prev:
            prev.activate()
    return out


def dev_draw(app, job):
    """Draw mode without a mouse: points (mm, projected onto the bone first)."""
    st = state()
    st.session = new_session(_body(app, job['bone']))
    s = st.session
    srv = server_for(s['bone'])
    for p in job['points']:
        q = srv.ask({'op': 'attach', 'pts': [p]})['pts'][0]
        add_ctrl(s['fit'], tuple(q))
        respline(s)
    _, total = arc_positions(s['Q'])
    return {'fit': len(s['fit']), 'samples': len(s['Q']), 'length_mm': total}


def dev_write(app, job):
    s = state().session
    out = write_back(app, s)
    # the spline that came out, read back from Fusion: still on the bone, and
    # through the session's fit points?
    root = _design(app).rootComponent
    sk = root.sketches.itemByName(out['spline'])
    spl = sk.sketchCurves.sketchFittedSplines.item(sk.sketchCurves.sketchFittedSplines.count - 1)
    pts = curve_polyline(spl, 0.01)
    srv = server_for(s['bone'])
    sub = pts[::max(1, len(pts) // 2000)]
    chk = srv.ask({'op': 'attach', 'pts': sub})['pts']
    off = sorted(_dist(p, q) for p, q in zip(chk, sub))
    dev = sorted(min(_dist(q, p) for p in s['Q']) for q in sub)
    out.update({'spline_off_bone_max_mm': off[-1],
                'spline_off_bone_p95_mm': off[int(0.95 * (len(off) - 1))],
                'fusion_vs_editor_max_mm': dev[-1]})
    return out


def dev_base(app, job):
    """Guide Base through build_base, as OK in the dialog runs it."""
    bone = _body(app, job['bone'])
    curves = _sketch_curves(app, job['sketch'])
    if job.get('spline_only'):
        curves = [c for c in curves if adsk.fusion.SketchFittedSpline.cast(c)][-1:]
    holes = cylinders([_body(app, n) for n in job.get('holes', [])])
    tubes = [(c, 'cyl') for c in cylinders([_body(app, n) for n in job.get('tubes_cyl', [])])]
    tubes += [(c, 'pent') for c in cylinders([_body(app, n) for n in job.get('tubes_pent', [])])]
    p = dict(DEFAULTS)
    p.update({'voxel': 1})
    p.update(job.get('params', {}))
    p['voxel_mm'] = float(VOXELS[int(p['voxel'])][0])
    body, res = build_base(app, None, bone, curves, holes, p, tubes)
    mb = adsk.fusion.MeshBody.cast(body)
    res['body'] = {'name': body.name, 'triangles': mb.mesh.triangleCount if mb else None,
                   'solid': bool(adsk.fusion.BRepBody.cast(body))}
    if mb:
        # MeshBody.volume is not trusted; the signed volume of its own triangles is
        c, t = mb.mesh.nodeCoordinatesAsDouble, mb.mesh.triangleNodeIndices
        v6 = 0.0
        for k in range(0, len(t), 3):
            a, b2, c2 = (c[3 * t[k + j]:3 * t[k + j] + 3] for j in range(3))
            v6 += (a[0] * (b2[1] * c2[2] - b2[2] * c2[1]) - a[1] * (b2[0] * c2[2] - b2[2] * c2[0])
                   + a[2] * (b2[0] * c2[1] - b2[1] * c2[0]))
        res['body']['signed_volume_mm3'] = v6 / 6.0 * CM ** 3
    res['holes'] = len(holes)
    res['report'] = base_report(p, res)
    return res


def dev_base_selftest(app, job):
    """Guide Base end to end in a scratch document, closed afterwards: the bone
    STL, the outline as a fitted spline, guide cylinders as bodies, OK as the
    dialog runs it, the result measured and captured. Refused unless Fusion is
    idle."""
    ui = app.userInterface
    if ui.activeCommand != 'SelectCommand':
        return {'refused': 'Fusion is busy (%s)' % ui.activeCommand}
    prev = app.activeDocument
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    out = {}
    try:
        des = _design(app)
        root = des.rootComponent
        base = root.features.baseFeatures.add()
        base.startEdit()
        try:
            root.meshBodies.add(job['bone_stl'], adsk.fusion.MeshUnits.MillimeterMeshUnit, base)
        finally:
            base.finishEdit()
        bone = root.meshBodies.item(root.meshBodies.count - 1)
        sk = root.sketches.add(root.xYConstructionPlane)
        coll = adsk.core.ObjectCollection.create()
        for q in job['curve']:
            coll.add(_pt(q))
        spl = sk.sketchCurves.sketchFittedSplines.add(coll)
        spl.isClosed = True
        tbm = adsk.fusion.TemporaryBRepManager.get()
        bf = root.features.baseFeatures.add()
        bf.startEdit()
        try:
            for c in job['tubes_cyl'] + job['tubes_pent']:
                root.bRepBodies.add(tbm.createCylinderOrCone(
                    _pt(c['p0']), c['r'] * MM, _pt(c['p1']), c['r'] * MM), bf)
        finally:
            bf.finishEdit()
        bodies = [root.bRepBodies.item(i) for i in range(root.bRepBodies.count)]
        ncyl = len(job['tubes_cyl'])
        tubes = [(c, 'cyl') for c in cylinders(bodies[:ncyl])]
        tubes += [(c, 'pent') for c in cylinders(bodies[ncyl:])]
        for b in bodies:
            b.isLightBulbOn = False
        p = dict(DEFAULTS)
        p.update(job.get('params', {}))
        p['voxel_mm'] = float(VOXELS[int(p['voxel'])][0])
        t0 = time.time()
        body, res = build_base(app, None, bone, [spl], [], p, tubes)
        out['secs_total'] = round(time.time() - t0, 1)
        brep = adsk.fusion.BRepBody.cast(body)
        out['body'] = {'name': body.name, 'solid': bool(brep and brep.isSolid),
                       'faces': brep.faces.count if brep else None,
                       'volume_mm3': brep.volume * CM ** 3 if brep else None}
        out['engine_volume_mm3'] = res['volume_mm3']
        out['solid'] = res.get('solid')
        out['decimation'] = res.get('decimation')
        out['tubes'] = res.get('tubes')
        out['report'] = base_report(p, res)
        bone.isLightBulbOn = True
        vp = app.activeViewport
        cam = vp.camera
        cam.isSmoothTransition = False
        cam.viewOrientation = adsk.core.ViewOrientations.IsoTopRightViewOrientation
        vp.camera = cam
        vp.fit()
        adsk.doEvents()
        if job.get('eye'):
            cam = vp.camera
            cam.isSmoothTransition = False
            cam.target = _pt(job['target'])
            cam.eye = _pt(job['eye'])
            cam.upVector = adsk.core.Vector3D.create(*job.get('up', (0, 0, 1)))
            if job.get('extent'):
                cam.viewExtents = float(job['extent']) * MM
            vp.camera = cam
            adsk.doEvents()
        vp.refresh()
        out['capture'] = bool(vp.saveAsImageFile(job['png'], 1100, 950))
    finally:
        doc.close(False)
        if prev:
            prev.activate()
    return out


def dev_convert_probe(app, job):
    """Scratch document: import an STL as a mesh body, Convert Mesh (faceted) with
    the given operation type, and list every body that exists afterwards."""
    ui = app.userInterface
    if ui.activeCommand != 'SelectCommand':
        return {'refused': 'Fusion is busy (%s)' % ui.activeCommand}
    prev = app.activeDocument
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    out = {}
    try:
        des = _design(app)
        root = des.rootComponent
        mb = import_mesh(des, root, job['stl'], 'probe', {})
        out['mesh'] = {'triangles': mb.mesh.triangleCount, 'closed': mb.isClosed}
        t0 = time.time()
        ci = root.features.meshConvertFeatures.createInput([mb])
        ci.meshConvertMethodType = adsk.fusion.MeshConvertMethodTypes.FacetedMeshConvertMethodType
        ci.meshConvertOperationType = (
            adsk.fusion.MeshConvertOperationTypes.BaseFeatureMeshConvertOperationType
            if job.get('op') == 'base' else
            adsk.fusion.MeshConvertOperationTypes.ParametricFeatureMeshConvertOperationType)
        feat = root.features.meshConvertFeatures.add(ci)
        out['secs'] = round(time.time() - t0, 1)
        out['feature_bodies'] = [(b.name, b.faces.count, b.volume * CM ** 3, b.isSolid)
                                 for b in feat.bodies] if feat else None
        out['all_bodies'] = [(b.name, b.faces.count, b.volume * CM ** 3, b.isSolid)
                             for b in root.bRepBodies]
    except Exception:
        out['error'] = traceback.format_exc()[-800:]
    finally:
        doc.close(False)
        if prev:
            prev.activate()
    return out


def dev_bodies(app, job):
    return [{'name': b.name, 'mesh': bool(adsk.fusion.MeshBody.cast(b)),
             'visible': b.isVisible} for b in _all_bodies(app)]


def dev_export_mesh(app, job):
    """Read-only: a named body's triangles (mm, float64 / int32) to job['path']_v/_f.bin."""
    b = _body(app, job['name'])
    if b is None:
        return {'error': 'no body %r' % job['name']}
    c, t = body_triangles(b)
    with open(job['path'] + '_v.bin', 'wb') as fh:
        array.array('d', [v * CM for v in c]).tofile(fh)
    with open(job['path'] + '_f.bin', 'wb') as fh:
        array.array('i', t).tofile(fh)
    return {'nodes': len(c) // 3, 'tris': len(t) // 3}


def dev_capture(app, job):
    """Viewport image looking at `target` (mm) from `direction`, `extent` mm wide."""
    vp = app.activeViewport
    cam = vp.camera
    t = job['target']
    d = job['direction']
    n = math.sqrt(sum(x * x for x in d))
    dist = float(job.get('distance', 200.0))
    cam.isSmoothTransition = False
    cam.cameraType = adsk.core.CameraTypes.OrthographicCameraType
    cam.target = _pt(t)
    cam.eye = _pt([t[k] + d[k] / n * dist for k in range(3)])
    cam.upVector = adsk.core.Vector3D.create(*job.get('up', (0, 0, 1)))
    cam.viewExtents = float(job.get('extent', 60.0)) * MM
    vp.camera = cam
    adsk.doEvents()
    vp.refresh()
    path = job['path']
    return {'saved': bool(vp.saveAsImageFile(path, int(job.get('width', 1000)),
                                             int(job.get('height', 900)))), 'path': path}


def dev_visible(app, job):
    """Show only the named bodies and sketches (test documents)."""
    keep = set(job.get('bodies', []))
    for b in _all_bodies(app):
        b.isLightBulbOn = b.name in keep
    root = _design(app).rootComponent
    for i in range(root.sketches.count):
        sk = root.sketches.item(i)
        sk.isLightBulbOn = sk.name in set(job.get('sketches', []))
    return {'shown': sorted(keep)}


def dev_open(app, job):
    """Open a command's dialog as its button does (it appears after this job)."""
    ui = app.userInterface
    ui.commandDefinitions.itemById(job['cmd_id']).execute()
    return {'opened': job['cmd_id']}


def dev_dialog(app, job):
    """Read the open Bone Curve dialog; optionally pick the bone and curve into it
    and turn on editing, through the same handler a click runs."""
    st = state()
    cmd = st.command
    if cmd is None:
        return {'error': 'no Bone Curve dialog open'}
    inputs = cmd.commandInputs
    if job.get('pick'):
        inputs.itemById('bone').addSelection(_body(app, job['bone']))
        for n in job.get('shapes', []):
            inputs.itemById('shapes').addSelection(_body(app, n))
        for n in job.get('shapes2', []):
            inputs.itemById('shapes2').addSelection(_body(app, n))
        if job.get('sketch'):
            for c in _sketch_curves(app, job['sketch']):
                inputs.itemById('curve').addSelection(c)
        for iid, v in job.get('set', {}).items():
            inputs.itemById(iid).value = float(v) * MM
        inputs.itemById('edit').value = True
        on_input_changed(app, cmd, 'edit')
    if job.get('unedit'):
        inputs.itemById('edit').value = False
        on_input_changed(app, cmd, 'edit')
    out = {}
    for i in range(inputs.count):
        it = inputs.item(i)
        v = None
        for attr in ('value', 'selectionCount', 'formattedText'):
            try:
                v = getattr(it, attr)
                break
            except Exception:
                continue
        out[it.id] = v if isinstance(v, (int, float, str, bool)) else str(v)
    out['session_fit'] = len(st.session['fit']) if st.session else 0
    groups = _design(app).rootComponent.customGraphicsGroups
    out['graphics'] = [(groups.item(i).id, groups.item(i).count) for i in range(groups.count)]
    out['enabled'] = {iid: inputs.itemById(iid).isEnabled for iid in ('bone', 'shapes', 'curve',
                                                                       'smooth', 'undo')}
    return out


def dev_terminate(app, job):
    app.userInterface.terminateActiveCommand()
    return {'terminated': True}


def dev_busy(app, job):
    """Read-only: the command Fusion has open ('SelectCommand' when idle)."""
    return {'active': app.userInterface.activeCommand,
            'doc': app.activeDocument.name if app.activeDocument else None}


def dev_shutdown(app, job):
    shutdown()
    return {'server': None}


def _brep():
    """bg_brep, re-imported with this module (the hot reload)."""
    if ADDIN_DIR not in sys.path:
        sys.path.insert(0, ADDIN_DIR)
    import bg_brep
    return importlib.reload(bg_brep)


def dev_brep_test(app, job):
    """Guide Base, BRep, on COPIES: the named bodies (exact BRep copies) and the
    named sketch's fitted spline are copied out of the active document into a NEW
    one, the guide is built there and measured, a picture is saved, and the new
    document is closed unsaved (left open with keep). The active document is
    only read. Refused unless Fusion is idle."""
    ui = app.userInterface
    if ui.activeCommand != 'SelectCommand':
        return {'refused': 'Fusion is busy (%s)' % ui.activeCommand}
    src = _design(app)
    prev = app.activeDocument
    tbm = adsk.fusion.TemporaryBRepManager.get()
    names = [job['bone']] + list(job.get('tubes_cyl', [])) + list(job.get('tubes_pent', []))
    copies = []
    for n in names[1 if job.get('bone_tsm') else 0:]:
        b = _body(app, n)
        if b is None or not adsk.fusion.BRepBody.cast(b):
            return {'error': 'no BRep body %r' % n}
        copies.append(tbm.copy(b))
    sk = src.rootComponent.sketches.itemByName(job['sketch'])
    spl = sk.sketchCurves.sketchFittedSplines.item(0)
    fit = [_v3(spl.fitPoints.item(i).worldGeometry) for i in range(spl.fitPoints.count)]
    closed = spl.isClosed
    if closed and len(fit) > 1 and _dist(fit[0], fit[-1]) < 1e-6:
        fit = fit[:-1]
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    out = {'source': prev.name, 'fit_points': len(fit)}
    try:
        des = _design(app)
        if des.designType != adsk.fusion.DesignTypes.ParametricDesignType:
            des.designType = adsk.fusion.DesignTypes.ParametricDesignType
        root = des.rootComponent
        base = root.features.baseFeatures.add()
        base.startEdit()
        try:
            for c in copies:
                root.bRepBodies.add(c, base)
        finally:
            base.finishEdit()
        bodies = [base.bodies.item(i) for i in range(base.bodies.count)]
        if job.get('bone_tsm'):
            # the bone as a Form built from a TSM file (a Mesh to Part remesh)
            with open(job['bone_tsm'], encoding='utf-8') as fh:
                text = fh.read()
            ff = root.features.formFeatures.add()
            ff.startEdit()
            ff.tSplineBodies.addByTSMDescription(text)
            ff.finishEdit()
            bodies = [ff.bodies.item(0)] + bodies
            out['bone_faces'] = bodies[0].faces.count
        for b, n in zip(bodies, names):
            b.name = n
        sk2 = root.sketches.add(root.xYConstructionPlane)
        sk2.name = job['sketch']
        coll = adsk.core.ObjectCollection.create()
        for q in fit:
            coll.add(sk2.modelToSketchSpace(adsk.core.Point3D.create(*q)))
        curve = sk2.sketchCurves.sketchFittedSplines.add(coll)
        curve.isClosed = closed
        p = dict(load_settings())
        p.update(job.get('params', {}))
        bone = bodies[0]
        nc = len(job.get('tubes_cyl', []))
        tubes = []
        for i, b in enumerate(bodies[1:]):
            cyl = cylinders([b])
            if cyl:
                tubes.append((cyl[0], 'cyl' if i < nc else 'pent', b.name, cylinder_centre(b)))
        t0 = time.time()
        guide, rep = guide_from(app, bone, [] if job.get('no_outline') else [curve], tubes,
                                [], p)
        out['report'] = rep
        out['secs'] = round(time.time() - t0, 1)
        for b in bodies[1:]:
            b.isLightBulbOn = False
        if job.get('capture'):
            if job['capture'].get('hide_bone'):
                bone.isLightBulbOn = False
                sk2.isVisible = False
            bb = guide.boundingBox
            c = [(getattr(bb.minPoint, k) + getattr(bb.maxPoint, k)) * 0.5 * CM for k in 'xyz']
            ext = max((getattr(bb.maxPoint, k) - getattr(bb.minPoint, k)) * CM for k in 'xyz')
            for i, d in enumerate(job['capture'].get('directions', [[1, 1, 1]])):
                path = job['capture']['path'] % i
                dev_capture(app, {'target': c, 'direction': d, 'extent': 1.3 * ext,
                                  'distance': 400.0, 'path': path,
                                  'up': job['capture'].get('up', (0, 0, 1))})
            out['captures'] = len(job['capture'].get('directions', [[1, 1, 1]]))
    except Exception:
        out['error'] = short_error(traceback.format_exc())
        log(traceback.format_exc())
        try:
            root = _design(app).rootComponent
            out['bodies'] = []
            for b in root.bRepBodies:
                bb = b.boundingBox
                out['bodies'].append([b.name, b.isSolid, b.faces.count,
                                      [round(v * CM, 2) for v in (bb.minPoint.x, bb.minPoint.y,
                                                                  bb.minPoint.z, bb.maxPoint.x,
                                                                  bb.maxPoint.y, bb.maxPoint.z)],
                                      round(b.volume * 1000.0, 1) if b.isSolid else None])
        except Exception:
            out['bodies_error'] = short_error(traceback.format_exc())
    finally:
        if not job.get('keep'):
            doc.close(False)
            if prev:
                prev.activate()
    return out


def dev_probe_split(app, job):
    """In a scratch (Untitled) document only: distance between the slab and the
    curtain, and which split call takes them."""
    doc = app.activeDocument
    if not doc.name.startswith('Untitled'):
        return {'refused': doc.name}
    comp = _design(app).rootComponent
    slab = comp.bRepBodies.itemByName('Guide base slab')
    band = comp.bRepBodies.itemByName('Guide base curtain')
    out = {}
    try:
        m = app.measureManager.measureMinimumDistance(slab, band)
        out['min_distance_mm'] = m.value * CM
    except Exception:
        out['distance_error'] = short_error(traceback.format_exc())
    # where the band crosses the slab's faces, per face
    tbm = adsk.fusion.TemporaryBRepManager.get()
    for i, f in enumerate(slab.faces):
        try:
            a = tbm.copy(f)
            b = tbm.copy(band)
            ok = tbm.booleanOperation(a, b, adsk.fusion.BooleanTypes.IntersectionBooleanType)
            out.setdefault('face_cross', []).append([i, f.geometry.surfaceType, bool(ok),
                                                     a.edges.count if a else None])
        except Exception:
            out.setdefault('face_cross', []).append([i, 'error', short_error(traceback.format_exc())])
    for label, tool, ext in (('body', band, False), ('face', band.faces.item(0), False),
                             ('body_ext', band, True)):
        try:
            sin = comp.features.splitBodyFeatures.createInput(slab, tool, ext)
            sf = comp.features.splitBodyFeatures.add(sin)
            out[label] = 'ok, %d bodies' % sf.bodies.count
            sf.deleteMe()
        except Exception:
            out[label] = short_error(traceback.format_exc())
    return out


def dev_probe_fillet(app, job):
    """In a scratch (Untitled) document only: which rim fillets build on the
    'Guide base' body. Each attempt is undone."""
    doc = app.activeDocument
    if not doc.name.startswith('Untitled'):
        return {'refused': doc.name}
    comp = _design(app).rootComponent
    g = comp.bRepBodies.itemByName('Guide base')
    bone = comp.bRepBodies.item(0)
    faces = []
    for f in g.faces:
        d = app.measureManager.measureMinimumDistance(f, bone).value * CM
        faces.append((d, f))
    out = {'faces': [[round(d, 3), f.geometry.surfaceType, f.edges.count,
                      round(f.area * 100, 1)] for d, f in faces]}
    out['edges'] = [[round(e.length * CM, 2), e.faces.count] for e in g.edges]
    top = max(faces, key=lambda t: t[0] if t[1].edges.count < 6 else -1)[1]
    def attempt(edges, r, chain):
        coll = adsk.core.ObjectCollection.create()
        for e in edges:
            coll.add(e)
        try:
            fin = comp.features.filletFeatures.createInput()
            fin.edgeSetInputs.addConstantRadiusEdgeSet(coll, adsk.core.ValueInput.createByReal(r * MM), chain)
            f = comp.features.filletFeatures.add(fin)
            f.deleteMe()
            return 'ok'
        except Exception as exc:
            return str(exc).splitlines()[-1][:90]
    tedges = [e for e in top.edges]
    for r in (1.0, 0.5, 0.25):
        out['rim_r%.2f_chain' % r] = attempt(tedges, r, True)
    out['per_edge_r1'] = [attempt([e], 1.0, False) for e in tedges]
    out['per_edge_r0.5'] = [attempt([e], 0.5, False) for e in tedges]
    return out


def dev_probe_rim(app, job):
    """Scratch (Untitled) documents only: the loops and edges of the Guide base's
    faces, and which single-edge rounds build. Each attempt is undone."""
    doc = app.activeDocument
    if not doc.name.startswith('Untitled'):
        return {'refused': doc.name}
    comp = _design(app).rootComponent
    g = comp.bRepBodies.itemByName('Guide base')
    out = {'faces': []}
    for fi, f in enumerate(g.faces):
        fr = {'i': fi, 'area': round(f.area * 100, 1), 'loops': []}
        for lp in f.loops:
            fr['loops'].append({'outer': lp.isOuter,
                                'edges': [[round(e.length * CM, 3), e.geometry.curveType,
                                           e.isTolerant, round(e.tolerance * CM, 5)]
                                          for e in lp.edges]})
        out['faces'].append(fr)

    def attempt(e, r):
        coll = adsk.core.ObjectCollection.create()
        coll.add(e)
        try:
            fin = comp.features.filletFeatures.createInput()
            fin.edgeSetInputs.addConstantRadiusEdgeSet(coll, adsk.core.ValueInput.createByReal(r * MM), False)
            ft = comp.features.filletFeatures.add(fin)
            ft.deleteMe()
            return 'ok'
        except Exception as exc:
            return str(exc).splitlines()[-1].strip()[:70]
    out['edges'] = []
    for ei, e in enumerate(g.edges):
        out['edges'].append([ei, round(e.length * CM, 3), [x.tempId for x in e.faces],
                             attempt(e, 0.2), attempt(e, 2.0) if e.length * CM > 3 else '-'])
    return out


def dev_probe_cyl(app, job):
    """Scratch (Untitled) documents only: cut 'Guide base slab' with an open
    cylinder standing through it at job['centre'] (mm) along job['axis'], then
    try 0.2 and 2 mm rounds on the cut edges."""
    doc = app.activeDocument
    if not doc.name.startswith('Untitled'):
        return {'refused': doc.name}
    comp = _design(app).rootComponent
    slab = comp.bRepBodies.itemByName('Guide base slab')
    tbm = adsk.fusion.TemporaryBRepManager.get()
    c = [v * MM for v in job['centre']]
    a = job['axis']
    r = float(job.get('radius', 12.0)) * MM
    h = float(job.get('half', 8.0)) * MM
    cyl = tbm.createCylinderOrCone(adsk.core.Point3D.create(c[0] - a[0] * h, c[1] - a[1] * h, c[2] - a[2] * h), r,
                                   adsk.core.Point3D.create(c[0] + a[0] * h, c[1] + a[1] * h, c[2] + a[2] * h), r)
    caps = [f for f in cyl.faces if adsk.core.Plane.cast(f.geometry)]
    tbm.deleteFaces(caps, False)
    base = comp.features.baseFeatures.add()
    base.startEdit()
    comp.bRepBodies.add(cyl, base)
    base.finishEdit()
    tool = base.bodies.item(0)
    out = {}
    before = set(b.entityToken for b in comp.bRepBodies)
    try:
        sin = comp.features.splitBodyFeatures.createInput(slab, tool, bool(job.get('extended')))
        sf = comp.features.splitBodyFeatures.add(sin)
        out['split'] = 'ok'
    except Exception:
        out['split'] = short_error(traceback.format_exc())
        return out
    pieces = [b for b in comp.bRepBodies if b.isSolid]
    small = min(pieces, key=lambda b: b.volume)
    out['piece_volume_mm3'] = round(small.volume * 1000, 1)
    out['piece_faces'] = [f.geometry.surfaceType for f in small.faces]
    cut = [e for e in small.edges if any(f.geometry.surfaceType == 1 for f in e.faces)]
    out['cut_edges'] = len(cut)
    for rr in (0.2, 2.0):
        res = []
        for e in cut:
            coll = adsk.core.ObjectCollection.create()
            coll.add(e)
            try:
                fin = comp.features.filletFeatures.createInput()
                fin.edgeSetInputs.addConstantRadiusEdgeSet(coll, adsk.core.ValueInput.createByReal(rr * MM), False)
                ft = comp.features.filletFeatures.add(fin)
                ft.deleteMe()
                res.append('ok')
            except Exception as exc:
                res.append(str(exc).splitlines()[-1].strip()[:60])
        out['r%.1f' % rr] = res
    return out


def dev_probe_thicken(app, job):
    """Scratch (Untitled) documents only: the slab's contact face alone as a
    surface, thickened by job['t'] mm away from the bone, then its edges rounded."""
    doc = app.activeDocument
    if not doc.name.startswith('Untitled'):
        return {'refused': doc.name}
    comp = _design(app).rootComponent
    slab = comp.bRepBodies.itemByName('Guide base')
    bone = comp.bRepBodies.item(0)
    tbm = adsk.fusion.TemporaryBRepManager.get()
    dist = [(app.measureManager.measureMinimumDistance(f, bone).value, i) for i, f in enumerate(slab.faces)]
    bi = min(dist)[1]
    sheet = tbm.copy(slab)
    tbm.deleteFaces([f for i, f in enumerate(sheet.faces) if i != bi], False)
    base = comp.features.baseFeatures.add()
    base.startEdit()
    comp.bRepBodies.add(sheet, base)
    base.finishEdit()
    sb = base.bodies.item(0)
    out = {'sheet_faces': sb.faces.count}
    f = sb.faces.item(0)
    ok, uv = f.evaluator.getParameterAtPoint(f.pointOnFace)
    ok, n = f.evaluator.getNormalAtParameter(uv)
    p = f.pointOnFace
    probe = adsk.core.Point3D.create(p.x + n.x * 0.05, p.y + n.y * 0.05, p.z + n.z * 0.05)
    away = bone.pointContainment(probe) != adsk.fusion.PointContainment.PointInsidePointContainment
    t = float(job.get('t', 3.5)) * MM * (1 if away else -1)
    faces = adsk.core.ObjectCollection.create()
    faces.add(f)
    t0 = time.time()
    try:
        tin = comp.features.thickenFeatures.createInput(faces, adsk.core.ValueInput.createByReal(t), False,
                                                         adsk.fusion.FeatureOperations.NewBodyFeatureOperation, False)
        th = comp.features.thickenFeatures.add(tin)
        out['thicken'] = 'ok %.1f s' % (time.time() - t0)
        solid = next(b for b in th.bodies if b.isSolid)
        out['solid_faces'] = solid.faces.count
        res = []
        for e in solid.edges:
            coll = adsk.core.ObjectCollection.create()
            coll.add(e)
            try:
                fin = comp.features.filletFeatures.createInput()
                fin.edgeSetInputs.addConstantRadiusEdgeSet(coll, adsk.core.ValueInput.createByReal(0.2), False)
                ft = comp.features.filletFeatures.add(fin)
                ft.deleteMe()
                res.append('ok')
            except Exception as exc:
                res.append(str(exc).splitlines()[-1].strip()[:50])
        out['edge_fillets_2mm'] = res
    except Exception:
        out['thicken'] = short_error(traceback.format_exc())
    return out


def dev_probe_gap(app, job):
    """Scratch document: a tube cut by a cylinder 'bone', its cut faces offset by
    +-gap, and the result measured (min distance, overlap volume) - which sign of
    Offset Faces makes the gap. Closed unsaved."""
    ui = app.userInterface
    if ui.activeCommand != 'SelectCommand':
        return {'refused': 'Fusion is busy (%s)' % ui.activeCommand}
    BB = _brep()
    prev = app.activeDocument
    doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    out = {}
    try:
        des = _design(app)
        des.designType = adsk.fusion.DesignTypes.ParametricDesignType
        comp = des.rootComponent
        tbm = adsk.fusion.TemporaryBRepManager.get()
        P = adsk.core.Point3D.create
        for sign in (-1.0, 1.0):
            bone_t = tbm.createCylinderOrCone(P(-5, 0, 0), 1.5, P(5, 0, 0), 1.5)
            tube_t = tbm.createCylinderOrCone(P(0.3, 0, 0.5), 0.35, P(0.3, 0, 4.0), 0.35)
            bodies, _ = BB.add_bodies(comp, [bone_t, tube_t], ['bone %+d' % sign, 'tube %+d' % sign])
            bone, tube = bodies
            BB.combine(comp, tube, [bone], adsk.fusion.FeatureOperations.CutFeatureOperation, True)
            tube = [b for b in comp.bRepBodies if b.name == 'tube %+d' % sign][0]
            faces = [f for f in tube.faces if bone.pointContainment(f.pointOnFace) ==
                     adsk.fusion.PointContainment.PointOnPointContainment]
            v0 = tube.volume
            oin = comp.features.offsetFacesFeatures.createInput(
                faces, adsk.core.ValueInput.createByReal(sign * 0.02))
            comp.features.offsetFacesFeatures.add(oin)
            tube = [b for b in comp.bRepBodies if b.name == 'tube %+d' % sign][0]
            md = app.measureManager.measureMinimumDistance(tube, bone)
            inter = tbm.copy(tube)
            try:
                tbm.booleanOperation(inter, tbm.copy(bone),
                                     adsk.fusion.BooleanTypes.IntersectionBooleanType)
                ov = inter.volume * 1000.0
            except Exception:
                ov = 0.0
            out['sign %+d' % sign] = {'faces': len(faces), 'vol_before_mm3': round(v0 * 1000, 3),
                                      'vol_after_mm3': round(tube.volume * 1000, 3),
                                      'min_dist_mm': round(md.value * CM, 4),
                                      'overlap_mm3': round(ov, 4)}
            bone.isLightBulbOn = False
            tube.isLightBulbOn = False
    except Exception:
        out['error'] = short_error(traceback.format_exc())
    finally:
        doc.close(False)
        if prev:
            prev.activate()
    return out


def dev_probe_chain(app, job):
    """Scratch (Untitled) documents only: chained rounds on the whole bottom and
    top loops of 'Guide base' (bottom = the face nearest the bone), each undone."""
    doc = app.activeDocument
    if not doc.name.startswith('Untitled'):
        return {'refused': doc.name}
    comp = _design(app).rootComponent
    g = comp.bRepBodies.itemByName('Guide base')
    bone = comp.bRepBodies.item(0)
    big = sorted(g.faces, key=lambda f: -f.area)[:3]
    d = sorted([(app.measureManager.measureMinimumDistance(f, bone).value, i) for i, f in enumerate(big)])
    bottom, top = big[d[0][1]], big[d[-1][1]]

    def chain(face_fn, r, then=None):
        f = face_fn()
        coll = adsk.core.ObjectCollection.create()
        for e in f.edges:
            coll.add(e)
        try:
            fin = comp.features.filletFeatures.createInput()
            fin.edgeSetInputs.addConstantRadiusEdgeSet(coll, adsk.core.ValueInput.createByReal(r * MM), True)
            ft = comp.features.filletFeatures.add(fin)
        except Exception as exc:
            return str(exc).splitlines()[-1].strip()[:50]
        res = 'ok'
        if then:
            res += ' + ' + then()
        ft.deleteMe()
        return res
    out = {}
    for r in job.get('bottom_r', [0.2, 0.3, 0.5, 1.0]):
        out['bottom %.1f' % r] = chain(lambda: bottom, r)
    for r in job.get('top_r', [2.0]):
        out['top %.1f' % r] = chain(lambda: top, r)
    out['top 2.0 then bottom 0.2'] = chain(lambda: top, 2.0, lambda: chain(
        lambda: min(comp.bRepBodies.itemByName('Guide base').faces,
                    key=lambda f: app.measureManager.measureMinimumDistance(f, bone).value), 0.2))
    return out


def _v3(p):
    return (p.x, p.y, p.z)


DEV_COMMANDS = {
    'brep_test': dev_brep_test,
    'probe_split': dev_probe_split,
    'probe_fillet': dev_probe_fillet,
    'probe_rim': dev_probe_rim,
    'probe_cyl': dev_probe_cyl,
    'probe_thicken': dev_probe_thicken,
    'probe_chain': dev_probe_chain,
    'probe_gap': dev_probe_gap,
    'ping': dev_ping,
    'new_test_doc': dev_new_test_doc,
    'close_doc': dev_close_doc,
    'session': dev_session,
    'drag': dev_drag,
    'click': dev_click,
    'selftest': dev_selftest,
    'busy': dev_busy,
    'draw': dev_draw,
    'write': dev_write,
    'base': dev_base,
    'base_selftest': dev_base_selftest,
    'convert_probe': dev_convert_probe,
    'bodies': dev_bodies,
    'export_mesh': dev_export_mesh,
    'shutdown': dev_shutdown,
    'capture': dev_capture,
    'visible': dev_visible,
    'open': dev_open,
    'dialog': dev_dialog,
    'terminate': dev_terminate,
}
