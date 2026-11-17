"""Guide Base as B-Rep bodies built from scratch. Imported by bg_core.

Fusion never splits, offsets or thickens the bone - on a scanned bone of 20 000
faces that is what takes minutes. Instead:

    1. the engine (bg_slab, py -3.10) samples the bone's triangles near the
       outline, projects the outline onto it along its normals, and fits a
       NURBS slab over the patch and a margin: the contact face is the bone
       offset by the clearance, the top that face offset by the thickness, the
       walls ruled - six faces sharing their edges exactly;
    2. Fusion assembles the slab (BRepBodyDefinition) and cuts it with a curtain
       standing along the bone normals through the projected outline (Outline
       Patch's band: a ruled surface between two control-point offsets of one
       spline), keeps the piece inside, and rounds the top rim;
    3. guide tubes - cylinders, or pentagonal prisms made as the intersection of
       five slabs - are built from scratch, trimmed by the bone, joined, and a
       hole is cut along every drill.

Parameters and reports in mm; geometry handed to Fusion in cm.
"""

import importlib
import importlib.util
import json
import math
import os
import sys
import time

import adsk.core
import adsk.fusion

CM = 10.0
MM = 0.1
ADDIN_DIR = os.path.dirname(os.path.realpath(__file__))


def outline_patch():
    """Outline Patch's op_core, from the add-in installed beside this one."""
    path = os.path.join(os.path.dirname(ADDIN_DIR), 'OutlinePatch', 'op_core.py')
    if not os.path.exists(path):
        raise ValueError('Guide Base uses Outline Patch for the cutting curtain: '
                         'install it (fusion-psd-suite, or its own install.ps1)')
    mod = sys.modules.get('op_core')
    if mod is not None and os.path.realpath(getattr(mod, '__file__', '')) == \
            os.path.realpath(path):
        return importlib.reload(mod)
    spec = importlib.util.spec_from_file_location('op_core', path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules['op_core'] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- #
# small helpers, cm tuples
# --------------------------------------------------------------------------- #

def _v(p):
    return (p.x, p.y, p.z)


def _pt(t):
    return adsk.core.Point3D.create(t[0], t[1], t[2])


def _cm(q):
    return (q[0] * MM, q[1] * MM, q[2] * MM)


def _vec(t):
    return adsk.core.Vector3D.create(t[0], t[1], t[2])


def _add(a, b, s=1.0):
    return (a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2])


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a, b):
    return (a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0])


def _unit(a):
    n = math.sqrt(_dot(a, a)) or 1.0
    return (a[0] / n, a[1] / n, a[2] / n)


def add_bodies(comp, tmps, names):
    """Temporary bodies into the design through one BaseFeature."""
    base = comp.features.baseFeatures.add()
    base.startEdit()
    try:
        for t in tmps:
            comp.bRepBodies.add(t, base)
    finally:
        base.finishEdit()
    out = [base.bodies.item(i) for i in range(base.bodies.count)]
    for b, n in zip(out, names):
        b.name = n
    return out, base


def fillet(comp, edges, r_cm, chain=False):
    if not edges or r_cm <= 0:
        return True, ''
    coll = adsk.core.ObjectCollection.create()
    for e in edges:
        coll.add(e)
    try:
        fin = comp.features.filletFeatures.createInput()
        fin.edgeSetInputs.addConstantRadiusEdgeSet(
            coll, adsk.core.ValueInput.createByReal(r_cm), chain)
        comp.features.filletFeatures.add(fin)
        return True, ''
    except Exception as exc:
        return False, str(exc).splitlines()[-1][:200]


def fillet_down(comp, edges_of, r_mm, chain=True, step=0.2, search=True):
    """The round at r_mm; if Fusion refuses it, the next smaller one in steps of
    `step` mm until one builds. (ok, radius built mm, first error, radii failed)."""
    tried, first = [], ''
    r = r_mm
    while r > 0.05:
        edges = edges_of()
        if not edges:
            return False, 0.0, 'edges not found', tried
        ok, msg = fillet(comp, edges, r * MM, chain)
        if ok:
            return True, r, first, tried
        tried.append(round(r, 2))
        first = first or msg
        if not search:
            break
        r = round(r - step, 6)
    return False, 0.0, first, tried


def run_tag():
    """Suffix that keeps this run's bodies apart from an earlier guide's in the
    same component (their names are how the build finds them again)."""
    return ' #%06x' % (int(time.time() * 1000) & 0xffffff)


def named(comp, name):
    return next(b for b in comp.bRepBodies if b.name == name)


def combine(comp, target, tools, op, keep):
    coll = adsk.core.ObjectCollection.create()
    for t in tools:
        coll.add(t)
    cin = comp.features.combineFeatures.createInput(target, coll)
    cin.operation = op
    cin.isKeepToolBodies = keep
    return comp.features.combineFeatures.add(cin)


def remove(comp, bodies):
    out = []
    for b in bodies:
        try:
            out.append(comp.features.removeFeatures.add(b))
        except Exception:
            pass
    return out


def inside(body, p_cm):
    return body.pointContainment(_pt(p_cm)) == \
        adsk.fusion.PointContainment.PointInsidePointContainment


# --------------------------------------------------------------------------- #
# the slab from the engine's faces
# --------------------------------------------------------------------------- #

def slab_body(topo):
    """BRepBodyDefinition: 8 vertices, 12 edges, 6 faces, each face's u x v
    outward and its loop counter-clockwise (laid out so by bg_slab)."""
    bd = adsk.fusion.BRepBodyDefinition.create()
    bd.doFullHealing = True
    lump = bd.lumpDefinitions.add()
    shell = lump.shellDefinitions.add()
    vd = {k: bd.createVertexDefinition(_pt(_cm(q))) for k, q in topo['vertices'].items()}
    ed = {}
    for k, e in topo['edges'].items():
        curve = adsk.core.NurbsCurve3D.createNonRational(
            [_pt(_cm(q)) for q in e['ctrl']], int(e['deg']), list(e['knots']), False)
        if curve is None:
            raise RuntimeError('edge %s: no curve' % k)
        ed[k] = bd.createEdgeDefinitionByCurve(vd[e['start']], vd[e['end']], curve)
    Open = adsk.core.NurbsSurfaceProperties.OpenNurbsSurface
    for f in topo['faces']:
        sp = f['surface']
        surf = adsk.core.NurbsSurface.create(
            int(sp['du']), int(sp['dv']), int(sp['nu']), int(sp['nv']),
            [_pt(_cm(q)) for q in sp['net']], list(sp['ku']), list(sp['kv']), [],
            Open, Open)
        if surf is None:
            raise RuntimeError('face %s: no surface' % f['name'])
        fd = shell.faceDefinitions.add(surf, False)
        loop = fd.loopDefinitions.add()
        for key, rev in f['loop']:
            loop.bRepCoEdgeDefinitions.add(ed[key], bool(rev))
    body = bd.createBody()
    info = list(bd.outcomeInfo or [])
    if body is None:
        raise RuntimeError('the slab would not assemble: %s' % '; '.join(info))
    return body, info


def band_face(lo, hi):
    """The cutting band as an EXPLICIT degree (3, 1) NURBS face between the two
    rails (clamped closed cubics on the same knots), assembled like the slab's
    walls: a closed face with one seam. A band from createRuledSurface left cut
    edges that no fillet would take, at any radius; the slab's own ruled walls,
    built this way, take 0.2 and 2 mm rounds."""
    a = [_pt(_cm(q)) for q in lo['ctrl']]
    b = [_pt(_cm(q)) for q in hi['ctrl']]
    ku = list(lo['knots'])
    n = len(a)
    net = []
    for i in range(n):
        net += [a[i], b[i]]
    Open = adsk.core.NurbsSurfaceProperties.OpenNurbsSurface
    surf = adsk.core.NurbsSurface.create(3, 1, n, 2, net, ku, [0.0, 0.0, 1.0, 1.0], [],
                                         Open, Open)
    if surf is None:
        raise RuntimeError('the cutting band surface would not build')
    bd = adsk.fusion.BRepBodyDefinition.create()
    bd.doFullHealing = True
    shell = bd.lumpDefinitions.add().shellDefinitions.add()
    va, vb = bd.createVertexDefinition(a[0]), bd.createVertexDefinition(b[0])
    ca = adsk.core.NurbsCurve3D.createNonRational(a, 3, ku, False)
    cb = adsk.core.NurbsCurve3D.createNonRational(b, 3, ku, False)
    ea = bd.createEdgeDefinitionByCurve(va, va, ca)
    eb = bd.createEdgeDefinitionByCurve(vb, vb, cb)
    es = bd.createEdgeDefinitionByCurve(va, vb, adsk.core.Line3D.create(a[0], b[0]))
    fd = shell.faceDefinitions.add(surf, False)
    loop = fd.loopDefinitions.add()
    for e, rev in ((ea, False), (es, False), (eb, True), (es, True)):
        loop.bRepCoEdgeDefinitions.add(e, rev)
    body = bd.createBody()
    if body is None:
        raise RuntimeError('the cutting band would not assemble: %s'
                           % '; '.join(bd.outcomeInfo or []))
    return body


def build_base(app, srv_ask, bone_paths, curves, p, rep, log, bone=None):
    """The base alone: engine slab -> Fusion solid -> curtain cut -> rim round.
    On failure every feature it made is deleted again, so a failed run leaves
    nothing behind in the design."""
    made = []
    build = _build_base if p.get('base_mode') == 'slab' else _build_base_surface
    try:
        return build(app, srv_ask, bone_paths, curves, p, rep, log, made, bone)
    except Exception:
        if p.get('_debug'):
            raise
        for f in reversed(made):
            try:
                f.deleteMe()
            except Exception:
                pass
        raise


def surface_body(r):
    """The engine's net as a one-face open surface body (u-major net; u x v along
    the bone's outward normal), its four boundary rows as the edges."""
    nu, nv = int(r['nu']), int(r['nv'])
    net = [_pt(_cm(q)) for q in r['net']]
    ku, kv = list(r['ku']), list(r['kv'])

    def at(i, j):
        return net[i * nv + j]
    Open = adsk.core.NurbsSurfaceProperties.OpenNurbsSurface
    surf = adsk.core.NurbsSurface.create(3, 3, nu, nv, net, ku, kv, [], Open, Open)
    if surf is None:
        raise RuntimeError('the bone surface would not build')
    bd = adsk.fusion.BRepBodyDefinition.create()
    bd.doFullHealing = True
    shell = bd.lumpDefinitions.add().shellDefinitions.add()
    c = [bd.createVertexDefinition(q) for q in (at(0, 0), at(nu - 1, 0),
                                                at(nu - 1, nv - 1), at(0, nv - 1))]
    rows = [([at(i, 0) for i in range(nu)], ku, 0, 1),
            ([at(nu - 1, j) for j in range(nv)], kv, 1, 2),
            ([at(i, nv - 1) for i in range(nu)], ku, 3, 2),
            ([at(0, j) for j in range(nv)], kv, 0, 3)]
    eds = [bd.createEdgeDefinitionByCurve(
        c[i0], c[i1], adsk.core.NurbsCurve3D.createNonRational(pts, 3, kn, False))
        for pts, kn, i0, i1 in rows]
    fd = shell.faceDefinitions.add(surf, False)
    loop = fd.loopDefinitions.add()
    for e, rev in zip(eds, (False, False, True, True)):
        loop.bRepCoEdgeDefinitions.add(e, rev)
    body = bd.createBody()
    if body is None:
        raise RuntimeError('the bone surface would not assemble: %s'
                           % '; '.join(bd.outcomeInfo or []))
    return body


def _face_at(comp, body, q_cm, tol=0.06):
    found = comp.findBRepUsingPoint(_pt(q_cm), adsk.fusion.BRepEntityTypes.BRepFaceEntityType,
                                    tol, False)
    for i in range(found.count):
        f = adsk.fusion.BRepFace.cast(found.item(i))
        if f and f.body.entityToken == body.entityToken:
            return f
    return None


def _build_base_surface(app, srv_ask, bone_paths, curves, p, rep, log, made, bone=None):
    """The base from ONE surface: the bone round the outline fitted as a single
    NURBS face, split by the outline itself (closest point) and the outside
    deleted; that face thickened by clearance + thickness, minus the same face
    thickened by the clearance; then the rim and contact rounds. One face in,
    so one rim edge at the top and one at the bottom."""
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.activeComponent
    raw = outline_patch().chain_outline(curves)
    pts_mm = [(q[0] * CM, q[1] * CM, q[2] * CM) for q in raw]
    c, t = float(p['clearance']), float(p['thickness'])
    t0 = time.time()
    r = srv_ask({'op': 'surface', 'coords': bone_paths['coords'], 'tris': bone_paths['tris'],
                 'pts': pts_mm, 'spacing': float(p.get('surface_spacing', 3.0)),
                 'margin': 4.0, 'height': c + t, 'open': not p.get('_no_ease'),
                 'rim': max(float(p.get('round', 0.0)), float(p.get('round_contact', 0.0)))})
    rep['engine_s'] = round(time.time() - t0, 1)
    rep['fit'] = r['report']
    log('engine surface %.1f s' % rep['engine_s'])

    # the closed rim edge's one vertex sits where the outline starts; a round
    # can fail where that lands, so a failed rim round is tried again with the
    # outline starting on another straight stretch, and on that last try a
    # radius that still fails steps down 0.2 mm at a time until one builds
    starts = (r['report'].get('starts') or [0])[:2]
    tag = run_tag()
    for attempt, start in enumerate(starts):
        mine = []
        last = attempt == len(starts) - 1
        guide = _base_attempt(app, comp, r, curves, p, rep, c, t, bone, start, mine, tag, last)
        made.extend(mine)
        rim_ok = rep.get('rim', {}).get('ok', True)
        rep['attempts'] = attempt + 1
        if rim_ok or last:
            break
        for f in reversed(mine):
            try:
                f.deleteMe()
            except Exception:
                pass
            made.remove(f)
    named(comp, 'Guide base surface' + tag).name = 'Guide base surface'
    rep['faces'] = guide.faces.count
    if adsk.fusion.BRepBody.cast(bone):
        rep['bone_gap'] = body_gap(guide, bone)
    return guide


def _base_attempt(app, comp, r, curves, p, rep, c, t, bone, start, made, tag, search):
    """Split, thicken, cut and round once, the outline starting at `start`.
    Bodies carry `tag` in their names while the guide is built."""
    t1 = time.time()
    SURF, BASE = 'Guide base surface' + tag, 'Guide base' + tag
    (sb,), bf = add_bodies(comp, [surface_body(r)], [SURF])
    made.append(bf)
    faces = adsk.core.ObjectCollection.create()
    faces.add(sb.faces.item(0))
    if r['report'].get('outline_eased'):
        # the outline with its too-tight convex bends opened for the rim round,
        # as a hidden sketch of its own; the user's sketch is not touched
        sk = comp.sketches.add(comp.xYConstructionPlane)
        sk.name = 'Guide base outline (eased)'
        made.append(sk)
        pts = adsk.core.ObjectCollection.create()
        pts_ = r['report']['outline_eased']
        for q in pts_[start:] + pts_[:start]:
            pts.add(sk.modelToSketchSpace(_pt(_cm(q))))
        sk.isComputeDeferred = True
        spl = sk.sketchCurves.sketchFittedSplines.add(pts)
        spl.isClosed = True
        sk.isComputeDeferred = False
        tool = spl
    elif len(curves) == 1:
        tool = curves[0]
    else:
        tool = adsk.core.ObjectCollection.create()
        for cv in curves:
            tool.add(cv)
    sin = comp.features.splitFaceFeatures.createInput(faces, tool, False)
    sin.setClosestPointSplitType()
    made.append(comp.features.splitFaceFeatures.add(sin))
    sb = named(comp, SURF)
    rep['split_faces'] = sb.faces.count
    outer = [f for f in sb.faces if f.loops.count > 1]
    if sb.faces.count != 2 or len(outer) != 1:
        raise RuntimeError('the outline did not split the surface into an inside and an '
                           'outside (%d faces)' % sb.faces.count)
    coll = adsk.core.ObjectCollection.create()
    coll.add(outer[0])
    made.append(comp.features.surfaceDeleteFaceFeatures.add(coll))
    sb = named(comp, SURF)
    face = sb.faces.item(0)
    rep['split_s'] = round(time.time() - t1, 1)

    # thicken away from the bone, whichever way the face points: a solid bone
    # says which side of the face is inside it; a mesh bone, the engine's normal
    q = _cm(r['report']['inside'])
    ok, uv = face.evaluator.getParameterAtPoint(_pt(q))
    ok, nf = face.evaluator.getNormalAtParameter(uv)
    nf = _unit(_v(nf))
    if adsk.fusion.BRepBody.cast(bone):
        IN = adsk.fusion.PointContainment.PointInsidePointContainment
        reads = {h: bone.pointContainment(_pt(_add(q, nf, h * MM))) == IN
                 for h in (-2.0, -1.0, 1.0, 2.0)}
        score = sum(reads[h] for h in (-2.0, -1.0)) - sum(reads[h] for h in (1.0, 2.0))
        rep['side_reads'] = {str(h): v for h, v in reads.items()}
        if score == 0:
            out_side = _dot(nf, _unit(tuple(r['report']['inside_normal']))) > 0
            rep['side_from'] = 'engine normal (containment undecided)'
        else:
            out_side = score > 0
            rep['side_from'] = 'bone containment'
    else:
        out_side = _dot(nf, _unit(tuple(r['report']['inside_normal']))) > 0
        rep['side_from'] = 'engine normal'
    sign = 1.0 if out_side else -1.0
    nb = nf if out_side else (-nf[0], -nf[1], -nf[2])

    t1 = time.time()

    def thicken(d_mm, name):
        before = set(b.entityToken for b in comp.bRepBodies)
        fc = adsk.core.ObjectCollection.create()
        fc.add(face)
        tin = comp.features.thickenFeatures.createInput(
            fc, adsk.core.ValueInput.createByReal(sign * d_mm * MM), False,
            adsk.fusion.FeatureOperations.NewBodyFeatureOperation, False)
        made.append(comp.features.thickenFeatures.add(tin))
        new = [b for b in comp.bRepBodies if b.isSolid and b.entityToken not in before]
        if len(new) != 1:
            raise RuntimeError('thickening by %.2f mm gave %d solids' % (d_mm, len(new)))
        new[0].name = name
        return new[0]

    # where the fitted face dips into the bone, the clearance starts from there:
    # both thickenings go up by the dip, so the gap is at least c and the base
    # stays t thick
    dip = float(r['report'].get('sunk_max_mm', 0.0)) if c > 0 else 0.0
    rep['clearance_added_mm'] = round(dip, 3)
    guide = thicken(c + dip + t, BASE)
    if c > 0:
        inner = thicken(c + dip, 'Guide base clearance' + tag)
        made.append(combine(comp, guide, [inner],
                            adsk.fusion.FeatureOperations.CutFeatureOperation, False))
        guide = named(comp, BASE)
    sb.isLightBulbOn = False
    for sk in [m for m in made if adsk.fusion.Sketch.cast(m)]:
        sk.isVisible = False
    rep['thicken_s'] = round(time.time() - t1, 1)

    t1 = time.time()
    if p.get('_export_rim'):
        # dev: the top rim edge as built, with both faces' normals along it
        f = _face_at(comp, guide, _add(q, nb, (c + dip + t) * MM))
        dump = {'faces': guide.faces.count, 'top_edges': f.edges.count, 'edges': []}
        for e in f.edges:
            ok, t0_, t1_ = e.evaluator.getParameterExtents()
            ok, pts = e.evaluator.getStrokes(t0_, t1_, 0.0005)
            other = [x for x in e.faces if x.tempId != f.tempId]
            rows = []
            for pt in pts:
                row = [pt.x * CM, pt.y * CM, pt.z * CM]
                for face_ in [f] + other[:1]:
                    ok1, uv_ = face_.evaluator.getParameterAtPoint(pt)
                    ok2, nn = face_.evaluator.getNormalAtParameter(uv_)
                    row += [nn.x, nn.y, nn.z]
                rows.append(row)
            dump['edges'].append({'n': len(rows), 'rows': rows,
                                  'other_face_geom': type(other[0].geometry).__name__ if other else None})
        dump['face_types'] = [type(x.geometry).__name__ + ':%d' % x.edges.count for x in guide.faces]
        with open(p['_export_rim'], 'w') as fh:
            json.dump(dump, fh)
    if p.get('_rim_probe'):
        # dev: which rim radii build on this edge, each tried and undone
        probe = {}
        for rr_mm in p['_rim_probe']:
            f = _face_at(comp, guide, _add(q, nb, (c + dip + t) * MM))
            coll = adsk.core.ObjectCollection.create()
            for e in f.edges:
                coll.add(e)
            try:
                fin = comp.features.filletFeatures.createInput()
                fin.edgeSetInputs.addConstantRadiusEdgeSet(
                    coll, adsk.core.ValueInput.createByReal(rr_mm * MM), True)
                ft = comp.features.filletFeatures.add(fin)
                probe[str(rr_mm)] = True
                ft.deleteMe()
            except Exception:
                probe[str(rr_mm)] = False
            guide = named(comp, BASE)
        rep['rim_probe'] = probe
    keys = (('top', 'round', 'rim', c + dip + t), ('bottom', 'round_contact', 'contact', c + dip))
    for key, pkey, label, h in keys:
        rr = float(p.get(pkey, 0.0)) * MM
        if rr <= 0:
            continue
        def edges_of():
            f = _face_at(comp, named(comp, BASE), _add(q, nb, h * MM))
            return [e for e in f.edges] if f else []
        n_edges = len(edges_of())
        ok, built, msg, tried = fillet_down(comp, edges_of, rr * CM, True, 0.2,
                                            search or key == 'bottom')
        rep[label] = {'ok': ok, 'msg': msg, 'radius_mm': rr * CM, 'built_mm': built,
                      'failed_mm': tried, 'edges': n_edges}
        guide = named(comp, BASE)
    rep['rim_s'] = round(time.time() - t1, 1)
    return guide


def body_gap(body, bone):
    """Measured, not assumed: the body's minimum distance to the bone (mm) and the
    volume they share (mm3)."""
    app = adsk.core.Application.get()
    tbm = adsk.fusion.TemporaryBRepManager.get()
    d = app.measureManager.measureMinimumDistance(body, bone).value * CM
    inter = tbm.copy(body)
    try:
        tbm.booleanOperation(inter, tbm.copy(bone),
                             adsk.fusion.BooleanTypes.IntersectionBooleanType)
        ov = inter.volume * 1000.0
    except Exception:
        ov = 0.0
    return {'min_mm': round(d, 4), 'overlap_mm3': round(ov, 4)}


def _build_base(app, srv_ask, bone_paths, curves, p, rep, log, made, bone=None):
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.activeComponent
    op = outline_patch()
    raw = op.chain_outline(curves)
    pts_mm = [(q[0] * CM, q[1] * CM, q[2] * CM) for q in raw]
    t0 = time.time()
    r = srv_ask({'op': 'slab', 'coords': bone_paths['coords'], 'tris': bone_paths['tris'],
                 'pts': pts_mm, 'clearance': float(p['clearance']),
                 'thickness': float(p['thickness']),
                 'rail_spacing': float(p.get('rail_spacing', 1.0)),
                 'rim': float(p.get('round', 0.0))})
    rep['engine_s'] = round(time.time() - t0, 1)
    rep['fit'] = r['report']
    log('engine slab %.1f s' % rep['engine_s'])

    t1 = time.time()
    tmp, info = slab_body(r['topo'])
    (slab,), sf = add_bodies(comp, [tmp], ['Guide base slab'])
    made.append(sf)
    rep['assemble_s'] = round(time.time() - t1, 1)
    if not slab.isSolid:
        raise RuntimeError('the slab is not a closed solid (%s)' % '; '.join(info))
    rep['slab_info'] = info
    if p.get('_stop_after_slab'):
        slab.name = 'Guide base'
        return slab

    t1 = time.time()
    cur = r['curtain']
    tbm = adsk.fusion.TemporaryBRepManager.get()
    band = band_face(cur['lo'], cur['hi'])
    (band_body,), bf = add_bodies(comp, [band], ['Guide base curtain'])
    made.append(bf)
    before = set(b.entityToken for b in comp.bRepBodies)
    try:
        sin = comp.features.splitBodyFeatures.createInput(slab, band_body, False)
        made.append(comp.features.splitBodyFeatures.add(sin))
        rep['split'] = 'band'
    except Exception:
        # Fusion has reported no intersection for a closed band that plainly
        # crosses the slab; extended, it splits
        sin = comp.features.splitBodyFeatures.createInput(slab, band_body, True)
        made.append(comp.features.splitBodyFeatures.add(sin))
        rep['split'] = 'band, extended'
    pieces = [b for b in comp.bRepBodies if b.isSolid and
              (b.entityToken not in before or b.entityToken == slab.entityToken)]
    keep_in = _cm(r['keep']['inside'])
    guide = next((b for b in pieces if inside(b, keep_in)), None)
    if guide is None:
        raise RuntimeError('the curtain did not cut the slab round the outline (%d pieces)'
                           % len(pieces))
    made.extend(remove(comp, [b for b in pieces if b.entityToken != guide.entityToken]
                       + [band_body]))
    for sk in [m for m in made if adsk.fusion.Sketch.cast(m)]:
        sk.isVisible = False
    rep['cut_s'] = round(time.time() - t1, 1)

    t1 = time.time()

    def face_at(key):
        q = _cm(r['keep'][key])
        found = comp.findBRepUsingPoint(_pt(q), adsk.fusion.BRepEntityTypes.BRepFaceEntityType,
                                        0.06, False)
        for i in range(found.count):
            f = adsk.fusion.BRepFace.cast(found.item(i))
            if f and f.body.entityToken == guide.entityToken:
                return f
        return None
    # the rim first, then the round at the bone
    for key, pkey, label in (('top', 'round', 'rim'), ('bottom', 'round_contact', 'contact')):
        rr = float(p.get(pkey, 0.0)) * MM
        if rr <= 0:
            continue
        f = face_at(key)
        if f is None:
            rep[label] = {'ok': False, 'msg': '%s face not found' % key}
            continue
        edges = [e for e in f.edges]
        ok, msg = fillet(comp, edges, rr, True)
        part = None
        if not ok:
            # edge by edge: whatever builds, builds
            done = 0
            for i in range(len(edges)):
                f = face_at(key)
                es = [e for e in f.edges] if f else []
                if i < len(es) and fillet(comp, [es[i]], rr, False)[0]:
                    done += 1
            part = '%d of %d edges' % (done, len(edges))
            ok = done == len(edges)
        rep[label] = {'ok': ok, 'msg': '' if ok else msg, 'radius_mm': rr * CM,
                      'edges': len(edges), 'by_edge': part}
    rep['rim_s'] = round(time.time() - t1, 1)
    guide.name = 'Guide base'
    return guide


# --------------------------------------------------------------------------- #
# guide tubes, from scratch
# --------------------------------------------------------------------------- #

def principal_direction(pts):
    n = len(pts)
    c = [sum(p[k] for p in pts) / n for k in range(3)]
    C = [[sum((p[i] - c[i]) * (p[j] - c[j]) for p in pts) / n for j in range(3)] for i in range(3)]
    v = (1.0, 0.7, 0.3)
    for _ in range(100):
        v = _unit(tuple(sum(C[i][j] * v[j] for j in range(3)) for i in range(3)))
    return v


def tube_frame(u, ref, spin_deg):
    e1 = _sub(ref, tuple(_dot(ref, u) * x for x in u))
    if math.sqrt(_dot(e1, e1)) < 1e-6:
        e1 = _cross(u, (1.0, 0.0, 0.0) if abs(u[0]) < 0.9 else (0.0, 1.0, 0.0))
    e1 = _unit(e1)
    e2 = _cross(u, e1)
    s = math.radians(spin_deg)
    return (tuple(math.cos(s) * a + math.sin(s) * b for a, b in zip(e1, e2)),
            tuple(-math.sin(s) * a + math.cos(s) * b for a, b in zip(e1, e2)))


def tube_tmp(kind, bottom, top, R, e1, e2):
    """Cylinder of radius R, or the pentagonal prism inscribed in it with a
    vertex along e1 (five slabs intersected), between bottom and top (cm)."""
    tbm = adsk.fusion.TemporaryBRepManager.get()
    if kind != 'pent':
        return tbm.createCylinderOrCone(_pt(bottom), R, _pt(top), R)
    u = _unit(_sub(top, bottom))
    L = math.sqrt(_dot(_sub(top, bottom), _sub(top, bottom)))
    mid = _add(bottom, _sub(top, bottom), 0.5)
    ap = R * math.cos(math.pi / 5)
    W = 4.0 * R
    body = None
    for k in range(5):
        th = math.pi / 5 + 2 * math.pi * k / 5
        nk = tuple(math.cos(th) * a + math.sin(th) * b for a, b in zip(e1, e2))
        box = adsk.core.OrientedBoundingBox3D.create(
            _pt(_add(mid, nk, ap - W / 2)), _vec(u), _vec(nk), L, W, W)
        slab = tbm.createBox(box)
        if body is None:
            body = slab
        else:
            tbm.booleanOperation(body, slab, adsk.fusion.BooleanTypes.IntersectionBooleanType)
    return body


def tube_edges(body, kind, top, u, R):
    out = []
    for e in body.edges:
        g = e.geometry
        if kind == 'pent':
            line = adsk.core.Line3D.cast(g)
            if line and abs(_dot(_sub(_v(line.startPoint), top), u)) < 1e-5 and \
                    abs(_dot(_sub(_v(line.endPoint), top), u)) < 1e-5:
                out.append(e)
        else:
            c = adsk.core.Circle3D.cast(g)
            if c and math.sqrt(_dot(_sub(_v(c.center), top), _sub(_v(c.center), top))) < 1e-5 \
                    and abs(c.radius - R) < 1e-5:
                out.append(e)
    return out


def _native_token(ent):
    n = getattr(ent, 'nativeObject', None)
    return (n or ent).entityToken


def axis_hits(comp, bone, org, a, srv_ask):
    """Every point (cm) where the axis line meets the bone. A solid bone: Fusion's
    own ray test, both ways (a ray reports each face once) - no engine, no load of
    the bone's triangles. A mesh bone: the engine's ray."""
    hits = []
    if adsk.fusion.BRepBody.cast(bone):
        token = _native_token(bone)
        for s in (1.0, -1.0):
            pts = adsk.core.ObjectCollection.create()
            found = comp.findBRepUsingRay(_pt(_add(org, a, -s * 100.0)), _vec(tuple(s * x for x in a)),
                                          adsk.fusion.BRepEntityTypes.BRepFaceEntityType,
                                          -1.0, False, pts)
            for i in range(found.count):
                f = adsk.fusion.BRepFace.cast(found.item(i))
                if f is not None and _native_token(f.body) == token:
                    hits.append(_v(pts.item(i)))
        return hits
    for s in (1.0, -1.0):
        o = _add(org, a, -s * 100.0)
        h = srv_ask({'op': 'ray', 'o': [x * CM for x in o], 'd': [s * x for x in a]})
        if h.get('hit'):
            hits.append(_cm(h['hit']))
    return hits


def bone_normal(comp, bone, q, srv_ask):
    """Outward unit normal of the bone at a point on it (cm)."""
    if adsk.fusion.BRepBody.cast(bone):
        found = comp.findBRepUsingPoint(_pt(q), adsk.fusion.BRepEntityTypes.BRepFaceEntityType,
                                        0.01, False)
        token = _native_token(bone)
        for i in range(found.count):
            f = adsk.fusion.BRepFace.cast(found.item(i))
            if f is not None and _native_token(f.body) == token:
                ok, uv = f.evaluator.getParameterAtPoint(_pt(q))
                ok, n = f.evaluator.getNormalAtParameter(uv)
                n = _unit(_v(n))
                if bone.pointContainment(_pt(_add(q, n, 0.01))) == \
                        adsk.fusion.PointContainment.PointInsidePointContainment:
                    n = (-n[0], -n[1], -n[2])
                return n
    return tuple(srv_ask({'op': 'attach', 'pts': [[x * CM for x in q]]})['nrm'][0])


def build_tubes(app, comp, srv_ask, bone, guide, tubes, outline_cm, p, rep):
    """Each tube stands where its axis enters the bone on the guide's side (the
    engine's ray against the bone's triangles), is rounded, and trimmed by the
    bone when the bone is a solid. With a base (guide) the tubes are joined to
    it; without one - no outline picked - each tube is its own body, and the
    side it stands on is the one nearest the picked cylinder itself."""
    if outline_cm:
        ref = principal_direction(outline_cm[::max(1, len(outline_cm) // 300)])
    else:
        nears = [t[3] for t in tubes]
        ref = principal_direction(nears) if len(nears) >= 2 else (0.0, 0.0, 1.0)
    built, rows = [], []
    for (org, axis, r), kind, name, near in tubes:
        a = _unit(axis)
        hits = axis_hits(comp, bone, org, a, srv_ask)
        if not hits:
            rows.append({'name': name, 'error': 'the axis does not meet the bone'})
            continue
        if outline_cm:
            hit = min(hits, key=lambda q: min(_dot(_sub(q, c), _sub(q, c))
                                              for c in outline_cm[::5]))
        else:
            hit = min(hits, key=lambda q: _dot(_sub(q, near), _sub(q, near)))
        n = bone_normal(comp, bone, hit, srv_ask)
        u = a if _dot(a, n) > 0 else (-a[0], -a[1], -a[2])
        pent = kind == 'pent'
        R = 0.5 * float(p['pent_d' if pent else 'guide_d']) * MM
        tilt = min(math.acos(min(1.0, abs(_dot(u, n)))), math.radians(75.0))
        below = R * math.tan(tilt) + 0.2
        top = _add(hit, u, float(p['pent_h' if pent else 'guide_h']) * MM)
        e1, e2 = tube_frame(u, ref, float(p.get('guide_spin', 0.0)))
        built.append({'name': name, 'kind': kind, 'top': top, 'u': u, 'R': R,
                      'tmp': tube_tmp(kind, _add(hit, u, -below), top, R, e1, e2),
                      'axis': (org, axis, r, kind), 'tilt_deg': round(math.degrees(tilt), 1)})
    if not built:
        rep['tubes'] = rows
        return
    pre = set(b.entityToken for b in comp.bRepBodies)
    # the picked bodies' names can repeat (and an earlier run's tubes may be
    # there): each tube is found again by a name of its own until the end
    tag = run_tag()
    for i, b in enumerate(built):
        b['body_name'] = 'Guide tube %s%s.%d' % (b['name'], tag, i)
    bodies, _ = add_bodies(comp, [b['tmp'] for b in built], [b['body_name'] for b in built])
    brep_bone = adsk.fusion.BRepBody.cast(bone) is not None
    for b, body in zip(built, bodies):
        row = {'name': b['name'], 'kind': b['kind'], 'tilt_deg': b['tilt_deg']}
        if b['kind'] == 'pent':
            # the five TOP edges, as on the hand-made guides; the long edges stay sharp
            rr = min(float(p['pent_top_round']) * MM, 0.45 * b['R'])
        else:
            rr = min(float(p['guide_top_round']) * MM, 0.999 * b['R'])
        ok, built_mm, msg, tried = fillet_down(
            comp, lambda b=b: tube_edges(named(comp, b['body_name']), b['kind'], b['top'],
                                         b['u'], b['R']), rr * CM, b['kind'] == 'pent')
        row['round'] = 'ok' if ok else msg
        if ok and tried:
            row['round_built_mm'] = built_mm
        if brep_bone:
            before = set(x.entityToken for x in comp.bRepBodies)
            combine(comp, body, [bone], adsk.fusion.FeatureOperations.CutFeatureOperation, True)
            gap = float(p.get('guide_clear', 0.0)) * MM
            if gap > 0:
                row['gap'] = pull_off_bone(comp, bone, b, before, gap)
        rows.append(row)
    tube_bodies = [b for b in comp.bRepBodies if b.isSolid and b.entityToken not in pre]
    for tb in tube_bodies:
        tb.name = tb.name.split(tag)[0]
    if tube_bodies and guide is not None:
        combine(comp, guide, tube_bodies, adsk.fusion.FeatureOperations.JoinFeatureOperation,
                False)
    rep['tubes'] = rows
    rep['tube_bodies'] = tube_bodies if guide is None else []
    if not brep_bone:
        rep['tubes_note'] = 'the bone is a mesh: tube bottoms were not trimmed by it'
    return [b['axis'] for b in built]


def pull_off_bone(comp, bone, b, before, gap):
    """The faces a tube took from the bone in the cut, moved back gap off it
    (Offset Faces, along their own normal - into the tube)."""
    bodies = [x for x in comp.bRepBodies if x.isSolid and x.name == b['body_name']]
    faces = []
    for tb in bodies:
        for f in tb.faces:
            if bone.pointContainment(f.pointOnFace) == \
                    adsk.fusion.PointContainment.PointOnPointContainment:
                faces.append(f)
    if not faces:
        return 'no face on the bone'
    try:
        oin = comp.features.offsetFacesFeatures.createInput(
            faces, adsk.core.ValueInput.createByReal(-gap))
        comp.features.offsetFacesFeatures.add(oin)
    except Exception as exc:
        return str(exc).splitlines()[-1].strip()[:120]
    return 'ok, %d faces, %s' % (len(faces), measured_gap(comp, bone, b['body_name']))


def measured_gap(comp, bone, name):
    """The tube's real distance to the bone after the offset, and any overlap."""
    app = adsk.core.Application.get()
    tbm = adsk.fusion.TemporaryBRepManager.get()
    out = []
    for tb in [x for x in comp.bRepBodies if x.isSolid and x.name == name]:
        d = app.measureManager.measureMinimumDistance(tb, bone).value * 10.0
        inter = tbm.copy(tb)
        try:
            tbm.booleanOperation(inter, tbm.copy(bone),
                                 adsk.fusion.BooleanTypes.IntersectionBooleanType)
            ov = inter.volume * 1000.0
        except Exception:
            ov = 0.0
        out.append('min gap %.3f mm, overlap %.3f mm3' % (d, ov))
    return '; '.join(out) or 'not measured'


# --------------------------------------------------------------------------- #
# the whole guide
# --------------------------------------------------------------------------- #

def build_guide(app, bone, curves, tubes, holes, p, srv_ask, bone_paths, log=None):
    """bone: a solid, Form or mesh body (only its triangles are sampled);
    curves: the closed outline; tubes [((origin, axis, r) mm, kind, name)];
    holes [(origin, axis, r)] mm; p in mm; srv_ask: the engine; bone_paths: the
    bone's triangle files. Returns (body, report)."""
    log = log or (lambda *_: None)
    t0 = time.time()
    to_cm = lambda c: (tuple(x * MM for x in c[0]), tuple(c[1]), c[2] * MM)
    tubes = [(to_cm(t[0]), t[1], t[2], _cm(t[3]) if len(t) > 3 and t[3] else
              to_cm(t[0])[0]) for t in tubes]
    holes = [to_cm(c) for c in holes]
    des = adsk.fusion.Design.cast(app.activeProduct)
    comp = des.activeComponent
    rep = {}
    guide = build_base(app, srv_ask, bone_paths, curves, p, rep, log, bone) if curves else None
    t1 = time.time()
    axes = []
    if tubes:
        outline = outline_patch().chain_outline(curves) if curves else []
        axes = build_tubes(app, comp, srv_ask, bone, guide, tubes, outline, p, rep) or []
    rep['tubes_s'] = round(time.time() - t1, 1)
    t1 = time.time()
    cut = list(axes) + list(holes)
    if cut:
        tbm = adsk.fusion.TemporaryBRepManager.get()
        tmps = []
        for c in cut:
            org, axis, r = c[:3]
            kind = c[3] if len(c) > 3 else None
            key = {'cyl': 'hole_d', 'pent': 'pent_hole_d'}.get(kind)
            rr = 0.5 * float(p[key]) * MM if key and float(p.get(key, 0)) > 0 else r
            a = _unit(axis)
            tmps.append(tbm.createCylinderOrCone(_pt(_add(org, a, -100.0)), rr,
                                                 _pt(_add(org, a, 100.0)), rr))
        hb, _ = add_bodies(comp, tmps, ['Guide hole'] * len(tmps))
        if guide is not None:
            combine(comp, guide, hb, adsk.fusion.FeatureOperations.CutFeatureOperation, False)
        else:
            # every tube on its own: each takes every hole, the hole bodies go after
            for tb in rep.get('tube_bodies', []):
                combine(comp, tb, hb, adsk.fusion.FeatureOperations.CutFeatureOperation, True)
            remove(comp, hb)
    rep['holes'] = len(cut)
    rep['holes_s'] = round(time.time() - t1, 1)
    if guide is None:
        bodies = rep.pop('tube_bodies', [])
        rep.update({'body': '%d guide tubes' % len(bodies),
                    'solid': all(b.isSolid for b in bodies),
                    'faces': sum(b.faces.count for b in bodies),
                    'volume_mm3': round(sum(b.volume for b in bodies) * 1000.0, 1),
                    'secs': round(time.time() - t0, 1), 'tubes_only': True})
        return (bodies[0] if bodies else None), rep
    rep.pop('tube_bodies', None)
    guide.name = 'Guide base'
    rep.update({'body': guide.name, 'solid': guide.isSolid, 'faces': guide.faces.count,
                'volume_mm3': round(guide.volume * 1000.0, 1),
                'secs': round(time.time() - t0, 1)})
    return guide, rep
