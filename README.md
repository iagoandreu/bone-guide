# Bone Guide

A Fusion 360 add-in for patient-specific guides and plates, in two steps
- a curve on the bone, then a base grown inside it - without lofts and
without fillets that fail.

Two commands, in *Solid > Create* and *Surface > Create*:

| Command | What it does |
|---|---|
| **Bone Curve** | A closed curve on the bone - made from the drill cylinders (their alpha shape cut with the bone) or from a sketch curve - edited like a sketch spline whose points slide **over** the bone. |
| **Guide Base** | The base inside a closed curve: thickness along the bone normal, wall along the bone normal, rim round, holes. |

## How to use

### Where the buttons are

Design workspace. Both buttons, **Bone Curve** and **Guide Base**, are placed in
two panels:

- **Solid** tab > **Create** panel
- **Surface** tab > **Create** panel

### External Python (both commands)

Fusion's own Python has no numpy, so the geometry runs in a separate Python 3
process. It needs **numpy**, **scipy**, **trimesh** and **pyvista**
(`pip install numpy scipy trimesh pyvista`). Guide Base with **Output** set to
**Solid body (faceted)** and a **Facet tolerance** above 0 (the defaults) also
needs **pymeshlab** (`pip install pymeshlab`). The add-in does not check for
pymeshlab in advance; without it that build fails at the reduction step.

The add-in uses the first interpreter that can `import numpy, scipy, trimesh, pyvista`,
tried in this order:

1. the path in the environment variable **`BONEGUIDE_PYTHON`** (set before Fusion
   starts);
2. the path written in **`python_path.txt`** beside the installed add-in
   (`%APPDATA%\Autodesk\Autodesk Fusion 360\API\AddIns\BoneGuide\python_path.txt`),
   one line, the full path to `python.exe`;
3. `%LOCALAPPDATA%\Programs\Python\Python310\python.exe`, then `Python311`,
   `Python312`, `Python313`, `Python39` in the same folder;
4. `py -3.10`.

If none works, the command stops with "No Python with numpy, scipy, trimesh and
pyvista was found". Once found, the interpreter is kept until Fusion restarts.

### Bone Curve

A closed curve that lies on the bone, which you edit by dragging its points over
the bone.

#### Before you start

- **The bone**: a solid body or a mesh body. A scan STL imported as a mesh body is
  fine.
- And one of:
  - **Shapes**: the bodies the curve should wrap - drill or wire cylinders, a slot
    block, a sleeve - or an alpha shape body that is already made. The curve is
    where their alpha shape meets the bone.
  - **A closed sketch curve** near the bone, to lay on it and edit.
  - Nothing: you then click points on the bone to draw a new curve.
- The external Python (see above).

#### Step by step

1. Click **Bone Curve**.
2. **Bone**: select the bone. The selection then moves on to **Shapes**.
3. Either select the shapes in **Shapes**, or click the **Or a curve** box and
   select the sketch curve(s), or leave both empty. If both are filled, the shapes
   are used.
4. With shapes: set **Margin**, **Closing radius** and **Closing**.
5. Tick **Edit curve**. The curve is made and drawn as a thin blue line with round
   fit points on the bone. The selection boxes are emptied while you edit (the
   selection highlight would hide the curve); unticking **Edit curve** puts the
   picks back. The first time for a bone, loading it into the engine takes a
   moment; from shapes the curve can take 10-30 s.
6. Edit with the left mouse button:
   - **Drag a fit point**: it slides over the bone under the cursor and the curve
     follows, every sample on the bone.
   - **Click the curve** between points: a fit point is added there, and you can
     drag it straight away.
   - **Shift+click a fit point**: it is removed (a closed curve keeps at least 3).
   - With no shapes and no curve picked: **click on the bone** to add a point. It
     goes where it lengthens the outline least.
7. **Smooth the curve** smooths the fit points along the curve over the
   **Smoothing** length and lays them back on the bone (their number stays).
   **Undo last change** steps back (up to 30 steps).
8. Click **OK** to write the curve. OK is enabled once the curve has 3 or more fit
   points.

| Input | Default | What it does |
|---|---|---|
| **Bone** | none | The bone the curve lies on. One solid or mesh body. |
| **Shapes** | empty | Drill cylinders and other bodies (solid or mesh), or an alpha shape. A body that is one plain cylinder is used as its exact cylinder; any other body by its triangles. |
| **Margin** | 2.4 mm | Grows every shape in **Shapes** by this much before the closing - the cannula and the wall round it on a drill cylinder. Use 0 for a picked alpha shape. |
| **Shapes, own margin** | empty | Bodies that take **Their margin** instead - a cutting guide next to the drills, say. A body picked in both boxes takes this margin. |
| **Their margin** | 1 mm | The margin of the bodies in **Shapes, own margin**. |
| **Closing radius** | 3.5 mm | The alpha radius: gaps between the shapes narrower than twice this are filled. 0 cuts the shapes with the bone as they are. |
| **Closing** | In planes across the drills (alpha shape) | **In planes across the drills (alpha shape)**: closes section by section in planes across the cylinders' mean axis; side tabs that wrap round the bone edge stay with the rest. **3D**: closes with a ball, the same in every direction. |
| **Join every shape** | on | A shape further from the others than the closing bridges ends up on its own island. On: each island is joined to the region the other shapes are in by a strip along the shortest path over the bone, as wide as the closing radius (at least 1 mm each side). Off: islands are left out and named in the text box. A drill that also cuts the main region is never bridged where it comes out of the far cortex. |
| **Or a curve** | empty | A closed sketch curve (one, or several that join) to lay on the bone and edit. A single fitted spline keeps its own fit points. |
| **Edit curve** | off | On: the curve is made and the mouse edits it. While it is on, the inputs above are locked. |
| **Influence distance** | 0 mm | 0: only the dragged point moves, as in a sketch spline. Above 0: fit points within this distance along the curve follow, less and less with distance (cosine fall-off). |
| **Fit point spacing** | 3 mm | Spacing of the fit points on a curve made from the shapes, or from sketch curves that are not one fitted spline. |
| **Smoothing** | 1 mm | Smoothing length used by **Smooth the curve**. |
| **Smooth the curve** | button | Smooths the fit points and lays them back on the bone. Only while editing. |
| **Undo last change** | button | Undoes the last drag, add, remove or smooth. Only while editing. |

The values are remembered for the next time (`last_settings.json` in the add-in
folder).

#### What you get

- If you picked **one fitted spline** and did not add or remove points: the fit
  points of that spline are moved in place, in its own sketch. Anything built on
  the spline follows.
- Otherwise: a new 3D sketch named **Bone curve** in the **root component** (on the
  XY plane, but the curve is fully 3D), holding one closed fitted spline through
  exactly the fit points you saw.
- The bone and the shapes are not changed.

#### Checking the result

- The text box shows the number of fit points and the curve length. From shapes it
  also shows the closing used, the length of the loop kept, how many other loops
  on the bone were left out (where drills come out of the far cortex) and the time
  taken. The loop kept is the outer edge of the region holding the most shapes.
- Every shape the curve does not simply hold is named: **joined** (it was on an
  island, with the length of the strip), **left out** (an island, with Join every
  shape off), or **does not reach the bone** (with how far short it stops even
  after its margin; more than 5 mm is reported as such). A drill that stops above
  the bone cannot be in the curve - lengthen it to the bone or raise its margin.
- The overlay is drawn 0.15 mm off the bone so it stays visible. The written
  spline passes through the fit points on the bone; between fit points a sketch
  spline can leave the bone slightly, more where the points are far apart on a
  strongly curved bone.
- Look at the curve from several sides before OK, especially where it crosses
  ridges or wraps an edge.
- Problems while making the curve (no bone picked, an open curve, "the shapes do
  not reach the bone", "the curve is not near the bone") appear in bold in the text
  box. Other errors go to `BoneGuide.log` in the add-in folder.

#### Where it fits

Use it after the bone is in Fusion (the segmented STL as a mesh body, or a solid)
and the drill or screw trajectories exist as cylinders. It gives the outline of the
guide or plate footprint on the bone, which Guide Base (or Outline Patch) then
turns into the bone-contacting base.

#### Known limits

- The mouse editing works only while **Edit curve** is ticked.
- Bodies are read in model space: a bone inside a moved component occurrence is not
  handled yet.
- From shapes, only the longest loop on the bone is kept.
- A picked spline is updated in place only when the number of fit points did not
  change; otherwise a new **Bone curve** sketch is made and the old spline stays.

### Guide Base

The guide base inside a closed outline on the bone, with optional guide tubes and
holes, built from Fusion bodies: one surface, thickened, rounded. No mesh.

#### Before you start

- **The bone as a solid or Form body.** A scan is a mesh: convert it first with
  Mesh to Part, **Organic Form (quad remesh)** (about 20 s for a femur). The tubes are
  trimmed against this body and the gaps are measured against it, so use the same
  bone body you will check the guide against.
- **A closed outline on the bone**: one or several sketch curves that join end to
  end, for example the result of Bone Curve.
- Optional: **cylinders** for holes and guide tubes - cylindrical faces, or solid
  bodies (every cylindrical face of a picked body is used).
- The external Python (see above).

#### Step by step

1. Click **Guide Base**.
2. **Bone**: select the bone body.
3. **Outline**: select the closed outline curve(s). Leave it empty to build only
   the guide tubes (each as its own body).
4. Optional: select cylinders in **Holes**, **Cylindrical guides** and
   **Pentagonal guides**.
5. If you use guides, open the **Guide tubes** group and check its values.
6. Set **Clearance**, **Thickness** and the two rounds.
7. Click **OK** and read the report at the end.

| Input | Default | What it does |
|---|---|---|
| **Bone** | none | The bone the guide sits on: one solid or Form body. |
| **Outline** | none | The closed curve on the bone: one or more sketch curves that join. |
| **Holes** | empty | Cylindrical faces or cylinder bodies. A hole is cut along each, at its own radius. |
| **Cylindrical guides** | empty | Drill or wire cylinders. A round guide tube stands on the bone along each, with the hole through it. |
| **Pentagonal guides** | empty | Drill or wire cylinders. A pentagonal guide tube along each. A cylinder picked in both boxes becomes a round tube. |
| **Clearance** | 0.2 mm | Gap between the bone and the base. |
| **Thickness** | 3.5 mm | Of the base, from the clearance outward (the top is clearance + thickness off the bone). |
| **Rim round** | 2 mm | Round on the top rim. |
| **Contact edge round** | 0.2 mm | Round on the rim at the bone. |

**Guide tubes** group (collapsed). The defaults were measured on a hand-made distal
femur cutting guide:

| Input | Default | What it does |
|---|---|---|
| **Round: outside diameter** | 7.0 mm | Round tube diameter. |
| **Round: hole diameter** | 3.5 mm | 0: the picked cylinder's own diameter. |
| **Round: height above the bone** | 12.1 mm | From where the axis enters the bone, along the axis, to the top. |
| **Round: top round** | 1.5 mm | Round on the top rim. |
| **Pentagon: outside diameter** | 6.0 mm | Circle the pentagon is inscribed in (flats 4.85 mm apart). |
| **Pentagon: hole diameter** | 2.0 mm | 0: the picked cylinder's own diameter. |
| **Pentagon: height above the bone** | 10.5 mm | As for round tubes. |
| **Pentagon: top edge round** | 1.0 mm | Round on the five top edges; the long edges stay sharp. |
| **Gap to the bone** | 0.2 mm | The tube's bottom stays this far off the bone. |
| **Pentagon rotation** | 0 deg | Turn about the axis. 0: a corner points along the outline's long direction. |

The values are remembered for the next time (`last_settings.json` in the add-in
folder).

#### What you get

In the active component:

- **Guide base**: one solid - base, tubes and holes - built by Fusion features
  (base feature, split face, delete face, two thickens, combine, fillets, tubes,
  combines), all in the timeline;
- **Guide base surface** (hidden): the single surface the base was thickened from;
- **Guide base outline (eased)** (hidden sketch): the outline the base was cut
  with (see below). Your outline sketch is not changed.

Without an outline: one **Guide tube <name>** body per guide.

#### Checking the result

The report gives the measured **gap between the base and the bone** and their
**overlap** (should be 0), and the same for every tube; whether each round was
built; the bends that were opened in the outline and by how much it moved; and the
time. Before you trust the part, check:

- the base gap is at least the clearance and the overlap is 0;
- every tube shows its round and the gap you asked for;
- the narrowest hole is the drill or wire size you need;
- the base in the viewport, especially side tabs that wrap the bone.

Errors appear in a message box. Details go to `BoneGuide.log` in the add-in folder.

#### Where it fits

Use it after Bone Curve (or any closed outline on the bone), once the drill or
screw trajectories exist as cylinders. It produces the guide body ready for a final
check and STL export.

#### Known limits

- The bone must be a solid or Form body; convert a mesh with Mesh to Part first.
- Bodies are read in model space: a bone inside a moved component occurrence is not
  handled yet.
- Sleeves (cannulas) as separate parts are not generated; tubes and holes are.
- Convex bends of the outline tighter than 1.25 x the rim round are opened (the
  outline moves there, by about 0.3 mm on the test case); the report says where.

## Bone Curve in detail

Pick the **bone** (a solid or mesh body), then one of:

- **Shapes**: drill cylinders and any other bodies (a slot block, a sleeve), or
  an alpha shape already made. The curve is where their alpha shape meets the
  bone, computed inside the tool:
  - **Margin** grows every shape (the cannula and the wall round it on a drill
    cylinder; 0 for a picked alpha shape);
  - **Closing radius** is the alpha radius - gaps narrower than twice it are
    filled (0: the shapes are cut with the bone as they are);
  - **Closing** *In planes across the drills* (default) closes section by section
    in planes across the cylinders' mean axis, as an alpha shape lofted through
    plane sections is built; *3D* closes with a ball.
  A cylinder body is used as its exact cylinder; other bodies by their triangles.
  The loop round the region holding the most shapes is kept; shapes on islands
  are joined to it (**Join every shape**) or named as left out, and shapes that do
  not reach the bone are named with the gap.
- **Or a curve**: a closed sketch curve, laid on the bone. A single fitted spline
  keeps its own fit points.
- Neither: each click on the bone adds a point.

Tick **Edit curve**. The curve is a closed spline through **fit points that lie on
the bone**, drawn like a sketch spline - a thin blue line and round points:

- **Drag a point**: it slides under the cursor over the bone and the spline
  follows, every sample laid on the bone. **Influence distance** 0 moves only that
  point (a sketch spline); above 0 the fit points within that distance along the
  curve follow with a cosine fall-off.
- **Click the curve** between points to add one there (and drag it straight
  away); **Shift+click** a point to remove it.
- **Smooth the curve** smooths the fit points along the curve and lays them back
  on the bone (their number stays); **Undo last change** goes back one step.

OK writes the curve: into the picked spline when the number of fit points did not
change (moved in place, so whatever uses the spline follows), otherwise into a new
3D sketch *Bone curve* through exactly these fit points.

While editing, the selection boxes are emptied (the selection highlight paints the
bone blue over the curve); unticking **Edit curve** puts the picks back.

Measured on a distal radius test case: the 16 r 3.6 drill margins of the guide
plus the slot block, closed at 3.5 mm, give one near-side loop matching the
case's own alpha shape (r 3.50 mm, 6 sections) cut with the bone, median 0.29 mm (3D) / 0.47 mm (in
planes) apart; the other loops are the drills' exits on the far cortex. Cutting
that alpha shape body with the bone directly takes 2 s; closing the margins
inside the tool 12-25 s.

## Guide Base in detail

The base is built the way it would be by hand in Fusion, from bodies:

1. **One surface.** The engine fits a single bicubic NURBS face to the bone over the
   outline and 4 mm round it (a least-conformal flat map of that patch, 3 mm control
   spacing, faired). On the test femur it lies within 0.15 mm of the bone.
2. **Cut by the outline.** Split Face with the outline (closest point), outer face
   deleted: one face, one closed boundary edge.
3. **Thickened.** That face thickened by clearance + thickness, minus the same face
   thickened by the clearance (Combine, cut). Where the fitted face dips into the
   bone, both thickenings go up by the dip, so the measured gap is never under the
   clearance and the base keeps its thickness. Result: three faces - top, contact,
   one side wall - and one closed edge at each rim.
4. **Rounded.** Rim round on the top edge, contact round on the bottom edge.

Two things decide whether the rim round builds, both measured on the test case:

- **A round cannot be larger than a convex bend it rolls round**, and Fusion fails
  unpredictably when the two are close (1.5 - 1.8 mm failed on a 1.77 mm bend; 2.0
  built in some runs and not in others). So convex bends tighter than 1.25 x the rim
  round are opened first: a morphological opening (erode, then dilate by the same
  disc) of the region inside the outline, on the surface's flat map. It changes only
  those bends: on the test outline the tightest went from 1.77 to 3.2 mm and the
  outline moved 0.30 mm at most (p95 0.12 mm).
- **A closed edge has one vertex, where the closed outline starts**, and a round
  fails if that vertex sits in a bend: with it in a 2.2 mm bend, every round from
  1.5 to 2.5 mm failed; started on the straightest stretch, every round from 0.5 to
  2.5 mm built. The outline is restarted there, and if the round still fails the base
  is rebuilt with the outline starting on the next straight stretch (up to three
  tries; the report says how many were needed).

Measured on an artificial femur (the bone a Mesh to Part quad-remesh Form, 433 faces):
three runs: gap to the bone 0.200 mm, overlap 0, rim 2 mm and contact 0.2 mm rounds built
at the first try every time, 12 tubes each 0.200 mm off the bone, 94 - 96 s in all
(fit 17 s, split 11 s, thicken 5 s, rounds 27 s, tubes 18 s).

## Install

```
powershell -ExecutionPolicy Bypass -File install.ps1
```

copies the add-in to Fusion's AddIns folder; then *Scripts and Add-Ins > Bone
Guide > Run* (it runs on startup afterwards). The geometry runs in a Python with
numpy, scipy, trimesh and pyvista (plus pymeshlab for the faceted solid output):
set `BONEGUIDE_PYTHON`, or put its path in `python_path.txt` next to the
installed add-in; otherwise it tries `%LOCALAPPDATA%\Programs\Python\Python310`
(to 313, and 39) and then `py -3.10` (see *External Python* above).

## Tests

- `py -3.10 tests/test_engine.py` - offline. On a flat bone: two drill margins
  closed at 5 mm meet the bone in one loop bridged by radius-5 arcs (waist
  sqrt(48) - 5, within 0.02 mm in 3D and 0.1 mm in planes); the spline keeps its
  fit points and stays on a circle through 12 of them to 0.002 mm. On a cylinder
  "bone" where every answer is known: contact and top radii, volume = outline area x radial shell,
  the wall radial all round a 240 deg wrap, tab thickness, round and hole
  volumes.
- Live, through the dev bridge (`tools/bg.py`, needs `devbridge.on`): a new test
  document (`new_test_doc`), an editing session, a drag without a mouse (the same
  code path the mouse takes), write-back, Guide Base, captures. The dev commands
  modify only documents they created.

## Limits

- The mouse itself can only be tried by hand; everything under it (drag, click
  to add, Shift+click to remove, overlay, write-back) runs through the dev
  bridge's `selftest`, in a scratch document it closes afterwards.
- Bodies are read in model space: a bone inside a moved occurrence is not yet
  handled.
- Output is a mesh. For a solid, Mesh to Part.
- Sleeves (cannulas) are not generated yet; holes are.
