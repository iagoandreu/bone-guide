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

The guide base inside a closed curve on the bone, optionally with guide tubes and
holes, built as one mesh and optionally converted to a faceted solid.

#### Before you start

- **The bone**: a solid body or a mesh body.
- **A closed outline on the bone**: one or several sketch curves that join end to
  end, for example the result of Bone Curve.
- Optional: **cylinders** for holes and guide tubes - cylindrical faces, or solid
  bodies (every cylindrical face of a picked body is used).
- The external Python (see above), with pymeshlab for the default solid output.

#### Step by step

1. Click **Guide Base**.
2. **Bone**: select the bone.
3. **Outline**: select the closed outline curve(s).
4. Optional: select cylinders in **Holes**, **Cylindrical guides** and
   **Pentagonal guides**.
5. If you use guides, open the **Guide tubes** group and check its values.
6. Set **Clearance**, **Thickness**, the rounds, **Resolution** and **Output**.
7. Click **OK**. A progress dialog shows the engine's current step and the
   elapsed time; you can cancel it. A build takes from tens of seconds to a few
   minutes, depending on size and resolution.
8. Read the report in the message box at the end.

| Input | Default | What it does |
|---|---|---|
| **Bone** | none | The bone the guide sits on. One solid or mesh body. |
| **Outline** | none | The closed curve on the bone: one or more sketch curves that join. |
| **Holes** | empty | Cylindrical faces or cylinder bodies. A hole is cut along each, at its own radius. |
| **Cylindrical guides** | empty | Drill or wire cylinders (faces or bodies). A round guide tube stands on the bone along each, with the hole through it. |
| **Pentagonal guides** | empty | Drill or wire cylinders. A pentagonal guide tube along each, with the hole through it. A cylinder picked in both boxes becomes a round tube. |
| **Clearance** | 0.2 mm | Gap between the bone and the contact face. |
| **Thickness** | 3.5 mm | Measured along the bone normal from the contact face: the top is the bone offset by clearance + thickness. |
| **Rim round** | 1 mm | Radius of the round between the top and the wall. A field blend: it cannot fail and never touches the contact face. |
| **Contact edge round** | 0 mm | Optional round where the wall meets the contact face. 0 keeps the contact face whole. |
| **Outline smoothing** | 0.5 mm | Smoothing along the outline before it is laid on the bone. 0 = as drawn. |
| **Resolution** | Medium - 0.15 mm voxels | **Fine - 0.10 mm voxels**, **Medium - 0.15 mm voxels** or **Draft - 0.25 mm voxels**. Finer is truer and slower. |
| **Output** | Solid body (faceted) | **Mesh body**: the contoured mesh as it is. **Solid body (faceted)**: the mesh is first reduced to within the facet tolerance, then converted with Convert Mesh (faceted, one face per triangle). |
| **Facet tolerance** | 0.05 mm | For the solid output: how far the reduced mesh may stray from the contoured surface, measured both ways. 0 skips the reduction. |

**Guide tubes** group (collapsed):

| Input | Default | What it does |
|---|---|---|
| **Outside diameter** | 5.2 mm | Tube diameter. A pentagonal tube is inscribed in it (5.2 mm: flats 4.21 mm apart). |
| **Height above the bone** | 10 mm | From where the axis enters the bone, along the axis, to the top. |
| **Gap to the bone** | 0.5 mm | The tube's bottom stays this far off the bone where it stands out of the base. |
| **Hole diameter** | 0 mm | 0: the picked cylinder's own diameter. |
| **Cylinder top round** | 1 mm | Round on the top rim of a cylindrical tube. |
| **Pentagon edge round** | 1 mm | Round on the five long edges of a pentagonal tube; its top is left sharp. |
| **Pentagon rotation** | 0 deg | Turn of the pentagon about its axis. 0: a corner points along the guide's length (the outline's long direction). |

The values are remembered for the next time (`last_settings.json` in the add-in
folder).

#### What you get

All in the **root component**:

- a mesh body named `Guide base t<thickness> c<clearance> r<rim round>`, for
  example `Guide base t3.50 c0.20 r1.00`. In a parametric design it sits in a base
  feature of the same name;
- with **Solid body (faceted)**: a Convert Mesh feature and a solid body with the
  same name; the mesh body is hidden.

The base, the tubes and the holes are one body. The bone, the outline and the
cylinders are not changed. The parameters and results are stored as attributes on
the body (attribute group `BoneGuide`).

#### Checking the result

The message box at the end reports:

- triangles, volume (mm3) and footprint (mm2);
- the contact face's median distance from the bone, against the clearance asked,
  and the highest point of the top;
- the outline length, the number of holes and the time;
- with guides: their number by type, the narrowest hole diameter, the range of
  tube top heights above the bone, and the smallest tube-to-bone gap against the
  gap asked;
- loose specks left out (pieces not connected to the main body), with their total
  volume and largest size;
- with the solid output: triangles before and after the reduction, the measured
  deviation against the facet tolerance, and the number of faces.

Before you trust the part, check:

- the contact distance equals the clearance you asked for;
- the narrowest hole diameter is the drill or wire size you need;
- the tube-to-bone gap is not smaller than asked;
- the loose specks: a large one means part of the design came out disconnected
  and was dropped;
- the base in the viewport, especially side tabs that wrap the bone, and the
  rounds.

Errors appear in a message box (for example "The outline must be closed",
"nothing lies inside the curve", "a guide axis does not meet the bone"). Details go
to `BoneGuide.log` in the add-in folder.

#### Where it fits

Use it after Bone Curve (or any closed outline on the bone), once the drill or
screw trajectories exist as cylinders. It produces the guide body - base, tubes and
holes - ready for a final check and STL export for printing.

#### Known limits

- The output is a mesh, or a faceted solid (one face per triangle). For a smooth
  solid, convert it with another tool (for example Mesh to Part).
- The solid output needs pymeshlab in the external Python.
- Bodies are read in model space: a bone inside a moved component occurrence is not
  handled yet.
- Sleeves (cannulas) as separate parts are not generated; tubes and holes are.
- Resolution is the voxel size of the field: details much smaller than a voxel are
  lost.

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

Pick the bone, the closed outline (one or several sketch curves that join) and,
optionally, **Holes**: cylindrical faces or cylinder bodies, each cut through at
its own radius.

The inputs, including the guide tubes and the output type, are listed in
*How to use > Guide Base* above.

The base is **one field expression**, contoured once:

- `d(x)`: signed distance to the bone (dense samples with interpolated normals).
- `w(x)`: signed distance to the **wall** - the ruled surface of bone normals
  standing on the outline.
- base = `{clearance < d < clearance + thickness}` and `{w < 0}`, the top/wall
  edge rounded with a circular blend (Quilez's round intersection), holes cut by
  `max(f, r - distance to axis)`.

So the contact face is the bone offset by the clearance exactly; the top is the
true offset at clearance + thickness (a concave spot rounds off instead of
folding); the wall stands perpendicular to the bone all round, so **side tabs
that wrap the bone keep their thickness** - which a loft or a plan-view
extrusion does not; and the round is a property of the field, so it cannot fail
and never touches the contact face. The result is a watertight mesh body in its
own base feature, converted to a faceted solid when **Output** asks for it.

Measured on a distal radius test case (the radius as the bone, offline scripts):

- From the footprint outline of the case's own guide, the base regenerates that
  guide: volume 2,227 vs 2,235 mm3, generated-to-guide p95 0.73 mm,
  guide-to-generated p95 0.58 mm,
  contact face 0.200 mm. 45 s at 0.15 mm.
- In Fusion, on the intersection contour after a 4 mm drag: 342 k triangles,
  contact 0.200 mm, 18 holes, 44 s.

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
