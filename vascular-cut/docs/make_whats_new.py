"""Generate docs/Whats-New.pdf.

The illustrated User-Guide.pdf has no source in this repository, so it cannot be
regenerated faithfully and is left alone. This produces a short companion covering
what changed since it was written, which is what a collaborator picking the tool up
from the share drive actually needs.

    python docs/make_whats_new.py
"""
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import LETTER
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    KeepTogether, ListFlowable, ListItem, PageBreak, Paragraph,
    SimpleDocTemplate, Spacer, Table, TableStyle,
)

OUT = Path(__file__).resolve().parent / "Whats-New.pdf"
INK = colors.HexColor("#1b2733")
DIM = colors.HexColor("#5b6875")
LINE = colors.HexColor("#d3dae1")
ACCENT = colors.HexColor("#c85a17")

ss = getSampleStyleSheet()
H1 = ParagraphStyle("H1", parent=ss["Heading1"], fontName="Helvetica-Bold",
                    fontSize=19, leading=23, textColor=INK, spaceAfter=4)
SUB = ParagraphStyle("SUB", parent=ss["Normal"], fontSize=10, leading=14,
                     textColor=DIM, spaceAfter=16)
H2 = ParagraphStyle("H2", parent=ss["Heading2"], fontName="Helvetica-Bold",
                    fontSize=13, leading=16, textColor=INK,
                    spaceBefore=16, spaceAfter=5)
H3 = ParagraphStyle("H3", parent=ss["Heading3"], fontName="Helvetica-Bold",
                    fontSize=10.5, leading=13, textColor=ACCENT,
                    spaceBefore=10, spaceAfter=3)
BODY = ParagraphStyle("BODY", parent=ss["Normal"], fontSize=9.8, leading=13.6,
                      textColor=INK, alignment=TA_LEFT, spaceAfter=7)
NOTE = ParagraphStyle("NOTE", parent=BODY, fontSize=9, leading=12.6, textColor=DIM)
CODE = ParagraphStyle("CODE", parent=BODY, fontName="Courier", fontSize=8.6,
                      leading=12, backColor=colors.HexColor("#f2f5f7"),
                      borderPadding=6, spaceBefore=3, spaceAfter=9)


def bullets(items):
    return ListFlowable(
        [ListItem(Paragraph(t, BODY), leftIndent=14, value="circle") for t in items],
        bulletType="bullet", start="circle", leftIndent=12,
        bulletFontSize=5, spaceAfter=6)


def table(rows, widths):
    t = Table(rows, colWidths=widths, hAlign="LEFT")
    t.setStyle(TableStyle([
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.6),
        ("TEXTCOLOR", (0, 0), (-1, -1), INK),
        ("LINEBELOW", (0, 0), (-1, 0), 0.7, LINE),
        ("LINEBELOW", (0, 1), (-1, -2), 0.3, colors.HexColor("#eef1f4")),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ("LEFTPADDING", (0, 0), (0, -1), 0),
    ]))
    return t


story = []
story.append(Paragraph("Vascular Cut Tool", H1))
story.append(Paragraph(
    "What changed since the illustrated User Guide. Read that first for the basics; "
    "this covers what is new. August 2026.", SUB))

# ── the arch preset ──────────────────────────────────────────────────────────
story.append(Paragraph("The aortic arch preset now knows where the arch is", H2))
story.append(Paragraph(
    "The arch preset used to place its two bounds by proportion along the model, "
    "which meant the same preset landed in a different anatomical place on every "
    "patient. Both bounds are now anchored on the vessels the region is named for.", BODY))
story.append(Spacer(1, 4))
story.append(table([
    ["Bound", "Placed at", "Why"],
    ["Ascending (proximal)", "brachiocephalic ostium − 10 mm",
     "matches hand-placed cuts: median −11.0 mm, IQR 1.05 mm"],
    ["Distal", "left subclavian ostium + 20 mm", "the existing rule, unchanged"],
], [1.5 * inch, 2.0 * inch, 2.6 * inch]))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "The 10 mm is measured, not chosen. Against thirty arches cut by hand, a fixed "
    "distance from the brachiocephalic ostium described where a person actually put "
    "the knife far better than a distance scaled by vessel calibre &mdash; "
    "coefficient of variation 0.27 against 0.59. Both bounds now carry a ladder of "
    "fallbacks, so a vessel that will not take a clean plane at the first offset is "
    "retried further along rather than refused: <b>30 of 30</b> arches land inside "
    "the intended band on the ascending bound, and 29 of 30 on both at once.", BODY))

story.append(Paragraph("A region boundary may cut its own aorta", H3))
story.append(Paragraph(
    "A region bound used to be refused where its plane crossed the aorta it belonged "
    "to &mdash; which is exactly what a region boundary is supposed to do. It now "
    "owns the vessels along its own route, so the bound cuts the aorta and stops at "
    "the branches, instead of declining because it found the aorta in the way.", BODY))
story.append(Paragraph(
    "Clearance from a branch ostium is applied to the distal bound only. The "
    "ascending bound is placed <i>relative to</i> the brachiocephalic ostium, so "
    "pushing it clear of that same ostium moved it away from the thing defining it.", NOTE))

story.append(Paragraph("Left carotid and left subclavian are no longer swapped", H3))
story.append(Paragraph(
    "On some arches the two left-sided vessels were named the wrong way round, which "
    "put the distal bound 20 mm past the wrong vessel. Their order along the arch is "
    "now checked and corrected &mdash; on the left only, because an aberrant right "
    "subclavian genuinely does reverse the order and must not be &ldquo;fixed&rdquo;.", BODY))

story.append(Paragraph("The export records the cut, not the request", H3))
story.append(Paragraph(
    "An exported region wrote down the station that was <i>asked for</i>. Where the "
    "tool had moved a bound to find a clean plane &mdash; which the fallback ladders "
    "now do routinely &mdash; the file disagreed with the geometry beside it, by up "
    "to 32 mm. It now records where each plane actually landed, keeping the request "
    "alongside when the two differ.", BODY))
story.append(PageBreak())

# ── cutting ──────────────────────────────────────────────────────────────────
story.append(Paragraph("Cuts that used to tear", H2))
story.append(Paragraph(
    "The biggest change in this release is not a feature. Cuts near a junction used "
    "to come out torn rather than sliced &mdash; the opening ragged, and a bite taken "
    "out of the vessel next door &mdash; and it was rarely obvious that it had "
    "happened. The cause was how the tool decided which vessel a piece of surface "
    "belonged to.", BODY))
story.append(Paragraph(
    "It used to give each face to the <i>nearest centerline</i>. That is wrong "
    "wherever a thin branch runs alongside a fat one: a trunk's axis lies a whole "
    "trunk radius deep, so a patch of aortic wall can sit closer to a small branch's "
    "centerline than to its own. A face is now given to whichever vessel it lies on "
    "the <b>wall of</b> &mdash; its distance compared against that vessel's own "
    "measured radius at that station.", BODY))
story.append(Spacer(1, 4))
story.append(table([
    ["Over 3,384 cuts on 66 aortas", "Before", "After"],
    ["Openings clean on all three checks", "1,653", "1,751"],
    ["Surface taken from the wrong vessel", "55,787 mm²", "26,377 mm²"],
], [3.0 * inch, 1.25 * inch, 1.25 * inch]))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "Clean means all three at once: the rim lies on the plane that made it, the "
    "result has exactly one more opening than it started with, and no surface was "
    "taken from a vessel the cut was not cutting. The last of those is the one nobody "
    "could see before it was measured &mdash; a cut can leak into the parent, still "
    "leave one connected surface, still leave one boundary loop, and still have the "
    "right Euler number.", NOTE))

story.append(Paragraph("It now refuses rather than hands back a bad opening", H3))
story.append(Paragraph(
    "Where a vessel has no clean cross-section left at the depth asked for, the tool "
    "says so and declines, instead of returning an opening that looks plausible in "
    "the window and is 15 mm out of plane. It also refuses a mesh and centerline that "
    "are not of the same anatomy, which is the mistake that is otherwise silent and "
    "total.", BODY))

story.append(Paragraph("The pull-back button", H3))
story.append(Paragraph(
    "&ldquo;Pull it back&rdquo; is the right advice and a useless amount of it. When "
    "an opening will not cut, the row now carries a button that finds <b>the deepest "
    "offset that actually works</b> and moves the slider there. The candidates are "
    "tested with the real rule &mdash; the rim must land on its plane, and the cut "
    "must take no more than a square millimetre of any other vessel &mdash; so the "
    "answer is one the tool can honour, not an estimate.", BODY))

story.append(Paragraph("Faster, and no longer freezing", H3))
story.append(Paragraph(
    "Each vessel's radius is measured once per model instead of on every cut, and "
    "read off the surface that is already there rather than by casting rays at it. "
    "Re-cuts resume from the last unchanged step rather than starting over, and the "
    "opening you moved is cut last so any slider can resume. On a five-opening model "
    "a slider commit takes about 0.6 s and dragging previews at roughly 36 ms a "
    "frame.", BODY))

# ── measurement ──────────────────────────────────────────────────────────────
story.append(Paragraph("Measuring a geometry", H2))
story.append(Paragraph(
    "The tool now measures as well as cuts. Both measurements run over the "
    "isolated region when there is one, and over the whole vessel otherwise.", BODY))

story.append(Paragraph("Diameter along the vessel", H3))
story.append(Paragraph(
    "The centerline is sampled every 2 mm and the mesh sectioned at each station. "
    "What is reported is the <b>equivalent-area diameter</b>, 2&#215;&#8730;(A/&#960;): "
    "the diameter of the circle with the same cross-sectional area. A vessel section "
    "is rarely a circle, so there is no single honest &ldquo;diameter&rdquo;. EAD is "
    "the standard vascular measure because area is what governs flow, and it does not "
    "depend on which way a caliper happens to be laid across an oval.", BODY))
story.append(Paragraph(
    "You get the maximum and where along the vessel it occurs, the minimum, the mean, "
    "both end diameters, the length and the number of stations.", BODY))

story.append(Paragraph("Neck angulation", H3))
story.append(Paragraph(
    "The circle that best fits the centerline over a region is found, and the angle "
    "reported is the arc it subtends: <b>&#952; = L / R</b>, centerline length divided "
    "by the fitted radius. For a seal zone that is neck angulation. Unlike a "
    "hand-drawn pair of tangents the answer does not depend on where the tangents are "
    "drawn, which is where inter-observer disagreement comes from.", BODY))
story.append(Paragraph(
    "Two diagnostics come with it and are worth reading before trusting the angle. "
    "<b>Fit RMS</b> is how far the centerline sits from that circle on average, and "
    "<b>planarity</b> is how far out of plane it wanders. &#952; = L/R only means "
    "something where the stretch really is close to a circular arc.", NOTE))
story.append(Spacer(1, 5))
story.append(table([
    ["Region", "Length", "Radius", "Angle", "Fit RMS"],
    ["Aortic neck", "58.0 mm", "123.0 mm", "27.0°", "0.12 mm"],
    ["Whole infrarenal aorta", "116.8 mm", "127.2 mm", "52.6°", "1.42 mm"],
], [1.9 * inch, 0.95 * inch, 0.95 * inch, 0.8 * inch, 0.9 * inch]))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "A worked example on one abdominal case. The neck fits a circle to a tenth of a "
    "millimetre, so 27&#176; is a number worth quoting. The full infrarenal stretch "
    "fits ten times worse, which is the honest signal that 117 mm of aorta is not one "
    "arc. Note also that a jagged centerline inflates arc length, and the angle "
    "inflates with it; leave the smoothing option on where the readout flags a "
    "staircase.", NOTE))

# ── regions ──────────────────────────────────────────────────────────────────
story.append(Paragraph("Regions that keep what is below them", H2))
story.append(Paragraph(
    "A region used to be two cut planes facing each other. Its distal plane stops "
    "just short of the end of the centerline, and the iliac ostia sit exactly at the "
    "bifurcation, so they fell outside every region: <b>Infrarenal aorta</b> reported "
    "no branches inside and cut the iliacs away rather than keeping them.", BODY))
story.append(Paragraph(
    "There is now a one-sided region: a single plane, with everything below it kept. "
    "<b>Infrarenal aorta with iliacs</b> uses it, so you can remove everything above "
    "the renals and still have the iliacs as ordinary openings, each with its own "
    "slider. The end slider is disabled while such a region is active, since there is "
    "no end to move.", BODY))
story.append(Paragraph(
    "<b>Fixed:</b> moving the start plane of such a region used to close it. The "
    "disabled end slider still held a number &mdash; zero in a fresh session, or the "
    "end of whichever two-plane region you used before &mdash; and the first nudge of "
    "the start slider read that number back. A distal plane appeared where there had "
    "been none, the iliacs fell outside the region, and it could not be undone, "
    "because the slider now driving that plane was still greyed out. The start plane "
    "moves freely now and the open end stays open. Reported by a collaborator who cut "
    "a 49-case cohort; their workaround was to use <b>Infrarenal aorta</b> instead, "
    "which is no longer necessary.", NOTE))

story.append(Paragraph("Trimming inside a region", H2))
story.append(Paragraph(
    "Trims taken while a region was isolated were recorded and then discarded, so "
    "nothing changed on screen. They apply now, and as many as the region contains. "
    "The menu is also scoped to what is in view: it used to offer trims for branches "
    "the region does not contain.", BODY))
story.append(Paragraph(
    "So the abdominal workflow is: isolate <b>Infrarenal aorta with iliacs</b>, then "
    "trim each common iliac before it splits. Both iliac pairs collapse to one "
    "opening each.", BODY))

story.append(Paragraph("Which side is left, and when the tool does not know", H2))
story.append(Paragraph(
    "Left and right are decided for a whole case by a single reading, taken from a "
    "witness in the anatomy: on an abdominal aorta the splenic and hepatic branches "
    "of the celiac trunk (the spleen is on the left, whatever the file says), on an "
    "arch the brachiocephalic trunk. Every side name on the case follows from that "
    "one reading, so if it is wrong they are all wrong together &mdash; renals and "
    "iliacs alike.", BODY))
story.append(Paragraph(
    "When no witness is found the tool falls back to the coordinate frame "
    "(+x = patient left) and carries on naming sides. That fallback was recorded only "
    "in the exported metadata, so on screen a guessed case looked exactly like a case "
    "the tool had actually read. It no longer does: the caption now says <b>LEFT AND "
    "RIGHT ARE NOT CONFIRMED</b> and explains that a mirrored export would swap every "
    "side name without the geometry being able to tell. In the last 59-case abdominal "
    "run, 11 cases were in this position.", BODY))
story.append(Paragraph(
    "Both the desktop window and the browser tool now also draw an <b>L / P / S "
    "orientation triad</b> in the corner of the view, so a side can be checked "
    "against something. Worth knowing when you do: an anterior view puts the "
    "patient's left on <i>your</i> right.", NOTE))

story.append(Paragraph("Aortic neck", H2))
story.append(Paragraph(
    "The neck runs from the lowest renal artery to where the aorta first exceeds that "
    "calibre. That threshold is now <b>10%</b>, the reporting convention, rather than "
    "the 20% used before. 20% reports a longer neck than a reader following the "
    "standard would measure, which matters against the 15 mm minimum for standard "
    "EVAR.", BODY))

story.append(PageBreak())

# ── the web version ──────────────────────────────────────────────────────────
story.append(Paragraph("The tool in a browser", H2))
story.append(Paragraph(
    "There is now a web version, for collaborators who would rather not install "
    "anything. It runs the same VTK pipeline as the desktop window, rendered on the "
    "server and streamed to the browser, so the cut planes stay live and draggable.", BODY))
story.append(Paragraph("python webapp/cut_app.py&nbsp;&nbsp;&nbsp;&#8594;&nbsp;&nbsp;"
                       "http://127.0.0.1:8080", CODE))
story.append(Paragraph(
    "It binds to 127.0.0.1 deliberately. The app has <b>no authentication</b>, and "
    "binding to every interface would put it on the institutional network and raise a "
    "firewall prompt. If you ever see that prompt when starting it, do not allow it: "
    "check the configuration first.", NOTE))
story.append(Paragraph(
    "One person at a time. The built-in server is single-session, so run your own "
    "copy rather than sharing one.", NOTE))
story.append(Paragraph(
    "Not yet in the web version: the trim menu, the pinned measurement markers, and "
    "file upload, which still takes typed paths. The desktop tool remains the "
    "complete one.", NOTE))

# ── installing ───────────────────────────────────────────────────────────────
story.append(Paragraph("Installing from the share drive", H2))
story.append(Paragraph(
    "Copy the zip anywhere and unzip it. Everything sits under a "
    "<font face='Courier'>vascular-cut</font> folder. Then either double-click "
    "<font face='Courier'>run.bat</font>, or:", BODY))
story.append(Paragraph(
    "python -m venv .venv<br/>"
    ".venv\\Scripts\\pip install -r requirements.txt<br/>"
    ".venv\\Scripts\\python branch_gui.py", CODE))
story.append(Paragraph(
    "Python 3.10, 3.11 or 3.12. Not 3.13: several dependencies have no wheels there "
    "and the install fails. Two sample aortas ship with it, one branched and one not, "
    "so there is something to open straight away.", NOTE))
story.append(Paragraph(
    "<font face='Courier'>run.bat</font> now finds conda when it is installed but not "
    "on PATH &mdash; Miniforge, Mambaforge, Miniconda and Anaconda, per-user or "
    "system-wide, falling back to whatever <font face='Courier'>CONDA_EXE</font> "
    "points at. Previously it reported that conda was missing on a machine that had "
    "it, which is a confusing way to be told to open a different prompt.", NOTE))
story.append(Paragraph(
    "The curvature metrics (&#948;K) additionally need the lab's CAMP source, which "
    "is vendored under <font face='Courier'>camp/</font>. See "
    "<font face='Courier'>camp/README.md</font> for its provenance and how to refresh "
    "it when the lab's version moves on.", NOTE))

# ── naming, and the preset that was hiding in the batch tool ─────────────────
story.append(Paragraph("Saying which vessel is which", H2))
story.append(Paragraph(
    "An opening's name box is now a list as well as a text box. Typing in it still "
    "relabels the opening; picking an artery from the list <b>declares</b> it, and "
    "the region presets are built from that. On an abdominal aorta a declared renal "
    "is what brings back Aortic neck, Infrarenal aorta and Infrarenal aorta with "
    "iliacs.", BODY))
story.append(Paragraph(
    "This is not a workaround for a weak rule. Across 47 abdominal scans the tool "
    "names no renal on 18 of them, and on 14 of those because <i>several</i> lateral "
    "branch pairs fit the renal pattern equally well &mdash; an accessory or polar "
    "renal is about 30% of people. Naming either would be a guess presented as a "
    "fact, so it correctly names neither. You can see which pair it is; the geometry "
    "cannot say. A declared name also displaces the automatic name of the same "
    "artery, so a boundary placed at “the lowest renal” cannot land on the "
    "vessel you were correcting.", NOTE))

story.append(Paragraph("The nub-in preset is now a button", H3))
story.append(Paragraph(
    "Cut every branch back to a short stub past the aortic wall and square off the "
    "aorta's own ends, in one action. The stub is 0.75 of each vessel's own radius "
    "&mdash; measured over 82 vessels in 30 cases, not chosen: a branch leaves its "
    "aorta as a funnel, a plane through a funnel does not come back flat, and so "
    "each vessel has a shortest nub it can take. Until now this existed only in the "
    "batch cohort tool.", BODY))
story.append(Paragraph(
    "A nub that will not cut cleanly is moved <i>outward</i> until it does, never "
    "inward, so a rescued nub is longer than asked for and never ends up back inside "
    "the aorta. Every rescue is reported, as is any aortic end too oblique for a "
    "horizontal plane, because neither is visible in the finished geometry.", NOTE))

story.append(Paragraph("Paths can be typed", H3))
story.append(Paragraph(
    "The mesh and centerline boxes accept a pasted path, like the output folder "
    "always has. Quotes and trailing spaces are trimmed, so a path copied from "
    "Explorer goes straight in, and a path that is not there says so immediately "
    "rather than at Load.", BODY))

story.append(Paragraph("A note on what is measured", H2))
story.append(Paragraph(
    "Both measurements depend on where the geometry was cut. Surface area, and "
    "therefore any curvature metric computed from it, changes with the cut extent; "
    "angulation and diameter change with which region is isolated. Record the region "
    "you used. The exported sidecar JSON does this for you: it carries every cut "
    "plane, offset, tilt, plane size and resulting diameter, plus the region's "
    "stations, so a measurement can be reproduced rather than approximated.", BODY))

# ── delta K: the reproducibility finding ─────────────────────────────────────
#
# The most consequential thing learned this cycle, and the easiest for a reader to
# get wrong. It goes in the shipped PDF rather than only in the repository, because
# the people who read a delta K off the tool are not the people who read CURVATURE.md.
story.append(Paragraph("How far to trust a &#948;K", H2))
story.append(Paragraph(
    "<b>These figures replace the ones in earlier copies of this document.</b> Two "
    "defects were found in how CAMP builds its patches and both are now fixed, so "
    "&#948;K is markedly steadier than the old table said. If you are holding a "
    "Whats-New printed before September 2026, discard its repeatability figures: they "
    "were measured with those defects present and are far more pessimistic than the "
    "tool now is.", NOTE))
story.append(Paragraph(
    "&#948;K is not a property of the surface alone. It is a property of the surface "
    "<i>and</i> of a partition that CAMP rebuilds from scratch on every run. The two "
    "defects were patches arriving in <b>two disconnected pieces</b> (a branch wall "
    "and the aortic wall beside it are close through the air and half a circumference "
    "apart along the surface) and patches <b>wrapped around</b> a junction, which no "
    "height field can describe and which R&#178; and &#951; both pass. Fixing them "
    "tightened the healthy-aorta band from 0.14&#8211;3.47 to 0.66&#8211;1.37.", BODY))
story.append(Spacer(1, 4))
story.append(table([
    ["30 aortas, 100k faces, no trim at all", "Spread of the normals", "CV"],
    ["Before either fix", "0.14 &#8211; 3.47", "1.04"],
    ["Connectivity fix only", "0.21 &#8211; 2.53", "0.78"],
    ["Both fixes (current)", "0.66 &#8211; 1.37", "0.23"],
], [3.0 * inch, 1.4 * inch, 1.1 * inch]))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "The trim is no longer doing structural work. It used to be needed to suppress "
    "those defects; with them gone, <b>trim 0 now beats the old reference trim</b>, "
    "which managed CV 0.30. Note the trim takes <i>both</i> tails, so the fraction "
    "shown is per tail.", NOTE))

story.append(Paragraph("One partition is one draw", H3))
story.append(Paragraph(
    "This is the part to carry into any comparison. &#948;K is carried by a small "
    "number of patches &mdash; typically about 22 of some 850 &mdash; sitting where "
    "curvature changes fastest, which on a branched aorta means the branch ostia. "
    "Half a per cent of the faces hold about 85% of the run-to-run variance. Where "
    "k-means happens to put a boundary in those few places decides the number, so "
    "<b>re-running the same mesh with a different seed moves &#948;K by about "
    "16%</b>.", BODY))
story.append(Paragraph(
    "The draws are independent, verified against disjoint halves, so averaging "
    "delivers what the arithmetic promises. The panel now prints the band beside the "
    "value at every appearance.", BODY))
story.append(Spacer(1, 4))
story.append(table([
    ["Partitions averaged", "1", "4", "8", "16"],
    ["Error on the mean", "&#177;16%", "&#177;8%", "&#177;6%", "&#177;4%"],
], [1.9 * inch, 0.85 * inch, 0.85 * inch, 0.85 * inch, 0.85 * inch]))
story.append(Spacer(1, 4))
story.append(Paragraph(
    "What that supports: &#948;K separates populations and large individual "
    "differences. It does <i>not</i> support fine ordering &mdash; across 30 cases "
    "spanning 0.68 to 2.5, neighbours sit about 6% apart, so at 8 partitions only one "
    "of 29 adjacent-in-rank pairs is resolvable at two standard errors.", NOTE))
story.append(Paragraph(
    "<b>Comparable only at a fixed face count.</b> &#948;K drifts about 21% across a "
    "fourfold change in density, because the minimum patch size is a count of "
    "<i>vertices</i>: the smallest allowed patch covers about 16 mm&#178; at 100k "
    "faces and 4 mm&#178; at 400k. Those are measurements at different physical "
    "scales, not repeat draws, and averaging does not touch it. Always compare at the "
    "same target. 100,000 faces is the reference, chosen across seven densities from "
    "25k to 400k.", NOTE))
story.append(Paragraph(
    "Patches whose curvature fit never converges are counted, always, and cannot be "
    "excluded. A patch that will not take a quadratic is precisely where the surface "
    "is least like one, which is what &#948;K exists to measure; dropping those does "
    "not clean the number up, it removes most of it.", NOTE))

story.append(Paragraph("Folded surfaces are flagged before the run", H3))
story.append(Paragraph(
    "The tool checks any surface it is given for <b>folds</b> &mdash; neighbouring "
    "faces meeting at more than 90&#176;, where the surface doubles back through "
    "itself &mdash; and says so before spending half a minute measuring it. That is "
    "geometry the vessel does not have, and it arrives in the segmentation rather "
    "than being introduced by the cut. 14 of the 30 reference meshes carry some.", BODY))
story.append(Paragraph(
    "They are reported, not corrected, and the reason is now better understood than "
    "when this document first said so. Fold and sharp-crease counts <i>do</i> track "
    "how much a case moves between runs, at &#961; +0.56 to +0.61 &mdash; an earlier "
    "figure of &#8722;0.07 here was measured before the two partition fixes and no "
    "longer holds. But repairing them was tried at four thresholds (90&#176;, "
    "75&#176;, 60&#176;, 45&#176;) and <b>every one made the error bars wider</b>, "
    "because the repair smooths exactly where curvature is sharpest, which is exactly "
    "where &#948;K comes from. The defects are a marker of sharp anatomy, not its "
    "cause. Fix a folded surface upstream if the value matters.", NOTE))

doc = SimpleDocTemplate(
    str(OUT), pagesize=LETTER,
    leftMargin=0.9 * inch, rightMargin=0.9 * inch,
    topMargin=0.85 * inch, bottomMargin=0.8 * inch,
    title="Vascular Cut Tool - What's New", author="Pocivavsek Lab")
doc.build(story)
print(f"wrote {OUT}  ({OUT.stat().st_size / 1024:.0f} KB)")
