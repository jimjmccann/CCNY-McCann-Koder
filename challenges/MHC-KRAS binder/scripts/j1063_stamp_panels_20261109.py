#!/usr/bin/env python3
"""j1063 stage 3/4 -- crop the raw renders to ONE COMMON BOX, stamp a TITLE BAND, re-measure.

    python3 -I j1063_stamp_panels_20261109.py <rawdir> <outdir>

⭐ METHOD AND CREDIT: `an earlier challenge's scripts/j575_stamp_panel_title_20261019.py` (the band ABOVE
   the grid, for a marker that a per-cell label can never carry), `TfR1 binder/scripts/
   j957_stamp_panels_20261102.py` and `competition5/<WORKDIR>/scripts/s94_stamp_arm_panels.py`
   (the fit/record/re-measure discipline and the wrapping legend). Local, PIL only.
⛔ Opens no structure and changes no geometry: this is presentation.

⭐⭐ WHY THE BAND EXISTS, AND WHY IT IS NOT DECORATION. `[M]` three times on 2026-10-19 a PyMOL
   per-cell labeller strips whatever is COMMON across the label set: `DES_s07_3_s05` drew as
   `s07_3_s05`. ⇒ A CONSTANT MARKER CAN NEVER SURVIVE INSIDE A PER-CELL LABEL, because within one
   panel it is by definition common [[figure-label-cannot-carry-a-constant-marker]]. The TIER WORD
   on these panels is exactly such a constant, and it is the single most important thing on the
   figure. ⇒ it goes in a band ABOVE the image area, where nothing can strip it.

⛔⛔ THE TIER WORDING IS EXACT AND IS NOT MINE TO SOFTEN. James McCann put all seven designs in the
   submission; that is a SCOPE decision and it does not change what any of them is.
     4 designs  VALIDATED                       -- clear the conjunction at 4 of 4 seeds, 20 draws
     3 designs  3 of 4 SEEDS, INSIDE NOISE, NOT PROMOTED   -- J1058 section 6 tiering
   ⛔ The word CANDIDATE is NOT used anywhere: in that pre-registered scheme it means a
   5-draw-only design and none of these are. ⛔ No panel may read validated, promoted or result
   unless its own manifest entry says VALIDATED.

⛔⛔ A STAMPED LINE THAT RUNS OFF THE RIGHT EDGE IS A **WRONG LABEL**, not a cosmetic problem:
   `[M]` the Challenge-4 first pass truncated the column name `rg_ratio` to `rg_`. So
   (a) every piece of text goes through `fit()`, which shrinks until the line fits the width
       actually available and RAISES below a minimum size rather than drawing something
       truncated; and
   (b) every drawn line is recorded with its text, font, size, x, baseline y and allotted width,
       and a SEPARATE PASS re-opens the FINISHED PNGs and re-measures every line against its own
       canvas. ⛔ That pass, not the drawing pass, is what reports 0 overflows.
   ⛔⛔ AND THE LEGEND FIT MUST INCLUDE THE CHIPS AND THE GAPS, NOT JUST THE TEXT: `[M]` a 6-key
   row passed a text-only fit and then ran 258 px off the right edge once the chips, gaps and
   separators were laid down. Here the overview carries NINE keys, so the row WRAPS over as many
   rows as it needs rather than overflowing.

⛔ THE CROP CAN ONLY REMOVE BORDER. The box is the union of the white-background content boxes of
   ALL EIGHT renders, so it cannot clip structure in any of them, and it is the SAME box for every
   panel so the panels stay pixel-comparable. Before/after sizes are printed.

⛔⛔ AND THE NUMBERS ARE NOT FREE TEXT. Every figure stamped here is read from the stage-2
   manifest, which carried it from stage 1, which asserted it cell-by-cell against the scored
   tables. Nothing is measured here and nothing is claimed that is not in that chain.
"""
import json
import os
import sys

from PIL import Image, ImageChops, ImageDraw, ImageFont

FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONT_B = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"

TARGET_GREY = (0xCC, 0xCC, 0xD0)
PEPTIDE_PURPLE = (0x7D, 0x55, 0xC7)      # CCNY purple
INK = (20, 20, 20)
INK2 = (55, 55, 55)
MUTED = (100, 100, 100)
RULE = (200, 200, 200)
TIER_OK = (25, 25, 25)
TIER_NOISE = (150, 75, 0)                # ⛔ the WORDS carry the tier; this only reinforces them

MARGIN = 56        # left/right text margin inside a canvas
PAD = 30           # crop pad, white border only
MINSIZE = 22       # ⛔ refuse to draw below this rather than truncate
EDGE_TOL = 2
MIN_CANVAS_W = 1900   # pad narrow crops out to this so the band has room; see the loop

# ⛔ TWO caveat lines on the panel, not five. The band must stay a band: at the honest crop width
#    of these renders a five-line caveat set made the band nearly as tall as the picture, and a
#    figure that is mostly text is a figure nobody reads. The rest is in the README beside them.
CAVEATS = [
    "These are VIEWS, not evidence. A maintained geometry is NEVER evidence of affinity, nothing "
    "here has been tested, and every coordinate is a model.",
    "Contact margin is peptide contacts above this design's OWN shard floor, set by cross-paired "
    "non-binders. Numbers are carried from the scored tables, never re-measured here.",
]

LABELS = []


def fit(dr, text, path, start, maxw, minsize=MINSIZE):
    """⛔ Shrink until the line fits `maxw`. RAISE below `minsize` rather than truncate."""
    size = start
    while size >= minsize:
        f = ImageFont.truetype(path, size)
        if dr.textlength(text, font=f) <= maxw:
            return f, size
        size -= 2
    raise SystemExit("LABEL WILL NOT FIT at >=%dpx in %dpx: %r" % (minsize, maxw, text))


def put(dr, img, canvas_w, xy, text, path, start, maxw, fill=INK):
    f, size = fit(dr, text, path, start, maxw)
    dr.text(xy, text, font=f, fill=fill)
    LABELS.append({"img": img, "text": text, "font": path, "size": size,
                   "x": xy[0], "y": xy[1], "allowed": maxw, "canvas_w": canvas_w,
                   "measured": int(dr.textlength(text, font=f))})
    return size


def put_block(dr, img, canvas_w, xy, text, path, size, maxw, fill=INK):
    """Draw `text` WORD-WRAPPED to `maxw` at a fixed size; returns the height consumed.

    ⛔⛔ WHY WRAPPING AND NOT ONLY SHRINKING. `[M]` these panels crop to 1097 px wide -- the
    complex is a narrow diagonal and the union box is honestly that narrow -- and an 88-character
    title does not fit on one line even at the 22 px refusal floor, so `fit()` correctly REFUSED
    and stopped the run. Shrinking further would have produced a label nobody can read on a
    printed panel, which is the failure [[read-printed-labels-not-bar-lengths]] records. ⇒ the
    line WRAPS. ⛔ It is never truncated and never drawn below the floor: a word longer than the
    whole available width still raises, because that is a genuinely unfittable label.
    Every produced line is recorded separately, so the re-measure pass checks each one.
    """
    # ⛔ A design name is ONE 38-character token (`arm_shard03_MB_nomotif_broad_2_model_3`) and
    #    wrapping cannot help it, so the SIZE shrinks until the longest single token fits, and
    #    only then does the line wrap. Below MINSIZE it still raises rather than truncate.
    while size > MINSIZE:
        f = ImageFont.truetype(path, size)
        if max(dr.textlength(w, font=f) for w in text.split(" ")) <= maxw:
            break
        size -= 2
    f = ImageFont.truetype(path, size)
    words, lines, cur = text.split(" "), [], ""
    for w in words:
        trial = (cur + " " + w).strip()
        if cur and dr.textlength(trial, font=f) > maxw:
            lines.append(cur)
            cur = w
        else:
            cur = trial
        if dr.textlength(cur, font=f) > maxw:
            raise SystemExit("SINGLE TOKEN WILL NOT FIT in %dpx at %dpx: %r" % (maxw, size, cur))
    if cur:
        lines.append(cur)
    x, y = xy
    for ln in lines:
        dr.text((x, y), ln, font=f, fill=fill)
        LABELS.append({"img": img, "text": ln, "font": path, "size": size, "x": x, "y": y,
                       "allowed": maxw, "canvas_w": canvas_w,
                       "measured": int(dr.textlength(ln, font=f))})
        y += int(size * 1.28)
    return y - xy[1]


def legend_rows(dr, img, canvas_w, x0, y, avail, keys, start=34, sw=34, gap=16, sep=44):
    """Colour chips + labels, wrapping over as many rows as the keys need.

    ⛔⛔ The fit measures CHIPS, GAPS AND SEPARATORS, not just the text -- a text-only fit is how a
    6-key row ran 258 px off the edge. ⛔ And a key that does not fit WRAPS; it never overflows,
    because a key running off the edge is a wrong label.
    ⛔ Identity is carried by the TEXT of each key, never by colour alone.
    """
    def wrap(f):
        rows, cur, w = [], [], 0.0
        for c, t in keys:
            kw = sw + gap + dr.textlength(t, font=f) + sep
            if cur and w + kw > avail:
                rows.append(cur)
                cur, w = [], 0.0
            cur.append((c, t))
            w += kw
        if cur:
            rows.append(cur)
        return rows

    size = start
    while size >= MINSIZE:
        f = ImageFont.truetype(FONT, size)
        rows = wrap(f)
        if all(sum(sw + gap + dr.textlength(t, font=f) + sep for _c, t in r) <= avail
               for r in rows):
            break
        size -= 2
    else:
        raise SystemExit("LEGEND WILL NOT FIT even wrapped at >=%dpx" % MINSIZE)

    f = ImageFont.truetype(FONT, size)
    yy = y
    for row in rows:
        x = x0
        for c, t in row:
            dr.rectangle([x, yy + 5, x + sw, yy + 5 + size - 8], fill=c, outline=(80, 80, 80))
            x += sw + gap
            dr.text((x, yy), t, font=f, fill=INK2)
            w = dr.textlength(t, font=f)
            LABELS.append({"img": img, "text": t, "font": FONT, "size": size, "x": x, "y": yy,
                           "allowed": canvas_w - MARGIN - x, "canvas_w": canvas_w,
                           "measured": int(w)})
            x += w + sep
        yy += size + 14
    return yy - y


def union_box(paths):
    """⛔ The union of every panel's non-white content box -> it cannot clip structure in ANY of
    them, and it is the same box for all, so the panels stay pixel-comparable."""
    box = None
    size = None
    for p in paths:
        im = Image.open(p).convert("RGB")
        size = im.size
        bb = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
        if bb is None:
            raise SystemExit("panel is entirely blank: %s" % p)
        box = bb if box is None else (min(box[0], bb[0]), min(box[1], bb[1]),
                                      max(box[2], bb[2]), max(box[3], bb[3]))
    return (max(0, box[0] - PAD), max(0, box[1] - PAD),
            min(size[0], box[2] + PAD), min(size[1], box[3] + PAD))


def tier_line(p):
    if p["tier"] == "VALIDATED":
        return ("VALIDATED - clears the pre-registered conjunction at %d of %d seeds over 20 draws"
                % (p["seeds_clear"], p["seeds_total"]), TIER_OK)
    return ("%d of %d SEEDS - INSIDE NOISE, NOT PROMOTED (J1058 s6 tiering); the %s limb failed at "
            "1 of %d seeds" % (p["seeds_clear"], p["seeds_total"], p["failed_limb"],
                               p["seeds_total"]), TIER_NOISE)


def main(rawdir, outdir):
    man = json.load(open(os.path.join(rawdir, "j1063_renders_20261109.json")))
    if not man.get("ok"):
        raise SystemExit("REFUSING: the stage-2 manifest does not carry ok=true")
    panels = man["panels"]
    designs = [p for p in panels if p["kind"] == "design"]
    os.makedirs(outdir, exist_ok=True)

    paths = [os.path.join(rawdir, p["png"]) for p in panels]
    for p in paths:
        if not os.path.exists(p):
            raise SystemExit("missing raw render: %s" % p)
    box = union_box(paths)
    print("== j1063 stage 3/4: stamp %d panels ==" % len(panels))
    print("   common crop box (union of all %d content boxes, pad %d): %s -> %dx%d"
          % (len(paths), PAD, box, box[2] - box[0], box[3] - box[1]))

    crops = {}
    for p in panels:
        im = Image.open(os.path.join(rawdir, p["png"])).convert("RGB")
        crops[p["png"]] = im.crop(box)
        print("   %-58s %dx%d -> %dx%d" % (p["png"], im.size[0], im.size[1],
                                           crops[p["png"]].size[0], crops[p["png"]].size[1]))

    nval = sum(1 for p in designs if p["tier"] == "VALIDATED")
    nnoise = len(designs) - nval
    written = []

    for p in panels:
        im = crops[p["png"]]
        # ⛔ THE CROP CAN ONLY REMOVE BORDER; THIS ONLY ADDS IT. These renders crop honestly to
        #    1097 px because the complex is a narrow diagonal, and at that width the band wrapped
        #    to nearly the height of the picture. Padding the canvas out to a minimum width and
        #    CENTRING the image gives the band room without clipping a single pixel of structure.
        #    Every panel is padded by the same amount, so they stay pixel-comparable.
        if im.size[0] < MIN_CANVAS_W:
            pad = Image.new("RGB", (MIN_CANVAS_W, im.size[1]), (255, 255, 255))
            pad.paste(im, ((MIN_CANVAS_W - im.size[0]) // 2, 0))
            im = pad
        W, Himg = im.size
        avail = W - 2 * MARGIN
        design = p["kind"] == "design"
        out_name = ("00_overview.png" if not design
                    else "%02d_%s.png" % (p["slot"], p["design"]))
        if out_name == p["png"]:
            raise SystemExit("REFUSING: output name equals input name for %s" % p["png"])

        # ---- the band is SIZED FROM ITS CONTENT, never a constant a longer caption overruns
        probe = ImageDraw.Draw(Image.new("RGB", (W, 10)))
        if design:
            keys = [(TARGET_GREY, "HLA-A*11:01 + beta-2 microglobulin (one shared copy)"),
                    (PEPTIDE_PURPLE, "KRAS G12D 10-mer VVVGADGVGK (chain C)"),
                    (tuple(int(p["color_hex"][i:i + 2], 16) for i in (0, 2, 4)),
                     "this design, chain D (%d aa)" % p["binder_len"])]
        else:
            keys = [(TARGET_GREY, "HLA-A*11:01 + b2m (ONE shared copy)"),
                    (PEPTIDE_PURPLE, "KRAS G12D peptide")]
            keys += [(tuple(int(d["color_hex"][i:i + 2], 16) for i in (0, 2, 4)),
                      "%d. %s (%s)" % (d["slot"], d["sub_name"],
                                       "VALIDATED" if d["tier"] == "VALIDATED" else "3/4 seeds"))
                     for d in designs]
        # ⛔ Type sizes SCALE WITH THE CANVAS, so a narrow crop gets proportionate type rather
        #    than type tuned to one cohort's width. Floored at MINSIZE by fit()/put_block().
        # Ratios follow the deposited precedent (a 58 px title on a ~2900 px canvas in
        # `j957_stamp_panels`, 54 on ~2950 in `s94`), so the band reads the same as the
        # already-deposited challenge figures rather than being retuned by eye here.
        S_TITLE = max(MINSIZE + 8, W // 46)
        S_TIER = max(MINSIZE + 4, W // 60)
        S_BODY = max(MINSIZE, W // 76)
        S_SMALL = max(MINSIZE, W // 88)
        S_CAP = max(MINSIZE, W // 92)

        # ⛔ THE BAND IS SIZED FROM ITS OWN CONTENT, measured on a throwaway canvas first, never
        #    a constant that a longer caption set silently overruns. The probe's labels are
        #    discarded so they cannot pollute the re-measure pass.
        def lay(d2, nm, yy, draw):
            if design:
                t, tcol = tier_line(p)
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "SLOT %d of %d   %s" % (p["slot"], len(designs), p["design"]),
                           FONT_B, S_TITLE, avail) + 16
                yy += draw(d2, nm, W, (MARGIN, yy), t, FONT_B, S_TIER, avail, fill=tcol) + 16
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "submission name %s  |  generated in the %s  |  arm %s  |  MPNN rank "
                           "%d  |  binder %d aa" % (p["sub_name"], p["origin_round"], p["arm"],
                                                    p["rank"], p["binder_len"]),
                           FONT, S_BODY, avail, fill=INK2) + 8
                sd = ("  |  sd of per-seed peptide contacts %s" % p["sd_pep_contacts"]
                      if p["sd_pep_contacts"] else "")
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "mean contact margin %s over %d seeds  |  mean binder-peptide ipTM "
                           "%s%s" % (p["margin_mean"], p["seeds_total"], p["iptm_mean"], sd),
                           FONT, S_BODY, avail, fill=INK2) + 8
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "Shown: fold shard %s, model_0, i.e. that seed's OWN TOP-RANKED draw "
                           "of 5, not a random one. Superposed HLA-on-HLA onto slot 1 at %.3f A "
                           "over %d of %d CA." % (p["shard"], p["fit_rms"], p["fit_n"],
                                                  man["hla_len"]),
                           FONT, S_SMALL, avail, fill=MUTED) + 16
            else:
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "CHALLENGE 5 SUBMISSION SET   %d designs vs KRAS G12D on "
                           "HLA-A*11:01, one colour per design" % len(designs),
                           FONT_B, S_TITLE, avail) + 16
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "%d VALIDATED (4 of 4 seeds, 20 draws)  +  %d at 3 of 4 SEEDS - "
                           "INSIDE NOISE, NOT PROMOTED. These are NOT %d equally supported "
                           "designs." % (nval, nnoise, len(designs)),
                           FONT_B, S_TIER, avail, fill=TIER_NOISE) + 16
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "All %d binders in ONE frame against ONE copy of the target: the "
                           "slot-1 reference, with every other design superposed HLA-on-HLA onto "
                           "it at %.3f-%.3f A over %d-%d of %d CA. Each design's colour here is "
                           "the colour of its own panel, so this image is the legend for the set."
                           % (len(designs), min(d["fit_rms"] for d in designs[1:]),
                              max(d["fit_rms"] for d in designs),
                              min(d["fit_n"] for d in designs[1:]), man["hla_len"],
                              man["hla_len"]), FONT, S_BODY, avail, fill=INK2) + 8
                yy += draw(d2, nm, W, (MARGIN, yy),
                           "The binders overlap because they are drawn in one frame; the peptide "
                           "is drawn as cartoon AND sticks so it stays readable underneath them.",
                           FONT, S_SMALL, avail, fill=MUTED) + 16
            yy += legend_rows(d2, nm, W, MARGIN, yy, avail, keys, start=S_SMALL) + 10
            for cap in CAVEATS:
                yy += draw(d2, nm, W, (MARGIN, yy), cap, FONT, S_CAP, avail, fill=INK2) + 6
            return yy

        band = lay(probe, "_probe", 26, put_block) + 26
        LABELS[:] = [x for x in LABELS if x["img"] != "_probe"]

        canvas = Image.new("RGB", (W, Himg + band), (255, 255, 255))
        canvas.paste(im, (0, band))
        dr = ImageDraw.Draw(canvas)
        y = lay(dr, out_name, 26, put_block)
        if y > band:
            raise SystemExit("REFUSING: band content (%d px) overruns the %d px band on %s"
                             % (y, band, out_name))
        # inset rule: a full-width rule would put ink on the canvas edge and make the independent
        # pixel check ambiguous about what touched the edge
        dr.line([(MARGIN, band - 3), (W - MARGIN, band - 3)], fill=RULE, width=3)

        op = os.path.join(outdir, out_name)
        canvas.save(op, "PNG", dpi=(man["dpi"], man["dpi"]))
        written.append((op, out_name))
        print("   stamped %-56s %dx%d  band %d" % (out_name, canvas.size[0], canvas.size[1], band))

    # ---- ⛔ THE SEPARATE PASS: re-open the FINISHED files and re-measure EVERY recorded line
    print("\n== re-measuring %d recorded lines on the FINISHED files ==" % len(LABELS))
    over = []
    sizes = {}
    for L in LABELS:
        if L["img"] not in sizes:
            sizes[L["img"]] = Image.open(os.path.join(outdir, L["img"])).size
        cw, ch = sizes[L["img"]]
        if cw != L["canvas_w"]:
            raise SystemExit("%s: recorded canvas width %d != actual %d"
                             % (L["img"], L["canvas_w"], cw))
        f = ImageFont.truetype(L["font"], L["size"])
        scratch = ImageDraw.Draw(Image.new("RGB", (10, 10)))
        bb = scratch.textbbox((L["x"], L["y"]), L["text"], font=f)
        why = []
        if bb[2] > cw - EDGE_TOL:
            why.append("ink right %d > canvas %d" % (bb[2], cw))
        if scratch.textlength(L["text"], font=f) > L["allowed"]:
            why.append("advance > allotted %d" % L["allowed"])
        if bb[0] < EDGE_TOL:
            why.append("ink left %d" % bb[0])
        if bb[3] > ch - EDGE_TOL:
            why.append("ink bottom %d > canvas %d" % (bb[3], ch))
        if L["size"] < MINSIZE:
            why.append("font %d below the %d refusal floor" % (L["size"], MINSIZE))
        if why:
            over.append({"img": L["img"], "text": L["text"][:70], "why": "; ".join(why)})

    print("   images covered: %d  (%s)" % (len(sizes), ", ".join(sorted(sizes))))
    for o in over:
        print("   OVERFLOW %-24s %s  <-- %s" % (o["img"], o["text"], o["why"]))
    print("\n   OVERFLOWING LABELS: %d of %d" % (len(over), len(LABELS)))

    # ---- independent pixel check, which does not trust the text metrics at all
    bad = 0
    for img in sorted(sizes):
        im = Image.open(os.path.join(outdir, img)).convert("RGB")
        w, h = im.size
        for name, bx in (("left", (0, 0, 4, h)), ("right", (w - 4, 0, w, h))):
            nz = sum(1 for px in im.crop(bx).getdata() if px != (255, 255, 255))
            if nz:
                bad += 1
                print("   %-24s %s edge: %d non-white px" % (img, name, nz))
    if bad == 0:
        print("   0 images with ink touching a left or right edge")

    lp = os.path.join(outdir, "j1063_labels_20261109.json")
    json.dump({"margin": MARGIN, "minsize": MINSIZE, "crop_box": list(box),
               "n_labels": len(LABELS), "overflows": over, "edge_touching": bad,
               "labels": LABELS}, open(lp, "w"), indent=1)
    print("   wrote %s" % lp)
    if over or bad:
        raise SystemExit("REFUSING TO CLAIM SUCCESS: %d overflowing labels, %d edge-touching"
                         % (len(over), bad))
    print("== j1063 stage 3 ok: %d stamped panels ==" % len(written))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
