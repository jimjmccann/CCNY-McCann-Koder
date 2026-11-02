#!/usr/bin/env python3
"""j957 stamp -- crop the j957 raw renders to ONE COMMON BOX, stamp a title band on each,
and build the combined 2x2 panel figure.

Local, PIL only. ⛔ Opens no structure and changes no geometry: this is presentation.

⭐ A LABEL GOES IN A BAND, never into the image area (j902's recorded lesson: an automatic
labeller that strips a common prefix and suffix cannot carry a constant marker).

⛔ DEFECT 3 -- A STAMPED LINE THAT RUNS OFF THE RIGHT EDGE IS A **WRONG LABEL**, not a
cosmetic problem: the first pass truncated the column name `rg_ratio` to `rg_`. So:
  (a) every single piece of text drawn here goes through `fit()`, which shrinks until the
      line fits the width actually available to it and RAISES below a minimum size rather
      than silently drawing something truncated; and
  (b) every drawn line is recorded to j957_labels_20261102.json with its text, font, size,
      x, baseline y and the width it was allowed, so a SEPARATE pass
      (j957_verify_labels_20261102.py) can re-open the finished PNGs and re-measure every
      line against its own canvas. That pass, not this one, is what reports 0 overflows.

⛔ THE CROP CAN ONLY REMOVE BORDER. The box is the union of the white-background content
boxes of all four renders, so it cannot clip structure in any of them, and it is the SAME
box for every panel so the panels stay pixel-comparable. Before/after sizes are printed.

usage: python3 -I j957_stamp_panels_20261102.py <renderdir>
"""
import json
import os
import sys

from PIL import Image, ImageChops, ImageDraw, ImageFont

FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONT_B = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"

# the session's own colours, which ARE the information content of the figure
PURPLE = (0x7D, 0x55, 0xC7)
GREY = (0xBF, 0xBF, 0xBF)
YELLOW = (0xF2, 0xC2, 0x30)
TEAL = (0x5F, 0xB3, 0xB3)

INK = (20, 20, 20)
INK2 = (45, 45, 45)
MUTED = (95, 95, 95)
RULE = (200, 200, 200)

PAD = 48          # crop pad, white border only
MARGIN = 48       # left/right text margin inside a canvas
MINSIZE = 20      # ⛔ refuse to draw below this rather than truncate

LEGEND = [(PURPLE, "purple = GIVEN protein, sequence-fixed (143 res)"),
          (GREY, "grey = BUILT N-terminal extension (27 res)"),
          (YELLOW, "yellow = fixed guidepost pin (resi 26)"),
          (TEAL, "teal = TfR1 target")]

CAPTION_1 = ("Ranks 1 and 3 come from the SAME batch (s01b0), so the three panels are NOT "
             "three independent draws.")
CAPTION_2 = "Engagement is not affinity, and every coordinate here is a model."
# ⭐ added 2026-11-02 after reading the render: at high surface transparency the built
# segment is visibly INSIDE the target's solvent-excluded surface. That is what a 2.1-2.6 A
# minimum separation looks like, and these are RAW generator outputs gated at a 2.0 A floor
# with no relaxation. Saying so on the figure is cheaper than a reader inferring packing.
CAPTION_3 = ("These are RAW generator backbones gated at a 2.0 A minimum separation and NOT "
             "relaxed; the built segment sits within the target's surface envelope at 2.1-2.6 A.")

LABELS = []       # every stamped line, for the independent re-measurement pass


def fit(dr, text, path, start, maxw, minsize=MINSIZE):
    """⛔ shrink until the line fits `maxw`; RAISE below `minsize` rather than truncate."""
    s = start
    while s >= minsize:
        f = ImageFont.truetype(path, s)
        if dr.textlength(text, font=f) <= maxw:
            return f
        s -= 2
    raise RuntimeError("LABEL WILL NOT FIT at >=%dpx in %dpx: %r" % (minsize, maxw, text))


def put(dr, img_name, canvas_w, xy, text, path, start, maxw, fill=INK):
    """draw one line, fitted, and RECORD it so a separate pass can re-measure it"""
    f = fit(dr, text, path, start, maxw)
    dr.text(xy, text, font=f, fill=fill)
    w = dr.textlength(text, font=f)
    LABELS.append({"image": img_name, "text": text, "font": path, "size": f.size,
                   "x": xy[0], "y": xy[1], "maxw": maxw, "measured_w": w,
                   "canvas_w": canvas_w})
    return f


def legend_row(dr, img_name, canvas_w, x0, y, avail, start, sw=30, gap=16, sep=40):
    """lay out the colour legend; the WHOLE row is fitted before anything is drawn"""
    s = start
    while s >= MINSIZE:
        f = ImageFont.truetype(FONT, s)
        total = 0
        for i, (_c, lab) in enumerate(LEGEND):
            total += sw + gap + dr.textlength(lab, font=f)
            if i < len(LEGEND) - 1:
                total += sep
        if total <= avail:
            break
        s -= 2
    else:
        raise RuntimeError("LEGEND ROW WILL NOT FIT at >=%dpx in %dpx" % (MINSIZE, avail))
    f = ImageFont.truetype(FONT, s)
    x = x0
    for c, lab in LEGEND:
        dr.rectangle([x, y + 4, x + sw, y + 4 + sw], fill=c, outline=(80, 80, 80))
        tx = x + sw + gap
        dr.text((tx, y), lab, font=f, fill=INK2)
        w = dr.textlength(lab, font=f)
        LABELS.append({"image": img_name, "text": lab, "font": FONT, "size": s,
                       "x": tx, "y": y, "maxw": canvas_w - MARGIN - tx,
                       "measured_w": w, "canvas_w": canvas_w})
        x = tx + w + sep
    return s


def union_box(paths, pad=PAD):
    """union of the white-background content boxes -> ONE crop box for every panel"""
    L = T = 10 ** 9
    R = B = 0
    size = None
    for p in paths:
        im = Image.open(p).convert("RGB")
        size = im.size
        bb = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).getbbox()
        if bb is None:
            continue
        L, T = min(L, bb[0]), min(T, bb[1])
        R, B = max(R, bb[2]), max(B, bb[3])
    L = max(0, L - pad); T = max(0, T - pad)
    R = min(size[0], R + pad); B = min(size[1], B + pad)
    return (L, T, R, B)


def main(d):
    man = json.load(open(os.path.join(d, "j957_renders_20261102.json")))
    panels = man["panels"]
    paths = [os.path.join(d, p["png"]) for p in panels]
    for p in paths:
        if not os.path.exists(p):
            raise RuntimeError("missing raw render: %s" % p)

    box = union_box(paths)
    print("== j957 stamp ==")
    print("   common crop box (union of all 4 content boxes, pad %d): %s -> %dx%d"
          % (PAD, box, box[2] - box[0], box[3] - box[1]))

    crops = {}
    for p in panels:
        im = Image.open(os.path.join(d, p["png"])).convert("RGB")
        c = im.crop(box)
        crops[p["png"]] = c
        print("   %-34s %dx%d -> %dx%d" % (p["png"], im.size[0], im.size[1],
                                           c.size[0], c.size[1]))

    written = []

    # ---------------- one labelled panel per render ----------------
    for p in panels:
        im = crops[p["png"]]
        W = im.size[0]
        avail = W - 2 * MARGIN
        band = 250
        out_name = p["png"].replace("j957_panel_", "j957_labelled_")
        canvas = Image.new("RGB", (W, im.size[1] + band), (255, 255, 255))
        canvas.paste(im, (0, band))
        dr = ImageDraw.Draw(canvas)

        if p["kind"] == "design":
            t = "RANK %d   %s" % (p["rank"], p["design"])
            put(dr, out_name, W, (MARGIN, 20), t, FONT_B, 58, avail)
            sub = ("batch %s   |   BUILT 27-residue N-terminal segment: %d heavy-atom "
                   "contact pairs under 4.5 A to TfR1" % (p["batch"], p["pairs"]))
            put(dr, out_name, W, (MARGIN, 92), sub, FONT, 38, avail, fill=INK2)
        else:
            t = "OVERVIEW   all three designs, one TfR1 surface, same camera"
            put(dr, out_name, W, (MARGIN, 20), t, FONT_B, 58, avail)
            sub = ("ranks 1/2/3 = %s / %s / %s built-segment contact pairs under 4.5 A; "
                   "batches s01b0 / s03b2 / s01b0"
                   % (panels[0]["pairs"], panels[1]["pairs"], panels[2]["pairs"]))
            put(dr, out_name, W, (MARGIN, 92), sub, FONT, 38, avail, fill=INK2)

        legend_row(dr, out_name, W, MARGIN, 150, avail, 34)
        put(dr, out_name, W, (MARGIN, 200), CAPTION_2, FONT, 32, avail, fill=MUTED)
        # inset to the margin: a full-width rule would put ink on the canvas edge and make
        # the verifier's independent pixel check ambiguous about what touched the edge
        dr.line([(MARGIN, band - 2), (W - MARGIN, band - 2)], fill=RULE, width=3)

        op = os.path.join(d, out_name)
        canvas.save(op, "PNG", dpi=(man["dpi"], man["dpi"]))
        written.append(op)
        print("   labelled %-34s %dx%d" % (out_name, canvas.size[0], canvas.size[1]))

    # ---------------- the combined 2x2 figure ----------------
    cell_w = 1600
    scale = cell_w / float(box[2] - box[0])
    cell_h = int(round((box[3] - box[1]) * scale))
    cap_h = 104                       # per-cell caption strip above each panel
    gut = 44
    marg = 60
    head = 300
    foot = 400   # +70: CAPTION_3 added a footer line (the raw-structure caveat)
    CW = marg * 2 + cell_w * 2 + gut
    CH = head + 2 * (cap_h + cell_h) + gut + foot
    comb_name = "j957_combined_nterm_top3_20261102.png"
    comb = Image.new("RGB", (CW, CH), (255, 255, 255))
    dr = ImageDraw.Draw(comb)
    avail = CW - 2 * marg

    put(dr, comb_name, CW, (marg, 34),
        "RFd3 N-terminal extension designs vs TfR1 -- corrected top 3 of 66 survivors",
        FONT_B, 76, avail)
    put(dr, comb_name, CW, (marg, 130),
        "Ranked on heavy-atom contact pairs under 4.5 A made by the BUILT 27-residue "
        "N-terminal segment. Common camera; the target is identical in every panel.",
        FONT, 40, avail, fill=INK2)
    legend_row(dr, comb_name, CW, marg, 200, avail, 40, sw=36, gap=18, sep=56)
    dr.line([(marg, head - 24), (CW - marg, head - 24)], fill=RULE, width=3)

    order = [panels[0], panels[1], panels[2], panels[3]]
    for i, p in enumerate(order):
        col, row = i % 2, i // 2
        x = marg + col * (cell_w + gut)
        y = head + row * (cap_h + cell_h + gut)
        if p["kind"] == "design":
            c1 = "RANK %d   %s" % (p["rank"], p["design"])
            c2 = ("batch %s   |   built-segment contact pairs < 4.5 A: %d"
                  % (p["batch"], p["pairs"]))
        else:
            c1 = "OVERVIEW   all three, one TfR1 surface"
            c2 = "same camera as the three panels above/left"
        put(dr, comb_name, CW, (x, y), c1, FONT_B, 46, cell_w)
        put(dr, comb_name, CW, (x, y + 56), c2, FONT, 36, cell_w, fill=INK2)
        comb.paste(crops[p["png"]].resize((cell_w, cell_h), Image.LANCZOS),
                   (x, y + cap_h))

    fy = head + 2 * (cap_h + cell_h) + gut + 20
    dr.line([(marg, fy - 26), (CW - marg, fy - 26)], fill=RULE, width=3)
    put(dr, comb_name, CW, (marg, fy), CAPTION_1, FONT_B, 44, avail)
    put(dr, comb_name, CW, (marg, fy + 66), CAPTION_2, FONT, 42, avail, fill=INK2)
    put(dr, comb_name, CW, (marg, fy + 128), CAPTION_3, FONT, 42, avail, fill=INK2)
    put(dr, comb_name, CW, (marg, fy + 190),
        "Purple is the protein GIVEN to the generator and sequence-fixed; grey is what the "
        "generator BUILT; yellow is the one fixed guidepost; teal is the TfR1 target surface.",
        FONT, 36, avail, fill=MUTED)
    put(dr, comb_name, CW, (marg, fy + 244),
        "Source session: j956_nterm_top3_corrected_20261102.pse. Surface transparency is a "
        "presentation choice; the cut-face interior colour is cosmetic, not data.",
        FONT, 32, avail, fill=MUTED)

    cp = os.path.join(d, comb_name)
    comb.save(cp, "PNG", dpi=(man["dpi"], man["dpi"]))
    written.append(cp)
    print("   combined %-34s %dx%d" % (comb_name, comb.size[0], comb.size[1]))

    lp = os.path.join(d, "j957_labels_20261102.json")
    json.dump({"margin": MARGIN, "minsize": MINSIZE, "labels": LABELS},
              open(lp, "w"), indent=1)
    print("\n   %d stamped lines recorded -> %s" % (len(LABELS), lp))
    for w in written:
        print("   wrote %s (%d bytes)" % (w, os.path.getsize(w)))


main(sys.argv[1])
