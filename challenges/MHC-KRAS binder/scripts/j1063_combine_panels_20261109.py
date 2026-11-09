#!/usr/bin/env python3
"""j1063 stage 4/4 -- combine the stamped panels into ONE submission figure.

    python3 -I j1063_combine_panels_20261109.py <stampdir> <out.png> [cols] [width]

⭐ METHOD AND CREDIT: `an earlier challenge's scripts/j582_combine_panels_20261019.py`. ⛔ ADDITIVE: it
   reads the stamped PNGs and writes a new file; it touches no earlier stage and no source image.

⭐ EACH INPUT ALREADY CARRIES ITS OWN STAMPED BAND, so every design's TIER WORD, numbers and
   colour key survive the combination. ⛔ The master header does NOT replace them -- that is the
   whole point of combining already-titled panels rather than bare renders: a reader who crops one
   cell out of this figure still gets a cell that says what it is.

⛔⛔ AND THE HEADER MUST STATE THE 4/3 SPLIT. Seven panels side by side is exactly the image that
   could be mistaken for seven equally supported designs, which is the one misreading this figure
   must not permit. James McCann putting all seven in the submission is a SCOPE decision and does not make
   the 3 that clear only 3 of 4 seeds into results.

⛔ Every header line goes through the same fit-and-record discipline as stage 3 and is
   RE-MEASURED on the FINISHED file. A stamped line that runs off the edge is a WRONG label.
⛔ Refuses to overwrite: a fresh `--out` per run (STOP S5).
"""
import json
import os
import sys

from PIL import Image, ImageDraw, ImageFont

FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONT_B = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
CCNY = (0x7D, 0x55, 0xC7)
INK2 = (55, 55, 55)
TIER_NOISE = (150, 75, 0)
RULE = (205, 205, 205)
MARGIN = 70
MINSIZE = 22
GAP = 40
EDGE_TOL = 2

LABELS = []


def fit(dr, text, path, start, maxw, minsize=MINSIZE):
    size = start
    while size >= minsize:
        f = ImageFont.truetype(path, size)
        if dr.textlength(text, font=f) <= maxw:
            return f, size
        size -= 2
    raise SystemExit("HEADER LINE WILL NOT FIT at >=%dpx in %dpx: %r" % (minsize, maxw, text))


def put(dr, img, canvas_w, xy, text, path, start, maxw, fill=INK2):
    f, size = fit(dr, text, path, start, maxw)
    dr.text(xy, text, font=f, fill=fill)
    LABELS.append({"img": img, "text": text, "font": path, "size": size, "x": xy[0], "y": xy[1],
                   "allowed": maxw, "canvas_w": canvas_w})
    return int(size * 1.45)


def main(stampdir, out, cols=4, width=6000):
    if os.path.exists(out):
        raise SystemExit("exists, refusing to overwrite: %s  (fresh --out per run, STOP S5)" % out)
    lab = json.load(open(os.path.join(stampdir, "j1063_labels_20261109.json")))
    if lab["overflows"] or lab["edge_touching"]:
        raise SystemExit("REFUSING: stage 3 reported %d overflowing labels and %d edge-touching "
                         "images; combining them would only multiply a wrong label."
                         % (len(lab["overflows"]), lab["edge_touching"]))

    names = sorted(n for n in os.listdir(stampdir) if n.endswith(".png"))
    if len(names) != 8:
        raise SystemExit("REFUSING: expected 8 stamped panels (overview + 7), found %d: %s"
                         % (len(names), names))
    if names[0] != "00_overview.png":
        raise SystemExit("REFUSING: the first cell must be the overview, found %s" % names[0])
    print("== j1063 stage 4/4: combine %d panels ==" % len(names))

    colw = (width - GAP * (cols - 1) - 2 * MARGIN) // cols
    imgs = []
    for n in names:
        im = Image.open(os.path.join(stampdir, n)).convert("RGB")
        # ⛔ trim the outer white margin first, or the grid inherits each panel's padding
        bb = im.convert("L").point(lambda v: 0 if v > 247 else 255).getbbox()
        if bb:
            im = im.crop((max(0, bb[0] - 10), max(0, bb[1] - 10),
                          min(im.size[0], bb[2] + 10), min(im.size[1], bb[3] + 10)))
        h = round(im.size[1] * colw / float(im.size[0]))
        imgs.append((n, im.resize((colw, h), Image.LANCZOS)))
        print("   %-50s -> %dx%d" % (n, colw, h))

    hpx = max(MINSIZE + 20, width // 58)
    out_name = os.path.basename(out)
    probe = ImageDraw.Draw(Image.new("RGB", (width, 10)))
    lines = [
        ("CCNY McCann / Koder -- Challenge 5 submission set: 7 de novo binders against the "
         "KRAS G12D 10-mer on HLA-A*11:01", FONT_B, hpx, CCNY),
        ("4 are VALIDATED (clear the pre-registered conjunction at 4 of 4 seeds over 20 draws). "
         "3 clear only 3 of 4 seeds and are INSIDE NOISE, NOT PROMOTED. Each panel states its own "
         "tier on its face.", FONT_B, int(hpx * 0.62), TIER_NOISE),
        ("One colour per design, consistent between each design's own panel and the overview, so "
         "the overview is the legend. Grey = HLA-A*11:01 + beta-2 microglobulin, ONE shared copy. "
         "CCNY purple = the KRAS G12D peptide.", FONT, int(hpx * 0.52), INK2),
        ("These are VIEWS, not evidence: a maintained geometry is NEVER evidence of affinity and "
         "nothing here has been tested. Panel 00 is the overview; panels 01-07 are one design "
         "each.", FONT, int(hpx * 0.52), INK2),
    ]
    band = MARGIN
    for t, f, s, _c in lines:
        band += put(probe, "_probe", width, (MARGIN, 0), t, f, s, width - 2 * MARGIN)
    LABELS[:] = [x for x in LABELS if x["img"] != "_probe"]
    band += 30

    rows = [imgs[i:i + cols] for i in range(0, len(imgs), cols)]
    rowh = [max(im.size[1] for _n, im in r) for r in rows]
    total = band + sum(rowh) + GAP * (len(rows) - 1) + MARGIN
    canvas = Image.new("RGB", (width, total), (255, 255, 255))
    dr = ImageDraw.Draw(canvas)
    y = MARGIN
    for t, f, s, c in lines:
        y += put(dr, out_name, width, (MARGIN, y), t, f, s, width - 2 * MARGIN, fill=c)
    dr.line([(MARGIN, band - 16), (width - MARGIN, band - 16)], fill=RULE, width=4)

    y = band
    for ri, row in enumerate(rows):
        for ci, (_n, im) in enumerate(row):
            canvas.paste(im, (MARGIN + ci * (colw + GAP), y))
        y += rowh[ri]
        if ri < len(rows) - 1:
            dr.line([(MARGIN, y + GAP // 2), (width - MARGIN, y + GAP // 2)],
                    fill=RULE, width=3)
            y += GAP
    canvas.save(out, "PNG", dpi=(400, 400))
    print("   wrote %s  (%dx%d, %d cells in %d cols)"
          % (out, canvas.size[0], canvas.size[1], len(imgs), cols))

    # ---- ⛔ re-measure the header on the FINISHED file
    over = []
    im = Image.open(out)
    for L in LABELS:
        f = ImageFont.truetype(L["font"], L["size"])
        bb = ImageDraw.Draw(Image.new("RGB", (10, 10))).textbbox((L["x"], L["y"]), L["text"],
                                                                 font=f)
        if bb[2] > im.size[0] - EDGE_TOL or bb[0] < EDGE_TOL or L["size"] < MINSIZE:
            over.append({"img": L["img"], "text": L["text"][:70]})
    print("   OVERFLOWING LABELS: %d of %d" % (len(over), len(LABELS)))
    w, h = im.size
    rgb = im.convert("RGB")
    bad = sum(1 for name, bx in (("left", (0, 0, 4, h)), ("right", (w - 4, 0, w, h)))
              if any(px != (255, 255, 255) for px in rgb.crop(bx).getdata()))
    print("   %d of 2 vertical edges carry ink" % bad)
    if over or bad:
        raise SystemExit("REFUSING TO CLAIM SUCCESS: %d overflows, %d edges with ink"
                         % (len(over), bad))
    print("== j1063 stage 4 ok ==")


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], int(a[2]) if len(a) > 2 else 4, int(a[3]) if len(a) > 3 else 6000)
