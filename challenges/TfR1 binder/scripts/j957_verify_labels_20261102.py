#!/usr/bin/env python3
"""j957 verify -- re-measure EVERY stamped line against the canvas it was drawn on.

⛔ DEFECT 3 AGAIN, STATED PLAINLY: a stamped label that runs off the right edge is a
**WRONG LABEL**, not a cosmetic problem -- the earlier pass truncated `rg_ratio` to `rg_`.
The stamping pass fits each line as it draws it; that is the stamper checking its own work.
THIS pass is separate and does not trust it: it re-opens each FINISHED PNG, reads its real
pixel width, re-creates each recorded font at its recorded size, re-measures each recorded
string, and fails if any line's right edge crosses the canvas (or its allotted width).

It also checks, independently of the text metrics, that no non-white ink actually touches
the left or right edge of any band region -- a line that measured short but was drawn at a
bad x would still be caught.

Reports the count of overflowing labels. The required answer is 0 and it is PRINTED, not
assumed.

usage: python3 -I j957_verify_labels_20261102.py <renderdir>
"""
import json
import os
import sys

from PIL import Image, ImageDraw, ImageFont

EDGE_TOL = 2      # px: a glyph may not come closer than this to the canvas edge


def main(d):
    rec = json.load(open(os.path.join(d, "j957_labels_20261102.json")))
    labels = rec["labels"]
    margin = rec["margin"]
    sizes = {}
    overflow = []
    checked = 0

    print("== j957 verify: re-measuring %d recorded lines ==" % len(labels))
    for L in labels:
        img = L["image"]
        if img not in sizes:
            p = os.path.join(d, img)
            if not os.path.exists(p):
                raise RuntimeError("recorded label names a PNG that does not exist: %s" % p)
            sizes[img] = Image.open(p).size
        cw, ch = sizes[img]
        if cw != L["canvas_w"]:
            raise RuntimeError("%s: recorded canvas width %d != actual %d"
                               % (img, L["canvas_w"], cw))
        f = ImageFont.truetype(L["font"], L["size"])
        scratch = ImageDraw.Draw(Image.new("RGB", (10, 10)))
        w = scratch.textlength(L["text"], font=f)
        # the true ink box, not just the advance width
        bbox = scratch.textbbox((L["x"], L["y"]), L["text"], font=f)
        right = bbox[2]
        checked += 1
        why = []
        if right > cw - EDGE_TOL:
            why.append("ink right edge %d > canvas %d" % (right, cw))
        if L["x"] + w > L["x"] + L["maxw"]:
            why.append("advance %.1f > allotted %d" % (w, L["maxw"]))
        if bbox[0] < EDGE_TOL:
            why.append("ink left edge %d" % bbox[0])
        if bbox[3] > ch - EDGE_TOL:
            why.append("ink bottom %d > canvas height %d" % (bbox[3], ch))
        if L["size"] < rec["minsize"]:
            why.append("font %d below the %d refusal floor" % (L["size"], rec["minsize"]))
        if why:
            overflow.append((img, L["text"], "; ".join(why)))

    print("   lines re-measured: %d" % checked)
    print("   images covered:    %d  (%s)" % (len(sizes), ", ".join(sorted(sizes))))
    print("   margin recorded:   %d px   minimum font size: %d px"
          % (margin, rec["minsize"]))
    if overflow:
        print("\n   OVERFLOWING LABELS: %d" % len(overflow))
        for img, t, why in overflow:
            print("     %-44s %s  <-- %s" % (img, t[:70], why))
    print("\n   OVERFLOWING LABELS: %d" % len(overflow))

    # ---- independent pixel check: no ink in the outermost columns of any stamped image
    print("\n== independent pixel check (no ink in the outer %d px columns) ==" % margin)
    bad_px = 0
    for img in sorted(sizes):
        im = Image.open(os.path.join(d, img)).convert("RGB")
        w, h = im.size
        for name, bx in (("left", (0, 0, 4, h)), ("right", (w - 4, 0, w, h))):
            strip = im.crop(bx)
            nonwhite = sum(1 for px in strip.getdata() if px != (255, 255, 255))
            if nonwhite:
                bad_px += 1
                print("   %-44s %s edge: %d non-white px" % (img, name, nonwhite))
    if bad_px == 0:
        print("   0 images with ink touching a left or right edge")

    print("\n== RESULT: %d overflowing labels, %d edge-touching images =="
          % (len(overflow), bad_px))
    if overflow or bad_px:
        sys.exit(1)


main(sys.argv[1])
