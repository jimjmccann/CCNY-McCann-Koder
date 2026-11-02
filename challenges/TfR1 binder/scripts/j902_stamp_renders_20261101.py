#!/usr/bin/env python3
"""j902 stamp -- autocrop the white margins off each render and stamp a TITLE BAND.

Local, PIL only. ⛔ Opens no structure and changes no geometry: this is presentation.

⭐ A LABEL GOES IN A TITLE BAND, not into the image area: an automatic labeller that strips the
common prefix and suffix cannot carry a constant marker. The band carries the design id, what
the colours mean, and the design's own measured numbers, so a PNG that is detached from its
surrounding text still says what it is.

⛔ The crop is computed from the WHITE BACKGROUND only and is applied with a fixed pad, so it can
only ever remove uniform-white border -- it cannot clip structure. The pre/post sizes are PRINTED
so a bad crop is visible rather than silent.

usage: python3 -I j902_stamp_renders_20261101.py <renderdir>
"""
import json
import os
import sys

from PIL import Image, ImageDraw, ImageFont

FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
FONT_B = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
PURPLE = (125, 85, 199)
GREY = (150, 150, 150)
YELLOW = (210, 165, 30)
TEAL = (95, 179, 179)
PAD = 40


def autocrop(im, pad=PAD):
    """remove uniform WHITE border only; returns (image, before, after)"""
    rgb = im.convert("RGB")
    bg = Image.new("RGB", rgb.size, (255, 255, 255))
    from PIL import ImageChops
    diff = ImageChops.difference(rgb, bg)
    bb = diff.getbbox()
    before = im.size
    if not bb:
        return im, before, before
    l, t, r, b = bb
    l = max(0, l - pad); t = max(0, t - pad)
    r = min(im.size[0], r + pad); b = min(im.size[1], b + pad)
    out = im.crop((l, t, r, b))
    return out, before, out.size


def fit(dr, text, path, start, maxw, minsize=18):
    """⛔ a stamped label that runs off the edge is a WRONG label, not a cosmetic issue -- the
    first pass truncated 'rg_ratio' to 'rg_'. Shrink until it fits, and REFUSE below minsize."""
    s = start
    while s > minsize:
        f = ImageFont.truetype(path, s)
        if dr.textlength(text, font=f) <= maxw:
            return f
        s -= 2
    return ImageFont.truetype(path, minsize)


def main(d):
    meta = {m["png"]: m for m in json.load(open(os.path.join(d, "j902_renders.json")))}
    outd = os.path.join(d, "labelled")
    os.makedirs(outd, exist_ok=True)
    print("== j902 stamp ==")
    for png in sorted(meta):
        p = os.path.join(d, png)
        if not os.path.exists(p):
            print("  SKIP %s (missing)" % png)
            continue
        m = meta[png]
        im = Image.open(p).convert("RGB")
        im, before, after = autocrop(im)
        W = im.size[0]
        band_h = 190 if m["design"] else 150
        canvas = Image.new("RGB", (W, im.size[1] + band_h), (255, 255, 255))
        canvas.paste(im, (0, band_h))
        dr = ImageDraw.Draw(canvas)
        fb = ImageFont.truetype(FONT_B, 46)
        fs = ImageFont.truetype(FONT, 34)
        fsm = ImageFont.truetype(FONT, 30)

        if m["design"]:
            t = "%s   (%s)" % (m["obj"], m["opt"])
            dr.text((PAD, 18), t, font=fit(dr, t, FONT_B, 46, W - 2 * PAD), fill=(20, 20, 20))
            dr.text((PAD, 76), m["why"],
                    font=fit(dr, m["why"], FONT, 30, W - 2 * PAD), fill=(90, 90, 90))
            line = ("max CA-CA step %s A   |   backbone min to target %s A   |   "
                    "min_cut %s   |   rg_ratio %s   |   res under 2.5 A: %s"
                    % (m["max_step"], m["bb_min"], m["min_cut"], m["rg_ratio"],
                       m["n_under_2p5"]))
            fline = fit(dr, line, FONT, 30, W - 2 * PAD)
            dr.text((PAD, 118), line, font=fline, fill=(40, 40, 40))
            swatches = [(PURPLE, "given 5XOX core (%d res)" % m["n_held"]),
                        (GREY, "RFd3 BUILT (%d res)" % m["n_built"]),
                        (YELLOW, "4 fixed pins"), (TEAL, "TfR1")]
            legend = "  ".join(s[1] for s in swatches)
            fleg = fit(dr, legend + "          ", FONT, 30, W - 2 * PAD)
            x = PAD
            for c, lab in swatches:
                dr.rectangle([x, 156, x + 28, 180], fill=c, outline=(80, 80, 80))
                dr.text((x + 38, 154), lab, font=fleg, fill=(40, 40, 40))
                x += 38 + int(dr.textlength(lab, font=fleg)) + 30
        else:
            t = "OVERVIEW -- 6 designs, 6 distinct batches, arms opt05 / opt01 / opt04"
            dr.text((PAD, 18), t, font=fit(dr, t, FONT_B, 46, W - 2 * PAD), fill=(20, 20, 20))
            dr.text((PAD, 76), "Not a binding result: no co-fold, no sequence design yet.",
                    font=fsm, fill=(90, 90, 90))
            swatches = [(PURPLE, "given 5XOX core"), (GREY, "RFd3 BUILT"),
                        (YELLOW, "fixed pins"), (TEAL, "TfR1")]
            x = PAD
            for c, lab in swatches:
                dr.rectangle([x, 116, x + 28, 140], fill=c, outline=(80, 80, 80))
                dr.text((x + 38, 114), lab, font=fsm, fill=(40, 40, 40))
                x += 38 + int(dr.textlength(lab, font=fsm)) + 46
        dr.line([(0, band_h - 2), (W, band_h - 2)], fill=(200, 200, 200), width=3)

        op = os.path.join(outd, png)
        canvas.save(op, "PNG")
        print("  %-22s %sx%s -> crop %sx%s -> %sx%s" %
              (png, before[0], before[1], after[0], after[1],
               canvas.size[0], canvas.size[1]))
    print("\n   labelled PNGs -> %s" % outd)


main(sys.argv[1])
