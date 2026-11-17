"""The fit-point images of the Bone Curve overlay (resources/points), drawn 8x
oversize and reduced, so the discs are antialiased.

    py -3.10 tools/make_points.py
"""

import os

from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(HERE, 'resources', 'points')
BLUE = (0, 110, 220, 255)
WHITE = (255, 255, 255, 255)


def disc(path, size, fill, ring, ring_px):
    k = 8
    im = Image.new('RGBA', (size * k, size * k), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    d.ellipse([0, 0, size * k - 1, size * k - 1], fill=ring)
    w = int(round(ring_px * k))
    d.ellipse([w, w, size * k - 1 - w, size * k - 1 - w], fill=fill)
    im.resize((size, size), Image.LANCZOS).save(path)


if __name__ == '__main__':
    os.makedirs(OUT, exist_ok=True)
    disc(os.path.join(OUT, 'fit.png'), 11, WHITE, BLUE, 2.0)
    disc(os.path.join(OUT, 'hot.png'), 13, BLUE, WHITE, 1.5)
    print('written to', OUT)
