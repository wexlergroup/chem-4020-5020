"""Figure style for the CHEM 4020/5020 notes.

Every figure in the notes is drawn by a function ``draw(c)`` that returns a Matplotlib figure. ``show(draw)``
calls it twice, once with the colors for the light page and once with the colors for the dark page, and the
book displays the version that matches the reader's theme.

    import sys
    sys.path.append("..")            # the book folder, which holds figstyle
    import figstyle as fs

    def draw(c):
        fig, ax = plt.subplots(figsize=(6, 4))     # inches; fs.WIDTH is the full width of the text
        ax.plot(x, y, color=c.blue)
        return fig

    fs.show(draw, alt="One sentence that says what the figure shows.")

Colors are named by role and belong to the argument ``c``: ``c.blue``, ``c.vermilion`` and ``c.green`` for data
series (in that order of use), ``c.gray`` for reference curves, ``c.ink`` for text and annotations, ``c.muted``
for secondary text and axes, ``c.fill`` for quiet boxes, ``c.ground`` for the page behind the figure, and
``c.solid``, ``c.liquid``, ``c.gas``, ``c.carbon``, ``c.oxygen`` for phases and atoms. The series colors are the
Okabe-Ito colors, which remain distinguishable with color-vision deficiency. The typeface is Inter (SIL Open
Font License, files in ``fonts/``).
"""
import base64
import io
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch

__all__ = ["LIGHT", "DARK", "WIDTH", "BASE", "LINE", "show", "panel", "arrow", "mix", "on"]

FONT = "Inter"
for _file in sorted((Path(__file__).resolve().parent / "fonts").glob("*.otf")):
    font_manager.fontManager.addfont(str(_file))

BASE = 9.5           # font size (pt) of axis labels and annotations
LINE = 2.4           # width (pt) of data lines
DPI = 200            # pixels per inch in the saved image
PX_PER_INCH = 96     # size on the page: a figure 6 in wide is shown 576 px wide
WIDTH = 7.0          # width (in) of the text column; a wider figure is shrunk to fit, and its text with it


def mix(color1, color2, t):
    """Color between ``color1`` (t = 0) and ``color2`` (t = 1)."""
    a, b = mcolors.to_rgb(color1), mcolors.to_rgb(color2)
    return mcolors.to_hex([x + (y - x) * t for x, y in zip(a, b)])


def on(color):
    """Near-black or white, whichever is easier to read on the given fill."""
    r, g, b = mcolors.to_rgb(color)
    return "#111111" if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.55 else "#ffffff"


@dataclass(frozen=True)
class Theme:
    mode: str
    ground: str
    ink: str
    muted: str
    rule: str
    fill: str
    blue: str
    vermilion: str
    green: str
    gray: str
    carbon: str

    @property
    def series(self):
        return [self.blue, self.vermilion, self.green, self.gray]

    solid = property(lambda self: self.blue)
    liquid = property(lambda self: self.green)
    gas = property(lambda self: self.vermilion)
    oxygen = property(lambda self: self.vermilion)

    def tint(self, color, t=0.6):
        """Paler version of a color, for fills (it moves toward the page color)."""
        return mix(color, self.ground, t)

    def rc(self):
        return {
            "font.family": FONT, "font.size": BASE,
            "axes.labelsize": BASE, "axes.titlesize": BASE, "legend.fontsize": BASE - 1,
            "xtick.labelsize": BASE - 1, "ytick.labelsize": BASE - 1,
            "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": f"{FONT}:italic",
            "mathtext.bf": f"{FONT}:bold", "mathtext.bfit": f"{FONT}:bold:italic",
            "mathtext.cal": FONT, "mathtext.sf": FONT,
            "text.color": self.ink, "axes.labelcolor": self.ink, "axes.edgecolor": self.muted,
            "xtick.color": self.muted, "ytick.color": self.muted,
            "xtick.labelcolor": self.ink, "ytick.labelcolor": self.ink,
            "axes.linewidth": 0.9, "xtick.major.width": 0.9, "ytick.major.width": 0.9,
            "xtick.major.size": 4, "ytick.major.size": 4, "xtick.direction": "out", "ytick.direction": "out",
            "axes.spines.top": False, "axes.spines.right": False,
            "axes.grid": False, "grid.color": self.rule, "grid.linewidth": 0.8, "axes.axisbelow": True,
            "axes.prop_cycle": mpl.cycler(color=self.series),
            "lines.linewidth": LINE, "lines.solid_capstyle": "round",
            "legend.frameon": False,
            "figure.facecolor": "none", "axes.facecolor": "none", "savefig.facecolor": "none",
            "savefig.dpi": DPI,
        }


LIGHT = Theme("light", ground="#ffffff", ink="#1c1917", muted="#78716c", rule="#e7e5e4", fill="#f1f0ee",
              blue="#0072B2", vermilion="#D55E00", green="#009E73", gray="#78716c", carbon="#57534e")
DARK = Theme("dark", ground="#1c1917", ink="#f5f5f4", muted="#a8a29e", rule="#3a3531", fill="#2a2623",
             blue="#56B4E9", vermilion="#F4915A", green="#3FCF9F", gray="#a8a29e", carbon="#a8a29e")


def panel(ax, letter, title=""):
    """Panel label: a bold letter and a short title above the top-left corner of the axes."""
    ax.set_title(rf"$\mathbf{{{letter}}}$   {title}", loc="left", fontsize=BASE + 0.5, pad=9)


def arrow(ax, start, end, color, lw=1.4, scale=13):
    """Straight arrow from ``start`` to ``end`` in data coordinates."""
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=scale, lw=lw, color=color,
                                 shrinkA=0, shrinkB=0, zorder=4))


def _png(draw, theme):
    with mpl.rc_context(theme.rc()):
        fig = draw(theme)
        buffer = io.BytesIO()
        fig.savefig(buffer, format="png", dpi=DPI, transparent=True)
        width = fig.get_figwidth()
        plt.close(fig)
    return base64.b64encode(buffer.getvalue()).decode("ascii"), width


def show(draw, alt=""):
    """Draw the figure for the light and the dark page and display the pair.

    ``draw`` takes one argument, the colors of the page (``LIGHT`` or ``DARK``), and returns the figure.
    ``alt`` is a sentence that describes the figure to a reader who cannot see it.
    """
    from html import escape
    from IPython.display import display

    (light, width), (dark, _) = _png(draw, LIGHT), _png(draw, DARK)
    image = ('<div class="{cls}"><img alt="{alt}" src="data:image/png;base64,{data}" '
             'style="width:{px}px;max-width:100%;height:auto;margin:0"></div>')
    px = round(width * PX_PER_INCH)
    html = (image.format(cls="dark:hidden", alt=escape(alt, quote=True), data=light, px=px)
            + image.format(cls="hidden dark:block", alt=escape(alt, quote=True), data=dark, px=px))
    display({"text/html": html, "image/png": light}, raw=True)
