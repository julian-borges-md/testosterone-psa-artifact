"""Figure style conforming to the Endocrine Society Digital Art Guidelines.

The guidelines this encodes, so a reader can check the code against them:

  Width      one column 3.25 in, two column 5.8 in. Submit at final publication
             size and do not scale. Height at most 9.25 in.
  Font       a standard font: Arial, Helvetica, Times, Symbol, Mathematical Pi
             or European Pi. Same font for every figure in the manuscript.
             Preferred size 8 point. Fonts embedded in the file.
  Sizes      do not vary letter size within a single figure.
  Panels     lowercase bold with parentheses, 12 point: (a), (b), (c).
  Lines      at least 0.5 to 1 point when the figure is at its column width.
  Colour     submit in original RGB.

Helvetica and Arial are not installed on every machine that builds these
figures. The family list below prefers them and falls back to their metric
compatible clones, Nimbus Sans for Helvetica and Liberation Sans for Arial, so
the rendered glyph metrics match whichever machine runs the build. Font type 42
embeds TrueType outlines in the PDF and EPS rather than matplotlib's default
Type 3, which some production systems reject.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ONE_COLUMN = 3.25
TWO_COLUMN = 5.8
MAX_HEIGHT = 9.25

BASE_PT = 8
PANEL_LABEL_PT = 12

# Fixed per exposure across every figure.
C = {'biot': '#1f4e79', 'tt': '#7f7f7f', 'shbg': '#c55a11',
     'bmi': '#548235', 'igf': '#7030a0'}

# Okabe and Ito, chosen so the two competing paths in the causal diagram stay
# distinguishable under deuteranopia and protanopia. A red against green pair
# does not, and this is a paper about men's health.
PATH_DETECTION = '#D55E00'
PATH_BIOLOGY = '#009E73'


def apply():
    plt.rcParams.update({
        'font.family': 'sans-serif',
        'font.sans-serif': ['Helvetica', 'Arial', 'Nimbus Sans',
                            'Liberation Sans', 'DejaVu Sans'],
        'font.size': BASE_PT,
        'axes.titlesize': BASE_PT,
        'axes.labelsize': BASE_PT,
        'xtick.labelsize': BASE_PT,
        'ytick.labelsize': BASE_PT,
        'legend.fontsize': BASE_PT,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'savefig.dpi': 600,
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
    })


def tint(hex_colour, fraction, over='#ffffff'):
    """The colour an alpha blend would have produced, as a solid value.

    EPS has no transparency: the PostScript backend renders a partially
    transparent artist opaque, so a band drawn with alpha becomes a solid block
    in the submitted vector file. Computing the blend here keeps the PNG, the
    PDF and the EPS identical.
    """
    def rgb(h):
        h = h.lstrip('#')
        return [int(h[i:i + 2], 16) for i in (0, 2, 4)]
    a, b = rgb(hex_colour), rgb(over)
    return '#%02x%02x%02x' % tuple(
        int(round(fraction * x + (1 - fraction) * y)) for x, y in zip(a, b))


def panel_label(ax, letter, x=-0.02, y=1.06, ha='right'):
    """Lowercase bold with parentheses, as the guidelines require.

    Placed clear of the panel title. A two line title needs a higher y and left
    alignment, otherwise the label and the title overprint.
    """
    ax.text(x, y, '(%s)' % letter, transform=ax.transAxes,
            fontsize=PANEL_LABEL_PT, fontweight='bold',
            ha=ha, va='bottom')


def check(fig, name):
    """Refuse to write a figure that breaches the published limits."""
    w, h = fig.get_size_inches()
    if w > TWO_COLUMN + 1e-6:
        raise SystemExit('%s is %.2f in wide, over the %.2f in two column '
                         'maximum' % (name, w, TWO_COLUMN))
    if h > MAX_HEIGHT:
        raise SystemExit('%s is %.2f in tall, over the %.2f in maximum'
                         % (name, h, MAX_HEIGHT))
    return fig
