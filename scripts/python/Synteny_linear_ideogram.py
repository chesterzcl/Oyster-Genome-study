import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mpath

# -----------------------------
# Function to draw smooth chromosome bars with fully rounded ends
# -----------------------------
def draw_chromosome_bar(ax, x0, x1, y, color, label, height=0.06):
    """
    Draw a chromosome bar with pill-shaped fully rounded ends.
    """
    width = x1 - x0
    rounding_size = height / 1.5  # This ensures full-round ends

    # Enforce pill shape with correct corner rounding
    bar = patches.FancyBboxPatch(
        (x0, y), width, height,
        boxstyle=patches.BoxStyle("Round", pad=0.0, rounding_size=rounding_size),
        facecolor=color,
        edgecolor="black",
        linewidth=1.0
    )
    ax.add_patch(bar)
    ax.text((x0 + x1) / 2, y + 0.08, label, ha='center', va='bottom', fontsize=10)

# -----------------------------
# Function to draw tapered curved chord between two aligned segments
# -----------------------------
def draw_tapered_chord(ax, top_start, top_end, y_top,
                       bottom_start, bottom_end, y_bottom,
                       bar_height=0.06,
                       taper_width=0.1, taper_height=0.2,
                       color='orchid', alpha=0.6, lw=0.5):
    """
    Draw a tapered ribbon aligned to the boundary edges of pill-shaped chromosome bars.
    """
    Path = mpath.Path

    # Adjust y to align to chromosome edge (bottom of top bar, top of bottom bar)
    y_top_edge = y_top                # bottom edge of top bar
    y_bottom_edge = y_bottom + bar_height  # top edge of bottom bar

    # Compute midpoints for curve
    mid_y = (y_top_edge + y_bottom_edge) / 2 + taper_height

    mid_x_left = (top_start + bottom_start) / 2 + taper_width * (bottom_start - top_start)
    mid_x_right = (top_end + bottom_end) / 2 + taper_width * (bottom_end - top_end)

    verts = [
        (top_start, y_top_edge),         # top-left
        (mid_x_left, mid_y),             # control point left
        (bottom_start, y_bottom_edge),   # bottom-left
        (bottom_end, y_bottom_edge),     # bottom-right
        (mid_x_right, mid_y),            # control point right
        (top_end, y_top_edge),           # top-right
        (top_start, y_top_edge)          # close path
    ]

    codes = [
        Path.MOVETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.CLOSEPOLY
    ]

    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor=color, edgecolor='black', lw=lw, alpha=alpha)
    ax.add_patch(patch)

# -----------------------------
# Example plot layout
# -----------------------------
fig, ax = plt.subplots(figsize=(10, 6))

# Y positions
y_cv = 0.85
y_mg = 0.4
bar_height = 0.06

# Chromosome coordinates
cv_start, cv_end = 1.0, 7.0
mg_start, mg_end = 1.5, 7.5

# Draw chromosome bars
draw_chromosome_bar(ax, cv_start, cv_end, y_cv, 'gray', "C. virginica", height=bar_height)
draw_chromosome_bar(ax, mg_start, mg_end, y_mg, 'skyblue', "M. gigas", height=bar_height)

# Draw properly aligned tapered synteny ribbons
draw_tapered_chord(ax, 1.2, 2.8, y_cv, 2.0, 3.5, y_mg,
                   taper_width=1, taper_height=0.2, color='violet')
draw_tapered_chord(ax, 2.9, 3.9, y_cv, 3.6, 4.4, y_mg,
                   taper_width=1, taper_height=0.2, color='violet')
draw_tapered_chord(ax, 4.0, 5.8, y_cv, 4.5, 6.2, y_mg,
                   taper_width=1, taper_height=0.2, color='mediumpurple')

# Final formatting
ax.set_xlim(0.5, 8.5)
ax.set_ylim(0.2, 1.05)
ax.axis('off')
plt.tight_layout()
plt.show()