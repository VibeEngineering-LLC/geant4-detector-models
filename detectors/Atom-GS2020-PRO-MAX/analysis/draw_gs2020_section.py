import sys
sys.stdout.reconfigure(encoding='utf-8')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
import matplotlib.patches as mpatches

fig, ax = plt.subplots(figsize=(11, 8))

# Data definition
# (name, color, points)
polygons = [
    ('Al крышка', '#b0b7c3', [(0,37.75),(29.5,37.75),(29.5,10.25),(29.1,10.25),(29.1,18.75),(28.5,18.75),(28.5,36.75),(0,36.75)]),
    ('Al стенка', '#b0b7c3', [(28.5,18.75),(28.9,18.75),(28.9,8.75),(29.5,8.75),(29.5,-18.75),(28.7,-18.75),(28.7,-17.25),(26.75,-17.25),(26.75,-15.25),(28.5,-15.25)]),
    ('Отражатель', '#f1e6a8', [(0,36.75),(28.5,36.75),(28.5,-15.25),(25.4,-15.25),(25.4,33.65),(0,33.65)]),
    ('NaI', '#8fd3f4', [(0,33.55),(25.4,33.55),(25.4,-17.25),(0,-17.25)]),
    ('Окно', '#c9a6e8', [(0,-17.25),(28.5,-17.25),(28.5,-20.25),(0,-20.25)]),
    ('Эпоксид', '#f4a261', [(28.9,18.75),(29.1,18.75),(29.1,10.25),(29.5,10.25),(29.5,8.75),(28.9,8.75)]),
    ('Эпоксид дно', '#f4a261', [(28.5,-17.25),(28.7,-17.25),(28.7,-18.75),(29.5,-18.75),(29.39,-19.2),(29.18,-19.62),(28.88,-19.97),(28.5,-20.25)])
]

# Draw polygons (mirrored)
for name, color, points in polygons:
    # Right side (r >= 0)
    poly_right = Polygon(points, closed=True, facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.9)
    ax.add_patch(poly_right)
    
    # Left side (r <= 0)
    points_left = [(-r, z) for r, z in points]
    poly_left = Polygon(points_left, closed=True, facecolor=color, edgecolor='black', linewidth=0.5, alpha=0.9)
    ax.add_patch(poly_left)

# SiPM box
rect_sipm = Rectangle((-28.5, -25), 57, 4.7, facecolor='none', edgecolor='green', linestyle='--', linewidth=1.5)
ax.add_patch(rect_sipm)
ax.text(0, -22.65, 'SiPM + плата FR4 + корпус (STL, схематично)', ha='center', va='center', fontsize=8, color='green')

# Dimensions
# Vertical 1: Total height
ax.annotate('', xy=(-34, 37.75), xytext=(-34, -20.25),
            arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
ax.text(-34, 8.75, '58,0 мм', ha='right', va='center', fontsize=9, color='red', rotation=90)

# Vertical 2: Crystal height
ax.annotate('', xy=(-31.5, 33.55), xytext=(-31.5, -17.25),
            arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
ax.text(-31.5, 8.15, '50,8 мм (кристалл)', ha='right', va='center', fontsize=9, color='red', rotation=90)

# Horizontal 1: Outer diameter
ax.annotate('', xy=(29.5, 40), xytext=(-29.5, 40),
            arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
ax.text(0, 40, 'Ø59,0 мм', ha='center', va='bottom', fontsize=9, color='red')

# Horizontal 2: Crystal diameter
ax.annotate('', xy=(25.4, 23), xytext=(-25.4, 23),
            arrowprops=dict(arrowstyle='<->', color='red', lw=1.5))
ax.text(0, 23, 'Ø50,8 мм', ha='center', va='bottom', fontsize=9, color='red')

# Callouts
callouts = [
    ('Al крышка 1,0 мм; G4_Al', (14, 37.2), 36),
    ('Отражатель MgO 3,1 мм; 1,8 г/см³', (27, 28), 28),
    ('NaI Ø50,8×50,8; 3,667 г/см³', (10, 10), 10),
    ('Эпоксидный шов 0,2–0,4 мм; 1,2 г/см³', (29.3, 14), 17),
    ('Al стенка 0,8–1,0 мм', (29.2, -8), -4),
    ('Окно SiO₂ Ø57×3; 2,2 г/см³', (20, -18.75), -14),
    ('Эпоксидный шов дна', (29.2, -19.2), -20)
]

for text, point, text_y in callouts:
    ax.annotate(text, xy=point, xytext=(36, text_y),
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.8),
                fontsize=9, ha='left', va='center')

# Axis settings
ax.set_xlim(-40, 80)
ax.set_ylim(-27, 44)
ax.set_aspect('equal')
ax.set_xlabel('r, мм', fontsize=12)
ax.set_ylabel('z, мм', fontsize=12)
ax.set_title('GS2020 + Capescint NaI 51×51: осевой разрез (GDML)', fontsize=14)
ax.grid(True, linestyle=':', alpha=0.5)

# Save
plt.savefig('<WORKDIR>/GEANT4/gs2020_section_dims.png', dpi=130, bbox_inches='tight')
