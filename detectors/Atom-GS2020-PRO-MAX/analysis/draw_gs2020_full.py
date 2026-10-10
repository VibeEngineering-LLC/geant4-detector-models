import sys
sys.stdout.reconfigure(encoding='utf-8')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Polygon

fig, ax = plt.subplots(figsize=(12, 11))
ax.set_aspect('equal')
ax.set_xlim(-47, 80)
ax.set_ylim(-47, 47)
ax.set_xlabel('x, мм')
ax.set_ylabel('z, мм')
ax.set_title('GS2020 в сборке (GDML, осевой разрез)')
ax.grid(True, linestyle=':', alpha=0.3)

# LAYER A
layer_a = [
    ('Плата FR4 №2', -30.5, 29.0, -33.9, -26.8, '#7fb77e', 0.35),
    ('Плата FR4 №1 (SiPM)', -29.3, 29.2, -23.2, -21.6, '#5fa35e', 0.35),
    ('Пена PE', -29.7, 29.8, 37.8, 39.8, '#fff2b3', 0.35),
    ('SiPM Si (корпуса)', -23.3, 23.2, -21.5, -20.2, '#d9534f', 0.8)
]
for name, x0, x1, z0, z1, color, alpha in layer_a:
    ax.add_patch(Rectangle((x0, z0), x1 - x0, z1 - z0, facecolor=color, alpha=alpha, edgecolor='gray', linestyle='--', zorder=1))

# LAYER B
def draw_half_profile(pts, color, zorder):
    x = [p[0] for p in pts]
    z = [p[1] for p in pts]
    ax.add_patch(Polygon(list(zip(x, z)), closed=True, facecolor=color, edgecolor='black', linewidth=0.6, zorder=zorder))
    ax.add_patch(Polygon(list(zip([-v for v in x], z)), closed=True, facecolor=color, edgecolor='black', linewidth=0.6, zorder=zorder))

# 0 Shell
shell_pts = [(29.75,-36.5),(31.5,-36.5),(31.5,41.5),(0,41.5),(0,39.75),(29.75,39.75)]
draw_half_profile(shell_pts, '#9aa5b1', 2)

# 0b Bottom cover
cover_pts = [(3.18,-41.5),(31.0,-41.5),(31.5,-41.0),(31.5,-36.5),(29.75,-36.5),(29.75,-30.5),(28.25,-30.5),(28.25,-32.99),(27.25,-32.99),(27.25,-36.5),(5.0,-36.5),(5.0,-31.4),(0,-31.4),(0,-32.5),(3.18,-32.5)]
draw_half_profile(cover_pts, '#7d8a99', 3)

# 0c Ring
ax.add_patch(Rectangle((22.43, -28.19), 29.82 - 22.43, -23.19 - (-28.19), facecolor='#b0b7c3', edgecolor='black', linewidth=0.6, zorder=4))
ax.add_patch(Rectangle((-29.82, -28.19), 29.82 - 22.43, -23.19 - (-28.19), facecolor='#b0b7c3', edgecolor='black', linewidth=0.6, zorder=4))

# 1 Al cap
cap_pts = [(0,37.75),(29.5,37.75),(29.5,10.25),(29.1,10.25),(29.1,18.75),(28.5,18.75),(28.5,36.75),(0,36.75)]
draw_half_profile(cap_pts, '#b0b7c3', 5)

# 2 Al wall
wall_pts = [(28.5,18.75),(28.9,18.75),(28.9,8.75),(29.5,8.75),(29.5,-18.75),(28.7,-18.75),(28.7,-17.25),(26.75,-17.25),(26.75,-15.25),(28.5,-15.25)]
draw_half_profile(wall_pts, '#b0b7c3', 6)

# 3 Reflector
refl_pts = [(0,36.75),(28.5,36.75),(28.5,-15.25),(25.4,-15.25),(25.4,33.65),(0,33.65)]
draw_half_profile(refl_pts, '#f1e6a8', 7)

# 4 NaI
nai_pts = [(0,33.55),(25.4,33.55),(25.4,-17.25),(0,-17.25)]
draw_half_profile(nai_pts, '#8fd3f4', 8)

# 5 Window
win_pts = [(0,-17.25),(28.5,-17.25),(28.5,-20.25),(0,-20.25)]
draw_half_profile(win_pts, '#c9a6e8', 9)

# 6 Epoxy
epoxy_pts = [(28.9,18.75),(29.1,18.75),(29.1,10.25),(29.5,10.25),(29.5,8.75),(28.9,8.75)]
draw_half_profile(epoxy_pts, '#f4a261', 10)

# 7 Epoxy bottom
epoxy_b_pts = [(28.5,-17.25),(28.7,-17.25),(28.7,-18.75),(29.5,-18.75),(29.39,-19.2),(29.18,-19.62),(28.88,-19.97),(28.5,-20.25)]
draw_half_profile(epoxy_b_pts, '#f4a261', 11)

# SCREWS
ax.add_patch(Rectangle((26.5, -27.9), 31.5 - 26.5, -23.2 - (-27.9), facecolor='#444', edgecolor='black', linewidth=0.6, zorder=10))
ax.add_patch(Rectangle((26.5, -35.4), 31.5 - 26.5, -30.6 - (-35.4), facecolor='#444', edgecolor='black', linewidth=0.6, zorder=10))
ax.add_patch(Rectangle((-17.8, -27.9), -12.3 - (-17.8), -23.2 - (-27.9), facecolor='none', edgecolor='#444', linestyle='--', linewidth=0.8, zorder=10))
ax.add_patch(Rectangle((-17.8, -35.4), -12.3 - (-17.8), -30.6 - (-35.4), facecolor='none', edgecolor='#444', linestyle='--', linewidth=0.8, zorder=10))

# CALLOUTS
callouts = [
    ('Пена PE 2 мм; 0,03 г/см³', (20, 38.8)),
    ('Отражатель MgO 3,1 мм; 1,8 г/см³', (27, 25)),
    ('Кожух Shell, Al: стенка и крышка 1,75 мм; Ø63×78', (31, 20)),
    ('NaI Ø50,8×50,8; 3,667 г/см³', (10, 10)),
    ('Al корпус кристалла; 0,8–1,0 мм', (29.2, 0)),
    ('Окно SiO₂ Ø57×3; 2,2 г/см³', (15, -18.7)),
    ('SiPM (Si) 1,3 мм + выводы Cu', (15, -20.8)),
    ('Плата FR4 №1, 1,6 мм; 1,86 г/см³', (15, -22.4)),
    ('Нажимное кольцо Al: 5 мм по z, r 22,4…29,8 мм', (23.5, -27.5)),
    ('Винты нерж. сталь: радиальные, 3 поз. по 120°, Ø4,7×5 мм', (29, -25.5)),
    ('Плата FR4 №2, 7,1 мм', (25, -30)),
    ('Нижняя крышка, Al: плита 5 мм, борт 2,5 мм', (28, -38))
]
for i, (text, target) in enumerate(callouts):
    y = 44 - i * 88 / 11
    ax.annotate(text, xy=target, xytext=(40, y), fontsize=9,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

# DIMENSIONS
# Vertical 83
ax.annotate('', xy=(-41, 41.5), xytext=(-41, -41.5), arrowprops=dict(arrowstyle='<->', color='red'))
ax.text(-42, 0, 'полная высота 83 мм', rotation=90, ha='right', va='center', color='red')

# Vertical 58
ax.annotate('', xy=(-35, 37.75), xytext=(-35, -20.25), arrowprops=dict(arrowstyle='<->', color='red'))
ax.text(-36, 8.75, 'кристалл+корпус 58 мм', rotation=90, ha='right', va='center', color='red')

# Horizontal 63
ax.annotate('', xy=(31.5, -44), xytext=(-31.5, -44), arrowprops=dict(arrowstyle='<->', color='red'))
ax.text(0, -44, 'Ø63 мм (кожух)', ha='center', va='top', color='red')

# Horizontal 59
ax.annotate('', xy=(29.5, 43.5), xytext=(-29.5, 43.5), arrowprops=dict(arrowstyle='<->', color='red'))
ax.text(0, 43.5, 'Ø59 мм (корпус кристалла)', ha='center', va='bottom', color='red')

# FOOTNOTE
fig.text(0.5, 0.01, 'Кожух, оптический блок и нижняя крышка — по профилю GDML/STL; платы, кольцо, SiPM и винты — габаритные прямоугольники; отверстие в крышке глухое (дно z −32,5), верх бобышки по грани z −31,4 условно.', ha='center', fontsize=8)

plt.savefig('<WORKDIR>/GEANT4/gs2020_full_assembly.png', dpi=130, bbox_inches='tight')
