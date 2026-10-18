# Réseau autoroutier de moindre coût pour le Cameroun.
#
# Question : si un planificateur devait relier les 23 plus grandes villes du Cameroun au plus
# faible coût de construction, quel réseau construirait-il ? Ressemble-t-il aux routes actuelles ?
# Quelles villes gagneraient ou perdraient en accès au marché ?
#
# Méthode inspirée du "réseau couvrant de moindre coût" de Faber (2014) :
#   1. carte des coûts de construction (pente du terrain, rivières, frontière)
#   2. chemin le moins cher entre chaque paire de villes (Dijkstra)
#   3. arbre couvrant de coût minimal sur ces coûts (Kruskal) = réseau de moindre coût
#   4. comparaison avec les routes actuelles et calcul de l'accès au marché
#
# Lancement : python run_analysis.py   (environ 2 minutes)

import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import numpy as np
import pandas as pd
import xarray as xr
from rasterio import features
from rasterio.transform import from_origin
from scipy import ndimage
from shapely.geometry import LineString

from src.graph_algorithms import dijkstra_grid, kruskal, trace_path

ROOT = Path(__file__).parent
plt.rcParams.update({"figure.dpi": 130})

# ------------------------------------------------------------------ paramètres
UPSAMPLE = 4            # grille du relief affinée de 10' à 2,5' (cases d'environ 4,6 km)
SLOPE_WEIGHT = 0.5      # coût supplémentaire par % de pente
WATER_COST = 25.0       # traverser une rivière ou un lac coûte 25 fois plus (ponts)
OFFROAD_RATIO = 5       # se déplacer hors route coûte 5 fois plus que sur une route
THETA = 1.0             # élasticité à la distance dans l'accès au marché

# ------------------------------------------------------------------ données
border = gpd.read_file(ROOT / "data/cameroon_boundary.geojson")
cities = gpd.read_file(ROOT / "data/cities.geojson").reset_index(drop=True)
roads = gpd.read_file(ROOT / "data/existing_roads_natural_earth.geojson")
rivers = gpd.read_file(ROOT / "data/rivers.geojson")
lakes = gpd.read_file(ROOT / "data/lakes.geojson")
topo = xr.open_dataset(ROOT / "data/topography_10arcmin_cameroon.nc")["topography"]
n = len(cities)

# relief : nord en haut, puis interpolation bilinéaire sur une grille 4 fois plus fine
z = ndimage.zoom(topo.values[::-1], UPSAMPLE, order=1)
n_rows, n_cols = z.shape
lat_top = float(topo.latitude.max())
lon_left = float(topo.longitude.min())
res = float(topo.longitude.diff("longitude")[0]) / UPSAMPLE          # taille d'une case en degrés
transform = from_origin(lon_left - res / 2, lat_top + res / 2, res, res)
cell_km_y = 111.32 * res                                              # 1 degré de latitude = 111,32 km
cell_km_x = 111.32 * res * math.cos(math.radians(6.5))                # latitude moyenne du pays : 6,5°

# ------------------------------------------------------------------ 1. CARTE DES COÛTS
# pente (%) = 100 x racine(pente nord-sud² + pente est-ouest²)
# pente d'une case = (altitude du voisin d'après - altitude du voisin d'avant) / (2 x taille de case)
# sur les bords, on utilise la case elle-même au lieu du voisin manquant.
def slope_percent(elev, dy_m, dx_m):
    rows, cols = elev.shape
    slope = np.zeros((rows, cols))
    for r in range(rows):
        for c in range(cols):
            if r == 0:
                gy = (elev[1, c] - elev[0, c]) / dy_m
            elif r == rows - 1:
                gy = (elev[r, c] - elev[r - 1, c]) / dy_m
            else:
                gy = (elev[r + 1, c] - elev[r - 1, c]) / (2 * dy_m)
            if c == 0:
                gx = (elev[r, 1] - elev[r, 0]) / dx_m
            elif c == cols - 1:
                gx = (elev[r, c] - elev[r, c - 1]) / dx_m
            else:
                gx = (elev[r, c + 1] - elev[r, c - 1]) / (2 * dx_m)
            slope[r, c] = 100 * math.sqrt(gx * gx + gy * gy)
    return slope


slope_pct = slope_percent(z, cell_km_y * 1000, cell_km_x * 1000)

# cases dans le pays, et cases d'eau (dessin des formes géographiques sur la grille)
border_shapes = []
for g in border.geometry:
    border_shapes.append((g, 1))
water_shapes = []
for g in list(rivers.geometry) + list(lakes.geometry):
    water_shapes.append((g, 1))
inside = features.rasterize(border_shapes, out_shape=z.shape, transform=transform, fill=0,
                            all_touched=False).astype(bool)
water = features.rasterize(water_shapes, out_shape=z.shape, transform=transform, fill=0,
                           all_touched=True).astype(bool)


# coût par km = 1 + 0,5 x pente + 25 si eau ; infini hors du pays (réseau national)
def construction_cost_grid(slope_weight, water_cost):
    grid = []
    for r in range(n_rows):
        line = []
        for c in range(n_cols):
            if not inside[r, c]:
                line.append(math.inf)
            else:
                value = 1 + slope_weight * slope_pct[r, c]
                if water[r, c]:
                    value = value + water_cost
                line.append(value)
        grid.append(line)
    return grid


cost = construction_cost_grid(SLOPE_WEIGHT, WATER_COST)


def to_cell(lon, lat):
    return int(round((lat_top - lat) / res)), int(round((lon - lon_left) / res))


def cell_to_lonlat(cell):
    return lon_left + cell[1] * res, lat_top - cell[0] * res


# case de chaque ville ; une ville tombée sur l'eau ou hors du pays est déplacée vers la case
# valide la plus proche
city_cell = []
for i in range(n):
    r, c = to_cell(cities.geometry[i].x, cities.geometry[i].y)
    if cost[r][c] == math.inf or water[r, c]:
        best, best_d = None, math.inf
        for rr in range(n_rows):
            for cc in range(n_cols):
                if cost[rr][cc] != math.inf and not water[rr, cc]:
                    d = (rr - r) ** 2 + (cc - c) ** 2
                    if d < best_d:
                        best, best_d = (rr, cc), d
        r, c = best
    city_cell.append((r, c))

# ------------------------------------------------------------------ 2. CHEMINS LES MOINS CHERS
# POUR chaque ville i : Dijkstra depuis i, puis lire le coût et le tracé vers chaque ville j > i
pair_cost = {}
pair_path = {}
for i in range(n):
    dist, prev = dijkstra_grid(cost, city_cell[i], cell_km_y, cell_km_x)
    for j in range(i + 1, n):
        rj, cj = city_cell[j]
        pair_cost[(i, j)] = dist[rj][cj]
        pair_path[(i, j)] = trace_path(prev, city_cell[j])

# ------------------------------------------------------------------ 3. ARBRES COUVRANTS
edges = []
for (i, j) in pair_cost:
    edges.append((pair_cost[(i, j)], i, j))
lcn_links = kruskal(n, edges)                       # réseau de moindre coût

geo_m = cities.to_crs(32633)                        # projection UTM 33N, en mètres
straight_edges = []
for i in range(n):
    for j in range(i + 1, n):
        straight_edges.append((geo_m.geometry[i].distance(geo_m.geometry[j]) / 1000, i, j))
straight_links = kruskal(n, straight_edges)         # arbre en lignes droites (référence naïve)

lcn_rows = []
for (i, j) in lcn_links:
    points = []
    for cell in pair_path[(i, j)]:
        points.append(cell_to_lonlat(cell))
    lcn_rows.append({"from": cities["name"][i], "to": cities["name"][j],
                     "construction_cost": pair_cost[(i, j)], "geometry": LineString(points)})
lcn_lines = gpd.GeoDataFrame(lcn_rows, crs=4326)
lcn_lines["length_km"] = lcn_lines.to_crs(32633).length / 1000

straight_rows = []
for (i, j) in straight_links:
    straight_rows.append({"from": cities["name"][i], "to": cities["name"][j],
                          "geometry": LineString([cities.geometry[i], cities.geometry[j]])})
straight_lines = gpd.GeoDataFrame(straight_rows, crs=4326)
straight_lines["length_km"] = straight_lines.to_crs(32633).length / 1000

lcn_lines.to_file(ROOT / "results/least_cost_network.geojson", driver="GeoJSON")
lcn_lines.drop(columns="geometry").sort_values("construction_cost").to_csv(
    ROOT / "results/least_cost_network_links.csv", index=False, float_format="%.1f")


# coût de construction d'un tracé : on avance par demi-case et on additionne le coût local
def construction_cost_of(lines):
    total = 0.0
    for geom in lines.geometry:
        line_m = gpd.GeoSeries([geom], crs=4326).to_crs(32633).iloc[0]
        length_km = line_m.length / 1000
        steps = max(int(length_km / (cell_km_x / 2)), 1)
        points_m = []
        for k in range(steps + 1):
            points_m.append(line_m.interpolate(k / steps, normalized=True))
        points = gpd.GeoSeries(points_m, crs=32633).to_crs(4326)
        values = []
        for p in points:
            r, c = to_cell(p.x, p.y)
            r = min(max(r, 0), n_rows - 1)
            c = min(max(c, 0), n_cols - 1)
            if cost[r][c] == math.inf:
                values.append(WATER_COST)
            else:
                values.append(cost[r][c])
        total = total + sum(values) / len(values) * length_km
    return total


roads_in = gpd.clip(roads, border)
lcn_length = lcn_lines["length_km"].sum()
straight_length = straight_lines["length_km"].sum()
roads_length = roads_in.to_crs(32633).length.sum() / 1000
lcn_build = construction_cost_of(lcn_lines)
straight_build = construction_cost_of(straight_lines)

# part du réseau située à moins de 10 km d'une route existante
road_buffer = roads_in.to_crs(32633).buffer(10_000).union_all()
lcn_m = lcn_lines.to_crs(32633)
straight_m = straight_lines.to_crs(32633)
share_lcn = lcn_m.intersection(road_buffer).length.sum() / lcn_m.length.sum()
share_straight = straight_m.intersection(road_buffer).length.sum() / straight_m.length.sum()

comparison = [
    {"Network": "Least-cost network (MST on least-cost paths)", "Links": len(lcn_lines), "Length (km)": lcn_length,
     "Construction cost (index)": lcn_build, "Cost per km": lcn_build / lcn_length,
     "Share within 10 km of existing road": share_lcn},
    {"Network": "Straight-line MST", "Links": len(straight_lines), "Length (km)": straight_length,
     "Construction cost (index)": straight_build, "Cost per km": straight_build / straight_length,
     "Share within 10 km of existing road": share_straight},
    {"Network": "Existing trunk roads (Natural Earth)", "Links": len(roads_in), "Length (km)": roads_length},
]
pd.DataFrame(comparison).to_csv(ROOT / "results/network_comparison.csv", index=False, float_format="%.3f")

# ------------------------------------------------------------------ 4. ACCÈS AU MARCHÉ
# ALGORITHME AccesAuMarche
#   coût de déplacement par km : 1 sur une case de route, 5 x (1 + 0,1 x pente) ailleurs
#   POUR chaque ville i FAIRE
#       T(i, j) <- coût du trajet le moins cher de i à j (Dijkstra)
#       MA(i) <- somme sur j différent de i de Population(j) / T(i, j)
def travel_cost_grid(road_lines):
    if len(road_lines) > 0:
        shapes = []
        for g in road_lines.geometry:
            shapes.append((g, 1))
        on_road = features.rasterize(shapes, out_shape=z.shape, transform=transform, fill=0,
                                     all_touched=True).astype(bool)
    else:
        on_road = np.zeros(z.shape, bool)
    grid = []
    for r in range(n_rows):
        line = []
        for c in range(n_cols):
            if not inside[r, c]:
                line.append(math.inf)
            elif on_road[r, c]:
                line.append(1.0)
            else:
                line.append(OFFROAD_RATIO * (1 + 0.1 * slope_pct[r, c]))
        grid.append(line)
    return grid


def market_access(grid):
    population = list(cities["pop_max"].astype(float))
    access = []
    for i in range(n):
        dist, _ = dijkstra_grid(grid, city_cell[i], cell_km_y, cell_km_x)
        total = 0.0
        for j in range(n):
            if j != i:
                rj, cj = city_cell[j]
                total = total + population[j] * dist[rj][cj] ** (-THETA)
        access.append(total)
    return access


scenarios = {"No network": gpd.GeoDataFrame(geometry=[], crs=4326), "Existing roads": roads_in,
             "Least-cost network": lcn_lines, "Straight-line MST": straight_lines}
access = {}
for name in scenarios:
    access[name] = market_access(travel_cost_grid(scenarios[name]))

ma_rows = []
for i in range(n):
    none, existing = access["No network"][i], access["Existing roads"][i]
    lcn, straight = access["Least-cost network"][i], access["Straight-line MST"][i]
    ma_rows.append({"name": cities["name"][i], "population": cities["pop_max"][i], "region": cities["adm1name"][i],
                    "No network": none, "Existing roads": existing, "Least-cost network": lcn,
                    "Straight-line MST": straight,
                    "Gain from existing roads (%)": 100 * (existing / none - 1),
                    "Gain from least-cost network (%)": 100 * (lcn / none - 1),
                    "Least-cost vs existing (%)": 100 * (lcn / existing - 1)})
pd.DataFrame(ma_rows).to_csv(ROOT / "results/market_access.csv", index=False, float_format="%.2f")

# moyenne pondérée par la population
weighted = {}
for name in scenarios:
    num, den = 0.0, 0.0
    for i in range(n):
        num = num + access[name][i] * cities["pop_max"][i]
        den = den + cities["pop_max"][i]
    weighted[name] = num / den
lengths = {"No network": 0, "Existing roads": roads_length, "Least-cost network": lcn_length,
           "Straight-line MST": straight_length}
summary_rows = []
for name in scenarios:
    gain = 100 * (weighted[name] / weighted["No network"] - 1)
    if lengths[name] > 0:
        per_1000 = gain / lengths[name] * 1000
    else:
        per_1000 = float("nan")
    summary_rows.append({"Scenario": name, "Population-weighted market access": weighted[name],
                         "Gain vs no network (%)": gain, "Network length (km)": lengths[name],
                         "MA gain per 1,000 km (%)": per_1000})
pd.DataFrame(summary_rows).to_csv(ROOT / "results/market_access_summary.csv", index=False, float_format="%.3f")

# ------------------------------------------------------------------ 5. SENSIBILITÉ AUX PARAMÈTRES
# Le réseau change-t-il si on modifie le poids de la pente ou le coût de l'eau ?
base_links = []
for link in lcn_links:
    base_links.append(tuple(sorted(link)))
sensitivity = []
for slope_w, water_c in [(0.0, 25), (0.25, 25), (1.0, 25), (0.5, 5), (0.5, 50)]:
    grid = construction_cost_grid(slope_w, water_c)
    test_edges = []
    for i in range(n):
        dist, _ = dijkstra_grid(grid, city_cell[i], cell_km_y, cell_km_x)
        for j in range(i + 1, n):
            rj, cj = city_cell[j]
            test_edges.append((dist[rj][cj], i, j))
    shared = 0
    for link in kruskal(n, test_edges):
        if tuple(sorted(link)) in base_links:
            shared = shared + 1
    sensitivity.append({"Slope weight": slope_w, "Water cost": water_c, "Links shared with baseline": shared,
                        "Total links": len(base_links)})
pd.DataFrame(sensitivity).to_csv(ROOT / "results/sensitivity.csv", index=False)

# ------------------------------------------------------------------ 6. CARTES
extent = [lon_left - res / 2, lon_left + res * (n_cols - 0.5), lat_top - res * (n_rows - 0.5), lat_top + res / 2]
relief = np.where(inside, z, np.nan)
fig, axes = plt.subplots(1, 2, figsize=(12, 8.5))
titles = ["Least-cost network vs existing trunk roads", "Straight-line MST (naive benchmark)"]
for k in range(2):
    ax = axes[k]
    ax.imshow(relief, extent=extent, cmap="terrain", vmin=-550, vmax=2100, alpha=0.75)
    border.boundary.plot(ax=ax, color="black", lw=0.8)
    rivers.plot(ax=ax, color="#3b8bc2", lw=0.6, alpha=0.7)
    ax.set_title(titles[k], fontsize=11)
roads_in.plot(ax=axes[0], color="#555", lw=1.2, ls="--", label="Existing trunk roads")
lcn_lines.plot(ax=axes[0], color="#c0392b", lw=2.2, label="Least-cost network")
straight_lines.plot(ax=axes[1], color="#8e44ad", lw=2.2, label="Straight-line MST")
for ax in axes:
    ax.scatter(cities.geometry.x, cities.geometry.y, s=np.sqrt(cities.pop_max) / 6, color="black",
               zorder=5, edgecolor="white", lw=0.5)
    for i in range(n):
        if cities["pop_max"][i] > 150_000 or cities["featurecla"][i] != "Populated place":
            ax.annotate(cities["name"][i], (cities.geometry[i].x, cities.geometry[i].y), fontsize=7,
                        xytext=(3, 3), textcoords="offset points")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_xlim(8.3, 16.4)
    ax.set_ylim(1.4, 13.3)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
fig.suptitle("A least-cost highway network for Cameroon's 23 largest cities", fontsize=13)
fig.tight_layout()
fig.savefig(ROOT / "figures/network_maps.png")
plt.close(fig)

log_cost = np.full(z.shape, np.nan)
for r in range(n_rows):
    for c in range(n_cols):
        if cost[r][c] != math.inf:
            log_cost[r, c] = math.log(1 + cost[r][c])
fig, ax = plt.subplots(figsize=(7, 8.5))
im = ax.imshow(log_cost, extent=extent, cmap="magma_r")
border.boundary.plot(ax=ax, color="black", lw=0.8)
lcn_lines.plot(ax=ax, color="#00b3b3", lw=2)
ax.scatter(cities.geometry.x, cities.geometry.y, s=12, color="black", zorder=5)
ax.set_xlim(8.3, 16.4)
ax.set_ylim(1.4, 13.3)
ax.set_xticks([])
ax.set_yticks([])
ax.set_xlabel("")
ax.set_ylabel("")
ax.set_title("Construction-cost surface (log) and least-cost network")
plt.colorbar(im, ax=ax, shrink=0.6, label="log(1 + cost per km)")
fig.tight_layout()
fig.savefig(ROOT / "figures/cost_surface.png")
plt.close(fig)

pairs = []
for row in ma_rows:
    pairs.append((row["Least-cost vs existing (%)"], row["name"]))
pairs.sort()                                        # de la plus forte perte au plus fort gain
names, values, bar_colors = [], [], []
for value, name in pairs:
    names.append(name)
    values.append(value)
    bar_colors.append("#2ca02c" if value >= 0 else "#c0392b")
fig, ax = plt.subplots(figsize=(8, 6.5))
ax.barh(names, values, color=bar_colors)
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("Market access under the least-cost network relative to existing roads (%)")
ax.set_title("Who would gain from a least-cost network?")
fig.tight_layout()
fig.savefig(ROOT / "figures/market_access_gains.png")
plt.close(fig)

print(pd.DataFrame(comparison).round(3).to_string())
print(pd.DataFrame(summary_rows).round(3).to_string())
print(pd.DataFrame(sensitivity).to_string())
