"""A least-cost national highway network for Cameroon.

Following the "least-cost spanning network" approach of Faber (2014), and the use of
terrain-based construction costs common in the economics of African road networks
(e.g. Müller-Crepon and co-authors), we ask: if a planner had to connect Cameroon's 23 largest
cities at minimum construction cost, which network would she build, how does it compare with the
existing trunk roads, and which cities gain market access?

Steps
 1. Construction-cost surface from terrain slope, water bodies and the national border.
 2. Least-cost paths (Dijkstra on the raster) between every pair of cities.
 3. Minimum spanning tree (MST) on those pairwise costs = least-cost network (LCN).
    A straight-line MST is computed as a naive benchmark.
 4. Comparison with existing trunk roads (Natural Earth) and market-access gains.

Usage: python run_analysis.py
"""
from itertools import combinations
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import geopandas as gpd
import networkx as nx
import numpy as np
import pandas as pd
import xarray as xr
from rasterio import features
from rasterio.transform import from_origin
from scipy import ndimage
from shapely.geometry import LineString
from skimage.graph import MCP_Geometric

ROOT = Path(__file__).parent
plt.rcParams.update({"figure.dpi": 130})

# ------------------------------------------------------------------ parameters
UPSAMPLE = 4            # 10 arc-min DEM resampled to 2.5 arc-min (about 4.6 km) for smoother paths
SLOPE_WEIGHT = 0.5      # extra cost per % of slope
WATER_COST = 25.0       # crossing a river or lake cell costs 25 times flat land (bridges)
OFFROAD_SPEED_RATIO = 5 # moving off the road network is 5 times costlier than on a road
THETA = 1.0             # distance elasticity in market access

# ------------------------------------------------------------------ data
border = gpd.read_file(ROOT / "data/cameroon_boundary.geojson")
cities = gpd.read_file(ROOT / "data/cities.geojson").reset_index(drop=True)
roads = gpd.read_file(ROOT / "data/existing_roads_natural_earth.geojson")
rivers = gpd.read_file(ROOT / "data/rivers.geojson")
lakes = gpd.read_file(ROOT / "data/lakes.geojson")
topo = xr.open_dataset(ROOT / "data/topography_10arcmin_cameroon.nc")["topography"]

# upsample DEM (bilinear), north-up orientation
z = ndimage.zoom(topo.values[::-1], UPSAMPLE, order=1)
lat0, lat1 = float(topo.latitude.max()), float(topo.latitude.min())
lon0 = float(topo.longitude.min())
res = float(topo.longitude.diff("longitude")[0]) / UPSAMPLE
transform = from_origin(lon0 - res / 2, lat0 + res / 2, res, res)
nrows, ncols = z.shape
lats = lat0 - res * np.arange(nrows)
mean_lat = np.deg2rad(6.5)
cell_km_y = 111.32 * res
cell_km_x = 111.32 * res * np.cos(mean_lat)

# slope (%) from the DEM
gy, gx = np.gradient(z, cell_km_y * 1000, cell_km_x * 1000)
slope_pct = 100 * np.hypot(gx, gy)

inside = features.rasterize([(g, 1) for g in border.geometry], out_shape=z.shape, transform=transform,
                            fill=0, all_touched=False).astype(bool)
water = features.rasterize([(g, 1) for g in list(rivers.geometry) + list(lakes.geometry)],
                           out_shape=z.shape, transform=transform, fill=0, all_touched=True).astype(bool)
cost = 1 + SLOPE_WEIGHT * slope_pct + WATER_COST * water
cost[~inside] = np.inf                                     # a domestic network stays in the country


def to_rc(lon, lat):
    c = int(round((lon - lon0) / res))
    r = int(round((lat0 - lat) / res))
    return r, c


cities["rc"] = [to_rc(p.x, p.y) for p in cities.geometry]
for i, (r, c) in enumerate(cities["rc"]):          # snap any city falling on an invalid cell
    if not np.isfinite(cost[r, c]) or water[r, c]:
        rr, cc = np.where(np.isfinite(cost) & ~water)
        k = np.argmin((rr - r) ** 2 + (cc - c) ** 2)
        cities.at[i, "rc"] = (rr[k], cc[k])


def rc_to_lonlat(path):
    return [(lon0 + c * res, lat0 - r * res) for r, c in path]


# ------------------------------------------------------------------ least-cost paths between all pairs
n = len(cities)
pair_cost, pair_path = {}, {}
for i in range(n):
    mcp = MCP_Geometric(cost, sampling=(cell_km_y, cell_km_x))
    cum, _ = mcp.find_costs([cities.at[i, "rc"]])
    for j in range(i + 1, n):
        pair_cost[(i, j)] = cum[cities.at[j, "rc"]]
        pair_path[(i, j)] = mcp.traceback(cities.at[j, "rc"])

G_lcp = nx.Graph()
for (i, j), c in pair_cost.items():
    G_lcp.add_edge(i, j, weight=c)
mst_lcp = nx.minimum_spanning_tree(G_lcp)

geo = cities.to_crs(32633)                                  # UTM 33N, metres
G_euc = nx.Graph()
for i, j in combinations(range(n), 2):
    G_euc.add_edge(i, j, weight=geo.geometry[i].distance(geo.geometry[j]) / 1000)
mst_euc = nx.minimum_spanning_tree(G_euc)

# geometries of both networks
lcn_lines = gpd.GeoDataFrame(
    [{"from": cities.at[i, "name"], "to": cities.at[j, "name"],
      "construction_cost": pair_cost[tuple(sorted((i, j)))],
      "geometry": LineString(rc_to_lonlat(pair_path[tuple(sorted((i, j)))]))}
     for i, j in mst_lcp.edges()], crs=4326)
lcn_lines["length_km"] = lcn_lines.to_crs(32633).length / 1000
euc_lines = gpd.GeoDataFrame(
    [{"from": cities.at[i, "name"], "to": cities.at[j, "name"],
      "geometry": LineString([cities.geometry[i], cities.geometry[j]])} for i, j in mst_euc.edges()], crs=4326)
euc_lines["length_km"] = euc_lines.to_crs(32633).length / 1000
lcn_lines.to_file(ROOT / "results/least_cost_network.geojson", driver="GeoJSON")
lcn_lines.drop(columns="geometry").sort_values("construction_cost").to_csv(
    ROOT / "results/least_cost_network_links.csv", index=False, float_format="%.1f")


def construction_cost_of(lines):
    """Integrate the cost surface along each line (sampled every half cell)."""
    total = 0.0
    for geom in lines.geometry:
        g = gpd.GeoSeries([geom], crs=4326).to_crs(32633).iloc[0]
        steps = max(int(g.length / 1000 / (cell_km_x / 2)), 1)
        pts = gpd.GeoSeries([g.interpolate(t, normalized=True) for t in np.linspace(0, 1, steps + 1)],
                            crs=32633).to_crs(4326)
        vals = []
        for p in pts:
            r, c = to_rc(p.x, p.y)
            r, c = min(max(r, 0), nrows - 1), min(max(c, 0), ncols - 1)
            vals.append(cost[r, c] if np.isfinite(cost[r, c]) else WATER_COST)
        total += np.mean(vals) * g.length / 1000
    return total


roads_in = gpd.clip(roads, border)
summary = pd.DataFrame([
    {"Network": "Least-cost network (MST on least-cost paths)", "Links": len(lcn_lines),
     "Length (km)": lcn_lines.length_km.sum(), "Construction cost (index)": construction_cost_of(lcn_lines)},
    {"Network": "Straight-line MST", "Links": len(euc_lines),
     "Length (km)": euc_lines.length_km.sum(), "Construction cost (index)": construction_cost_of(euc_lines)},
    {"Network": "Existing trunk roads (Natural Earth)", "Links": len(roads_in),
     "Length (km)": roads_in.to_crs(32633).length.sum() / 1000, "Construction cost (index)": np.nan},
]).set_index("Network")
summary["Cost per km"] = summary["Construction cost (index)"] / summary["Length (km)"]

# How much of the least-cost network follows today's roads?
buf = roads_in.to_crs(32633).buffer(10_000).union_all()
lcn_m = lcn_lines.to_crs(32633)
share_on_roads = lcn_m.intersection(buf).length.sum() / lcn_m.length.sum()
euc_m = euc_lines.to_crs(32633)
share_euc_on_roads = euc_m.intersection(buf).length.sum() / euc_m.length.sum()
summary.loc["Least-cost network (MST on least-cost paths)", "Share within 10 km of existing road"] = share_on_roads
summary.loc["Straight-line MST", "Share within 10 km of existing road"] = share_euc_on_roads
summary.to_csv(ROOT / "results/network_comparison.csv", float_format="%.3f")

# ------------------------------------------------------------------ market access
def travel_cost_surface(road_lines):
    """Travel cost per km: 1 on a road cell, OFFROAD_SPEED_RATIO elsewhere, plus terrain."""
    on_road = features.rasterize([(g, 1) for g in road_lines.geometry], out_shape=z.shape,
                                 transform=transform, fill=0, all_touched=True).astype(bool) \
        if len(road_lines) else np.zeros(z.shape, bool)
    t = np.where(on_road, 1.0, OFFROAD_SPEED_RATIO * (1 + 0.1 * slope_pct))
    t[~inside] = np.inf
    return t


def market_access(tcost):
    pop = cities["pop_max"].values.astype(float)
    ma = np.zeros(n)
    for i in range(n):
        mcp = MCP_Geometric(tcost, sampling=(cell_km_y, cell_km_x))
        cum, _ = mcp.find_costs([cities.at[i, "rc"]])
        tij = np.array([cum[cities.at[j, "rc"]] for j in range(n)])
        mask = np.arange(n) != i
        ma[i] = np.sum(pop[mask] * tij[mask] ** (-THETA))
    return ma


scen = {
    "No network": gpd.GeoDataFrame(geometry=[], crs=4326),
    "Existing roads": roads_in,
    "Least-cost network": lcn_lines,
    "Straight-line MST": euc_lines,
}
ma = pd.DataFrame({k: market_access(travel_cost_surface(v)) for k, v in scen.items()},
                  index=cities["name"])
ma["Gain from existing roads (%)"] = 100 * (ma["Existing roads"] / ma["No network"] - 1)
ma["Gain from least-cost network (%)"] = 100 * (ma["Least-cost network"] / ma["No network"] - 1)
ma["Least-cost vs existing (%)"] = 100 * (ma["Least-cost network"] / ma["Existing roads"] - 1)
ma.insert(0, "population", cities["pop_max"].values)
ma.insert(1, "region", cities["adm1name"].values)
ma.to_csv(ROOT / "results/market_access.csv", float_format="%.2f")

agg = pd.Series({
    k: np.average(ma[k], weights=ma["population"]) for k in scen
}, name="Population-weighted market access")
agg_df = agg.to_frame()
agg_df["Gain vs no network (%)"] = 100 * (agg / agg["No network"] - 1)
agg_df["Network length (km)"] = [0, summary.iloc[2]["Length (km)"], summary.iloc[0]["Length (km)"],
                                 summary.iloc[1]["Length (km)"]]
agg_df["MA gain per 1,000 km (%)"] = agg_df["Gain vs no network (%)"] / agg_df["Network length (km)"] * 1000
agg_df.to_csv(ROOT / "results/market_access_summary.csv", float_format="%.3f")

# ------------------------------------------------------------------ sensitivity: does the LCN change with the cost weights?
def lcn_edges(slope_w, water_c):
    c2 = 1 + slope_w * slope_pct + water_c * water
    c2[~inside] = np.inf
    g = nx.Graph()
    for i in range(n):
        mcp = MCP_Geometric(c2, sampling=(cell_km_y, cell_km_x))
        cum, _ = mcp.find_costs([cities.at[i, "rc"]])
        for j in range(i + 1, n):
            g.add_edge(i, j, weight=cum[cities.at[j, "rc"]])
    return {tuple(sorted(e)) for e in nx.minimum_spanning_tree(g).edges()}


base = {tuple(sorted(e)) for e in mst_lcp.edges()}
sens = []
for sw, wc in [(0.0, 25), (0.25, 25), (1.0, 25), (0.5, 5), (0.5, 50)]:
    e = lcn_edges(sw, wc)
    sens.append({"Slope weight": sw, "Water cost": wc, "Links shared with baseline": len(e & base),
                 "Total links": len(base)})
pd.DataFrame(sens).to_csv(ROOT / "results/sensitivity.csv", index=False)

# ------------------------------------------------------------------ maps
extent = [lon0 - res / 2, lon0 + res * (ncols - 0.5), lat0 - res * (nrows - 0.5), lat0 + res / 2]
zz = np.where(inside, z, np.nan)
fig, axes = plt.subplots(1, 2, figsize=(12, 8.5))
for ax, title in zip(axes, ["Least-cost network vs existing trunk roads", "Straight-line MST (naive benchmark)"]):
    ax.imshow(zz, extent=extent, cmap="terrain", vmin=-550, vmax=2100, alpha=0.75)
    border.boundary.plot(ax=ax, color="black", lw=0.8)
    rivers.plot(ax=ax, color="#3b8bc2", lw=0.6, alpha=0.7)
    ax.set_xlim(8.3, 16.4)
    ax.set_ylim(1.4, 13.3)
    ax.set_title(title, fontsize=11)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
roads_in.plot(ax=axes[0], color="#555", lw=1.2, ls="--", label="Existing trunk roads")
lcn_lines.plot(ax=axes[0], color="#c0392b", lw=2.2, label="Least-cost network")
euc_lines.plot(ax=axes[1], color="#8e44ad", lw=2.2, label="Straight-line MST")
for ax in axes:
    ax.scatter(cities.geometry.x, cities.geometry.y, s=np.sqrt(cities.pop_max) / 6, color="black",
               zorder=5, edgecolor="white", lw=0.5)
    for _, row in cities.iterrows():
        if row.pop_max > 150_000 or row.featurecla != "Populated place":
            ax.annotate(row["name"], (row.geometry.x, row.geometry.y), fontsize=7, xytext=(3, 3),
                        textcoords="offset points")
    ax.legend(loc="lower right", fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
fig.suptitle("A least-cost highway network for Cameroon's 23 largest cities", fontsize=13)
fig.tight_layout()
fig.savefig(ROOT / "figures/network_maps.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(7, 8.5))
im = ax.imshow(np.where(inside, np.log1p(cost), np.nan), extent=extent, cmap="magma_r")
border.boundary.plot(ax=ax, color="black", lw=0.8)
lcn_lines.plot(ax=ax, color="#00b3b3", lw=2)
ax.set_xlabel("")
ax.set_ylabel("")
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

top = ma.sort_values("Least-cost vs existing (%)")
fig, ax = plt.subplots(figsize=(8, 6.5))
ax.barh(top.index, top["Least-cost vs existing (%)"],
        color=np.where(top["Least-cost vs existing (%)"] >= 0, "#2ca02c", "#c0392b"))
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("Market access under the least-cost network relative to existing roads (%)")
ax.set_title("Who would gain from a least-cost network?")
fig.tight_layout()
fig.savefig(ROOT / "figures/market_access_gains.png")
plt.close(fig)

pd.set_option("display.width", 200)
print(summary.round(3)); print(agg_df.round(3)); print(pd.DataFrame(sens))
print(ma.round(1).sort_values("Least-cost vs existing (%)"))
print(lcn_lines.drop(columns="geometry").round(1))
