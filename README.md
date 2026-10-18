# A Least-Cost Highway Network for Cameroon

If a planner had to connect Cameroon's 23 largest cities **at minimum construction cost**, which
network would she build? How does it compare with the trunk roads that exist today, and which cities
would gain or lose market access?

The method follows the **least-cost spanning network** used by Faber (2014) to build a counterfactual
of China's highway system, combined with the terrain-based construction costs common in the spatial
economics of African road networks (e.g. the work of Carl Müller-Crepon and co-authors).

![Maps](figures/network_maps.png)

## Method

1. **Construction-cost surface** on a 2.5 arc-minute grid (about 4.6 km):
   `cost per km = 1 + 0.5 × slope (%) + 25 × water`, where slope comes from a global 10 arc-minute
   topography grid and water from major rivers and lakes (bridges are expensive). Cells outside
   Cameroon are excluded: the network must be domestic.
2. **Least-cost paths** between all 253 pairs of cities (Dijkstra on the raster, `skimage.graph`).
3. **Minimum spanning tree** on the pairwise least-cost values: the cheapest network that connects every
   city, the **least-cost network (LCN)**. A straight-line MST serves as a naive benchmark.
4. **Market access** of each city, `MA_i = Σ_j Pop_j × T_ij^(-1)`, where `T_ij` is the travel cost
   on a surface where moving on a road costs 1 per km and moving off-road costs 5 per km (more on
   slopes). Computed with no network, today's trunk roads, the LCN and the straight-line MST.

## Results

### The cheapest network

| Network | Links | Length | Construction cost (index) | Cost per km | Share within 10 km of an existing road |
|---|---|---|---|---|---|
| **Least-cost network** | 22 | 2,288 km | **2,987** | **1.31** | **63%** |
| Straight-line MST | 22 | 1,983 km | 5,233 | 2.64 | 74% |
| Existing trunk roads | 247 segments | 10,068 km | | | |

- **Terrain matters a lot:** the straight-line network is 13% shorter but would cost **75% more** to
  build, because it climbs straight across the Western Highlands and the Adamawa plateau and crosses
  rivers where they are wide. The LCN avoids them: it reaches Kumbo from Bamenda rather than straight
  over the highlands from Foumban, and its paths bend along valleys instead of following straight lines.
- **Two thirds of the least-cost network lies within 10 km of today's trunk roads**: Cameroon's
  main corridors broadly follow the cheapest routes, especially in the West (Douala, Bafoussam,
  Bamenda) and the North (Ngaoundéré, Garoua, Maroua).
- The result is robust to the cost weights: changing the slope weight from 0 to 1 or the water cost from
  5 to 50 keeps at least 19 of the 22 links (`results/sensitivity.csv`).

![Cost surface](figures/cost_surface.png)

### Market access

| Scenario | Population-weighted MA gain vs no network | Network length | Gain per 1,000 km |
|---|---|---|---|
| Existing trunk roads | +371% | 10,068 km | 37 points |
| **Least-cost network** | **+275%** | **2,288 km** | **120 points** |
| Straight-line MST | +278% | 1,983 km | 140 points |

- The LCN delivers **three quarters of the market-access gain of today's network with less than a
  quarter of its length**: a cheap backbone captures most of the value.
- **But a pure cost-minimising tree misallocates links:** it never connects Douala and Yaoundé directly
  (the tree goes Douala, Edéa, Kribi, Ebolowa, Mbalmayo, Yaoundé), because a spanning tree only
  minimises total cost and ignores traffic. As a result, Yaoundé (-37%), Edéa (-39%) and Bafia (-52%)
  would have much lower market access than with today's roads, while Kumba and Ebolowa
  would gain slightly.

![Market access](figures/market_access_gains.png)

**Takeaway:** the existing trunk network is consistent with least-cost engineering logic, and the
places where it departs from the LCN (the Douala to Yaoundé corridor) are exactly where traffic,
not construction cost, justifies extra links. This is why Faber uses the LCN as an *instrument*:
it predicts which places get connected for cost reasons, independent of their economic potential.

## Limits

10 arc-minute terrain smooths mountains (Mount Cameroon appears at under 2,000 m); only major rivers
are included; city populations are Natural Earth estimates; travel costs are stylised, not
calibrated on real speeds. Natural extensions: SRTM 90 m terrain, OpenStreetMap roads with surface
type, a traffic-weighted network (Steiner tree or gravity-weighted links), and use of the LCN as an
instrument for road access in a regression on local development (night lights, census data).

## Run it

```bash
pip install -r requirements.txt
python run_analysis.py      # about 10 seconds
```

```
data/cities.geojson                         23 cities above 50,000 inhabitants (Natural Earth)
data/cameroon_boundary.geojson              national border (Natural Earth)
data/existing_roads_natural_earth.geojson   trunk roads (Natural Earth 1:10m)
data/rivers.geojson, data/lakes.geojson     water bodies (Natural Earth 1:10m)
data/topography_10arcmin_cameroon.nc        topography (fatiando-data, ETOPO1 at 10 arc-minutes)
results/                                     network links, comparisons, market access, sensitivity
figures/                                     maps and charts
```

**Reference:** Faber, B. (2014). Trade Integration, Market Size, and Industrialization: Evidence
from China's National Trunk Highway System. *Review of Economic Studies*, 81(3), 1046-1070.

**Author:** Klaus KUEKOU
