# Algorithmes de graphes écrits à la main : Dijkstra sur une grille et arbre couvrant minimal (Kruskal).

import heapq
import math

# Les 8 voisins d'une case : (déplacement en ligne, déplacement en colonne)
NEIGHBOURS = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


# ------------------------------------------------------------------------------
# ALGORITHME DijkstraSurGrille
# ENTREES : cost[r][c] = coût par km pour traverser la case (infini = interdit)
#           depart     = case de départ (ligne, colonne)
#           dy, dx     = taille d'une case en km (verticale, horizontale)
# SORTIES : dist[r][c] = coût total du chemin le moins cher depuis le départ
#           prev[r][c] = case précédente sur ce chemin (pour reconstruire le tracé)
# DEBUT
#   dist <- infini partout ; dist[depart] <- 0
#   file <- {(0, depart)}                      (file de priorité : le plus petit coût sort en premier)
#   TANT QUE file non vide FAIRE
#       (d, case) <- retirer le plus petit de file
#       SI d > dist[case] ALORS continuer      (ancienne entrée déjà améliorée)
#       POUR chaque voisin v de case FAIRE
#           longueur <- dy, dx ou la diagonale racine(dy² + dx²)
#           pas <- longueur x (cost[case] + cost[v]) / 2
#           SI dist[case] + pas < dist[v] ALORS
#               dist[v] <- dist[case] + pas ; prev[v] <- case ; ajouter (dist[v], v) dans file
#   FIN TANT QUE
# FIN
# ------------------------------------------------------------------------------
def dijkstra_grid(cost, start, dy, dx):
    n_rows = len(cost)
    n_cols = len(cost[0])
    diagonal = math.sqrt(dy * dy + dx * dx)

    dist = []
    prev = []
    for r in range(n_rows):
        dist.append([math.inf] * n_cols)
        prev.append([None] * n_cols)

    r0, c0 = start
    dist[r0][c0] = 0.0
    queue = [(0.0, r0, c0)]

    while len(queue) > 0:
        d, r, c = heapq.heappop(queue)
        if d > dist[r][c]:
            continue
        for step_r, step_c in NEIGHBOURS:
            nr = r + step_r
            nc = c + step_c
            if nr < 0 or nr >= n_rows or nc < 0 or nc >= n_cols:
                continue
            if cost[nr][nc] == math.inf:
                continue
            if step_r != 0 and step_c != 0:
                length = diagonal
            elif step_r != 0:
                length = dy
            else:
                length = dx
            new_d = d + length * (cost[r][c] + cost[nr][nc]) / 2
            if new_d < dist[nr][nc]:
                dist[nr][nc] = new_d
                prev[nr][nc] = (r, c)
                heapq.heappush(queue, (new_d, nr, nc))
    return dist, prev


def trace_path(prev, end):
    """Remonte les cases précédentes depuis l'arrivée jusqu'au départ."""
    path = [end]
    cell = end
    while prev[cell[0]][cell[1]] is not None:
        cell = prev[cell[0]][cell[1]]
        path.append(cell)
    path.reverse()
    return path


# ------------------------------------------------------------------------------
# ALGORITHME Kruskal (arbre couvrant de coût minimal)
# ENTREE : liste d'arêtes (coût, i, j) entre n villes
# DEBUT
#   trier les arêtes par coût croissant
#   chaque ville forme son propre groupe
#   POUR chaque arête (coût, i, j) dans l'ordre FAIRE
#       SI i et j ne sont pas déjà dans le même groupe ALORS
#           garder l'arête ; fusionner les deux groupes
#   RETOURNER les n - 1 arêtes gardées
# FIN
# On relie ainsi toutes les villes au coût total le plus faible, sans jamais faire de boucle.
# ------------------------------------------------------------------------------
def kruskal(n, edges):
    edges = sorted(edges)
    group = list(range(n))                          # group[i] = représentant du groupe de i

    def find(i):
        while group[i] != i:
            i = group[i]
        return i

    kept = []
    for weight, i, j in edges:
        gi = find(i)
        gj = find(j)
        if gi != gj:
            kept.append((i, j))
            group[gj] = gi
        if len(kept) == n - 1:
            break
    return kept
