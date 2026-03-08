### Dependencies ###

# Standard library
import ast
import itertools
import json
import math
import os
import warnings
import random
import time
import argparse
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.interpolate import PchipInterpolator
from scipy.optimize import LinearConstraint, curve_fit, minimize
from scipy.sparse.csgraph import shortest_path
from scipy.stats import skew, kurtosis
from openpyxl import Workbook, load_workbook
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = REPO_ROOT / "data"
FIGURES_DIR = REPO_ROOT / "figures"
SENSITIVITY_DIR = DATA_DIR / "sensitivity_analysis"

for d in (DATA_DIR, FIGURES_DIR, SENSITIVITY_DIR):
    d.mkdir(parents=True, exist_ok=True)



### Main functions ###

# Notch-Delta simulation
def simulate(M, t_final, dt, signalling_labels=None, k=2, h=2, Ka=0.1, Kr=0.001, nu=1, randomQ=True):

    if not randomQ:
        np.random.seed(42)

    matrix = M
    
    if signalling_labels is not None:
        mask = np.zeros(M.shape[0], bool)
        mask[signalling_labels] = True
        M_eff = M.copy()
        M_eff[~mask, :] = 0
        M_eff[:, ~mask] = 0
    else:
        M_eff = M
        mask = np.ones(M.shape[0], bool)

    N = M.shape[0]
    y = np.random.rand(2*N) * 0.1
    steps = int(t_final / dt)
    for _ in range(steps):
        n, d = y[:N], y[N:]
        avg = M_eff.dot(d)
        dn = avg**k / (Ka**k + avg**k) - n
        dd = nu * (Kr**h / (Kr**h + n**h) - d)
        y[:N] += dt * dn
        y[N:] += dt * dd
        if signalling_labels is not None:
            y[N:][~mask] = 0

    notch = y[:N].copy()
    delta = y[N:].copy()

    if signalling_labels is not None:
        notch[~mask] = 0
        delta[~mask] = 0
        
    return notch, delta

# Auxiliary functions for eps-Collier model
def long_range_M(M, cent_list, prot_len, eps, labels):

    M = np.asarray(M)
    N = M.shape[0]

    cent = np.array(cent_list, float)

    dx = cent[:, None, 0] - cent[None, :, 0]
    dy = cent[:, None, 1] - cent[None, :, 1]
    dist = np.sqrt(dx*dx + dy*dy)

    M_long = (dist <= prot_len).astype(float)
    M_long[M != 0] = 0

    np.fill_diagonal(M_long, 0)

    labels_out = np.setdiff1d(np.arange(N), labels)
    M_long[labels_out, :] = 0
    M_long[:, labels_out] = 0

    M_combined = (1 - eps) * M + eps * M_long

    return M_combined
def min_long_range_distance(M, centroids, labels):
    M = np.asarray(M)
    cent = np.array(centroids, float)
    labels = np.array(labels, int)

    cent_lab = cent[labels]
    dist = np.linalg.norm(cent_lab[:, None, :] - cent_lab[None, :, :], axis=-1)

    M_lab = M[np.ix_(labels, labels)]

    np.fill_diagonal(dist, np.inf)

    dist_non_nb = np.where(M_lab == 0, dist, np.inf)

    return dist_non_nb.min()

# Layer weight computation
def compute_omega_k(omega_func, wing_region, Lmax=None, method='simpson', num_simpson=101):
    
    gap = gap_dict[wing_region]
    n   = n_dict[wing_region]

    def omega_cut(z):
        z = np.asarray(z)
        vals = omega_func(z)
        if Lmax is not None:
            mask = (z >= 0) & (z <= Lmax)
            vals = np.where(mask, vals, 0.0)
        return vals

    if method == 'midpoint':
        mid_z = (np.arange(n) + 0.5) * gap
        return omega_cut(mid_z) * gap

    elif method == 'trapezoid':
        z_edges = np.linspace(0, n * gap, n + 1)
        y = omega_cut(z_edges)
        return (gap / 2) * (y[:-1] + y[1:])

    elif method == 'simpson':
        if num_simpson % 2 == 0:
            raise ValueError("num_simpson must be odd for Simpson's rule")
        omega_k = np.empty(n)
        for k in range(n):
            a, b = k * gap, (k + 1) * gap
            zs = np.linspace(a, b, num_simpson)
            ys = omega_cut(zs)
            h  = (b - a) / (num_simpson - 1)
            omega_k[k] = (h / 3) * (
                ys[0]
                + ys[-1]
                + 4 * ys[1:-1:2].sum()
                + 2 * ys[2:-1:2].sum()
            )
        return omega_k

    else:
        raise ValueError(f"Unknown method '{method}'")

# Load centroid data
def load_centroids(wing_region):
    path = DATA_DIR / f"data_cell_geometry_{wing_region}.csv"
    df = pd.read_csv(path)
    max_frame = df["frame"].max()
    centroid_lists = []
    for frame in range(max_frame + 1):
        subdf = df[df["frame"] == frame]
        arr = subdf[["centroid_x", "centroid_y"]].to_numpy()
        centroid_lists.append(arr)
    return centroid_lists

# Load apical area data
def load_areas(wing_region):
    path = DATA_DIR / f"data_cell_geometry_{wing_region}.csv"
    df = pd.read_csv(path)
    subdf = df[df["frame"] == 0]
    area_list = subdf["area"].to_numpy()
    return area_list

# Color function
def shade(color, factor):
    rgb = np.array(mcolors.to_rgb(color))
    return mcolors.to_hex(rgb*factor + (1-factor))

# SOP spacing plots
def fancy_plot(spacing_dict, Lmax_list, weight_type, wing_regions,
               saveQ=False, fig_size=(5,4), xlim=(0,25), ylim=(0.95, 2.35), smooth=0,
               degenplotQ=False, errorbarQ=False, colors = [shade("#1f77b4", 0.6),shade("#1f77b4", 1),shade("#1f77b4", 0.4)],
               title='sop_spacing', x_title='Depth (Âµm)',
               legend_loc='lower right', meanQ=False, meancolor='black', mergeQ=False):

    plt.figure(figsize=fig_size)
    labels = [f"WD {i+1}" for i in range(len(wing_regions))]
    x = np.array(Lmax_list)

    all_y = []

    for color, label, region in zip(colors, labels, wing_regions):
        data   = spacing_dict[region]
        y      = np.array([d for d, vr, flag in data])
        ystd   = np.array([vr for d, vr, flag in data])
        flags  = np.array([flag for d, vr, flag in data])

        inf_idx = np.where(~np.isfinite(y))[0]
        if inf_idx.size:
            flags[inf_idx] = True
            for i in inf_idx:
                if i > 0:
                    y[i] = y[i-1]
                else:
                    y[i] = y[np.isfinite(y)][0]

        if meanQ or mergeQ:
            all_y.append(y.copy())

        if mergeQ:
            continue

        if errorbarQ:
            plt.errorbar(x, y, yerr=ystd, fmt='none',
                         ecolor=color, alpha=0.5, capsize=3)

        if smooth > 0 and len(x) > 1:
            xs    = np.linspace(x.min(), x.max(), 300)
            pchip = PchipInterpolator(x, y)
            ys    = pchip(xs)
            plt.plot(xs, ys, color=color, linewidth=3, label=label)
        else:
            plt.plot(x, y, color=color, linewidth=3, label=label)

        if degenplotQ and flags.any():
            idxs = np.where(flags)[0]
            runs = np.split(idxs, np.where(np.diff(idxs) > 1)[0]+1)
            for run in runs:
                xr = x[run]; yr = y[run]
                if len(xr) > 1:
                    plt.plot(xr, yr, '-', color='white', linewidth=4, zorder=5)
                    plt.plot(xr, yr, '--', color='red',   linewidth=2, zorder=6)
                else:
                    plt.plot(xr, yr, 'o', color='white', markersize=8, zorder=5)
                    plt.plot(xr, yr, 'o', color='red',   markersize=5, zorder=6)

    if (meanQ or mergeQ) and len(all_y) > 0:
        Y = np.vstack(all_y)
        ymean = Y.mean(axis=0)
        ysd   = Y.std(axis=0)

        if smooth > 0 and len(x) > 1:
            xs    = np.linspace(x.min(), x.max(), 300)
            pchip_m = PchipInterpolator(x, ymean)
            ys_m    = pchip_m(xs)
            plt.plot(xs, ys_m, color=meancolor, linewidth=4, label='Mean')

            if mergeQ:
                pchip_u = PchipInterpolator(x, ymean + ysd)
                pchip_l = PchipInterpolator(x, ymean - ysd)
                plt.fill_between(xs, pchip_l(xs), pchip_u(xs), color=meancolor, alpha=0.2, linewidth=0)
        else:
            plt.plot(x, ymean, color=meancolor, linewidth=4, label='Mean')
            if mergeQ:
                plt.fill_between(x, ymean - ysd, ymean + ysd, color=meancolor, alpha=0.2, linewidth=0)

    plt.xlabel(x_title, fontsize=22)
    plt.ylabel('Mean SOP spacing', fontsize=22)
    plt.tick_params(labelsize=18)
    plt.xlim(xlim)
    plt.ylim(ylim); plt.grid(True, linestyle='--', alpha=0.6)

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    if degenplotQ and not mergeQ:
        deg_h = Line2D([0], [0], color='red', linestyle='--', linewidth=2)
        handles.append(deg_h); labels.append('Degen.')
    ax.legend(handles, labels, fontsize=12, frameon=True, loc=legend_loc) \
      .get_frame().set_edgecolor('black')
    ax.set_position([0.10, 0.15, 0.85, 0.75])

    ax = plt.gca()
    handles, labels = ax.get_legend_handles_labels()
    if degenplotQ and not mergeQ:
        deg_h = Line2D([0], [0], color='red', linestyle='--', linewidth=2)
        handles.append(deg_h); labels.append('Degen.')

    handles = []
    leg = ax.legend(
        handles, labels,
        fontsize=12,
        frameon=False,
        loc=legend_loc
    )

    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_alpha(0.8)
    leg.set_zorder(10)
    leg.get_frame().set_edgecolor('black')

    if saveQ:
        plt.savefig(FIGURES_DIR / f'{title}_{weight_type}.pdf',
                    bbox_inches='tight', transparent=True)
    plt.show()

# Plots
def plot_layer_graph(delta, centroids, A_layers, framenumber=0,
                     wing_region='',
                     sop_idx=None, adj0=None,
                     draw_rect=False, height=None, y_shift=0.0,
                     show_labels=True, figsize=(7,7),
                     node_size=300, label_font_size=8,
                     edge_color='black', edge_width=1.0,
                     show_other_layers=True,
                     other_edge_color='lightblue', other_edge_width=0.5,
                     saveQ=False):

    pos = centroids[framenumber]
    adj = A_layers[framenumber]
    valid = np.nonzero(adj.sum(0) > 0)[0]

    fig, ax = plt.subplots(figsize=figsize)

    if show_other_layers:
        union_adj = np.zeros_like(adj, dtype=bool)
        for k, m in enumerate(A_layers):
            if k == framenumber:
                continue
            union_adj |= (m > 0)

        r_u, c_u = np.nonzero(union_adj)
        mask_u   = (r_u < c_u) & np.isin(r_u, valid) & np.isin(c_u, valid)
        segs_u = [((pos[i,0], pos[i,1]), (pos[j,0], pos[j,1]))
                  for i, j in zip(r_u[mask_u], c_u[mask_u])]

        if segs_u:
            ax.add_collection(LineCollection(
                segs_u,
                colors=other_edge_color,
                linewidths=other_edge_width,
                zorder=0
            ))

    r, c = np.nonzero(adj)
    mask_f = (r < c) & np.isin(r, valid) & np.isin(c, valid)
    segs_f = [((pos[i,0], pos[i,1]), (pos[j,0], pos[j,1]))
              for i, j in zip(r[mask_f], c[mask_f])]

    if segs_f:
        ax.add_collection(LineCollection(
            segs_f,
            colors=edge_color,
            linewidths=edge_width,
            zorder=1
        ))

    pts = pos[valid]
    cols = np.zeros((len(valid), 3))
    cols[:, 0] = delta[valid]
    ax.scatter(
        pts[:, 0], pts[:, 1],
        c=cols, s=node_size,
        edgecolors='none', zorder=2
    )

    if show_labels:
        for i in valid:
            ax.text(
                pos[i, 0], pos[i, 1], str(i),
                color='white', ha='center', va='center',
                fontsize=label_font_size, zorder=3
            )

    if draw_rect:
        x, y = pos[:,0], pos[:,1]
        x_min, x_max = x.min(), x.max()
        y_min, y_max = y.min(), y.max()
        y0 = y_min + y_shift * ((y_max - height) - y_min)

        rect = Rectangle(
            (x_min, y0),
            x_max - x_min, height,
            edgecolor='magenta', facecolor='none',
            linewidth=2, zorder=4
        )
        ax.add_patch(rect)

        graph = (adj0 > 0).astype(int)
        full_dist = shortest_path(csgraph=graph, directed=False, unweighted=True)
        d = band_min_sop_distance_at_shift(
            cents=np.array(centroids[0]),
            sop_idx=sop_idx,
            full_dist=full_dist,
            height=height,
            ys=y_shift
        )
        print(f"Mean SOP spacing at band shift={y_shift:.2f}: {d:.3f}")

    ax.set_aspect('equal')
    ax.axis('off')

    if saveQ:
        plt.savefig(
            FIGURES_DIR / f'graph_plot_{wing_region}.pdf',
            dpi=300, bbox_inches='tight', transparent=False
        )
    plt.show()

# SOP distance function (shift)
def band_min_sop_distance_at_shift(cents, sop_idx, full_dist, height=50, ys=0.):

    x = cents[:,0]
    y = cents[:,1]
    x_min, x_max = x.min(), x.max()
    y_min, y_max = y.min(), y.max()

    y0 = y_min + ys * ((y_max - height) - y_min)

    in_rect = np.where(
        (x >= x_min) & (x <= x_max) &
        (y >= y0)    & (y <= y0 + height)
    )[0]
    sops = [i for i in sop_idx if i in in_rect]
    if len(sops) < 2:
        return [np.nan, np.nan]

    dm = full_dist[np.ix_(sops, sops)]
    np.fill_diagonal(dm, np.inf)
    mins = np.min(dm, axis=1) - 1
    finite = mins[np.isfinite(mins)]
    
    return [finite.mean(), finite.std()] if finite.size else [np.nan, np.nan]

# SOP distance function (average)
def band_avg_min_sop_distance(cents, sop_idx, adj0, height, y_shift_steps):

    cents = np.asarray(cents)
    graph = (adj0 > 0).astype(int)
    full_dist = shortest_path(csgraph=graph, directed=False, unweighted=True)

    d_list = []
    v_list = []
    for ys in np.linspace(0, 1, y_shift_steps):
        d, v = band_min_sop_distance_at_shift(
            cents, sop_idx, full_dist, height, ys
        )
        if not np.isnan(d):
            d_list.append(d)
            v_list.append(v)
            
    d_list = [x for x in d_list if x <= 6.5]

    return [float(np.mean(d_list)), float(np.mean(v_list))] if d_list else [np.inf, np.inf]

#  SOP distance function (master)
def compute_band_distance(wing_region,
                          weight_list=None,
                          omega_func=None,
                          t_final=1000.,
                          dt=0.1,
                          Lmax=None,
                          height=80,
                          y_shift_steps=20,
                          sim_number=1,
                          quad_method='simpson',
                          alpha=0.,
                          degen_T=1.0,
                          plotQ=False, graphsaveQ=False, normalQ=False, str_type='centroid',
                          epsmodelQ=False, eps=0., prot_len=0., marker_type='delta',
                          k=2, h=2, Ka=0.1, Kr=0.001, nu=1,
                          show_labels=True, show_other_layers=True,
                          randomQ=True,
                          **quad_kwargs):
    
    if weight_list is None:
        if omega_func is None:
            raise ValueError("Either weight_list or omega_func must be provided")
        if normalQ: # option to normalise the weights
            weight_list_base = compute_omega_k(omega_func, wing_region, Lmax=Lmax, method=quad_method, **quad_kwargs)
            weight_list_all = compute_omega_k(omega_func, wing_region, Lmax=32, method=quad_method, **quad_kwargs)
            weight_list = np.array([i*sum(weight_list_base)/sum(weight_list_all) for i in weight_list_base])
        else:
            weight_list = compute_omega_k(omega_func, wing_region, Lmax=Lmax, method=quad_method, **quad_kwargs)
    weight_list0 = weight_list
    
    thresholds = [0.1, 0.2, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
    results = {th: [] for th in thresholds}
    results_v = {th: [] for th in thresholds}
    
    if str_type == 'centroid':
        A_list = straight_adjacency(A_dict, centroids_dict, alpha)[wing_region]
    
    cents   = np.array(centroids_dict[wing_region][0])
    adj0    = A_dict[wing_region][0]
    labels  = signalling_labels_dict[wing_region]
    
    sop_pairs = []
    bimodal_list = []

    for sim in range(sim_number):
        wk    = weight_list
        if str_type == 'prune':
            straight_data = straight_prune(A_dict, centroids_dict, alpha)
            A_list, prune_idx = straight_data[0][wing_region], straight_data[1][wing_region]
            wk_new = [0 if i in prune_idx else val for i,val in enumerate(wk)]
            wk_new = [v * sum(wk) / sum(wk_new) for v in wk_new]
            wk = wk_new

        if epsmodelQ:
            dist_min = min_long_range_distance(A_dict[wing_region][0], centroids_dict[wing_region][0], signalling_labels_dict[wing_region])
            M = long_range_M(A_dict[wing_region][0], centroids_dict[wing_region][0], (1+prot_len)*(dist_min-1e-6), eps, signalling_labels_dict[wing_region])
        else:
            M = sum(w*A for w,A in zip(wk, A_list))
            
        notch, delta = simulate(M, t_final, dt, labels, k=k, h=h, Ka=Ka, Kr=Kr, nu=nu, randomQ=randomQ)

        bimodal_list.append(SOP_bimodality(delta[signalling_labels_dict[wing_region]], notch[signalling_labels_dict[wing_region]]))
        
        for th in thresholds:
            sop_idx = np.where(delta>th)[0]
            d, v = band_avg_min_sop_distance(cents, sop_idx, adj0, height, y_shift_steps)
            results[th].append(d)
            results_v[th].append(v)

        if plotQ:
            if marker_type == 'notch':
                marker = notch
            elif marker_type == 'delta':
                marker = delta
            plot_layer_graph(marker, centroids_dict[wing_region], A_dict[wing_region], wing_region=wing_region, framenumber=0, saveQ=graphsaveQ,
                             show_labels=show_labels, show_other_layers=show_other_layers)
            
        sop_pairs.append(np.triu(adj0[np.ix_(sop_idx, sop_idx)] > 0, k=1).sum())
    
    avg_results = {}
    for th, ds in results.items():
        finite = [x for x in ds if np.isfinite(x)]
        avg_results[th] = np.mean(finite) if finite else np.inf
        if th==0.1:
            ds0=ds
    avg_results_v = {}
    for th, ds in results_v.items():
        finite = [x for x in ds if np.isfinite(x)]
        avg_results_v[th] = np.mean(finite) if finite else np.inf

    sop_pairs_mean = np.mean(sop_pairs)
    degenQ = sop_pairs_mean > degen_T or (sum(np.isinf(x) and x > 0 for x in ds0) / len(ds0)) > 0.5
    avg_bimodal = np.mean(bimodal_list)

    return avg_results, avg_results_v, degenQ, avg_bimodal, notch, delta

# Straightening model
def straight_adjacency(A_dict, centroids_dict, alpha):
    
    straight = {}
    for region, layers in A_dict.items():
        A0 = layers[0].astype(float)
        if alpha >= 1.0:
            straight[region] = [A0.copy() for _ in layers]
            continue

        pos0 = np.asarray(centroids_dict[region][0])
        new_layers = [A0.copy()]

        for k in range(1, len(layers)):
            Ak   = layers[k].astype(float)
            posk = np.asarray(centroids_dict[region][k])
            disp = np.linalg.norm(posk - pos0, axis=1)

            bin0 = A0 > 0
            bink = Ak > 0

            ex_i, ex_j = np.where(bink & ~bin0)
            disp_ex    = [max(disp[i],disp[j]) for i,j in zip(ex_i,ex_j) if i<j]
            mi_i, mi_j = np.where(bin0 & ~bink)
            disp_mi    = [max(disp[i],disp[j]) for i,j in zip(mi_i,mi_j) if i<j]

            tau_rm  = np.quantile(disp_ex, 1 - alpha) if disp_ex else np.inf
            tau_add = np.quantile(disp_mi, alpha)     if disp_mi else -np.inf

            Anew = Ak.copy()

            for i,j in zip(ex_i, ex_j):
                if i<j and max(disp[i],disp[j]) >= tau_rm:
                    Anew[i,j] = Anew[j,i] = 0.0
            for i,j in zip(mi_i, mi_j):
                if i<j and max(disp[i],disp[j]) <= tau_add:
                    Anew[i,j] = Anew[j,i] = A0[i,j]

            new_layers.append(Anew)

        straight[region] = new_layers

    return straight

# Straightening plots
def plot_straightening_nonapical(A_str_dict, region, alphas, saveQ=False):

    A0 = A_str_dict[0.0][region][0] > 0
    deg_apical = A0.sum(axis=1)
    mask       = (deg_apical > 0)
    deg_apical = deg_apical[mask]
    deg_apical = deg_apical[deg_apical <= 12]

    deg_nonap_list = []
    for alpha in alphas:
        layers = A_str_dict[alpha][region]
        nonap = np.zeros_like(A0, dtype=bool)
        for Ak in layers[1:]:
            nonap |= (Ak > 0) & ~A0
        deg_nonap = nonap.sum(axis=1)[mask]
        deg_nonap = deg_nonap[deg_nonap <= 12]
        deg_nonap_list.append(deg_nonap)

    xs = np.arange(len(alphas)+1)

    fig, ax = plt.subplots(figsize=(8,3), constrained_layout=False)

    bp0 = ax.boxplot(
        [deg_apical], positions=[xs[0]], widths=0.6,
        patch_artist=True, manage_ticks=False
    )
    for box in bp0['boxes']:
        box.set(facecolor='#F06B68', edgecolor='black')
    for part in ('whiskers','caps','medians','fliers'):
        for line in bp0[part]:
            line.set(color='black', linewidth=1)

    bp1 = ax.boxplot(
        deg_nonap_list, positions=xs[1:], widths=0.6,
        patch_artist=True, manage_ticks=False
    )
    for box in bp1['boxes']:
        box.set(facecolor='#1f77b4', edgecolor='black', alpha=.7)
    for part in ('whiskers','caps','medians','fliers'):
        for line in bp1[part]:
            line.set(color='black', linewidth=1)

    labels = ['apical'] + [f"{int(100*alpha)}" for alpha in alphas]
    ax.set_xticks(xs)
    ax.set_xticklabels(labels, rotation=45)

    ax.set_xlabel('Straightening percentage (%)')
    ax.set_ylabel('Number of neighbours')
    ax.legend(
        [bp0["boxes"][0], bp1["boxes"][0]],
        ['Apical neighbours', 'Non-apical neighbours'],
        loc='upper right'
    )
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.set_position([0.10, 0.15, 0.85, 0.75])
    plt.title(f'WD {wd_dict[region]}')
    if saveQ:
        plt.savefig(FIGURES_DIR / f'straight_neighbour_{region}.pdf', bbox_inches='tight', transparent=True)
    plt.show()

    return deg_nonap_list

# Sensitivity analysis functions
def best_fit(L, best_array=[1., 2., 2., 1.]):
    nums, bools = L
    return float('inf') if any(bools) else np.linalg.norm(np.array(nums) - np.array(best_array))
def plot_heatmap(
    data,
    wd_label='wd_1',
    allQ=False,
    h_range=None,
    nu_range=None,
    logQ=False,
    show_label=False,
    figsize_auto=True,
    tick_step=None,
    saveQ=False,
    filename=None,
    tick_fontsize=12,
    title_name=''
):

    if allQ:
        hs_all  = sorted({h for w, h, n in data})
        nus_all = sorted({n for w, h, n in data})
        title = "all wing discs"
    else:
        hs_all  = sorted({h for w, h, n in data if w == wd_label})
        nus_all = sorted({n for w, h, n in data if w == wd_label})
        title = wd_label

    hs  = [h for h in hs_all  if h_range  is None or (h_range[0] <= h <= h_range[1])]
    nus = [n for n in nus_all if nu_range is None or (nu_range[0] <= n <= nu_range[1])]

    if allQ:
        M = np.empty((len(hs), len(nus)))
        M[:] = np.nan
        for i, h in enumerate(hs):
            for j, nu in enumerate(nus):
                vals = [data.get((w, h, nu), np.nan) for w in wing_discs]
                vals = np.array(vals, float)
                vals = vals[np.isfinite(vals)]
                M[i, j] = np.nan if vals.size == 0 else np.mean(vals)
    else:
        M = np.array([[data.get((wd_label, h, nu), np.nan) for nu in nus] for h in hs])

    M = np.ma.masked_where(~np.isfinite(M), M)
    cmap = plt.cm.viridis.copy()
    cmap.set_bad(color='white')

    if figsize_auto:
        w = max(4, len(nus) * 0.6)
        h_size = max(4, len(hs) * 0.6)
        plt.figure(figsize=(w, h_size))
    else:
        plt.figure()

    finite_vals = M.compressed()
    if logQ and np.any(finite_vals > 0):
        pos_vals = finite_vals[finite_vals > 0]
        norm = LogNorm(vmin=pos_vals.min(), vmax=pos_vals.max())
        im = plt.imshow(M, origin='upper', cmap=cmap, norm=norm)
    else:
        im = plt.imshow(M, origin='upper', cmap=cmap)

    def choose_tick_positions(n):
        if tick_step is not None:
            return list(range(0, n, tick_step))
        step = max(1, int(np.ceil(n / 10)))
        return list(range(0, n, step))

    xticks = choose_tick_positions(len(nus))
    yticks = choose_tick_positions(len(hs))
    plt.xticks(xticks, [f"{nus[i]:.2g}" for i in xticks], fontsize=tick_fontsize)
    plt.yticks(yticks, [f"{hs[i]:.2g}" for i in yticks], fontsize=tick_fontsize)
    plt.xlabel(r'$\nu$', fontsize=tick_fontsize)
    plt.ylabel(r'$h$', fontsize=tick_fontsize)
    plt.title(title, fontsize=tick_fontsize + 2)
    if show_label:
        for i in range(len(hs)):
            for j in range(len(nus)):
                if np.isfinite(M[i, j]):
                    plt.text(j, i, f'{M[i,j]:.2f}', ha='center', va='center', fontsize=tick_fontsize - 1)
    cbar = plt.colorbar(im)
    cbar.ax.tick_params(labelsize=tick_fontsize)
    plt.tight_layout()

    if saveQ:
        if filename is None:
            filename = title.replace(" ", "_") + ".pdf"
        plt.savefig(FIGURES_DIR / f"{title_name}_{filename}", transparent=True, bbox_inches='tight')
    plt.show()
    
def SOP_bimodality(delta, notch):
    corr_score = -np.corrcoef(delta, notch)[0,1]
    delta_bin = delta > np.median(delta)
    notch_bin = notch > np.median(notch)
    xor_score = np.mean(delta_bin ^ notch_bin)
    return corr_score
    
def _fmt_eta(seconds):
    seconds = max(0, int(round(seconds)))
    return f"{seconds//60:02d}:{seconds%60:02d}"

def longrun(Kr, h_list = np.arange(2, 8.1, 1), nu_list = [10**i for i in np.arange(-1., 0.1, 0.2)], sim_number = 50, threshold=0.1, k=2, Ka=0.1, dt=0.1):
    
    Lmax_list = np.linspace(0.5, 25, 2)
    alpha_list = np.linspace(0., 1., 2)
    degen_T = 1.
    normalQ = False    
    heatmap_dict = {}
    heatmap_dict_bimodal = {}
    start_time = time.time()
    total_iters = len(h_list) * len(nu_list)
    completed = 0
    
    for h_i in h_list:
        for nu_i in nu_list:
            spacing_dict_exp = {region: [] for region in wing_regions}
            bimodal_dict_exp = {region: [] for region in wing_regions}
            it = 0
            for Lmax in Lmax_list:
                for region in wing_regions:
                    d, vr, degenQ, avg_bimodal = compute_band_distance(
                        region,
                        omega_func=omega_exp,
                        Lmax=Lmax,
                        sim_number=sim_number,
                        quad_method='simpson',
                        height=heights_dict[region],
                        degen_T=degen_T,
                        normalQ=normalQ,
                        t_final=t_final,
                        y_shift_steps=20,
                        k=k, h=h_i, Ka=Ka, Kr=Kr, nu=nu_i, dt=dt
                    )
                    d = d[threshold]
                    vr = vr[threshold]
                    spacing_dict_exp[region].append([d, vr, degenQ])
                    bimodal_dict_exp[region].append(avg_bimodal)
                    it += 1
                    print(f'{it}/{len(Lmax_list)*len(wing_regions)}')#, end='\r')

            spacing_dict_exp_straight = {region: [] for region in wing_regions}
            bimodal_dict_exp_straight = {region: [] for region in wing_regions}
            it = 0
            for alpha in alpha_list:
                for region in wing_regions:
                    d, vr, degenQ, avg_bimodal = compute_band_distance(
                        region,
                        omega_func=omega_exp,
                        Lmax=25,
                        alpha=alpha,
                        sim_number=sim_number,
                        quad_method='simpson',
                        height=heights_dict[region],
                        degen_T=degen_T,
                        normalQ=normalQ,
                        t_final=t_final,
                        y_shift_steps=20,
                        str_type='centroid',
                        k=k, h=h_i, Ka=Ka, Kr=Kr, nu=nu_i, dt=dt
                    )
                    d = d[threshold]
                    vr = vr[threshold]
                    spacing_dict_exp_straight[region].append([d, vr, degenQ])
                    bimodal_dict_exp_straight[region].append(avg_bimodal)
                    it += 1
                    print(f'{it}/{len(alpha_list)*len(wing_regions)}')#, end='\r')

            for region in wing_regions:
                ll_dist = [[spacing_dict_exp[region][0][0],spacing_dict_exp[region][1][0],spacing_dict_exp_straight[region][0][0],spacing_dict_exp_straight[region][1][0]]]
                ll_bool = [[spacing_dict_exp[region][0][2],spacing_dict_exp[region][1][2],spacing_dict_exp_straight[region][0][2],spacing_dict_exp_straight[region][1][2]]]
                ll_all = ll_dist+ll_bool
                heatmap_dict[(region, h_i, nu_i)] = ll_all
                
                heatmap_dict_bimodal[(region, h_i, nu_i)] = [bimodal_dict_exp[region], bimodal_dict_exp_straight[region]]
    
            completed += 1
            elapsed = time.time() - start_time
            avg_time = elapsed / completed
            remaining = avg_time * (total_iters - completed)
            print(f"Running simulations for h={h_i}, nu={nu_i} "
                  f"({completed}/{total_iters}, "
                  f"elapsed {elapsed/60:.1f} min, "
                  f"~{remaining/60:.1f} min remaining)")
    
            print(f" â†’ Finished h={h_i}, nu={nu_i}\n")#, end='\r')
    
    total_time = (time.time() - start_time) / 60
    print(f"Total estimated runtime: {total_time:.1f} minutes")

    return heatmap_dict, heatmap_dict_bimodal

def runsens(Kr=10**-3, h_list=np.arange(2, 8.1, 1), nu_list = np.arange(0.075, 0.11, 0.0025), threshold=0.1, sim_number=2, reso='SMALL', sname='', bimodQ=False,
            k=2, Ka=0.1, dt=0.1):
    
    ttname1 = f'heatmap_spacing_{reso}_{sname}'
    ttname2 = f'heatmap_bimodal_{reso}_{sname}'
    heatmap_dict, heatmap_dict_bimodal = longrun(Kr, h_list=h_list, nu_list=nu_list, sim_number=sim_number, threshold=threshold, k=k, Ka=Ka, dt=dt)
    heatmap_error_dict = {i: best_fit(j) for i, j in heatmap_dict.items()}
    plot_heatmap(
        heatmap_error_dict,
        allQ=True,
        h_range=(2, 10),
        nu_range=(min(nu_list), 1000),
        logQ=False,
        tick_fontsize=10,
        saveQ=True,
        title_name=ttname1
    )

    d_to_save = {str(k): v for k, v in heatmap_error_dict.items()}
    with open(SENSITIVITY_DIR / f"{ttname1}.txt", "w", encoding="utf-8") as f:
        json.dump(d_to_save, f, indent=2)
    
    if bimodQ:
    
        heatmap_error_dict_bimodal = {}
        for key, val in heatmap_dict_bimodal.items():
            dist = np.linalg.norm(np.array([val[0][1], val[1][0]]) - 1)
            heatmap_error_dict_bimodal[key] = dist
        plot_heatmap(
            heatmap_error_dict_bimodal,
            allQ=True,
            h_range=(2, 10),
            nu_range=(0.05, 1000),
            logQ=False,
            tick_fontsize=10,
            saveQ=True,
            title_name=ttname2
        )
    
        d_to_save2 = {str(k): v for k, v in heatmap_error_dict_bimodal.items()}
        with open(SENSITIVITY_DIR / f"{ttname2}.txt", "w", encoding="utf-8") as f:
            json.dump(d_to_save2, f, indent=2)

def parse_key(k_str: str):

    try:
        return ast.literal_eval(k_str)
    except Exception:
        s = k_str.strip()
        s = re.sub(r"\barray\(\s*([^)]+?)\s*\)", r"\1", s)
        s = re.sub(r"\bnp\.float64\(\s*([^)]+?)\s*\)", r"\1", s)
        s = re.sub(r"\bfloat\(\s*([^)]+?)\s*\)", r"\1", s)
        s = re.sub(r"\bDecimal\(\s*([^)]+?)\s*\)", r"\1", s)
        return ast.literal_eval(s)

def plothm(path, cbar_range=None, title=None, debug_bad_keys=True, stitle='plot', saveQ=False):
    path = Path(path)
    if not path.is_absolute():
        candidate = REPO_ROOT / path
        if candidate.exists():
            path = candidate
    with path.open("r", encoding="utf-8") as f:
        d = json.loads(f.read())
    vals = defaultdict(list)
    for k_str, v in d.items():
        try:
            wd, k2, k3 = parse_key(k_str)
        except Exception:
            if debug_bad_keys:
                print("BAD KEY (could not parse):", k_str)
            raise
        vals[(float(k2), float(k3))].append(float(v))
    rows = [
        {"k2": k2, "k3": k3, "mean": float(np.mean(vs)), "n": len(vs)}
        for (k2, k3), vs in vals.items()
    ]
    df = pd.DataFrame(rows)
    bad = df[df["n"] != 3]
    if len(bad):
        print("Warning: some (k2,k3) pairs are missing wd entries:")
        print(bad.sort_values(["k2", "k3"]).to_string(index=False))
    k2_sorted = np.sort(df["k2"].unique())
    k3_sorted = np.sort(df["k3"].unique())
    grid = np.full((len(k2_sorted), len(k3_sorted)), np.nan, float)
    k2_to_i = {k2: i for i, k2 in enumerate(k2_sorted)}
    k3_to_j = {k3: j for j, k3 in enumerate(k3_sorted)}
    for _, r in df.iterrows():
        grid[k2_to_i[r["k2"]], k3_to_j[r["k3"]]] = r["mean"]
    plt.figure(figsize=(10, 5))
    imshow_kwargs = dict(origin="lower", aspect="auto")
    if cbar_range is not None:
        vmin, vmax = cbar_range
        imshow_kwargs.update(vmin=vmin, vmax=vmax)
    im = plt.imshow(grid, **imshow_kwargs)
    plt.colorbar(im, label="Mean value over wd_1-wd_3")
    plt.xticks(
        np.arange(len(k3_sorted)),
        [f"{x:.4g}" for x in k3_sorted],
        rotation=45,
        ha="right",
    )
    plt.yticks(
        np.arange(len(k2_sorted)),
        [f"{y:.4g}" for y in k2_sorted],
    )
    plt.xlabel("Third key (k3)")
    plt.ylabel("Second key (k2)")
    plt.title(title if title is not None else "Heatmap: mean over wd_1-wd_3")
    plt.tight_layout()
    if saveQ:
        plt.savefig(FIGURES_DIR / f'sens_{stitle}.pdf', bbox_inches='tight', transparent=True)
    plt.show()

# Distance heights
def compute_avg_count(cents, h1, n_shifts=20):
    y = cents[:, 1]
    y_min, y_max = y.min(), y.max()
    counts = []
    for s in np.linspace(0, 1, n_shifts):
        y0 = y_min + s * ((y_max - h1) - y_min)
        inside = (y >= y0) & (y <= y0 + h1)
        counts.append(inside.sum())
    return np.mean(counts)



### Parameter setup and data loading ###

# Notch-Delta parameters
k, h = 2, 2
Ka, Kr = 0.1, 0.001
nu = 1.
t_final = 1000.
dt = 0.1

# Regionâ€specific parameters
wing_regions = ['wd_1', 'wd_2', 'wd_3', 'wd_2_mbs', 'wd_3_mbs', 'wd_8_mbs']
wing_discs = ["wd_1", "wd_2", "wd_3"]
wd_dict = {'wd_1': '1', 'wd_2': '2', 'wd_3': '3', 'wd_2_mbs': '2mbs', 'wd_3_mbs': '3mbs', 'wd_8_mbs': '8mbs'}
gap_dict = {'wd_1': 0.5, 'wd_2': 0.3, 'wd_3': 0.5, 'wd_2_mbs': 0.3, 'wd_3_mbs': 0.3, 'wd_8_mbs': 0.3}
n_dict = {'wd_1': 50, 'wd_2': 105, 'wd_3': 60, 'wd_2_mbs': 80, 'wd_3_mbs': 80, 'wd_8_mbs': 91}
signalling_labels_dict = {
    'wd_1': np.array([23,25,26,27,28,29,34,36,37,41,44,45,47,48,51,52,54,55,
                              56,57,58,60,63,64,65,67,68,69,71,73,74,76,79,80,82,
                              83,84,86,88,89,91,92,93,95,96,97,98,99,103,106,121,
                              124,125,130,139], dtype=int),
    'wd_2': np.array([5,8,12,13,14,16,17,18,19,20,21,22,28,29,30,33,34,37,
                              38,39,40,41,42,43,46,48,49,51,52], dtype=int),
    'wd_3': np.array([20,22,27,28,29,30,31,34,35,36,37,38,39,42,43,45,48,49,
                              51,52,53,54,55,56,57,58,59,60,62,64,65,66,67,68,69,
                              70,72,73,76,77,78,94,96,101], dtype=int),
    'wd_2_mbs': np.array([1, 5, 6, 7, 8, 9, 10, 13, 14, 17, 21, 22, 23, 24, 25, 27,
                          29, 30, 31, 32, 35, 37, 41, 45, 48, 49, 50, 51, 60, 71, 72,
                          78, 79, 99], dtype=int),
    'wd_3_mbs': np.array([0, 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 17, 18, 19, 20, 23, 25,
                          26, 27, 28, 30, 31, 32, 33, 34, 36, 44, 45, 53, 58, 60, 61,
                          68, 69, 71, 72, 78, 83, 87, 93, 94], dtype=int),           ############## +100?????
    'wd_8_mbs': np.array([5, 13, 15, 16, 17, 18, 21, 29, 30, 31, 32, 34, 37, 38, 39, 41,
                          42, 43, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 62,
                          64, 65, 68, 73, 74, 75, 76, 85, 98, 99, 125], dtype=int)
}

# Notch intensity data
notch_data = np.array([
    78.882, 70.702, 62.452, 55.169, 49.589, 44.339, 40.863, 38.124,
    35.518, 33.322, 31.955, 31.507, 30.938, 30.339, 29.484, 29.244,
    29.084, 28.970, 28.560, 28.366, 28.093, 28.167, 27.262, 26.740,
    26.212, 26.267, 27.035, 26.509, 26.545, 25.619, 26.036, 26.059,
    25.479, 25.472, 25.191, 25.506, 24.864, 24.734, 24.919, 25.314,
    26.016, 25.677, 25.618, 25.282, 25.482, 25.690, 25.258, 25.165,
    25.084, 25.022, 25.243, 25.178, 25.359, 25.294, 25.470, 25.763,
    25.659, 25.891, 25.820, 25.924, 26.172, 26.732, 26.190, 26.571
])

# Data loading
height_list = [85., 170., 105., 105.9, 89., 81.5]
heights_dict = dict(zip(wing_regions, height_list))

# Dictionary creation
path_dict = {
    r: DATA_DIR / f'adjacency_matrices_{r}.xlsx'
    for r in wing_regions
}
sheets_dict = {
    r: pd.read_excel(path_dict[r], sheet_name=None, header=None)
    for r in wing_regions
}
A_dict = {
    r: [sheet.values.astype(float) for sheet in sheets_dict[r].values()]
    for r in wing_regions
}
centroids_dict = {
    r: load_centroids(r)
    for r in wing_regions
}
area_apical_dict = {
    r: load_areas(r)
    for r in wing_regions
}
diam_apical_dict = {
    r: 2 * np.sqrt(area_apical_dict[r] / np.pi)
    for r in wing_regions
}
signalling_labels_apical_dict = {
    region: [i for i in labels if A_dict[region][0][i].sum() > 0]
    for region, labels in signalling_labels_dict.items()
}
apical_neighbours_dict = {
    region: {
        i: [j for j in labels if A_dict[region][0][i, j] > 0]
        for i in labels
    }
    for region, labels in signalling_labels_dict.items()
}
nonapical_neighbours_dict = {
    region: {
        i: [j for j in labels if A_dict[region][0][i, j] == 0]
        for i in labels
    }
    for region, labels in signalling_labels_dict.items()
}

# Folder creation
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# Signalling function and heights
def omega_exp(z):
    return np.exp(-0.4 * z) + 0.5
def omega_cnt(z):
    return 1.5
def omega_lin(z):
    return -1.5/25 * z + 1.5
def omega_exp0(z):
    return np.exp(-0.1 * z) + 0.5
omega_map = {
    "exp":  omega_exp,
    "cnt":  omega_cnt,
    "lin":  omega_lin,
    "exp0": omega_exp0,
}
def height_set(diam_list, hs=0.04):
    return np.mean(diam_list[diam_list != 0]) * hs

