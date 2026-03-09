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
import re
from functools import lru_cache
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
def compute_omega_k(
    omega_func,
    wing_region,
    Lmax=None,
    method='simpson',
    num_simpson=101,
    gap=None,
    n_layers=None,
):

    if gap is None:
        gap = get_region_gap(wing_region)
    if n_layers is None:
        n_layers = get_region_n_layers(wing_region)
    n = int(n_layers)

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

# Load apical area data

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
                          signalling_labels=None,
                          adjacency_layers=None,
                          centroid_layers=None,
                          **quad_kwargs):

    if weight_list is None:
        if omega_func is None:
            raise ValueError("Either weight_list or omega_func must be provided")
        if normalQ:
            weight_list_base = compute_omega_k(
                omega_func, wing_region, Lmax=Lmax, method=quad_method, **quad_kwargs
            )
            weight_list_all = compute_omega_k(
                omega_func, wing_region, Lmax=32, method=quad_method, **quad_kwargs
            )
            weight_list = np.array([
                i * sum(weight_list_base) / sum(weight_list_all)
                for i in weight_list_base
            ])
        else:
            weight_list = compute_omega_k(
                omega_func, wing_region, Lmax=Lmax, method=quad_method, **quad_kwargs
            )

    if adjacency_layers is None:
        adjacency_layers = get_adjacency_layers(wing_region)
    else:
        adjacency_layers = [np.asarray(a, dtype=float) for a in adjacency_layers]

    if centroid_layers is None:
        centroid_layers = get_centroids_layers(wing_region)
    else:
        centroid_layers = [np.asarray(c, dtype=float) for c in centroid_layers]

    if signalling_labels is None:
        labels = get_signalling_labels(wing_region)
    else:
        labels = np.asarray(signalling_labels, dtype=int)

    thresholds = [0.1, 0.2, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
    results = {th: [] for th in thresholds}
    results_v = {th: [] for th in thresholds}

    if str_type == 'centroid':
        A_list = straight_adjacency(
            {wing_region: adjacency_layers},
            {wing_region: centroid_layers},
            alpha,
        )[wing_region]
    else:
        A_list = adjacency_layers

    cents = np.array(centroid_layers[0])
    adj0 = adjacency_layers[0]

    sop_pairs = []
    bimodal_list = []

    for _ in range(sim_number):
        wk = weight_list

        if str_type == 'prune':
            if 'straight_prune' not in globals():
                raise NameError(
                    "str_type='prune' requires a straight_prune(...) helper, "
                    "which is not defined in msm_model_test.py"
                )
            straight_data = straight_prune(
                {wing_region: adjacency_layers},
                {wing_region: centroid_layers},
                alpha,
            )
            A_list, prune_idx = straight_data[0][wing_region], straight_data[1][wing_region]
            wk_new = [0 if i in prune_idx else val for i, val in enumerate(wk)]
            wk_new = [v * sum(wk) / sum(wk_new) for v in wk_new]
            wk = wk_new

        if epsmodelQ:
            dist_min = min_long_range_distance(adjacency_layers[0], centroid_layers[0], labels)
            M = long_range_M(
                adjacency_layers[0],
                centroid_layers[0],
                (1 + prot_len) * (dist_min - 1e-6),
                eps,
                labels,
            )
        else:
            M = sum(w * A for w, A in zip(wk, A_list))

        notch, delta = simulate(
            M, t_final, dt, labels, k=k, h=h, Ka=Ka, Kr=Kr, nu=nu, randomQ=randomQ
        )

        bimodal_list.append(SOP_bimodality(delta[labels], notch[labels]))

        for th in thresholds:
            sop_idx = np.where(delta > th)[0]
            d, v = band_avg_min_sop_distance(cents, sop_idx, adj0, height, y_shift_steps)
            results[th].append(d)
            results_v[th].append(v)

        if plotQ:
            marker = notch if marker_type == 'notch' else delta
            plot_layer_graph(
                marker,
                centroid_layers,
                adjacency_layers,
                wing_region=wing_region,
                framenumber=0,
                saveQ=graphsaveQ,
                show_labels=show_labels,
                show_other_layers=show_other_layers,
            )

        sop_pairs.append(np.triu(adj0[np.ix_(sop_idx, sop_idx)] > 0, k=1).sum())

    avg_results = {}
    for th, ds in results.items():
        finite = [x for x in ds if np.isfinite(x)]
        avg_results[th] = np.mean(finite) if finite else np.inf
        if th == 0.1:
            ds0 = ds

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
    plt.title(f'WD {get_region_label(region)}')
    if saveQ:
        plt.savefig(FIGURES_DIR / f'straight_neighbour_{region}.pdf', bbox_inches='tight', transparent=True)
    plt.show()

    return deg_nonap_list

# Sensitivity analysis functions
    
def SOP_bimodality(delta, notch):
    corr_score = -np.corrcoef(delta, notch)[0,1]
    delta_bin = delta > np.median(delta)
    notch_bin = notch > np.median(notch)
    xor_score = np.mean(delta_bin ^ notch_bin)
    return corr_score
    





# Distance heights



### Parameter setup and data loading ###

# Notch-Delta parameters
k, h = 2, 2
Ka, Kr = 0.1, 0.001
nu = 1.
t_final = 1000.
dt = 0.1

# Region-specific parameters are externalized to data files.
# This makes the model reusable for arbitrary datasets.

_REGION_METADATA_PATH = DATA_DIR / "wing_region_metadata.csv"
SIGNALLING_LABELS_DIR = DATA_DIR / "signalling_labels"  # legacy fallback only


def _coerce_region(name):
    return str(name).strip()


def _parse_signalling_labels_cell(value):
    if value is None:
        return np.array([], dtype=int)
    s = str(value).strip()
    if not s or s.lower() == 'nan':
        return np.array([], dtype=int)
    if s.startswith('[') and s.endswith(']'):
        vals = json.loads(s)
        return np.asarray(vals, dtype=int)
    parts = [p for p in re.split(r'[;,\s]+', s) if p]
    return np.asarray([int(p) for p in parts], dtype=int)


@lru_cache(maxsize=1)
def _load_region_metadata_table():
    if _REGION_METADATA_PATH.exists():
        df = pd.read_csv(_REGION_METADATA_PATH)
        if "wing_region" not in df.columns:
            raise ValueError("wing_region_metadata.csv must include a 'wing_region' column")
        df["wing_region"] = df["wing_region"].astype(str)
        return df.set_index("wing_region", drop=False)
    cols = ["wing_region", "wing_label", "gap", "n_layers", "default_height", "is_wing_disc"]
    return pd.DataFrame(columns=cols).set_index("wing_region", drop=False)


def list_wing_regions():
    names = set()
    for p in (DATA_DIR / "adjacency_matrices").glob("adjacency_matrices_*.xlsx"):
        names.add(p.stem.replace("adjacency_matrices_", ""))
    for p in (DATA_DIR / "cell_geometry").glob("cell_geometry_*.csv"):
        names.add(p.stem.replace("cell_geometry_", ""))
    if not names:
        names = set(_load_region_metadata_table().index.tolist())
    return sorted(names)


def list_wing_discs(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    meta = _load_region_metadata_table()
    out = []
    for region in wing_regions:
        region = _coerce_region(region)
        if region in meta.index and "is_wing_disc" in meta.columns:
            flag = meta.loc[region, "is_wing_disc"]
            if str(flag).strip().lower() in {"1", "true", "yes"}:
                out.append(region)
        elif region.startswith("wd_") and "mbs" not in region:
            out.append(region)
    return out


def get_region_metadata(wing_region):
    wing_region = _coerce_region(wing_region)
    meta = _load_region_metadata_table()
    if wing_region not in meta.index:
        raise KeyError(
            f"Region '{wing_region}' missing from metadata. "
            f"Add it to {_REGION_METADATA_PATH}."
        )
    row = meta.loc[wing_region]
    return row.to_dict()


def get_region_gap(wing_region):
    row = get_region_metadata(wing_region)
    if pd.isna(row.get("gap")):
        raise ValueError(f"Missing gap for region '{wing_region}'")
    return float(row["gap"])


def get_region_n_layers(wing_region):
    row = get_region_metadata(wing_region)
    if not pd.isna(row.get("n_layers")):
        return int(row["n_layers"])
    return len(get_adjacency_layers(wing_region))


def get_region_height(wing_region):
    row = get_region_metadata(wing_region)
    if not pd.isna(row.get("default_height")):
        return float(row["default_height"])
    diam = get_diam_apical(wing_region)
    return float(height_set(diam, hs=0.04))


def get_region_label(wing_region):
    row = get_region_metadata(wing_region)
    label = row.get("wing_label")
    return wing_region if pd.isna(label) else str(label)


def build_gap_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_region_gap(r) for r in wing_regions}


def build_n_layers_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_region_n_layers(r) for r in wing_regions}


def build_default_height_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_region_height(r) for r in wing_regions}


def build_wd_label_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_region_label(r) for r in wing_regions}


@lru_cache(maxsize=None)
def _get_adjacency_layers_cached(wing_region):
    wing_region = _coerce_region(wing_region)
    path = DATA_DIR / "adjacency_matrices" / f"adjacency_matrices_{wing_region}.xlsx"
    sheets = pd.read_excel(path, sheet_name=None, header=None)
    return tuple(sheet.values.astype(float) for sheet in sheets.values())


def get_adjacency_layers(wing_region):
    return [arr.copy() for arr in _get_adjacency_layers_cached(_coerce_region(wing_region))]


def build_adjacency_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_adjacency_layers(r) for r in wing_regions}


@lru_cache(maxsize=None)
def _get_centroids_layers_cached(wing_region):
    wing_region = _coerce_region(wing_region)
    path = DATA_DIR / "cell_geometry" / f"cell_geometry_{wing_region}.csv"
    df = pd.read_csv(path)
    max_frame = int(df["frame"].max())
    out = []
    for frame in range(max_frame + 1):
        subdf = df[df["frame"] == frame]
        out.append(subdf[["centroid_x", "centroid_y"]].to_numpy(dtype=float))
    return tuple(out)


def get_centroids_layers(wing_region):
    return [arr.copy() for arr in _get_centroids_layers_cached(_coerce_region(wing_region))]


def build_centroids_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_centroids_layers(r) for r in wing_regions}


@lru_cache(maxsize=None)
def _get_area_apical_cached(wing_region):
    wing_region = _coerce_region(wing_region)
    path = DATA_DIR / "cell_geometry" / f"cell_geometry_{wing_region}.csv"
    df = pd.read_csv(path)
    subdf = df[df["frame"] == 0]
    return subdf["area"].to_numpy(dtype=float)


def get_area_apical(wing_region):
    return _get_area_apical_cached(_coerce_region(wing_region)).copy()


def build_area_apical_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_area_apical(r) for r in wing_regions}


def get_diam_apical(wing_region):
    area = get_area_apical(wing_region)
    return 2 * np.sqrt(area / np.pi)


def build_diam_apical_dict(wing_regions=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    return {r: get_diam_apical(r) for r in wing_regions}


def get_signalling_labels(wing_region, labels_dir=SIGNALLING_LABELS_DIR):
    wing_region = _coerce_region(wing_region)

    # Primary source: signalling_labels column in wing_region_metadata.csv
    meta = _load_region_metadata_table()
    if wing_region in meta.index and "signalling_labels" in meta.columns:
        parsed = _parse_signalling_labels_cell(meta.loc[wing_region, "signalling_labels"])
        if parsed.size > 0:
            return parsed

    # Legacy fallback: data/signalling_labels/<wing_region>.json
    path = labels_dir / f"{wing_region}.json"
    if path.exists():
        vals = json.loads(path.read_text(encoding="utf-8"))
        return np.asarray(vals, dtype=int)

    raise FileNotFoundError(
        f"No signalling labels found for '{wing_region}'. "
        "Expected wing_region_metadata.csv column 'signalling_labels' or a legacy JSON file."
    )


def load_signalling_labels_dict(wing_regions=None, labels_dir=SIGNALLING_LABELS_DIR):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    out = {}
    for r in wing_regions:
        out[r] = get_signalling_labels(r, labels_dir=labels_dir)
    return out




def get_signalling_labels_apical(wing_region, signalling_labels=None, adjacency_layers=None):
    if signalling_labels is None:
        signalling_labels = get_signalling_labels(wing_region)
    if adjacency_layers is None:
        adjacency_layers = get_adjacency_layers(wing_region)
    adj0 = adjacency_layers[0]
    return [int(i) for i in signalling_labels if adj0[int(i)].sum() > 0]


def get_apical_neighbours(wing_region, signalling_labels=None, adjacency_layers=None):
    if signalling_labels is None:
        signalling_labels = get_signalling_labels(wing_region)
    if adjacency_layers is None:
        adjacency_layers = get_adjacency_layers(wing_region)
    adj0 = adjacency_layers[0]
    out = {}
    for i in signalling_labels:
        i = int(i)
        out[i] = [int(j) for j in signalling_labels if adj0[i, int(j)] > 0]
    return out


def get_nonapical_neighbours(wing_region, signalling_labels=None, adjacency_layers=None):
    if signalling_labels is None:
        signalling_labels = get_signalling_labels(wing_region)
    if adjacency_layers is None:
        adjacency_layers = get_adjacency_layers(wing_region)
    adj0 = adjacency_layers[0]
    out = {}
    for i in signalling_labels:
        i = int(i)
        out[i] = [int(j) for j in signalling_labels if adj0[i, int(j)] == 0]
    return out


def build_signalling_labels_apical_dict(wing_regions=None, signalling_labels_dict=None, adjacency_dict=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    if signalling_labels_dict is None:
        signalling_labels_dict = load_signalling_labels_dict(wing_regions)
    out = {}
    for r in wing_regions:
        adj_layers = adjacency_dict[r] if adjacency_dict and r in adjacency_dict else None
        out[r] = get_signalling_labels_apical(r, signalling_labels=signalling_labels_dict[r], adjacency_layers=adj_layers)
    return out


def build_apical_neighbours_dict(wing_regions=None, signalling_labels_dict=None, adjacency_dict=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    if signalling_labels_dict is None:
        signalling_labels_dict = load_signalling_labels_dict(wing_regions)
    out = {}
    for r in wing_regions:
        adj_layers = adjacency_dict[r] if adjacency_dict and r in adjacency_dict else None
        out[r] = get_apical_neighbours(r, signalling_labels=signalling_labels_dict[r], adjacency_layers=adj_layers)
    return out


def build_nonapical_neighbours_dict(wing_regions=None, signalling_labels_dict=None, adjacency_dict=None):
    if wing_regions is None:
        wing_regions = list_wing_regions()
    if signalling_labels_dict is None:
        signalling_labels_dict = load_signalling_labels_dict(wing_regions)
    out = {}
    for r in wing_regions:
        adj_layers = adjacency_dict[r] if adjacency_dict and r in adjacency_dict else None
        out[r] = get_nonapical_neighbours(r, signalling_labels=signalling_labels_dict[r], adjacency_layers=adj_layers)
    return out


# Signalling functions
FIGURES_DIR.mkdir(parents=True, exist_ok=True)


def omega_exp(z):
    return np.exp(-0.4 * z) + 0.5


def omega_cnt(z):
    return 1.5


def omega_lin(z):
    return -1.5 / 25 * z + 1.5


def omega_exp0(z):
    return np.exp(-0.1 * z) + 0.5


omega_map = {
    "exp": omega_exp,
    "cnt": omega_cnt,
    "lin": omega_lin,
    "exp0": omega_exp0,
}


def height_set(diam_list, hs=0.04):
    diam_list = np.asarray(diam_list, dtype=float)
    valid = diam_list[diam_list != 0]
    return np.mean(valid) * hs
