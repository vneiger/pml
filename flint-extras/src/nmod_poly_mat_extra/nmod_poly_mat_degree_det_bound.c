/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of PML.

    PML is free software: you can redistribute it and/or modify it under
    the terms of the GNU General Public License version 2.0 (GPL-2.0-or-later)
    as published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version. See
    <https://www.gnu.org/licenses/>.
*/

#include <flint/flint.h>
#include <flint/nmod_poly.h>
#include <flint/nmod_poly_mat.h>

#include "nmod_poly_mat_forms.h"

/* marks forbidden pairs and unreached jobs; only used in comparisons, never
 * in arithmetic */
#define _ASSIGN_INF WORD_MAX

/* Maximum weight assignment of k >= 1 "workers" to pairwise distinct "jobs"
 * among N >= k, with weights w[a*N + b] (worker a, job b), where
 * _ASSIGN_INF marks a forbidden pair. Returns 0 if every assignment of all
 * workers uses a forbidden pair. Otherwise returns 1, sets total to the
 * maximum total weight, and if assign is not NULL, sets assign[a] to the job
 * of worker a in an optimal assignment.
 *
 * Successive shortest augmenting paths with node potentials (Edmonds-Karp
 * 1972, Tomizawa 1971), on the costs c = -w: workers are inserted one at a
 * time; for each of them, a Dijkstra-like search on the reduced costs
 * c(a,b) - u[a] - v[b] (nonnegative, zero on matched pairs) finds a
 * shortest alternating path to a free job, then the potentials are updated
 * so that reduced costs stay nonnegative, and the matching is augmented
 * along the path. Each search costs O(k N), hence O(k^2 N) in total.
 *
 * Size of the integers. Let B be the maximum of |w| over allowed pairs. At
 * all times, v[b] = 0 if job b is free and v[b] <= 0 otherwise, and a worker
 * a matched to job b has u[a] = c(a,b) - v[b]; a worker not inserted yet has
 * its initial potential, in [-B, B]. In the search from worker a0, the final
 * distance of a job b is dist[b] = L(b) - u[a0] - v[b], where L(b) is the
 * minimum, over alternating paths from a0 to b, of the cost of their
 * unmatched pairs minus the cost of their matched pairs; such a path has at
 * most 2k-1 pairs, so |L(b)| <= (2k-1) B. The update sets v[b] to L(b) -
 * L(end) for the fixed jobs b, hence -(4k-2) B <= v[b] <= 0 at all times, and
 * -B <= u[a] <= (4k-1) B. The final distances are in [0, D], where D =
 * dist[end] = L(end) - u[a0] <= 2k B. As a result, all quantities computed
 * below, including the partial sums in the tentative distances da + c(a,b) -
 * u[a] - v[b], are at most 6k B in absolute value. Thus, if B <= WORD_MAX /
 * (8k) (which the documentation requires from the caller), there is no
 * overflow and all these quantities are smaller than _ASSIGN_INF. */
static int _max_weight_assignment(slong * total, slong * assign,
                                  const slong * w, slong k, slong N)
{
    slong * u = flint_malloc(k * sizeof(slong));        /* worker potentials */
    slong * v = flint_calloc(N, sizeof(slong));         /* job potentials */
    slong * job_of = flint_malloc(k * sizeof(slong));   /* matched job, or -1 */
    slong * worker_of = flint_malloc(N * sizeof(slong));/* matched worker, or -1 */
    slong * dist = flint_malloc(N * sizeof(slong));     /* tentative distances */
    slong * pred = flint_malloc(N * sizeof(slong));     /* worker preceding job */
    char * fixed = flint_malloc(N);                     /* distance is final */
    int ok = 1;

    for (slong b = 0; b < N; b++)
        worker_of[b] = -1;

    /* initial potentials: u[a] = min_b c(a,b) = -max_b w(a,b), so that all
     * reduced costs are nonnegative */
    for (slong a = 0; ok && a < k; a++)
    {
        slong best = _ASSIGN_INF;
        for (slong b = 0; b < N; b++)
            if (w[a*N + b] != _ASSIGN_INF && -w[a*N + b] < best)
                best = -w[a*N + b];
        if (best == _ASSIGN_INF)
            ok = 0;  /* worker without allowed job */
        u[a] = best;
        job_of[a] = -1;
    }

    for (slong a0 = 0; ok && a0 < k; a0++)
    {
        for (slong b = 0; b < N; b++)
        {
            dist[b] = _ASSIGN_INF;
            fixed[b] = 0;
        }

        /* Dijkstra from worker a0; da is the (final) distance of worker a */
        slong a = a0, da = 0, end = -1;
        while (1)
        {
            for (slong b = 0; b < N; b++)
            {
                if (fixed[b] || w[a*N + b] == _ASSIGN_INF)
                    continue;
                const slong nd = da - w[a*N + b] - u[a] - v[b];
                if (nd < dist[b])
                {
                    dist[b] = nd;
                    pred[b] = a;
                }
            }
            slong bmin = -1, dmin = _ASSIGN_INF;
            for (slong b = 0; b < N; b++)
                if (!fixed[b] && dist[b] < dmin)
                {
                    dmin = dist[b];
                    bmin = b;
                }
            if (bmin == -1)
            {
                ok = 0;  /* no augmenting path: no complete assignment */
                break;
            }
            fixed[bmin] = 1;
            if (worker_of[bmin] == -1)
            {
                end = bmin;
                break;
            }
            /* the matched pair has reduced cost zero */
            a = worker_of[bmin];
            da = dmin;
        }
        if (!ok)
            break;

        /* potentials: shift by the distances, relative to D = dist[end] */
        const slong D = dist[end];
        u[a0] += D;
        for (slong b = 0; b < N; b++)
            if (fixed[b] && b != end)
            {
                u[worker_of[b]] += D - dist[b];
                v[b] -= D - dist[b];
            }

        /* augment along the path ending at job `end` */
        for (slong b = end; ; )
        {
            const slong pa = pred[b];
            const slong prev_job = job_of[pa];
            job_of[pa] = b;
            worker_of[b] = pa;
            if (pa == a0)
                break;
            b = prev_job;
        }
    }

    if (ok)
    {
        *total = 0;
        for (slong a = 0; a < k; a++)
            *total += w[a*N + job_of[a]];
        if (assign)
            for (slong a = 0; a < k; a++)
                assign[a] = job_of[a];
    }

    flint_free(fixed);
    flint_free(pred);
    flint_free(dist);
    flint_free(worker_of);
    flint_free(job_of);
    flint_free(v);
    flint_free(u);
    return ok;
}

slong nmod_poly_mat_degree_det_bound(slong * assign,
                                     const nmod_poly_mat_t pmat,
                                     const slong * shift,
                                     orientation_t orient)
{
    const slong m = pmat->r;
    const slong n = pmat->c;
    const int tall = (m >= n);        /* workers: columns if tall, else rows */
    const slong k = tall ? n : m;
    const slong N = tall ? m : n;
    const int rowwise = (orient == ROW_LOWER || orient == ROW_UPPER);
    const slong slen = rowwise ? n : m;

    if (k == 0)
        return 0;

    /* min(shift), for the value returned when there is no assignment */
    slong min_shift = 0;
    if (shift)
    {
        min_shift = shift[0];
        for (slong i = 1; i < slen; i++)
            min_shift = FLINT_MIN(min_shift, shift[i]);
    }

    slong * w = flint_malloc(k * N * sizeof(slong));
    for (slong i = 0; i < m; i++)
        for (slong j = 0; j < n; j++)
        {
            const slong d = nmod_poly_degree(nmod_poly_mat_entry(pmat, i, j));
            slong wij = _ASSIGN_INF;
            if (d >= 0)
                wij = d + (shift ? (rowwise ? shift[j] : shift[i]) : 0);
            if (tall)
                w[j*N + i] = wij;
            else
                w[i*N + j] = wij;
        }

    slong res;
    if (! _max_weight_assignment(&res, assign, w, k, N))
        res = k * min_shift - 1;  /* below the weight of any assignment */

    flint_free(w);
    return res;
}

#undef _ASSIGN_INF
