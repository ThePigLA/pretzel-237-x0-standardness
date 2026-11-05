#!/usr/bin/env python3
"""Diagnostic skeletons and candidate gap joins; no PD or topology admission.

Uses NumPy, SciPy and Pillow from the analysis runtime. All skeleton branches,
small fragments and competing gap candidates are retained for human audit.
No candidate connection is applied automatically.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy import ndimage
from PIL import Image, ImageDraw


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def thin(mask):
    s = np.pad(mask.astype(np.uint8), 1)
    iterations = 0
    for iterations in range(1, 201):
        changed = 0
        for phase in (0, 1):
            p = [s[:-2, 1:-1], s[:-2, 2:], s[1:-1, 2:], s[2:, 2:],
                 s[2:, 1:-1], s[2:, :-2], s[1:-1, :-2], s[:-2, :-2]]
            count = sum(p)
            transitions = sum((p[i] == 0) & (p[(i + 1) % 8] == 1) for i in range(8))
            n, e, south, w = p[0], p[2], p[4], p[6]
            if phase == 0:
                corner = (n * e * south == 0) & (e * south * w == 0)
            else:
                corner = (n * e * w == 0) & (n * south * w == 0)
            delete = (s[1:-1, 1:-1] != 0) & (count >= 2) & (count <= 6) & (transitions == 1) & corner
            changed += int(delete.sum())
            s[1:-1, 1:-1][delete] = 0
        if not changed:
            break
    if iterations == 200 and changed:
        raise ValueError('Skeleton thinning did not converge')
    return s[1:-1, 1:-1].astype(bool), iterations


def make_graph(skeleton):
    points = {tuple(map(int, p)) for p in np.argwhere(skeleton)}
    graph = {}
    for y, x in points:
        adj = []
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if not (dx or dy) or (y + dy, x + dx) not in points:
                    continue
                # A diagonal redundant with an orthogonal two-edge route is
                # omitted; otherwise a digital right-angle creates a triangle.
                if dx and dy and ((y + dy, x) in points or (y, x + dx) in points):
                    continue
                adj.append((y + dy, x + dx))
        graph[y, x] = sorted(adj)
    return graph


def chains(graph):
    seen = set()
    result = []
    def walk(a, b):
        path = [a, b]
        seen.add(tuple(sorted((a, b))))
        while len(graph[b]) == 2:
            c = next(p for p in graph[b] if p != a)
            edge = tuple(sorted((b, c)))
            if edge in seen:
                break
            seen.add(edge)
            path.append(c)
            a, b = b, c
        result.append(path)
    for a in sorted(graph):
        if len(graph[a]) != 2:
            for b in graph[a]:
                if tuple(sorted((a, b))) not in seen:
                    walk(a, b)
    for a in sorted(graph):
        for b in graph[a]:
            if tuple(sorted((a, b))) not in seen:
                walk(a, b)
    return result


def analyse(path, output, name, scale, threshold):
    started = time.monotonic()
    rgb = np.asarray(Image.open(path).convert('RGB'))
    ink = np.max(255 - rgb, axis=-1) >= threshold
    _, before_components = ndimage.label(ink, structure=np.ones((3, 3)))
    skel, iterations = thin(ink)
    graph = make_graph(skel)
    traces = chains(graph)
    fragments, endpoints = [], []
    for i, points in enumerate(traces):
        length = sum(math.dist(a, b) for a, b in zip(points, points[1:])) / scale
        fragments.append({'id': i, 'length_points': length,
                          'coordinates_pdf': [[x / scale, y / scale] for y, x in points],
                          'endpoint_degrees': [len(graph[points[0]]), len(graph[points[-1]])],
                          'tiny_fragment_under_2_points': length < 2})
        for end in (0, 1):
            q = points if end == 0 else points[::-1]
            if len(graph[q[0]]) != 1:
                continue
            # This is a directional estimate, not a certified tangent.
            k = min(len(q) - 1, max(1, int(6 * scale)))
            dy, dx = q[0][0] - q[k][0], q[0][1] - q[k][1]
            norm = math.hypot(dx, dy)
            if not norm:
                continue
            endpoints.append({'id': len(endpoints), 'fragment': i, 'end': end,
                              'position_pdf': [q[0][1] / scale, q[0][0] / scale],
                              'outward_direction_estimate': [dx / norm, dy / norm]})
    candidates = []
    for i, a in enumerate(endpoints):
        for b in endpoints[i + 1:]:
            if a['fragment'] == b['fragment']:
                continue
            ax, ay = a['position_pdf']; bx, by = b['position_pdf']
            d = math.hypot(bx - ax, by - ay)
            if d == 0 or d > 24:
                continue
            direction = ((bx - ax) / d, (by - ay) / d)
            align_a = sum(u*v for u, v in zip(a['outward_direction_estimate'], direction))
            align_b = -sum(u*v for u, v in zip(b['outward_direction_estimate'], direction))
            if min(align_a, align_b) < 0.5:
                continue
            candidates.append({'endpoints': [a['id'], b['id']], 'distance_points': d,
                               'endpoint_alignments': [align_a, align_b],
                               'heuristic_score': d + 20*(2-align_a-align_b)})
    candidates.sort(key=lambda c: c['heuristic_score'])
    by_endpoint = {e['id']: [] for e in endpoints}
    for c in candidates:
        for e in c['endpoints']:
            by_endpoint[e].append(c)
    for e in endpoints:
        e['candidate_count'] = len(by_endpoint[e['id']])
    report = {'status': 'DIAGNOSTIC_ONLY_NO_CONNECTIONS_APPLIED',
              'source_layer': str(path), 'source_sha256': sha(path),
              'colour': name, 'scale_pixels_per_point': scale,
              'ink_threshold': threshold, 'ink_pixels': int(ink.sum()),
              'ink_components_eight_neighbour': int(before_components),
              'skeleton_pixels': int(skel.sum()), 'thinning_iterations': iterations,
              'skeleton_degree_histogram': dict(Counter(len(v) for v in graph.values())),
              'junction_pixels_pdf': [[x/scale, y/scale] for (y, x), n in graph.items() if len(n) > 2],
              'fragments': fragments, 'endpoints': endpoints, 'gap_candidates': candidates,
              'gap_candidates_applied': False, 'pd_generated': False,
              'crossing_over_under_assigned': False,
              'annotations_and_incidental_marks_not_removed': True,
              'wall_seconds': time.monotonic() - started}
    (output / f'{name}_fragments.json').write_text(json.dumps(report, separators=(',', ':')) + '\n')
    canvas = Image.fromarray(rgb.copy());draw = ImageDraw.Draw(canvas)
    for e in endpoints:
        x, y = [v*scale for v in e['position_pdf']]
        draw.ellipse((x-4,y-4,x+4,y+4), outline=(255,0,100), width=1)
        draw.text((x+4,y-5), str(e['id']), fill=(0,0,0))
    for y, x in graph:
        if len(graph[y, x]) > 2:
            draw.rectangle((x-2,y-2,x+2,y+2), fill=(0,0,0))
    canvas.save(output / f'{name}_endpoints.png')
    return {k:report[k] for k in ('colour','ink_components_eight_neighbour','skeleton_pixels',
                                   'skeleton_degree_histogram','wall_seconds')} | {
        'fragments':len(fragments), 'endpoints':len(endpoints), 'gap_candidates':len(candidates),
        'endpoints_without_candidate':sum(not e['candidate_count'] for e in endpoints),
        'endpoints_with_competing_candidates':sum(e['candidate_count']>1 for e in endpoints)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('layers',type=Path)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--scale',type=float,default=2)
    parser.add_argument('--threshold',type=int,default=20)
    parser.add_argument('--colours',nargs='+',default=['blue','cyan','purple','red','yellow'])
    args=parser.parse_args()
    if args.out.exists():raise ValueError('Use a new output directory')
    args.out.mkdir(parents=True)
    rows=[]
    for colour in args.colours:
        row=analyse(args.layers/f'layer_{colour}.png',args.out,colour,args.scale,args.threshold)
        rows.append(row)
        print(json.dumps(row),flush=True)
    report={'status':'DIAGNOSTIC_ONLY_PD_TRANSCRIPTION_PENDING','layers':rows,
            'no_gap_connections_applied':True,'script_sha256':sha(__file__)}
    (args.out/'CENTERLINE_DIAGNOSTIC.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
