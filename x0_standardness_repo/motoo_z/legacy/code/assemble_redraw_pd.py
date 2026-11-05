#!/usr/bin/env python3
"""Extract an auditable crossing record from guided closed centre-line traces.

The source image, not a target invariant, determines over/under choices.
Unresolved choices remain explicit and prevent PD admission.  This tool does
not certify an isotopy from a different diagram.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage

from extract_diagrams import make_pd, linking_matrix, validate_with_libraries


COLOURS = ['red', 'salmon', 'green', 'blue', 'navy', 'yellow', 'purple', 'cyan']
PALETTE = {
    'salmon': (0.964706, 0.325490, 0.325490),
    'navy': (0.011765, 0.274495, 0.560769),
    'yellow': (0.819608, 0.827451, 0),
    'purple': (0.509804, 0.129412, 0.545098),
    'cyan': (0.090196, 0.850980, 1),
    'blue': (0, 0.478431, 1),
    'red': (0.823529, 0, 0),
    'green': (0.494102, 0.776471, 0.211749),
}
DOTTED = {'salmon', 'green', 'navy'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def save(p, value):
    Path(p).write_text(json.dumps(value, indent=2) + '\n')


def traces_from(files):
    result = {}
    for filename in files:
        value = json.loads(filename.read_text())
        parts = value.get('components', value.get('traces', value)) if isinstance(value, dict) else value
        if isinstance(parts, dict):
            parts = [dict(v, colour=k) if isinstance(v, dict) else {'colour': k, 'points_pdf': v}
                     for k, v in parts.items() if k in COLOURS]
        for c in parts:
            colour = c.get('colour', c.get('color'))
            if colour not in COLOURS:
                continue
            points = c.get('points_pdf', c.get('points'))
            if points is None:
                raise ValueError(f'No points for {colour}')
            p = []
            for q in points:
                q = list(map(float, q))
                if not p or q != p[-1]:
                    p.append(q)
            if len(p) < 4:
                raise ValueError(f'Trace {colour} is too short')
            if p[-1] != p[0]:
                p.append(p[0])
            result[colour] = dict(c, points_pdf=p, source_trace_file=str(filename),
                                  source_trace_sha256=sha(filename))
    return result


def det(a, b):
    return a[0]*b[1]-a[1]*b[0]


def intersections(traces):
    segments, grid = [], defaultdict(list)
    cell = 12.0
    lengths = {}
    for colour in COLOURS:
        if colour not in traces:
            continue
        points = traces[colour]['points_pdf']
        lengths[colour] = len(points)-1
        for n, (a, b) in enumerate(zip(points, points[1:])):
            v = [b[0]-a[0], b[1]-a[1]]
            norm = math.hypot(*v)
            if not norm:
                continue
            s = {'colour': colour, 'n': n, 'a': a, 'b': b,
                 'v': v, 'tangent': [x/norm for x in v]}
            i = len(segments)
            segments.append(s)
            for x in range(math.floor(min(a[0], b[0])/cell), math.floor(max(a[0], b[0])/cell)+1):
                for y in range(math.floor(min(a[1], b[1])/cell), math.floor(max(a[1], b[1])/cell)+1):
                    grid[x, y].append(i)
    checked, hits = set(), []
    for bucket in grid.values():
        for ii, i in enumerate(bucket):
            a = segments[i]
            for j in bucket[ii+1:]:
                key = tuple(sorted((i, j)))
                if key in checked:
                    continue
                checked.add(key)
                b = segments[j]
                if a['colour'] == b['colour']:
                    delta = abs(a['n']-b['n'])
                    if delta in (0, 1, lengths[a['colour']]-1):
                        continue
                denominator = det(a['v'], b['v'])
                if abs(denominator) < 1e-10:
                    continue
                offset = [b['a'][k]-a['a'][k] for k in (0, 1)]
                t, u = det(offset, b['v'])/denominator, det(offset, a['v'])/denominator
                # Half-open segments count a vertex intersection once.
                if not (-1e-9 <= t < 1-1e-9 and -1e-9 <= u < 1-1e-9):
                    continue
                point = [a['a'][k]+t*a['v'][k] for k in (0, 1)]
                visits = [{'colour': s['colour'], 'segment': s['n'], 'fraction': w,
                           'order': s['n']+w, 'tangent_pdf': s['tangent']}
                          for s, w in ((a, t), (b, u))]
                hits.append({'position_pdf': point, 'visits': visits,
                             'angle_sine': abs(det(a['tangent'], b['tangent']))})
    hits.sort(key=lambda c: (c['position_pdf'][1], c['position_pdf'][0]))
    for n, c in enumerate(hits):
        c['id'] = n
    return hits, len(segments)


class InkEvidence:
    def __init__(self, layers, scale):
        self.scale = scale
        self.layers = {}
        for colour in COLOURS:
            p = layers/f'layer_{colour}.png'
            a = np.asarray(Image.open(p).convert('RGB'))
            ink = np.max(255-a, axis=-1) > 20
            self.layers[colour] = ndimage.distance_transform_edt(~ink)/scale
        self.original = np.asarray(Image.open(layers/'original_page.png').convert('RGB'))/255
        self.height, self.width = self.original.shape[:2]

    def distance(self, colour, point):
        return float(ndimage.map_coordinates(self.layers[colour],
                     [[point[1]*self.scale], [point[0]*self.scale]], order=1, mode='nearest')[0])

    def visible_distance(self, colour, point, radius=8):
        """Distance to a source pixel with this colour (allow antialiasing)."""
        x, y = [v*self.scale for v in point]
        r = int(math.ceil(radius*self.scale))
        x0, x1 = max(0, int(x)-r), min(self.width, int(x)+r+1)
        y0, y1 = max(0, int(y)-r), min(self.height, int(y)+r+1)
        a = self.original[y0:y1, x0:x1]
        v = 1-np.array(PALETTE[colour])
        q = 1-a
        alpha = np.sum(q*v, axis=2)/np.sum(v*v)
        residual = np.linalg.norm(q-alpha[..., None]*v, axis=2)
        pixels = np.argwhere((alpha > .12) & (alpha < 1.1) & (residual < .065))
        if not len(pixels):
            return None
        return float(np.min(np.hypot(pixels[:, 0]+y0-y, pixels[:, 1]+x0-x))/self.scale)

    def decide(self, crossing):
        a, b = crossing['visits']
        p = crossing['position_pdf']
        colours = [a['colour'], b['colour']]
        evidence = {'isolated_layer_distance': [self.distance(k, p) for k in colours],
                    'original_visible_distance': [self.visible_distance(k, p) for k in colours]}
        neighbourhoods = []
        for v in (a, b):
            neighbourhoods.append([self.distance(v['colour'],
                [p[k]+offset*v['tangent_pdf'][k] for k in (0,1)])
                for offset in (-4.5,-3,-1.5,1.5,3,4.5)])
        evidence['branch_neighbourhood_ink_distances'] = neighbourhoods
        if colours[0] == colours[1]:
            means = [sum(row)/len(row) for row in neighbourhoods]
            best = 0 if means[0] < means[1] else 1
            if means[best] <= .75 and means[1-best]-means[best] > .8:
                return best, dict(evidence, method='self_crossing_branch_gap_requires_visual_audit')
            return None, dict(evidence, method='self_crossing_requires_trace_gap_or_review')
        ds = [10 if d is None else d for d in evidence['original_visible_distance']]
        best = 0 if ds[0] < ds[1] else 1
        if ds[best] <= 1.6 and ds[1-best]-ds[best] >= 1.0:
            return best, dict(evidence, method='source_visible_colour_gap')
        raw = evidence['isolated_layer_distance']
        best = 0 if raw[0] < raw[1] else 1
        if raw[best] <= 1.0 and raw[1-best]-raw[best] >= 2.0:
            return best, dict(evidence, method='isolated_layer_gap_with_source_review_pending')
        return None, dict(evidence, method='manual_review_required')


def annotate(source, traces, crossings, output, scale):
    im = Image.open(source).convert('RGB')
    d = ImageDraw.Draw(im)
    try:
        font = ImageFont.truetype('/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf', 15)
    except OSError:
        font = ImageFont.load_default()
    for c in crossings:
        x, y = [v*scale for v in c['position_pdf']]
        color = (200, 0, 0) if c['over_visit'] is None else (0, 0, 0)
        d.ellipse((x-4, y-4, x+4, y+4), outline=color, width=1)
        d.text((x+5, y-16), str(c['id']), fill=color, font=font, stroke_width=1, stroke_fill=(255,255,255))
    im.save(output/'crossings_numbered.png')
    overlay = Image.open(source).convert('RGB')
    d = ImageDraw.Draw(overlay)
    for colour, trace in traces.items():
        points = [tuple(v*scale for v in p) for p in trace['points_pdf']]
        d.line(points, fill=tuple(int(x*255) for x in PALETTE[colour]), width=2)
    overlay.save(output/'traces_overlay.png')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--traces', nargs='+', required=True, type=Path)
    p.add_argument('--layers', required=True, type=Path)
    p.add_argument('--figure', required=True, type=int, choices=(1,2))
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--scale', type=float, default=2)
    p.add_argument('--overrides', type=Path)
    p.add_argument('--link-python', type=Path)
    args = p.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    traces = traces_from(args.traces)
    cross, nsegments = intersections(traces)
    ink = InkEvidence(args.layers, args.scale)
    overrides = json.loads(args.overrides.read_text()) if args.overrides else {}
    overrides = overrides.get('crossings', overrides)
    for c in cross:
        over, evidence = ink.decide(c)
        choice = overrides.get(str(c['id']))
        if choice is not None:
            if not isinstance(choice, dict) or not choice.get('reason'):
                raise ValueError('Manual crossing choices require a reason')
            if math.dist(choice['position_pdf'], c['position_pdf']) > 1.0:
                raise ValueError('Manual crossing override no longer matches position')
            over = choice.get('over_visit')
            if over not in (0, 1):
                wanted = choice['over_colour']
                matches = [i for i,v in enumerate(c['visits']) if v['colour'] == wanted]
                if len(matches) != 1:
                    raise ValueError('Self-crossing override needs over_visit')
                over = matches[0]
            evidence['manual_override'] = choice
        c['over_visit'], c['evidence'] = over, evidence
    unknown = [c['id'] for c in cross if c['over_visit'] is None]
    report = {'schema': 'motoo.guided_redraw_crossings.v1', 'figure': args.figure,
              'source_layers': str(args.layers), 'source_original_png_sha256': sha(args.layers/'original_page.png'),
              'trace_files': [{'path':str(f),'sha256':sha(f)} for f in args.traces],
              'component_colours': list(traces), 'segments': nsegments,
              'crossing_count': len(cross), 'unresolved_crossings': unknown, 'crossings': cross,
              'status': 'CROSSING_RECORD_REQUIRES_VISUAL_AUDIT'}
    save(args.out/'crossing_record.json', report)
    annotate(args.layers/'original_page.png', traces, cross, args.out, args.scale)
    print(json.dumps({'components': len(traces), 'crossings': len(cross), 'unresolved': unknown}), flush=True)
    if unknown or set(traces) != set(COLOURS):
        return
    imported_crossings = []
    for c in cross:
        o, u = c['visits'][c['over_visit']], c['visits'][1-c['over_visit']]
        imported_crossings.append({'id':c['id'], 'position_pdf':c['position_pdf'],
            'over_component':COLOURS.index(o['colour']), 'under_component':COLOURS.index(u['colour']),
            'over_order': [o['order']], 'under_order':[u['order']],
            'over_tangent_pdf':o['tangent_pdf'], 'under_tangent_pdf':u['tangent_pdf'],
            'sign':1 if -det(o['tangent_pdf'],u['tangent_pdf']) > 0 else -1})
    pd, arcs, visits = make_pd(imported_crossings, len(COLOURS))
    matrix, writhe = linking_matrix(imported_crossings, len(COLOURS))
    components = [{'index':i, 'id':colour, 'name':colour, 'dotted':colour in DOTTED,
                   'handle_index':1 if colour in DOTTED else 2,
                   'framing':None if colour in DOTTED else (-1 if colour == 'blue' else 0),
                   'arc_labels':arcs[i], 'crossing_visits':visits[i], 'self_writhe':writhe[i]}
                  for i,colour in enumerate(COLOURS)]
    result = {'schema':'motoo.guided_redraw_link.v1', 'figure':args.figure,
              'pd_code':pd, 'crossing_count':len(pd), 'component_count':len(components),
              'components':components, 'crossings':imported_crossings, 'linking_matrix':matrix,
              'crossing_record_sha256':sha(args.out/'crossing_record.json'),
              'source_diagram_isotopy_independently_verified':False}
    save(args.out/'candidate_link.json', result)
    if args.link_python:
        validate_with_libraries(result, args.link_python)
        save(args.out/'native_checked_link.json', result)
        print(json.dumps({'native_PD_validation':result['pd_validation'], 'linking_matrix':matrix}), flush=True)


if __name__ == '__main__':
    main()
