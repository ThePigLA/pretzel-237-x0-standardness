# Certified unmarked contraction prefix

The closed, source-oriented control has been reduced from **103,744 to
44,728 pentachora** by **3,307 simplicial edge contractions**.  An independent
standard-library replayer passed every step.  The marked master was only read;
it was not changed.

This endpoint is an **unmarked oriented PL descendant**.  No deck action,
ramification sphere, twin sphere, puncture, peripheral chart, or other marking
is transferred to the endpoint by this certificate.

## Files and exact bindings

| File | Role |
| --- | --- |
| `source_facets.json.gz` | Exact source facets and orientations, bound to the marked master and its native admission |
| `certificate.json.gz` | Ordered 3,307-step contraction trace, endpoints, local hashes, and theorem premise |
| `endpoint.json.gz` | Exact endpoint facets and orientations, retaining the original vertex labels |
| `endpoint_regina.json.gz` | Native four-dimensional simplex gluing data generated independently from the replayed endpoint |
| `producer_summary.json` | Producer counts, resource use, and interpretation of its stopping event |
| `independent_replay.json` | Independent replay receipt and exact file hashes |

Source facet SHA-256:

`65de21dd16c5c11fa880e6f0c80a3cbf25a4780ba3d04ba7a6ab838ef37e5587`

Endpoint facet SHA-256:

`496c7a5cb6f5e8526db056bbd79e293252b181d86cbe0d7d7889681ceecc4d98`

Native endpoint file SHA-256:

`4155432efdffb657934814e7faa5c3ad4902f304357005dfef7787d56d459f24`

The endpoint has

\[
f=(1003,25367,91458,111820,44728),\qquad \chi=2.
\]

The replay checks connectedness, closed codimension-one incidence, the entire
oriented top-dimensional boundary chain, and the native gluing permutations.

## Why every step preserves PL type

The source premise is a closed PL four-manifold.  It is supplied by the
source-bound native admission of the immutable control, including recognition
of its vertex links as PL three-spheres.  A homology-sphere test would not be a
substitute for that premise.

For an edge \(uv\), the certificate requires the full simplicial link condition

\[
\operatorname{lk}(u)\cap\operatorname{lk}(v)
=\operatorname{lk}(uv)
\]

in every face dimension.  The empty face is automatic because this is an actual
edge with nonempty incident top cells; all nonempty faces are checked explicitly.

Put \(L=\operatorname{lk}(u)\).  Since the current complex is a closed PL
four-manifold, \(L\) is a PL three-sphere.  The closed antistar
\(A=\operatorname{ast}(v,L)\) is a PL three-ball.  The contraction replaces

\[
u*L \quad\text{by}\quad v*A.
\]

Both are PL four-balls, and both have boundary exactly \(L\).  The full link
condition excludes extra identifications with the retained complex.  Replacing
one ball by the other along the identical boundary is therefore a PL
homeomorphism relative to the retained complex.  The exact equality of the
oriented boundary chains fixes the orientation of this replacement.  Induction
preserves the closed PL-manifold premise and the global orientation at every
subsequent step.

This is the local proof in Eran Nevo, *Higher Minors and Van Kampen's
Obstruction*, Theorem 1.4 and Section 5.1,
[arXiv:math/0602531](https://arxiv.org/pdf/math/0602531),
published in *Mathematica Scandinavica* 101(2), 161–176 (2007),
[doi:10.7146/math.scand.a-15037](https://doi.org/10.7146/math.scand.a-15037).
The antistar assertion is a fact about a PL three-sphere; no four-dimensional
Schoenflies assertion is needed.

## Independence and scope of the finite checks

The producer is `diagram/contract_unmarked4.py`.  The replayer is
`diagram/replay_unmarked4.py`; it does not import the producer or reuse the
producer's link cache.  For every recorded step it reconstructs the full links,
the complete vertex star, the antistar cone, and all faces of the new ball.  It
checks the old and new oriented boundary chains are identical and checks every
new face against the actual retained global complex.  The replay performed
3,307 complete link checks and 1,190,001 full-face overlap checks.

It then compares every endpoint facet and orientation sign with the supplied
endpoint and constructs native gluing data by matching literal codimension-one
faces.  The native endpoint is suitable for subsequent invariant calculations
or a separate simplification search.  Any later PL standardization must include
this contraction prefix in its ancestry and independently certify its remaining
steps.

Neither simple connectivity nor standardness of the endpoint is asserted by
this contraction replay.  Those are separate mathematical checks.

## Budget and stopping event

| Stage | Wall seconds | CPU seconds | Peak resident bytes |
| --- | ---: | ---: | ---: |
| Producer | 42.6285 | 42.6262 | 349,491,200 |
| Independent replay | 42.3475 | 42.3256 | 399,519,744 |

The producer was capped at 55 seconds for search/preparation, with output time
reserved within the requested 60-second first checkpoint.  Both processes
enforced an 8 GiB address-space limit.  Each contraction strictly reduced the
number of pentachora; it never exceeded the original 103,744.

The producer stopped when its local priority queue was exhausted after 272,938
candidate attempts.  The original certificate's legacy string
`no_remaining_admissible_candidate` denotes this queue event.  It does **not**
certify that every remaining edge is noncontractible: candidate requeueing was
local, and a fresh global minimum scan was not performed.  The producer code
now names this same event `local_priority_queue_exhausted`.

## Recheck

From the project parent directory:

```sh
python j_source_v1/diagram/replay_unmarked4.py
```

The replayer uses only the Python standard library and checks the exact current
source/master bindings before accepting the trace.
