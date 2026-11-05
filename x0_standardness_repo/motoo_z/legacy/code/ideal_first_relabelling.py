"""Lawful cache-friendly relabelling of a valid ideal Regina triangulation.

No cache is edited and no vertex type is presumed. Vertex types come from
Regina's completed skeleton computation. randomiseLabelling() retains Regina's
already established topological properties. The returned explicit isomorphism
is sufficient to verify the relabelling independent of the random sequence.

For best performance after this relabelling, do moves whose removed stars
exclude simplex zero. Its local vertex zero then remains the ideal vertex.
"""
import regina


def ideal_first(tri, excluded_simplices=(), max_attempts=1000):
    """Relabel in place; return an explicit source-to-destination isomorphism.

The first simplex in the destination avoids ``excluded_simplices`` (indices
in the original input), and its local vertex zero is known ideal. Do not pass
stale face/vertex objects across this routine; map a local face through the
returned isomorphism instead.
"""
    n=tri.size()
    if n == 0:
        raise ValueError("empty triangulation")
    flags=[[p.vertex(v).isIdeal() for v in range(5)] for p in tri.simplices()]
    if not any(any(row) for row in flags):
        raise ValueError("no ideal vertex was recognized")
    excluded=set(excluded_simplices)
    allowed=[p not in excluded for p in range(n)]
    if not any(allowed[p] and any(flags[p]) for p in range(n)):
        raise ValueError("no ideal anchor outside the excluded star")
    total=regina.Isomorphism4.identity(n)
    attempts=0
    while not (flags[0][0] and allowed[0]):
        if attempts >= max_attempts:
            raise RuntimeError("failed to find an ideal-first relabelling")
        iso=tri.randomiseLabelling()
        new_flags=[[False]*5 for _ in range(n)]
        new_allowed=[False]*n
        for p in range(n):
            image=iso.simpImage(p);perm=iso.facetPerm(p)
            new_allowed[image]=allowed[p]
            for v in range(5):
                new_flags[image][perm[v]]=flags[p][v]
        flags,allowed=new_flags,new_allowed
        new_total=regina.Isomorphism4(n)
        for p in range(n):
            image=total.simpImage(p)
            new_total.setSimpImage(p,iso.simpImage(image))
            new_total.setFacetPerm(p,iso.facetPerm(image)*total.facetPerm(p))
        total=new_total
        attempts+=1
    return total,attempts


def encode_isomorphism(iso):
    return [[iso.simpImage(p),[iso.facetPerm(p)[v] for v in range(5)]] for p in range(iso.size())]


def decode_isomorphism(data):
    iso=regina.Isomorphism4(len(data))
    for p,(image,perm) in enumerate(data):
        iso.setSimpImage(p,image)
        iso.setFacetPerm(p,regina.Perm5(*perm))
    if sorted(iso.simpImage(p) for p in range(iso.size())) != list(range(iso.size())):
        raise ValueError("not a bijection of simplices")
    return iso
