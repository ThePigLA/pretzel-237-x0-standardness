"""Compare two marked links (archive component order): decorated signature + isometry matching
components in order, with every meridian and longitude sent to +-itself; decorations must agree."""
import json, sys, snappy
A = json.load(open(sys.argv[1])); B = json.load(open(sys.argv[2]))
da = A['decorations']; db = B.get('decorations') or [c.get('decoration', c.get('framing')) for c in B['components']]
da = [None if x == 'x' else x for x in da]; db = [None if x == 'x' else x for x in db]
print('decorations', da, db, 'equal' if da == db else 'DIFFER')
MA = snappy.Link(A['pd_code']).exterior(); MB = snappy.Link(B['pd_code']).exterior()
sa = MA.isometry_signature(of_link=True); sb = MB.isometry_signature(of_link=True)
print('signatures equal:', sa == sb, ' volumes %.6f %.6f' % (MA.volume(), MB.volume()))
good = 0
for iso in MA.is_isometric_to(MB, return_isometries=True):
    if iso.cusp_images() != list(range(MA.num_cusps())):
        continue
    if all(abs(m[0, 0]) == 1 and m[0, 1] == 0 and m[1, 0] == 0 and m[1, 1] == m[0, 0] for m in iso.cusp_maps()):
        good += 1
print('isometries preserving component order with meridians/longitudes -> +-themselves:', good)
print('CERTIFIED_SAME_MARKED_LINK' if (sa == sb and good and da == db) else 'NOT_CERTIFIED')
