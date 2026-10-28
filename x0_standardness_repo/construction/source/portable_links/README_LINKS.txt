X0 VERTEX-LINK CERTIFICATE

Result: every vertex link of the frozen 103,744-pentachoron X0 cover master is a
PL 3-sphere. This establishes that this master is a closed combinatorial
4-manifold. The certificate does not reduce the four-manifold to S4.

The 4,310 links are certified by 529,460 three-dimensional bistellar flips.
Each flip has its full-star and missing-face preconditions checked, and each
link ends as the boundary of a four-simplex. No manifold recognition software
or manifold premise is used for these three-dimensional reductions.

The X0/X2 cap replacement leaves 3,732 labelled vertex links exactly equal.
Their explicit X2 flip lists were copied only after exact equality of the
complete sorted tetrahedron lists. All 578 changed links were freshly reduced,
using 31,938 new flips, then locally replayed. Finally, a separate fresh process
reconstructed every X0 link from the frozen master and replayed all 529,460
moves, including the reused lists. This final replay is self-contained: neither
the X2 master nor the X2 link certificate is needed to verify X0.

Files:
- x0_vertex_links.json.gz: complete self-contained certificate.
- check_links.py: byte-identical copy of the reviewer's inspected producer and
  full-precondition replay, retained as the elementary reference checker.
- replay_x0_links.py: source/hash/schema/coverage checks plus a fresh full replay;
  its reference checker hash is pinned in code.
- produce_x0_links.py: bounded production and exact-equality reuse provenance.
- production_receipt.json: historical production result; full replay was still
  pending at that stage.
- replay_receipt.json and replay.log: completed fresh all-links PASS evidence.

Recommended fresh replay from the X0 package root (Python standard library):

python3 portable_links/replay_x0_links.py --receipt /path/to/fresh/x0_links_replay.json

The wrapper defaults to candidate/marked_candidate.json.gz and the two local
portable_links certificate/checker files. Explicit --master, --certificate,
--checker and --receipt paths are supported. It applies a 4 GiB address-space
ceiling and a 120 CPU-second ceiling using the Unix resource module. It refuses
Python optimization, since the reference checker uses assertions.

The elementary checker itself is also portable standard Python:

python3 portable_links/check_links.py verify candidate/marked_candidate.json.gz portable_links/x0_vertex_links.json.gz

Use ordinary Python without -O or PYTHONOPTIMIZE. This elementary command checks
all actual moves and vertices but omits the wrapper's resource and file-binding
checks; the package manifest supplies the file identity.

Measured production: 27.304 CPU seconds, 27.311 wall seconds, 489,627,648 bytes
peak RSS. Fresh complete replay: 19.391 CPU seconds, 19.397 wall seconds,
665,632,768 bytes peak RSS. Executed on Linux; the explicit combinatorial flip
sequences do not depend on operating-system random-number engines.

Frozen master SHA256:
251d7ad38235869fa268e9ee6ec9da2f397cf1ca71db267d0bf0ca9b434c5e4c
Certificate SHA256:
ae285aa6b2dcbc4e23584a5480fd6a0e2ee916a56179b8d39cd1b8ff15d79c47
Reference checker SHA256:
aa6a86a25e6a79975cc3d90076e2861e2e407ec216a903b26331ac9e0bb2638c
