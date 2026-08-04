# EventMixing

ROOT/C++ tools to perform **event mixing** for combinatorial background estimation in ALICE (O2 framework) invariant-mass and/or femtoscopy analyses. The code builds mixed-event (or same-event, angle-rotated) pairs from candidate trees produced by O2 derived-data tasks, and fills the resulting invariant-mass spectra together with QA histograms.

Two analyses are currently supported, sharing a common core:

- **Li4** (`include/li4`) — mixing of ³He and hadron (proton) candidates, used for the search of the `Li-4` (³He–p) system.
- **H4L** (`include/h4l`) — mixing of ⁴He and pion candidates, used for hypernuclei (H4Λ-like, ⁴He–π) analyses.

The main executable currently wired up in `src/` is `mixingLi4.cxx`.

## Repository layout

```
.
├── build.cpp              # ROOT macro to compile the analysis with ACLiC
├── config/
│   └── configMixingLi4.yml   # Example YAML configuration for the Li4 mixing
├── include/
│   ├── core/               # Shared utilities (physics helpers, tree I/O, generic mixing types)
│   │   ├── candidates.hh
│   │   ├── indexTableUtils.hh
│   │   ├── mixing.hh
│   │   ├── physics.hh
│   │   └── treeUtils.hh
│   ├── li4/                 # ³He–hadron (Li4) candidate classes, mixer and QA histograms
│   │   ├── histograms.hh
│   │   ├── li4candidates.hh
│   │   ├── mixing.hh
│   │   └── selections.h
│   └── h4l/                 # ⁴He–pion (H4L) candidate classes and mixer
│       ├── h4lcandidates.hh
│       └── mixing.hh
└── src/
    └── mixingLi4.cxx        # Entry-point macro/executable for the Li4 event mixing
```

## What the code does

`mixingLi4.cxx` drives the workflow:

1. **Merge** (optional, `doMerge` in the config) — reads the O2 derived-data trees (`O2he3hadtable` for candidates, `O2he3hadmult` for collision multiplicity/vertex info) from possibly multiple directories in the input file and merges them into single trees (`output/inputCands.root`, `output/inputColls.root`) via `treeUtils::treeMerging`.
2. **Load** — fills in-memory vectors of `He3Candidate`, `HadCandidate` and `CollisionCandidate` objects from the merged trees, grouping hadron indices into per-collision brackets (`mixing::fillParticlesFromTree`), optionally applying preliminary quality/PID cuts.
3. **Mix** — a `Mixer` object combines candidates from different (or the same, rotated) collisions using one of two strategies:
   - `kEvent` (`mixingStrategy: 0`) — standard event mixing: pair a candidate with candidates from other, similar (vertex-z / multiplicity binned) collisions.
   - `kRotation` (`mixingStrategy: 1`) — same-event background estimate obtained by randomly rotating the azimuthal angle of one of the daughters.
4. **Output** — writes the mixed pairs to a `MixedTree` and saves QA histograms (candidate `p_T` before/after mixing, like-sign/unlike-sign invariant-mass spectra, comparison canvases) to the output ROOT file.

Binning for the event-mixing pools (vertex-z and multiplicity/centrality) is handled by `HistVertexMultiplicity` in `include/core/indexTableUtils.hh`.

## Dependencies

- [ROOT](https://root.cern/) (tested with ROOT 6.32) with `RDataFrame`, `Math/Vector4D`, `TTree`, etc.
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) for reading the YAML configuration.
- A C++17-capable compiler (the code is compiled through ROOT's ACLiC).

## Building

The project is built with ROOT's ACLiC via the `build.cpp` macro, which compiles `src/mixingLi4.cxx` into the `build/` directory and links against `yaml-cpp`.

```bash
root -l -b 'build.cpp("fast")'   # incremental build
root -l -b 'build.cpp("force")'  # force a full recompilation
```

`build.cpp` currently points to a `yaml-cpp` installation under `/home/galucia/yaml-cpp`; update the `gSystem->AddIncludePath(...)` and `gSystem->Load(...)` calls to match your local `yaml-cpp` include/library paths before building.

## Configuration

The mixing job is configured through a YAML file, e.g. `config/configMixingLi4.yml`:

```yaml
inputFileName: "/path/to/input_same.root"   # O2 derived-data file (candidates + collisions)
outputFileName: "/path/to/output_mixing.root"

inputCands: "output/inputCands.root"        # intermediate merged candidate tree
inputColls: "output/inputColls.root"        # intermediate merged collision tree

doMerge: true          # merge per-directory trees from inputFileName before mixing
mixingStrategy: 0       # 0: event mixing, 1: angle/rotation mixing
mixingDepth: 10          # number of mixed partners per candidate
randomSeed: 256          # seed for the random-number generator
is23: false               # dataset-specific flag (affects preliminary cuts)
applyCuts: false            # apply preliminary PID/DCA/quality cuts before mixing
```

## Running

Once built, load the macro and call the entry point from a ROOT session:

```bash
root -l
root [0] .x build.cpp
root [1] mixingLi4("config/configMixingLi4.yml")
```

This produces the configured `outputFileName` ROOT file containing:

- `MixedTree` — the tree of mixed (or rotated) ³He–hadron pairs, with the same branch layout as the input candidates.
- `HistogramsQA/` — a directory with QA histograms and comparison canvases (candidate kinematics and invariant-mass spectra before/after mixing, like-sign and unlike-sign).

## Status

This is an active analysis codebase for ALICE nuclear/hypernuclear physics studies; interfaces (especially the H4L mixer and the generic mixing utilities in `include/core/mixing.hh`) are still evolving, and some paths/parameters in `build.cpp` and the example config are specific to the author's local setup.
