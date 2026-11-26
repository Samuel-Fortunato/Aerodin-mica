#!/usr/bin/env bash
set -euo pipefail

# Run from the case directory, independent of the caller's working directory.
cd "$(dirname "$0")"

# Start from a clean mesh/solution. Remove processor directories if present.
rm -rf constant/polyMesh processor*
foamListTimes -rm 2>/dev/null || true

blockMesh | tee log.blockMesh

# Parallel meshing
decomposePar -copyZero | tee log.decomposePar
mpirun -np 6 --bind-to core --map-by core snappyHexMesh -parallel -overwrite | tee log.snappyHexMesh
reconstructParMesh -constant -mergeTol 1e-6 | tee log.reconstructParMesh.log

# Serial post-processing of the mesh
extrudeMesh | tee log.extrudeMesh
createPatch -overwrite | tee log.createPatch
checkMesh -allTopology -allGeometry | tee log.checkMesh

# Re-decompose
rm -rf processor*
cp -r 0.orig 0
decomposePar | tee log.decomposePar.final

# Parallel run
mpirun -np 6 --bind-to core --map-by core pimpleFoam -parallel | tee log.pimpleFoam
reconstructPar | tee log.reconstructPar
