# SmlmTools

A set of tools for processing point cloud based superresolution/single molecule localization microscopy, including but not limited to 
- point cloud to image conversion
- fiducial tracking
- cross-channel alignment

Code Coverage [![codecov](https://codecov.io/gh/bencardoen/SmlmTools.jl/branch/master/graph/badge.svg?token=qFQ3PGsBBY)](https://codecov.io/gh/bencardoen/SmlmTools.jl)

Automated testing [![CircleCI](https://dl.circleci.com/status-badge/img/gh/bencardoen/SmlmTools.jl/tree/master.svg?style=svg&circle-token=51454c475b36421e7f42be42ebcf3dea1b77c483)](https://dl.circleci.com/status-badge/redirect/gh/bencardoen/SmlmTools.jl/tree/master)

## Installation
```bash
git clone https://github.com/bencardoen/SmlmTools.jl.git
cd SmlmTools.jl
julia --project=.
```
or adding as a package
```bash
julia -e 'using Pkg; Pkg.add(url="https://github.com/bencardoen/SmlmTools.jl.git")'
```

## Usage
### 2-channel alignment
Let F and S be the file names of the GSD bin files of either channel:
```bash
julia --project=. scripts/align.jl --f [F] --s [S] --outdir [mydirectory] --colocalize --align --segment
```
Adding --colocalize runs colocalization metrics.

### Point cloud to image
If you have a 3D point cloud in coords3d (Nx3), project image using nm_per_px.
```julia
image, smoothed, density,  beadmask = project_image(coords3d, nm_per_px; mx=nothing, remove_bead=false, log_scale=true, σnm=10)
```

### Using on cluster
- Log in
#### Get resources
```bash
cd scratch
cp /project/rrg-hamarneh/data/nabi-robert-ivan/software/smlmtools/SmlmTools.sif .
salloc --mem=32GB --account=rrg-hamarneh --cpus-per-task=4 --time=3:00:00
```
#### Run
```bash
export SINGULARITY_BINDPATH="/scratch/$USER,$SLURM_TMPDIR"
./SmlmTools.sif -f testdata/1C8PTRF_3_1_Cav_647.bin -s testdata/1C8PTRF_3_1_PTRF_568.bin -t GSD -p 10 -a -n 10 --outdir . -c
```
This will produce colocalization for 7 metrics, 2D image projections and 3D alignment.


#### See also
- [Colocalization](https://github.com/bencardoen/Colocalization.jl)
