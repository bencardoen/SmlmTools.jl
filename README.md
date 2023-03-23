# SmlmTools

A set of tools for processing point cloud based superresolution/single molecule localization microscopy, including but not limited to 
- point cloud to image conversion
- fiducial detection/tracking
- cross-channel alignment
- temporal drift correction

While microscopes often have built-in or software solutions to do this for you, they may not be perfect or you may no longer have access to the microscope.
You can use this software when:
- you have a point cloud dataset
- no access to microscope or raw data
- you're not happy with the current correction
- correction is in 2D, not 3D
- ...

Code Coverage [![codecov](https://codecov.io/gh/bencardoen/SmlmTools.jl/branch/master/graph/badge.svg?token=qFQ3PGsBBY)](https://codecov.io/gh/bencardoen/SmlmTools.jl)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7632321.svg)](https://doi.org/10.5281/zenodo.7632321)

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

## Algorithm
You can review the documentation of the code, but in short this is in plain English what will happen:
- Load the dataset (3D point clouds) for both channels
- Find up to k fiducials (default 2) per channel
  - Use the fact that fiducials in SMLM continuously emit, so find peaks in the spatial density distribution (those would show up as 'bright points')
- Pair the fiducials across channels
- If the distance between the closest pair is > than a threshold (default 400nm), refuse to continue. You can increase this, but be careful, you could be aligning a fiducial close to a bright part to the cell if you do : k = 5, distance = 1000nm, and there are only 2 fiducials.
- For each channel
  - Look at the mean location of the fiducial over time
  - Use this offset over time to correct all points in this channel
  - The trajectory will be plotted for you so you can inspect it
- Once corrected across time
  - Find the linear translation (offset) between the time-corrected fiducials, and use this align the channels

To allow you to inspect each stage, plots are saved with before/after data and for example the temporal trajectory in XYZ of the beads.

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
This example runs on Compute Canada
```bash
# Change to local scratch directory
cd scratch
# Copy the singularity image
cp /project/rrg-hamarneh/data/nabi-robert-ivan/software/smlmtools/SmlmTools.sif .
# Set executable
chmod u+x SmlmToools.sif
# Request resources to run from SLURM
salloc --mem=32GB --account=rrg-hamarneh --cpus-per-task=4 --time=3:00:00
```
#### Run
```bash
# Make sure Singularity has rw access
export SINGULARITY_BINDPATH="/scratch/$USER,$SLURM_TMPDIR"
# Execute
./SmlmTools.sif -f testdata/1C8PTRF_3_1_Cav_647.bin -s testdata/1C8PTRF_3_1_PTRF_568.bin -t GSD -p 10 -a -n 10 --outdir . -c
```
This will produce colocalization for 7 metrics, 2D image projections and 3D alignment.


### Cite
If you find this useful, consider citing
```bibtext
@software{cardoen_ben_2023_7632321,
  author       = {Cardoen, Ben},
  title        = {{SmlmTools: A Julia package for computational 
                   methods for single molecule localization /
                   superresolution microscopy}},
  month        = feb,
  year         = 2023,
  note         = {https://github.com/bencardoen/SmlmTools.jl},
  publisher    = {Zenodo},
  version      = {0.1},
  doi          = {10.5281/zenodo.7632321},
  url          = {https://doi.org/10.5281/zenodo.7632321}
}
``` 

#### See also
- [Colocalization](https://github.com/bencardoen/Colocalization.jl)
- [DataCurator](https://github.com/bencardoen/DataCurator.jl) This package allows you to run the above without any code, as a part of your pre or postprocessing workflows
