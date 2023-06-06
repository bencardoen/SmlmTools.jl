# Copyright (C) 2018-2023 Ben Cardoen bcardoen@sfu.ca
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


module SmlmTools

import Base.Filesystem
using Distributions
using Statistics
using ImageFiltering
using ImageMorphology
using Printf
using Match
using HypothesisTests
using PyCall
using LaTeXStrings
using Images, DataFrames, CSV, LinearAlgebra
import Glob
using Logging
import Random
using ProgressMeter
using LoggingExtras
using Dates
using Base.Threads
# using Plots
# @pyimport scipy
# @pyimport scipy.stats as st

# keys = ["x", "y", "z", "dx", "dy", "dz", "amplitude", "frame_number"]
export parseSRFile, onlineVariance, estimateError, getValues, getTimes, getMagnitude, freqCount, nmz, sample_mean, detect_bead, logz,
project_image, adjust_to_first, align_using_time_mean, minpair, readfile, roi_pts, findbeads, track_sample_mean, beadcoords, summarize_colocalization, align, nmz


function sample_mean(xs; reps=100, seed=0)
    if seed != 0
        Random.seed!(seed)
    end
    nrows = size(xs, 1)
    sel=Int(floor(nrows/2))
    rowindices = range(1, nrows) |> collect
    ms = zeros(reps, 3)
    for r in 1:reps
        sampled_indices=sample(rowindices, sel; replace=true)
        sm = mean(xs[sampled_indices, :], dims=1)
        ms[r, :] .= sm[1, 1:3]
    end
    return ms
end

function findbeads(image, nbeads=1)
    ccs = label_components(tomask(image))
    N = maximum(ccs)
    vals = zeros(N)
    for (ic, ci) in enumerate(component_indices(ccs)[2:end])
        vals[ic] = maximum(image[ci])
    end
    si=sortperm(vals)
    beads = aszero(image)
    cind = component_indices(ccs)[2:end]
    for i in 1:nbeads
        m = si[end-i+1]
        beads[cind[m]] .= 1
    end
    return tomask(beads)
end


function tomask(xs)
	ys = copy(xs)
	ys[ys .> 0] .= 1
	return ys
end


function aszero(xs)
	return zeros(eltype(xs), size(xs))
end

"""
    detect_bead(xs, ys, nm_per_px, beads=1)

    Detect two closest fiducials in xs, ys, considering up to `beads`.

    Returns the bead masks, the union of the masks, the projection images with bead removed and without removal, as well as the indices of the beads and their minimmum distance 
"""
function detect_bead(coordsc1, coordsc2, nm_per_px, beads=1; σ=10, maxdistance=500)
    qx = max(maximum(coordsc1), maximum(coordsc2))
    @debug "Maximum range = $qx"
    i1, i2 = [project_image(c, nm_per_px; mx=qx, remove_bead=true, log_scale=false, σnm=σ, beads=beads) for c in [coordsc1, coordsc2]]
	r1, r2 = [project_image(c, nm_per_px; mx=qx, remove_bead=false, log_scale=false, σnm=σ, beads=beads) for c in [coordsc1, coordsc2]]
    d1, d2 = i1[end], i2[end]
    @showprogress for i in 1:beads
        # @info "Testing bead count $i"
        m1i = findbeads(d1, i)
        m2i = findbeads(d2, i)
        p1, p2, indices, distance = minpair(m1i, m2i)
        @info "Nearest pair is $(indices) with distance $(distance) px = $(distance*nm_per_px) nm"
        if distance*nm_per_px < maxdistance
            @info "Found configuration of $i beads with 1 pair close enough"
            # @warn "Removeme"
            un = dilate(dilate(p1)) .+ dilate(dilate(p2))
            return d1, d2, un, i1, i2, r1, r2, indices, distance 
        end
    end
    m1i = findbeads(d1, beads)
    m2i = findbeads(d2, beads)
    p1, p2, indices, distance = minpair(m1i, m2i)
    @info "Nearest pair is $(indices) with distance $(distance) px = $(distance*nm_per_px) nm"
	un = p1 .+ p2
    return d1, d2, un, i1, i2, r1, r2, indices, distance 
end


"""
    minpair(xs, ys)

For two binary masks, find the pair of components that are closest to each other.
Return the masks of the two components, the indices of the components, and the distance between the components.
"""
function minpair(xs, ys)
    if (sum(xs) < 1) || (sum(ys) < 1)
        @warn "No fiducials in one or both channels !!!!"
        return zeros(eltype(xs), size(xs)), zeros(eltype(ys), size(ys)), (-1, -1), Inf
    end
    cssx = Images.label_components(xs)
    ctrx = Images.component_centroids(cssx)[2:end]
    xin = Images.component_indices(cssx)[2:end]
    cssy = Images.label_components(ys)
    ctry = Images.component_centroids(cssy)[2:end]
    yin = Images.component_indices(cssy)[2:end]
    # @info ctrx
    # @info ctry
    md = Inf
    mindices = -1, -1
    for (i, ctr_x) in enumerate(ctrx)  
        for (j, ctr_y) in enumerate(ctry)
            x1, y1 = ctr_x
            x2, y2 = ctr_y
            d = sqrt((x1-x2)^2 + (y1-y2)^2)
            if d < md
                md = d
                mindices = i, j
            end
        end
    end
    xsi = zeros(size(xs))
    ysj = zeros(size(ys))
    i, j = mindices
    xsi[xin[i]] .= 1
    ysj[yin[j]] .= 1
    return xsi, ysj, mindices, md
end


function nmz(xs)
    return xs ./ maximum(xs)
end

function readfile(filename, args)
    fext = split(filename, ".")[end]
    @match fext begin
        "bin" => readgsd(filename, args)
        "ascii" => readgsd(filename, args, false)
        "csv" => readcsv(filename, args)
        _ => throw(ArgumentError("Unknown file type $fext"))
    end
end

function readgsd(filename, args, binary=true)
    @debug "Reading binary GSD format with $filename $args"
    s=pyimport("smlmvis.gsdreader")
    r=s.GSDReader(filename, binary=binary)
    pts = copy(r.points)
    pts[:,1:2] .*= args["gsd_nmpx"]    
    meta=copy(r.values)
    return pts, meta
end

function readcsv(filename, args)
    @debug "Reading CSV data"
    C1 = CSV.read(filename, DataFrame)
    if args["type"] != "thunderstorm"
        @debug "Reading generic CSV with x,y,z,id,frame"
        return Matrix(C1[:,1:3]), Matrix(C1[:,4:end])
    else
        @debug "Reading Thunderstorm CSV with x [nm], y [nm], id, frame"
        cols = names(C1)
        if "z [nm]" in cols
            @info "Loading Z data from Thunderstorm"
            pts, meta = Matrix(C1[:,["x [nm]", "y [nm]", "z [nm]"]]), Matrix(C1[:,["id", "frame"]])
        else
            @warn "No Z data in Thunderstorm CSV, using 0"
            pts, meta = Matrix(C1[:,["x [nm]", "y [nm]"]]), Matrix(C1[:,["id", "frame"]])
            pts=hcat(pts, zeros(size(pts,1)))
        end
        return pts, meta       
    end
end
"""
	align(first, second; outdir=".",  nm_per_px=10, σ=10, gsd_nmpx=159.9, maxframe=20000, interval=4000)

	Align 2 point clouds using their metadata, encoded in GSD bin files (first and second are filenames)

	Saves plots of the correction, 2D projected images, as well as CSV files of the correct 3D point clouds.

	First, second are eiter GSD superres format or CSV, in the 2nd case the first 3 columns should be X, Y, Z in nm, and time as 5th column.
    If `type` is set to `thunderstorm`, will read columns `x [nm]` and `y [nm]` as well as `frame` and `id`. 
    In this case a zero 3rd dimension is used.
"""
function align(first, second; outdir=".",  nm_per_px=10, σ=10, gsd_nmpx=159.9, maxframe=nothing, interval=4000, type="gsd", maxbeaddistancenm=300, maxbeads=2)
	fext = split(first, ".")[end]
    @info "Loading $first and $second"
    if isnothing(maxframe)
        @info "Maximum frame not specified, using all frames"
    end
    @info "Using $nm_per_px nm per pixel, max bead distance nm=$maxbeaddistancenm nm, max beads=$maxbeads"
    if ! (fext in ["ascii", "bin", "csv"])
        @error "Unsupported files : should be CSV or GSD bin/ascii"
    end
    args = Dict("gsd_nmpx"=>gsd_nmpx, "type"=>type)
    C1_pts, C1_meta_all = readfile(first, args)
    C2_pts, C2_meta_all = readfile(second, args)
    @info "Detecting bead ..."
	mx = max(maximum(C1_pts[:, 1:2]), maximum(C2_pts[:, 1:2]))
    bd = detect_bead(C1_pts[:, 1:2], C2_pts[:, 1:2], nm_per_px, maxbeads, maxdistance=maxbeaddistancenm, σ=σ)

	_, _, _, _, _, i1, i2, _, dist =bd
    if dist*nm_per_px > maxbeaddistancenm
        @error "Nearest Beads in both channels are too far apart ($dist nm > $(maxbeaddistancenm)) !"
        throw(ArgumentError("Beads are too far apart ($(dist*nm_per_px) nm > $(maxbeaddistancenm)) !"))
    end
	Images.save(joinpath(outdir, "C1_notaligned.tif"), N0f16.(nmz(i1[2])))
	Images.save(joinpath(outdir, "C2_notaligned.tif"), N0f16.(nmz(i2[2])))
    a, b = bd[4], bd[5]
    (minx, miny), (maxx, maxy) = beadcoords(bd[3])

    X = (minx-1)*nm_per_px, (maxx+1)*nm_per_px
    Y = (miny-1)*nm_per_px, (maxy+1)*nm_per_px

    @info "Estimated Fiducal location (nm): X $X Y $Y"
    C1_fiducial, C1_meta = roi_pts(C1_pts, X, Y, C1_meta_all)
    C2_fiducial, C2_meta = roi_pts(C2_pts, X, Y, C2_meta_all)

    # q=Plots.scatter(C1_fiducial[:,1], C1_fiducial[:,2] , alpha=.25, markersize=2, color=:blue, label="C1")
    # Plots.scatter!(C2_fiducial[:,1], C2_fiducial[:,2] , alpha=.25, markersize=2, color=:red, label="C2")
    # q=Plots.plot(q, dpi=300, size=(800, 600), xlabel="X axis nm", ylabel="Y axis nm", title="Bead location in XY")
	# _f = joinpath(outdir, "bead.svg")
    # @info "Saving fiducal XY plot in: $(_f)"
    # Plots.savefig(_f)
    if isnothing(maxframe)
        Max1 = maximum(C1_meta[:,2])
    else
        Max1 = maxframe
    end
    if isnothing(maxframe)
        Max2 = maximum(C2_meta[:,2])
    else
        Max2 = maxframe
    end
    # MAXFRAME=maxframe
    INTERVAL=interval
    # @debug "Change maxframe to be dynamic"
    @info "Max frame set to $Max1 $Max2"
    @info "Computing trajectory of bead over time"
    ### Find the sampled mean location of each fiducial over time
    sms_c, sts_c = track_sample_mean(C2_fiducial, C2_meta, INTERVAL, Max1)
    sms_p, sts_p = track_sample_mean(C1_fiducial, C1_meta, INTERVAL, Max2)
    ### Compute the relative offset over time compared to t=0
    offset_cav = adjust_to_first(sms_c)
    offset_ptrf = adjust_to_first(sms_p)
    #
    # q=Plots.scatter(offset_cav[:,1], markershape=:cross, offset_cav[:,2], alpha=.5, markersize=5,label="C2", color=:reds, marker_z=range(1, size(offset_cav, 1)))
    # Plots.scatter!(offset_ptrf[:,1], offset_ptrf[:,2], markershape=:xcross, alpha=.5, dpi=150, markersize=5, label="C1", color=:blues, marker_z=range(1, size(offset_ptrf, 1)), xlabel="Fiducial trajectory over time (X, nm)", ylabel="Fiducial trajectory over time (Y, nm)")
    # Plots.savefig(joinpath(outdir, "bead_trajectory.svg"))
    # @debug "Saving trajectory plot to $(joinpath(outdir, "bead_trajectory.svg"))"

    ### Compute the translation between T=1 for each channel
    offset_translate = sms_p[1, :] .- sms_c[1, :]

    @info "Offset between channels at time 0: X/Y/Z in nm:  $(offset_translate)"
    ### Aligned the locations using the timed offset
    ptrf_aligned_timed = vcat(align_using_time_mean(C1_fiducial, offset_ptrf, C1_meta)...)
    cav_aligned_timed = vcat(align_using_time_mean(C2_fiducial, offset_cav, C2_meta)...)
    pc = copy(ptrf_aligned_timed)
    pc[:, 3] .-= offset_translate[3]
    pc[:, 1] .-= offset_translate[1]
    pc[:, 2] .-= offset_translate[2]

    # @debug "Plotting aligned fiducials"
    # q=Plots.scatter(pc[:,1], pc[:,2] , alpha=.25, markersize=2, color=:blue, label="C1 fiducial aligned")
    # Plots.scatter!(cav_aligned_timed[:,1], cav_aligned_timed[:,2] , alpha=.25, markersize=2, color=:red, label="C2 fiducial aligned")
    # q=Plots.plot(q, dpi=300, size=(800, 600), xlabel="X axis nm", ylabel="Y axis nm", title="Bead location in XY after alignment")
    # @debug "Saving fiducal XY plot in: $(joinpath(outdir, "bead_aligned.svg"))"
    # Plots.savefig(joinpath(outdir, "bead_aligned.svg"))

    @info "Aligning full channels"
    ptrf_aligned_time_full = vcat(align_using_time_mean(C1_pts, offset_ptrf, C1_meta_all)...)
    cav_aligned_time_full = vcat(align_using_time_mean(C2_pts, offset_cav, C2_meta_all)...)
    aligned_ptrf = copy(ptrf_aligned_time_full)
    aligned_ptrf[:,3] .-= offset_translate[3]
    aligned_ptrf[:,1] .-= offset_translate[1]
    aligned_ptrf[:,2] .-= offset_translate[2]

    # Plots.scatter(aligned_ptrf[:, 1], aligned_ptrf[:, 2], alpha=.125, markersize=2, color=:red, label="C1 aligned")
    # Plots.scatter!(cav_aligned_time_full[:, 1], dpi=300, size=(800, 600), cav_aligned_time_full[:, 2], alpha=.125, markersize=2, color=:blue, label="C2 aligned", xlabel="X nm", ylabel="Y nm")
    # @debug "Saving scatterplot of aligned channels to $(joinpath(outdir, "aligned_xy_both_channels.svg"))"
    # Plots.savefig(joinpath(outdir, "aligned_xy_both_channels.svg"))

    MX=max(maximum(aligned_ptrf[:,1:2]), maximum(cav_aligned_time_full[:,1:2]))
	MX=max(MX, mx)
    C1P = project_image(aligned_ptrf, nm_per_px; mx=MX, remove_bead=false, log_scale=false, σnm=σ)
    C2P = project_image(cav_aligned_time_full, nm_per_px; mx=MX, remove_bead=false, log_scale=false, σnm=σ)
	

    c1f = joinpath(outdir, "C1.tif")
    c2f = joinpath(outdir, "C2.tif")
    @debug "Saving projection 2D images to $(c1f) and $(c2f)"
    Images.save(c1f, N0f16.(nmz(C1P[2])))
    Images.save(c2f, N0f16.(nmz(C2P[2])))
    
	@debug "Saving points to CSV"
    pc = copy(aligned_ptrf)
    @info "Min $(minimum(pc[:,1:2])) - Max $(maximum(pc[:,1:2])) nm"
    df1a = DataFrame(xnm=pc[:,1], ynm=pc[:,2], znm=pc[:,3])
    df1a[!,:originalfile] = first
    df2a = DataFrame(xnm=qc[:,1], ynm=qc[:,2], znm=qc[:,3])
    df2a[!,:originalfile] = second
	CSV.write(joinpath(outdir, "aligned_c1.csv"), df1a)
	qc = copy(cav_aligned_time_full)
	CSV.write(joinpath(outdir, "aligned_c2.csv"), df2a)
    @debug "Done"
    @info "Writing to VTU"
    @info "Using filenames $(joinpath(outdir, "aligned_c1.vtu")) and $(joinpath(outdir, "aligned_c2.vtu"))"
	writetovtu(joinpath(outdir, "aligned_c1.vtu"), aligned_ptrf, aligned_ptrf[:,3:3])
	writetovtu(joinpath(outdir, "aligned_c2.vtu"), cav_aligned_time_full, cav_aligned_time_full[:,3:3])
	return aligned_ptrf, cav_aligned_time_full, C1P, C2P, bd
end


function writetovtu(fname, pts, meta)
    # @info "Writing $(size(pts)) to $(fname)"
    # @warn "Todo -- replace with https://juliapackages.com/p/writevtk"
    # s = pyimport("smlmvis.vtuwriter")
    # s.VtuWriter(fname, pts, meta)
end

"""
	logz(x)

	isinf(log(x)) ? 0 : log(x)
"""
function logz(x)
    if x == 0
        return 0
    else
        return log(x)
    end
end

"""
	project_image(coords3d, nm_per_px; mx=nothing, remove_bead=false, log_scale=true, σnm=10)

	Project the coordinates (Nx3) onto a 2D pixel grid (of nm_per_px).
	Returns the image, smoothed image, density and bead mask.

	Assumes the bead is identifiable by repeated blinks.
"""
function project_image(coords3d, nm_per_px; mx=nothing, remove_bead=false, log_scale=true, σnm=10, beads=1)
    C1X, C1Y = coords3d[:, 1], coords3d[:, 2]
    if isnothing(mx)
        @info "Finding pixel coordinates .."
        mx=max(maximum(C1X), maximum(C1Y))
    end
    N = Int(ceil(mx / nm_per_px))
    image = zeros(N, N)
    @info "Creating 2D images $N x $N pixels with $(nm_per_px) nm / px"
    @debug "Creating density map"
    npoints = 0
    for (cx, cy) in zip(C1X, C1Y)
        if cx <= 0 || cy <= 0
            npoints += 1
            continue
        end
        xc = Int(round(cx/nm_per_px))
        yc = Int(round(cy/nm_per_px))
        xc = max(xc, 1)
        yc = max(yc, 1)
        xc = min(xc, N)
        yc = min(yc, N)
        image[xc, yc] += 1
        # image[Int(round(cx/nm_per_px)), Int(round(cy/nm_per_px))] += 1
    end
    if npoints > 0
        @warn "You have $npoints negative coordinates!!"
    end
    dense = copy(image)
    beadmask = aszero(image)
    if remove_bead
		@info "Detecting bead"
		beadmask = findbeads(image, beads)
        beadmask = tomask(ImageMorphology.dilate(beadmask))
        image[beadmask .>= 1] .= 0
    end
    if log_scale
        @debug "Log transform of density"
        image[image .> 0] .= log.(image[image .> 0])
    end
    σ= σnm/nm_per_px
    @debug "Gaussian σ $(σnm) → $(σ) px"
	if length(logz.(image[image .> 0])[:]) != 0
	    qz=quantile(logz.(image[image .> 0])[:], 0.95)
	    image[logz.(image) .>= qz] .= exp(qz)
	else
		@warn "Image empty without beads ... --> bead detection too aggresive?"
	end
    # qz=quantile(logz.(image[image .> 0])[:], 0.95)
    # image[logz.(image) .>= qz] .= exp(qz)
    smooth =  ImageFiltering.imfilter(image, ImageFiltering.Kernel.gaussian((σ, σ)))
    return image, smooth, dense, beadmask
end


function adjust_to_first(xs)
    ys = zeros(size(xs)...)
    for row in 1:size(xs, 1)
        ys[row, :] = xs[row, :] .- xs[1, :]
    end
    return ys
end

"""
	align_using_time_mean(pts_unaligned, offset, metadata)

	Return aligned pts based on `offset`, an array indexed by time (`metadata[:,2]`).
	Iow `x, y, z = offset[i]` is the correction to be applied at time i to pts_unaligned[:, metadata[:,2].==i]
"""
function align_using_time_mean(pts_unaligned, offset, metadata)
    pts=copy(pts_unaligned)
    met=copy(metadata)
    PT=[]
    @showprogress for t in 1:size(offset, 1)
        mk= (met[:,2] .== t)
        _p = pts[mk,:]
        if length(_p) < 1
            continue
        end
        if any(isnan.(_p))
            @warn "Nan!"
            continue
        end
        for row in 1:size(_p, 1)
            _p[row, :] = _p[row, :] .- offset[t, :]
        end
        if any(isnan.(_p))
            @warn "Nan! offset"
            continue
        end
        push!(PT, _p)
    end
    pts = PT
    return pts
end

"""
	roi_pts(points, X, Y, metadata)

	Select a 2D ROI defined by xmin, xmax = X, ymin, ymax = Y.
	Apply the same mask to corresponding rows in metadata.
"""
function roi_pts(pts, X, Y, meta)
    m=(pts[:,1] .> X[1]) .& (pts[:,1] .< X[2]) .& (pts[:,2] .> Y[1]) .& (pts[:,2] .< Y[2])
    return copy(pts[m, :]), copy(meta[m,:])
end
#
# function track_sample_mean(pts, meta, interval, maxframe)
#     span = Int(floor(interval/2))
#     m=minimum(meta[:, 2])
#     M=min(maximum(meta[:, 2]), maxframe)-1
#     @debug "Min time $m to max time $M"
#     PTS = copy(pts)
#     iters = range(Int(m+span),Int(M-span)) |> collect
#     N = length(iters)
#     ms = zeros(N, 3)
#     ts = zeros(N, 1)
#     # @info "Results size $(size(ms))"
#     @showprogress for (i,t) in enumerate(range(Int(m+span),Int(M-span)))
#         mks=(meta[:,2] .> t-span) .& (meta[:,2] .< t+span)
#         xs = PTS[mks,:]
#         sx = sample_mean(xs; seed=i+1)
#         _m = mean(sx, dims=1)
#         ms[i, :] .= _m[1, 1:3]
#         ts[i] = t
#     end
#     return ms, ts
# end
function track_sample_mean(pts, meta, interval, maxframe)
    span = Int(floor(interval/2))
    m=minimum(meta[:, 2])
    M=min(maximum(meta[:, 2]), maxframe)-1
    @info "Min time $m to max time $M"
    PTS = copy(pts)
    iters = range(Int(m+span),Int(M-span)) |> collect
    N = length(iters)
    ms = zeros(N, 3)
    ts = zeros(N, 1)
    # @info "Results size $(size(ms))"
    p = Progress(N)
    @threads for i in 1:N
        t = iters[i]
        mks=(meta[:,2] .> t-span) .& (meta[:,2] .< t+span)
        xs = @view PTS[mks,:]
        sx = sample_mean(xs; seed=i+1)
        _m = mean(sx, dims=1)
        ms[i, :] .= _m[1, 1:3]
        ts[i] = t
        next!(p)
    end
    return ms, ts
end








function beadcoords(bm)
    if maximum(bm) < 2
        @error "No valid mask"
        throw(ArgumentError("Invalid mask -- bead detection failed, try increasing minimum distance or number of beads"))
    end
    bmask = tomask(bm)
    ccs = Images.label_components(bmask)
    bmask = aszero(bmask)
    for (bx, ind) in zip(component_boxes(ccs), component_indices(ccs))
        if maximum(bm[ind]) == 2
            return bx
        end
    end
	@warn "No overlapping bead mask found!"
    return nothing
end

function parseSRFile(filename)
    f = open(filename)
    lines = readlines(f)
    content = []
    keys = ["x", "y", "z", "dx", "dy", "dz", "amplitude", "frame_number"]
    mM = Dict{String, Tuple{Float64, Float64}}([k=>(Inf, -Inf) for k in keys])
    @printf "File %s  has %d lines\n" filename length(lines)
    for (index, line) in enumerate(lines)
        row = split(line)
        d = Dict{String,Float64}()
        for (k,v) in zip(keys, row)
            value::Float64 = parse(Float64, v)
            d[k] = value
            m, M = mM[k]
            m = min(m, value)
            M = max(m, value)
            mM[k] = (m, M)
        end
        push!(content, d)
    end
    close(f)
    return content, mM
end

"""
Online mean and variance calculation based on
B. P. Welford (1962)."Note on a method for calculating corrected sums of squares and products".
Technometrics 4(3):419–420.
"""
function onlineVariance(data)
    n::Int64 = 0
    M::Float64, S::Float64 = 0.0, 0.0
    for (index::Int64, value::Float64) in enumerate(data)
        n = index
        Mo::Float64 = M
        M += (value - M) / index
        S += (value-M)*(value-Mo)
    end
    return M, S / n
end

function estimateError(db)
    lkeys = ["dx", "dy", "dz"]
    n::Int64 = 0
    M, S = [zeros(Float64, length(lkeys)) for _ in 1:2]
    for (index, row) in enumerate(db)
        n = index
        for (i, k) in enumerate(lkeys)
            value::Float64 = row[k]
            Mo = M[i]
            M[i] += ((value - M[i])  / index)
            S[i] += ((value - M[i])*(value - Mo))
        end
    end
    return M, map(x->x/n, S)
end

function getValues(db, key, f)
    return [f(row[key]) for row in db]
end

function windowValues(db, key, lower, upper)
    return [row[key] for row in db if row[key] >= lower && row[key] <= upper]
end

function getTimes(db)
    return getValues(db, "frame_number", x->x)
end

function getMagnitude(db)
    return getValues(db, "amplitude", x->x)
end

function freqCount(values)
    sd = Dict{Int64, Int64}()
    for v in values
        if haskey(sd, v)
            sd[v]+=1
        else
            sd[v]=1
        end
    end
    return sort(collect(sd), by=x->x[1])
end

function binValues(db)
    keys = ["x", "y", "z", "dx", "dy", "dz", "amplitude"]
    results = Dict{Int64, Any}()
    for row in db
        fn::Int64 = row["frame_number"]
        if !haskey(results, fn)
            results[fn] = Dict(k=>Array{Float64, 1}() for k in keys)
        end
        for key in keys
            value::Float64 = row[key]
            push!(results[fn][key], value)
        end
    end
    return sort(collect(results), by=x->x[1])
end

function unzip(lst)
    Nt = length(lst[1])
    lsts = [[] for _ in 1:Nt]
    for el in lst
        for index in 1:Nt
            push!(lsts[index],el[index])
        end
    end
    return lsts
end

function unzipr(lst)
    return reduce( ([],[]), lst) do x,y; for i in 1:length(y) push!(x[i], y[i]) end;x;end;
end

function normality(X)
    st=pyimport("scipy.stats")
    return st.anderson(X)
end

function evaluateNormality(X)
    r, l, p = normality(X)
    return Dict(pvalue=> (threshold>r) for (threshold, pvalue) in zip(l,p))
end

function evaluateNormalityVerbose(X)
    normal = evaluateNormality(X)
    found = false
    for n in keys(normal)
        if normal[n]
            @printf "Normal distribution for p value %.2f\n" normal[n]
            found = true
        end
    end
    if !found
        println("Values are not distributed on a normal distribution.")
    end
end


end # module
