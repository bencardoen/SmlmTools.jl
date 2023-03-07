
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
using SmlmTools
using CSV, Statistics
import Glob
using Logging
using Dates
using ArgParse
using LoggingExtras
using SmlmTools
using Colocalization
using Images

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--first", "-f"
            help = "Filename of first channel"
            arg_type = String
            required = true
		"--type", "-t"
            help = "Type of file to process, defaults to GSD"
            arg_type = String
            default = "GSD"
		"--second", "-s"
            help = "Filename of second channel"
            arg_type = String
            required = true
		"--outdir", "-o"
            help = "output folder"
            arg_type = String
            required = true
		"--colocalize", "-c"
            help = "Set to true to do colocalization (default = false). If align is true, runs on the 2D projects produced after alignment completes. If not, expects 2 image files."
            action = :store_true
            default = false
		"--segment", "-g"
            help = "If active, segment the 2D images before colocalization. Not all coloc methods work without segmentation."
            action = :store_true
            default = false
		"--align", "-a"
            help = "Set to true to do alignment (default = true). Needs a SRM/SMLM type format, e.g. GSD's bin files"
            action = :store_true
            default = false
		"--nmpx", "-n"
			help = "Each pixel is n x n nm wide, default 10"
			arg_type = Float64
			default = 10.0
		"--precision", "-p"
			help = "Precision of system, σ used in creating 2D image of 3D point clouds, default 10"
			arg_type = Float64
			default = 10.0
		"--windowsize", "-w"
			help = "Windowsize used in colocalization (default = 3, should be odd >= 3)"
			arg_type = Int64
			default = 3
		"--beads", "-b"
			help = "Maximum expected nr of fiducials (default 2)"
			arg_type = Int64
			default = 2
		"--maxdistancebeads", "-m"
			help = "The maximum centroid to centroid distance between closest beads that is still acceptable, default 300nm"
			arg_type = Float64
			default = 300.0
    end

    return parse_args(s)
end

function runalign()
    date_format = "yyyy-mm-dd HH:MM:SS"
    timestamp_logger(logger) = TransformerLogger(logger) do log
      merge(log, (; message = "$(Dates.format(now(), date_format)) $(basename(log.file)):$(log.line): $(log.message)"))
    end
    ConsoleLogger(stdout, Logging.Info) |> timestamp_logger |> global_logger
    parsed_args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in parsed_args
        @info "  $arg  =>  $val"
    end
	first = parsed_args["first"]
	second = parsed_args["second"]
	outdir = parsed_args["outdir"]
	nm_per_px = parsed_args["nmpx"]
	σ = parsed_args["precision"]
	do_align =  parsed_args["align"]
	do_coloc = parsed_args["colocalize"]
	maxdistancebeads = parsed_args["maxdistancebeads"]
	res = nothing
	if do_align
		@info "Running alignment"
		res = align(first, second; type=parsed_args["type"], nm_per_px=nm_per_px, outdir=outdir, σ=σ, maxbeaddistancenm=maxdistancebeads, maxbeads=parsed_args["beads"])
	end
	if do_coloc
		if do_align
			_, _, C1P, C2P, _ = res
			_, C1, _, _ = C1P
			_, C2, _, _ = C2P
		else
			C1 = Images.load(first)
			C2 = Images.load(second)
		end
		if parsed_args["segment"]
			@info "Segmenting"
			C1 = segment(C1)
			C2 = segment(C2)
		end
		@info "Running colocalization"
		results = colocalize_all(C1, C2; windowsize=parsed_args["windowsize"])
		for k in keys(results)
			Images.save(joinpath(outdir, "$k.tif"), N0f16.(nmz(abs.(results[k]))))
		end
		df=summarize_colocalization(results, first, second)
		CSV.write(joinpath(outdir, "colocalization_results.csv"), df)
		@info "Done"
	end
end


runalign()
