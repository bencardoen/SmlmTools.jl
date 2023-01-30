export LOCALPKG=/opt/SmlmTools.jl
    #/opt/julia/julia-1.6.2/bin
export JLMJV=1.7
export JLV=$JLMJV.1
export PATH=/opt/julia/julia-$JLV/bin:$PATH
export JULIA_DEPOT_PATH=/opt/juliadepot
export GKSwstype=100
julia --project=/opt/SmlmTools.jl --sysimage=/opt/SmlmTools.jl/sys_img.so /opt/SmlmTools.jl/scripts/align.jl "$@"
