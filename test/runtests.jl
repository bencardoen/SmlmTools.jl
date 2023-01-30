using SmlmTools
using Test
using PyCall
using Random
using Statistics

@testset "SmlmTools.jl" begin
    @testset "OV" begin
        @test onlineVariance([1 2 3]) == (2.0, 0.6666666666666666)
    end

    @testset "Full run" begin
        @info pwd()
        t=mktempdir()
        first="1.csv"
        second="2.csv"
        d = pwd()
        if ! (isfile(first) && isfile(second))
            @error "Failed finding $first $second in $(pwd())"
            if isfile(joinpath("test", first))
                @info "changing path"
                cd("test")
            end
        end
        r=align(first, second; outdir=t)
        if pwd() != d
            @info "resetting path to $(d)"
            cd(d)
        end
        @test isdir(t)
        rm(t;recursive=true)
        @test !isdir(t)
    end


    @testset "sm" begin
        xs = rand(100, 3) .* 100
        for _ in 1:10
            sx = sample_mean(xs)
            m = mean(sx)
            @test 40 < m < 60
        end
        xs = rand(100, 3) .* 100
        mx=sample_mean(xs; seed=42)
        my=sample_mean(xs; seed=42)
        @test all(mx .== my)
    end

    @testset "tools" begin
        @test iszero(logz(0))
        @test logz(ℯ) == log(ℯ)
        @test maximum(nmz([1 2])) <= 1.0
        xs = zeros(2, 3)
        xs[1, :] .= [1,1,1]
        xs[2, :] .= [2,2,2]
        ys = adjust_to_first(xs)
        @test sum(ys) == 3
    end

    @testset "project" begin
        Random.seed!(42)
        C1 = 10 .+ rand(100000, 3) .* 10000
        C2 = 10 .+ rand(100000, 3) .* 10000
        rs = project_image(C1, 10)
        bd=detect_bead(C1, C2, 10)[3]
        @test maximum(bd) == 2
    end

    @testset "track" begin
        C1 = 10 .+ rand(100000, 3) .* 10000
        N=size(C1, 1)
        MT = rand(N, 2) * 10000
        MT[:,2] .= Int.(round.(MT[:, 2]))
        MT[:,2] .= sort(MT[:,2])
        sms_c, sts_c = track_sample_mean(C1, MT, 100, 1000)
        size(sms_c, 1) == 900
        isapprox(mean(mean(sms_c, dims=1))/1000, 5.0, atol=.1)
    end
end
