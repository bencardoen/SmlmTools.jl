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
        # Random.seed!(42)
        # Cb1 = 100 .+ rand(100, 3) .* 50
        # Cb12 = 20 .+ rand(100, 3) .* 50
        # Cb2 = 105 .+ rand(100, 3) .* 60
        # Cb22 = 10 .+ rand(100, 3) .* 50
        # C1p = vcat(Cb1, Cb12)
        # C2p = vcat(Cb2, Cb22)
        # C1 = C1p
        # C2 = C2p
        # N=size(C1, 1)
        # MT = rand(N, 2) * 10000
        # MT[:,2] .= Int.(round.(MT[:, 2]))
        # MT[:,2] .= sort(MT[:,2])
        # using DataFrames, CSV
        # vs = hcat([C1, MT]...)
        # CSV.write("1p.csv", DataFrame(vs, :auto))
        # vs = hcat([C2, MT]...)
        # CSV.write("2p.csv", DataFrame(vs, :auto))
        first="1p.csv"
        second="2p.csv"
        if !isfile(first)
            first="test/$first"
            @test isfile(first)
        end
        if !isfile(second)
            first="test/$second"
            @test isfile(second)
        end
        r=align(first, second; outdir=t, type="thunderstorm")
        @test isdir(t)
        rm(t;recursive=true)
        @test !isdir(t)
    end

    @testset "bmp" begin
        xs = zeros(10, 10)
        xs[2:3,2:3] .= 1
        ys = zeros(10, 10)
        ys[10, 10] = 1
        ys[4:5, 4:5] .=1
        xi, yi, ind, d = minpair(xs, ys)
        @test ind == (1, 1)
        @test d == sqrt(8)
        @test sum(xi) == sum(xs)
        @test sum(yi) <= sum(ys)
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
        # C1 = 10 .+ rand(100000, 3) .* 10000
        # C2 = 10 .+ rand(100000, 3) .* 10000
        Cb1 = 100 .+ rand(100, 3) .* 50
        Cb12 = 20 .+ rand(100, 3) .* 50
        Cb2 = 105 .+ rand(100, 3) .* 60
        Cb22 = 10 .+ rand(100, 3) .* 50


        C1 = vcat(Cb1, Cb12)
        C2 = vcat(Cb2, Cb22)
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
