function pmap_io(f, io, n, chunksize = 10000)
    np = nprocs()
    results = Array(Tuple{String, Float64, Int64}, n)
    i = 1
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while !eof(io)
                        snp = readline(io)
                        idx = nextidx()
                        idx%chunksize == 0 && println("$(Dates.format(now(), "HH:MM")): working on SNP #$(idx)")
                        results[idx] = remotecall_fetch(f, p, snp)
                    end
                end
            end
        end
    end
    results
end

if haskey(ENV, "LSB_HOSTS")
  addprocs(split(ENV["LSB_HOSTS"]))
else
  addprocs(nprocs() - 1)
end

filein = ARGS[1]
logB = parse(ARGS[2])
@eval @everywhere logB = $logB

chunksize = length(ARGS) > 2 ? ARGS[3] : 10000

using Distributions, DataFrames
@everywhere using Distributions, DataFrames
@everywhere include("aspu_utils.jl")

nsnp, ntraits = makecov(filein)

@everywhere begin
  estv = readdlm("vcov_aspu.txt", ',')
  mvn = MvNormal(convert(Matrix{Float32}, estv))
  thisrun = Aspurun(logB, mvn, 3)
  runvals = Aspuvals(
    zeros(UInt32,2, ceil(Int64, 10^logB)),
    Matrix{Float64}(9, ceil(Int64, 10^logB)),
    zeros(9),
    zeros(size(estv,2)),
    ones(9)
  )
end
gc()

insnp = open(filein)
readline(insnp)
println("Setup complete")

tic()
res = pmap_io(x->runsnp!(x, thisrun, runvals), insnp, nsnp, chunksize)

println("Job completed in $(round(toq()/3600,5)) hours")

writedlm("aspu_results.txt", res)
exit()
















#
