function pmap_io(f, io, n, chunksize = 10000)
    np = nprocs()  # determine the number of processes available
    results = Array(Tuple{String,Float64},n)
    i = 1
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while !eof(io)
                        snp = readline(io)
                        idx = nextidx()
                        idx%chunksize == 0 && println("$(Dates.format(now(), "HH:MM")): working on SNP #$(idx): $(snp[1])")
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
B = ARGS[2]

@eval @everywhere B = $B
nsnp, ntraits = makecov(filein)

@everywhere begin
  using Distributions, DataFrames
  include("aspu_functions.jl")
  estv = readdlm("vcov_aspu.txt", ',')
  mvn = MvNormal(convert(Matrix{Float64}, estv))
  zb = Matrix{Float32}(9, ceil(Int64, 10^B)) #could be unsigned
  thisrun = aspurun(B, mvn, 3)
end

resio = open("aspu_results.txt", "w+")
insnp = open(filein)
readline(insnp)
println("Setup complete")

tic()
res = pmap_io(x->runsnp!(x,thisrun,zb), insnp, nsnp) #r
println("Job completed in $(round(toq()/3600,3)) hours")

exit()