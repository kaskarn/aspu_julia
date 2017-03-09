if haskey(ENV, "LSB_HOSTS")
  addprocs(split(ENV["LSB_HOSTS"]))
else
  addprocs(nprocs() - 1)
end

filein = ARGS[1]
using Distributions, DataFrames
include("aspu_functions.jl")

nsnp, snpnames, ntraits = makecov(filein)
@everywhere begin
  using Distributions, DataFrames
  include("aspu_functions.jl")
  B = 8 #will be input
  chunksize = 100000
  estv = readdlm("vcov_aspu.txt", ',')
  mvn = MvNormal(convert(Matrix{Float64}, estv))
  zb = Matrix{Float32}(9, ceil(Int64, 10^B)) #could be unsigned
  thisrun = aspurun(B, mvn, 3)
end

resio = open("test.txt", "w+")
insnp = open(filein)
readline(insnp)

chunkn = ceil(Int, nsnp/chunksize)
for i in 1:chunkn
  println("Chunk $i started at $(Dates.format(now(), "HH:MM"))")
  res = pmap(x->runsnp!(x,thisrun,zb), take(eachline(insnp),chunksize)) #replace with chunksize
  println("Chunk completed at $(Dates.format(now(), "HH:MM"))")
  for j in eachindex(res)
    join(resio, (res[j]..., "\n"), ",", "")
    (j%2000 == 0) && flush(resio)
  end
end
