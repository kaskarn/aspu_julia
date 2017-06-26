using DataFrames, Distributions

#Start workers
if haskey(ENV, "LSB_HOSTS")
  addprocs(split(ENV["LSB_HOSTS"])[2:end]) #keep master on own slot
else
  addprocs(Sys.CPU_CORES-1)
end

@everywhere using Distributions
@everywhere include("aspu_utils.jl")
@everywhere using aspu
using aspu_utils

#Parse options
# ARGS = ["--filein aspu_geno.csv --logB 7"] #example for testing
# ARGS = ["--filein aspu_geno.csv --replicate reptest_1e7.csv"]

inp = parse_aspu(ARGS)
isarg(a, d=inp) = in(a, keys(d))

println("\nInputs:")
display(inp)

logB = isarg("replicate") ? log10(nlines(inp["replicate"])) : parse(inp["logB"])
float_t = inp["floatsize"] == "32" ? Float32 : Float64

if isarg("nullsim")
  pown = parse(Int64, inp["nullsim"])
  estv = readdlm(inp["incov"], ',', float_t)
  mvn_p = MvNormal(estv)
  in_tstats = Array(float_t, pown, size(estv, 1))
  for i in 1:pown
    rand!(mvn_p, view(in_tstats,i, :))
  end
  snpnames = ["sim_$i" for i=1:size(in_tstats, 1)]
  fcov = makecov(in_tstats, "nullsim_$(inp["outcov"])")
else
  snpnames, in_tstats = readtstats(inp["filein"], float_t)
  isarg("replicate") || (fcov = isarg("incov") ? inp["incov"] : makecov(in_tstats, inp["outcov"]))
end

# @everywhere using aspu
@eval @everywhere logB = $logB
@eval @everywhere floatsize = $(inp["floatsize"])
@everywhere float_t = floatsize == "32" ? Float32 : Float64

#Setup aspu objects
@sync @everywhere runvals = Aspuvals(
    zeros(UInt32,2, round(Int64, 10^logB)),
    Matrix{float_t}(9, round(Int64, 10^logB)),
    ones(9)
)

if haskey(inp, "replicate")
  @everywhere estv = eye(float_t, 2)
  @eval @everywhere repfile = $(inp["replicate"])
else
  @eval @everywhere fcov = $fcov
  @everywhere estv = readdlm(fcov, ',', float_t)
end
@everywhere mvn = MvNormal(estv)
@sync @everywhere thisrun = Aspurun(logB, mvn, 3)


#Setup replication
isarg("replicate") && @everywhere rep_setup!(runvals.zb, runvals.rnk, repfile)
isarg("keepzb") && @everywhere rep_setup!(runvals.zb, runvals.rnk, thisrun.mvn)
println("\n\n$(dtnow()): Preprocessing complete")

if isarg("norun"); println("All ready to go! Remove norun option to launch for good\n"); exit(); end

#Chunk input (to do: maybe send directly to workers instead. 10K seems sweet spot for aspu on killdevil)
tstats_chunks = chunkify(in_tstats, nworkers()*10000)


#Run models
println("$(dtnow()): ASPU started on $(nworkers()) workers")
if haskey(inp, "keepzb") || haskey(inp, "replicate")
  @time out = pmap_msg(x->runsnp_rep!(x,thisrun,runvals), tstats_chunks)
else
  @time out = pmap_msg(x->runsnp!(x,thisrun,runvals), tstats_chunks)
end
println("$(dtnow()): ASPU run finished. Writing to file...")


#Write results
fout = open(inp["fileout"], "w")
join(fout, ["snpid", "aspu_p", "gamma"], ',')
snp_i = 1
for res in out
  for i in 1:size(res,1)
    write(fout, "\n$(snpnames[snp_i]),$(join(res[i], ','))")
    snp_i = snp_i+1
  end
end
println("\n$(dtnow()): Writing complete.\n\nJob done.\n\n")
