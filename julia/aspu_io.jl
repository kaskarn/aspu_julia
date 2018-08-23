using Distributions, Distributed, CSV, DelimitedFiles, Random, ClusterManagers

include("/proj/epi/CVDGeneNas/antoine/bin/aspu_julia/julia/aspu_utils_io.jl")
using Main.aspu_utils
mylog = open("aspu_log.txt", "a")

#Parse options
#for testing:
if length(ARGS) == 0
  tmarg = ["--filein /proj/epi/CVDGeneNas/antoine/bin/aspu_julia/tests/inputs/testfin.txt --incov /proj/epi/CVDGeneNas/antoine/bin/aspu_julia/tests/inputs/testcov.txt --logB 6 --ncpu 1 --floatsize 64"]
  inp = parse_aspu(tmarg)
else
  inp = parse_aspu(ARGS) 
end
isarg(a, d=inp) = in(a, keys(d))

#Start workers
np=parse(Int, inp["ncpu"])
if (np > 1)
  if haskey(ENV, "LSB_HOSTS")
    addprocs(split(ENV["LSB_HOSTS"])[2:end]) #keep master on own slot
  elseif haskey(ENV, "SLURM_JOB_NODELIST")
    addprocs(SlurmManager(np),topology=:master_worker)
  end
end

@everywhere include("/proj/epi/CVDGeneNas/antoine/bin/aspu_julia/julia/aspu_utils_alt.jl")
@everywhere using Distributions, Random, DelimitedFiles
@everywhere using Main.aspu_module
println("\nInputs:")
display(inp)
println("\n")
logB = isarg("replicate") ? log10(nlines(inp["replicate"])) : parse(Int, inp["logB"])
float_t = inp["floatsize"] == "32" ? Float32 : Float64
pows = Array{Int64, 1}()

if isarg("pows")
  for i in eval(parse(Int, inp["pows"]))
    x = collect(i)
    pows=[pows..., x...]
  end
else
  pows=collect(1:9)
end

if isarg("nullsim")
  pown = parse(Int64, inp["nullsim"])
  estv = readdlm(inp["incov"], ',', float_t)
  mvn_p = MvNormal(estv)
  in_tstats = Array{float_t}(pown, size(estv, 1))
  for i in 1:pown
    rand!(mvn_p, view(in_tstats,i, :))
  end
  snpnames = ["sim_$i" for i=1:size(in_tstats, 1)]
  fcov = makecov(in_tstats, "nullsim_$(inp["outcov"])")
else
  snpnames, in_tstats = readtstats(inp["filein"], float_t)
  isarg("replicate") || (fcov = isarg("incov") ? inp["incov"] : makecov(in_tstats, inp["outcov"]))
end

ntraits = size(in_tstats, 2)
@eval @everywhere ntraits = $ntraits
@eval @everywhere logB = $logB
@eval @everywhere logB = Int(logB)
@eval @everywhere pows = $pows
@eval @everywhere floatsize = $(inp["floatsize"])
@everywhere float_t = floatsize == "32" ? Float32 : Float64

#Setup aspu objects
@everywhere B0 = min(10^7, Int(floor((10^logB))))
@everywhere runvals = Aspuvals{float_t}(
    Int(floor(logB)),
    zeros(Int64, 2, Int(round(B0, digits = 0))), #rank
    zeros(Int64, 2, Int(round(B0, digits = 0))), #rank static,
    zeros(Int64, min(logB,7) - 2, 2, round(Int64, B0)), #all ranks,
    
    ones(length(pows)), #pvals
    Matrix{float_t}(undef, ntraits, Int(round(B0, digits = 0))), #zb
    
    Array{float_t}(undef, length(pows)), #lowmin
    Array{float_t}(undef, length(pows), B0), #A0
    Array{float_t}(undef, length(pows), B0), #Astk
)

if haskey(inp, "replicate")
  @everywhere estv = eye(float_t, 2)
  @eval @everywhere repfile = $(inp["replicate"])
else
  @eval @everywhere fcov = $fcov
  @everywhere estv = readdlm(fcov, ',', float_t)
end

@everywhere mvn = MvNormal(estv)
@sync @everywhere thisrun = Aspurun(logB, mvn, 3, pows)
@sync @everywhere init_spus!(runvals, pows, mvn, Int(10^logB))

#Setup replication
isarg("replicate") && @everywhere rep_setup!(runvals.zb, runvals.rnk, repfile)
isarg("keepzb") && @everywhere rep_setup!(runvals.zb, runvals.rnk, thisrun.mvn)

#Chunk input (to do: maybe send directly to workers instead. 10K seems sweet spot for aspu on killdevil)
tstats_chunks = chunkify(in_tstats, min(size(in_tstats, 1), nworkers()*10000));
snpnames_chunks = chunkify(snpnames, min(size(in_tstats, 1), nworkers()*10000));

zi = vec(tstats_chunks[1]);

if isarg("norun")
  @time aspu_iter!(pows, 6, zi, mvn, runvals)
  println("")
  @time aspu_iter!(pows, 6, zi, mvn, runvals)
  println("")
  @time aspu_first!(pows, 6, zi, mvn, runvals)
  @time aspu_first!(pows, 6, zi, mvn, runvals)
  @time runsnp!(zi, thisrun, runvals)
  @time runsnp!(zi, thisrun, runvals)
end

printlog(mylog, "\n\n$(dtnow()): Preprocessing complete\n")
if isarg("norun"); print("All ready to go! Remove norun option to launch for good\n"); exit(); end

#Prepare output file
fout = open(inp["fileout"], "w")
write(fout, "chunk,index,snpid,aspu_p")
for i in pows
  write(fout, ",pval_$(i)")
end
write(fout, ",gamma")
flush(fout)

#Run models
printlog(mylog, "$(dtnow()): ASPU started on $(nworkers()) workers\n")
if haskey(inp, "keepzb") || haskey(inp, "replicate")
  @time pmap_msg(x->runsnp_rep!(x,thisrun,runvals), tstats_chunks, mylog, fout, snpnames_chunks, np)
else
  @time pmap_msg(x->runsnp!(x,thisrun,runvals), tstats_chunks, mylog, fout, snpnames_chunks, np)
end
printlog(mylog, "$(dtnow()): ASPU run finished. Job done...\n")
# printlog(mylog, "$(dtnow()): ASPU run finished. Writing to file...\n")

# snp_i = 1
# for res in out
  # for i in 1:length(res)
    # write(fout, "\n$(snpnames[snp_i]),$(join(res[i], ','))")
    # write(fout, "\n$(snpnames[snp_i]),$(res[i][1]),$(join(res[i][2], ',')),$(res[i][3])")
    # snp_i = snp_i+1
  # end
# end
# print("\n$(dtnow()): Writing complete.\n\nJob done.\n\n")

# part=ENV["SLURM_JOB_NODELIST"][1]
    # snode=strip(ENV["SLURM_JOB_NODELIST"],[part, '[', ']'])
    # cpu_string=split(ENV["SLURM_JOB_CPUS_PER_NODE"],',')
    # scpu=Array{Int64,1}()
    # for i in cpu_string
      # if in('x', i)
        # tmp = split(i,'(')
        # n = parse(Int,tmp[1])
        # for x in 1:parse(Int,chop(tmp[2], head=1, tail=1))
          # push!(scpu, n)
        # end
      # else
        # push!(scpu, parse(Int, i))
      # end
    # end

    # machlist = Array{String,1}()
    # local j = 0
    # for i in split(snode, ',')
      # if in('-', i)
        # nstart, nend = split(i, '-')
        # k = nstart
        # while parse(Int,k) <= parse(Int, nend)
          # j = j+1
          # length(k)<length(nstart) && (k=string("0",k))
          # for x in 1:scpu[j]
            # push!(machlist, string(part,k))
          # end
          # k = string(parse(Int, k) + 1)
        # end
      # else
        # j = j+1
        # for x in 1:scpu[j]
          ## addprocs([string("c",k)])
          # push!(machlist, string(part,i))
        # end
      # end
    # end
    # printlog(mylog, "Assigned cores:")
    # printlog(mylog, machlist)
    ## print(machlist)
    # addprocs(machlist[2:end])
  # end