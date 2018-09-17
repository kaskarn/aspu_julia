using Distributions, Distributed, CSV, DelimitedFiles, Random, ClusterManagers, Dates

sdir = ARGS[1]
sargs = ARGS[2:end]
include("$(sdir)/julia/aspu_utils_io.jl")


using Main.aspu_utils
isarg(a, d=inp) = in(a, keys(d))

#Parse options
#for testing:
if length(ARGS) == 0
  tmarg = ["--filein /proj/epi/CVDGeneNas/antoine/bin/aspu_julia/tests/inputs/testfin.txt --incov /proj/epi/CVDGeneNas/antoine/bin/aspu_julia/tests/inputs/testcov.txt --logB 6 --ncpu 2 --floatsize 64"]
  inp = parse_aspu(tmarg)
else
  inp = parse_aspu(sargs) 
end

#Start workers
np=parse(Int, inp["ncpu"])
if (np > 1)
  if haskey(ENV, "LSB_HOSTS")
    addprocs(split(ENV["LSB_HOSTS"])[2:end]) #keep master on own slot?
  elseif haskey(ENV, "SLURM_JOB_NODELIST")
    addprocs(SlurmManager(np))
  else
    addprocs(np)
  end
end

#Load modules & packages
@eval @everywhere sdir = $sdir
@everywhere include("$(sdir)/julia/aspu_utils_alt.jl")
@everywhere using Distributions, Random, DelimitedFiles
@everywhere using Main.aspu_module

println("\nInputs:"); display(inp); println("\n")

#Set program parameters
logB = isarg("replicate") ? log10(nlines(inp["replicate"])) : parse(Int, inp["logB"])
float_t = inp["floatsize"] == "32" ? Float32 : Float64
pows=collect(1:9)

#Make covariance matrix
fin=open(inp["filein"],"r")
headline=split(readline(fin), ',')
# isarg("replicate") || (fcov = isarg("incov") ? inp["incov"] : makecov(in_tstats, inp["outcov"]))
fcov = inp["incov"]

#Create constants
ntraits = length(headline)-1
@eval @everywhere ntraits = $ntraits
@eval @everywhere logB = $logB
@eval @everywhere logB = Int(logB)
@eval @everywhere pows = $pows
@eval @everywhere floatsize = $(inp["floatsize"])
@everywhere float_t = floatsize == "32" ? Float32 : Float64

#Setup aspu objects
## aspuvals
@everywhere B0 = min(10^MINB, Int(floor((10^logB))))
@everywhere const runvals = Aspuvals{float_t}(
    Int(floor(logB)),
    zeros(Int64, 2, Int(round(B0, digits = 0))), #rank
    zeros(Int64, 2, Int(round(B0, digits = 0))), #rank static,
    zeros(Int64, min(logB,MINB) - 2, 2, round(Int64, B0)), #all ranks,
    
    ones(length(pows)), #pvals
    Matrix{float_t}(undef, ntraits, Int(round(B0, digits = 0))), #zb
    
    Array{float_t}(undef, length(pows)), #lowmin
    Array{float_t}(undef, length(pows), B0), #A0
    Array{float_t}(undef, length(pows), B0), #Astk
)
## aspurun
@eval @everywhere fcov = $fcov
@everywhere estv = readdlm(fcov, ',', float_t)

@everywhere mvn = MvNormal(estv)
@sync @everywhere thisrun = Aspurun(logB, mvn, 3, pows)
@sync @everywhere init_spus!(runvals, pows, mvn, Int(10^logB))

#Preprocessing complete
println("\n\n$(dtnow()): Preprocessing complete")
if isarg("norun"); print("All ready to go! Remove norun option to launch for good\n"); exit(); end

#Output function and file
fout = open(inp["fileout"], "w")
write(fout, "snpid,aspu_p")
for i in pows
  write(fout, ",pval_$(i)")
end
write(fout, ",gamma", '\n')
flush(fout)

out = runsnp!(readline(fin), thisrun, runvals);

#Setup channel
const buffer_s = 2*np;
const jobs = RemoteChannel(()->Channel{String}(buffer_s));
const results = RemoteChannel(()->Channel{typeof(out)}(buffer_s));

#Dispatchers
@everywhere function do_work(jobs, results)
  while true
    snp = take!(jobs)
    out = runsnp!(snp, thisrun, runvals)
    put!(results, out)
  end
end;
function make_jobs(n, io)
  for i in 1:n
    put!(jobs, readline(io))
  end
end;

#Start workers
close(fin); fin=open(inp["filein"],"r"); readline(fin);
@async make_jobs(buffer_s, fin)
for p in workers()
  remote_do(do_work, p, jobs, results)
end;

#rock and roll
@time @sync begin
  #bulk
  while true
    for k in 1:10^5
      eof(fin) && break
      put!(jobs, readline(fin))
      write_aspu(fout, take!(results))
    end
    eof(fin) && break
    println("$(Dates.format(now(), "Yud_HhMM:SS")) -- 100,000 SNPs processed")
    flush(fout)
    flush(stdout)
  end
  sleep(60)
  #flush buffer
  for i in 1:buffer_s
    write_aspu(fout, take!(results))
  end
end

flush(fout)
close(fout)

## all done!!




#old parser
# if isarg("pows")
  # for i in eval(parse(Int, inp["pows"]))
    # x = collect(i)
    # pows=[pows..., x...]
  # end
# else


#Run models
# printlog(mylog, "$(dtnow()): ASPU started on $(nworkers()) workers\n")
# @time pmap_msg(x->runsnp!(x,thisrun,runvals), tstats_chunks, mylog, fout, snpnames_chunks, nprocs())
# @time pmap(x->runsnp!(x,thisrun,runvals),tstats_chunks)

# fin=open(inp["filein"],"r")
# readline(fin)
# @time runsnp!(readline(fin),thisrun,runvals)
# @time pmap(x->runsnp!(x,thisrun,runvals), eachline(fin))
# @time map(x->runsnp!(x,thisrun,runvals), eachline(fin))


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