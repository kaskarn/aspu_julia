## fix timeout issue in ClusterManagers package
module aspu_utils

using Distributions, Dates, CSV, DelimitedFiles, Distributed, Random

export
  readtstats, makecov, chunkify,
  pmap_msg, dtnow, parse_aspu,
  nline, printlog

dtnow() = Dates.format(now(), "Yud_HhMM")

function printlog(io, msg)
  print(msg)
  (io != stdout) && (print(io, msg))
end

function chunkify(mat, n)
  r = size(mat,1)
  chunksize = vcat(0, cumsum(fill(div(r, n), n) .+ (1:n .<= (r % n))))
  s_inds = [ar[1]:ar[2] for ar in zip(chunksize[1:n] .+ 1, chunksize[2:end])]
  tstats_chunks = [mat[ind,:] for ind in s_inds]
  tstats_chunks
end

function makecov(mat, outname="vcov_aspu.txt")
  minZ = [minimum(map(x->2*cdf(Normal(),-abs(x)), mat[i,:])) for i in 1:size(mat,1)]
  #nullsnps = minZ[:] .> 0.05/(size(mat,1)/size(mat,2))
  nullsnps = minZ .> 10e-5
  estv = cor(mat[nullsnps,:])
#  estv = cov(mat[nullsnps,:])
  writedlm(outname, estv, ',')
  outname
end

function readtstats(infile, T::DataType)
  # in_t = readtable(infile, header=true)
  in_t = CSV.read(infile)
  snpnames = copy(in_t[:,1])
  tstats = convert(Matrix{T}, in_t[:,2:length(in_t)])
  snpnames, tstats
end
 
function pmap_msg(f, lst, logf, fout, snpnames, np)
    n = length(lst)
    i = 1
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > n
                            break
                        end
                        if (idx % (div(n,100))) == 0 
                          printlog(logf, "$(dtnow()): Chunk $(idx) of $(n) started!\n")
                          flush(logf)
                          flush(fout)
                        end
                        res=remotecall_fetch(f, p, lst[idx])
                        for j in 1:length(res)
                          write(fout, "\n$(idx),$(j),$(snpnames[idx][j]),$(res[j][1]),$(join(res[j][2], ',')),$(res[j][3])")
                        end
                    end
                end
            end
        end
    end
end

function pathcheck(p)
  isfile(p) || print("ERROR: file $(p) not found")
  1-isfile(p)
end

nand(a,b) = !(a & b)
checkargs(inp, a, b, f) = f(in(a,keys(inp)), in(b,keys(inp)))
argnor(inp, a, b) = !(in(a,keys(inp) | in(b,keys(inp))))
function argexcl(inp, a, b)
  iserr = in(a,keys(inp)) & in(b,keys(inp))
  iserr && print("ERROR: Cannot use both --$(a) and --$(b)")
  iserr
end

function argeither(inp, a, b)
  iserr = !xor(in(a,keys(inp)), in(b,keys(inp)))
  iserr && print("ERROR: Need either --$(a) or --$(b)")
  iserr
end

function argdep(inp, a, b)
  iserr = in(a, keys(inp)) & !(in(b, keys(inp)))
  iserr && print("ERROR: --$(a) must be used with --$(b)")
  iserr
end

function parse_aspu(argsin)
  tbeg = dtnow()
  incmd = join(argsin, " ")
  pinputs = Dict()
  [get!(pinputs, split(i)[1], size(split(i),1) > 1 ? split(i)[2] : true) for i in split(incmd, "--")[2:end]]
  isarg(a) = in(a, keys(pinputs))
  #Check errors
  haserr = 0
  if isarg("floatsize") && !in(pinputs["floatsize"], ["32", "64"])
    print("ERROR: --floatsize must be 32 or 64\n")
    haserr += 1
  end
  if isarg("nullsim") && !(typeof(parse(Int, pinputs["nullsim"])) <: Int)
    print("ERROR: --nullsim must be an integer number")
    haserr += 1
  end
  haserr += argeither(pinputs, "nullsim", "filein")
  haserr += argexcl(pinputs, "incov", "outcov")
  haserr += argexcl(pinputs, "replicate", "keezb")
  haserr += argexcl(pinputs, "replicate", "logB")
  haserr += argexcl(pinputs, "replicate", "incov")
  haserr += argdep(pinputs, "replicate", "filein")
  isarg("filein") && (haserr += pathcheck(pinputs["filein"]))
  isarg("incov") && (haserr += pathcheck(pinputs["incov"]))
  if haserr > 0;print("There were $(haserr) command-line error(s)\n");exit();end
  #Set defaults
  checkargs(pinputs, "replicate", "logB", |) || get!(pinputs, "logB", "5")
  checkargs(pinputs, "outcov", "replicate", |) || get!(pinputs, "outcov", "vcov_aspu_$(tbeg).txt")
  isarg("floatsize") || get!(pinputs, "floatsize", "32")
  if isarg("replicate")
    dfltname = "aspu_results_$(tbeg)_replicate.csv"
  else
    logB = parse(Int, pinputs["logB"])
    bl, br = round(10^logB/10^floor(logB),digits=2), floor(Int64,logB)
    dfltname = "aspu_results_$(tbeg)_$(bl)E$(br).csv"
  end
  isarg("fileout") || get!(pinputs, "fileout", dfltname)
  return pinputs
end

function nlines(filein)
  zbrep = open(filein, "r")
  zb_n = 0
  while(!eof(zbrep))
    readline(zbrep)
    zb_n += 1
  end
  zb_n
end

end
