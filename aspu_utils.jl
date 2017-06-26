module aspu_utils

using DataFrames, Distributions

export
  readtstats, makecov, chunkify,
  pmap_msg, dtnow, parse_aspu,
  nlines

dtnow() = Dates.format(now(), "Yud_HhMM")

function chunkify(mat, n)
  r = size(mat,1)
  chunksize = vcat(0, cumsum(fill(div(r, n), n) .+ (1:n .<= (r % n))))
  s_inds = [ar[1]:ar[2] for ar in zip(chunksize[1:n]+1, chunksize[2:end])]
  tstats_chunks = [mat[ind,:] for ind in s_inds]
  tstats_chunks
end

function makecov(mat, outname="vcov_aspu.txt")
  minZ = minimum(2*cdf(Normal(), -abs(mat)), 2)
  nullsnps = minZ[:] .> 0.05/(size(mat,1)/size(mat,2))
  estv = cor(mat[nullsnps,:])
  writecsv(outname, estv)
  outname
end

function readtstats(infile, T::DataType)
  in_t = readtable(infile, header=true)
  snpnames = copy(in_t[:,1])
  tstats = convert(Matrix{T}, in_t[:,2:length(in_t)])
  snpnames, tstats
end

function pmap_msg(f, lst)
    np = nprocs()
    n = length(lst)
    results = Vector{Any}(n)
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
                        (idx % (div(n,100))) == 0 && (println("$(dtnow()): Chunk $(idx) of $(n) started!"))
                        results[idx] = remotecall_fetch(f, p, lst[idx])
                    end
                end
            end
        end
    end
    results
end

function pathcheck(p)
  isfile(p) || println("ERROR: file $(p) not found")
  1-isfile(p)
end
nand(a,b) = !(a & b)
checkargs(inp, a, b, f) = f(in(a,keys(inp)), in(b,keys(inp)))
argnor(inp, a, b) = !(in(a,keys(inp) | in(b,keys(inp))))

function argexcl(inp, a, b)
  iserr = in(a,keys(inp)) & in(b,keys(inp))
  iserr && println("ERROR: Cannot use both --$(a) and --$(b)")
  iserr
end
function argeither(inp, a, b)
  iserr = !(in(a,keys(inp)) $ in(b,keys(inp)))
  iserr && println("ERROR: Need either --$(a) or --$(b)")
  iserr
end
function argdep(inp, a, b)
  iserr = in(a, keys(inp)) & !(in(b, keys(inp)))
  iserr && println("ERROR: --$(a) must be used with --$(b)")
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
    println("ERROR: --floatsize must be 32 or 64\n")
    haserr += 1
  end
  if isarg("nullsim") && !(typeof(parse(pinputs["nullsim"])) <: Int)
    println("ERROR: --nullsim must be an integer number")
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
  if haserr > 0;println("There were $(haserr) command-line error(s)\n");exit();end

  #Set defaults
  checkargs(pinputs, "replicate", "logB", |) || get!(pinputs, "logB", "5")
  checkargs(pinputs, "outcov", "replicate", |) || get!(pinputs, "outcov", "vcov_aspu_$(tbeg).txt")
  isarg("floatsize") || get!(pinputs, "floatsize", "32")

  if isarg("replicate")
    dfltname = "aspu_results_$(tbeg)_replicate.csv"
  else
    logB = parse(pinputs["logB"])
    bl, br = round(10^logB/10^floor(logB),2), floor(Int64,logB)
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

end #end of module



module aspu

using Distributions

export
  Aspuvals, Aspurun,
  runsnp!, runsnp_rep!, getspu!, aspu!,
  rep_setup!

type Aspuvals
  rnk::Array{UInt32, 2}
  zb::Array
  pval::Vector{UInt32}
end

type Aspurun
  logB::Int64
  mvn::MvNormal
  p0::Int64
end

# modified from: https://github.com/JuliaLang/julia/issues/939
# http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Julia
function InsertionSort!{T<:Real}(A::AbstractArray{T, 1},
  order::AbstractArray{UInt32, 1}, ii=1, jj=length(A))

    for i = ii+1 : jj
        j = i - 1
        temp  = A[i]
        itemp = order[i]

        while true
            if j == ii-1
                break
            end
            if A[j] <= temp
                break
            end
            A[j+1] = A[j]
            order[j+1] = order[j]
            j -= 1
        end

        A[j+1] = temp
        order[j+1] = itemp
    end  # i
    return A
end # function InsertionSort!

# modified from: https://github.com/JuliaLang/julia/issues/939
# http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Julia
function quicksort!{T<:Real}(A::AbstractArray{T, 1},
  order::AbstractArray{UInt32, 1}, i=1, j=length(A))

    if j > i
      if  j - i <= 10
        # Insertion sort for small groups is faster than Quicksort
        InsertionSort!(A,order, i,j)
        return A
      end
      pivot = A[ div(i+j,2) ]
      left, right = i, j
      while left <= right
        while A[left] < pivot
            left += 1
        end
        while A[right] > pivot
            right -= 1
        end
        if left <= right
            A[left], A[right] = A[right], A[left]
            order[left], order[right] = order[right], order[left]
            left += 1
            right -= 1
        end
      end  # left <= right
      quicksort!(A,order, i, right)
      quicksort!(A,order, left, j)
    end  # j > i
    return A
end # function quicksort!


function getspu!{T<:Real}(spu::AbstractArray{T, 1},
  z::AbstractArray{T, 1}, n::Int64)

  for i = 1:8
    @inbounds spu[i] = z[1]^i
  end
  @inbounds spu[9] = abs(z[1])
  for j = 2:n
    for i = 1:8
      @inbounds spu[i] += z[j]^i
    end
    @inbounds spu[9] < abs(z[j]) && (spu[9] = abs(z[j]))
  end
  for i = 1:8
    @inbounds spu[i] = abs(spu[i])
  end
  0
end

function getspu{T<:Real}(z::Vector{T}, n::Int64)
  tmpspu = Array(T, 9)
  getspu!(tmpspu,z,n)
  tmpspu
end

function calc_spus!{T<:Real}(x::Aspuvals, t_in::Vector{T},
  mvn::MvNormal, B::Int64)

  fill!(x.pval, 1)
  n = length(mvn)
  zi_spu = getspu(t_in, n)

  tmval = view(x.zb, 1:n, 1:B)
  rand!(mvn, tmval)
  firstline = tmval[:,1]
  for i in 2:B
    zbnow = view(x.zb, :, i-1)
    getspu!(zbnow, tmval[:,i], n)
    for j = 1:9
      zbnow[j] > zi_spu[j] && (x.pval[j] += 1)
    end
  end
  getspu!(view(x.zb, :, B), firstline, n)
  x.pval += (x.zb[:, B] .> zi_spu)

  aspu_gamma = sortperm(x.pval)[1]
  minp = x.pval[aspu_gamma]/(B+1)
  minp, aspu_gamma
end

function rank_spus!{T<:Real}(rnk::Array{UInt32, 2},
  zb::Array{T,2}, B::Int64)

  rnk1_v = view(rnk,1,1:B)
  # rnk2_v = view(rnk,2,1:B) slower than getindex()
  # zbnow = view(zb,1, 1:B) slower than getindex()

  for i in 1:B
    rnk1_v[i] = i
  end
  quicksort!(zb[1,1:B], rnk1_v)
  for (i, val) in enumerate(rnk1_v)
    rnk[2,val] = i
    rnk1_v[i] = i
  end

  for i in 2:8
    quicksort!(zb[i,1:B], rnk1_v)
    for (j, val) in enumerate(rnk1_v)
      rnk[2,val] < j && (rnk[2,val] = j)
      rnk1_v[j] = j
    end
  end

  quicksort!(zb[9,1:B], rnk1_v)
  for (j, val) in enumerate(rnk1_v)
    rnk[2,val] < j && (rnk[2,val] = j)
  end

  0
end

function aspu!{T<:Real}(B::Int64, t_in::Array{T},
  mvn::MvNormal, x::Aspuvals)

  @inbounds @fastmath minp, aspu_gamma = calc_spus!(x, t_in, mvn, B)
  @inbounds rank_spus!(x.rnk, x.zb, B)

  aspu_p = 1
  @simd for i in 1:B
    @inbounds ((1 + B - x.rnk[2,i])/B) <= minp && (aspu_p += 1)
  end
  aspu_p/(B+1), aspu_gamma
end

parse32(x) = parse(Float32, x)
function parsesnp(snpnow)
  #parse line from file. Must start with snp name.
  tm = split(snpnow, ',')
  n = length(tm)
  return tm[1], broadcast(parse32, tm[2:n])
end

function runsnp!{T<:Real}(zi::Array{T, 1}, r::Aspurun,
  x::Aspuvals, bmax = Inf)

  p0, mvn = r.p0, r.mvn
  logB = min(r.logB, bmax)
  aspu = 0
  aspu_gamma = 0
  while (p0 <= logB) && (aspu < 15/(10^(p0)))
    p0 > logB && (p0 = logB)
    aspu, aspu_gamma = aspu!(round(Int64,10^p0), zi, mvn, x)
    p0 += 1
  end
  return aspu, aspu_gamma
end

function runsnp!(snp::AbstractString, r::Aspurun, x::Aspuvals,
  bmax = Inf)
  snpname, zi = parsesnp(snp)
  aspu, aspu_gamma = runsnp!(zi, r, x, bmax)
  return snpname, aspu, aspu_gamma
end

function runsnp!{T<:Real}(snpmat::Matrix{T}, r::Aspurun,
  x::Aspuvals, bmax = Inf)
  return [runsnp!(snpmat[i,:], r, x, bmax) for i in 1:size(snpmat, 1)]
end



function aspu!{T<:Real}(B::Int64, t_in::Array{T}, x::Aspuvals)
  fill!(x.pval, 1)
  zi_spu = getspu(t_in, length(t_in))
  for i in 1:B
    for j in 1:9
      @inbounds x.zb[j,i] > zi_spu[j] && (x.pval[j] += 1)
    end
  end
  aspu_gamma = sortperm(x.pval)[1]
  minp = x.pval[aspu_gamma]/(B+1)

  aspu_p = 1
  @simd for i in 1:B
    @inbounds ((1 + B - x.rnk[2,i])/B) <= minp && (aspu_p += 1)
  end
  aspu_p/(B+1), aspu_gamma
end
function runsnp_rep!{T<:Real}(zi::AbstractArray{T, 1}, r::Aspurun,
  x::Aspuvals, bmax = Inf)
  return aspu!(round(Int64, 10^r.logB), zi, x)
end
function runsnp_rep!{T<:Real}(snpmat::AbstractArray{T, 2}, r::Aspurun,
  x::Aspuvals, bmax = Inf)
  return [runsnp_rep!(snpmat[i,:], r, x, bmax) for i in 1:size(snpmat, 1)]
end

function readtrans!{T<:Real, S<:AbstractString}(ar::Array{T, 2}, fname::S)
  #read in z0 to replicate
  filein = open(fname, "r")
  tm = split(readline(filein), ',')
  tmlen = length(tm)
  for j in 1:length(tm)
    ar[j, 1] = parse(T, tm[j])
  end

  i = 1
  while(!eof(filein))
    i += 1
    tm = split(readline(filein), ',')
    for j in 1:tmlen
      ar[j, i] = parse(T, tm[j])
    end
  end
  tmlen
end

function rep_setup!{T<:Real, S<:AbstractString}(ar::Array{T, 2},
  rk::Array{UInt32, 2}, fname::S)
  ntraits = readtrans!(ar, fname)
  calc_once!(ar, rk, ntraits)
end

function rep_setup!{T<:Real}(ar::Array{T, 2}, rk::Array{UInt32, 2}, mvn::MvNormal)
  rand!(mvn, view(ar, 1:length(mvn), :))
  calc_once!(ar, rk, length(mvn))
end

function calc_once!{T<:Real}(ar::Array{T, 2}, rk::Array{UInt32, 2}, n::Int)
  #get z0 spus, and ranks
  tmval = view(ar, 1:n, :)
  B = size(ar, 2)
  firstline = tmval[:,1]
  for i in 2:B
    zbnow = view(ar, :, i-1)
    getspu!(zbnow, tmval[:,i], n)
  end
  getspu!(view(ar, :, B), firstline, n)
  @inbounds rank_spus!(rk, ar, B)
end

end #end of module





#
