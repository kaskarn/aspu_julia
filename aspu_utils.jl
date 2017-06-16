type Aspuvals
  rnk::Array{UInt32, 2}
  zb::Array
  randspu::Vector
  randval::Vector
  pval::Vector{Int32}
end

type Aspurun
  logB::Real
  mvn::MvNormal
  p0::Int64
end

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
  writecsv("vcov_aspu.txt", estv)
  outname
end

function readtstats(infile, T::DataType)
  in_t = readtable(infile, header=true)
  snpnames = copy(in_t[:,1])
  tstats = convert(Matrix{T}, in_t[:,2:length(in_t)])
  snpnames, tstats
end

# modified from: https://github.com/JuliaLang/julia/issues/939
# http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Julia
function InsertionSort!(A, order, ii=1, jj=length(A))
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
function quicksort!(A, order, i=1, j=length(A))
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


function getspu!(spu, z, n::Int64)
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

  n = length(mvn)
  zi_spu = getspu(t_in, n)

  tmval = view(x.zb, 1:n, :)
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
  x.pval += (x.zb[:,B] .> zi_spu)

  aspu_gamma = sortperm(x.pval)[1]
  minp = x.pval[aspu_gamma]/(B+1)
  minp, aspu_gamma
end

function rank_spus!{T<:Real}(rnk::Array{UInt32, 2},
  zb::Array{T,2}, B::Int)

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

function aspu!{T<:Real}(B::Int64, t_in::Array{T}, mvn, x::Aspuvals)
  fill!(x.pval, 1)

  @inbounds @fastmath minp, aspu_gamma = calc_spus!(x, t_in, mvn, B)
  @inbounds rank_spus!(x.rnk, x.zb, B)

  aspu = 1
  @simd for i in 1:B
    @inbounds ((1 + B - x.rnk[2,i])/B) < minp && (aspu += 1)
  end
  aspu/(B+1), aspu_gamma, minp
end

parse32(x) = parse(Float32, x)
function parsesnp(snpnow)
  #parse line from file. Must start with snp name.
  tm = split(snpnow, ',')
  n = length(tm)
  return tm[1], broadcast(parse32, tm[2:n])
end

function runsnp!{T<:Real}(zi::Array{T, 1}, r::Aspurun, x::Aspuvals, bmax = Inf)
  p0, mvn = r.p0, r.mvn
  logB = min(r.logB, bmax)
  aspu = 0
  aspu_gamma = 0
  while (p0 <= logB) && (aspu < 15/(10^(p0)))
    p0 > logB && (p0 = logB)
    aspu, aspu_gamma = aspu!(ceil(Int,10^p0), zi, mvn, x)
    p0 += 1
  end
  return aspu, aspu_gamma
end

function runsnp!(snp::AbstractString, r::Aspurun, x::Aspuvals, bmax = Inf)
  snpname, zi = parsesnp(snp)
  aspu, aspu_gamma = runsnp!(zi, r, x, bmax)
  return snpname, aspu, aspu_gamma
end

function runsnp!{T<:Real}(snpmat::Matrix{T}, r::Aspurun, x::Aspuvals, bmax = Inf)
  return [runsnp!(snpmat[i,:],r, x, bmax) for i in 1:size(snpmat,1)]
end

# function runsnp!(snp::AbstractString,r::Aspurun, x::Aspuvals, bmax = Inf)
#   p0, mvn = r.p0, r.mvn
#   logB = min(r.logB, bmax)
#   snpname, zi = parsesnp(snp)
#   aspu = 0
#   aspu_gamma = 0
#   while (p0 <= logB) && (aspu < 15/(10^(p0)))
#     p0 > logB && (p0 = logB)
#     aspu, aspu_gamma = aspu!(ceil(Int,10^p0), zi, mvn, x)
#     p0 += 1
#   end
#   return snpname, aspu, aspu_gamma
# end












#
