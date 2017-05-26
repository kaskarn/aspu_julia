type Aspuvals
  rnk::Array{UInt32, 2}
  zb::Array{Float32}
  randspu::Vector{Float32}
  randval::Vector{Float32}
  pval::Vector{Int32}
end

type Aspurun
  logB::Real
  mvn::MvNormal
  p0::Int64
end


function makecov(infile)
  # Compute covariance matrix and write to disk
  in_res = readtable(infile, header=true)
  nsnp = size(in_res)[1]
  in_mat = convert(Matrix{Float64}, in_res[:,2:length(in_res)])
  ntraits = size(in_mat)[2]

  minZ = minimum(2*cdf(Normal(), -abs(in_mat)), 2)
  nullsnps = minZ[:] .> 0.05/(nsnp/ntraits)
  estv = cor(in_mat[nullsnps,:])
  writecsv("vcov_aspu.txt", estv)
  return nsnp, ntraits
end

# Written by rShekhtman
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

# Written by rShekhtman
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

function getspu(z::Vector{Float32}, n::Int64)
  tmpspu = Array(Float32, 9)
  getspu!(tmpspu,z,n)
  tmpspu
end

function calc_spus!(x::Aspuvals, t_in::Vector{Float32},
  mvn::MvNormal, B::Int64)

  n = size(x.randval,1)
  zi_spu = getspu(t_in, n)

  tmval = view(x.zb, 1:5, :)
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

function rank_spus!(rnk::Array{UInt32, 2},
  zb::Array{Float32,2}, B::Int)

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

function aspu!(B::Int64, t_in, mvn, x::Aspuvals)
  fill!(x.pval, 1)

  @time @inbounds @fastmath minp, aspu_gamma = calc_spus!(x, t_in, mvn, B)
  @time @inbounds rank_spus!(x.rnk, x.zb, B)

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

function runsnp!(snp,r::Aspurun,x::Aspuvals, bmax = Inf)
  p0, mvn = r.p0, r.mvn
  logB = min(r.logB, bmax)
  snpname, zi = parsesnp(snp)
  aspu = 0
  aspu_gamma = 0
  while (p0 <= logB) && (aspu < 15/(10^(p0)))
    p0 > logB && (p0 = logB)
    aspu, aspu_gamma = aspu!(ceil(Int,10^p0), zi, mvn, x)
    p0 += 1
  end
  return snpname, aspu, aspu_gamma
end
