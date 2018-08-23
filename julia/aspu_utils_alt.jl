module aspu_module
using Distributions, Random

#Defines Aspuvals and Aspurun types
export
  Aspuvals, Aspurun,
  runsnp!, runsnp_rep!, getspu!, aspu!,
  rep_setup!, calc_spus!, rank_spus!,
  getspu, getspu!, init_spus!, calc_spus_first!,
  calc_spus_iter!, CHUNKN, MAXN, aspu_first!, aspu_iter!, aspu!
  
const CHUNKN = 10000
const MAXN = Int64(1e6)


struct Aspuvals{T<:AbstractFloat}
  maxb::Int64
  rnk::Array{Int64, 2}
  rnk_s::Array{Int64, 2}
  rnk_all::Array{Int64, 3}
  
  pval::Vector{Int64}
  zb::Array{T, 2}
  
  lowmin::Array{T, 1}
  A0::Array{T, 2}
  # A::Array{T, 3}
  # A_top::Array{T, 3}
  Astk::Array{T, 2}
  
  # R::Array{Int64, 2}
  # DUPS::Array{Bool, 2}
end

struct Aspurun
  logB::Int64
  mvn::MvNormal
  p0::Int64
  pows::Array{Int64,1}
end


#Sorting functions 

# modified from: https://github.com/JuliaLang/julia/issues/939
# http://rosettacode.org/wiki/Sorting_algorithms/Quicksort#Julia
function InsertionSort!(A::AbstractArray{T, 1}, order::AbstractArray{Int64, 1}, ii=1, jj=length(A)) where {T<:Real}
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

function quicksort!(A::AbstractArray{T, 1}, order::AbstractArray{Int64, 1}, i=1, j=length(A)) where {T<:Real}
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

#SPU calculating function (optimal memory access pattern)
# function getspu!{T<:Real}(spu::AbstractArray{T,1}, pows::Array{Int64,1}, z::Vector{T}, n::Int64)
function getspu!(spu, pows, z, n)
  for i in eachindex(pows)
    @inbounds pows[i] < 9 && (spu[i] = z[1]^(pows[i]))
    @inbounds pows[i] == 9 && (spu[i] = abs(z[1]))
  end
  for j = 2:n
    for i in eachindex(pows)
      @inbounds (pows[i] < 9) && (spu[i] += z[j]^(pows[i]))
      @inbounds (pows[i] == 9) && (abs(z[j]) > spu[i]) && (spu[i] = abs(z[j]))      
    end
  end
  for i = eachindex(pows)
    @inbounds spu[i] = abs(spu[i])
  end
end
function getspu(pows::Array{Int64, 1}, z::Vector{T}, n::Int64) where {T<:Real}
  tmpspu = Array{T}(undef, length(pows))
  getspu!(tmpspu,pows,z,n)
  tmpspu
end



#Function ranking SPUs to create aSPU
function rank_spus!(rnk::AbstractArray{Int64, 2}, zb::Array{T,2}, B = size(zb, 2)) where {T<:Real} 
  rnk1_v = view(rnk,1,1:B)
  for i in 1:B
    rnk1_v[i] = i
  end
  quicksort!(zb[1,1:B], rnk1_v)
  for (i, val) in enumerate(rnk1_v)
    rnk[2,val] = i
    rnk1_v[i] = i
  end
  for i in 2:(size(zb,1)-1)
    quicksort!(zb[i,1:B], rnk1_v)
    for (j, val) in enumerate(rnk1_v)
      rnk[2,val] < j && (rnk[2,val] = j)
      rnk1_v[j] = j
    end
  end
  quicksort!(zb[size(rnk,1),1:B], rnk1_v)
  for (j, val) in enumerate(rnk1_v)
    rnk[2,val] < j && (rnk[2,val] = j)
  end
  0 
end

#Init
function init_spus!(x::Aspuvals{T}, pows::Array{Int64, 1}, mvn::MvNormal, B::Int64) where {T<:Real}
  B0 = min(10^7, Int(floor(B)))

  n = length(mvn)
  # zb = view(x.zb, :, 1:B0)
  zb = x.zb
  rand!(mvn, zb)
  for i in 1:B0
    zbnow = view(zb, :, i)
    destnow = view(x.A0, :, i)
    getspu!(destnow, pows, zbnow, n)
  end
  
  #x.lowmin = Array{T}(undef, length(pows))
  for i in eachindex(pows)
    tmp_r = sortperm(x.A0[i, :])
    top_r = tmp_r[Int(B0-CHUNKN+1):Int(B0)]
    x.lowmin[i] = minimum(x.A0[i, top_r])
    # x.A_top[i,:,:] = x.A0[:, top_r]
    # for j in 1:size(x.R,2)
      # x.R[i,j] = j
    # end
  end
  
  logB0 = Int(floor(log10(B0)))
  for i in 3:min(logB0, 7)
    tmp = copy(x.A0[:,1:10^i])
    tmpr = view(x.rnk_all, i-2, :, :)
    rank_spus!(tmpr, tmp)
  end
  rank_spus!(x.rnk_s, x.A0)
end

#SPU calculating functions
## Storing all zb
function calc_spus!(x::Aspuvals{T}, pows::Array{Int64, 1}, t_in::Vector{T}, B::Int64) where {T<:Real}
  fill!(x.pval, 1)
  n = length(mvn)
  
  zi_spu = getspu(pows, t_in, n)
  tmval = zeros(T, n, B)
  zb = view(x.zb, 1:B, :)
  rand!(mvn, zb)
  
  for i in 1:B
    zbnow = view(zb, :, i)
    getspu!(zbnow, pows, tmval, n)
    for j in eachindex(pows)
      zbnow[j] > zi_spu[j] && (x.pval[j] += 1)
    end
  end
  
  aspu_gamma = sortperm(x.pval)[1]
  minp = x.pval[aspu_gamma]/(B+1)
  
  minp, x.pval, aspu_gamma
end

function calc_spus_first!(x::Aspuvals{T}, pows::Array{Int64, 1}, t_in::Vector{T}, B) where {T<:Real}
  fill!(x.pval, 1)
  zi_spu = getspu(pows, t_in, length(t_in))
  for i in 1:B
    for j in eachindex(pows)
      x.A0[j, i] > zi_spu[j] && (x.pval[j] += 1)
    end
  end
  aspu_gamma = sortperm(x.pval)[1]
  minp = x.pval[aspu_gamma]/(B+1)
  minp, x.pval./(B+1), aspu_gamma
end

## Storing subset of zb
function calc_spus_iter!(x::Aspuvals{T}, pows::Array{Int64, 1}, t_in::Vector{T}, mvn::MvNormal, B::Int64) where {T<:Real}
  fill!(x.pval, 1)
  
  n = length(mvn)
  # tmval = zeros(T, length(mvn))

  tmspu = zeros(T, length(pows))
  k = 0
  
  B0 = min(size(x.A0,2), Int(floor(B/10)))
  
  zi_spu = getspu(pows, t_in, n)
  for i in 1:B0
    for j in eachindex(pows)
      x.A0[j, i] > zi_spu[j] && (x.pval[j] += 1)
      if x.A0[j, i] > x.lowmin[j]
        k = k + 1
        # x.Astk[:, k] = x.A0[:, i]
        for z in eachindex(pows)
          x.Astk[z, k] = x.A0[z, i]
        end
      end
    end
  end
  
  chunks = div(B,B0)-1
  for chk in chunks
    zbnow = view(x.zb, :, 1:B0)
    rand!(mvn, zbnow)
    for i in 1:B0
      keepsim = true
      # rand!(mvn, tmval)
      getspu!(tmspu, pows, zbnow[:,i], n)

      for j in eachindex(pows)
        tmspu[j] > zi_spu[j] && (x.pval[j] += 1)
        if tmspu[j] > x.lowmin[j] && keepsim
          keepsim = false
          k = k + 1
          # x.Astk[:, k] = tmspu
          for z in eachindex(pows)
            x.Astk[z, k] = tmspu[z]
          end
        end
      end
    end
  end
  min_iter, aspu_gamma = findmin(x.pval)
  minp = min_iter/(B+1)
  
  minp, x.pval./(B+1), aspu_gamma, k
end



#aSPU function
## If storing all zb
function aspu_first!(pows::Array{Int64, 1}, logB::Int64, t_in::Vector{T}, mvn::MvNormal, x::Aspuvals{T}) where {T<:Real}
  B = 10^logB
  @inbounds @fastmath minp, pvals, aspu_gamma = calc_spus_first!(x, pows, t_in, B)
  aspu_p = 1
  ind = logB - 2
  @simd for i in 1:B
    @inbounds ((1 + B - x.rnk_all[ind,2,i])/B) <= minp && (aspu_p += 1)
  end
  aspu_p/(B+1), pvals, aspu_gamma
end


## If storing subset of zb
function aspu_iter!(pows::Array{Int64, 1}, logB::Int64, t_in::Vector{T}, mvn::MvNormal, x::Aspuvals{T}) where {T<:Real}
  B = Int(10^logB)
  @inbounds @fastmath minp, pvals, aspu_gamma, krows = calc_spus_iter!(x, pows, t_in, mvn, B)
  # n_uniq = sum(.!x.DUPS)
  @inbounds rank_spus!(x.rnk, x.Astk, krows)
  aspu_p = 1
  @simd for i in 1:krows
    @inbounds ((1 + krows - x.rnk[2,i])/B) <= minp && (aspu_p += 1)
  end
  aspu_p/(B+1), pvals, aspu_gamma
end

#Run aSPU on vector
function runsnp!(zi::Vector{T}, r::Aspurun, x::Aspuvals{T}, bmax = Inf) where {T<:Real}
  p0, mvn, pows = r.p0, r.mvn, r.pows
  logB = Int64(floor(min(r.logB, bmax)))
  aspu = 0
  aspu_gamma = 0
  pvals = Array{Float64}(undef, length(pows))
  
  while (p0 <= min(logB,7)) && (aspu < 15/(10^(p0)))
    p0 > logB && (p0 = logB)
    aspu, pvals, aspu_gamma = aspu_first!(pows, p0, zi, mvn, x)
    p0 += 1
  end
  
  while (p0 <= logB && (aspu < 15/10^(p0)))
    p0 > logB && (p0 = logB)
    aspu, pvals, aspu_gamma = aspu_iter!(pows, p0, zi, mvn, x)
    p0 += 1
  end
  return aspu, pvals, aspu_gamma
end

#Run aSPU on Matrix
function runsnp!(snpmat::Matrix{T}, r::Aspurun, x::Aspuvals, bmax = Inf) where {T<:Real}
  return [runsnp!(vec(snpmat[i,:]), r, x, bmax) for i in 1:size(snpmat, 1)]
end


#### OLD #####


function calc_spus!(x::Aspuvals{T}, pows::Array{Int64, 1}, t_in::Vector{T}, mvn::MvNormal, B::Int64) where {T<:Real}
  fill!(x.pval, 1)
  n = length(mvn)
  zi_spu = getspu(pows, t_in, n)
  tmval = copy(t_in)
  for i in 1:B
    zbnow = view(x.zb, :, i)
    rand!(mvn, tmval)
    getspu!(zbnow, pows, tmval, n)
    for j in eachindex(pows)
      zbnow[j] > zi_spu[j] && (x.pval[j] += 1)
    end
  end
  #getspu!(view(x.zb, :, B), pows, firstline, n)
  #x.pval += (x.zb[:, B] .> zi_spu)
  aspu_gamma = sortperm(x.pval)[1]
  minp = x.pval[aspu_gamma]/(B+1)
  minp, x.pval, aspu_gamma
end


function aspu!(pows::Array{Int64, 1}, B::Int64, t_in::Array{T}, mvn::MvNormal, x::Aspuvals{T}) where {T<:Real}
  @inbounds @fastmath minp, pvals, aspu_gamma = calc_spus!(x, pows, t_in, mvn, B)
  @inbounds rank_spus!(x.rnk, x.zb, B)
  aspu_p = 1
  @simd for i in 1:B
    @inbounds ((1 + B - x.rnk[2,i])/B) <= minp && (aspu_p += 1)
  end
  aspu_p/(B+1), pvals./(B+1), aspu_gamma
end


end











