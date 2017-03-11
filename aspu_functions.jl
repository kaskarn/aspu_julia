type aspurun
  B::Int64
  mvn::MvNormal
  p0::Int
end

#read significance threshold, and input filename
function makecov(infile)
  # Compute covariance matrix and write to disk
  in_res = readtable(infile, header=true)
  nsnp = size(in_res)[1]
  in_mat = convert(Matrix{Float32}, in_res[:,2:length(in_res)])
  ntraits = size(in_mat)[2]

  minZ = minimum(2*cdf(Normal(), -abs(in_mat)), 2)
  nullsnps = minZ[:] .> 0.05/(nsnp/ntraits)
  estv = cor(in_mat[nullsnps,:])
  writecsv("vcov_aspu.txt", estv)
  return nsnp, ntraits
end

function parsesnp(snpnow)
  #parse line from file. Must start with snp name.
  tm = split(snpnow, ',')
  n = length(tm)
  return tm[1], broadcast(parse, tm[2:n])
end

function getspu!(spu::Vector{Float32}, z, n::Int64)
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
end

function getspu(z, n::Int64)
  spu = Array(Float32, 9)
  getspu!(spu,z,n)
  return spu
end

function aspu!(B, Zi, mvn, zb, n)
  Zi_spu = getspu(Zi,n)
  spu = zeros(Float32, 9)
  pval = ones(9)
  tm = zeros(Float32, n)
  rnk = Array(Float64, 2, B)
  for i = 1:B
    rand!(mvn, tm)
    getspu!(spu, tm, n)
    zb[:,i] = spu
    for j = 1:9
      zb[j,i] > Zi_spu[j] && (pval[j] += 1)
    end
  end
  minp = minimum(pval)/(B+1)
  rnk[1,:] = tiedrank(zb[1,1:B])
  for i = 2:9
    rnk[2,:] = tiedrank(zb[i,1:B])
    for j = 1:B
      rnk[2,j] > rnk[1,j] && (rnk[1,j] = rnk[2,j])
    end
  end
  for i = 1:B
    rnk[1,i] = (1 + B - rnk[1,i])/B
  end
  aspu = 1
  for i = 1:B
    rnk[1,i] < minp && (aspu += 1)
  end
  aspu/(B+1) #,zb[1,:],minp
end

function runsnp!(snp,r::aspurun,zb)
  snpname, Zi = parsesnp(snp)
  p0 = r.p0 - 1
  aspu = 0
  while (p0 < r.B) && (aspu < 15/(10^p0))
    p0 += 1
    p0 > r.B && (p0 = r.B)
    aspu = aspu!(ceil(Int,10^p0), Zi, r.mvn, zb, length(r.mvn.μ))
  end
  return snpname, aspu
end

function runsnp_steps!(snp,r::aspurun,zb)
  snpname, Zi = parsesnp(snp)
  p0 = r.p0
  n = length(r.mvn.μ)
  aspu = aspu!(10^p0, Zi, r.mvn, zb, n)
  while (p0 < r.B) && (aspu < 15)
    p0 += 1
    p0 > r.B && (p0 = r.B)
    aspu += getaspu!(ceil(Int,10^p0)-10^ceil(Int, p0-1), Zi, r.mvn, zb, n) - 1
  end
  return snpname, aspu/(10^p0+1)
end
