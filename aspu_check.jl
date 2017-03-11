#### FREEZE ####
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
function aspu_Rtrans!(B, Zi, Z0, zb, n)
  Zi_spu = getspu(Zi,n)
  spu = zeros(Float32, 9)
  pval = ones(9)
  tm = zeros(Float32, n)
  for i = 1:size(Z0,2)
    rand!(mvn, tm)
    getspu!(spu, Z0[:,i], n)
    zb[:,i] = spu
    for j = 1:9
      zb[j,i] > Zi_spu[j] && (pval[j] += 1)
    end
  end
  minp = minimum(pval)/(B+1)
  p0 = Array(Float64,9,B)
  for i in 1:9
    p0[i,:] = (1+B-tiedrank(zb[i,1:B]))/B
  end
  p0m = minimum(p0, 1)
  aspu = (1+sum(p0m .< minp))/(B+1)
  #aspu,p0,p0m,minp
  aspu
end

### FREEZE ####
function aspu_check(B, Zi, Z0, n)
  zb = Array(Float32, 9, ceil(Int, B))
  Zi_spu = getspu(Zi,n)
  spu = zeros(Float32, 9)
  pval = ones(9)
  tm = zeros(Float32, n)
  rnk = Array(Float64, 2, B)
  for i = 1:size(Z0,2)
    getspu!(spu, Z0[:,i], n)
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

R"""
screenaSPU = function(Zi, B, Z0){
   pval = numeric(9)
   for(k in 1:8){
     z0b = abs(rowSums(Z0^k))
     pval[k] = (1+sum(z0b > abs(sum(Zi^k))))/(B+1)
   }
   ## SPU(max)
   z0 = apply(abs(Z0), 1, max)
   pval[9] = (1+sum(z0 > max(abs(Zi))))/(B+1)
   ## aSPU

   p1m = min(pval)
   p0 = matrix(NA, B, 9)
   for(k in 1:8){
     zb = abs(rowSums(Z0^k))
     p0[,k] = (1+B-rank(abs(zb)))/B
   }
   zb = apply(abs(Z0), 1, max)
   p0[,9] = (1+B-rank(zb))/B

   p0m = apply(p0, 1, min)
   k = sum(p0m < min(pval))
   pval[10] = (1+sum(p0m < min(pval) ))/(B+1)
   names(pval) <- c(paste('SPU', 1:9, sep=''), 'aSPU')
   return(pval[10])
}
"""
