using Distributions, DataFrames, Base.Dates
include("./aspu_utils.jl")

#Parse options
#ARGS = ["--filein aspu_test.csv"] #for testing

incmd = join(ARGS, " ")
pinputs = Dict()
[get!(pinputs, split(i)[1], split(i)[2]) for i in split(incmd, "--")[2:end]]

#Set Floating point width
haskey(pinputs, "floatsize") || get!(pinputs, "floatsize", "32")
if !in(pinputs["floatsize"], ["32", "64"]); println("Error: floatsize must be 32 or 64\n"); exit(); end
float_t = pinputs["floatsize"] == "32" ? Float32 : Float64

#Set input file path
if !haskey(pinputs, "filein"); println("Error: no input file defined with filein\n"); exit(); end
if !isfile(pinputs["filein"]); println("Error: input file not found\n"); exit(); end
snpnames, in_tstats = readtstats(pinputs["filein"], float_t)

#Set covariance
haskey(pinputs, "incov") && if haskey(pinputs, "outcov"); println("Error: use either incov or outcov\n"); exit(); end
haskey(pinputs, "incov") || haskey(pinputs, "outcov") || get!(pinputs, "outcov", "vcov_aspu.txt")
fcov = haskey(pinputs, "incov") ? pinputs["incov"] : makecov(in_tstats, pinputs["outcov"])

#Set logB
haskey(pinputs, "logB") || get!(pinputs, "logB", "5")
logB = parse(pinputs["logB"])

#Set outname
outname = haskey(pinputs, "fileout") ? pinputs["fileout"] : "aspu_results_$(Dates.format(now(), "Yud_HhMM")).csv"

println("\nInputs:\n")
@show pinputs

#Start workers
if haskey(ENV, "LSB_HOSTS")
  addprocs(split(ENV["LSB_HOSTS"])[2:end])
else
  addprocs(Sys.CPU_CORES-1)
end

@everywhere using Distributions, DataFrames
@everywhere include("./aspu_utils.jl")

@eval @everywhere fcov = $fcov
@eval @everywhere logB = $logB
@eval @everywhere floatsize = $(pinputs["floatsize"])
@everywhere float_t = floatsize == "32" ? Float32 : Float64

#Setup aspu objects
@sync @everywhere begin
  estv = readdlm(fcov, ',', float_t)
  mvn = MvNormal(estv)
  thisrun = Aspurun(logB, mvn, 3)
  runvals = Aspuvals(
    zeros(UInt32,2, ceil(Int64, 10^logB)),
    Matrix{float_t}(9, ceil(Int64, 10^logB)),
    zeros(float_t, 9),
    zeros(float_t, size(estv,2)),
    ones(9)
  )
end

if haskey(pinputs, "norun"); prinln("All ready to go! Remove norun option to launch for good\n"); exit(); end

#Chunk input (to do: maybe send directly to workers instead)
tstats_chunks = chunkify(in_tstats, nworkers()*10)

#Run models
println("$(Dates.format(now(), "Yud_HhMM")): ASPU started on $(nworkers()) workers")
@time out = pmap(x->runsnp!(x,thisrun,runvals), tstats_chunks)
println("$(Dates.format(now(), "Yud_HhMM")): ASPU run finished. Writing to file...")

fout = open(outname, "w")
join(fout, ["snpid", "aspu_p", "gamma"], ',')
snp_i = 1
for res in out
  for i in 1:size(res,1)
    write(fout, "\n$(snpnames[snp_i]),$(join(res[i], ','))")
    snp_i = snp_i+1
  end
end
println("\nWriting complete.\n\nJob done.\n\n")

exit()
