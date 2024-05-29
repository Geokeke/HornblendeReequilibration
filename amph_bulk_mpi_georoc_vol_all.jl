using MPI, StatGeochem, DelimitedFiles, DataFrames, CSV

# Absolute paths to perplex resources
# const perplexdir = joinpath(resourcepath,"perplex-6.8.7")
const perplexdir = joinpath(resourcepath,"perplex-stable")
const scratchdir = "/dartfs-hpc/scratch/f004dmk/geothermo/K_solution_georoc_vol_all/" # Location of directory to store output files
const resultdir = "./results/geothermo/hbl_bulk_georoc_vol_all/"
# Attempt to install perplex, if not already extant
if !isfile(joinpath(perplexdir,"vertex"))
    # Make sure resourcepath exists
    run(`mkdir -p $resourcepath`)

    # Try to compile PerpleX from source; if that fails, try to download linux binaries
    try
        # Check if there is a fortran compiler
        run(`gfortran -v`)

        # Download Perplex v6.8.7 -- known to work with interface used here
        file = download("https://storage.googleapis.com/statgeochem/perplex-6.8.7-source.zip", joinpath(resourcepath,"perplex-stable.zip"))

        # # For a more updated perplex version, you might also try
        # file = download("https://petrol.natur.cuni.cz/~ondro/perplex-sources-stable.zip", joinpath(resourcepath,"perplex-stable.zip"))

        run(`unzip -u $file -d $resourcepath`) # Extract
        system("cd $perplexdir; make") # Compile
    catch
        @warn "Failed to compile from source, trying precompiled linux binaries instead"
        run(`mkdir -p $perplexdir`)
        file = download("https://petrol.natur.cuni.cz/~ondro/Perple_X_6.8.7_Linux_64_gfortran.tar.gz","perplex-6.8.7-linux.tar.gz")
        run(`tar -xzf $file -C $perplexdir`)
    end
end

if !isdir(resultdir)
    run(`mkdir -p $resultdir`)
end
if !isdir(scratchdir)
    run(`mkdir -p $scratchdir`)
end

# Start MPI
MPI.Init()

# Run perplex for a given sample and save results
function perplex_geothermo(data::Dict, allmajors::Matrix, elements, P_range, T_surf, sample::Int, sample_set)
    # Emphasis on phases from Green (2016), but only those likely to occur in igneous rocks. Use with hp11ver.dat or newer
    G_solution_phases = "Augite(G)\nOpx(JH)\ncAmph(G)\noAmph(DP)\nO(JH)\nSp(JH)\nGrt(JH)\nMica(W)\nBio(TCC)\nF\n"
    G_excludes ="ged\nfanth\ngl\n" # Use solutions instead of endmembers
    fsp_model = "feldspar_B"; G_solution_phases *= fsp_model*"\n"
    fsp_model = ""

    # # Also exclude zircon since we'll calculate it ourselves instead (and also sphene if using Ayers sphene)
    # excludes *= "zrc\nsph\n"

    # Melt model of Holland, Green, & Powell 2018 (doi: 10.1093/petrology/egy048)
    # Includes Ti, Cr, and Fe3, and covers "Peridotites through to Granites", from 0.001 to 70kbar!
    melt_model = "melt(HGP)"
    thermo_dataset="hp633ver.dat"; dataset_uppercase=false

    # Starting composition
    composition = allmajors[sample,:]

    # # Emphasis on phases from Holland and Powell -- all phases can be used with hp02ver.dat.
    # HP_solution_phases = "Omph(HP)\nOpx(HP)\nGlTrTsPg\nAnth\nO(HP)\nSp(HP)\nGt(HP)\nfeldspar_B\nMica(CF)\nBio(TCC)\nChl(HP)\nCtd(HP)\nSapp(HP)\nSt(HP)\nIlHm(A)\nDo(HP)\nT\nB\nF\n"
    # HP_excludes = ""

    # # Emphasis on magmatic phases from Holland and Powell, plus good Amph
    # solution_phases = "Omph(HP)\nOpx(HP)\ncAmph(G)\nAnth\nO(HP)\nSp(HP)\nGt(HP)\nMica(CF)\nBio(TCC)\nChl(HP)\nIlHm(A)\nF\n"
    # excludes = ""

    # Modified from Jeff Moyen's Excel spreadsheet, replacing feldspar
    K_solution_phases = "feldspar_B\nSp(HGP)\nGt(HGP)\nO(HGP)\nOpx(HGP)\nCpx(HGP)\nCrd(HGP)\nBi(HGP)\nMica(W)\nEp(HP11)\ncAmph(G)\nIlm(WPH)\nChl(W)\n"
    excludes = ""

    # ---
    # Configure (run build and vertex)
    @time perplex_configure_geotherm(perplexdir, scratchdir, composition, elements, P_range, T_surf, 0.04,
            dataset=thermo_dataset,
            excludes=excludes,
            solution_phases=melt_model*"\n"*K_solution_phases,
            index=sample)

    ## ---
    modes = perplex_query_modes(perplexdir, scratchdir, index=sample)

    calc = perplex_query_system(perplexdir, scratchdir, index=sample)
    melts = perplex_query_phase(perplexdir, scratchdir, "melt(HGP)", index=sample)
    hbl = perplex_query_phase(perplexdir, scratchdir, "cAmph(G)", index=sample)

    # find system bulk composition
    sys_comp = Dict{String, Array{Float64}}()
    for e in elements
        sys_comp[e] = calc[e]
    end

    # find hbl composition
    hbl_comp = Dict{String, Array{Float64}}()
    for e in elements
        hbl_comp[e] = hbl[e]
    end

    # query solid
    bulk_list = modes["elements"]
    # deleteat!(sol_list, findall(x->(x == "P(bar)" || x ==  "T(K)" || x ==  "melt(HGP)"),sol_list))
    deleteat!(bulk_list, findall(x->(x == "P(bar)" || x ==  "T(K)"), bulk_list))
    bulk_modes = Dict{String, Array{Float64}}()
    for b in bulk_list
        bulk_modes[b] = modes[b]
    end

    tot_vol = vec(nansum(Matrix(DataFrame(bulk_modes)), dims = 2))

    melts_dic = Dict{String, Array{Float64}}()
    for e in elements
        melts_dic[e] = melts[e]
    end
    melts_dic["melt_vol"] = modes[melt_model]
    melts_dic["T(K)"] = calc["T(K)"]
    melts_dic["P(bar)"] = calc["P(bar)"]
    renormalize!(melts_dic,elements,total=100)

    solid_dic = Dict{String, Array{Float64}}()
    solid_dic["solid_vol"] = tot_vol .- melts_dic["melt_vol"]
    for e in elements
        solid_dic[e] = (calc[e] .- (melts_dic[e] .* melts_dic["melt_vol"]./100)) ./ (solid_dic["solid_vol"]./100)
    end
    solid_dic["T(K)"] = calc["T(K)"]
    solid_dic["P(bar)"] = calc["P(bar)"]
    renormalize!(solid_dic,elements,total=100)

    melts_df = DataFrame(melts_dic)
    solid_df = DataFrame(solid_dic)
    sys_df = DataFrame(sys_comp)
    hbl_df = DataFrame(hbl_comp)

    # Remove data that bulk vol is not greater than 99 or not less than 101
    good_vol = 99 .<= tot_vol .<= 101
    melts_df = melts_df[good_vol, :]
    solid_df = solid_df[good_vol, :]
    sys_df = sys_df[good_vol, :]
    hbl_df = hbl_df[good_vol, :]

    # # Select melts vol >= 30 and <= 70
    # extra_vol = 10 .<= melts_df[!, "melt_vol"] .<= 90
    # melts_df = melts_df[extra_vol, :]
    # solid_df = solid_df[extra_vol, :]

    ## --- 
    melt_cols = vcat(["T(K)", "P(bar)", "melt_vol"], elements)
    melt_sel = melts_df[!, melt_cols]
    sol_cols = vcat(["solid_vol",], elements)
    sol_sel = solid_df[!, sol_cols]

    col_re = ["T(K)", "P(bar)", "melt_vol", "melt_SiO2", "melt_TiO2", "melt_Al2O3", "melt_FeO", "melt_MgO", "melt_CaO", "melt_Na2O", "melt_K2O", "melt_H2O", "melt_CO2", "solid_vol", "solid_SiO2", "solid_TiO2", "solid_Al2O3", "solid_FeO", "solid_MgO", "solid_CaO", "solid_Na2O", "solid_K2O", "solid_H2O", "solid_CO2"]
    bulk_ele = "bulk_" .* names(sys_df)
    hbl_ele = "hbl_" .* names(hbl_df)
    col_new = vcat(col_re, hbl_ele, bulk_ele)


    df = hcat(melt_sel, sol_sel, hbl_df, sys_df,makeunique=true)
    rename!(df, col_new)
    insertcols!(df, 1, :Row => 1:nrow(df))

    CSV.write(joinpath(resultdir, string(sample) * ".csv"), df)
end

# Main function handling parallelization
function mpiperplex()
    # Get world size (number of MPI processes) and world rank (# of this process)
    comm = MPI.COMM_WORLD
    world_size = MPI.Comm_size(comm)
    world_rank = MPI.Comm_rank(comm)
    nworkers = world_size - 1

    cd(@__DIR__)

    # data = importdataset("exp_cpx_hbl.csv",',', importas=:Dict)
    data = importdataset("arc_vol_georoc.csv",',', importas=:Dict)
    # Major elements to include
    elements = ["SiO2", "TiO2", "Al2O3", "FeO", "MgO", "CaO", "Na2O", "K2O", "H2O", "CO2"] # No Cr

    allmajors = hcat((data[x] for x in elements)...)
    # # Renormalize to 100%
    # majors_tot = sum(eachcol(allmajors))
    # allmajors = allmajors ./ majors_tot .* 100    # hasallmajors = vec(sum(isnan.(allmajors), dims=2) .== 0)

    # # Input parameters
    # indexes = findall(hasallmajors)
    dpdz = 2818 * 9.8 / 10^5 *10^3
    maxdepth = 70 # km
    P_range = (dpdz, maxdepth*dpdz) # Pressure range to explore, bar (roughly 1-100 km depth)
    T_surf = 873.15
    # T_range = (673.15, 1273.15) # Temperature range to explore, K
    melts_solids = DataFrame()

    ntasks = length(data["SiO2"])

    hasdata = trues(ntasks)
    same_bulk_idx = Vector{Vector{Int64}}()
    for m in elements
        hasdata .&= .!isnan.(data[m])
    end
    hashes = hash.(eachrow(allmajors)) 
    for i in findall(hasdata)
        push!(same_bulk_idx, findall(hashes .== hashes[i]))
        hasdata[i] &= hashes[i] ∉ view(hashes,1:i-1) 
    end

    same_bulk_idx = collect(Set(same_bulk_idx))

    if (world_rank == 0)
        # Some variables used only on the root node
        wrs = fill(0, nworkers)
        requests = fill(MPI.REQUEST_NULL, nworkers)

        # Listen for task requests from the worker nodes
        for i = 1:nworkers
            requests[i] = MPI.Irecv!(view(wrs,i), comm, source=i, tag=0)
        end
        nr = 0
        # Once any worker asks for a new task, send next task to that worker and keep listening
        for sample in findall(hasdata)
            nr = MPI.Waitany(requests) #+ 1 # Returns rank of next ready task
            MPI.Send(sample, comm, dest=nr, tag=1) # Tell them to work on `sample`
            requests[nr] = MPI.Irecv!(view(wrs,nr), comm, source=nr, tag=0)
        end

        # Wait for all workers to complete, then send the stop signal
        MPI.Waitall(requests)
        sample = -1
        for i = 1:nworkers
            requests[i] = MPI.Isend(sample, comm, dest=i, tag=1)
        end

        # MPI.Waitall(requests)

    else
        # Ask root node for new task
        request = MPI.Isend(world_rank, comm, dest=0, tag=0)
        # Get our new assignment
        sample = MPI.Recv(Int64, comm, source=0, tag=1)

        while sample > 0
            @info "Running sample $sample on rank $world_rank"
            ############################# Do work ##############################
            try
                melt_solid = perplex_geothermo(data, allmajors, elements, P_range, T_surf, sample, same_bulk_idx)
                # melts_solids = vcat(melts_solids, melt_solid)
                # @info "Calculation completed: Index $sample"
            catch
                @warn "PerpleX calculations failed for sample $sample on rank $world_rank"
            end
            ######################## Done doing work ###########################

            # Ask root node for new task
            request = MPI.Isend(world_rank, comm, dest=0, tag=0)
            # Get our new assignment
            sample = MPI.Recv(Int64, comm, source=0, tag=1)
        end
        CSV.write(joinpath(resultdir, "hbl_bulk_georoc_vol_" * string(world_rank) * ".csv"), melts_solids)
    end
end

mpiperplex()

MPI.Finalize()
