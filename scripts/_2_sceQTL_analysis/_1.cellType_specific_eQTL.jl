println("loading Julia packages..")
using Pkg
using HDF5
# using DataFramesMeta
using Logging
using VCFTools
using Suppressor
using DataFrames
using CSV
using MixedModels
using Distributions
using Statistics
using CategoricalArrays
using ProgressMeter
using IJulia
using Base.Threads
using ArgParse
using Debugger

# Function to parse command-line arguments
function parse_arguments()
    s = ArgParseSettings()

    @add_arg_table s begin
        "file_path";
        help = "Anndata object including var and obs for analysis";
        arg_type = String;
        required = true;

        "faFai_path";
        help = "Fai index file of fasta file used for alignment of single cell data";
        arg_type = String;
        required = true;
        
        "vcf_path";
        help = "Genotype VCF file";
        arg_type = String;
        required = true;
        
		"formula_str";
        help = "fomular for association";
        arg_type = String;
        required = true;
        
		"colName_celltype";
        help = "name of column in Obs for celltype";
        arg_type = String;
        required = true;
        default = "final_celltype";
        
        "num_threads";
        help = "Number of threads used in mixed model running";
        arg_type = Int;  # Changed to Int for the number of threads
        required = false;
        default = 20;
        
        "output_file_path";
        help = "The name of the output file";
        arg_type = String;
        required = false;
        default = "sceQTL.result.tsv";
        
        "pad";
        help = "Padding area upstream from TSS and downstream from TES";
        arg_type = Float64;
        required = false;
        default = 1000000.0;
        
    end

    args = parse_args(s)
    return args
end

# Main function
function main()
    args = parse_arguments()
	
    # Access and use the arguments
    file_path = args["file_path"]
    faFai_path = args["faFai_path"]
    vcf_path = args["vcf_path"]
    num_threads = args["num_threads"]
    output_file_path = args["output_file_path"]
    pad = args["pad"]
    colName_celltype = args["colName_celltype"]
	# Choose family
    
    formula_str = args["formula_str"]
    model_formula = eval(Meta.parse(formula_str))
    # formula = eval(Meta.parse(formula_str))
    # Set environment variable
    ENV["JULIA_NUM_THREADS"] = string(num_threads)

    # Print the arguments to confirm they are correctly parsed
    println("file_path: ", file_path)
    println("faFai_path: ", faFai_path)
    println("vcf_path: ", vcf_path)
    println("num_threads: ", num_threads)
    println("output_file_path: ", output_file_path)
    println("pad: ", pad)
    println("colName_celltype: ", colName_celltype)
    println("formula_str: ", formula_str)

# Functions
ENV["JULIA_NUM_THREADS"] = Threads.nthreads()
println("num of threads to be used : "  , Threads.nthreads())

function fit_multiple_models(snplist, geno_gene, df_obs, formula_str, geneName)
    result_forGene_list = Vector{Any}(undef, size(snplist)[1])
    # println("start fit_multiple_models")
    @threads for snpIndex in 1:size(snplist)[1]
        snpName = snplist[snpIndex,"snp"][1]
		sampleVCF = DataFrame(snp = geno_gene[:,snpIndex], sampleid = sampleName)
		sampleVCF.snp = convert(Vector{Int64}, sampleVCF.snp)
		sampleVCF[!, "snp"] = float.(sampleVCF[!, "snp"])
		assocDf = innerjoin(df_obs,sampleVCF, on= :sampleid)
        result_forGene_list[snpIndex] = suppress_specific_warning("Fixed-effects matrix is rank deficient", () -> doAssoc(assocDf, snpName, formula_str, geneName))
    end
    # return result_forGene
    return result_forGene = reduce(vcat, result_forGene_list)
end

function doAssoc(assocDf, snpName, formula_str, geneName)
	result_for_snp_list = Vector{Any}(undef, length(celltypeList))
	#println("start assoc")
    @threads for celltypeIndex in 1:length(celltypeList)
	    celltype = celltypeList[celltypeIndex]
	    #println(celltype)
	    
		assocDf_celltype = assocDf[assocDf.assocCelltype .== celltype,:]
		
		# println(describe(assocDf_celltype))
		try
			m1 = fit(MixedModel, eval(Meta.parse(formula_str)), assocDf_celltype, Poisson(), fast = true)
			fixed_effects = coef(m1)[2]
			std_errors = stderror(m1)[2]
			z_scores = fixed_effects ./ std_errors
			p_values = 2 .* cdf.(Normal(0, 1), -abs.(z_scores))
			
			genotypeMat = makeGenotype(assocDf_celltype)
			alt0 = genotypeMat[genotypeMat.snp .== 0,"count"][1]
			alt1 = genotypeMat[genotypeMat.snp .== 1,"count"][1]
			alt2 = genotypeMat[genotypeMat.snp .== 2,"count"][1]
			
			#println("makeGenotype done")
			
			result_for_snp_list[celltypeIndex] = DataFrame(
			snp = snpName,
			gene = geneName,
			coef = fixed_effects,
			se = std_errors,
			zscroe =z_scores,
			pvale =p_values,
			alt0 = alt0,
			alt1 = alt1,
			alt2 = alt2,
			af =  (alt1+ (2*alt2)) / ((alt0 + alt1 + alt2)*2),
			cellType = celltype
			)
		catch e
			result_for_snp_list[celltypeIndex] = DataFrame(
			snp = snpName,
			gene = geneName,
			coef = "NaN",
			se = "NaN",
			zscroe ="NaN",
			pvale ="NaN",
			alt0 = "NaN",
			alt1 = "NaN",
			alt2 = "NaN",
			af =  "NaN",
			cellType = celltype)
		end
		#println(result_for_snp_list[celltypeIndex])
	end
	# return result_for_snp
	# return result_for_snp = reduce(vcat, result_for_snp_list)
	return result_for_snp = reduce(vcat, result_for_snp_list)
	#CSV.write(csv_file_path, result_for_snp, append=true, header=false, delim='\t')
end

function makeGenotype(assocDf)
	#println("start makeGenotype")
	genotypeMat = combine(groupby(unique(assocDf[:,["sampleid","snp"]]), :snp), nrow => :count)
	for snp_value in 0:2
		exists = any(row -> row.snp == snp_value, eachrow(genotypeMat))
		# Append a new row if the snp value is not present
		if !exists
		    new_row = DataFrame(snp=[snp_value], count=[0])
		    genotypeMat = vcat(genotypeMat, new_row)
		end
	end
	return genotypeMat
end

function read_dataset(dataset)
    if isa(dataset, HDF5.Group) && "categories" in keys(dataset) && "codes" in keys(dataset)
        # Extract categories and codes for categorical data
        categories = dataset["categories"][:]
        codes = dataset["codes"][:]
        
        # Ensure codes are within bounds
        max_code = maximum(codes)
        if max_code >= length(categories)
            error("Some codes are out of bounds for the categories array")
        end
        
        return CategoricalArray(categories[codes .+ 1])  # HDF5 codes are 0-based, so we adjust to 1-based indexing
    else
        # Handle regular dataset
        return dataset[:]
    end
end

function suppress_specific_warning(warning_text::String, func)
    # Capture stderr output to suppress specific warnings
    open("/dev/null", "w") do null
        redirect_stderr(null)  # Suppress stderr output (which includes warnings)
        
        result = func()
        
        redirect_stderr(stdout)  # Restore stderr output
        
        return result
    end
end

# Read chromosome length
faFai = CSV.read(faFai_path, DataFrame; header=false)
faFai = faFai[[occursin(r"chr\d+", val) for val in faFai.Column1],:]
faFai.Column1 .= replace.(faFai.Column1, "chr" => "")
if eltype(faFai.Column1) != Int
    # Convert to Int
    faFai.Column1 = parse.(Int, faFai.Column1)
end

# Open the HDF5 file
adata = h5open(file_path, "r")
obs_group = adata["obs"]

# Initialize dictionary to hold data
data = Dict{String, Any}()

# Extract data from the 'obs' group
for name in keys(obs_group)
    dataset = obs_group[name]
    data[name] = read_dataset(dataset)
end

# Create DataFrame
df_obs = DataFrame(data)

# Optionally, set index if '_index' is present
if "_index" in keys(obs_group)
    df_obs.index = obs_group["_index"][:]
end

df_obs[!,"assocCelltype"] = df_obs[!,colName_celltype]
celltypeList = unique(df_obs[!,"assocCelltype"] )

# Gene name ---------------------------------------------------------------------------------
df_var = adata["var"]
df_var = DataFrame(Dict(k => read(adata["var"][k]) for k in keys(adata["var"])))
if eltype(df_var.chromosome) != Int
    # Convert to Int
    df_var.chromosome = parse.(Int, df_var.chromosome)
end
# df_var.end = parse.(Int, df_var.end)
# df_var.start = parse.(Int, df_var.start)

# size of dimension -----------------------------------------------------------------------
geneNum = size(df_var)[1]
cellNum = size(df_obs)[1]

# vcf ---------------------------------------------------------------------------
println("Loading variants from VCF file...")
@time A = convert_gt(Int8, vcf_path; model = :additive, impute = false, save_snp_info = true, center = false, scale = false)
geno = DataFrame(transpose(A[1]), :auto)
println("A total number of variants :",size(geno)[1] )

vcf_snpInfo = DataFrame(
    chr = A[3],
    pos = A[4],
    snp = A[5])
sampleName = A[2]
vcf_snpInfo.chr = parse.(Int, vcf_snpInfo.chr)
if eltype(vcf_snpInfo.chr) != Int
    # Convert to Int
    vcf_snpInfo.chr = parse.(vcf_snpInfo, df_var.chr)
end



# get genes ----------------- For loop ----------------- 
A = nothing
data = nothing

result_forAllGene=DataFrame()

@showprogress for gene_index in 1:size(df_var)[1]
	result_forGene=DataFrame()
	geneName = df_var[gene_index,"features"]
	geneChr = df_var[gene_index,"chromosome"]
	geneStart = max(df_var[gene_index,"start"] - pad,1)
	geneEnd = min(df_var[gene_index,"end"] + pad, faFai[faFai.Column1 .== geneChr,"Column2"][1])
	
	cellStart=(gene_index-1)*cellNum+1
	cellEnd=(gene_index-1)*cellNum+ cellNum
	df_obs[!,"exp"]= Int.(adata["X"]["data"][cellStart:cellEnd])
	
	# extract vcf
	indices = findall((vcf_snpInfo.chr .== geneChr) .& (vcf_snpInfo.pos .> geneStart) .& (vcf_snpInfo.pos .< geneEnd))
	snplist= vcf_snpInfo[indices,:]
	geno_gene = permutedims(Matrix(geno[indices, :]))
	
	result_forGene = fit_multiple_models(snplist, geno_gene, df_obs, formula_str, geneName)
	
	#result_forAllGene = vcat(result_forAllGene,result_forGene)
	
	CSV.write(output_file_path, result_forGene, append=true, header=false, delim='\t')
end

end

main()
