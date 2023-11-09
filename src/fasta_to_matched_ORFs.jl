# take SARS-CoV_2 .fasta files, extract all possible open reading frames, and match these to NCBI reference ORFs

#input paths to folders with .fasta files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"] 
#input path to reference ORFs
reference_path = "data/reference/coding_sequences.fasta" 

#load dependencies
using DataFrames
using CSV

# read in reference ORF fasta
function ref_readin(file_path::AbstractString)
    orf_names = String[]
    sequences = String[]
    current_sequence = ""

    open(file_path) do file
        current_orf_name = ""
        for line in eachline(file)
            if occursin(">", line)
                if current_orf_name != ""
                    push!(orf_names, current_orf_name)
                    push!(sequences, current_sequence)
                    current_sequence = ""
                end
                current_orf_name = split(strip(line), '>')[2]
            else
                current_sequence *= strip(line)
            end
        end
        if current_orf_name != ""
            push!(orf_names, current_orf_name)
            push!(sequences, current_sequence)
        end
    end

    df = DataFrame(
        Reference_orf_name = orf_names,
        Sequence = sequences
    )
    df = hcat(df,
    DataFrame(reduce(vcat, permutedims.(split.(df.Reference_orf_name, '|'))),
    [:Gene_accession_and_position, :ORF_name]))
    select!(df, Not(:Reference_orf_name))
    df = hcat(df,
    DataFrame(reduce(vcat, permutedims.(split.(df.Gene_accession_and_position, ':'))),
    [:Accession, :Position]))
    select!(df, Not(:Gene_accession_and_position))
    df = hcat(df,
    DataFrame(reduce(vcat, permutedims.(split.(df.Position, ".."))),
    [:Start_Position, :End_Position]))
    select!(df, Not(:Position))
    select!(df, [:Accession, :ORF_name, :Start_Position, :End_Position, :Sequence])
    df.Start_Position = parse.(Int64, df.Start_Position)
    df.End_Position = parse.(Int64, df.End_Position)
    return df
end

# read in fasta data
function fasta_readin(filename::AbstractString)
    fasta_dict = Dict{String, String}()
    current_header = ""
    current_sequence = ""

    open(filename) do file
        for line in eachline(file)
            if startswith(line, ">")
                if current_header != ""
                    fasta_dict[current_header] = current_sequence
                    current_sequence = ""
                end
                current_header = strip(line[2:end])
            else
                current_sequence *= strip(line)
            end
        end
    end

    if current_header != ""
        fasta_dict[current_header] = current_sequence
    end

    return fasta_dict
end

# find possible ORFs
function find_orfs(fasta_dict::Dict{String, String})
    orfs_df = DataFrame(Sequence_Name = String[], Start_Position = Int[], End_Position = Int[], ORF_Length = Int[], Sequence = String[])

    for (sequence_name, sequence) in fasta_dict
        orfs = find_orfs_in_sequence(sequence, sequence_name)
        for orf in orfs
            push!(orfs_df, orf)
        end
    end

    return orfs_df
end

function find_orfs_in_sequence(sequence::String, sequence_name::String)
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    min_orf_length = 90
    orfs = []

    for i in 1:length(sequence) - 2
        codon = sequence[i:i+2]
        if codon == start_codon
            orf_start = i
            for j in i+3:3:length(sequence)-2
                codon = sequence[j:j+2]
                if codon in stop_codons
                    orf_end = j + 2
                    orf_length = orf_end - orf_start + 1
                    if orf_length >= min_orf_length
                        orf_sequence = sequence[orf_start:orf_end]
                        push!(orfs, (sequence_name, orf_start, orf_end, orf_length, orf_sequence))
                    end
                    break
                end
            end
        end
    end

    return orfs
end

# match ORFs to reference ORFs
function match_orfs(orf_df::DataFrame, ref_orf_df::DataFrame)
    result_df = DataFrame(ORF_name = String[], Sequence_Name = String[], Start_Position = Int[], End_Position = Int[], Reference_Sequence = String[], Matched_Sequence = String[])

    for seq_name in unique(orf_df.Sequence_Name)
        orf_df_subset = filter(row -> row.Sequence_Name == seq_name, orf_df)
        for i in 1:size(ref_orf_df, 1)
            ref_start = ref_orf_df[i, :Start_Position]
            ref_end = ref_orf_df[i, :End_Position]

            best_match_idx = 0
            min_distance = typemax(Int)
            for j in 1:size(orf_df_subset, 1)
                orf_start = orf_df_subset[j, :Start_Position]
                orf_end = orf_df_subset[j, :End_Position]

                if orf_start >= ref_start && orf_end <= ref_end
                    distance = min(abs(orf_start - ref_start), abs(orf_end - ref_end))
                    if distance < min_distance
                        min_distance = distance
                        best_match_idx = j
                    end
                end
            end

            if best_match_idx != 0
                push!(result_df, (ref_orf_df[i, :ORF_name], seq_name, orf_df_subset[best_match_idx, :Start_Position], orf_df_subset[best_match_idx, :End_Position], ref_orf_df[i, :Sequence], orf_df_subset[best_match_idx, :Sequence]))
            end
        end
    end

    return result_df
end

# read in reference ORFs
ref_orfs = ref_readin(reference_path)
println("Reading in reference open reading frames from: $reference_path")

# run above functions
function find_match_orfs(folder::String)
    fasta_file = joinpath(folder, "sequences.fa") #find sample .fasta files in folder
    println("Running fasta_readin() function on folder: $folder")
    df = fasta_readin(fasta_file) #readin fasta file

    println("Running find_orfs() function on folder: $folder")
    orfs_result = find_orfs(df) #find possible ORFs in fasta file
    println("Running match_orfs() function on folder: $folder")
    match_orfs_result = match_orfs(orfs_result, ref_orfs) #find best matches to reference ORFs by length

    csv_path = joinpath(folder, "matched_orfs.csv") 
    println("Writing matched_orfs.csv to folder: $folder")
    CSV.write(csv_path, match_orfs_result) #write results to a CSV in sample folder

    return match_orfs_result
end

# run on folder_list folders
for folder in folder_list
    find_match_orfs(folder)
end