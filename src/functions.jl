# somewhere to put functions as they're worked on

# convert Times.in to Dates.in 
using DataFrames
using Dates
using DelimitedFiles
using CSV

function convert_to_date(folder::String)
    times_file = joinpath(folder, "Times.in")
    dates_file = joinpath(folder, "Dates.in")
    if !isfile(times_file)
        error("File 'Times.in' not found in the folder.")
    end

    times = readdlm(times_file, '\t', Int)
    dates = Date[]
    for time in times
        push!(dates, Date(1900, 1, 1) + Day(time))
    end

    df = DataFrame(Dates = dates)
    CSV.write(dates_file, df)

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

# read in fasta data
using DataFrames

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


# translate dictionary ORFs
function translate_orfs(fasta_dict::Dict{String, Vector{String}})
    translated_orfs = Dict{String, Vector{String}}()

    codon_table = Dict(
        "TTT" => "F", "TTC" => "F", "TTA" => "L", "TTG" => "L",
        "TCT" => "S", "TCC" => "S", "TCA" => "S", "TCG" => "S",
        "TAT" => "Y", "TAC" => "Y", "TAA" => "*", "TAG" => "*",
        "TGT" => "C", "TGC" => "C", "TGA" => "*", "TGG" => "W",
        "CTT" => "L", "CTC" => "L", "CTA" => "L", "CTG" => "L",
        "CCT" => "P", "CCC" => "P", "CCA" => "P", "CCG" => "P",
        "CAT" => "H", "CAC" => "H", "CAA" => "Q", "CAG" => "Q",
        "CGT" => "R", "CGC" => "R", "CGA" => "R", "CGG" => "R",
        "ATT" => "I", "ATC" => "I", "ATA" => "I", "ATG" => "M",
        "ACT" => "T", "ACC" => "T", "ACA" => "T", "ACG" => "T",
        "AAT" => "N", "AAC" => "N", "AAA" => "K", "AAG" => "K",
        "AGT" => "S", "AGC" => "S", "AGA" => "R", "AGG" => "R",
        "GTT" => "V", "GTC" => "V", "GTA" => "V", "GTG" => "V",
        "GCT" => "A", "GCC" => "A", "GCA" => "A", "GCG" => "A",
        "GAT" => "D", "GAC" => "D", "GAA" => "E", "GAG" => "E",
        "GGT" => "G", "GGC" => "G", "GGA" => "G", "GGG" => "G"
    )

    for (header, sequences) in fasta_dict
        translated_sequences = Vector{String}()
        for seq in sequences
            codons = [seq[i:i+2] for i in 1:3:length(seq)-2]
            translated_sequence = join([get(codon_table, codon, "X") for codon in codons], "")
            push!(translated_sequences, translated_sequence)
        end
        translated_orfs[header] = translated_sequences
    end

    return translated_orfs
end

# read in reference FASTAs function
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
    DataFrame(reduce(vcat, permutedims.(split.(df.Position, '-'))),
    [:Start_Position, :End_Position]))
    select!(df, Not(:Position))
    select!(df, [:Accession, :ORF_name, :Start_Position, :End_Position, :Sequence])
    df.Start_Position = parse.(Int64, df.Start_Position)
    df.End_Position = parse.(Int64, df.End_Position)
    return df
end


#first go at an orf-matching function
function match_orfs(orf_df::DataFrame, ref_orf_df::DataFrame)
    result_df = DataFrame(ORF_name = String[], Sequence_Name = String[], Start_Position = Int[], End_Position = Int[], Matched_Sequence = String[])

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
                push!(result_df, (ref_orf_df[i, :ORF_name], seq_name, orf_df_subset[best_match_idx, :Start_Position], orf_df_subset[best_match_idx, :End_Position], orf_df_subset[best_match_idx, :Sequence]))
            end
        end
    end

    return result_df
end


#test
kemp_fasta = fasta_readin("data/Kemp/Sequences.fa")
kemp_orfs = find_orfs(kemp_fasta)
ref_orfs = ref_readin("data/reference/coding_sequences.fasta")
kemp_matched_orfs = match_orfs(kemp_orfs, ref_orfs)

