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
using DataFrames

function find_orfs(fasta_dict::Dict{String, String})
    orfs_df = DataFrame(Sequence_Name = String[], Start_Position = Int[], End_Position = Int[], ORF_Length = Int[], ORF_Sequence = String[])

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


# find open reading grames & return as a list
function find_orfs(fasta_dict::Dict{String, String})
    orfs_dict = Dict{String, Vector{String}}()

    for (name, sequence) in fasta_dict
        orfs = find_long_orfs(sequence)
        if !isempty(orfs)
            orfs_dict[name] = orfs
        end
    end

    return orfs_dict
end

function find_long_orfs(sequence::String)
    orfs = Vector{String}()
    for frame in 1:3
        for i in 1:3:length(sequence) - 2
            codon = sequence[i:i+2]
            if codon == "ATG"
                orf = codon
                j = i + 3
                while j <= length(sequence) - 2
                    current_codon = sequence[j:j+2]
                    if current_codon in ["TAA", "TAG", "TGA"]
                        if length(orf) >= 90
                            push!(orfs, orf)
                        end
                        break
                    else
                        orf *= current_codon
                        j += 3
                    end
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

#first go at an orf-matching function
using DataFrames
using BioSequences
using BioAlignments

function match_orfs(orfs::Dict{String, Vector{String}}, reference_orfs::Dict{String, String})
    result_df = DataFrame(Reference_orf_name = String[], Sample_name = String[], Sequence = String[])

    for (ref_name, ref_sequence) in reference_orfs
        for (seq_name, orfs_list) in orfs
            best_score = 0
            best_match = ""
            best_seq = ""
            for orf in orfs_list
                global_alignment = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), orf, ref_sequence)
                score = BioAlignments.score(global_alignment)
                if score >= best_score
                    best_score = score
                    best_match = seq_name
                    best_seq = orf
                end
            end
            if best_score >= 0.8 # Adjust thresholds as needed
                push!(result_df, (ref_name, best_match, best_seq))
            end
        end
    end

    return result_df
end


#test
kemp_fasta = fasta_readin("data/Kemp/Sequences.fa")
kemp_orfs = find_orfs(kemp_fasta)
ref_orfs = fasta_readin("data/reference/coding_sequences.fasta")
#kemp_matched_orfs = match_orfs(kemp_orfs, ref_orfs)

# chop and change column names to be more manageable
function split_orfs_df(df::DataFrame)
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
    [:Start, :End]))
    select!(df, Not(:Position))
    select!(df, [:Accession, :Sample_name, :ORF_name, :Start, :End, :Sequence])
    return df
end

#kemp_clean_orfs = split_orfs_df(kemp_matched_orfs)