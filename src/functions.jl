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

# find open reading grames & return as a list

function find_orfs(fasta_dict::Dict{String, String})
    all_orfs = Dict{String, Vector{String}}()

    function find_orfs_in_seq(header, seq)
        orfs = Vector{String}()
        for frame in [seq, reverse(seq)]
            for i in 1:3:length(frame) - 2
                codon = frame[i:i + 2]
                if codon == "ATG"
                    current_orf = "ATG"
                    j = i + 3
                    while j <= length(frame) - 2
                        current_codon = frame[j:j + 2]
                        if current_codon in ["TAA", "TAG", "TGA"]
                            push!(orfs, current_orf * current_codon)
                            break
                        else
                            current_orf *= current_codon
                            j += 3
                        end
                    end
                end
            end
        end
        return orfs
    end

    for (header, sequence) in fasta_dict
        orfs = find_orfs_in_seq(header, sequence)
        all_orfs[header] = orfs
    end

    return all_orfs
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


# example with Kemp data
kemp_fasta = fasta_readin("data/Kemp/Sequences.fa")
kemp_orfs = find_orfs(kemp_fasta)
kemp_aa_orfs = translate_orfs(kemp_orfs)

# example with CAMP000427 data
CAMP000427_fasta = fasta_readin("data/CAMP000427/Sequences.fa")
CAMP000427_orfs = find_orfs(CAMP000427_fasta)
CAMP000427_aa_orfs = translate_orfs(CAMP000427_orfs)

#first go at an orf-matching function
using DataFrames
using BioSequences
using BioAlignments

function match_orfs(orfs::Dict{String, Vector{String}}, reference_orfs::Dict{String, String})
    result_df = DataFrame(Reference_ORF_Name = String[], Sequence_Name = String[], Sequence = String[])

    for (ref_name, ref_sequence) in reference_orfs
        for (orfs_name, orfs_list) in orfs
            for orf in orfs_list
                global_alignment = BioAlignments.pairalign(BioAlignments.LevenshteinDistance(), orf, ref_sequence)
                score = BioAlignments.score(global_alignment)
                if score >= 0.8 * length(ref_sequence) # Adjust the threshold as needed
                    push!(result_df, (ref_name, orfs_name, orf))
                end
            end
        end
    end

    return result_df
end

#test
ref_orfs = fasta_readin("data/reference/coding_sequences.fasta")
match_orfs(kemp_orfs, ref_orfs)