# convert Times.in to Dates.in and join dates to matched_orfs.csv

#input paths to folders with Dates.in and matched_orfs.csv files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
                "data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
                "data/CAMP004884", "data/CAMP007136", "data/Kemp"] 

#run fasta_to_matched_ORFs.jl script
# include("src/fasta_to_matched_ORFs.jl")

# load dependencies
using DataFrames
using Dates
using DelimitedFiles
using CSV
using FASTX

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

# match sequence names to date
function sequence_dates(folder_path::String)
    # Construct the file paths
    dates_file_path = joinpath(folder_path, "Dates.in")
    sequences_file_path = joinpath(folder_path, "Sequences.fa")

    # Check if the files exist
    if isfile(dates_file_path) && isfile(sequences_file_path)
        # Read the dates from the file
        dates = []
        open(dates_file_path, "r") do file
            # Skip the first line
            readline(file)
            for line in eachline(file)
                push!(dates, strip(line))
            end
        end

        # Read the sequences and sequence names from the FASTA file
        sequence_names = []
        sequence_dates = []
        open(sequences_file_path, "r") do file
            current_sequence = ""
            current_header = ""
            for line in eachline(file)
                if startswith(line, '>')
                    if current_header != ""
                        push!(sequence_names, current_header)
                        push!(sequence_dates, dates[findfirst(x -> x == current_header, sequence_names)])
                    end
                    current_header = split(strip(line[2:end]), ' ')[1]
                else
                    current_sequence *= strip(line)
                end
            end
            push!(sequence_names, current_header)
            push!(sequence_dates, dates[findfirst(x -> x == current_header, sequence_names)])
        end

        # Create a DataFrame with sequence names, sequences, and corresponding dates
        result_df = DataFrame(
            Sequence_Name = sequence_names,
            Date = sequence_dates
        )

        # Save the DataFrame as a CSV
        output_file_path = joinpath(folder_path, "sequence_dates.csv")
        CSV.write(output_file_path, result_df)

        return result_df
    else
        println("Files not found.")
        return DataFrame()
    end
end

# join dates to matched_orfs.csv and save
function join_dates(folder_path::String)
    # Construct the file paths
    matched_orfs_file_path = joinpath(folder_path, "matched_orfs.csv")
    sequence_dates_file_path = joinpath(folder_path, "sequence_dates.csv")

    # Check if the files exist
    if isfile(matched_orfs_file_path) && isfile(sequence_dates_file_path)
        # Read the data from "matched_orfs.csv" and "sequence_dates.csv"
        matched_orfs_df = CSV.File(matched_orfs_file_path) |> DataFrame
        sequence_dates_df = CSV.File(sequence_dates_file_path) |> DataFrame

        # Join the dataframes by the "Sequence_Name" column
        joined_df = leftjoin(matched_orfs_df, sequence_dates_df, on=:Sequence_Name)

        # Save the resultant joined dataframe as "matched_orfs.csv"
        output_file_path = joinpath(folder_path, "matched_orfs.csv")
        CSV.write(output_file_path, joined_df)

        return joined_df
    else
        println("Files not found.")
        return DataFrame()
    end
end

# Main function
function times_to_date_join(folder)
    println("Running convert_to_date() function on folder: $folder")
    convert_to_date(folder)
    println("Running sequence_dates() function on folder: $folder")
    sequence_dates(folder)
    println("Running join_dates() function on folder: $folder")
    join_dates(folder)
end

# run on folder_list folders
for folder in folder_list
    times_to_date_join(folder)
end