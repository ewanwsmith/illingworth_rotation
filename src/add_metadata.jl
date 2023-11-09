using CSV
using DataFrames
using FilePaths

include("fasta_to_matched_ORFS.jl")

function join_metadata(folder::String)
    # Check if the folder contains the required files
    matched_orfs_path = joinpath(folder, "matched_orfs.csv")
    metadata_path = joinpath(folder, "metadata.csv")

    if !isfile(matched_orfs_path) || !isfile(metadata_path)
        println("Warning: One or both files not found in the folder.")
        return
    end

    # Read the data into dataframes
    matched_orfs_df = CSV.read(matched_orfs_path, DataFrame)
    metadata_df = CSV.read(metadata_path, DataFrame)

    # Perform the join operation
    joined_df = rightjoin(matched_orfs_df, metadata_df, on = :Sequence_Name)

    # Save the joined dataframe in the folder as matched_orfs.csv
    joined_df_path = joinpath(folder, "matched_orfs.csv")
    CSV.write(joined_df_path, joined_df)

    println("ORF dataframe and Metadata dataframe have been joined and saved as matched_orfs.csv in $folder.")
end

for folder in folder_list
    join_metadata(folder)
end