# #input paths to folders with Variant_list.out files
folder_list = ["data/CAMP000427", "data/CAMP001339", "data/CAMP001490",
"data/CAMP001523", "data/CAMP002274", "data/CAMP003468",
"data/CAMP004884", "data/CAMP007136", "data/Kemp"] 

# load dependencies
using DataFrames
using CSV
using FileIO

function variant_readin(folder::AbstractString)
    # Construct the file path for Variant_list.out
    variant_file_path = joinpath(folder, "Variant_list.out")

    # Check if the file exists
    if isfile(variant_file_path)
        # Read the data from Variant_list.out with space delimiter
        data = readdlm(variant_file_path, ' ', header=false)

        # Create a DataFrame
        df = DataFrame(Position = data[:, 1], Original_Base = data[:, 2], Variant_Base = data[:, 3])

        # Save the DataFrame as a CSV file
        csv_file_path = joinpath(folder, "Variant_list.csv")
        CSV.write(csv_file_path, df)

        println("DataFrame saved as CSV file: $csv_file_path")
    else
        println("Error: Variant_list.out not found in the specified folder.")
    end
end

# run on folder_list folders
for folder in folder_list
    variant_readin(folder)
end