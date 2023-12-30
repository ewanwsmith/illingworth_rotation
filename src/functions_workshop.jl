# somewhere to put functions as they're worked on

#load functions
include("/Users/ewansmith/Documents/PhD/Rotation 1 - Illingworth/illingworth_rotation/src/functions.jl")

function true_to_T(df::DataFrame)
    # Replace "true" with "T" in Original_Base column
    df[!, :Original_Base] .= replace.(df[!, :Original_Base], "true" => "T")
    
    # Replace "true" with "T" in Variant_Base column
    df[!, :Variant_Base] .= replace.(df[!, :Variant_Base], "true" => "T")
    
    return df
end


folder = "data/CAMP007136"

df = readin_consensus(folder)
var = readin_variants(folder)
df = pull_proteins(df, positions_df)
df = locate_variants(var, df)
df = true_to_T(df)