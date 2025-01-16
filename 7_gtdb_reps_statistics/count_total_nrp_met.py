import pandas as pd
import csv

# Load the reformatted CSV file
file_path = 'All_genomes_AS7_parsed_products_with_taxonomy.csv'  # Replace with your file path

# Read the CSV file
data = pd.read_csv(file_path, sep=',', quoting=csv.QUOTE_MINIMAL)

# Ensure 'Parsed_Products' column is treated as a list of strings
data['Parsed_Products'] = data['Parsed_Products'].apply(lambda x: x.strip('[]').split(', ') if pd.notna(x) else [])

# Extract phylum from `gtdb_taxonomy`
def extract_phylum(taxonomy):
    if pd.isna(taxonomy):
        return None
    for part in taxonomy.split(';'):
        if part.startswith(' p__'):
            return part.split('__')[1].strip()
    return None

# Add a column for phylum
data['Phylum'] = data['gtdb_taxonomy'].apply(extract_phylum)

# Filter data to include only bacterial genomes
bacterial_data = data[data['gtdb_taxonomy'].str.contains('d__Bacteria', na=False)]

# Calculate the total number of NRP-metallophore regions, excluding NRPS-like and NRPS-independent-siderophore
total_nrp_metallophore_regions = bacterial_data['Parsed_Products'].apply(
    lambda x: sum("NRP-metallophore" in product and "NRPS-like" not in product and "NRPS-independent-siderophore" not in product for product in x)
).sum()

# Calculate the total number of NRPS regions, excluding NRPS-like and NRPS-independent-siderophore
total_nrps_regions = bacterial_data['Parsed_Products'].apply(
    lambda x: sum("NRPS" in product and "NRPS-like" not in product and "NRPS-independent-siderophore" not in product for product in x)
).sum()

# Calculate the percentage of NRPS regions that are NRP-metallophore
percentage_nrp_metallophore_of_nrps = (total_nrp_metallophore_regions / total_nrps_regions) * 100

# Print the results
print(f"Total NRP-metallophore regions: {total_nrp_metallophore_regions}")
print(f"Total NRPS regions: {total_nrps_regions}")
print(f"Percentage of NRPS regions that are NRP-metallophores: {percentage_nrp_metallophore_of_nrps:.2f}%")