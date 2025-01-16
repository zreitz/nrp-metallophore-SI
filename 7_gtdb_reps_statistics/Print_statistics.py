import pandas as pd
import csv

# Load the reformatted CSV file
file_path = 'All_genomes_AS7_parsed_products_with_taxonomy.csv'  # Replace with your file path
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

# Calculate statistics
total_bacterial_genomes = len(bacterial_data)
bacterial_with_bgc = bacterial_data['Parsed_Products'].apply(lambda x: len(x) > 0).sum()
bacterial_with_nrps = bacterial_data['Parsed_Products'].apply(lambda x: any("NRPS" in product for product in x)).sum()
bacterial_with_nrp_metallophore = bacterial_data['Parsed_Products'].apply(lambda x: any("NRP-metallophore" in product for product in x)).sum()

# Calculate statistics by phylum
phylum_stats = bacterial_data.groupby('Phylum').agg(
    total_genomes=('Phylum', 'size'),
    with_bgc=('Parsed_Products', lambda x: (x.apply(lambda y: len(y) > 0)).sum()),
    with_nrps=('Parsed_Products', lambda x: (x.apply(lambda y: any("NRPS" in product for product in y))).sum()),
    with_nrp_metallophore=('Parsed_Products', lambda x: (x.apply(lambda y: any("NRP-metallophore" in product for product in y))).sum()),
    total_nrp_metallophore=('Parsed_Products', lambda x: sum(product.count("NRP-metallophore") for products in x for product in products))
).reset_index()

# Sort phyla by total number of genomes in descending order
phylum_stats = phylum_stats.sort_values(by='total_genomes', ascending=False)

# Print results
print(f"Total bacterial genomes analyzed: {total_bacterial_genomes}")
print(f"Bacterial genomes with at least one BGC: {bacterial_with_bgc}")
print(f"Bacterial genomes with at least one NRPS: {bacterial_with_nrps}")
print(f"Bacterial genomes with at least one NRP-metallophore: {bacterial_with_nrp_metallophore}")

print("\nPhylum statistics sorted by abundance:")
print("Phylum,total_genomes, with_bgc, with_nrps, with_nrp_metallophore, total_nrp_metallophore")
for _, row in phylum_stats.iterrows():
    print(f"{row['Phylum']}, {row['total_genomes']}, {row['with_bgc']}, {row['with_nrps']}, {row['with_nrp_metallophore']}, {row['total_nrp_metallophore']}")

# Calculate the number of genomes with at least one NRP-metallophore for each phylum
nrp_metallophore_counts = bacterial_data[bacterial_data['Parsed_Products'].apply(lambda x: any("NRP-metallophore" in product for product in x))]['Phylum'].value_counts()

# Calculate the total number of bacterial genomes with at least one NRP-metallophore
total_nrp_metallophore = nrp_metallophore_counts.sum()

# Select the six most abundant phyla in NRP-metallophore
top_six_phyla = nrp_metallophore_counts.head(6)

# Calculate the percentage of all NRP-metallophores in each of these phyla
top_six_phyla_percentage = (top_six_phyla / total_nrp_metallophore) * 100

# Calculate the percentage of genomes with at least one NRP-metallophore in each of these phyla
top_six_phyla_genome_percentage = (top_six_phyla / phylum_stats.set_index('Phylum').loc[top_six_phyla.index, 'total_genomes']) * 100

# Print the results
print("\nSix most abundant phyla in NRP-metallophore and their percentages:")
for phylum, count in top_six_phyla.items():
    percentage = top_six_phyla_percentage[phylum]
    genome_percentage = top_six_phyla_genome_percentage[phylum]
    print(f"{phylum}: {count} genomes, {percentage:.2f}% of all NRP-metallophores, {genome_percentage:.2f}% of genomes in this phylum have at least one NRP-metallophore")
