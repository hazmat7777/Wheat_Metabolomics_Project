import pandas as pd

def create_simple_taxonomy_table(blast_file, output_file='ESV_simple_taxonomy.csv'):
    """Create a simple taxonomy table with just ESV and taxonomic levels"""
    
    # Read BLAST results - your actual format
    df = pd.read_csv(blast_file, sep='\t', header=None,
                     names=['qacc', 'sacc', 'stitle', 'pident', 'length',
                            'mismatch', 'gapopen', 'qstart', 'qend',
                            'sstart', 'send', 'evalue', 'bitscore', 'staxid'])
    
    # Add ESV column (same as qacc in your case)
    df['ESV'] = df['qacc']
    
    def parse_silva_taxonomy(title):
        if pd.isna(title):
            return ['Unknown'] * 7
        
        title = str(title)
        
        # SILVA format is: accession Kingdom;Phylum;Class;Order;Family;Genus;Species
        # Split by space and take taxonomy part (everything after first space)
        parts = title.split(' ', 1)
        if len(parts) < 2:
            return ['Unknown'] * 7
        
        taxonomy_string = parts[1]  # This should be "Bacteria;Verrucomicrobiota;..."
        tax_levels = taxonomy_string.split(';')
        tax_levels = [level.strip() if level.strip() else 'Unknown'
                     for level in tax_levels]
        
        while len(tax_levels) < 7:
            tax_levels.append('Unknown')
        
        return tax_levels[:7]
    
    # Get best hit per ESV
    best_hits = df.groupby('ESV').first().reset_index()
    
    # Create simple taxonomy table
    taxonomy_rows = []
    for _, row in best_hits.iterrows():  # Fixed the syntax here
        print(f"Processing ESV: {row['ESV']}")
        print(f"Title: {row['stitle']}")  # Debug: see what titles look like
        tax_levels = parse_silva_taxonomy(row['stitle'])
        print(f"Parsed taxonomy: {tax_levels}")  # Debug: see parsed result
        taxonomy_rows.append([row['ESV']] + tax_levels)
    
    simple_df = pd.DataFrame(taxonomy_rows,
                           columns=['ESV', 'Kingdom', 'Phylum', 'Class',
                                   'Order', 'Family', 'Genus', 'Species'])
    simple_df.set_index('ESV', inplace=True)
    simple_df.to_csv(output_file)
    
    print(f"Simple taxonomy table saved to: {output_file}")
    print(f"Shape: {simple_df.shape}")
    print("\nFirst few rows:")
    print(simple_df.head())
    
    return simple_df

# Usage
if __name__ == "__main__":
    # Create taxonomy table
    simple_table = create_simple_taxonomy_table(
        blast_file='../results/ESV_taxonomy2.txt',
        output_file='../results/ESV_tax_table.csv'  # Changed to .csv extension
    )
    