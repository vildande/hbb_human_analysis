from Bio import SeqIO


def find_variation_details(record):
    for feature in record.features:
        if feature.type == "variation":
            position = int(feature.location.start)
            original_nucleotide = feature.qualifiers.get('replace', [None])[0]
            return position, original_nucleotide if original_nucleotide == None else original_nucleotide.upper()
    return None, None


def analyze_mutation(sequence, position, original_nucleotide):
    if position is None:
        print("No variation position found in the record.")
        return

    actual_nucleotide = sequence[position]
    print(f"Nucleotide at position {position + 1}: {actual_nucleotide}")

    if actual_nucleotide == original_nucleotide:
        print(f"No change: Nucleotide is not mutated")
    else:
        print(f"Mutation detected: Found '{actual_nucleotide}' instead of expected '{original_nucleotide}'.")


def main():
    filename = "hbb_human.gb"
    with open(filename, "r") as file:
        record = SeqIO.read(file, "genbank")
        
    # Homozygous-29 A->G mutation detected at the promoter region of beta globin gene on a black patient with beta-thalassemia major
    position, original_nucleotide = find_variation_details(record)
    analyze_mutation(record.seq, position, original_nucleotide)


if __name__ == "__main__":
    main()