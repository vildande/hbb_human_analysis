from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
import numpy as np


def main():
    filename = "hbb_human.gb"
    with open(filename, "r") as file:
        record = SeqIO.read(file, "genbank")

    print(f"ID: {record.id}")
    print(f"Name: {record.name}")
    print(f"Description: {record.description}")

    feature_types = []
    gc_contents = []


    for feature in record.features:
        if feature.type in ["regulatory", "exon", "CDS"]:
            if feature.type == "exon" and "number" in feature.qualifiers:
                feature_label = f"{feature.type} {feature.qualifiers['number'][0]}"
            else:
                feature_label = feature.type
            print("========================================")

            print(f"Feature Type: {feature_label}")
            print(f"Location: {feature.location}")

            seq = feature.extract(record.seq)
            print(f"\nExtracted Sequence: {seq}")

            gc_content = gc_fraction(seq)
            print(f"\nGC Content: {gc_content}\n")


            feature_types.append(feature_label)
            gc_contents.append(gc_content)

            if feature.type == "CDS":
                print(f"mRNA Sequence: {seq.transcribe()}")
                print(f"\nProtein Sequence: {seq.translate()}")
            

    plt.figure(figsize=(10, 6))
    plt.bar(feature_types, gc_contents, color='springgreen', align='center', width=0.4)
    plt.xlabel('Feature Type')
    plt.xticks(rotation=45)
    plt.ylabel('GC Content')
    plt.ylim(0, 1.0)
    plt.yticks(np.arange(0, 1.1, 0.1))
    plt.title('GC Content by Feature Type in HBB Gene')
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()