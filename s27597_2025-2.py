#!/usr/bin/env python3
from Bio import Entrez, SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import os

class NCBIRetriever:
    def __init__(self, email, api_key):
        Entrez.email = email
        Entrez.api_key = api_key
        Entrez.tool = "BioScriptEx10"
    
    def search(self, taxid):
        try:
            handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
            records = Entrez.read(handle)
            self.organism = records[0]["ScientificName"]
            print(f"Organism: {self.organism}")
            search_term = f"txid{taxid}[Organism]"
            handle = Entrez.esearch(db="nucleotide", term=search_term, usehistory="y", retmax=0)
            res = Entrez.read(handle)
            self.count = int(res["Count"])
            self.webenv = res["WebEnv"]
            self.query_key = res["QueryKey"]
            return self.count
        except Exception as e:
            print(f"Search error: {e}")
            return 0

    def fetch_records(self, max_records=1000):
        batch_size = 500
        all_records = []
        for start in range(0, max_records, batch_size):
            handle = Entrez.efetch(
                db="nucleotide", rettype="gb", retmode="text",
                retstart=start, retmax=min(batch_size, max_records - start),
                webenv=self.webenv, query_key=self.query_key)
            all_records += list(SeqIO.parse(handle, "gb"))
        return all_records

def filter_records(records, min_len, max_len):
    return [r for r in records if min_len <= len(r.seq) <= max_len]

def save_csv(records, filename):
    data = [{
        "accession": r.id,
        "length": len(r.seq),
        "description": r.description
    } for r in records]
    df = pd.DataFrame(data)
    df.sort_values(by="length", ascending=False, inplace=True)
    df.to_csv(filename, index=False)
    return df

def plot_lengths(df, png_filename):
    df_sorted = df.sort_values("length", ascending=False)
    plt.figure(figsize=(12, 6))
    plt.plot(df_sorted["accession"], df_sorted["length"], marker='o')
    plt.xticks(rotation=90, fontsize=6)
    plt.xlabel("Accession")
    plt.ylabel("Sequence Length")
    plt.title("Sequence Lengths (sorted)")
    plt.tight_layout()
    plt.savefig(png_filename)
    plt.close()

def main():
    email = input("Enter NCBI email: ")
    api_key = input("Enter NCBI API key: ")
    taxid = input("Enter TaxID: ")
    min_len = int(input("Min sequence length: "))
    max_len = int(input("Max sequence length: "))

    retriever = NCBIRetriever(email, api_key)
    count = retriever.search(taxid)
    if count == 0:
        print("No results.")
        return

    print(f"Fetching up to 1000 records out of {count} available...")
    raw_records = retriever.fetch_records(max_records=1000)
    filtered = filter_records(raw_records, min_len, max_len)
    
    if not filtered:
        print("No records match length criteria.")
        return

    csv_name = f"taxid_{taxid}_filtered.csv"
    png_name = f"taxid_{taxid}_plot.png"
    df = save_csv(filtered, csv_name)
    plot_lengths(df, png_name)

    print(f"Saved {len(filtered)} records to {csv_name}")
    print(f"Plot saved as {png_name}")

if __name__ == "__main__":
    main()
