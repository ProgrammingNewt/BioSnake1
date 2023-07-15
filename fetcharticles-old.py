from Bio import Entrez

# Always tell NCBI who you are
Entrez.email = "your.email@example.com"

# Use the 'esearch' function to search the 'pubmed' database for articles
# that match the term "disease environmental factors" and have "free full text" available.
handle = Entrez.esearch(db="pubmed", term="disease environmental factors free full text[sb]")

# Read the resulting record and print the IDs of the matching articles
record = Entrez.read(handle)

print("PubMed IDs of matching articles with 'free full text':")
print(record["IdList"])

handle.close()
