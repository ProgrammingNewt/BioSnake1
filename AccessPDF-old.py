
# Import requests and BeautifulSoup libraries
import requests
from bs4 import BeautifulSoup

# Define the base url and the page url
base_url = "https://pubmed.ncbi.nlm.nih.gov/"
page_url = "https://pubmed.ncbi.nlm.nih.gov/?term=disease+environmental+factors&filter=simsearch2.ffrft"

# Define a function to get the pdf link from an article link
def get_pdf_link(article_link):
    # Get the html content of the article link
    article_html = requests.get(article_link).text
    # Parse the html with BeautifulSoup
    article_soup = BeautifulSoup(article_html, "html.parser")
    # Find the link to the PMC article
    pmc_link = article_soup.find("a", {"data-ga-action": "PMC"})["href"]
    # Get the html content of the PMC link
    pmc_html = requests.get(pmc_link).text
    # Parse the html with BeautifulSoup
    pmc_soup = BeautifulSoup(pmc_html, "html.parser")
    # Find the link to the pdf file
    pdf_link = pmc_soup.find("a", {"data-ga-action": "PDF"})["href"]
    # Return the pdf link
    return pdf_link

# Define a function to download a pdf file from a pdf link
def download_pdf(pdf_link):
    # Get the file name from the pdf link
    file_name = pdf_link.split("/")[-1]
    # Get the binary content of the pdf link
    pdf_content = requests.get(pdf_link).content
    # Open a file with the file name and write mode
    with open(file_name, "wb") as f:
        # Write the binary content to the file
        f.write(pdf_content)
        # Print a message that the file is downloaded
        print(f"Downloaded {file_name}")

# Get the html content of the page url
page_html = requests.get(page_url).text
# tests if page_html works 
# print(page_html)
# Parse the html with BeautifulSoup
page_soup = BeautifulSoup(page_html, "html.parser")
print(page_soup)

# Find all the links to the articles on the page
article_links = page_soup.find_all("a", {"class": "docsum-title"})
# print(article_links)
# Loop through each article link
for article_link in article_links:
    # Get the full article link by joining it with the base url
    full_article_link = base_url + article_link["href"]
    # Get the pdf link from the full article link using the function defined above
    pdf_link = get_pdf_link(full_article_link)
    # Download the pdf file from the pdf link using the function defined above
    download_pdf(pdf_link)
