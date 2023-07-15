from Bio import Entrez
from xml.etree import ElementTree as ET
from xml.dom.minidom import parseString
import requests
from bs4 import BeautifulSoup
from selenium import webdriver

def pretty_print_xml(xml_string):
    xml = parseString(xml_string)
    pretty_xml = xml.toprettyxml()
    print(pretty_xml)
    
def search(query):
    Entrez.email = 'codingnewt@gmail.com'
    handle = Entrez.esearch(db='pubmed', 
                            sort='relevance', 
                            retmax='1',  # adjust this to retrieve more or fewer articles
                            retmode='xml', 
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list):
    ids = ','.join(id_list)
    handle = Entrez.efetch(db='pubmed',
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    return results



def fetch_links(id_list):
    print(id_list)
    link_handle = Entrez.elink(dbfrom="pubmed", id=id_list, linkname="pubmed_pubmed", cmd="llinks")

    
    link_results = ET.fromstring(link_handle.read())  # read the response as XML
    print("link_result", link_results)
    link_handle.close()

    links = []
    for id_url_list in link_results.findall(".//IdUrlList"):
        for id_url_set in id_url_list.findall(".//IdUrlSet"):
            for obj_url in id_url_set.findall(".//ObjUrl"):
                url = obj_url.find(".//Url")
                if url is not None:
                    if "www.ncbi.nlm.nih.gov/pmc" in url.text:
                        links.append(url.text)
                        print("url:", url.text)
    

    return links
    


''' def main():
    results = search('disease environmental factors')
    id_list = results['IdList']
    papers = fetch_details(id_list)
    links = fetch_links(id_list)
    
    for i, paper in enumerate(papers['PubmedArticle']):
        print(f"Title: {paper['MedlineCitation']['Article']['ArticleTitle']}")
        print(f"Abstract: {paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]}")
        if 'free full text' in paper['MedlineCitation']['Article']['ELocationID'][0]:
            print(f"Full text link: {links['LinkSet'][i]['IdUrlList'][0]['IdUrlSet'][0]['ObjUrl'][0]['Url']}")
        print("\n")
'''
def main():
    results = search('disease environmental factors')
    id_list = results['IdList']
    papers = fetch_details(id_list)
    links = fetch_links(id_list)
    print("id_list", id_list)
    #print(results)
    
    for i, paper in enumerate(papers['PubmedArticle']):
        print(f"Title: {paper['MedlineCitation']['Article']['ArticleTitle']}")
        try:
            #print(f"Abstract: {paper['MedlineCitation']['Article']['Abstract']['AbstractText'][0]}")
            if 'free full text' in paper['MedlineCitation']['Article']['ELocationID'][0]:
                print(f"Full text link: {links[i]}")
        except KeyError:
            print("No ELocationID found for this article.")
        print("\n")

from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service



def get_matching_url_selenium(page_url, link_text):
    # Set up the driver (this example uses Chrome; adjust for your browser)
    driver = webdriver.Chrome(service=Service(r'C:\Windows\chromedriver.exe'))
    #driver = webdriver.Chrome()  
    driver.get(page_url)
    
    # Find all links and check their text
    #links = driver.find_elements_by_tag_name("a")
    links = driver.find_elements(By.TAG_NAME, 'a')
    for link in links:
        #print(link.text, "-", link.get_attribute("href"))
        if link_text in link.text:
            url = link.get_attribute("href")
            driver.quit()
            return url

    # Close the driver if no match was found
    driver.quit()
    return None


def download_pdf(url, destination):
    response = requests.get(url)

    with open(destination, 'wb') as output_file:
        output_file.write(response.content)


def get_matching_url(page_url, link_text):
    response = requests.get(page_url)
    soup = BeautifulSoup(response.text, 'html.parser')

    for link in soup.find_all('a'):
        print(link)
        if link.string == link_text:
            return link.get('href')

    return None  # No matching link was found

if __name__ == '__main__':
    main()
    
handle = Entrez.elink(dbfrom="pubmed", db="pmc", linkname="pubmed_pmc", id="['30388784']")
link_results = ET.fromstring(handle.read())  # read the response as XML
handle.close()

for db in link_results.findall(".//LinkSetDb"):
    for link in db.findall(".//Id"):
        pmc_id = [link.text]
        print("PMC id:", pmc_id)


page_url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC" + pmc_id[0] + "/"  # replace with your target URL
print(page_url)

'''
link_text = "PDF"  # replace with the text you're searching for
url = get_matching_url_selenium(page_url, link_text)

if url:
    print(f"Found URL: {url}")
else:
    print("No matching URL found.")

pdf_folder = "pdf_files/"

output_path = pdf_folder + pmc_id[0] + "_article.html"  # replace with where you want to save the PDF

download_pdf(url, output_path)
'''




                        

