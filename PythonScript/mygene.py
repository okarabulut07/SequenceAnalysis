import time
from bs4 import BeautifulSoup
import requests

with open('de_genes_new.txt') as f:
    lines = f.readlines()
    stripped = [s.strip() for s in lines]
for elem in stripped:
    if "." not in elem:
        #print(elem)
        url = "https://www.ncbi.nlm.nih.gov/gene/?term=" + elem
        # print(url)
        req = requests.get(url)
        time.sleep(0.5)
        soup = BeautifulSoup(req.text, "html.parser")
        # a = soup.select_one("a[href*=NCBI]")
        # print(a)
        # print(url)
        # for ultag in soup.find_all('ul', {'class': 'list-inline'}):
        #     print(ultag)
        try:
            target = soup.find('div', class_='ncbi-doc-id').find('ul')
            print(elem, target.text.split(":")[1].strip())
        except KeyError:
            print("Something went wrong about " + elem)