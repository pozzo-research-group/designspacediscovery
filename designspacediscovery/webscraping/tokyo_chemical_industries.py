"""
Scrape tci ids off of tci website

https://www.tcichemicals.com/US/en

"""
import urllib
import bs4
from bs4 import BeautifulSoup
import time
import re


def scrape_tci_productids(url_base):
    """
    Retrieves TCI product numbers for a particular class of chemicals. These product numbers can then be parsed to pubchem CIDs using the pubchem identifier exchange service. 

    Parameters:
    ----------
    url_base (str): base URL to scrape from. Should be in format 'www.tcichemicals.com/lots of junk/page='. This function fills in page nos

    Returns:
    -------
    tci_ids, cas_nos (sets): TCI product id and CAS no sets.
    """
    more_pages = True
    page_num = 0
    
    
    tci_ids = set()
    cas_nos = set()
    
    while more_pages:
        print('Scraping page number ', page_num)
        url = url_base+str(page_num) # tci indexes from 0
        req = urllib.request.Request(url, headers={'User-Agent':"Magic Browser"})
        resp = urllib.request.urlopen(req)
        bytespage = resp.read()
        page = bytespage.decode('UTF-8')
        pagelines = page.split('\n')

        # find the lines with the product ids in them
        match_strings = []
        for line in pagelines:
            m = re.search('data-product-code1=', line)
            if m is not None:
                match_strings.append(m.string)

        #parse the product ids and cas nos
        more_pages = eval_more_pages(page)
        
        time.sleep(0.25)

        for match in match_strings:
            prod = parse_tci_productdiv(match)

            try:
                tci_ids.add(prod['tci_id'])
                #print('appended tci ', prod['tci_id'])
                cas_nos.add(prod['cas_no'])
            except:
                print('Error with extracting cas or prod no.')
                
        page_num += 1
        
    return tci_ids, cas_nos


def parse_tci_productdiv(div):
    """
    regex to get product id out of their divs 
    """

    for chunk in div.split(' '):
        if re.match('data-product-code1', chunk):
            parts = chunk.split('=')
            pid = parts[-1]
            pid = re.sub('[^a-zA-Z0-9]', '', pid)

        if re.match('data-casNo', chunk):
            parts = chunk.split('=')
            cas = parts[-1]
            cas = re.sub('[^a-zA-Z0-9-]', '', cas)
            
    return {'tci_id':pid, 'cas_no':cas}

def eval_more_pages(page):
    """
    see if more pages to review. True if so, F otherwise
    """
    soup = BeautifulSoup(page, features='lxml')
    pagenope = soup.find_all("li", class_='pagination-next disabled')
    
    if len(pagenope) > 0:
        return False
    else:
        return True