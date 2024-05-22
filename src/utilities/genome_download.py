import os
import subprocess
import requests
from bs4 import BeautifulSoup
import gzip
from tqdm import tqdm
from urllib.request import urlopen
from urllib.parse import urlencode
from collections import defaultdict
from Bio import SeqIO
import pandas as pd

# from dysfunctional_dectector.src.utilities.logging import *

# input:target organism search url output:list of organism url
def get_search_result(organism_name):
    url = f"https://www.genome.jp/dbget-bin/www_bfind_sub?mode=bfind&max_hit=1000&locale=en&serv=kegg&dbkey=genome&keywords={organism_name}&page=1"
    response = requests.get(url)
    soup = BeautifulSoup(response.content, "html.parser")
    link_element = [a['href'] for a in soup.find_all("a", href=lambda href: href and "/entry/gn:" in href)]
    if link_element:
        return link_element
    else:
        raise ValueError("Failed to get genomes of the organism from KEGG database, please check the name of the organism.")

# input: list of organism url result and download folder result:genome compressed file downloaded
def download_genome(url_list, temp_folder):
    path_dict = defaultdict(list)
    download_result = defaultdict(list)
    download_folder = f"{temp_folder}/accessory_genome"
    os.makedirs(download_folder, exist_ok=True)

    ready_list = url_list[:5] if len(url_list) >= 5 else url_list

    for url in tqdm(ready_list):
        compl_link = f"https://www.genome.jp{url}"
        response = requests.get(compl_link)
        soup = BeautifulSoup(response.content, "html.parser")
        link_element_faa = [a['href'] for a in soup.find_all("a", href=lambda href: href and "https://ftp.ncbi.nlm.nih.gov/genomes/all" in href)]
        #print(link_element_faa)
        if not link_element_faa:
            continue
        
        genome_link = link_element_faa[0]
        response_faa = requests.get(genome_link)
        soup_faa = BeautifulSoup(response_faa.content, "html.parser")
        # parent_url = 
        # for link in soup_faa.find_all('a'):
        #     href = link.get('href')
        #     if href and not href.startswith(base_url):
        #         absolute_url = urljoin(url, href)
        for link in soup_faa.find_all('a'):
            href = link.get('href')
            if href and href.endswith('_translated_cds.faa.gz'):
                absolute_url = href

        #'/entry/GCA_000011965.2_ASM1196v2_translated_cds.faa.gz'
        #https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/011/965/GCA_000011965.2_ASM1196v2/GCA_000011965.2_ASM1196v2_translated_cds.faa.gz
        if not absolute_url:
            continue
        absolute_url = genome_link+"/"+absolute_url
        try:
            subprocess.run(['wget', '-P', download_folder, absolute_url], check=True)
        except subprocess.CalledProcessError as e:
            print(f'Genome download fails: {e}')
            continue

        gz_file = absolute_url.split("/")[-1]
        faa_path = f"{download_folder}/{gz_file}"
        with gzip.open(faa_path, "rb") as file_in:
            extracted_file_name = ".".join(gz_file.split(".")[:-1])
            out_path = f"{download_folder}/{extracted_file_name}"
            with open(out_path, "wb") as file_out:
                file_out.write(file_in.read())

        link_element_gbk = [a['href'] for a in soup.find_all("a", href=lambda href: href and "https://www.ncbi.nlm.nih.gov/nuccore" in href)]
        
        if not link_element_gbk:
            continue
        
        access_num = link_element_gbk[0].split('/')[-1]
        base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        params = {
            "db": "nucleotide",
            "id": access_num,
            "rettype": "gb",
            "retmode": "text"
        }

        download_url = f"{base_url}?{urlencode(params)}"
        os.makedirs(f"{download_folder}/{access_num}", exist_ok=True)
        try:
            gb_path = f"{download_folder}/{access_num}/{access_num}.gb"
            with urlopen(download_url) as response, open(gb_path, "wb") as file:
                file.write(response.read())
        except Exception as e:
            print(f'Genome download fails: {e}')
            continue

        gbk_path = f"{download_folder}/{access_num}/{access_num}.gbk"
        os.rename(gb_path, gbk_path)
        ###write the path of gbk and faa file to one tsv file
        record = SeqIO.read(gbk_path,"genbank")
        genome_id = record.id
        path_dict["genome Name"].append(genome_id)
        path_dict["protein file"].append(out_path)
        path_dict["gbk file"].append(gbk_path)
    df=pd.DataFrame(path_dict)
    table_path = f"{download_folder}/genome_info_table.tsv"
    df.to_csv(table_path, sep='\t',index=False)
    return table_path

# link_list = get_search_result("ruegeria")
# result = download_genome(link_list,"/mnt/storage3/yfdai/download/script/op/temp_file")
# print(result)

