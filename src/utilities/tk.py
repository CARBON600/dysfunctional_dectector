from os.path import exists,dirname
import os
import logging
from Bio import SeqIO
from bioservices.kegg import KEGG
import pandas as pd
import requests
from bs4 import BeautifulSoup
kegg_api = KEGG()


def check(ofile,cmd,name,dry_run=False,LOGGER=logging.getLogger('dysfunction_logger')):
    executed_cmd=[]
    if exists(ofile):
        LOGGER.debug(f'{name} existed')
    else:
        if not dry_run and not exists(dirname(ofile)):
            os.makedirs(dirname(ofile))
        LOGGER.debug(f'{name} run following command: \n {cmd}')
    if (not dry_run) and (not exists(ofile)):
        executed_cmd.append(cmd)
    elif dry_run:
        LOGGER.debug(f" '{cmd}' would not run.")
    return executed_cmd
        

def parse_gbk(gbk_file):
    info = {}
    records = list(SeqIO.parse(gbk_file, 'genbank'))
    for contig in records:
        all_fea = [_ for _ in contig.features if 'locus_tag' in _.qualifiers]
        for fea in all_fea:
            locus_tag = fea.qualifiers['locus_tag'][0]
            start = fea.location.start.real
            end = fea.location.end.real
            strand = fea.location.strand
            contig_name = contig.id
            info[locus_tag] = [contig_name, start, end, strand]
    return info

def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d : len(iter) + 1])
    return n_iter

def output_file(opath,df):
    if opath.endswith('.xlsx'):
        df.to_excel(opath)
    elif opath.endswith('.tsv') or opath.endswith('.tab'):
        df.to_csv(opath,sep='\t',index=1)
    elif opath.endswith('.csv'):
        df.to_csv(opath,sep=',',index=1)
    else:
        raise IOError(f'Unknown suffix for {opath}')

import re
from ete3 import Tree

brite_file = '/mnt/home-db/pub/protein_db/kegg/v20230301/brite/ko00002.keg'
def get_koid(row):
    return re.findall(r'K\d+',row)

def parse_name(row):
    idx = row[0]
    row = row[1:]
    if get_koid(row):
        return get_koid(row)[0].strip()
    if '<b>' in row or idx in 'AB':
        name = row.split('<b>')[-1].split(' ',1)[-1].replace('</b>','').strip()
        return name
    elif idx == 'C':
        name = row.split(' ',1)[-1].strip()
        return name
    else:
        name = row.split(' ',1)[-1].strip()
        #print(row)
        mid = name.split(' ')[0]
        return mid

level2depth = {'A':4,'B':3,'C':2,'D':1}
def parse_brite_file(brite_file):
    basic_tre = Tree()
    hier_infos = open(brite_file).read().strip().split('\n')
    br_name = hier_infos[0].split('\t')[-1]
    # if '+' in br_name:
    #     br_name = suppl_dict.get(br, 'MISSING')

    parent_node = ''
    for row in hier_infos:
        if row[0] not in 'ABCD':
            continue
        name = parse_name(row).strip()
        if not name:continue
        if row[0] == 'A':
            parent_node = basic_tre.add_child(name=name)
            parent_node.add_feature('level',row[0])
            continue
        ## trace back
        diff_num = level2depth[row[0]] - level2depth[parent_node.level] + 1
        if diff_num>=1:
            for _ in range(diff_num):
                parent_node = parent_node.up
        if row[0]!='D':
            parent_node = parent_node.add_child(name = name)
            parent_node.add_feature('level',row[0])
        else:
            _ = parent_node.add_child(name = name)
            _.add_feature('level',row[0])

    return basic_tre

# module_brite = parse_brite_file(brite_file)

def ko2gename(ko_list):
    ko2name={}
    for ko_num in ko_list:
        url = f"https://www.genome.jp/dbget-bin/www_btit?ko+-n+{ko_num}"
        response = requests.get(url)
        soup = BeautifulSoup(response.content, "html.parser")
        a_tag = soup.find('a', href=f"/entry/ko:{ko_num}")
        row_info = a_tag.parent.text.strip()
        start_index = row_info.find(f"ko:{ko_num}") + len(f"ko:{ko_num}")
        end_index = row_info.find("[")
        gename = row_info[start_index:end_index].strip().split(";", 1)[0]
        ko2name[ko_num]=gename
    return ko2name
