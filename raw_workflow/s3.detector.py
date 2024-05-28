"""

Use refined kegg annotation to calculate the completeness of each module and output some tables for further visualizations.


input: refined annotations
output: 
list of genes/locus/modules needed to be visualized


"""

import pandas as pd
from os.path import join,dirname,realpath
import pickle,click,os
from tqdm import tqdm
from collections import defaultdict
from dysfunctional_dectector.src.utilities.logging import *


ko2ec_file = '/mnt/home-db/pub/protein_db/kegg/v20230301/link/ko2ec'
ko2pathway_file = '/home-user/thliao/db/protein_db/kegg/v20230301/link/ko2pathway'
ko2module_file = '/home-user/thliao/db/protein_db/kegg/v20230301/link/ko2module'
def parse_linkfile(infile):
    other2name = {}
    ko2others = defaultdict(list)
    for row in open(infile).read().strip().split('\n')[1:]:
        rows = row.split('\t')
        ko,others = rows[0],rows[1]
        ko,others = ko.split(':')[1],others.split(':')[1]
        ko2others[ko].append(others)
        if len(rows)>2:
            other2name[others] = rows[2]
    ko2others = {ko:';'.join(others) for ko,others in ko2others.items()}
    
    return ko2others,other2name
ko2ec,ec2name = parse_linkfile(ko2ec_file)
ko2pathway,p2name = parse_linkfile(ko2pathway_file)
ko2module,m2name = parse_linkfile(ko2module_file)


def load_pickle(dictionary_location):
    pickle_in = open(dictionary_location,"rb")
    dictionary = pickle.load(pickle_in)
    return dictionary

def read_module2ko():
    data_folder = join(dirname(dirname(realpath(__file__))),'data')
    # Import all modules from dictionaries
    regular_modules = load_pickle(join(data_folder,"01.KEGG_Regular_Module_Information.pickle"))
    bifurcation_modules = load_pickle(join(data_folder,"02.KEGG_Bifurcating_Module_Information.pickle"))
    structural_modules = load_pickle(join(data_folder,"03.KEGG_Structural_Module_Information.pickle"))
    return regular_modules,bifurcation_modules,structural_modules

def main(indir,odir,gid,with_unannotated=False):
    refined_ko_infodf_path = f'{indir}/{gid}_refined_ko_info.tsv'
    ko_infodf = pd.read_csv(refined_ko_infodf_path,sep='\t',index_col=0)
    refined_ko_bindf_path = f'{indir}/{gid}_refined_ko_bin.tsv'
    ko_bindf = pd.read_csv(refined_ko_bindf_path,sep='\t',index_col=0)    
    subko_df = ko_infodf.loc[ko_infodf[gid]!='no KEGG-annotated',:]

    regular_modules,bifurcation_modules,structural_modules = read_module2ko()

    mentioned_kos = list(subko_df.index)
    related_modules = list(set([m for ko in mentioned_kos for m in ko2module.get(ko,'').split(';') if m]))
    logger.debug(f"Identified {len(related_modules)} modules in {gid} annotations...")
    logger.debug(f"Ouputing module pattern into {odir}/module_tables/")
    os.makedirs(odir,exist_ok=True)
    os.makedirs(f"{odir}/module_tables/",exist_ok=True)
    module2info = defaultdict(dict)
    for m in tqdm(related_modules):
        module_name = m2name[m]
        module2info[m]['Name'] = module_name
        if m not in regular_modules:
            continue
        steps_ko = regular_modules[m]
        ko2step = {}
        ordered_kos = []
        for _,step in steps_ko.items():
            step_name = f"step {_}"
            remaining_kos = [_ for _ in step if _ in subko_df.index]
            if len(remaining_kos)==0:
                remaining_kos = step
            ordered_kos.extend(remaining_kos)
            for k in remaining_kos:
                ko2step[k] = step_name
        module_ko_df = ko_infodf.T.reindex(columns=ordered_kos)
        module_ko_df.loc['STEP INFO',:] = [ko2step[_] for _ in ordered_kos]
        
        module2info[m]['Pattern'] = module_ko_df

        dfout = f"{odir}/module_tables/{m}_info.tsv"
        module_ko_df.to_csv(dfout,sep='\t',index=1)
        module_ko_df = ko_bindf.T.reindex(columns=ordered_kos)
        module_ko_df.loc['STEP INFO',:] = [ko2step[_] for _ in ordered_kos]        
        dfout = f"{odir}/module_tables/{m}_val.tsv"
        module_ko_df.to_csv(dfout,sep='\t',index=1)        


# parse args
@click.group()
@click.option('--odir','-o',help="output directory")
@click.pass_context
def cli(ctx,odir):
    ctx.ensure_object(dict)
    ctx.obj['odir'] = realpath(odir)

@cli.command()
@click.option('--genome','-gid',help="Genome ID")
@click.option('--indir','-i',help="Input directory which is also a output directory of s2.")
@click.pass_context
def workflow(ctx,genome,indir):
    outputdir = ctx.obj['odir']
    indir = realpath(indir)
    genome = genome
    main(indir,outputdir,genome)

######### RUN
if __name__ == '__main__':
    cli()
