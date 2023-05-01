import pandas as pd

def vj_parsing(
    sequences, 
    imgt_correction = True,
    cdr3_column = "junction_aa",
    v_gene_column = "v_call",
    j_gene_column = "j_call"
    ):
    
    # Only retain beta chain sequences.
    sequences = sequences[sequences[j_gene_column].str.contains('TRBJ') &
                          sequences[v_gene_column].str.contains('TRBV')]
    
    # Format the V/J genes
    sequences['J_family'] = (sequences[j_gene_column].str.replace('TRBJ', '')
                             .str.split('-').str[0].str.zfill(2))
    sequences['V_family'] = (sequences[v_gene_column].str.replace('TRBV', '')
                             .str.split('-').str[0].str.zfill(2))
    sequences['J_gene'] = sequences[j_gene_column].apply(
        lambda gene_name: '-'.join([sub_gene.zfill(2) for sub_gene in
                                    gene_name.replace('TRBJ', '').split('-')]))
    sequences['V_gene'] = sequences[v_gene_column].apply(
        lambda gene_name: '-'.join([sub_gene.zfill(2)
                                    if sub_gene not in 'ABC'
                                    else sub_gene for sub_gene in
                                    gene_name.replace('TRBV', '').split('-')]))
    if imgt_correction:
        # Get all IMGT genes and families.
        imgt_v_genes, imgt_j_genes, imgt_v_families, imgt_j_families =\
            _get_imgt_genes()

        # Only keep sequences with IMGT families.
        sequences = sequences[(sequences['J_family'].isin(imgt_j_families)) &
                              (sequences['V_family'].isin(imgt_v_families))]

        # Replace non-IMGT genes by IMGT genes.
        gene_map = {'20': '20-01',
                    '21': '21-01',
                    '22': '22-01',
                    '23': '23-01',
                    '24': '24-01',
                    '25': '25-01',
                    '29': '29-01',
                    '0A': 'A',
                    '0B': 'B',
                    '0C': 'C'}
        sequences['J_gene'] = sequences['J_gene'].apply(
            lambda x: x if x in imgt_j_genes
            else x.split('-')[0].zfill(2))
        sequences['V_gene'] = sequences['V_gene'].apply(
            lambda x: x if x in imgt_v_genes
            else gene_map.get(x.split('-')[0].zfill(2),
                              x.split('-')[0].zfill(2)))

    target_cols = {'V_gene':'TRBV', 
                   'J_gene':'TRBJ',
                   'V_family':'TRBV', 
                   'J_family':'TRBJ'}
    
    for column in target_cols.keys():
        sequences[column] = target_cols[column] + sequences[column]

    return sequences


def _get_imgt_genes():
    # http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=human&group=TRBV
    v_genes_families = pd.read_csv('./data/imgt_gene_parsing/v_genes.tsv', sep='\t', dtype=str)
    v_families = set(v_genes_families['V_family'])
    v_genes = set(v_genes_families['V_gene'])

    # http://www.imgt.org/IMGTrepertoire/index.php?section=LocusGenes&repertoire=genetable&species=human&group=TRBJ
    j_genes_families = pd.read_csv('./data/imgt_gene_parsing/j_genes.tsv', sep='\t', dtype=str)
    j_families = set(j_genes_families['J_family'])
    j_genes = set(j_genes_families['J_gene'])

    return v_genes, j_genes, v_families, j_families


def merge_results(
        original: pd.DataFrame,
        predictions: pd.DataFrame,
        clusters: pd.DataFrame
        ):
    
    clusters = clusters.merge(right = original, on = ["junction_aa", "v_call"])
    clusters = vj_parsing(clusters)
    clusters["tcr_id"] = clusters["V_gene"] + "_" + clusters["junction_aa"] + "_" + clusters["J_gene"]
    
    predictions["tcr_id"] = predictions["TRBV_gene"] + "_" + predictions["CDR3_beta"] + "_" + predictions["TRBJ_gene"]
    merged = clusters.merge(right = predictions, on = "tcr_id", how = "outer")
    
    return merged
    
def create_edgelist_vgene(clusters):
    '''
    Create tab-separated edgelist of edges with HD = 1, from a set of sequences.    
    '''
    
    # Set makes sure there are no dupes
    tcrs = [(clusters.iloc[i]["junction_aa"], 
             clusters.iloc[i]["v_call"]) for i in range(len(clusters))]
    # tcrs = set(tcrs)
    
    # Hashing
    cdr3hash = dict()
    for tcr in tcrs:
        cdr = tcr[0]
        for hash in (cdr[::2], cdr[1::2]):
            if hash not in cdr3hash:
                cdr3hash[hash] = set()
            cdr3hash[hash].add(tcr)
            
    # Generate network
    edgelist = set()
    for hash in cdr3hash:
        if len(cdr3hash[hash]) >= 1:
            for tcr1 in cdr3hash[hash]:
                for tcr2 in cdr3hash[hash]:
                    if tcr1 != tcr2:
                        if tcr1[0] <= tcr2[0]:
                            if sum(ch1 != ch2 for ch1, ch2 in zip(tcr1[0], tcr2[0])) <= 1:
                                edgelist.add((tcr1,tcr2))
    
    edges = pd.DataFrame(edgelist, columns=["source", "target"])
    for col in edges:
        edges[col] = edges[col].apply(lambda x: "_".join(x))

    return edges