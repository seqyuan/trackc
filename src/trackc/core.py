from pyscenic.aucell import aucell, aucell4r
from typing import Sequence, Type
from ctxcore.genesig import GeneSignature
from GSVA import gsva, gmt_to_dataframe
import pandas as pd
import os,sys
basedir = os.path.abspath(os.path.dirname(__file__))

def genesets2GeneSig(df: pd.DataFrame) -> Sequence[Type[GeneSignature]]:
    """
    Conver dataframe to GeneSig for AUCell. 

    :param df: A dataframe with columns ["name", "member", "description"].
               name        member     description
               signature1  gene1      signature1 description
               signature1  gene2      signature1 description
               signature2  gene4      signature2 description
    :return: GeneSignature list.
    """

    names = list(df["name"].unique())
    GeneSigs = [0]*len(names)
    
    for i, v in enumerate(names):
        genes = df[df["name"]==v]["member"].to_list()
        GeneSigs[i] = GeneSignature(name=v, gene2weight=genes)
        
    return(GeneSigs)

def metabolism_sigs(resources="KEGG") -> pd.DataFrame:
    """
    Get a set of default metabolism signature from KEGG or REACTOME. 
    :param resources: KEGG or REACTOME, default KEGG
    """
    file = "KEGG_metabolism_nc.gmt"
    if resources=="REACTOME":
        file = "REACTOME_metabolism.gmt"
    df = gmt_to_dataframe(os.path.join(basedir,'resources', file))

    return df

def sigc_auc(
		gene_exp, 
		GeneSigsDf: pd.DataFrame, 
		num_workers=4) -> pd.DataFrame:

    GeneSigs = genesets2GeneSig(GeneSigsDf)

    auc_mxt = aucell4r(
        gene_exp,
        GeneSigs,
        0.05,
        noweights=False,
        normalize=False,
        num_workers=num_workers)

    return auc_mxt

def sigc_score(
		ex_mtx: pd.DataFrame, 
		GeneSigs: pd.DataFrame, 
		method="AUCell", 
		num_workers=4) -> pd.DataFrame:
    """
    Get a set of signature score of a given gene expression matrix. 

    :param ex_mtx: The expression profile matrix. 
                   The rows should correspond to different cells, 
                   the columns to different genes (n_cells x n_genes).
    :param GeneSigs: A dataframe with columns ["name", "member", "description"].
                       name        member     description
                       signature1  gene1      signature1 description
                       signature1  gene2      signature1 description
                       signature2  gene4      signature2 description
    :param method: sinature score method [AUCell, GSVA, ssGSEA, ...] (default: AUCell).
    :param num_workers: The number of cores to use in AUCell method (default: 4).
    :return: A dataframe with cell signature score (n_cells x n_signatures).
    """

    sig_score_mat = pd.DataFrame()
    if method == "AUCell":
        sig_score_mat = sigc_auc(ex_mtx, GeneSigs, num_workers=num_workers)
    else:
        sys.stderr("Method {0} not support!".format(method))
        sys.exit(1)
    return sig_score_mat


