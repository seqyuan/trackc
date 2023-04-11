import cooler
import pandas as pd
import numpy as np
from typing import Union, List, Sequence, Tuple, Collection, Optional
from .._utils import GenomeRegion

class RegionsCMat:
    def __init__(self, cmat: np.array, row_regions: pd.DataFrame, col_regions: pd.DataFrame):
        self.cmat = cmat
        self.row_regions = row_regions
        self.col_regions = col_regions

def get_regions_contact_mat(
        coolMat: cooler.Cooler, 
        balance: bool = False,
        row_regions: Union[Sequence[str], str, None] = None,
        col_regions: Union[Sequence[str], str, None] = None, 
        ) -> RegionsCMat:
    """\
    Extract a set of regions matrix from the cool format Hi-C matrix.

    The extracted matrix will splice intra and Inter region interaction according to 
        the given order and direction of the regions.

    Parameters
    ----------
    coolMat 
        ``cooler.Cooler``: cool format Hi-C matrix (https://github.com/open2c/cooler)

    balance
        ``bool``: The ``'balance'`` parameters of ``coolMat.matrix(balance=False).fetch('chr6:119940450-123940450')``

    row_regions
        ``chrom region`` list: or ``chrom region`` or None. 
        The subset matrix row genome regions
        eg. ``"chr6:1000000-2000000"``, eg. ``["chr6:1000000-2000000", "chr3:5000000-4000000", 'chr5']``
        The start can be larger than the end (eg. ``"chr6:2000000-1000000"``), 
            which means you want to get the reverse region contact matrix

    col_regions
        ``chrom region`` list: or ``chrom region`` or None. 
        The subset matrix col genome regions, default is ``None``, which means the sample region as ``row_regions``
        
    Returns:
        np.array: row_regions and col_regions contact matrix
    """
    # -------
    row_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in row_regions])
    if col_regions == None:
        col_GenomeRegions = row_GenomeRegions.copy()
    else:
        col_GenomeRegions = pd.concat([GenomeRegion(i).GenomeRegion2df() for i in col_regions])
    
    # ------
    region_mat_dic = {}
    for _, row_row in row_GenomeRegions.iterrows():
        for _, col_row in row_GenomeRegions.iterrows():
            row_col_region_cmat = coolMat.matrix(balance=balance).fetch(row_row['region4coolFetch'], col_row['region4coolFetch'])
            if row_row['isReverse'] == True:
                row_col_region_cmat = np.flip(row_col_region_cmat, 0)
            if col_row['isReverse'] == True:
                row_col_region_cmat = np.flip(row_col_region_cmat, 1)

            region_mat_dic["{0}__{1}".format(row_row["region4coolFetch"], col_row["region4coolFetch"])] = row_col_region_cmat

    vstack_list = [None]*row_GenomeRegions.shape[0]
    hstack_list = [None]*col_GenomeRegions.shape[0]
    
    for i, row_region in enumerate(row_GenomeRegions["region4coolFetch"].to_list()):
        for ii, col_region in enumerate(col_GenomeRegions["region4coolFetch"].to_list()):
            hstack_list[ii] = region_mat_dic["{0}__{1}".format(row_region, col_region)]
        
        vstack_list[i] = np.hstack(tuple(hstack_list))
                          
    cMat = np.vstack(tuple(vstack_list))

    #### 
    row_GenomeRegions['cbins'] = 0
    col_GenomeRegions['cbins'] = 0

    for i in row_GenomeRegions.index:
        for ii in col_GenomeRegions.index:
            cmat_shape = region_mat_dic["{0}__{1}".format(
                row_GenomeRegions.loc[i, "region4coolFetch"], 
                col_GenomeRegions.loc[ii, "region4coolFetch"])].shape
            
            row_GenomeRegions.loc[i, "cbins"] = cmat_shape[0]
            col_GenomeRegions.loc[i, "cbins"] = cmat_shape[1]
            
    r_l_regions_cMat = RegionsCMat(cmat=cMat, row_regions=row_GenomeRegions, col_regions=col_GenomeRegions)
    
    return r_l_regions_cMat
        

