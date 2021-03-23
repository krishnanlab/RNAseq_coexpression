import numpy as np
import time
import argparse
import ast


parser = argparse.ArgumentParser()
parser.add_argument('-dir1',
                    default = '../data/small_data/BioGrid/', 
                    type = str,
                    help = 'The directory the intial edgelist sits in')
parser.add_argument('-edgelist',
                    default = 'BIOGRID-ORGANISM-Homo_sapiens-3.4.160.mitab_cleaned.edgelist', 
                    type = str,
                    help = 'The filenmae of the edgelist') 
parser.add_argument('-dir2',
                    default = '../data/large_data/BioGrid/', 
                    type = str,
                    help = 'The directory to save output') 
parser.add_argument('-binary', default = False, type = ast.literal_eval, 
	                help = 'Make a three column weighted edgelist into a binary adjaceny matrix') 
parser.add_argument('-diags', default = 0.0, type = float, 
	                help = 'What to make the diagonals')               
args     = parser.parse_args()
dir1     = args.dir1  
edgelist = args.edgelist
dir2     = args.dir2          
binary   = args.binary
diags    = args.diags

def read_edgelist(edgelist_file):
    '''
    This was a script from Anna, and Chris changed it a bit
    1. got rid of making it indexed by one
    2. added part to make into an array
    
    Reads an edgelist and returns an adjacency dictionary and corresponding
    IDs.
    Parameters
    ----------
        edgelist_file : str
            Path to a tab-separated file with lines formatted like
            gene_1 \\t gene_2 \\t edge_weight.
    Returns
    -------
        data : list(dict)
            A list of feature vector dicts where elements are accessed by node
            ID.
        IDs : list(str)
            Gene IDs corresponding to each row of `data`.
    '''
    data = []
    ID_map = {}
    current_node = 0

    print('Reading edgelist file {}...'.format(edgelist_file))
    with open(edgelist_file, 'r') as f:
        for line in f:
            if len(line.split('\t')) == 3:
                g1, g2, weight = line.split('\t')
                if binary == True:
                    weight = 1.0
                    txt1 = 'Binary'
                else:
                    weight = float(weight)
                    txt1 = 'Weighted'
            elif len(line.split('\t')) == 2:
                g1, g2 = line.split('\t')
                weight = 1.0
                txt1 = 'Binary'
            else:
                print('Not correct number of columns')

            if g1 not in ID_map:
                ID_map[g1] = current_node
                # dict is indexed by ID_map + 1 because
                # LibLINEAR segfaults when first feature is 0
                data.append({ID_map[g1]: 1})
                current_node += 1
            if g2 not in ID_map:
                ID_map[g2] = current_node
                data.append({ID_map[g2]: 1})
                current_node += 1

            # TODO: check for duplicate edges and print warning

            data[ID_map[g1]][ID_map[g2]] = weight
            data[ID_map[g2]][ID_map[g1]] = weight

    IDs = sorted(ID_map, key=ID_map.get)

    return data, IDs, txt1
    
def make_array(data,IDs):
    mysize = len(IDs)
    myarr = np.zeros((mysize,mysize))
    for idx, item in enumerate(data):
        col_inds = list(data[idx].keys())
        col_weights = list(data[idx].values())
        myarr[idx,col_inds] = col_weights
    np.fill_diagonal(myarr, diags)
    return myarr
    
def save_array(myarr,IDs,savepath_arr,savepath_nodelist):
    np.savetxt(savepath_arr,myarr)
    np.savetxt(savepath_nodelist,IDs,fmt='%s')
    
if __name__ == '__main__':
    tic0 = time.time()
    edgelist_file = dir1+edgelist
    print('Making data dict and ID list')
    tic = time.time()
    data,IDs,txt1 = read_edgelist(edgelist_file)
    print('The function took',time.time()-tic,'seconds to run')
    print('There are',len(IDs),'nodes in the graph')
    print('Making the array')
    tic = time.time()
    myarr = make_array(data,IDs)
    print('The function took',time.time()-tic,'seconds to run')
    print('Saving the arrays')
    savepath_arr = dir2 + edgelist + '_AdjMatrix_%s_Diags%.2f.txt'%(txt1,diags)
    savepath_nodelist = dir2 + edgelist + '_AdjMatrix_%s_Diags%.2f.nodelist'%(txt1,diags)
    save_array(myarr,IDs,savepath_arr,savepath_nodelist)
    print('The function took',time.time()-tic,'seconds to run')
    print('The script took',time.time()-tic0,'seconds to run')
    
