import h5py
from SIR_parameters import SIRparameters, set_atrophy, set_netw
from SIR_genes import SIRgenes
from sir_iterator import sir_iterator
from write_async import write_async


def sir_main(clear, opt=None):
    """
    Main function to run SIR simulation for gene pairs.

    Parameters:
        clear (str): List of clearance gene x region string arrays.
        opt (dict): Optional parameters with default values.
    """
    if opt is None:
        opt = {}

    # Setting defaults for optional parameters
    opt.setdefault('risk', None)
    opt.setdefault('out', '../SIR_simulator_gene_corrs/results_3/gene_corrs.csv')
    opt.setdefault('parc', 'S132')
    opt.setdefault('dbg', False)
    opt.setdefault('vis', False)
    opt.setdefault('null', 'none')

    # Load workspace variables from HDF5 file
    ws_file = f'data/workspace_{opt["parc"]}.h5'

    with h5py.File(ws_file, 'r') as f:
        gene_expr = f['gene_expr'][:]
        ifod_len_35 = f['ifod_len_35'][:]
        ifod_den_35 = f['ifod_den_35'][:]
        ROIsize = f['ROIsize'][:]
        bgs = f['bgs'][:]
        cobre = f['cobre'][:]
        hcpep = f['hcpep'][:]
        stages = f['stages'][:]

    # Initialize parameter object and set attributes
    params = SIR_parameters(opt)
    params = set_atrophy(params, bgs, cobre, hcpep, stages)
    params = set_netw(params, ifod_len_35, ifod_den_35, ROIsize)

    # Initialize gene expression object
    genes = SIRgenes(opt, clear, gene_expr)

    # Run SIR simulation for each gene pair
    gene_corrs = sir_iterator(params, genes)

    # Output the results
    if opt['dbg']:
        print(gene_corrs)
    else:
        write_async(gene_corrs, opt['out'])