## This fork tries to lift the code to support TensorFlow >2.5

I fact it seams to be just outdated and abandoned.
The new way to do this is possibly something like that:

```
import anndata
import scanpy
import scvi
import scipy
import os
import gzip
import subprocess

def gzip_files_in_directory(path):
    # Make sure the path exists
    if os.path.exists(path):
        # Run the gzip command on all files in the directory
        subprocess.run(f"gzip {path}/*", shell=True, check=True)
        print(f"All files in {path} have been gzipped.")
    else:
        print(f"The specified directory {path} does not exist.")

# Usage example:

def impute_expression_data( adata, output_dir ):
    # Convert to CSR format
    adata.X = scipy.sparse.csr_matrix(adata.X)
    # Set up the model
    scvi.model.SCVI.setup_anndata(adata)
    model = scvi.model.SCVI(adata)
    # Train the model
    model.train()
    # Get the imputed expression values for cells
    imputed_expression = scipy.sparse.csr_matrix(model.get_normalized_expression())
    adata.X = imputed_expression  # Replace with imputed data
    os.makedirs(output_dir, exist_ok=True)
    # Write the sparse matrix to a compressed MTX file
    scipy.io.mmwrite(os.path.join(output_dir, 'matrix.mtx'), imputed_expression)
    # Extract barcodes (cell identifiers)
    barcodes = adata.obs_names.values
    # Write barcodes to a TSV file
    with open(os.path.join(output_dir, 'barcodes.tsv'), 'w') as f:
        for barcode in barcodes:
            f.write(f"{barcode}\n")
            
    # Extract feature names (gene identifiers)
    features = adata.var_names.values
    # Write features to a TSV file
    with open(os.path.join(output_dir, 'features.tsv'), 'w') as f:
        for feature in features:
            f.write(f"{feature}\n")
    
    gzip_files_in_directory( output_dir )


adata = scanpy.read_10x_mtx('YourDataFolder/filtered_feature_bc_matrix/')

impute_expression_data( data2, 'YourDataFolder_IMPUTED/filtered_feature_bc_matrix/')
```


## Deep count autoencoder for denoising scRNA-seq data

A deep count autoencoder network to denoise scRNA-seq data and remove the dropout effect by taking the count structure, overdispersed nature and sparsity of the data into account using a deep autoencoder with zero-inflated negative binomial (ZINB) loss function.

See our [manuscript](https://www.nature.com/articles/s41467-018-07931-2) and [tutorial](https://nbviewer.ipython.org/github/theislab/dca/blob/master/tutorial.ipynb) for more details.

### Installation

#### pip

For a traditional Python installation of the count autoencoder and the required packages, use

```
$ pip install dca
```

#### conda

Another approach for installing count autoencoder and the required packages is to use [Conda](https://conda.io/docs/) (most easily obtained via the [Miniconda Python distribution](https://conda.io/miniconda.html)). Afterwards run the following commands.

```
$ conda install -c bioconda dca
```

### Usage

You can run the autoencoder from the command line:

`dca matrix.csv results`

where `matrix.csv` is a CSV/TSV-formatted raw count matrix with genes in rows and cells in columns. Cell and gene labels are mandatory. 

### Results

Output folder contains the main output file (representing the mean parameter of ZINB distribution) as well as some additional matrices in TSV format:

- `mean.tsv` is the main output of the method which represents the mean parameter of the ZINB distribution. This file has the same dimensions as the input file (except that the zero-expression genes or cells are excluded). It is formatted as a `gene x cell` matrix. Additionally, `mean_norm.tsv` file contains the library size-normalized expressions of each cell and gene. See `normalize_total` function from [Scanpy](https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.normalize_total.html) for the details about the default library size normalization method used in DCA.

- `pi.tsv` and `dispersion.tsv` files represent dropout probabilities and dispersion for each cell and gene. Matrix dimensions are same as `mean.tsv` and the input file.

- `reduced.tsv` file contains the hidden representation of each cell (in a 32-dimensional space by default), which denotes the activations of bottleneck neurons.

Use `-h` option to see all available parameters and defaults.

### Hyperparameter optimization

You can run the autoencoder with `--hyper` option to perform hyperparameter search.
