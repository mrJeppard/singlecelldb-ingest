"""
Author Duncan McColl dmccoll@ucsc.edu

Main script pulls all EBI expression atlas data into scanpy objects. If you change this code you'll likely
break the ebi ingest script.

To run, source an environment built from requirements.txt and change the TMPDIR and DATADIR globals
below the import statements. Currently the script has no arguments.

`python anndataize_ebi.py`

The data_ingest_pipeline function builds a single scanpy object from an accession number.

The all_accession_and_descriptions function returns two lists containing all the expression accessions
and their "experiment description".

-------------------------

A single scanpy object is built from the data_ingest_pipeline. That anndata object has the following fields.

.obs["sc3-K-*"] where * is all the values of K that were run.
.uns["sc3_preferred_cluster"] = str a sc3-K-* identifier in .obs
.uns["sc3_cluster_solutions"] = list (all the names of the sc3 cluster solutions)
.uns["species"] = str
.uns["idf_data_identity"] = dictionary with keys and then the line of text associated with them.
.uns["ebi_view_data_url"] = str
.uns["short_description"] = str

"""

from io import BytesIO
import os
import pandas as pd
import numpy as np
import requests
import scanpy as sc
import sys
from shutil import rmtree
import urllib.request
import zipfile

# Temp directory that stores unzipped .mtx files, both the files and the directory
# are removed after use.
TMPDIR = "/projects/sysbio/users/cellAtlas/data/EBI-starter-scanpyObjs/tmp"
# This directory will be filled with ebi anndata objects.
DATADIR = "/projects/sysbio/users/cellAtlas/data/EBI-starter-scanpyObjs"

# This is the url that gets all of the accession ids (used for endpoint generation) and
# experimental descriptions.
EBI_DATAENDPOINT = "https://www.ebi.ac.uk/gxa/sc/json/experiments/"

# These are tab filetypes that can be accessed through ebi (they are used to build the urls in tab_file_url)
# Could add type checking to functions but not currently used.
VALID_TAB_FILETYPES = {"sdrf", "idf", "experiment-design", "cluster"}

VALID_ZIP_FILETYPES = {"marker-genes", "quantification-filtered"}


def view_data_url(expression_accession):
    """Creates the url to view the data on the ebi's single cell expression atlas browser"""
    url = "https://www.ebi.ac.uk/gxa/sc/experiments/%s/Results" \
          % expression_accession
    return url


def tab_file_url(expression_accession, filetype):
    """Generate a url for access of a tab file of the EBI's filetypes"""

    url = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download?fileType=%s&accessKey=" \
          % (expression_accession, filetype)
    return url


def mtx_zip_url(expression_accession):
    """Make the url to grab the zipped mtx experssion data from the ebi's accessionI_id"""
    url = "https://www.ebi.ac.uk/gxa/sc/experiment/%s/download/zip?fileType=quantification-filtered&accessKey=" \
        % expression_accession
    return url


def open_url(url):
    return urllib.request.urlopen(url)


def mtx_to_scanpy(url, tmpdir):
    """
    Build a scanpy object with the .mtx files from an ebi scAtlas endpoint.
    Caution, creates a tmp directory in ./tmp then deletes the
    contents right afterwards.
    :param url: single cell expression atlas url
    :param tmpdir: path for temporary write of *mtx* contents.
    :return: scanpy object with .X and .obs_names and var_names filled.
    """
    opened_url = open_url(url)
    # Dump the zipped files into a temporary folder.
    zipped = zipfile.ZipFile(BytesIO(opened_url.read()))
    os.mkdir(tmpdir)
    zipped.extractall(tmpdir)

    # Grab the three filenames for access.
    filenames = zipped.namelist()

    # Get all the ./tmp file names
    mtxfile = [f for f in filenames if f.endswith(".mtx")][0]
    mtxfile = os.path.join(tmpdir, mtxfile)
    colsfile = [f for f in filenames if f.endswith(".mtx_cols")][0]
    colsfile = os.path.join(tmpdir, colsfile)
    rowsfile = [f for f in filenames if f.endswith(".mtx_rows")][0]
    rowsfile = os.path.join(tmpdir, rowsfile)

    # Fill the anndata object.
    anndata = sc.read_mtx(mtxfile).transpose()
    anndata.obs_names = pd.read_csv(colsfile, header=None, sep="\t")[1]
    anndata.var_names = pd.read_csv(rowsfile, header=None, sep="\t")[1]

    # Empty out the tmp dir...
    rmtree(tmpdir)
    return anndata


def add_sc3_clusters_as_obs(url, anndata):
    """
    Add the sc3 clusters and the 'preferred' cluster to the anndata object.
    :param url: 
    :param anndata: 
    :return: modifies anndata inplace
        anndata.obs["sc3-K-*"] where * includes all different values for K.
        list(str) anndata.uns["sc3_cluster_solutions"] 
            The names 
        str or None anndata.uns["sc3_preferred_cluster"]   
    """
    # Grab data and wrangle
    clusters = pd.read_csv(url, sep="\t")
    new_column_descript = [("sc3-K-%s" % str(k)) for k in clusters["K"]]
    preferred_cluster = np.array(new_column_descript)[clusters["sel.K"].values]

    if len(preferred_cluster) == 0:
        preferred_cluster = None
    else:
        preferred_cluster = preferred_cluster[0]

    sc3df = clusters[clusters.columns[2:]].transpose()
    sc3df.columns = new_column_descript
    
    # Add each cluster solution to the anndata obj.
    for column in sc3df.columns:
        anndata.obs[column] = sc3df[column]

    # Add unstructured data.
    anndata.uns["sc3_cluster_solutions"] = new_column_descript
    anndata.uns["sc3_preferred_cluster"] = preferred_cluster

    return None


def add_exp_attrs(url, anndata):
    """
    Adds 'experimental-design' data to the scanpy object.
    All column of the experimental-design are copied over to anndata.obs
    :param url:
    :param anndata:
    :return: modifies anndata.obs inplace.
    """

    expr_attrs = pd.read_csv(url, index_col=0, sep="\t")
    expr_attrs.head()
    for column in expr_attrs.columns:
        anndata.obs[column] = expr_attrs[column]

    return None


def add_id_data(url, anndata):
    """Adds 'idf' data to the anndata object. This is information about the data's identity.
    It is added to anndata.uns["idf_data_identity"]" as a dictionary."""

    # Ugly hack for reading these into a dict. Some bizarre string encoding thing
    # That pandas apparently deals with.
    df = pd.read_csv(url, sep="\t")

    idf_dict = {}
    #resource = open_url(url).read()
    #charset = chardet.detect(resource)["encoding"]
    for row in df.index:
        rlist = df.loc[row].tolist()
        rlist = [r for r in rlist if isinstance(r, str)]
        keyval = " ".join(rlist).strip().split(" ", 1)
        try:
            idf_dict[keyval[0]] = keyval[1]
        except IndexError:
            pass

    anndata.uns["idf_data_identity"] = idf_dict


def add_view_data_url(url, anndata):
    """Adds the url that can be ued to observe the data from their website."""
    anndata.uns["view_data_url"] = url


def add_srdf_data(url, anndata):
    raise NotImplementedError


def data_ingest_pipeline(accession_id, experiment_desc, tmpdir, species, add_idf=True):
    """Build a scanpy anndata object from 'quantification-filtered' matrix, 'idf' notes,
    sc3 clusters, and 'experiment-design' details from their prospective endpoints."""
    # Grab  the zipped mtx file and put it in an anndata object.
    url = mtx_zip_url(accession_id)
    anndata = mtx_to_scanpy(url, tmpdir)

    # Add the sc3 clusters as .obs
    url = tab_file_url(accession_id, "cluster")
    add_sc3_clusters_as_obs(url, anndata)

    # And the attributes that come along with the experiment to the scanpy object
    url = tab_file_url(accession_id, "experiment-design")
    add_exp_attrs(url, anndata)

    if add_idf:
        #add unstructured details about the experiment
        url = tab_file_url(accession_id, "idf")
        add_id_data(url, anndata)

    # Add the url to view the data on ebi's single cell expression atlas website.
    url = view_data_url(accession_id)
    add_view_data_url(url, anndata)

    anndata.uns["short_description"] = experiment_desc
    anndata.uns["species"] = species

    return anndata


def write_anndata(expressionAccession, datadir, filenamepostfix, anndata):
    filename = "%s%s" % (expressionAccession, filenamepostfix)
    filepath = os.path.join(datadir, filename)
    anndata.write(filepath)


def all_high_level(ebi_dataendpoint=EBI_DATAENDPOINT):
    """
    Returns accession ids, experimental description, and species for each of the
    datasets available in the ebi's scea.
    :param ebi_dataendpoint: the string url needed to get the data.
    :return: (list, list, list) accessions, descriptions, species
    """
    # Grab json and access array.
    expressionMeta = requests.get(ebi_dataendpoint).json()["aaData"]
    # Unpack meaningful contents.
    accessions = [dict["experimentAccession"] for dict in expressionMeta]
    descriptions = [dict["experimentDescription"] for dict in expressionMeta]
    species = [dict["species"] for dict in expressionMeta]

    return accessions, descriptions, species


def add_long_description(accession_id, anndata):
    try:

        long_description = anndata.uns["idf_data_identity"]["Experiment Description"]

        anndata.uns["long_description"] = long_description
    except KeyError:
        print("there was a key error for %s" % accession_id)

    return None


def main():
    sc.settings.verbosity = 0
    # The path to write out temporary zipped files from the mtx endpoint.
    tmpdir = TMPDIR

    # The path to dump each of the generated scanpy objects to.
    datadir = DATADIR

    # This is the string that is used for naming scanpy h5ad files.
    filenamepostfix = "_ebi.expdesign.idf.expr.h5ad"

    accessions, descriptions, all_species = all_high_level()
    
    try:
        os.mkdir(datadir)
    except FileExistsError:
        pass

    for accession, description, species in zip(accessions, descriptions, all_species):
        print(accession)
        anndata = data_ingest_pipeline(accession, description, tmpdir, species)
        try:
            write_anndata(accession, datadir, filenamepostfix, anndata)
        except UnicodeEncodeError:
            print("Couldn't write %s because of unicode error" % accession)
            print("Trying without the idf file....")
            anndata = data_ingest_pipeline(accession, description, tmpdir, species, add_idf=False)
            write_anndata(accession, datadir, filenamepostfix, anndata)


if __name__ == "__main__":
    sys.exit(main())
