#!/usr/bin/env python3
"""
this module provides functions to search ENA for samples by sample ID, which are used in the 
Jupyter notebook pf8genre.ipynb in the same folder.  

It also provides function for ENA data download, which can be used interactively or through a 
command line interface.  

To use the module, make sure all dependencies are installed. A file with dependencies is provide in this folder. To install dependencies, run the 
following command:
pip install requirements.txt 

The jupyter notebook provides usage examples. To run the command-line tool, consult the help 
output by running:

python3 ENA_data_helper.py download --help

"""
import pandas as pd
import requests
import argparse
import csv
import time
import pathlib
from io import StringIO
import urllib.request

def args_parser():     
    parser = argparse.ArgumentParser(
        description = "ENA_data_helper.py: helper tool for downloading FASTQ from ENA") 

    # Currently, there is only one command available (download). 
    # The rest of the functionality of this module is currently meant to be used in 
    # interactive Jupyter notebooks but could be accessed via sub-commands in the CLI too.  
    subparsers = parser.add_subparsers()
    subparsers.required = True
    subparsers.dest = 'command'
    
    download_parser = subparsers.add_parser(
        'download', help='download command')
    download_parser.set_defaults( func = cli_download_fastqs )

    download_parser.add_argument(
        '--insdc_manifest','-i',
        action = 'store',
        required = True,
        metavar = 'FILE', 
        type = str,
        help='path to the insdc manifest data file.')

    download_parser.add_argument(
        '--output_folder','-o',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to a directory (must not exist yet) where downloaded files are deposited')

    download_parser.add_argument(
        '--top3',
        action = 'store_true',
        required = False,
        help = 'when used, only the data for the first 3 rows in the data file are downloaded. Use this first to ensure that the output is what you are expecting before proceeding to a full download' )

    download_parser.add_argument(
        '--skip_errors',
        action = 'store_true',
        help ='when used, rows of data that would otherwise throw an exception are just skipped and a message is printed'
    )

    download_parser.add_argument(
        '--download_attempts',
        action = 'store',
        default= 3,
        required = False, 
        type = int,
        help='Optional number of attempts to try FTP downloads. Defaults to 3')

    return parser

def create_ena_data_frame_from_bioproject(bioproject:str):
    """
    Takes a BioProject id, and returns a new dataframe with ENA run data
    
    Args:
        bioproject (str): INSDC Bioproject

    Returns:
        pandas.DataFrame: A new DataFrame that contains the input sample IDs as well as the  
        following data fields from ENA:
        run_accession,study_accession,sample_accession,sample_title,submitted_ftp

    """  
    fields = ['study_accession','sample_accession','sample_title','experiment_accession','run_accession','submitted_ftp']
    query = _build_ena_bioproject_query(bioproject=bioproject)
    ena_data = _parse_ena_response( _search_ena(query , return_fields=fields, limit=0) )

    return ena_data

def create_ena_data_frame_from_samples(df:pd.DataFrame, sample_id_col_name:str='sample', chunk_size:int=50):
    """
    Takes a DataFrame with a sample ID column (default name 'sample') and searches 
    ENA by sample ID ('sample_title' field in ENA). Returns a new dataframe with sample ID 
    and ENA data. A single sample may have multiple ENA records, hence sample ID is not 
    unique in this DataFrame.  

    Args:
        df (pandas.DataFrame): A pandas DataFrame. Only needs a single column for sample IDs.  
            By default, this is expected to be called 'sample' but can be changed. Other data  
            is ignored.  
        sample_id_col_name (str): Name of the sample ID column, defaults to 'sample'
        chunk_size (int): number of samples queried in one request. If the total number of 
            samples IDs is larger than this, the ENA query will be run in chunks.  
            Defaults to 50.  

    Returns:
        pandas.DataFrame: A new DataFrame that contains the input sample IDs as well as the  
        following data fields from ENA:
        'run_accession','center_name','library_strategy','sample_accession','fastq_ftp','submitted_ftp'

    """
    try:
        sample_ids=df[sample_id_col_name].tolist()
    except KeyError:
        raise ValueError(
            f'DataFrame is missing sample ID column, expected name {sample_id_col_name} ' +
            '(can be changed with param "sample_id_col_name")')
    
    fields = ['run_accession','experiment_accession', 'study_accession', 'sample_title','read_count','center_name','library_strategy','sample_accession']
    
    # split the list of sample IDs into chunks of max_chunk_size IDs each in order 
    # to avoid an exception due to URLs becoming too long
    sample_id_chunks = [sample_ids[i:i + chunk_size] for i in range(0, len(sample_ids), chunk_size)]
    
    ena_dfs=[]
    for sample_ids_chunk in sample_id_chunks:
        ena_query = _build_ena_sample_query(sample_ids_chunk)
        ena_response = _search_ena(query=ena_query, return_fields=fields)
        ena_dfs.append(_parse_ena_response(ena_response))
    ena_data = pd.concat(ena_dfs)
    
    # rename the sample_title column to match the original column name for the sample ID field
    ena_data.rename(columns={'sample_title': sample_id_col_name}, inplace=True)
    return ena_data


def create_ena_data_frame_from_run_accessions(df:pd.DataFrame, run_acc_col_name:str='run_accession', chunk_size:int=50):
    """
    Takes a DataFrame with a run_accession column (default name 'run_accession') and searches 
    ENA by run acc ('run_accession' field in ENA). Returns a new dataframe with run accs and 
    ftp download locations. 

    Args:
        df (pandas.DataFrame): A pandas DataFrame. Only needs a single column for run accs.  
        run_acc_col_name (str): Name of the run acc column, defaults to 'run_accession'
        chunk_size (int): number of samples queried in one request. If the total number of 
            run accs IDs is larger than this, the ENA query will be run in chunks.  
            Defaults to 50.  

    Returns:
        pandas.DataFrame: A new DataFrame that contains the input run accs IDs as well as the  
        following data fields from ENA:
        'fastq_ftp',

    """
    try:
        run_accs=df[run_acc_col_name].tolist()
    except KeyError:
        raise ValueError(
            f'DataFrame is missing run accession column, expected name {run_acc_col_name} ' +
            '(can be changed with param "run_acc_col_name")')
    
    fields = ['run_accession','fastq_ftp']
    
    # split the list of accs into chunks of max_chunk_size IDs each in order 
    # to avoid an exception due to URLs becoming too long
    run_acc_chunks = [run_accs[i:i + chunk_size] for i in range(0, len(run_accs), chunk_size)]
    
    ena_dfs=[]
    for run_acc_chunk in run_acc_chunks:
        ena_query = _build_ena_run_acc_query(run_acc_chunk)
        ena_response = _search_ena(query=ena_query, return_fields=fields)
        ena_dfs.append(_parse_ena_response(ena_response))
    ena_data = pd.concat(ena_dfs)
    
    return ena_data


def _search_ena(query:str, return_fields:list, limit:int=10000):
    """
    Run a search with the ENA API, querying by sample ID and return the search results response  
    in raw text form. 
    For details, consult the ENA API v2 documentation at 
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ

    Args:
        samples (list): a list of sample IDs to search with
        return_fields (list): a list of the field names to return from ENA
        limit (int): limit on number of results to return, defaults to 10000
        
    Returns:
        Raw text response from ENA
    """
    ENA_SEARCH_BASE_URL='https://www.ebi.ac.uk/ena/portal/api/search'
    search_params = {
        'result': 'read_run',
        'query': query,
        'fields': ','.join(return_fields),
        'limit': limit
    }
    request = requests.get(ENA_SEARCH_BASE_URL, params=search_params, timeout=15)
    request.raise_for_status() # throws exception if bad status returned
    return request.text

def _parse_ena_response(response_text:str):
    """
    Parse the response from a ENA search query into a pandas.DataFrame

    Args:
        response_text (str): response text from ENA API query
        
    Returns:
        pandas.DataFrame, see _search_ena for column names (ENA fieldnames)
    """
    return pd.read_csv(StringIO(response_text), sep="\t")
    
def _build_ena_sample_query(samples:list):
    """
    Create the query string for the ENA search request.  
    The format of the query string is detailed here:
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ

    Args:
        samples (list): list of sample IDs
    """
    return '(' + ' OR '.join([f'sample_title="{sample}"' for sample in samples ]) + ')'

def _build_ena_run_acc_query(run_accs:list):
    """
    Create the query string for the ENA search request.  
    The format of the query string is detailed here:
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ

    Args:
        run_accs (list): list of run accessions
    """
    return '(' + ' OR '.join([f'run_accession="{run_acc}"' for run_acc in run_accs ]) + ')'


def _build_ena_bioproject_query(bioproject:str):
    """
    Create the query string for the ENA search request.  
    The format of the query string is detailed here:
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ

    Args:
        samples (list): list of sample IDs
    """
    return "study_accession=" + bioproject


def align_ena_results_with_sample_data_genre_pf8(sample_data:pd.DataFrame, ena_result:pd.DataFrame, genre_panel_map:dict, sample_id_col_name:str='sample', skip_errors:bool=False):
    """
    Align ENA results, where one sample is expected to have more than one "run accession", 
    with sample data (which has one row per sample).  
    Data is assigned to Pf8 or GenRe based on the information in column 'center_name'.  
    A sanity check id performed against column 'library_strategy' (Pf8 data should be WGS, GenRe 
    should be AMPLICON). The GenRe panel is extracted from auxilliary dataframe that maps
    run accessions to panel (this information is not reliably found in the ENA records, so 
    must be provided externally). 
    
    NOTE
    Unlike other functions in this module, this function is only applicable to the GenRe/Pf8 
    malaria dataset, but it can be used as a starting point to create a merge function for other 
    data.  
    In this case, we expect every sample to have one "run accession" for each of the three 
    AmpSeq panels (GRC1, GRC2, SPEC) in GenRe plus one for Pf8 WGS.
    We create 4 pairs of new columns accordingly:  
        - INSDC_GenRe_GRC1, INSDC_GenRe_GRC1_readcount
        - INSDC_GenRe_GRC2, INSDC_GenRe_GRC2_readcount,
        - INSDC_GenRe_SPEC, INSDC_GenRe_SPEC_readcount,
        - INSDC_Pf8, INSDC_Pf8_readcount
    A version of the input ENA data frame is also returned that:
        - is filtered to only include retained samples; and 
        - has an additional library_name column, formed from the sample name and panel name (where 
          "WGS" is used as the panel name for Pf8)
        
    Args:
        ena_result (pandas.DataFrame): a DataFrame of ENA search results, expected to contain results for 
            4 run accessions per sample ID.  
            Must have columns:
                -sample (or alternative name given via sample_id_col_name parameter)
                -run_accession
                -center_name
                -library_strategy
        sample_data (pandas.DataFrame): a DataFrame of the original sample data, where each sample ID is unique 
            and into which the ENA data is to be merged
            sample_id_col_name (str): Name of the sample ID column, defaults to 'sample'
        skip_errors (bool): if True, failed dat sanity checks do not throw errors, data is just skipped
            defaults to False. 
  
    Returns:
        pandas.DataFrame: a new DataFrame based on sample_data with extra columns as shown above. 
        pandas.DataFrame: a new DataFrame based on ena_result, restricted to entries with sample in the first dataframe 
    """
    if not sample_id_col_name in sample_data or not sample_id_col_name in ena_result:
        raise ValueError('both DataFrames need to have a column "sample"')
    
    new_df = sample_data.copy()
    new_df['INSDC_Pf8'] = None
    new_df['INSDC_Pf8_readcount'] = None
    new_df['INSDC_GenRe_GRC1'] = None
    new_df['INSDC_GenRe_GRC1_readcount'] = None
    new_df['INSDC_GenRe_GRC2'] = None
    new_df['INSDC_GenRe_GRC2_readcount'] = None
    new_df['INSDC_GenRe_SPEC'] = None
    new_df['INSDC_GenRe_SPEC_readcount'] = None
  
    for _,row in ena_result.iterrows():
        try:
            sample = row[sample_id_col_name]
            run_accession = row['run_accession']
            center = row['center_name']
            strategy = row['library_strategy']
            readcount = row['read_count']
        except KeyError as e:
            raise ValueError(f'ENA result DataFrame is missing required column: {e}')
        
        run_accession_type = None
        if 'Wellcome Sanger' in center:
            if strategy!='WGS':
                msg=f'This row of ENA search results is expected to contain a Pf8 WGS sample but the library strategy is not "WGS": {row}'
                if skip_errors:
                    print(msg + ' - skip_errors active, skipping this row')
                else:
                    raise ValueError(msg)
            run_accession_type='Pf8'
            
        elif 'GenRe-Mekong' in center:
            if strategy!='AMPLICON':
                msg=f'This row of ENA search results is expected to contain a GenRe AMPLICON sample but the library strategy is not "AMPLICON": {row}'
                if skip_errors:
                    print(msg + ' - skip_errors active, skipping this row')
                else:
                    raise ValueError(msg)
                
            if genre_panel_map[run_accession]:
                run_accession_type = "GenRe_" + genre_panel_map[run_accession]
            else:
                msg=f'could not extract primer panel from GenRe data row: {row}'
                if skip_errors:
                    print(msg + ' - skip_errors active, skipping this row')
                else:
                    raise ValueError(msg)
                
        if not run_accession_type:
            msg=f'could not assign ENA result to Pf8 or GenRe {row}'
            if skip_errors:
                print(msg + ' - skip_errors active, skipping this row')
            else:
                raise ValueError(msg)
        
        accession_col = 'INSDC_'+run_accession_type
        readcount_col = accession_col+'_readcount'
        if new_df.loc[new_df['sample'] == sample, accession_col].any():
            raise ValueError(f'More than one run accessions found for sample {sample}, field {accession_col}')
        else:
            new_df.loc[new_df['sample'] == sample, accession_col]=run_accession
            # ENA read counts are summed over both ends in paired end. So need to divide by two to get 
            # number of read pairs (which is almost always what we want)
            new_df.loc[new_df['sample'] == sample, readcount_col]=readcount / 2
 
    new_df = new_df[
        new_df[['INSDC_GenRe_GRC1','INSDC_GenRe_GRC2','INSDC_GenRe_SPEC','INSDC_Pf8']].notnull().all(1)
    ]
    new_ena_df = ena_result[ena_result["sample"].isin(new_df["sample"])]
    new_ena_df["library_name"] = new_ena_df["sample"] + "_" + new_ena_df["run_accession"].map(lambda x: genre_panel_map.get(x, "WGS"))

    return new_df,new_ena_df


def download_all_fastqs(outdir, data:pd.DataFrame=None, data_file_path:str=None, num_tries:int=3,skip_errors:bool=False,top3:bool=False):
    """
    Download the FASTQ files from a table of FTP URLs, which can be provided as a DataFrame or a path to a 
    csv file. The file must contains a column for the ENA run accession and one column each for the FASTQ 
    files for read1 and read2.  
    The run accession column is not used for the download, it is just passed through to the manifest file 
    that is created as a list of local file paths for the respective FASTQ files.  
    
    The output directory will be created and must not exist already.  

    Args:
        data (pandas.DataFrame): DataFrame with run_accessions 
        data_file_path (str): path to a csv file with the run_accessions 
        outdir (str): path to output dir (must not exist yet)
        skip_errors (bool): if True, skip download errors and continue with next time, don't throw exception. 
        top3 (bool): if True, only process the first 3 rows of the data (useful for testing)
        
    Returns:
        True on success
        
        If option create_manifest in use, creates a CSV file in outdir with the following fields:
        'run_accession','ftp_url_read_1','ftp_url_read_2','read_1_file','read_2_file'

    """
    if not isinstance(data, pd.DataFrame) and not data_file_path:
        raise ValueError('must provide either "data" or "data_file_path" parameter')
    if isinstance(data, pd.DataFrame) and data_file_path:
        raise ValueError('must provide either "data" or "data_file_path" parameter, not both')
    if data_file_path:
        data = pd.read_csv(data_file_path)

    ena_locs_df = create_ena_data_frame_from_run_accessions( data )
 
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
       
    n_rows = len(ena_locs_df)
    i = 0
    for _, row in ena_locs_df.iterrows():
        i+=1
        if top3 and i>3:
            print("option top3 in use: stopping after three rows processed")
            break

        fastq_ftp = row["fastq_ftp"]
        if not ';' in fastq_ftp:
            raise ValueError(f'FASTQ FTP field in this row does not contain 2 links: {row}')
    
        fq_url_1, fq_url_2 = fastq_ftp.split(';')
        fq_url_1.strip()
        if not fq_url_1.startswith('ftp://'):
                fq_url_1 = 'ftp://' + fq_url_1
        fq_url_2.strip()
        if not fq_url_2.startswith('ftp://'):
                fq_url_2 = 'ftp://' + fq_url_2

        _download_fastq_file( fq_url_1, outdir, num_tries=num_tries, skip_errors=skip_errors)
        _download_fastq_file( fq_url_2, outdir, num_tries=num_tries, skip_errors=skip_errors)
        print(f'downloaded FASTQ {i} of {n_rows}')
               
    return True

def _download_fastq_file( remote_ftp_url:str, dir, num_tries:int=3, skip_errors:bool=False):
    """
    Downloads a single FASTQ file from a remote URL path into a local file

    Args:
        remote_ftp_url (str): URL of the remote file
        dir (PosixPath): path to directory into which the file is downloaded
        num_tries (int): number of times download should be tried in case of errors
        skip_errors (bool): if True, skip download errors and continue with next time, don't throw exception. 

    Returns:
        PosixPath of locally downloaded file
    """
    if not remote_ftp_url.startswith('ftp://'):
        remote_ftp_url = 'ftp://'+remote_ftp_url
    
    file_name = remote_ftp_url.split('/')[-1]
    local_path = pathlib.Path(dir) / file_name
    local_path = local_path.resolve()
    
    attempt = 1
    last_error = None
    while not local_path.exists() and attempt <= num_tries:
        try:
            urllib.request.urlretrieve(remote_ftp_url, local_path)
        except Exception as e:
            print(f'download attempt {attempt} of {num_tries} failed for URL {remote_ftp_url}')
            last_error = e
            attempt += 1
            time.sleep(5)
            
    if not local_path.exists():
        msg = f'Failed to download {remote_ftp_url} after {attempt} attempts. Last error raised: {last_error}'
        if skip_errors:
            print("\n".join([msg,'   Option skip_errors in use, continue with next download']))
            local_path='DOWNLOAD FAILED'
        else:
            raise Exception(msg)
    else:
        # give some feedback if this took more than one attempt
        if attempt>1:
            print(f'download attempt {attempt} of {num_tries} succeeded for URL {remote_ftp_url}')
    
    return str(local_path)

def cli_download_fastqs(args):
    """
    CLI function to run download_all_fastqs
    """
    download_all_fastqs(
        outdir=args.output_folder,
        data_file_path=args.insdc_manifest,
        num_tries=args.download_attempts,
        skip_errors=args.skip_errors,
        top3=args.top3
    )
    return 0

def main():
    
    args = args_parser().parse_args()
    return args.func(args) 

if __name__ == "__main__":
    exit(main())
