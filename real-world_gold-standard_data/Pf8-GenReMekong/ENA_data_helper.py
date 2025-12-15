"""
This script downloads FASTQ files for the Pf8/GenRe Mekong datasets provided in this resource.  

To run the script, make sure all dependencies are installed. A file with dependencies is provide in this folder. To install dependencies, run the 
following command:
pip install requirements.txt 

The script will download FASTQ files for ALL samples in the datafile by default. 
Before going ahead with the full download, it is recommended to run the script with 
the --top3 option first. This will run the downloader on the top 3 rows of the data 
file only and provide an opportunity to check the output before downloading all 
data. 

The script creates an output directory into which it writes two manifest files, one 
for Pf8 and one for GenRe Mekong data. The manifest files contain sample and accession IDs 
as well as paths to the FASTQ files. For GenRe data, each sample has three rows, one for each 
of the three primer panels GRC1, GRC2 and SPEC. For Pf8, each sample only contains a single row 
of data.  
One sub-directory is created for Pf8 and one for GenRe, which contains the downloaded FASTQ files 
that are references in the manifest files.  

Run the script with:
python ENA_featch_fastq_pf8_genre.py --data DATAFILE

For a list of options run:
python ENA_featch_fastq_pf8_genre.py -h

"""
import pandas as pd
import requests
import argparse
import csv
import pathlib
from io import StringIO
import urllib.request

def args_parser():     
    parser = argparse.ArgumentParser(
        description = "ENA_fetch_fastq_pf8_genre.py: helper tool for downloading FASTQ from ENA") 

    parser.add_argument(
        '--data','-d',
        action = 'store',
        required = True,
        metavar = 'FILE', 
        type = str,
        help='path to the dataset csv file. The file must contain a column "sample"')

    parser.add_argument(
        '--out','-o',
        action = 'store',
        required = True,
        metavar = 'DIR', 
        type = str,
        help='path to a directory where results will be written into subdirectories. Must not exist yet')

    parser.add_argument(
        '--top3',
        action = 'store_true',
        help = 'when used, only the data for the first 3 rows in the data file are downloaded. Use this first to ensure that the output is what you are expecting before proceeding to a full download' )

    parser.add_argument(
        '--skip_errors',
        action = 'store_true',
        help ='when used, rows of data that would otherwise throw and exception are just skipped and a message is printed'
    )

    parser.add_argument(
        '--download_attempts',
        action = 'store',
        default= 3,
        required = False, 
        type = int,
        help='Optional number of attempts to try FTP downloads. Defaults to 3')


    return parser


def create_ena_data_frame(df:pd.DataFrame, sample_id_col_name:str='sample', chunk_size:int=50):
    """
    Takes a DataFrame with a sample ID column (default name 'sample') and searches 
    ENA by sample ID ('sample_title' field in ENA). Returns a new dataframe with sample ID 
    and ENA data. A single sample may have multiple ENA records, hence sample ID is not 
    unique in this DataFrame.  

    Args:
        df (pandas.DataFrame): A pandas DataFrame. Only needs a single column for sample IDs.  
            By default, this is expected to eb called 'sample' but can be changed. Other data  
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
    
    fields = ['sample_title','run_accession','center_name','library_strategy','sample_accession','fastq_ftp','submitted_ftp']
    
    # split the list of sample IDs into chunks of max_chunk_size IDs each in order 
    # to avoid an exception due to URLs becoming too long
    sample_id_chunks = [sample_ids[i:i + chunk_size] for i in range(0, len(sample_ids), chunk_size)]
    
    ena_dfs=[]
    for sample_ids_chunk in sample_id_chunks:
        ena_response = _search_ena(samples=sample_ids_chunk, return_fields=fields)
        ena_dfs.append(_parse_ena_response(ena_response))
    ena_data = pd.concat(ena_dfs)
    
    # rename the sample_title column to match the original column name for the sample ID field
    ena_data.rename(columns={'sample_title': sample_id_col_name}, inplace=True)
    return ena_data

def _search_ena(samples:list, return_fields:list, limit:int=10000):
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
        'query': _build_ena_query(samples),
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
    
def _build_ena_query(samples:list):
    """
    Create the query string for the ENA search request.  
    The format of the query string is detailed here:
    https://docs.google.com/document/d/1CwoY84MuZ3SdKYocqssumghBF88PWxUZ

    Args:
        samples (list): list of sample IDs
    """
    return '(' + ' OR '.join([f'sample_title="{sample}"' for sample in samples ]) + ')'


def merge_ena_results_into_sample_data_genre_pf8(sample_data:pd.DataFrame, ena_result:pd.DataFrame,  sample_id_col_name:str='sample', skip_errors:bool=False, include_download_link:bool=True):
    """
    Merges ENA results, where one sample is expected to have more than one "run accession", 
    into a DataFrame where sample ID is unique (index).  
    Data is assigned to Pf8 or GenRe based on the information in column 'center_name'.  
    A sanity check id performed against column 'library_strategy' (Pf8 data should be WGS, GenRe 
    should be AMPLICON). The GenRe panel is extracted from the submitted cram file names in 
    'submitted_ftp'. If download links are to be included, two mate FASTQ files are expected in 
    field 'fastq_ftp'.  
    
    NOTE
    Unlike other functions in this module, this function is only applicable to the GenRe/Pf8 
    malaria dataset, but it can be used as a starting point to create a merge function for other 
    data.  
    In this case, we expect every sample to have one "run accession" for each of the three 
    AmpSeq panels (GRC1, GRC2, SPEC) in GenRe plus one for Pf8 WGS.
    We create 4 new columns accordingly:  
        - ENA_acc_GenRe_GRC1
        - ENA_acc_GenRe_GRC2
        - ENA_acc_GenRe_SPEC
        - ENA_acc_Pf8
        
    If include_download_link is used (default), the following 8 additional columns are added, 
    which contain the FTP download links for the FASTQ files (mate 1 and 2) corresponding to the ENA 
    run accession IDs:
        - GenRe_GRC1_ENA_FASTQ_FTP_1 / GenRe_GRC1_ENA_FASTQ_FTP_2
        - GenRe_GRC2_ENA_FASTQ_FTP_1 / GenRe_GRC2_ENA_FASTQ_FTP_2
        - GenRe_SPEC_ENA_FASTQ_FTP_1 / GenRe_SPEC_ENA_FASTQ_FTP_2
        - Pf8_ENA_FASTQ_FTP_1 / Pf8_ENA_FASTQ_FTP_2

    Args:
        ena_result (pandas.DataFrame): a DataFrame of ENA search results, expected to contain results for 
            4 run accessions per sample ID.  
            Must have columns:
                -sample (or alternative name given via sample_id_col_name parameter)
                -run_accession
                -center_name
                -library_strategy
                -submitted_ftp
                -fastq_ftp (if include_download_link in use)
        sample_data (pandas.DataFrame): a DataFrame of the original sample data, where each sample ID is unique 
            and into which the ENA data is to be merged
            sample_id_col_name (str): Name of the sample ID column, defaults to 'sample'
        skip_errors (bool): if True, failed dat sanity checks do not throw errors, data is just skipped
            defaults to False. 
        include_download_link (bool): if True (default), adds a column with the FTP download link for every 
            run accession column

    Returns:
        pandas.DataFrame): a new DataFrame based on sample_data with extra columns as shown above. 
    
    """
    if not sample_id_col_name in sample_data or not sample_id_col_name in ena_result:
        raise ValueError('both DataFrames need to have a column "sample"')
    
    new_df = sample_data.copy().set_index('sample')
    new_df['ENA_acc_Pf8'] = None
    new_df['ENA_acc_GenRe_GRC1'] = None
    new_df['ENA_acc_GenRe_GRC2'] = None
    new_df['ENA_acc_GenRe_SPEC'] = None

    if include_download_link:
        for i in [1,2]:
            new_df['GenRe_GRC1_ENA_FASTQ_FTP_'+str(i)] = None
            new_df['GenRe_GRC2_ENA_FASTQ_FTP_'+str(i)] = None
            new_df['GenRe_SPEC_ENA_FASTQ_FTP_'+str(i)] = None
            new_df['Pf8_ENA_FASTQ_FTP_'+str(i)] = None

    for _,row in ena_result.iterrows():
        try:
            sample = row[sample_id_col_name]
            run_accession = row['run_accession']
            center = row['center_name']
            strategy = row['library_strategy']
            submitted_ftp = row['submitted_ftp']
            if include_download_link:
                fastq_ftp = row['fastq_ftp']
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
            if 'GRC1' in submitted_ftp:
                run_accession_type='GenRe_GRC1'
            elif 'GRC2' in submitted_ftp:
                run_accession_type='GenRe_GRC2'
            elif 'SPEC' in submitted_ftp:
                run_accession_type= 'GenRe_SPEC'
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
        
        accession_col = 'ENA_acc_'+run_accession_type
        if new_df.loc[[sample], accession_col].any():
            raise ValueError(f'More than one run accessions found for sample {sample}, field {accession_col}')
        else:
            new_df.loc[sample, accession_col]=run_accession
            if include_download_link:
                if not ';' in fastq_ftp:
                    raise ValueError(f'FASTQ FTP field in this row does not contain 2 links: {row}')
                for i, ftp_url in enumerate(fastq_ftp.split(';')):
                    ftp_url.strip()
                    if not ftp_url.startswith('ftp://'):
                        ftp_url = 'ftp://'+ftp_url
                    ftp_col = run_accession_type+'_ENA_FASTQ_FTP_'+str(i+1)
                    new_df.loc[sample, ftp_col]=ftp_url

    return new_df


def download_all_fastqs(outdir, data:pd.DataFrame=None, data_file_path:str=None, create_manifest:bool=True, num_tries:int=3, ftp_url_read_1_col:str='ftp_url_read_1',ftp_url_read_2_col:str='ftp_url_read_2',run_accession_col:str='run_accession'):
    """
    Download the FASTQ files from a table of FTP URLs, which can be provided as a DataFrame or a path to a 
    csv file. The file must contains a column for the ENA run accession and one column each for the FASTQ 
    files for read1 and read2.  
    The run accession column is not used for the download, it is just passed through to the manifest file 
    that is created as a list of local file paths for the respective FASTQ files.  
    
    The output directory will be created and must not exist already.  

    Args:
        data (pandas.DataFrame): DataFrame with run_accession and read1/2 FTP URLs to download from
        data_file_path (str): path to a csv file with the run_accession and FTP URLs 
        
        ftp_url_read_1_col (str): name of the column in {data} that contains the read 1 FTP URLs,  
            defaults to 'ftp_url_read_1'
        ftp_url_read_2_col (str): name of the column in {data} that contains the read 2 FTP URLs,  
            defaults to 'ftp_url_read_2'
        run_accession_col (str): name of the column in {data} that contains run accession IDs, 
            defaults to 'run_accession'
        outdir (str): path to output dir (must not exist yet)
        create_manifest (bool): if True, a manifest file is created in {outdir}
        
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
    
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    if create_manifest:
        manifest_path = outdir / 'manifest.csv'
        fields = 'run_accession','ftp_url_read_1','ftp_url_read_2','read_1_file','read_2_file'
        manifest_fh = open(manifest_path, 'w', newline='')
        manifest_writer = csv.DictWriter(manifest_fh, fieldnames=fields)
        manifest_writer.writeheader()
    
    expected_cols = [run_accession_col, ftp_url_read_1_col, ftp_url_read_2_col]
    if not set(expected_cols).issubset(set(data.columns.values)):
        raise ValueError(f'ENA data table does not have all of the expected column headers: {expected_cols}')

    n_rows = len(data)
    i = 0
    for _, row in data.iterrows():
        i+=1
        run_accession = row[run_accession_col]
        ftp_url_read_1 = row[ftp_url_read_1_col]
        ftp_url_read_2 = row[ftp_url_read_2_col]
        read_1_file = _download_fastq_file( ftp_url_read_1, outdir, num_tries=num_tries)
        read_2_file = _download_fastq_file( ftp_url_read_2, outdir, num_tries=num_tries)
        print(f'downloaded FASTQ pair {i} of {n_rows}')
            
        if create_manifest:
            manifest_writer.writerow({
                'run_accession': run_accession,
                'ftp_url_read_1': ftp_url_read_1,
                'ftp_url_read_2': ftp_url_read_2,
                'read_1_file': read_1_file,
                'read_2_file': read_2_file
            })
        
    if create_manifest:
        manifest_fh.close()
    
    return True

def _download_fastq_file( remote_ftp_url:str, dir, num_tries:int=3):
    """
    Downloads a single FASTQ file from a remote URL path into a local file

    Args:
        remote_ftp_url (str): URL of the remote file
        dir (PosixPath): path to directory into which the file is downloaded
        num_tries (int): number of times download should be tried in case of errors

    Returns:
        PosixPath of locally downloaded file
    """
    if not remote_ftp_url.startswith('ftp://'):
        remote_ftp_url = 'ftp://'+remote_ftp_url
    
    file_name = remote_ftp_url.split('/')[-1]
    local_path = pathlib.Path(dir) / file_name
    
    attempt = 1
    last_error = None
    while not local_path.exists() and attempt <= num_tries:
        try:
            urllib.request.urlretrieve(remote_ftp_url, local_path)
        except Exception as e:
            print(f'download attempt {attempt} of {num_tries} failed for URL {remote_ftp_url}')
            last_error = e
            attempt += 1
            
    if not local_path.exists():
        raise Exception(f'Failed to download {remote_ftp_url} after {attempt} attempts. Last error raised: {last_error}')
    
    return str(local_path)

def main():
    
    args = args_parser().parse_args()
    
    outdir = args.out
    pathlib.Path(outdir).mkdir(parents=True, exist_ok=False)
    samples = parse_sample_ids(args.data)
    ena_data = search_ena(samples)
    pf8_data, genre_data = split_search_result_by_resource(ena_data, skip_errors=args.skip_errors)
    
    download_attempts = args.download_attempts
    for name, data in dict(zip(['Pf8','GenRe_Mekong'],[pf8_data, genre_data])).items():
        download_fastq(name, data, outdir, create_manifest=True, num_tries=download_attempts)
    
if __name__ == "__main__":
    exit(main())