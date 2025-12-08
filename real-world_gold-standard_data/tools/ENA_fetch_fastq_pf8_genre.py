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

import requests
import argparse
import csv
import pathlib
import urllib.request

ENA_SEARCH_BASE_URL='https://www.ebi.ac.uk/ena/portal/api/search'

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

def parse_sample_ids(data_file:str):
    """
    Parse sample IDs from data file

    Args:
        data_file (str): path to a csv data file
    """
    samples = []
    with open(data_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        if not reader.fieldnames or not 'sample' in reader.fieldnames:
            raise ValueError(f"data file {data_file} is missing required header 'sample' (all lower case)")
        for row in reader:
            samples.append(row['sample'])
    return samples

def search_ena(samples:list):
    """
    Run a search with the ENA API, querying by sample ID and return the search results

    Args:
        samples (list): a list of sample IDs to search with
        
    Returns:
        search result data as a list of dicts
    """
    search_params = {
        'result': 'read_run',
        'query': _build_ena_query(samples),
        'fields': ','.join(['sample_title','center_name','library_strategy','sample_accession','fastq_ftp','submitted_ftp','run_accession']),
        'limit': 1000
    }
    request = requests.get(ENA_SEARCH_BASE_URL, params=search_params, timeout=10)
    request.raise_for_status() # throws exception if bad status returned
    return _parse_ena_response( request.text)

def _parse_ena_response(response_text:str):
    """
    Parse the response from a ENA search query into a list of dicts.  
    

    Args:
        response_text (str): response text from ENA API query
        
    Returns:
        list of dicts:
            [ 
                {'sample_title': TITLE, 'center_name': CENTER, 'fastq_ftp': FASTQ,'submitted_ftp': SUBMITTED},{...}
            ]
    """
    result = []
    response_lines = response_text.splitlines()
    # should receive tab delimited text - this would fail if the API should change output format
    reader = csv.DictReader(response_lines, delimiter="\t")
    for row in reader:
        result.append(row)
    return result
    
def _build_ena_query(samples:list):
    """
    Create the query string for the ENA search request

    Args:
        samples (list): list of sample IDs
    """
    return '(' + ' OR '.join([f'sample_title="{sample}"' for sample in samples ]) + ')'

def split_search_result_by_resource(ena_data:list, skip_errors:bool=True):
    """
    Splits the ENA search result into a Pf8 and a GenRe Mekong set.  
    The split is done on content of the field 'center_name' and a sanity 
    check is performed to ensure that the 'library_strategy' matches the center, 
    i.e. 'WGS' for 'Wellcome Sanger' (Pf8) and 'AMPLICON' for 'GenRe-Mekong'
    
    For GenRe data, we also add a field 'panel', which is derived from the name of the 
    submitted cram file and reflects the primer panel 'GRC1', 'GRC2' or 'SPEC'
    
    NOTE: This is a convenience method that is very specific to the Pf8/GenRe data, where a 
    single sample ID points to data from two resources. For a more generic version of this 
    tool, this method needs to be removed or made optional and any filtering of datasets needs 
    to be done by the user before FASTQ data can be retrieved.  
    
    Args:
        ena_data (list): ENA search result list of dicts (from search_ena())
        skip_errors (bool): If True, don't throw exceptions, just skip the offending rows
        
    Returns:
        Two lists of dicts in the same format as the input:
            pf8_data, genre_data
    """
    pf8_data = []
    genre_data = []
    for row in ena_data:
        try:
            center = row['center_name']
            strategy = row['library_strategy']
            submitted_ftp = row['submitted_ftp']
        except KeyError as e:
            raise ValueError(f'An expected header was missing from the search result data: {e}')
        if 'Wellcome Sanger' in center:
            if strategy!='WGS':
                msg=f'This row of ENA search results is expected to contain a Pf8 WGS sample but the library strategy is not "WGS": {row}'
                if skip_errors:
                    print(msg + ' - skip_errors active, skipping this row')
                else:
                    raise ValueError(msg)
            pf8_data.append(row)
        elif 'GenRe-Mekong' in center:
            if strategy!='AMPLICON':
                msg=f'This row of ENA search results is expected to contain a GenRe AMPLICON sample but the library strategy is not "AMPLICON": {row}'
                if skip_errors:
                    print(msg + ' - skip_errors active, skipping this row')
                else:
                    raise ValueError(msg)
            if 'GRC1' in submitted_ftp:
                row['panel'] = 'GRC1'
            elif 'GRC2' in submitted_ftp:
                row['panel'] = 'GRC2'
            if 'SPEC' in submitted_ftp:
                row['panel'] = 'SPEC'
            else:
                msg=f'could not extract primer panel from GenRe data row: {row}'
                if skip_errors:
                    print(msg + ' - skip_errors active, skipping this row')
                else:
                    raise ValueError(msg)
            genre_data.append(row)
    return pf8_data, genre_data

def download_all_fastqs(name, data, outdir, create_manifest:bool=True, num_tries:int=3):
    """
    Download the FASTQ files from ENA for a list of dicts of data. FASTQ files are downloaded 
    into a new sub-directory called {name}_fastq.  
    If create_manifest is used, a manifest file is created in {outdir}, which contains sample IDs and 
    read1/read2 fastq file paths to the download directory.  

    Args:
        name (str): Name of the dataset (Pf8 or GenRe_Mekong)
        data (list(dict)): sample data
        outdir (str): path to output dir
        create_manifest (bool): if True, a manifest file is created in {outdir}
        
    Returns:
        True on success
        
        If option create_manifest in use, creates a CSV file in outdir with the following fields:
        'sample','run_accession','center_name','library_strategy','panel','read_1_file','read_2_file'
        
        All field values are taken from the ENA search results, except for the two file paths, which point 
        to the locally downloaded FASTQ files.  
    """
    # create a subdirectory for downloading the data
    fastq_dir = pathlib.Path(outdir) / (name + '_fastq')
    fastq_dir.mkdir(parents=True, exist_ok=False)
    
    if create_manifest:
        manifest_path = pathlib.Path(outdir) / 'manifest.csv'
        fields = ['sample','run_accession','center_name','library_strategy','panel','read_1_file','read_2_file']
        manifest_fh = open(manifest_path, 'w', newline='')
        manifest_writer = csv.DictWriter(manifest_fh, fieldnames=fields)
        manifest_writer.writeheader()
    
    for row in data:
        try:
            fastq_ftps = row['fastq_ftp']
            sample = row['sample']
            run_accession = row['run_accession']
            center_name = row['center_name']
            library_strategy = row['library_strategy']
            if 'panel' in row:
                panel = row['panel']
            else:
                panel = 'N/A'
        except KeyError as e:
            raise ValueError(f'search result is missing required fields: {e} - this should not happen and indicates a bug in the script')
        
        local_paths = []
        for fastq_ftp in fastq_ftps.split(';'):
            fastq_ftp.strip()
            local_path = _download_fastq_file( fastq_ftp, fastq_dir, num_tries=num_tries)
            local_paths.append(str(local_path.resolve()))
            
        if create_manifest:
            manifest_writer.writerow({
                'sample': sample,
                'run_accession': run_accession,
                'center_name': center_name,
                'library_strategy': library_strategy,
                'panel': panel,
                'read_1_file': local_paths[0],
                'read_2_file': local_paths[1]
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
    local_path = dir / file_name
    
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
    
    return local_path

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