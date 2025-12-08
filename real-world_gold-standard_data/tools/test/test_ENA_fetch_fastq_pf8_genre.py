import pytest
from pathlib import Path
from urllib.error import URLError
from ..ENA_fetch_fastq_pf8_genre import parse_sample_ids, _build_ena_query, search_ena, _parse_ena_response, split_search_result_by_resource, _download_fastq_file, download_all_fastqs

TEST_FILES_DIR = Path(__file__).parent.resolve() / 'test_data'

def test_parse_sample_ids():
    test_data = TEST_FILES_DIR / 'Pf8_genre_test_data.csv'
    samples = parse_sample_ids(test_data)
    assert samples, 'could parse samples from file'
    assert isinstance(samples, list)
    assert len(samples) == 3, 'there are 3 rows in the test data'
    assert samples[0] == 'RCN13568', 'first sample ID from file is correct'

def test__build_ena_query():
    query = _build_ena_query(['abc123','def3456'])
    assert query=='(sample_title="abc123" OR sample_title="def3456")'
    
def test__parse_ena_response():
    pass
    
def test_search_ena():
    """
    This test will run a query against the ENA search API using three
    real-world Pf8/GenRe sample IDs (="sample_title" in ENA).  
    If this test fails, make sure it is not an issue with your network connections.  
    """
    # use three existing sample IDs from Pf8 (sample titles in ENA)
    samples = ['RCN13568','RCN13560','RCN15107']
    result = search_ena(samples)
    assert result, 'should receive some results'
    assert isinstance(result, list)
    assert isinstance(result[0], dict), 'result is a list of dicts'
    assert 'sample_title' in result[0], 'result has header "sample_title"'
    assert 'center_name' in result[0], 'result has header "center_name"'
    assert 'fastq_ftp' in result[0], 'result has header "fastq_ftp"'
    assert 'submitted_ftp' in result[0], 'result has header "submitted_ftp"'
    assert 'library_strategy' in result[0], 'result has header "library_strategy"'
    assert 'run_accession' in result[0],  'result has header "run_accession"'
    
    assert ';' in result[0]['fastq_ftp'], 'there should be two file paths, separated by ;'
    assert '_1.fastq' in result[0]['fastq_ftp'], 'there should be a path to the read 1 file'
    assert '_2.fastq' in result[0]['fastq_ftp'], 'there should be a path to the read 2 file'

def test_split_search_result_by_resource():
    search_data = [
        {
            'sample_title': 'RCN13560', # = sample ID
            'center_name': 'Wellcome Sanger Institute;WSI', 
            'library_strategy': 'WGS', 
            'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15626087/RCN13560.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15626087/RCN13560.cram.crai'
        }, 
        {
            'sample_title': 'RCN15107', 
            'center_name': 'The GenRe-Mekong Project;GenRe-Mekong', 
            'library_strategy': 'AMPLICON',
            'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_SPEC_29632.cram;ftp.sra.ebi.ac.../ERR14392568/RCN15107_SPEC_29632.cram.crai'
        }
    ]
    pf8_data, genre_data = split_search_result_by_resource(search_data)
    
    assert len(pf8_data) == 1, 'one of the rows is PF8/WSI data'
    assert len(genre_data) == 1, 'one of the rows is GenRe data'
    
    assert pf8_data[0]['sample_title'] == 'RCN13560','Pf8 row is the correct sample'
    assert genre_data[0]['sample_title'] == 'RCN15107','Pf8 row is the correct sample'

    assert 'panel' not in pf8_data[0], 'Pf8 data should not have an additional field for "panel"'
    assert 'panel' in genre_data[0], 'GenRe data has an additional field for primer "panel"'
    assert genre_data[0]['panel']=='SPEC', 'panel name was correctly extracted from submitted cram file name'
    
    # rows from WSI that are not "WGS" type and rows from GenRe that are not AMPLICON should raise exceptions
    search_data=[{
        'sample_title': '1234567', # = sample ID
        'center_name': 'Wellcome Sanger Institute;WSI', 
        'library_strategy': 'some unexpected strategy', 
        'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15626087/RCN13560.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR156/ERR15626087/RCN13560.cram.crai'
    }]
    
    # Pf8 (WSI) data should be WGS library_strategy
    with pytest.raises(ValueError, match=r'Pf8 WGS .+ library strategy'):
        pf8_data, genre_data = split_search_result_by_resource(search_data, skip_errors= False)
        
    search_data=[{
        'sample_title': '1234567', # = sample ID
        'center_name': 'The GenRe-Mekong Project;GenRe-Mekong', 
        'library_strategy': 'some unexpected strategy', 
        'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_SPEC_29632.cram;ftp.sra.ebi.ac.../ERR14392568/RCN15107_SPEC_29632.cram.crai'
    }]
    
    # GenRe data should be AMPLICON library_strategy
    with pytest.raises(ValueError, match=r'GenRe AMPLICON .+ library strategy'):
        pf8_data, genre_data = split_search_result_by_resource(search_data, skip_errors= False)
        
    # test GenRe where the submitted cram file doesn't contain the primer panel name SPEC, GRC1 or GRC2
    search_data=[{
        'sample_title': 'RCN15107', 
        'center_name': 'The GenRe-Mekong Project;GenRe-Mekong', 
        'library_strategy': 'AMPLICON',
        'submitted_ftp': 'ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_XXXXX_29632.cram;ftp.sra.ebi.ac.../ERR14392568/RCN15107_XXXXX_29632.cram.crai'
    }]
    
    # GenRe data must have panel name in submitted_ftp cram file name
    with pytest.raises(ValueError, match=r'could not extract primer panel'):
        pf8_data, genre_data = split_search_result_by_resource(search_data, skip_errors= False)
        
    # with "skip_errors" (default), no error is raised
    assert split_search_result_by_resource(search_data, skip_errors= True)

def test__download_fastq_file(tmp_path):
    # NOTE: this is a real SRA FTP path and relies on a file existing on the 
    # third-party resource. The test may fail because the file was removed remotely
    remote_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/003/ERR14388603/ERR14388603_1.fastq.gz'
    local_path = _download_fastq_file( remote_url, tmp_path)
    assert local_path, 'method returns the local file path'
    assert local_path.exists(), 'the file has been created'
    
    # Force an error by providing a non-existing URL
    remote_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/003/ERR14388603/THIS_DOES_NOT_EXIST.fastq.gz'
    with pytest.raises(Exception, match='Failed to download'):
        _download_fastq_file( remote_url, tmp_path)

def test_download_all_fastqs(tmp_path):
    data = [
        {
            'run_accession': 'ERR15626087',
            'sample': 'RCN15107', 
            'center_name': 'The GenRe-Mekong Project;GenRe-Mekong', 
            'library_strategy': 'AMPLICON', 
            'fastq_ftp': 'ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/068/ERR14392568/ERR14392568_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/068/ERR14392568/ERR14392568_2.fastq.gz',
            'panel': 'SPEC'}
    ]
    name = 'GenRe'
    download_all_fastqs(name= name, data=data,outdir=tmp_path)
    
    expected_manifest_path = tmp_path / 'manifest.csv'
    assert expected_manifest_path.exists(), 'a manifest file has been created'
    
    expected_fastq_dir_path = tmp_path /  (name + '_fastq')
    assert expected_fastq_dir_path.exists(), 'a dir has been created for the FASTQ files'
    
    
