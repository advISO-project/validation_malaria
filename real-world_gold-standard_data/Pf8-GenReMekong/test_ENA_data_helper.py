import pytest
from urllib.error import URLError
import gzip
import pathlib
import pandas as pd
from ENA_data_helper import create_ena_data_frame, _build_ena_query, _search_ena, _parse_ena_response, merge_ena_results_into_sample_data_genre_pf8, _download_fastq_file, download_all_fastqs

def test_create_ena_data_frame():
    data = pd.DataFrame(
        {'sample':['RCN13568','RCN13560','RCN15107'],'some_other_data': [1,2,3] }
    )
    ena_data = create_ena_data_frame(data, sample_id_col_name='sample')
    assert isinstance(ena_data, pd.DataFrame)
    assert len(ena_data)==12, 'the sample IDs are GenRe and Pf8 so should have 3 GenRe and 1 Pf8 result each = 12 in total'

    expected_col_names = ['sample','run_accession','center_name','library_strategy','sample_accession','fastq_ftp','submitted_ftp']
    assert sorted(list(ena_data.columns.values)) == sorted(expected_col_names), 'DataFrame has correct column (field) names'
    
    # assert that the fastq_ftp field contains a list of two mate files
    fastq_ftp_result1 = ena_data['fastq_ftp'][0]
    assert ';' in fastq_ftp_result1, 'there should be two file paths, separated by ;'
    assert '_1.fastq' in fastq_ftp_result1, 'there should be a path to the read 1 file'
    assert '_2.fastq' in fastq_ftp_result1, 'there should be a path to the read 2 file'
    
def test_create_ena_data_frame_large():
    """
    Special case: test with a large number of sample IDs, which would throw a 
    TTPError: 414 Client Error: Request-URI Too Large for url
    This tests the query in chunks behaviour
    """
    # a live search is performed, so these need to be real ENA sample IDs
    # There are 60 IDs, which should trigger a search in 2 chunks
    sample_ids = ['RCN12025', 'RCN12026', 'RCN12028', 'RCN12031', 'RCN12032',
        'RCN12033', 'RCN12035', 'RCN12036', 'RCN12040', 'RCN12041',
        'RCN12043', 'RCN12044', 'RCN12051', 'RCN12058', 'RCN12059',
        'RCN12060', 'RCN12070', 'RCN12071', 'RCN12072', 'RCN12075',
        'RCN12076', 'RCN12079', 'RCN12080', 'RCN12082', 'RCN12083',
        'RCN12084', 'RCN12100', 'RCN12103', 'RCN12106', 'RCN12111',
        'RCN12112', 'RCN12115', 'RCN12116', 'RCN12676', 'RCN12688',
        'RCN12689', 'RCN12693', 'RCN12694', 'RCN12696', 'RCN12698',
        'RCN12700', 'RCN12704', 'RCN12710', 'RCN12712', 'RCN12716',
        'RCN12720', 'RCN12722', 'RCN12725', 'RCN12726', 'RCN12732',
        'RCN12733', 'RCN12734', 'RCN12736', 'RCN12737', 'RCN12738',
        'RCN12740', 'RCN12741', 'RCN12742', 'RCN12744', 'RCN12745']
    data = pd.DataFrame(
        {'sample':sample_ids}
    )
    # use a small chunk_size (default is 50) to trigger 3 searches
    ena_data = create_ena_data_frame(data, chunk_size=20)
    assert sorted(ena_data['sample'].unique()) == sorted(sample_ids), 'ENA results obtained for all sample IDs provided'


def test__build_ena_query():
    query = _build_ena_query(['abc123','def3456'])
    assert query=='(sample_title="abc123" OR sample_title="def3456")'
    
def test__parse_ena_response():
    response_text = "\n".join([
        'run_accession\tsample_title\tcenter_name\tlibrary_strategy\tsample_accession\tfastq_ftp\tsubmitted_ftp',
        'ERR14392568\tRCN15107\tThe GenRe-Mekong Project;GenRe-Mekong\tAMPLICON\tSAMEA117705075\tftp.sra.ebi.ac.uk/vol1/fastq/ERR143/068/ERR14392568/ERR14392568_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/068/ERR14392568/ERR14392568_2.fastq.gz\tftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_SPEC_29632.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_SPEC_29632.cram.crai',
        'ERR14390721\tRCN13560\tThe GenRe-Mekong Project;GenRe-Mekong\tAMPLICON\tSAMEA117704460\tftp.sra.ebi.ac.uk/vol1/fastq/ERR143/021/ERR14390721/ERR14390721_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/021/ERR14390721/ERR14390721_2.fastq.gz\tftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14390721/RCN13560_GRC1_27488.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14390721/RCN13560_GRC1_27488.cram.crai'
    ])
    ena_data = _parse_ena_response(response_text)
    assert isinstance(ena_data, pd.DataFrame), 'correctly parsed into a DataFrame'
    assert len(ena_data) == 2,'2 rows of data'
    assert 'sample_title' in ena_data, 'expected column "sample_title" exists'
    
def test__search_ena():
    """
    This test will run a query against the ENA search API using three
    real-world Pf8/GenRe sample IDs (="sample_title" in ENA).  
    If this test fails, make sure it is not an issue with your network connections.  
    """
    # use three existing sample IDs from Pf8 (sample titles in ENA)
    samples = ['RCN13568','RCN13560','RCN15107']
    fields = ['sample_title','run_accession','center_name','library_strategy','sample_accession','fastq_ftp','submitted_ftp']

    result = _search_ena(samples=samples,return_fields=fields)
    assert result, 'should receive some results'
    assert isinstance(result, str), 'returns raw text'
    assert "\n" in result, 'result string contains line breaks (mutliple rows)'
    assert "\t" in result, 'result string contains tab as delimiter'
    assert "sample_title" in result, 'contains sample_title column name'
    assert 'RCN13568' in result, 'contains an expected sample ID (title)'

def test_merge_ena_results_into_sample_data_genre_pf8():
    # mock ENA result DataFrame, created after a real-world query but modified for brevity. 
    ena_result = pd.DataFrame(
        {
            'run_accession': ['ERR14392568', 'ERR14390721', 'ERR14390722','ERR15626087'], 
            'sample': ['RCN15107', 'RCN13560', 'RCN13560','RCN13560'], 
            'center_name': ['GenRe-Mekong', 'GenRe-Mekong', 'GenRe-Mekong','Wellcome Sanger Institute;WSI'], 
            'library_strategy': ['AMPLICON', 'AMPLICON', 'AMPLICON','WGS'], 
            'fastq_ftp': ['some1/read1.fastq.gz;some1/read2.fastq.gz','some2/read1.fastq.gz;some2/read2.fastq.gz','some3/read1.fastq.gz;some3/read2.fastq.gz','some4/read1.fastq.gz;some4/read2.fastq.gz'],
            'submitted_ftp': ['some1/sample_SPEC.cram','some2/sample_GRC1.cram','some3/sample_GRC2.cram','some4/sample_wgs.cram']
        }
    )
        
    # mock sample data DataFrame with a sample column + 1 other column that isn't used 
    sample_data = pd.DataFrame(
        {'sample':['RCN15107', 'RCN13560'],'some_other_data': [1,2] }
    )

    merged_data = merge_ena_results_into_sample_data_genre_pf8(sample_data=sample_data, ena_result=ena_result)

    assert len(merged_data) == 2, 'the merged data has two rows, exactly as input sample data'
    assert merged_data.at['RCN15107','ENA_acc_GenRe_SPEC']=='ERR14392568','sample RCN15107 GenRE SPEC accession is correctly inferred from data'
    assert merged_data.at['RCN13560','ENA_acc_Pf8']=='ERR15626087','sample RCN13560 Pf8 accession is correctly inferred from data'
    assert merged_data.at['RCN13560','GenRe_GRC2_ENA_FASTQ_FTP_1']=='ftp://some3/read1.fastq.gz','sample RCN13560 Pf8 read 1 FTP link is correctly extracted from ENA data and protocol prepended'
    assert merged_data.at['RCN13560','GenRe_GRC2_ENA_FASTQ_FTP_2']=='ftp://some3/read2.fastq.gz','sample RCN13560 Pf8 read 1 FTP link is correctly extracted from ENA data and protocol prepended'

def test__download_fastq_file(tmp_path):
    # NOTE: this is a real SRA FTP path and relies on a file existing on the 
    # third-party resource. The test may fail because the file was removed remotely
    remote_ftp_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/003/ERR14388603/ERR14388603_1.fastq.gz'
    local_path = _download_fastq_file( remote_ftp_url, tmp_path)
    assert local_path, 'method returns the local file path'
    assert pathlib.Path(local_path).exists(), 'the file has been created'
    
    # Force an error by providing a non-existing URL
    remote_ftp_url = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/003/ERR14388603/THIS_DOES_NOT_EXIST.fastq.gz'
    with pytest.raises(Exception, match='Failed to download'):
        _download_fastq_file( remote_ftp_url, tmp_path)

def test_download_all_fastqs(tmp_path):
    data = pd.DataFrame(
        {
            'run_accession': [
                'ERR14388605',
                'ERR14388608'],
            'ftp_url_read_1': [
                'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/005/ERR14388605/ERR14388605_1.fastq.gz',
                'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/008/ERR14388608/ERR14388608_1.fastq.gz'
                ],
            'ftp_url_read_2': [
                'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/005/ERR14388605/ERR14388605_2.fastq.gz',
                'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/008/ERR14388608/ERR14388608_2.fastq.gz'
            ]
        }
    )
    download_all_fastqs(data=data,outdir=tmp_path)
    
    expected_manifest_path = tmp_path / 'manifest.csv'
    assert expected_manifest_path.exists(), 'a manifest file has been created'
    with open(expected_manifest_path,'r') as f:
        manifest_lines = list(f)
        assert manifest_lines[0].startswith('run_accession,'), 'first line of manifest file starts with "run_accession" field and looks like csv'
        assert manifest_lines[1].startswith('ERR14388605,'), 'first line of manifest data starts with expected run accession ID'
    
    assert (tmp_path / 'ERR14388605_1.fastq.gz').exists(), 'a FASTQ file was downloaded'
    assert (tmp_path / 'ERR14388608_1.fastq.gz').exists(), 'a FASTQ file was downloaded'
    assert (tmp_path / 'ERR14388605_2.fastq.gz').exists(), 'a FASTQ file was downloaded'
    assert (tmp_path / 'ERR14388608_2.fastq.gz').exists(), 'a FASTQ file was downloaded'

    
    with open(tmp_path / 'ERR14388605_1.fastq.gz', 'rb') as f:
        assert f.read(2) == b'\x1f\x8b', 'the file is gzipped'
    with gzip.open(tmp_path / 'ERR14388605_1.fastq.gz', 'rb') as f:
        file_content = f.read()
        assert file_content.decode("utf-8").startswith('@'), 'first line if FASTQ file starts with "@"'

    data.rename(columns={"run_accession": "accession"}, inplace=True)
    with pytest.raises(ValueError, match='data is missing columns'):
        download_all_fastqs(data=data,outdir=tmp_path)
        
    assert download_all_fastqs(data=data,outdir=tmp_path, run_accession_col='accession'),'non-default col name correctly applied'