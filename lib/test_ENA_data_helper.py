import pathlib
from urllib.error import URLError

import pandas as pd
import pytest

import ENA_data_helper as edh


def test_create_ena_data_frame_from_samples(monkeypatch):
    data = pd.DataFrame(
        {"sample": ["RCN13568", "RCN13560", "RCN15107"], "some_other_data": [1, 2, 3]}
    )
    response_text = "\n".join(
        [
            "run_accession\texperiment_accession\tstudy_accession\tsample_title\tread_count\tcenter_name\tlibrary_strategy\tsample_accession",
            "RUN1\tEXP1\tSTUDY1\tRCN13568\t100\tCenterA\tAMPLICON\tS1",
            "RUN2\tEXP2\tSTUDY2\tRCN13560\t120\tCenterB\tWGS\tS2",
            "RUN3\tEXP3\tSTUDY3\tRCN15107\t140\tCenterC\tAMPLICON\tS3",
        ]
    )

    monkeypatch.setattr(edh, "_search_ena", lambda query, return_fields, limit=10000: response_text)

    ena_data = edh.create_ena_data_frame_from_samples(data, sample_id_col_name="sample")
    assert isinstance(ena_data, pd.DataFrame)
    assert len(ena_data) == 3

    expected_col_names = [
        "sample",
        "run_accession",
        "experiment_accession",
        "study_accession",
        "read_count",
        "center_name",
        "library_strategy",
        "sample_accession",
    ]
    assert sorted(list(ena_data.columns.values)) == sorted(expected_col_names), (
        "DataFrame has correct column (field) names"
    )


def test_create_ena_data_frame_from_samples_large(monkeypatch):
    """
    Special case: test with a large number of sample IDs, which would throw a
    HTTPError: 414 Client Error: Request-URI Too Large for url.
    This tests the query-in-chunks behaviour.
    """
    sample_ids = [
        "RCN12025",
        "RCN12026",
        "RCN12028",
        "RCN12031",
        "RCN12032",
        "RCN12033",
        "RCN12035",
        "RCN12036",
        "RCN12040",
        "RCN12041",
        "RCN12043",
        "RCN12044",
        "RCN12051",
        "RCN12058",
        "RCN12059",
        "RCN12060",
        "RCN12070",
        "RCN12071",
        "RCN12072",
        "RCN12075",
        "RCN12076",
        "RCN12079",
        "RCN12080",
        "RCN12082",
        "RCN12083",
        "RCN12084",
        "RCN12100",
        "RCN12103",
        "RCN12106",
        "RCN12111",
        "RCN12112",
        "RCN12115",
        "RCN12116",
        "RCN12676",
        "RCN12688",
        "RCN12689",
        "RCN12693",
        "RCN12694",
        "RCN12696",
        "RCN12698",
        "RCN12700",
        "RCN12704",
        "RCN12710",
        "RCN12712",
        "RCN12716",
        "RCN12720",
        "RCN12722",
        "RCN12725",
        "RCN12726",
        "RCN12732",
        "RCN12733",
        "RCN12734",
        "RCN12736",
        "RCN12737",
        "RCN12738",
        "RCN12740",
        "RCN12741",
        "RCN12742",
        "RCN12744",
        "RCN12745",
    ]
    data = pd.DataFrame({"sample": sample_ids})
    calls = []

    def fake_search(query, return_fields, limit=10000):
        calls.append(query)
        return (
            "run_accession\texperiment_accession\tstudy_accession\tsample_title\tread_count\tcenter_name\tlibrary_strategy\tsample_accession\n"
            "RUNX\tEXPX\tSTUDYX\tRCN12025\t10\tCenter\tAMPLICON\tSAMEA1\n"
        )

    monkeypatch.setattr(edh, "_search_ena", fake_search)

    ena_data = edh.create_ena_data_frame_from_samples(data, chunk_size=20)
    assert len(calls) == 3, "60 samples with chunk_size=20 should perform 3 ENA searches"
    assert len(ena_data) == 3


def test__build_ena_sample_query():
    query = edh._build_ena_sample_query(["abc123", "def3456"])
    assert query == '(sample_title="abc123" OR sample_title="def3456")'


def test__parse_ena_response():
    response_text = "\n".join(
        [
            "run_accession\tsample_title\tcenter_name\tlibrary_strategy\tsample_accession\tfastq_ftp\tsubmitted_ftp",
            "ERR14392568\tRCN15107\tThe GenRe-Mekong Project;GenRe-Mekong\tAMPLICON\tSAMEA117705075\tftp.sra.ebi.ac.uk/vol1/fastq/ERR143/068/ERR14392568/ERR14392568_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/068/ERR14392568/ERR14392568_2.fastq.gz\tftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_SPEC_29632.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14392568/RCN15107_SPEC_29632.cram.crai",
            "ERR14390721\tRCN13560\tThe GenRe-Mekong Project;GenRe-Mekong\tAMPLICON\tSAMEA117704460\tftp.sra.ebi.ac.uk/vol1/fastq/ERR143/021/ERR14390721/ERR14390721_1.fastq.gz;ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/021/ERR14390721/ERR14390721_2.fastq.gz\tftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14390721/RCN13560_GRC1_27488.cram;ftp.sra.ebi.ac.uk/vol1/run/ERR143/ERR14390721/RCN13560_GRC1_27488.cram.crai",
        ]
    )
    ena_data = edh._parse_ena_response(response_text)
    assert isinstance(ena_data, pd.DataFrame), "correctly parsed into a DataFrame"
    assert len(ena_data) == 2, "2 rows of data"
    assert "sample_title" in ena_data, 'expected column "sample_title" exists'


def test__search_ena(monkeypatch):
    class FakeResponse:
        text = "sample_title\trun_accession\nRCN13568\tRUN1\n"

        @staticmethod
        def raise_for_status():
            return None

    def fake_get(url, params, timeout):
        assert "ena/portal/api/search" in url
        assert params["result"] == "read_run"
        assert params["query"] == '(sample_title="RCN13568")'
        assert params["fields"] == "sample_title,run_accession"
        assert timeout == 15
        return FakeResponse()

    monkeypatch.setattr(edh.requests, "get", fake_get)

    result = edh._search_ena(
        query='(sample_title="RCN13568")',
        return_fields=["sample_title", "run_accession"],
    )
    assert result, "should receive some results"
    assert isinstance(result, str), "returns raw text"
    assert "\n" in result, "result string contains line breaks (multiple rows)"
    assert "\t" in result, "result string contains tab as delimiter"
    assert "sample_title" in result, "contains sample_title column name"
    assert "RCN13568" in result, "contains an expected sample ID (title)"


def test_align_ena_results_with_sample_data_genre_pf8():
    ena_result = pd.DataFrame(
        {
            "run_accession": [
                "ERR14392568",
                "ERR14390721",
                "ERR14390722",
                "ERR14390723",
                "ERR15626087",
            ],
            "sample": ["RCN15107", "RCN13560", "RCN13560", "RCN13560", "RCN13560"],
            "center_name": [
                "The GenRe-Mekong Project;GenRe-Mekong",
                "The GenRe-Mekong Project;GenRe-Mekong",
                "The GenRe-Mekong Project;GenRe-Mekong",
                "The GenRe-Mekong Project;GenRe-Mekong",
                "Wellcome Sanger Institute;WSI",
            ],
            "library_strategy": ["AMPLICON", "AMPLICON", "AMPLICON", "AMPLICON", "WGS"],
            "read_count": [10000, 12000, 14000, 16000, 20000],
        }
    )
    sample_data = pd.DataFrame({"sample": ["RCN15107", "RCN13560"], "some_other_data": [1, 2]})
    genre_panel_map = {
        "ERR14392568": "SPEC",
        "ERR14390721": "GRC1",
        "ERR14390722": "GRC2",
        "ERR14390723": "SPEC",
    }

    merged_data, filtered_ena = edh.align_ena_results_with_sample_data_genre_pf8(
        sample_data=sample_data, ena_result=ena_result, genre_panel_map=genre_panel_map
    )

    assert len(merged_data) == 1, (
        "only RCN13560 has all expected run accessions and should be retained"
    )
    assert merged_data.iloc[0]["sample"] == "RCN13560"
    assert merged_data.iloc[0]["INSDC_GenRe_GRC1"] == "ERR14390721"
    assert merged_data.iloc[0]["INSDC_GenRe_GRC2"] == "ERR14390722"
    assert merged_data.iloc[0]["INSDC_GenRe_SPEC"] == "ERR14390723"
    assert merged_data.iloc[0]["INSDC_Pf8"] == "ERR15626087"
    assert merged_data.iloc[0]["INSDC_GenRe_GRC1_readcount"] == 6000
    assert "library_name" in filtered_ena.columns


def test__download_fastq_file(tmp_path, monkeypatch):
    def fake_urlretrieve(remote, local):
        pathlib.Path(local).write_bytes(b"@fake\n")
        return local, None

    monkeypatch.setattr(edh.urllib.request, "urlretrieve", fake_urlretrieve)

    remote_ftp_url = (
        "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/003/ERR14388603/ERR14388603_1.fastq.gz"
    )
    local_path = edh._download_fastq_file(remote_ftp_url, tmp_path, num_tries=1)
    assert local_path, "method returns the local file path"
    assert pathlib.Path(local_path).exists(), "the file has been created"

    def fake_urlretrieve_raise(remote, local):
        raise URLError("network problem")

    monkeypatch.setattr(edh.urllib.request, "urlretrieve", fake_urlretrieve_raise)
    remote_ftp_url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR143/003/ERR14388603/THIS_DOES_NOT_EXIST.fastq.gz"
    with pytest.raises(Exception, match="Failed to download"):
        edh._download_fastq_file(remote_ftp_url, tmp_path, num_tries=1)


def test_download_all_fastqs_with_mocked_ena_lookup(tmp_path, monkeypatch):
    run_data = pd.DataFrame({"run_accession": ["ERR111", "ERR222"]})
    ena_locations = pd.DataFrame(
        {
            "run_accession": ["ERR111", "ERR222"],
            "fastq_ftp": [
                "ftp.sra.ebi.ac.uk/r1_1.fastq.gz;ftp.sra.ebi.ac.uk/r1_2.fastq.gz",
                "ftp.sra.ebi.ac.uk/r2_1.fastq.gz;ftp.sra.ebi.ac.uk/r2_2.fastq.gz",
            ],
        }
    )
    seen_downloads = []

    def fake_create_ena_data_frame_from_run_accessions(df, run_acc_col_name="run_accession", chunk_size=50):
        assert run_acc_col_name == "run_accession"
        assert chunk_size == 50
        assert list(df["run_accession"]) == ["ERR111", "ERR222"]
        return ena_locations

    def fake_download(url, outdir, num_tries=3, skip_errors=False):
        seen_downloads.append((url, str(outdir), num_tries, skip_errors))
        return str(pathlib.Path(outdir) / url.split("/")[-1])

    monkeypatch.setattr(edh, "create_ena_data_frame_from_run_accessions", fake_create_ena_data_frame_from_run_accessions)
    monkeypatch.setattr(edh, "_download_fastq_file", fake_download)

    assert edh.download_all_fastqs(data=run_data, outdir=tmp_path)
    assert len(seen_downloads) == 4, "two FASTQ URLs per run accession should be downloaded"
    assert seen_downloads[0][0] == "ftp://ftp.sra.ebi.ac.uk/r1_1.fastq.gz"
    assert seen_downloads[-1][0] == "ftp://ftp.sra.ebi.ac.uk/r2_2.fastq.gz"


def test_download_all_fastqs_with_missing_data_raises(tmp_path):
    with pytest.raises(ValueError, match='must provide either "data" or "data_file_path" parameter'):
        edh.download_all_fastqs(outdir=tmp_path)


def test_cli_download_fastqs_uses_output_folder(monkeypatch, tmp_path):
    called_with = {}

    def fake_download_all_fastqs(outdir, data_file_path, num_tries, skip_errors, top3):
        called_with["outdir"] = outdir
        called_with["data_file_path"] = data_file_path
        called_with["num_tries"] = num_tries
        called_with["skip_errors"] = skip_errors
        called_with["top3"] = top3
        return True

    monkeypatch.setattr(edh, "download_all_fastqs", fake_download_all_fastqs)

    class Args:
        output_folder = str(tmp_path)
        insdc_manifest = "manifest.csv"
        download_attempts = 5
        skip_errors = True
        top3 = False

    rc = edh.cli_download_fastqs(Args())
    assert rc == 0
    assert called_with == {
        "outdir": str(tmp_path),
        "data_file_path": "manifest.csv",
        "num_tries": 5,
        "skip_errors": True,
        "top3": False,
    }
