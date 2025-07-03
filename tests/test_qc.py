from pathlib import Path

import pytest

from vibecheck.src.qc import (
    check_barcodes,
    check_parse_float_fraction,
    check_query_file,
    check_threads,
    check_tree,
    setup_outdir,
    setup_tempdir,
)


def test_parse_fraction_valid():
    value = 0.3
    got = check_parse_float_fraction(0.3, "foo")
    assert got == value


def test_parse_fraction_value_int():
    value = 30
    got = check_parse_float_fraction(value, "foo")
    assert got == 0.3


def test_parse_fraction_invalid():
    value = 1204123
    with pytest.raises(SystemExit) as e:
        check_parse_float_fraction(value, "foo")
    assert e.value.code != 0


def test_parse_fraction_negative():
    value = -0.3
    with pytest.raises(SystemExit) as e:
        check_parse_float_fraction(value, "foo")
    assert e.value.code != 0


def test_check_threads_valid():
    got = check_threads(1, 4)
    assert got == 1


def test_check_threads_invalid():
    with pytest.raises(SystemExit) as e:
        check_threads(10, 3)
    assert e.value.code != 0


def test_check_threads_negative():
    with pytest.raises(SystemExit) as e:
        check_threads(-42, 10)
    assert e.value.code != 0


def test_check_query_file_single_file(tmp_path):
    file = tmp_path / "query.fasta"
    file.touch()  # Create a dummy file
    result, _ = check_query_file([str(file)])
    assert result == file


def test_check_query_file_multiple_files(tmp_path):
    fileA = tmp_path / "file1.fasta"
    fileB = tmp_path / "file2.fasta"
    fileA.touch()
    fileB.touch()

    with pytest.raises(SystemExit) as e:
        check_query_file([str(fileA), str(fileB)])
    assert e.value.code == -5


def test_check_query_file_nonexistent_file():
    with pytest.raises(SystemExit) as e:
        check_query_file(["nonexistent.fasta"])
    assert e.value.code == -2


def test_check_query_file_empty_list():
    with pytest.raises(SystemExit) as e:
        check_query_file([])
    assert e.value.code == -3


def test_check_query_file_relative_path():
    file = "tests/example_fasta/ERR037738.fasta"
    result, _ = check_query_file([file])
    assert result == Path(file).absolute()


def test_check_query_file_absolute_path():
    file = Path("tests/example_fasta/ERR037738.fasta").absolute()
    result, _ = check_query_file([str(file)])
    assert result == Path(file)


def test_check_query_file_tilde_path():
    file = Path("tests/example_fasta/ERR037738.fasta")
    home = Path.home()
    file_home_relative = "~/" / file.absolute().relative_to(home.expanduser())
    result, _ = check_query_file([str(file_home_relative)])
    assert result == Path(file_home_relative).expanduser()


def test_check_query_file_detect_fasta(tmp_path):
    file = tmp_path / "query.fasta"
    file.touch()
    result, use_usher = check_query_file([str(file)])
    assert result == file
    assert use_usher


def test_check_query_file_dectect_fastq(tmp_path):
    fileA = tmp_path / "read1.fastq.gz"
    fileB = tmp_path / "read2.fastq.gz"
    fileA.touch()
    fileB.touch()

    result, use_usher = check_query_file([str(fileA), str(fileB)])
    assert result == (fileA, fileB)
    assert not use_usher


def test_check_query_file_dectect_single_fastq(tmp_path):
    fileA = tmp_path / "read1.fastq.gz"
    fileA.touch()

    result, use_usher = check_query_file([str(fileA)])
    assert result == (fileA, None)
    assert not use_usher


def test_check_query_file_dectect_fq(tmp_path):
    fileA = tmp_path / "read1.fq.gz"
    fileB = tmp_path / "read2.fq.gz"
    fileA.touch()
    fileB.touch()

    result, use_usher = check_query_file([str(fileA), str(fileB)])
    assert result == (fileA, fileB)
    assert not use_usher


def test_check_query_file_error_mixed_type(tmp_path):
    fileA = tmp_path / "read1.fastq.gz"
    fileB = tmp_path / "reads.fasta"
    fileA.touch()
    fileB.touch()

    with pytest.raises(SystemExit) as e:
        result, use_usher = check_query_file([str(fileA), str(fileB)])
    assert e.value.code == -4


def test_check_query_file_error_too_many_fastqs(tmp_path):
    fileA = tmp_path / "read1.fastq.gz"
    fileB = tmp_path / "read2.fastq.gz"
    fileC = tmp_path / "read3.fastq.gz"
    fileA.touch()
    fileB.touch()
    fileC.touch()

    with pytest.raises(SystemExit) as e:
        result, use_usher = check_query_file([str(fileA), str(fileB), str(fileC)])
    assert e.value.code == -6


def test_check_query_file_error_missing_fastq(tmp_path):
    fileA = tmp_path / "read1.fastq.gz"
    fileB = tmp_path / "read2.fastq.gz"
    fileA.touch()

    with pytest.raises(SystemExit) as e:
        result, use_usher = check_query_file([str(fileA), str(fileB)])
    assert e.value.code == -2


def test_check_tree_valid_file(tmp_path):
    tree = tmp_path / "tree.pb"
    tree.touch()  # Create a dummy file
    result = check_tree(str(tree))
    assert result == tree


def test_check_tree_nonexistent_file():
    with pytest.raises(SystemExit) as e:
        check_tree("nonexistent.pb")
    assert e.value.code == -4


def test_check_barcodes_valid_file(tmp_path):
    tree = tmp_path / "barcodes.feather"
    tree.touch()  # Create a dummy file
    result = check_barcodes(str(tree))
    assert result == tree


def test_check_barcodes_nonexistent_file():
    with pytest.raises(SystemExit) as e:
        check_barcodes("nonexistent.feather")
    assert e.value.code != 0


def test_setup_outdir_existing_dir(tmp_path):
    result = setup_outdir(str(tmp_path), str(tmp_path))
    assert result == tmp_path


def test_setup_outdir_new_dir(tmp_path):
    new_dir = tmp_path / "new_dir"
    result = setup_outdir(str(new_dir), str(tmp_path))
    assert result == new_dir
    assert new_dir.exists()


def test_setup_outdir_invalid_dir():
    with pytest.raises(SystemExit) as e:
        setup_outdir("/invalid_dir", "/cwd")
    assert e.value.code == -5


def test_setup_outdir_no_outdir(tmp_path):
    result = setup_outdir("", str(tmp_path))
    assert result == Path(tmp_path)


def test_setup_tempdir_no_temp(tmp_path):
    result = setup_tempdir("", tmp_path, True)
    assert result == tmp_path


def test_setup_tempdir_new_tempdir(tmp_path):
    tempdir = tmp_path / "temp_dir"
    result = setup_tempdir(str(tempdir), tmp_path, False)
    assert result == tempdir
    assert tempdir.exists()


def test_setup_tempdir_invalid_tempdir():
    with pytest.raises(SystemExit) as e:
        setup_tempdir("/invalid_tempdir", Path("/outdir"), False)
    assert e.value.code == -6


def test_setup_tempdir_writable_tempdir(tmp_path):
    tempdir = setup_tempdir("", tmp_path, False)
    test_file = tempdir / "test.txt"
    test_file.write_text("This is a test")
    assert test_file.exists()
    assert test_file.read_text() == "This is a test"


def test_setup_tempdir_relative_path(tmp_path):
    tempdir = Path("tests/example_fasta/")
    home = Path.home()
    tempdir_home_relative = "~/" / tempdir.absolute().relative_to(home.expanduser())
    result = setup_tempdir(str(tempdir_home_relative), tempdir_home_relative, False)
    assert result == Path(tempdir_home_relative).expanduser()


def test_check_lineage_aliases(tmp_path):
    alias_file = tmp_path / "lineage_aliases.csv"
    alias_file.write_text(
        "alias,lineage\n"
        "a,1\nb,2\nc,3\n"
    )

    got = check_lineage_aliases(str(alias_file))
    want = {"a": "1", "b": "2", "c": "3"}
    assert got == want

def test_check_lineage_aliases_exits_with_nonexistant_file(tmp_path):
    alias_file = tmp_path / "lineage_aliases.csv"

    with pytest.raises(SystemExit) as e:
        check_lineage_aliases( str(alias_file) )
    assert e.value.code != 0
