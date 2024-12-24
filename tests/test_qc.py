from pathlib import Path

import pytest

from phonebook.src.qc import check_query_file, check_tree, setup_outdir, setup_tempdir


def test_check_query_file_single_file(tmp_path):
    file = tmp_path / "query.fasta"
    file.touch()  # Create a dummy file
    result = check_query_file([str(file)])
    assert result == file


def test_check_query_file_multiple_files():
    with pytest.raises(SystemExit) as e:
        check_query_file(["file1.fasta", "file2.fasta"])
    assert e.value.code == -1


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
    result = check_query_file([file])
    assert result == Path(file)


def test_check_query_file_absolute_path():
    # Please ignore my absolute path...
    file = "/Users/natem/Documents/Code/Killington/tests/example_fasta/ERR037738.fasta"
    result = check_query_file([file])
    assert result == Path(file)


def test_check_query_file_tilde_path():
    # I don't want any comments about my file structure. I'm not sure how else to test
    # this.
    file = "~/Documents/Code/Killington/tests/example_fasta/ERR037738.fasta"
    result = check_query_file([file])
    assert result == Path(file).expanduser()


def test_check_tree_valid_file(tmp_path):
    tree = tmp_path / "tree.pb"
    tree.touch()  # Create a dummy file
    result = check_tree(str(tree))
    assert result == tree


def test_check_tree_nonexistent_file():
    with pytest.raises(SystemExit) as e:
        check_tree("nonexistent.pb")
    assert e.value.code == -4


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
    tempdir = "~/Documents/Code/Killington/tests/example_fasta/"
    result = setup_tempdir(str(tempdir), Path(tempdir), False)
    assert result == Path(tempdir).expanduser()
