import argparse

import pytest

from tsarina import cli_common


def test_split_csv_trims_and_drops_empty():
    assert cli_common.split_csv(" A , ,B,, C ") == ["A", "B", "C"]
    assert cli_common.split_csv("") == []


def test_parse_lengths_ok_and_error():
    assert cli_common.parse_lengths("8,9, 10") == (8, 9, 10)
    with pytest.raises(argparse.ArgumentTypeError):
        cli_common.parse_lengths("8,x")


def test_iedb_cedar_args_map_to_path_dests():
    p = argparse.ArgumentParser()
    cli_common.add_iedb_cedar_args(p)
    assert p.parse_args(["--iedb", "a"]).iedb_path == "a"
    assert p.parse_args(["--cedar", "c"]).cedar_path == "c"
    # The historical *-path spellings are gone (single canonical spelling).
    with pytest.raises(SystemExit):
        p.parse_args(["--iedb-path", "b"])


def test_add_predictor_arg_choices_and_default():
    p = argparse.ArgumentParser()
    cli_common.add_predictor_arg(p, context="testing")
    assert p.parse_args([]).predictor == "mhcflurry"
    with pytest.raises(SystemExit):
        p.parse_args(["--predictor", "nope"])


def test_warn_ignored_flags(capsys):
    cli_common.warn_ignored_flags("--healthy-tissue", ["--allele", "--format"])
    err = capsys.readouterr().err
    assert "--healthy-tissue" in err and "--allele" in err and "--format" in err
    # No output when nothing is ignored.
    cli_common.warn_ignored_flags("--healthy-tissue", [])
    assert capsys.readouterr().err == ""
