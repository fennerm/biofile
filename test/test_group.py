"""Test suite for biofile.group"""
from pathlib import Path

from pytest import (
        fixture,
        raises,
        )

from fmbiopy.fmpaths import as_paths

from biofile.file import *
from biofile.group import *


@fixture()
def diff_prefix(dat):
    return BiofileGroup(dat['tiny']['diff_prefix'], Fasta)


@fixture()
def fwd_reads(dat):
    return BiofileGroup(dat['tiny']['fwd_reads'], FwdFastq)


@fixture()
def rev_reads(dat):
    return BiofileGroup(dat['tiny']['rev_reads'], RevFastq)


@fixture()
def assemblies(dat):
    return BiofileGroup(dat['tiny']['assemblies'], Fasta)


class TestBiofileGroup(object):

    def test_empty_input_raises_value_err(self):
        with raises(ValueError):
            BiofileGroup([], filetype=Fasta)

    def test_access_returns_path(self, dat, assemblies):
        fasta_paths = dat['tiny']['assemblies']
        expect = fasta_paths[0]
        assert assemblies[0] == expect

    def test_length_method(self, dat, assemblies):
        fasta_paths = dat['tiny']['assemblies']
        actual = len(assemblies)
        expect = len(fasta_paths)
        assert actual == expect

    def test_single_path_input_raises_type_err(self):
        with raises(TypeError):
            BiofileGroup(Path('foo.fa'), filetype=Fasta)

    def test_biofilegroups_can_be_zipped(self, assemblies, fwd_reads):
        max_index = max(len(assemblies), len(fwd_reads))
        for fa, reads, i in zip(assemblies, fwd_reads, range(max_index)):
            assert fa == assemblies[i]
            assert reads == fwd_reads[i]

    def test_length_nonexist_doesnt_raise_error(self):
        paths = as_paths(['foo.fa', 'bar.fa'])
        len(BiofileGroup(paths, filetype=Fasta))

    def test_equality_operator(self, assemblies, diff_prefix):
        assert assemblies != diff_prefix

    def test_different_extensions_raises_value_err(self):
        with raises(FileExtensionsNotSameError):
            paths = as_paths(['a.fa', 'b.fasta', 'c.fa'])
            BiofileGroup(paths, filetype=Fasta)

    def test_optional_param_are_passed_to_biofile(self, dat):
        with raises(GzipStatusError):
            BiofileGroup(
                    dat['tiny']['fwd_reads'],
                    filetype=Fastq,
                    gzipped=True)
