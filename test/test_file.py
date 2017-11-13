"""Test suite for biofile.file"""
from pathlib import Path

from typing import (
        Callable,
        Type,
        )

from pytest import (
        fixture,
        raises,
        )

from fmbiopy.fmclass import list_classes
from fmbiopy.fmpaths import (add_suffix)

from biofile.file import *

@fixture(scope='module')
def instance_of(example_file: Callable[[Type[Biofile], str], Path]):
    def _make_test_instance(
            cls: Type[Biofile],
            size: str)-> Biofile:

        input_example = example_file(cls, size)

        return cls(input_example)
    return _make_test_instance


@fixture(
        scope='module',
        params=list_classes(
            'file',
            'biofile',
            of_type=['Biofile']))
def biofiles(request):
    return request.param


@fixture(scope='module')
def inst_biofiles(instance_of, biofiles):
    """Return instances of all Biofile types"""
    return instance_of(biofiles, 'tiny')


class TestBiofile(object):
    """Test `Biofile` class"""

    def test_extension_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'extensions')

    def test_input_type_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'input_type')

    def test_name_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'name')

    def test_path_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'path')

    def test_gzipped_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'gzipped')

    def test_empty_input_raises_value_err(self, biofiles):
        with raises(TypeError):
            biofiles(Path(''))

    def test_undeclared_gzip_raises_gzip_error(self, dat):
        with raises(GzipStatusError):
            Biofile(dat['tiny']['zipped_fwd_reads'][0]).validate()

    def test_list_input_raises_type_error(self, example_file, biofiles):
        with raises(AttributeError):
            biofiles([example_file(biofiles.input_type, 'tiny')])

    def test_if_files_dont_exist_raise_err(self):
        with raises(FileNotFoundError):
            Biofile(Path('i_dont_exist.fa')).validate()

    def test_empty_files_raises_err(self, dat):
        bf = Biofile(dat['tiny']['empty_reads'][0])
        with raises(EmptyFileError):
            bf.validate()

    def test_possibly_empty_prevents_error(self, dat):
        bf = Biofile(dat['tiny']['empty_reads'][0], possibly_empty=True)
        assert bf.validate()

    def test_incorrect_extension_raises_extension_err(self, biofiles, dat):
        read_path = dat['tiny']['fwd_reads'][0]
        incorrect_suffix = add_suffix(read_path, '.foobar')

        if biofiles.extensions != ['ANY']:
            with raises(FileExtensionError):
                biofiles(incorrect_suffix)
