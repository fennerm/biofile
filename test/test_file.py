"""Test suite for biofile.file"""

from typing import (
        Callable,
        Type,
        )

from plumbum import LocalPath
from pytest import (
        fixture,
        raises,
        )

from fmbiopy.fmclass import list_classes

from biofile.file import *

@fixture(scope='module')
def instance_of(example_file: Callable[[Type[Biofile], str], LocalPath]):
    def _make_test_instance(
            cls: Type[Biofile],
            size: str)-> Biofile:

        input_example = example_file(cls, size)

        return cls(str(input_example))
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

    def test_gzipped_set(self, inst_biofiles):
        assert hasattr(inst_biofiles, 'gzipped')

    def test_empty_input_raises_err(self, biofiles):
        with raises(TypeError):
            biofiles()

    def test_path_methods(self, inst_biofiles):
        inst_biofiles.suffix
        inst_biofiles.is_dir()

    def test_undeclared_gzip_raises_gzip_error(self, dat):
        with raises(GzipStatusError):
            Biofile(str(dat['tiny']['zipped_fwd_reads'][0]))

    def test_list_input_raises_type_error(self, example_file, biofiles):
        with raises(AttributeError):
            biofiles([example_file(biofiles.input_type, 'tiny')])

    def test_if_files_dont_exist_raise_err(self):
        with raises(FileNotFoundError):
            Biofile('i_dont_exist.fa')

    def test_empty_files_raises_err(self, dat):
        with raises(EmptyFileError):
            Biofile(str(dat['tiny']['empty_reads'][0]))

    def test_possibly_empty_prevents_error(self, dat):
        Biofile(str(dat['tiny']['empty_reads'][0]), possibly_empty=True)

    def test_determine_gzip(self, dat):
        Fastq(dat['tiny']['zipped_fwd_reads'][0], determine_gzip=True)

    def test_incorrect_extension_raises_extension_err(
            self,
            gen_tmp,
            biofiles):

        incorrect_suffix = gen_tmp(empty=False, suffix='.foobar')

        if biofiles.extensions != ['ANY']:
            with raises(FileExtensionError):
                biofiles(incorrect_suffix)
