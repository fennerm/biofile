"""BiofileGroup classes"""
from collections.abc import Sized
from typing import (
        List,
        Sequence,
        Type,
        )

from fmbiopy.fmcheck import all_equal
from plumbum import LocalPath

from biofile.file import Biofile


class BiofileGroup(Sized):
    """Superclass for storing and validating groups of bioinformatics files.

    All parameters except files are passed to `Biofile` to initialize each
    individual file.

    Attributes
    ----------
    gzipped : bool
        Same as parameter
    names : List[str]
        The unique IDs of each file.
    paths : List[str]
        The paths of the files
    cls : Type[Biofile]
        The stored `Biofile` class
    """
    input_type: List[str] = ['ANY*']
    extensions = ['ANY']

    def __init__(
            self,
            paths: Sequence[LocalPath],
            filetype: Type[Biofile],
            *args,
            **kwargs,
            ) -> None:
        """Initialize

        Parameters
        ----------
        files
            The list of files to be stored in the group.
        filetype
            The stored filetype
        """

        # Store paramaters
        self.paths = paths
        self.filetype = filetype

        self.biofiles = self._initialize_biofiles(*args, **kwargs)
        self.names = [f.name for f in self.biofiles]
        self.validate()

    def validate(self)-> bool:
        """Validation function to be used upon attempted access

        It is only called upon the first access attempt
        """
        self._checkpaths_not_none()
        self._check_extensions_same()
        return True

    def __getitem__(self, item) -> str:
        """Get a file from the group"""
        return self.paths[item]

    def __eq__(self, other) -> bool:
        """Test for BiofileGroup is equal to another"""

        # Cannot be equal if lengths are different
        if len(self) != len(other):
            return False

        # Test that all paths are the same
        for me, you in zip(self, other):  # type: ignore
            if me != you:
                return False

        return True

    def _initialize_biofiles(self, *args, **kwargs) -> List[Biofile]:
        """Initalize a set of biofiles for the input file list"""
        return [self.filetype(p, *args, **kwargs) for p in self.paths]

    def __len__(self) -> int:
        """Length of file list"""
        return len(self.paths)

    def _checkpaths_not_none(self) -> bool:
        """Check paths not an empty list"""
        if not self.paths:
            raise ValueError('Empty paths in BiofileGroup')
        return True

    def _check_extensions_same(self) -> bool:
        """Check that the stored file extensions are all the same"""
        extensions = [f.extension for f in self.biofiles]
        if not all_equal(extensions):
            raise FileExtensionsNotSameError(self.paths)
        return True

#
#
#
# class Bowtie2Indices(BiofileGroup):
#     """Biofile class for holding Bowtie2 .bt2 fasta index files"""
#     input_type = ['bt2']
#     accepted_extension = ['bt2']
#
#     def __init__(self, *args, **kwargs):
#         """Initalize class"""
#         super().__init__()
#         self._names = None
#         self._name = self._get_name()
#
#     @property
#     def name(self) -> Sequence[str]:
#         """Simple getter function to retrieve the paths of all index
#            files."""
#
#         return self.paths
#
#     def _get_name(self) -> str:
#         """Get the prefixes of the bowtie2 indices"""
#
#         # Reverse bowtie2 indices need three extensions removed, the rest
#         # just need two
#         path = os.path.basename(self.paths[0])
#         if '.rev.' in path:
#             return fmpaths.remove_suffix(path, 3)[0]
#         return fmpaths.remove_suffix(path, 2)[0]
#
#     def _check_extensions(self):
#         assert False
#
#     def __getitem__(self, item) -> str:
#         """Get an index prefix from the list """
#         if not self.validated:
#             self.validate()
#
#         return self._index_prefixes[item]
#
#
# """
# module TypeVar which groups different types of fasta index classes
# into a single type
# """
# FastaIndexGroup = Union[SamtoolsFAIndex, Bowtie2Indices]
#
#

# -----------------------------------------------------------------------------
# Exceptions
# -----------------------------------------------------------------------------

class BiofileGroupValidationError(Exception):
    """Basic exception for errors raised by `BiofileGroup` validation checks

    Parameters
    ----------
    name
        The filenames which caused the error

    Attributes
    ----------
    msg
        The formatted error message
    """
    def __init__(self, names: Sequence[LocalPath] = None) -> None:
        self.names = names
        super().__init__()

    def _formatted_filenames(self) -> str:
        """Construct the part of the error message which lists the file"""
        return '\n'.join(['Files:', '\t' + ', '.join(str(self.names))])


class FileExtensionsNotSameError(BiofileGroupValidationError):
    """Exception raised when files do not have the expected file extension"""

    def _err_description(self) -> str:
        return "File extensions are not all equal"
