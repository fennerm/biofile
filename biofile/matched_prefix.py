"""MatchedPrefixGroup classes"""

from typing import (
        List,
        Sequence,
        )

from fmbiopy.fmcheck import all_equal
from fmbiopy.fmpaths import prefix

from biofile.group import (
        BiofileGroup,
        BiofileGroupValidationError,
        )


class MatchedPrefixGroup(object):
    """Stores groups of matched BiofileGroups with same prefix

    Prefixes are checked for equality. Groups are checked for equal length.

    Parameters
    ----------
    groups
        List of matched BiofileGroups

    Attributes
    ----------
    groups : List[BiofileGroup]
    """
    input_type = ['ANY*']
    extensions = ['ANY']

    def __init__(self, groups: List[BiofileGroup]) -> None:
        self.groups = groups
        self.validated = False
        self.validate()

    def validate(self)-> bool:
        """Validate the `BiofileGroup`s"""
        self.__check_files_not_same()
        self.__check_lengths_match()
        self.__check_same_file_prefix()
        if not self.validated:
            for group in self.groups:
                group.validate()
        self.validated = True
        return True

    def __check_files_not_same(self) -> None:
        """Check that none of the filegroups are the exact same"""
        for i, group1 in enumerate(self.groups):
            for j, group2 in enumerate(self.groups):
                if i != j and group1 == group2:
                    raise DuplicateFilegroupError(self.groups)

    def __check_lengths_match(self) -> None:
        group_lengths = [len(g) for g in self.groups]
        if not all_equal(group_lengths):
            raise GroupLengthError(self.groups)

    def __check_same_file_prefix(self) -> None:
        """Check that the stored BiofileGroups all have the same prefixes"""
        group_paths = [group._paths for group in self.groups]
        prefixes = []
        for paths in group_paths:
            prefixes.append([prefix(path) for path in paths])
        if not all_equal(prefixes):
            raise PrefixMatchError(self.groups)

    def __len__(self) -> int:
        """Length of the `MatchedPrefixGroup`"""
        return len(self.groups[0])

    def __getitem__(self, item) -> List[str]:
        """Index the `MatchedPrefixGroup`"""
        return [g[item] for g in self.groups]
#
#
# class PairedFastqGroup(MatchedPrefixGroup):
#     """Stores two groups of paired FastqGroup files
#
#     Parameters
#     ----------
#     forward_fastq, reverse_fastq
#         Fastq objects to be paired
#
#     Attributes
#     ----------
#     forward_fastq - FastqGroup, reverse_fastq - FastqGroup
#         Same as Parameters
#
#     """
#
#     input_type = ['fastq', 'fastq']
#
#
# class IndexedFastaGroup(object):
#     """Represents a Fasta file grouped with its indices"""
#     def __init__(
#             self,
#             fasta: FastaGroup,
#             *indices: FastaIndexGroup) -> None:
#         self.fasta = fasta
#         self.indices = indices
#
#     def __getitem__(self, item) -> List[str]:
#         index_items = [index[item] for index in self.indices]
#         return [self.fasta[item]] + index_items

# -----------------------------------------------------------------------------
# Exceptions
# -----------------------------------------------------------------------------


class MatchedPrefixGroupValidationError(BiofileGroupValidationError):
    """Basic exception for errors raised by `MatchedPrefixGroup` validation

    Parameters
    ----------
    groups
        Nested list of filenames involved in the error.
    """
    def __init__(self, groups: Sequence[BiofileGroup] = None) -> None:
        self.groups = groups
        super().__init__()

    def _formatted_filenames(self) -> str:
        """Format the file group names for printing"""
        formatted_filenames = 'Groups:\n'

        for group in self.groups:
            filename_str = ', '.join(as_strs(group._paths))
            formatted_filenames += ''.join(['[', filename_str, ']\n'])
        return formatted_filenames


class DuplicateFilegroupError(MatchedPrefixGroupValidationError):
    """Exception raised when a `BiofileGroup` is matched with itself"""
    pass


class PrefixMatchError(MatchedPrefixGroupValidationError):
    """Exception raised the files do not share a prefix"""
    pass


class GroupLengthError(MatchedPrefixGroupValidationError):
    """Exception raised `MatchedPrefixGroups`s don't have equal lengths"""
    pass
