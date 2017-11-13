"""Classes for storing and validating various types of bioinformatics files.

Most common use case will be for pipeline scripts, which benefit from strict
validation of files between analysis steps.

Alternatively the classes can simply be used to validate groups of files for
function arguments etc. E.g if `BiofileGroup`[list_of_fasta_files, runs without
error, you have some confidence that these files can be passed to analysis
modules for processing without error. However file validation is targeted
mostly at catching user error - e.g passing a `Sam` file rather than a `Bam`
file, it does almost no internal checks for file format correctness.
Consequently, errors introduced to files during analysis will not be caught.

Usage
-----
Fastq(<fastq_file>)
Sam(<sam_file>)
BiofileGroup([file, ...], type=<type>)

Notes
-----
Files need not necessarily exist upon initialization, but they must exist upon
the first access attempt. In this way, a file group class can act as a promise
that files will exist later.

of the file before the first dot:
    <unique_ID>.<processing_step>.<...>.<extension>
E.g FA_SC.trim.map.sam
This could be changed with a little work, but it would probably be easier to
just rename the files.

FileGroups are designed to be initialized and then used without editing. Once
a group has been initialized, it is not recommended to attempt to change the
files it maps to.
"""
from pathlib import (
        Path,
        PurePath,
        )
from typing import (
        List,
        )

from fmbiopy.fmpaths import (
        add_suffix,
        is_empty,
        )


class Biofile():
    """Superclass for storing and validating bioinformatics files.

    Classes of more specific filetypes inherit the majority of their attributes
    and validation methods from it.


    Attributes
    ----------
    input_type : str
        (class) Filetype which is stored
    extensions : List[str]
        (class) List of acceptable input file extensions
    gzipped : bool
        Same as parameter
    """
    input_type: str = 'ANY'
    extensions: List[str] = ['ANY']

    def __init__(
            self,
            path: Path,
            gzipped: bool = False,
            possibly_empty: bool = False) -> None:
        """Initialize

        Parameters
        ----------
        path
            File to be stored in the group (need not exist)
        gzipped, Optional
            If True, file should be gzipped and have the .gz extension. Detecting
            that files are gzipped is obviously trivial, this flag is just used to
            ensure that the output function knows the gzip state. If this flag is
            unset but the files are gzipped, this indicates that an unintended zip
            step has occured or vice versa.
        possibly_empty, Optional
            If True, empty files will not cause validation errors.
        """

        # Store paramaters
        self._path = path
        self.gzipped = gzipped
        self.validated = False

        # True if it is acceptable for the files to be empty
        self.possibly_empty = possibly_empty

        # Get the actual file extension
        self.extension = self._get_extension()

        self.validate()


    @property
    def path(self) -> Path:
        """The path to the stored file"""
        if not self.validated:
            self.validate()
        return Path(self._path)

    @property
    def name(self) -> str:
        """The basename of the stored file"""
        if not self.validated:
            self.validate()
        return self._path.name

    def validate(self) -> bool:
        """Validation function to be used upon attempted access

        It is only called upon the first access attempt

        Returns
        -------
        bool
            True if successful
        """
        self._check_not_dir()
        self._check_file_not_empty()
        self._check_gzip()
        self._check_extension()
        self.validated = True
        return True

    def _check_not_dir(self) -> bool:
        """Check that input is not a directory

        Returns
        -------
        bool
            True if successful

        Raises
        ------
        TypeError
            If path is a directory
        """
        if Path(self._path).is_dir():
            raise TypeError('File cannot be a directory')
        return True

    def _get_extension(self) -> str:
        """Get the file extension

        Returns
        -------
        str
            If gzipped then it returns the two part extension E.g fq.gz.
            Otherwise just returns the final extension.
        """

        # If gzipped we return the two part extension E.g 'fq.gz'
        if self.gzipped:
            return ''.join(self._path.suffixes[-2:])
        return self._path.suffix


    def _check_file_not_empty(self)-> None:
        """Check that the file has contents"""
        if not self.possibly_empty:
            if is_empty(self._path):
                raise EmptyFileError(self._path)

    def _check_extension(self) -> bool:
        """Check that the file extension matches the accepted extensions

        The superclass should not have accepted extensions, but subclasses will
        use this method for extension validation.
        """
        if self.extensions != ['ANY']:
            # Extension check is not caps sensitive
            if self.extension.lower() not in self.extensions:
                raise FileExtensionError(self._path)
        return True

    def _check_gzip(self) -> bool:
        """Check that the gzipped flag parameter is correct"""

        if self.gzipped:
            if '.gz' not in self.extension:
                raise GzipStatusError(self._path)
        else:
            if 'gz' in self.extension:
                raise GzipStatusError(self._path)
        return True

    def __eq__(self, other) -> bool:
        """Test for BiofileGroup is equal to another"""
        return self.path == other.path


class Fastq(Biofile):
    """Biofile class for holding .fastq files."""

    input_type = 'fastq'
    extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']


class FwdFastq(Fastq):
    """Biofile class for holding forward pairs in paired .fastqs"""
    pass


class RevFastq(Fastq):
    """Biofile class for holding reverse pairs in paired .fastqs"""
    pass


class UnpairedFastq(Fastq):
    """Biofile class for holding unpaired .fastq files"""
    pass


class Fasta(Biofile):
    """Biofile for class holding .fasta files."""
    input_type = 'fasta'
    extensions = ['.fasta', '.fa', 'mfa', 'fna']


class SamtoolsFAIndex(Biofile):
    """Biofile class for holding samtools .fai files"""
    input_type = 'fai'
    extensions = ['.fai']


class Sam(Biofile):
    """Biofile class for holding .sam files"""

    input_type = 'sam'
    extensions = ['.sam']


class Bam(Biofile):
    """Biofile class for holding .bam files"""

    input_type = 'bam'
    extensions = ['.bam']


class Gzipped(Biofile):
    """A gzipped file"""
    input_type = 'gz'
    extensions = ['.gz']

    def __init__(self, *args, **kwargs):
        """Initialize"""
        super().__init__(*args, **kwargs, gzipped=True)

    def _check_gzip(self)-> bool:
        pass

    def _check_extension(self) -> bool:
        """Check that the file extension matches the accepted extensions

        The superclass should not have accepted extensions, but subclasses will
        use this method for extension validation.
        """
        if '.gz' not in self.extension.lower():
            raise FileExtensionError(self._path)
        return True


class Unzipped(Biofile):
    """Any file which is not gzipped"""
    input_type = 'ANY'
    extensions = ['ANY']

    def __init__(self, *args, **kwargs):
        """Initialize"""
        super().__init__(*args, **kwargs, gzipped=False)


class CentrifugeDB(Biofile):
    """A centrifuge database file"""
    input_type = 'centrifugedb'
    extensions = ['ANY']

    def __init__(self, *args, **kwargs)-> None:
        """Initialize"""
        # A CentrifugeDB biofile is in fact a reference to multiple files. We
        # store them here
        try:
            path = args[0]
        except KeyError:
            path = kwargs['path']

        self._idx: List[Path] = []
        for i in range(1, 4):
            suffix = ''.join(['.', str(i), '.cf'])
            self._idx.append(add_suffix(path, suffix))
        super().__init__(*args, **kwargs)


    def _check_file_not_empty(self):
        """Check that the file has contents"""
        if not self.possibly_empty:
            for path in self._idx:
                if is_empty(path):
                    raise EmptyFileError(self._path)


class TSV(Biofile):
    """A .tsv file"""
    input_type = 'tsv'
    extensions = ['.tsv']


class CentrifugeOutput(TSV):
    """A centrifuge output file"""
    pass


class Hits(TSV):
    """A hits file compatible with blobtools"""
    pass


class Txt(Biofile):
    """A .txt file"""
    input_type = 'txt'
    extensions = ['.txt']


class Html(Biofile):
    """A html file"""
    input_type = 'html'
    extensions = ['.html']


class FastQCReport(Html):
    """A FastQC report"""
    pass


class Zip(Biofile):
    """A zip file"""
    input_type = 'zip'
    extensions = ['.zip']


class Hist(Txt):
    """A histogram text file"""
    pass


class Adapters(Fasta):
    """A fasta file with adapter sequences"""
    pass


# -----------------------------------------------------------------------------
# Exceptions
# -----------------------------------------------------------------------------


class BiofileValidationError(Exception):
    """Basic exception for errors raised by Biofile validation checks

    Parameters
    ----------
    name
        The filename which caused the error

    Attributes
    ----------
    msg
        The formatted error message
    """
    def __init__(self, name: PurePath = None)-> None:
        self.name = name
        self.msg = self._construct_msg()
        super().__init__(self.msg)

    def _formatted_filename(self) -> str:
        """Construct the part of the error message which lists the file"""
        return '\n'.join(['File:', '\t' + str(self.name)])

    def _err_description(self) -> str:
        """Construct the error description to output"""
        return ''

    def _construct_msg(self) -> str:
        """Construct the combined error message"""
        return '\n'.join([
            self._formatted_filename(), self._err_description()])


class EmptyFileError(BiofileValidationError):
    """Exception raised when one or more of the stored files are empty"""

    def _err_description(self) -> str:
        return "File is empty but `possibly_empty` is False"


class GzipStatusError(BiofileValidationError):
    """Exception raised when one or more of the stored files are empty"""

    def _err_description(self) -> str:
        return "Gzip status of files does not match the gzip argument"


class FileExtensionError(BiofileValidationError):
    """Exception raised when files do not have the expected file extension"""

    def _err_description(self) -> str:
        return "Unexpected file extension"
