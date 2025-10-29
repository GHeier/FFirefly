"""
Data loaders for FFirefly plotting package.
Handles loading different field types and text files.
"""

from typing import Union, Tuple, Any, Optional
from pathlib import Path
import numpy as np
import pandas as pd
from dataclasses import dataclass

try:
    import firefly as fly
except ImportError:
    fly = None


@dataclass
class FieldData:
    """Container for field data and metadata."""

    data: Any  # Field object or raw data
    field_type: str  # 'Field_R', 'Field_C', 'Field_RM', 'Field_CM', 'text', 'csv'
    dimension: int
    mesh: list
    domain: np.ndarray
    w_points: np.ndarray
    is_complex: bool
    is_matrix: bool
    filename: str

    def __call__(self, *args, **kwargs):
        """Allow calling the field directly."""
        return self.data(*args, **kwargs)

    @property
    def has_frequency(self) -> bool:
        """Check if field has frequency data."""
        return len(self.w_points) > 0


def detect_file_type(filename: str) -> str:
    """Detect file type from extension and content.

    Args:
        filename: Path to file

    Returns:
        File type string: 'h5', 'hdf5', 'dat', 'txt', 'csv'
    """
    path = Path(filename)
    ext = path.suffix.lower()

    if ext in ['.h5', '.hdf5']:
        return 'h5'
    elif ext in ['.dat', '.txt']:
        return 'dat'
    elif ext == '.csv':
        return 'csv'
    else:
        # Try to detect from content
        try:
            with open(filename, 'r') as f:
                first_line = f.readline()
                if ',' in first_line:
                    return 'csv'
                else:
                    return 'dat'
        except:
            return 'unknown'


def load_text_file(filename: str) -> FieldData:
    """Load data from text file (dat, txt, csv).

    Args:
        filename: Path to text file

    Returns:
        FieldData object with loaded data
    """
    # Detect if file has header
    with open(filename, 'r') as f:
        first_line = f.readline().strip()
        try:
            float(first_line.split()[0])
            header = None
        except (ValueError, IndexError):
            header = 0

    # Read file
    df = pd.read_csv(filename, sep=r'\s+', engine='python', header=header)

    # Remove comment column if present
    columns = df.columns.tolist()
    if '#' in columns:
        columns.remove('#')
        df = df[columns]

    # Convert to numpy arrays
    data = df.to_numpy()

    return FieldData(
        data=data,
        field_type='text',
        dimension=data.shape[1] - 1,  # Assume last column is values
        mesh=[],
        domain=np.array([]),
        w_points=np.array([]),
        is_complex=False,
        is_matrix=False,
        filename=filename
    )


def load_field_h5(filename: str) -> FieldData:
    """Load field from HDF5 file.

    Args:
        filename: Path to HDF5 file

    Returns:
        FieldData object with loaded field

    Raises:
        RuntimeError: If firefly module not available
        RuntimeError: If unable to determine field type
    """
    if fly is None:
        raise RuntimeError("firefly module not available. Cannot load HDF5 files.")

    # Try to load as different field types
    field = None
    field_type = None
    is_complex = False
    is_matrix = False

    # Try Field_R (real scalar)
    try:
        field = fly.Field_R(filename)
        field_type = 'Field_R'
        is_complex = False
        is_matrix = False
    except:
        pass

    # Try Field_C (complex scalar)
    if field is None:
        try:
            field = fly.Field_C(filename)
            field_type = 'Field_C'
            is_complex = True
            is_matrix = False
        except:
            pass

    # Try Field_RM (real matrix)
    if field is None:
        try:
            field = fly.Field_RM(filename)
            field_type = 'Field_RM'
            is_complex = False
            is_matrix = True
        except:
            pass

    # Try Field_CM (complex matrix)
    if field is None:
        try:
            field = fly.Field_CM(filename)
            field_type = 'Field_CM'
            is_complex = True
            is_matrix = True
        except:
            pass

    if field is None:
        raise RuntimeError(f"Unable to load {filename} as any field type")

    return FieldData(
        data=field,
        field_type=field_type,
        dimension=field.dimension,
        mesh=field.mesh,
        domain=np.array(field.domain) if field.domain else np.array([]),
        w_points=field.w_points if hasattr(field, 'w_points') else np.array([]),
        is_complex=is_complex,
        is_matrix=is_matrix,
        filename=filename
    )


def load(filename: str, force_type: Optional[str] = None) -> FieldData:
    """Load data from file with automatic type detection.

    Args:
        filename: Path to file
        force_type: Force specific field type ('Field_R', 'Field_C', etc.)

    Returns:
        FieldData object
    """
    file_type = detect_file_type(filename)

    if force_type is not None:
        # Force specific field type
        if fly is None:
            raise RuntimeError("firefly module not available")

        if force_type == 'Field_R':
            field = fly.Field_R(filename)
            return FieldData(
                data=field, field_type='Field_R', dimension=field.dimension,
                mesh=field.mesh, domain=np.array(field.domain),
                w_points=field.w_points, is_complex=False, is_matrix=False,
                filename=filename
            )
        elif force_type == 'Field_C':
            field = fly.Field_C(filename)
            return FieldData(
                data=field, field_type='Field_C', dimension=field.dimension,
                mesh=field.mesh, domain=np.array(field.domain),
                w_points=field.w_points, is_complex=True, is_matrix=False,
                filename=filename
            )
        elif force_type == 'Field_RM':
            field = fly.Field_RM(filename)
            return FieldData(
                data=field, field_type='Field_RM', dimension=field.dimension,
                mesh=field.mesh, domain=np.array(field.domain),
                w_points=field.w_points, is_complex=False, is_matrix=True,
                filename=filename
            )
        elif force_type == 'Field_CM':
            field = fly.Field_CM(filename)
            return FieldData(
                data=field, field_type='Field_CM', dimension=field.dimension,
                mesh=field.mesh, domain=np.array(field.domain),
                w_points=field.w_points, is_complex=True, is_matrix=True,
                filename=filename
            )
        else:
            raise ValueError(f"Unknown force_type: {force_type}")

    # Auto-detect
    if file_type == 'h5':
        return load_field_h5(filename)
    else:
        return load_text_file(filename)


def get_label_from_filename(filename: str) -> str:
    """Extract a clean label from filename.

    Args:
        filename: Path to file

    Returns:
        Clean label string
    """
    path = Path(filename)
    name = path.stem

    # Remove common extensions
    for ext in ['.dat', '.csv', '.txt', '.h5', '.hdf5']:
        name = name.replace(ext, '')

    # Replace underscores with spaces
    name = name.replace('_', ' ')

    return name
