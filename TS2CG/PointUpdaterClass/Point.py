from pathlib import Path
import numpy as np
import shutil
from datetime import datetime
import logging
from typing import Optional, Dict, Union, List

logger = logging.getLogger(__name__)

class Point:
    """
    A class representing a membrane structure with inclusions and exclusions.
    Loads and manages data from a point folder containing membrane definitions.

    Example:
        >>> membrane = Point("path/to/point/folder")
        >>> # Access membrane properties
        >>> print(membrane.inner.coordinates)
        >>> print(membrane.outer.mean_curvature)
        >>> # Work with inclusions and exclusions
        >>> print(membrane.inclusions.count)
        >>> membrane.inclusions.add_protein(type_id=1, point_id=42)
        >>> membrane.exclusions.add_point(point_id=50, radius=2.0)
    """

    class Membrane:
        """Represents a membrane layer with associated properties."""

        def __init__(self, data: np.ndarray):
            """Initialize membrane layer from raw data."""
            self._process_data(data)

        def _process_data(self, data: np.ndarray):
            """Convert raw data into accessible properties."""
            # Basic properties
            self.ids = data[0].astype(int)
            self.domain_ids = data[1].astype(int)
            self.area = data[2]

            # Coordinates and vectors
            self.coordinates = data[3:6].T  # X, Y, Z
            self.normals = data[6:9].T      # Nx, Ny, Nz
            self.principal_vectors = {
                'p1': data[9:12].T,   # P1x, P1y, P1z
                'p2': data[12:15].T    # P2x, P2y, P2z
            }
            self.curvature = {
                'c1': data[15],        # First principal curvature
                'c2': data[16]         # Second principal curvature
            }

        @property
        def mean_curvature(self) -> np.ndarray:
            """Calculate mean curvature for all points."""
            return (self.curvature['c1'] + self.curvature['c2']) / 2

        @property
        def gaussian_curvature(self) -> np.ndarray:
            """Calculate Gaussian curvature for all points."""
            return self.curvature['c1'] * self.curvature['c2']

        def get_points_by_domain(self, domain_id: int) -> np.ndarray:
            """Get coordinates of all points in a specific domain."""
            mask = self.domain_ids == domain_id
            return self.coordinates[mask]

    class Inclusion:
        """Manages protein inclusions in the membrane."""

        def __init__(self, data: Optional[np.ndarray] = None):
            self.points = []
            if data is not None:
                self._process_data(data)

        def _process_data(self, data: np.ndarray):
            """Process inclusion data."""
            if data is None or len(data) == 0:
                return

            for i in range(data.shape[1]):
                point = {
                    'id': int(data[0, i]),
                    'type_id': int(data[1, i]),
                    'point_id': int(data[2, i]),
                    'orientation': data[3:6, i]
                }
                self.points.append(point)

        @property
        def count(self) -> int:
            """Number of protein inclusions."""
            return len(self.points)

        def get_all(self) -> List[int]:
            """Get all points with protein inclusions."""
            return [p['point_id'] for p in self.points]

        def get_by_type(self, type_id: int) -> List[dict]:
            """Get all inclusions of a specific type."""
            return [p for p in self.points if p['type_id'] == type_id]

        def add_protein(self, type_id: int, point_id: int,
                       orientation: Optional[np.ndarray] = None):
            """
            Add a protein inclusion.

            Args:
                type_id: Type identifier for the protein
                point_id: Point ID where protein should be placed
                orientation: Vector specifying protein orientation
            """
            if orientation is None:
                orientation = np.array([0, 0, 1])

            point = {
                'id': len(self.points),
                'type_id': type_id,
                'point_id': point_id,
                'orientation': orientation
            }
            self.points.append(point)

    class Exclusion:
        """Manages lipid exclusions in the membrane."""

        def __init__(self, data: Optional[np.ndarray] = None):
            self.points = []
            if data is not None:
                self._process_data(data)

        def _process_data(self, data: np.ndarray):
            """Process exclusion data."""
            if data is None or len(data) == 0:
                return

            for i in range(data.shape[1]):
                point = {
                    'id': int(data[0, i]),
                    'type_id': int(data[1, i]),
                    'radius': float(data[2, i])
                }
                self.points.append(point)

        @property
        def count(self) -> int:
            """Number of exclusion points."""
            return len(self.points)

        def get_all(self) -> List[int]:
            """Get all excluded points."""
            return [p['point_id'] for p in self.points]

        def add_pore(self, point_id: int, radius: float = 1.0):
            """
            Add a pore in the lipid membrane.

            Args:
                point_id: Point ID where lipids should be excluded
                radius: Radius of exclusion zone
            """
            point = {
                'id': len(self.points),
                'type_id': 0,
                'point_id': point_id,
                'radius': radius
            }
            self.points.append(point)

    def __init__(self, path: Union[str, Path]):
        """
        Initialize point class from a point folder.

        Args:
            path: Path to the point folder
        """
        self.path = Path(path)
        if not self.path.exists():
            raise FileNotFoundError(f"Point folder not found: {self.path}")

        self._load_data()

    def _load_data(self):
        """Load all data from the point folder."""
        try:
            # Load membrane data
            inner_data = self._load_membrane_file("InnerBM.dat")
            outer_data = self._load_membrane_file("OuterBM.dat")
            self.inner = self.Membrane(inner_data)
            self.outer = self.Membrane(outer_data)

            # Load inclusions and exclusions
            inc_data = self._load_modification_file("IncData.dat")
            exc_data = self._load_modification_file("ExcData.dat")
            self.inclusions = self.Inclusion(inc_data)
            self.exclusions = self.Exclusion(exc_data)

        except Exception as e:
            logger.error("Failed to load membrane data", exc_info=True)
            raise

    def _load_membrane_file(self, filename: str) -> np.ndarray:
        """Load membrane definition file."""
        filepath = self.path / filename

        # Read the first few lines to check for Box information
        with open(filepath) as f:
            first_lines = [next(f) for _ in range(4)]

        # Store box dimensions if this is OuterBM.dat
        if "OuterBM" in filename:
            self.box = self._parse_box_line(first_lines[0])
            skiprows = 4  # Skip box line and other headers
        else:
            skiprows = 3  # Skip only standard headers

        return np.loadtxt(filepath, skiprows=skiprows).T

    def _parse_box_line(self, line: str) -> tuple:
        """Parse box dimensions from header line."""
        parts = line.split()
        return [float(x) for x in parts[1:4]]

    def _load_modification_file(self, filename: str) -> Optional[np.ndarray]:
        """Load modification (inclusion/exclusion) file."""
        try:
            return np.loadtxt(self.path / filename, skiprows=2).T
        except (ValueError, FileNotFoundError):
            return None

    def save(self, backup: bool = True):
        """
        Save membrane structure back to files.

        Args:
            backup: Whether to create a backup before saving
        """
        if backup:
            self._create_backup()

        self._save_membranes()
        self._save_modifications()

    def _create_backup(self):
        """Create a backup of the point folder."""
        backup_path = self.path.parent / f"#{self.path.name}#"
        shutil.copytree(self.path, backup_path)
        logger.info(f"Created backup at: {backup_path}")

    def _save_membranes(self):
        """Save membrane data to files."""
        for name, membrane in [("InnerBM.dat", self.inner), ("OuterBM.dat", self.outer)]:
            data = np.zeros((17, len(membrane.ids)))  # 17 columns as per format
            data[0] = membrane.ids
            data[1] = membrane.domain_ids
            data[2] = membrane.area
            data[3:6] = membrane.coordinates.T
            data[6:9] = membrane.normals.T
            data[9:12] = membrane.principal_vectors['p1'].T
            data[12:15] = membrane.principal_vectors['p2'].T
            data[15] = membrane.curvature['c1']
            data[16] = membrane.curvature['c2']

            # Create header
            headers = []

            # Add box dimensions for OuterBM.dat
            if name == "OuterBM.dat":
                headers.append(f"Box     {self.box[0]:.3f}     {self.box[1]:.3f}     {self.box[2]:.3f}")

            headers.extend([
                f"< Point NoPoints     {len(membrane.ids)}>",
                "< id domain_id area X Y Z Nx Ny Nz P1x P1y P1z P2x P2y P2z C1 C2  >",
                f"< {'Outer' if 'Outer' in name else 'Inner'} >"
            ])

            header = '\n'.join(headers)

            np.savetxt(
                self.path / name,
                data.T,
                header=header,
                comments='',
                fmt = ['%10d', '%4d', '%9.3f'] + ['%9.3f']*3 + ['%7.3f']*11
            )

    def _save_modifications(self):
        """Save inclusion and exclusion data to files."""
        # Save inclusions
        if self.inclusions.points:
            data = np.zeros((6, len(self.inclusions.points)))
            for i, point in enumerate(self.inclusions.points):
                data[0:3, i] = [point['id'], point['type_id'], point['point_id']]
                data[3:6, i] = point['orientation']

            header = f"< Inclusion NoInc {len(self.inclusions.points)} >\n"
            header += "< id typeid pointid lx ly lz >"

            np.savetxt(
                self.path / "IncData.dat",
                data.T,
                header=header,
                comments='',
                fmt=['%12d', '%12d', '%12d', '%8.3f', '%8.3f', '%8.3f']
            )

        # Save exclusions
        if self.exclusions.points:
            data = np.zeros((3, len(self.exclusions.points)))
            for i, point in enumerate(self.exclusions.points):
                data[0:3, i] = [point['id'], point['type_id'], point['radius']]
                print([point['id'], point['type_id'], point['radius']])

            header = f"< Exclusion NoExc {len(self.exclusions.points)} >\n"
            header += "< id typeid radius >"

            np.savetxt(
                self.path / "ExcData.dat",
                data.T,
                header=header,
                comments='',
                fmt=['%12d', '%12d', '%12d']
            )

    def update_domains(self, layer: str = "both", domain_ids: Optional[np.ndarray] = None):
        """
        Update domain assignments for membrane layers.

        Args:
            layer: Which layer to update ("both", "inner", or "outer")
            domain_ids: New domain assignments
        """
        if domain_ids is None:
            return

        if layer == "both":
            self.inner.domain_ids = domain_ids
            self.outer.domain_ids = domain_ids
        elif layer == "inner":
            self.inner.domain_ids = domain_ids
        elif layer == "outer":
            self.outer.domain_ids = domain_ids
        else:
            raise ValueError(f"Invalid layer: {layer}")
