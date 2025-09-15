import ffirefly.config as cfg
import numpy as np
from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator

class Field:
    def __init__(self, points, values, dimension=3, is_complex=False, is_vector=False):
        print("Initializing field")
        self.points = np.asarray(points)
        self.dimension = dimension
        self.is_complex = is_complex
        self.is_vector = is_vector

        if is_complex:
            self.values = np.asarray(values, dtype=np.complex128)
            rvals = self.values.real
            ivals = self.values.imag
            
            self._setup_interpolators(rvals, ivals)
        else:
            self.values = np.asarray(values)
            self._setup_interpolators(self.values)

        self.grid_mins = np.array([axis[0] for axis in self.interpolator_real.grid])
        self.grid_maxs = np.array([axis[-1] for axis in self.interpolator_real.grid])
        self.margin = 0.01

        print("Field initialized")

    def _is_mesh(self):
        """
        Determine if the dataset is on a structured mesh.
        """
        unique_coords = [np.unique(self.points[:, i]) for i in range(self.points.shape[1])]
        meshgrid = np.array(np.meshgrid(*unique_coords, indexing='ij')).T.reshape(-1, self.points.shape[1])
        return True
        return len(meshgrid) == len(self.points) and np.allclose(
            np.sort(meshgrid, axis=0), np.sort(self.points, axis=0)
        )

    def _setup_interpolators(self, real_values, imag_values=None):
        """
        Setup interpolators based on mesh or scattered data.
        """
        if self._is_mesh():
            print("Dataset detected as a structured mesh. Using RegularGridInterpolator.")
            grid_coords = [np.unique(self.points[:, i]) for i in range(self.points.shape[1])]
            grid_shape = tuple(len(coord) for coord in grid_coords)
            reshaped_real = real_values.reshape(grid_shape)
            self.interpolator_real = RegularGridInterpolator(grid_coords, reshaped_real, method='linear', bounds_error=False, fill_value=None)
            
            if imag_values is not None:
                reshaped_imag = imag_values.reshape(grid_shape)
                self.interpolator_imag = RegularGridInterpolator(grid_coords, reshaped_imag, method='linear', bounds_error=False, fill_value=None)
        else:
            print("Dataset detected as scattered. Using LinearNDInterpolator.")
            self.interpolator_real = LinearNDInterpolator(self.points, real_values)
            if imag_values is not None:
                self.interpolator_imag = LinearNDInterpolator(self.points, imag_values)

    def __call__(self, *coords):
        query_point = np.array(coords).reshape(1, -1)
        if isinstance(self.interpolator_real, RegularGridInterpolator):
            query_point = self.margin_clamp(query_point)
        
        real_part = self.interpolator_real(query_point)
        if self.is_complex:
            imag_part = self.interpolator_imag(query_point)
            result = (real_part + 1j * imag_part) if not self.is_vector else np.concatenate((real_part, imag_part), axis=None)
        else:
            result = real_part.reshape(-1) if self.is_vector else real_part.item()
        
        return result

    def margin_clamp(self, query_point):
        if np.any(query_point < (self.grid_mins - self.margin)) or np.any(query_point > (self.grid_maxs + self.margin)):
            print("\nQuery point:", query_point)
            print("Grid mins:", self.grid_mins)
            print("Grid maxs:", self.grid_maxs)
            raise ValueError("One or more query points exceed the allowed margin of " + str(self.margin))
        
        # Clamp query points to valid bounds with margin
        return np.clip(query_point, self.grid_mins - self.margin, self.grid_maxs + self.margin)



# Utility Functions

def load_field_from_file(filename, dimension=3, is_complex=False, is_vector=False):
    data = np.loadtxt(filename)
    points = data[:, :dimension]
    if not is_complex and not is_vector:
        values = data[:, dimension]
    elif is_complex and not is_vector:
        values = data[:, dimension] + 1j * data[:, dimension + 1]
    elif not is_complex and is_vector:
        values = data[:, dimension:]
    else:
        values = data[:, dimension:2*dimension] + 1j * data[:, 2*dimension:]

    return Field(points, values, dimension, is_complex, is_vector)


def empty_init():
    return Field([], [])


def load_field_from_data(points, values, dimension=3, is_complex=False, is_vector=False):
    return Field(points, values, dimension, is_complex, is_vector)
