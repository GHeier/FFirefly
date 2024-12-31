import numpy as np
from scipy.interpolate import LinearNDInterpolator


class Field:
    def __init__(self, points, values, dimension=3, is_complex=False, is_vector=False):
        self.points = np.asarray(points)
        self.dimension = dimension
        self.is_complex = is_complex
        self.is_vector = is_vector

        if is_complex:
            self.values = np.asarray(values, dtype=np.complex128)
            self.interpolator_real = LinearNDInterpolator(self.points, self.values.real)
            self.interpolator_imag = LinearNDInterpolator(self.points, self.values.imag)
        else:
            self.values = np.asarray(values)
            self.interpolator = LinearNDInterpolator(self.points, self.values)

    def __call__(self, *coords):
        query_point = np.array(coords).reshape(1, -1)  # Ensure correct shape

        if self.is_complex:
            real_part = self.interpolator_real(query_point)
            imag_part = self.interpolator_imag(query_point)
            result = real_part + 1j * imag_part
        else:
            result = self.interpolator(query_point)

        if self.is_vector:
            # Return as a flat array if it's a vector field
            return result.reshape(-1)

        if result.size == 1:
            return result.item()  # Return scalar if single value

        return result

def load_field_from_file(filename, dimension=3, is_complex=False, is_vector=False):
    data = np.loadtxt(filename)
    points = data[:, :dimension]
    values = data[:, dimension]
    return Field(points, values, dimension, is_complex, is_vector)

class test:
    def __init__(self):
        pass
        #self.field = load_field_from_file('data.txt', dimension=3, is_complex=True, is_vector=False)

