"""
Vector class used for SSF and PSF Calculation

This module defines the Vector datatype that's used to 
store position (x,y,z) information for the atoms and residues
in an MD Simulation. Documentation is provided for reference.
"""

import numpy as np

class Vector:
    """
    Represents a 3D vector with x, y, and z components.

    This class defines a data type 'Vector' that represents a 3D vector with
    x, y, and z components, assuming the vector originates at the origin (0,0,0).
    Expands the numpy.array object.

    Attributes
    ----------
    x : float
        The x-component of the vector.
    y : float
        The y-component of the vector.
    z : float
        The z-component of the vector.
    
    """
    def __init__(self, x: float, y: float, z: float) -> None:
        """
        Initialize a Vector instance.

        Parameters
        ----------
        x : float
            The x-component of the vector.
        y : float
            The y-component of the vector.
        z : float
            The z-component of the vector.
        """
        self.x = x
        self.y = y
        self.z = z
        self.components = np.array([x, y, z])

    def __add__(self, other: 'Vector') -> 'Vector':
        """
        Add two vectors element-wise.

        Parameters
        ----------
        other : Vector
            Another Vector to add to this vector.

        Returns
        -------
        Vector
            A new Vector representing the element-wise sum of the two vectors.
        """
        return Vector(*(self.components + other.components))

    def __sub__(self, other: 'Vector') -> 'Vector':
        """
        Subtract another vector element-wise from this vector.

        Parameters
        ----------
        other : Vector
            Another Vector to subtract from this vector.

        Returns
        -------
        Vector
            A new Vector representing the element-wise difference of the two vectors.
        """
        return Vector(*(self.components - other.components))

    def __eq__(self, other: 'Vector') -> bool:
        """
        Check if two Vectors are equal.

        Parameters
        ----------
        other : Vector
            Vector to check equality to.

        Returns
        -------
        bool
            True if vectors are equal, False otherwise.
        """
        return np.array_equal(self.components, other.components)
        
    def calculate_cross_product(self, b: 'Vector') -> 'Vector':
        """
        Calculate the cross product of two vectors.

        Calculates the cross product of two vectors, which is the unit vector that is 
        perpendicular to both vectors.

        Parameters
        ----------
        b : Vector
            The second vector to calculate the cross product with.

        Returns
        -------
        Vector
            The resulting vector from the cross product.
        """
        c = np.cross(self.components, b.components)
        return Vector(*c)

    def calculate_dot_product(self, b: 'Vector') -> float:
        """
        Calculate the dot product of this vector with vector b.

        Parameters
        ----------
        b : Vector
            The second vector to calculate the dot product with.

        Returns
        -------
        float
            The dot product of the two vectors.
        """
        return float(np.dot(self.components, b.components))

    def magnitude(self) -> float:
        """
        Calculate the magnitude (length) of the vector.

        Returns
        -------
        float
            The magnitude (length) of the vector.
        """
        return np.linalg.norm(self.components)
        
    def calculate_projection(self, b: 'Vector') -> 'Vector':
        """
        Calculate the projection of this vector onto vector b.

        Parameters
        ----------
        b : Vector
            The vector to project onto.

        Returns
        -------
        Vector
            The resulting vector from the projection.
        """
        b_np = b.components
        projected_vector = (np.dot(self.components, b_np) / np.dot(b_np, b_np)) * b_np
        return Vector(*projected_vector)

    def __str__(self) -> str:
        """
        Redefine the string representation of the Vector.

        Redefines the output of print(Vector()) to display the x, y, z attributes.

        Returns
        -------
        str
            The string representation of the vector.
        """
        x, y, z = self.components
        return "[ " + str(x) + "\n  " + str(y) + "\n  " + str(z) + " ]"
    
    def scale(self, a: float) -> 'Vector':
        """
        Scale the vector by a scalar a.

        Parameters
        ----------
        a : float
            The scalar to scale the vector by.

        Returns
        -------
        Vector
            The scaled vector.
        """
        scaled_vector = a * self.components
        return Vector(*scaled_vector)
    
if __name__ == "__main__":
    assert (Vector(1,2,3) + Vector(3,2,1) == Vector(4,4,4))
    assert (Vector(1,2,3).y == 2)
    assert (Vector(1,2,3) - Vector(-1,0,4) == Vector(2,2,-1))
    assert (Vector(1,0,0).calculate_cross_product(Vector(0,1,0)) == Vector(0,0,1))
    assert (Vector(1,0,0).calculate_cross_product(Vector(0,0,0)) == Vector(0,0,0))
    assert (Vector(1,0,0).magnitude() == 1)
    assert (Vector(0,3,4).magnitude() == 5)
    assert (Vector(-1,-2,-2).magnitude() == 3)
    assert (Vector(3,1,0).calculate_projection(Vector(1,0,0)) == Vector(3,0,0))
    assert (Vector(1,2,3).calculate_projection(Vector(0,0,0)) == Vector(0,0,0))