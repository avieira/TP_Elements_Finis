export ReferenceInterval, Point1D

"""
    Point1D

    Structure reprenseting a point in a one dimensional domain.
    Here: simply a float.
"""
const Point1D = Float64

"""
    ReferenceInterval

    Structure storing the 2 vertices composing the reference interval for the finite element computations.
    By default, the reference interval is [-1,1].

    # Arguments
    - a, b: Point1D - extremal vertices of the reference interval.
"""
struct ReferenceInterval
    a::Point1D
    b::Point1D
    ReferenceInterval() = new(-1.0, 1.0)
    ReferenceInterval(a::Point1D, b::Point1D) = new(a, b)
end

