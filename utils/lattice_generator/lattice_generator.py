import numpy as np

def detect_collision(rhombi_a, rhombi_b):

    '''
    (from wikipedia)
    The separating axis theorem (SAT) says that:
    Two convex objects do not overlap if there exists a line (called axis)
    onto which the two objects' projections do not overlap.

    The separating axis theorem can be applied for fast collision detection
    between polygon meshes.
    The set of separating axis is given by all face normals of both particles and all cross_products of those face normals.
    '''


    sep_axis = calculate_separating_axis(rhombus_a.facenormals, rhombus_b.facenormals)
