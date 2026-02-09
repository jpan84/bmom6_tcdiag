import numpy as np
import uxarray as ux

FILI = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0014-12-26-00000.nc'
FILO = '/glade/derecho/scratch/jpan/jpan_tcfields1/b.e23.BMOM.ne120np4_sx0.66av1.aqua.production.250417_ctrl.cam.h1i.0014-12-26-00000.nc.div200.nc'
UGRD = '/glade/p/cesmdata/inputdata/share/scripgrids/ne120np4_pentagons_100310.nc'

def check_topology(uxds):
    grid = uxds.uxgrid
    n_face = grid.n_face
    
    print(f"--- Grid Topology Report ---")
    print(f"Total Faces: {n_face}")
    
    # 1. Check for 'Holes' in connectivity
    # edge_face_connectivity (n_edge, 2)
    # If the second column is -1, it's a boundary edge.
    efc = grid.edge_face_connectivity
    boundary_edges = np.sum(efc[:, 1] == -1)
    print(f"Boundary Edges (edges with only 1 face): {boundary_edges}")
    
    # 2. Check Face-Node Connectivity
    # If this has fill values, the Green-Gauss loop skips those nodes.
    fnc = grid.face_node_connectivity
    missing_nodes = np.sum(fnc == -1)
    print(f"Missing Nodes in face-node map: {missing_nodes}")

    # 3. Simulate the _check_face_on_boundary logic
    # A face is a "boundary" if any of its edges are boundary edges.
    # We find edges for each face, then check if those edges are boundaries.
    fec = grid.face_edge_connectivity
    # Create a mask of boundary edges
    is_boundary_edge = (efc[:, 1] == -1)
    
    # Map that mask back to faces
    # If any edge connected to a face is a boundary, the face is a boundary.
    face_is_boundary = np.any(is_boundary_edge[fec], axis=1)
    n_boundary_faces = np.sum(face_is_boundary)
    
    print(f"Faces marked as 'Boundary': {n_boundary_faces}")
    print(f"Faces marked as 'Internal': {n_face - n_boundary_faces}")
    
    if n_boundary_faces > (n_face * 0.9):
        print("\nRESULT: Severe Topology Disconnection.")
        print("The gradient function will skip these faces and return NaN.")
    
    # 4. Check Coordinates (The 'Denominator' Check)
    # The function needs 3D Cartesian coords
    if hasattr(grid, 'face_x'):
        print(f"\nFace X sample: {grid.face_x[:3].values}")
        if np.all(np.abs(grid.face_x.values) < 1e-10):
             print("WARNING: Cartesian coordinates appear to be zero/uninitialized.")

# Usage
uxds = ux.open_dataset(UGRD, FILI)
check_topology(uxds)
