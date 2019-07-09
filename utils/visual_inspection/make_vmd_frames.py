import numpy as np 
import glob


# Set lenghts
alpha=(60*np.pi)/180.0;
Lx=1.0;
Ly=Lx;
Lz=0.1*Lx;
h=Lx*np.sin(alpha);
a_x = Lx*np.cos(alpha);

rax = Lx;
ray = 0.0;
raz = 0.0;

rbx = a_x;
rby = h;
rbz = 0.0;

rcx = 0.0;
rcy = 0.0;
rcz = Lz;

dist_center_pframe = np.array([
    np.array([(- rax - rbx - rcx)/2,(- ray - rby - rcy)/2, (- raz - rbz - rcz)/2]),
    np.array([(  rax - rbx - rcx)/2,(  ray - rby - rcy)/2, (  raz - rbz - rcz)/2]),
    np.array([(  rax + rbx - rcx)/2,(  ray + rby - rcy)/2, (  raz + rbz - rcz)/2]),
    np.array([(- rax + rbx - rcx)/2,(- ray + rby - rcy)/2, (- raz + rbz - rcz)/2]),
    np.array([(- rax - rbx + rcx)/2,(- ray - rby + rcy)/2, (- raz - rbz + rcz)/2]),
    np.array([(+ rax - rbx + rcx)/2,(+ ray - rby + rcy)/2, (+ raz - rbz + rcz)/2]),
    np.array([(+ rax + rbx + rcx)/2,(+ ray + rby + rcy)/2, (+ raz + rbz + rcz)/2]),
    np.array([(- rax + rbx + rcx)/2,(- ray + rby + rcy)/2, (- raz + rbz + rcz)/2])])

def rotation(phi):
    Rot_mat = np.zeros((3,3))
    Rot_mat[0,0] = np.cos(phi)
    Rot_mat[0,1] = np.sin(phi)
    Rot_mat[0,2] = 0

    Rot_mat[1,0] = -np.sin(phi)
    Rot_mat[1,1] = np.cos(phi)
    Rot_mat[1,2] = 0

    Rot_mat[2,0] = 0
    Rot_mat[2,1] = 0
    Rot_mat[2,2] = 1

    return Rot_mat

def calculate_edge_points(pos_i, orient_i):
    dist_center_lframe = np.array([ np.dot(dist_center_pframe, rotation(phi)) for phi in orient_i])
    edge_pos = np.array([ pi + dist_i for pi, dist_i in zip(pos_i, dist_center_lframe)])
    edge_pos = np.reshape(edge_pos, (len(pos_i)*len(dist_center_pframe), 3))

    return edge_pos

checkpoints= glob.glob("Box*.bin")
check_point_values = np.sort([ int(point.split("_")[-1].split(".")[0]) for point in checkpoints ])

with open("pos.xyz", 'a') as f:
    for val in check_point_values: 
        pos_i = np.fromfile("positions_{}.bin".format(val))
        pos_i = np.reshape(pos_i, (-1,3))

        orient_i = np.fromfile("orientations_{}.bin".format(val))
        orient_i = np.reshape(orient_i, (-1,5))[:,4]

        N=len(pos_i)
        Nconst=1000
        f.write("{}\n".format(Nconst*8))
        f.write("Particles of frame {}\n".format(val))

        " returns a (N*N_edge_points,3) array "
        edge_pos = calculate_edge_points(pos_i, orient_i)
        rest_pos = np.zeros((8000-len(edge_pos),3))
        edge_pos = np.concatenate((edge_pos,rest_pos), axis=0)

        [ f.write("cb       {}   {}   {}\n".format(elem[0], elem[1], elem[2])) for elem in edge_pos]


