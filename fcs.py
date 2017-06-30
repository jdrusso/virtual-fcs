#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Define a spot size (in nanometers)
spot_radius = 2
spotX = 2.5
spotY = 2.5
spotZ = 2.5

# Step size when iterating through frames
STEP = 20

def main():
    print("Running")

    # Load trajectory file
    t = md.load("min.xtc", top="min.gro")
    print("Imported trajectories")

    # Strip out the SRD atoms
    non_srd_atoms = [a.index for a in t.topology.atoms if not a.name == 'SRD']
    t.restrict_atoms(non_srd_atoms)


    # TODO: A better algorithm is probably to sort by x,y,z

    # This will be a list of detections per timestep
    detections = []

    # More efficient algorithm:
    # For each timestep:
    #   Filter only atoms with x,y,z within the detection volume
    #   Filter only residues with ALL their atoms in the previous list
    # Using filter function: On each timestep, look at the list of residues. Filter out
    #   any residue that doesn't have ALL its atoms in the detection area.
    # May be better (but way more computation?) to use a method that tracks each
    #   residue wherever it is, so that photobleaching etc can be incorporated

    # Iterate through each timestep (frame, in Gromacs terms)
    for frame_index in range(0,len(t), STEP):
        print("Processing frame %d out of %d" % (frame_index/STEP, len(t)/STEP))

        detected = 0


        #Iterate through each residue
        for res_index in range(1, t[frame_index].topology.n_residues):

            # Generate a list of the indices of all atoms in this residue
            atom_indices = t[frame_index].topology.select("resSeq %d" % res_index)


            # Iterate through each of these atoms, to make sure the entire residue is
            #   within the detection area
            in_detection_area = True
            centers = []
            for atom_index in atom_indices:

                #Get the XYZ coordinates of this atom
                atom = t.xyz[frame_index, atom_index,:]

                #Determine whether the atom is within the detection area
                in_detection_area = check_in_detection_volume(atom)

                # If the atom is within the detection area, add its center
                if in_detection_area:
                    centers.append(atom)
                # If it isn't, don't bother with the rest of the atoms in the residue.
                else:
                    break

            # If the entire residue is not within the detection area,
            #   go to the next residue.
            if not in_detection_area:
                continue

            # Average centers of atoms in residue.
            #   axis = 0 keyword means it'll take a list of (x,y,z) tuples and return
            #   a tuple of (avg_x, avg_y, avg_z)
            #   TODO: May just want to use this to decide if a particle is in the
            #    detection area
            center = np.mean(centers, axis=0)


            # At this point, I have the index of a residue that is entirely within
            #   the detection volume.
            detected += generate_detection(center)
            print("\r%d detections on frame %d             " % (detected, frame_index), end="")

        detections += [detected]
        print("")

    np.savetxt("dataout.dat", detections)

    return detections


# Returns true if a given tuple of (x,y,z) coordinates are within the detection
#   volume, otherwise false.
def check_in_detection_volume(coords):

    # Unpack coords tuple into x, y, z
    x, y, z = coords

    # Get magnitude of distance to the spot center
    distance = (x - spotX)**2 + (y - spotY)**2

    # Check if the distance is within the spot radius
    in_detection_area = distance <= spot_radius**2

    return in_detection_area


def generate_detection(center):

    #Radial and axial std. dev.s of the Gaussian beam profile
    w_xy = 1
    w_z = 1

    k = w_z/w_xy

    # Absorption probability
    epsilon = .5

    # Determine the probability of a detection
    # See Dix et al. (2006), J. Phys. Chem.
    probability = \
        1 * np.exp(
        -( (center[0] - spotX)**2 + (center[1] - spotY)**2 + ((center[2] - spotZ)/k)**2)
        / 2 * w_xy**2 )

    num = np.random.random()

    # Detection!
    if num < probability:
        return 1

    # No detection
    else:
        return 0

if __name__ == "__main__":
    detections = main()
    plt.plot(detections)
    plt.show()
    print(detections)
