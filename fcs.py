#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

# Define a spot size (in nanometers)
spot_radius = 2
spotX = 2.5
spotY = 2.5
spotZ = 2.5

# Base intensity of a particle
#   TODO: This is poorly defined here. This is just a scaling constant for the
#   intensity of a particle that's being scaled based on its position relative
#   to the beam.
INTENSITY = 1

# Gaussian beam profile parameters

#   Radial and axial std. dev.s of the Gaussian beam profile
w_xy = 1
w_z = 1
k = w_z/w_xy

# Absorption probability
epsilon = .5

# Step size when iterating through frames
STEP = 20

def main():
    print("Running")

    # Load trajectory file
    t = md.load("min.xtc", top="min.gro")
    print("Imported trajectories")

    # Select only phosphorous atoms in phosphate head group.
    #   This is a bit of a simplification, but it should significantly reduce the
    #   amount of atoms to iterate over if we're only considering the phosphorous
    #   at the center of the phosphate group.
    #   Error from this would be on the order of the bond lengths, so roughly
    #   1.5 angstrom.
    phosphorous_atoms = [a.index for a in t.topology.atoms if a.element.symbol == 'P']
    t.restrict_atoms(phosphorous_atoms)

    # This will be a list of detections per timestep
    detections = []

    # Iterate through each timestep (frame, in Gromacs terms)
    for frame_index in range(0,len(t), STEP):
        print("Processing frame %d out of %d" % (frame_index/STEP, len(t)/STEP))

        detected = 0

        # Iterate through each residue
        #   TODO: May want to track each residue independently. This can be
        #   accomplished by just adding another dimension to the detections array,
        #   and appending each residue to its own list within detections.
        for residue in t.topology.residues:

            # Do analysis if atom is in detection volume
            if not check_in_detection_volume(t, frame_index, residue):
                continue

            detected += generate_detection(t, frame_index, residue)
            print("\r%d detections on frame %d\r" % (detected, frame_index), end="")

        detections += [detected]
        print("")

    np.savetxt("dataout.dat", detections)

    return detections

#TODO: passing t is here bad practice, I think
# Returns true if a residue is within the detection volume, otherwise false.
def check_in_detection_volume(t, frame_index, residue):

    # For now, pay attention to all atoms, regardless of whether or not they're
    #   in the detection volume.
    #   TODO: Will this throw an unreachable exception?
    return True

    x, y, z = t.xyz[frame_index, residue._atoms[0].index]

    # Get magnitude of distance to the spot center
    distance = (x - spotX)**2 + (y - spotY)**2

    # Check if the distance is within the spot radius
    in_detection_area = distance <= spot_radius**2

    return in_detection_area

#TODO: passing t is here bad practice, I think
def generate_detection(t, frame_index, residue):

    # Get coordinates of residue (more correctly, of the P atom)
    x, y, z = t.xyz[frame_index, residue._atoms[0].index]

    # Calculate contribution to intensity from an atom, based on the Gaussian
    #   profile of the incident beam and the particle's position.
    intensity = \
        INTENSITY * np.exp(
        -( (x - spotX)**2 + (y - spotY)**2 + ((z - spotZ)/k)**2)
        / 2 * w_xy**2 )

    return intensity



if __name__ == "__main__":
    detections = main()
    plt.plot(detections)
    plt.show()
    print(detections)
