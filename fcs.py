#!/usr/bin/env python3
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt

######################## Detection area parameters ###########################
# Radius of detection area (in nanometers) - used for determining if a particle
#   is in the detection volume
spot_radius = 5

# Coordinates of detection area center (in nanometers)
spotX = 10.0
spotY = 10.0
spotZ = 10.0
##############################################################################

###################### Gaussian beam profile parameters ######################
#   Radial and axial std. dev.s of the Gaussian beam profile
w_xy = 2
w_z = 2
k = w_z/w_xy
##############################################################################

# Step size when iterating through data frames
STEP = 1

# Scaling constant for the intensity of a fluorescing particle
INTENSITY = 1000

# Percentage of tagged particles
SAMPLING_RATIO = .50

def main():

    # Gather list of intensities at each timestep for each lipid
    intensities = analyze_intensities()

    # Plot intensities and their autocorrelations
    plot(intensities)
    plt.show()

def plot(intensities):

    # Ignore data from particles that never entered the detection area
    # intensities = [x for x in intensities if np.mean(x) > 5]

    ################## Plot intensity traces #############################
    # Plot the intensities of detections associated with each lipid
    plt.subplot(221)
    for _data in intensities:
        plt.plot(np.arange(0,len(_data), 1), _data, marker = 'o', markersize=.5)
    plt.xlabel("Time")
    plt.ylabel("Intensity")
    plt.title("Intensity Trace")

    # Plot sum of intensities (i.e. what a detector would see)
    plt.subplot(222)
    summed_data = np.sum(intensities, axis=0)
    plt.plot(np.arange(0,len(summed_data), 1), summed_data, marker = 'o', markersize=.5)
    plt.xlabel("Time")
    plt.ylabel("Intensity")
    plt.title("Summed Intensity Trace")


    ################## Plot autocorrelations ##############################
    # Individual
    plt.subplot(223)
    for _data in intensities:
        plt.acorr(_data, maxlags=2000, usevlines=False, linestyle='-', marker="None")
    plt.xlabel("Time lag")
    plt.title("Autocorrelation")
    plt.xlim([0,2000])

    #Summed
    plt.subplot(224)
    plt.acorr(summed_data, maxlags=2000, usevlines=False, linestyle='-', marker="None")
    plt.xlabel("Time lag")
    plt.title("Autocorrelation of Summed Intensities")
    plt.xlim([0,2000])

def analyze_intensities():
    print("Beginning analaysis.")

    print("Importing trajectories")
    # Load trajectory file
    t = md.load("20nm/run.xtc", top="20nm/system.gro")
    print("Trajectories successfully imported.")

    print("Timestep for plotted data will be %f picoseconds" % (t.timestep * STEP))

    # Select only phosphorous atoms in phosphate head group.
    #   This is a bit of a simplification, but it should significantly reduce the
    #   amount of atoms to iterate over if we're only considering the phosphorous
    #   at the center of the phosphate group.
    #   Error from this would be on the order of the bond lengths, so roughly
    #   1.5 angstrom.
    print("Starting with %d atoms" % t.topology.n_atoms)

    phosphorous_atoms = [a.index for a in t.topology.atoms if a.element.symbol == 'P']
    t.atom_slice(phosphorous_atoms, inplace=True)

    print("Reduced to %d phosphorous atoms" % t.topology.n_atoms)

    # Reduce to the sampling ratio * number of phosphorous atoms
    num_sampled = int(t.topology.n_atoms * SAMPLING_RATIO)

    # Randomly select the sampled atoms
    sampled = np.random.choice([a.index for a in t.topology.atoms], num_sampled, replace=False)
    t.atom_slice(sampled, inplace=True)

    print("Reduced to %d \"tagged\" phosphorous atoms" % t.topology.n_atoms)

    # This will be a list of lists of detections for each residue per timestep
    #   Use a list comprehension here instead of [[]]*<some number> to deep copy the
    #   empty list.
    detections = [ [] for x in t.topology.residues ]

    # Iterate through each timestep (frame, in Gromacs terms)
    for frame_index in range(0,len(t), STEP):
        print("\rProcessing frame %d out of %d" % (frame_index/STEP, len(t)/STEP), end="\r")

        # Initialize this to -1, since it increments at the beginning of the residue loop
        idx = -1

        # Iterate through each residue
        #   TODO: May want to track each residue independently. This can be
        #   accomplished by just adding another dimension to the detections array,
        #   and appending each residue to its own list within detections.
        for residue in t.topology.residues:

            # Increment this here, so that it'll increment even if we skip the residue
            idx += 1

            # Do analysis if atom is in detection volume
            if not check_in_detection_volume(t, frame_index, residue):
                print("Not in detection volume, skipping.")
                continue

            detected = generate_detection(t, frame_index, residue)
            detections[idx].append(detected)

    # np.savetxt("dataout.dat", detections)

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
    main()
