#!/bin/env python

# script to calculate the power spectrum of a snapshot stored in the cosmosim hdf5 format
# Takes command line options: 
#     input filenames 
#     particle type (not relevant here, leave as 0)
#     mesh grid size
#     folding parameter (max number of times to fold)
#     output file name
# ifold parameter should be a positive integer (preferably a power of 2) after nmesh and before outfile

from __future__ import print_function

import sys

from nbodykit import CurrentMPIComm
import nbodykit.source.catalog.array
import nbodykit.algorithms

import numpy as np
import h5py


comm = CurrentMPIComm.get()
comm_size = comm.Get_size()
comm_rank = comm.Get_rank()



print(comm_size,comm_rank)



def read_snapshot(fnames, part_type):
    """
    Read particle positions and masses from a simulation 
    snapshot. Distributes files across MPI tasks.
    """
    
    # Find number of files per snapshot
    if comm_rank == 0:
        print(fnames % 0)
        snap = h5py.File(fnames % 0,"r")
        nfiles = snap["Header"].attrs["NumFilesPerSnapshot"]
        boxsize = snap["Header"].attrs["BoxSize"]
        mass_value = snap["Header"].attrs["ParticleMass.Matter"]
        snap.close()
    else:
        nfiles    = None
        boxsize   = None
        mass_value = None
    nfiles    = comm.bcast(nfiles)
    boxsize   = comm.bcast(boxsize)
    mass_value = comm.bcast(mass_value)

    # Read the files
    pos  = []
    mass = []
    for ifile in range(nfiles):
        if ifile % comm_size == comm_rank:
            print("Task %d reading file %d" % (comm_rank, ifile))
            snap = h5py.File(fnames % ifile,"r")
            group = snap["Matter"]
            pos.append(group["Position"][...])
                # Make mass array using value from header
            mass.append(np.ones(pos[-1].shape[0], dtype=np.float64)*mass_value)
            snap.close()

    # Combine arrays from separate files
    if len(pos) > 0:
        pos = np.concatenate(pos)
    else:
        pos = None
    if len(mass) > 0:
        mass = np.concatenate(mass)
    else:
        mass = None

    return pos, mass, boxsize


def main():
    
    # Get command line parameters
    args = {}
    if comm_rank == 0:
        args["fnames"]    = sys.argv[1]
        args["part_type"] = int(sys.argv[2])
        args["nmesh"]     = int(sys.argv[3])
        args["ifold"]     = int(sys.argv[4])
        args["outfile"]   = sys.argv[5]
    args = comm.bcast(args)
    # Read the particle data
    pos, mass, boxsize = read_snapshot(args["fnames"], args["part_type"])

    fold_param = 1
    max_prev_k = 0.

    with open(args["outfile"],"w") as outf:
                outf.write("%16s %16s\n" % ("k", "P(k)"))

    while (fold_param<=args["ifold"]):
    
        #Fold the data into a smaller box
        pos = pos%(boxsize/fold_param)
        # Make an nbodykit Catalogue object with the positions
        cat = nbodykit.source.catalog.array.ArrayCatalog({"Position":pos, "Mass":mass})
    
        # Convert to mesh and calculate density field
        # Folded so boxsize is reduced
        mesh = cat.to_mesh(Nmesh=args["nmesh"], BoxSize=(boxsize/fold_param), compensated=True, position="Position", weight="Mass", window='tsc')

        # Calculate power spectrum
        if comm_rank == 0:
            print("Calculating power spectrum, folding parameter: %d" % fold_param)
        
        pk = nbodykit.algorithms.FFTPower(mesh, mode='1d').power
        
        max_k = np.max(pk["k"])
        
        # Write out the results
        if comm_rank == 0:
            print("Writing file: %s" % args["outfile"])
            if fold_param!=args["ifold"]:
                with open(args["outfile"],"a") as outf:
                    for k, p in zip(pk["k"], pk["power"].real):
                        if ((k > (max_prev_k*0.75)) & (k < (max_k*0.75))): 
                            # power spectrum has an inaccuracy of 1% above around 75% of the grid nyquist frequency so ignore these and use pk from more folded grid
                            outf.write("%16.8e %16.8e\n" % (k, p*fold_param**3)) #need to multiply power by fold_param^3 to rescale power
            else:
                with open(args["outfile"],"a") as outf:
                    for k, p in zip(pk["k"], pk["power"].real):
                        if (k > (max_prev_k*0.75)): 
                            # might as well include all possible pk from the most folded grid
                            outf.write("%16.8e %16.8e\n" % (k, p*fold_param**3)) #need to multiply power by fold_param^3 to rescale power
        fold_param *= 2
        max_prev_k = np.max(pk["k"])

if __name__ == "__main__":
    main()
