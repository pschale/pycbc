import h5py
import argparse
import numpy
import math
from pycbc.events import newsnr

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--triggers-file', required=True,
                    help='HDF file containing the triggers from pycbc_inspiral.')
parser.add_argument('--output-blips', required=True,
                    help='Name for txt file where to write the list of blips.'
                         'If the file already exists, will append new blips at the end.')
parser.add_argument('--output-rejected', required=True,
                    help='Name for txt file where to write the list of rejected glitches.'
                         'If the file already exists, will append new glitches at the end.')
opts = parser.parse_args()

def group_glitches(sorted_list, indices):
    one_glitch = [ sorted_list[0] ]
    glitches = []
    for idx in range( 1, len( sorted_list) ):
        if idx in indices:
            glitches.append( one_glitch )
            one_glitch = [ sorted_list[idx] ]
        else:
            one_glitch.append( sorted_list[idx] )
            if idx == len( sorted_list ) - 1:
                glitches.append( one_glitch )
    return glitches

pycbcinspiral = opts.triggers_file
blips_list = opts.output_blips
rejected_list = opts.output_rejected

triggers = h5py.File(pycbcinspiral, 'r')

ifo = triggers.keys()[0]
times = triggers[ ifo + '/end_time' ]
snr = triggers[ ifo + '/snr' ]
chisq = triggers[ ifo + '/chisq' ]
DOF = 2 * numpy.array( triggers[ ifo + '/chisq_dof' ] ) - 2
red_chisq = chisq / DOF
new_snr = newsnr( snr, red_chisq )

# Times are not sorted: get the indices to sort the times and
# sort all the lists with it

t_indices = numpy.argsort(times)
sorted_times = numpy.sort(times)
sorted_snr = [snr[idx] for idx in t_indices]
sorted_chisq = [chisq[idx] for idx in t_indices]
sorted_redchisq = [red_chisq[idx] for idx in t_indices]
sorted_newsnr = [new_snr[idx] for idx in t_indices]

# Group the glitches based on the time separation between
# different groups of triggers
# I chose 0.3 after several tests because sometimes there are
# double blips very close to each other that would be missed
# later with the length restriction if grouped as the same glitch
st = sorted_times[0]
st_idx = 0
glitch = [st]
glitch_st = []
glitch_st_idx = []
all_glitches = []
for t_idx in range( 1, len(sorted_times) ):
    if sorted_times[ t_idx ] - sorted_times[ t_idx - 1 ] < 0.3:
        glitch.append( sorted_times[t_idx] )
        # If it is the last time, add the glitch to the list
        if t_idx == len(sorted_times) - 1:
            all_glitches.append( glitch )
            glitch_st.append( st )
            glitch_st_idx.append( st_idx )
    else:
        all_glitches.append( glitch )
        glitch_st.append( st )
        glitch_st_idx.append( st_idx )
        st = sorted_times[t_idx]
        st_idx = t_idx
        glitch = [ st ]
# And group the snr, chisq, red_chisq, and new_snr of the separate glitches
# using the indices obtained from the times
all_glitches_snr = group_glitches( sorted_snr, glitch_st_idx )
all_glitches_chisq = group_glitches( sorted_chisq, glitch_st_idx )
all_glitches_redchisq = group_glitches( sorted_redchisq, glitch_st_idx )
all_glitches_newsnr = group_glitches( sorted_newsnr, glitch_st_idx )

# Discard glitches that are shorter than 10 ms and longer than 0.1 s
# Even if the blip is 10 ms long, there might be less triggers
# After some tests, it turns out we miss some blips by setting a lower cutoff
# Therefore, set only upper cutoff. The newSNR cutoff later will get rid of
# blips shorter than 10 ms
blip_glitches = []
longer_glitches = []
for t_glitch in all_glitches:
    if ( t_glitch[-1] - t_glitch[0] ) <= 0.1:
        blip_glitches.append( t_glitch )
     # Sometimes a blip glitch has only one trigger!
    elif len( t_glitch ) == 1:
        blip_glitches.append( t_glitch )
    else:
        longer_glitches.append( t_glitch )
# Get the indices of the blips to get their snr, chisq, newsnr
blip_idx = [ all_glitches.index( blip ) for blip in blip_glitches ]
blip_glitches_snr = [ all_glitches_snr[idx] for idx in blip_idx ]
blip_glitches_chisq = [ all_glitches_chisq[idx] for idx in blip_idx ]
blip_glitches_redchisq = [ all_glitches_redchisq[idx] for idx in blip_idx ]
blip_glitches_newsnr = [ all_glitches_newsnr[idx] for idx in blip_idx ]
# Same for the longer glitches
longer_idx = [ all_glitches.index( longer ) for longer in longer_glitches ]
longer_glitches_snr = [ all_glitches_snr[idx] for idx in longer_idx ]
longer_glitches_chisq = [ all_glitches_chisq[idx] for idx in longer_idx ]
longer_glitches_redchisq = [ all_glitches_redchisq[idx] for idx in longer_idx ]
longer_glitches_newsnr = [ all_glitches_newsnr[idx] for idx in longer_idx ]

# List of blips (time, snr, newsnr, chisq, red_chisq)
# We choose the maximum newsnr of the blips
blips = []
for idx in range(len(blip_glitches)):
    max_idx = blip_glitches_newsnr[idx].index( max( blip_glitches_newsnr[idx] ) )
    blips.append( [ blip_glitches[idx][max_idx], blip_glitches_snr[idx][max_idx],
                    blip_glitches_newsnr[idx][max_idx], blip_glitches_chisq[idx][max_idx],
                    blip_glitches_redchisq[idx][max_idx] ] )
# List of longer glitches (time, snr, newsnr, chisq, red_chisq)
# We choose the maximum newsnr of the glitches
rejected = []
for idx in range(len(longer_glitches)):
    max_idx = longer_glitches_newsnr[idx].index( max( longer_glitches_newsnr[idx] ) )
    rejected.append( [ longer_glitches[idx][max_idx], longer_glitches_snr[idx][max_idx],
                    longer_glitches_newsnr[idx][max_idx], longer_glitches_chisq[idx][max_idx],
                    longer_glitches_redchisq[idx][max_idx] ] )
# Remove blips with SNR > 100, SNR < 7.5 and newSNR < 5
# The newSNR cutoff can be set in pycbc_inspiral, but I keep it here
# just in case
# Save the rejected glitches in a separate list
del_blips = []
for blip_idx in range( len( blips ) ):
    if blips[blip_idx][1] > 150:
        del_blips.append( blip_idx )
    elif blips[blip_idx][1] < 7.5:
        del_blips.append( blip_idx )
    elif blips[blip_idx][2] < 5:
        del_blips.append( blip_idx )
    # Remove also blips with chisq > 2500
    elif blips[blip_idx][3] > 2500:
        del_blips.append( blip_idx )
[ rejected.append( blips[del_idx] ) for del_idx in del_blips ]
del_blips.reverse()
for del_idx in del_blips:
    del blips[del_idx]

# Write blip glitches in output blips (time, snr, newsnr, chisq, red_chisq)
# and longer + rejected glitches in output rejected
with open(blips_list,'a') as output:
    [ output.write( '%f %f %f %f %f \n' % ( blips[idx][0], blips[idx][1], blips[idx][2], blips[idx][3], blips[idx][4] ) ) for idx in range( len( blips ) ) ]
with open(rejected_list,'a') as output:
    [ output.write( '%f %f %f %f %f \n' % ( rejected[idx][0], rejected[idx][1], rejected[idx][2], rejected[idx][3], rejected[idx][4] ) ) for idx in range( len( rejected ) ) ]
