import h5py
import argparse
import numpy
import math
import pycbc.io
from pycbc.events import newsnr

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--triggers-file', required=True,
                    help='HDF file containing the triggers from pycbc_inspiral.')
parser.add_argument('--veto-file', required=True,
                    help='Veto definer file.')
parser.add_argument('--bank-file', required=True,
                    help='HDF template-bank file.')
parser.add_argument('--detector', required=True,
                    help='H1 or L1.')
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
vetoes = opts.veto_file
tmplt_bank = opts.bank_file
ifo = opts.detector
blips_list = opts.output_blips
rejected_list = opts.output_rejected

segment_name = 'CUMULATIVE_CAT_12H'
filter_func = None
triggers = pycbc.io.SingleDetTriggers(pycbcinspiral, tmplt_bank, vetoes,
                                    segment_name, filter_func, ifo)

times = triggers.end_time
snr = triggers.snr
red_chisq = triggers.rchisq
new_snr = triggers.newsnr

# Times are not sorted: get the indices to sort the times and
# sort all the lists with it

t_indices = numpy.argsort(times)
sorted_times = numpy.sort(times)
sorted_snr = [snr[idx] for idx in t_indices]
sorted_redchisq = [red_chisq[idx] for idx in t_indices]
sorted_newsnr = [new_snr[idx] for idx in t_indices]

# Group the glitches based on the time separation between
# different groups of triggers
st = sorted_times[0]
st_idx = 0
glitch = [st]
glitch_st = []
glitch_st_idx = []
glitches = []
for t_idx in range( 1, len(sorted_times) ):
    if sorted_times[ t_idx ] - sorted_times[ t_idx - 1 ] < 0.1:
        glitch.append( sorted_times[t_idx] )
        # If it is the last time, add the glitch to the list
        if t_idx == len(sorted_times) - 1:
            glitches.append( glitch )
            glitch_st.append( st )
            glitch_st_idx.append( st_idx )
    else:
        glitches.append( glitch )
        glitch_st.append( st )
        glitch_st_idx.append( st_idx )
        st = sorted_times[t_idx]
        st_idx = t_idx
        glitch = [ st ]
# And group the snr, chisq, red_chisq, and new_snr of the separate glitches
# using the indices obtained from the times
glitches_snr = group_glitches( sorted_snr, glitch_st_idx )
glitches_redchisq = group_glitches( sorted_redchisq, glitch_st_idx )
glitches_newsnr = group_glitches( sorted_newsnr, glitch_st_idx )

# List of blips (time, snr, newsnr, red_chisq)
# We choose the maximum newsnr of the blips
blips = []
rejected = []
for idx in range(len(glitches)):
    newsnr_idx = glitches_newsnr[idx].index(max(glitches_newsnr[idx]))
    # Reject glitches that have many triggers with newsnr smaller than 5
    if numpy.median(glitches_newsnr[idx]) < 5:
        rejected.append( [ glitches[idx][newsnr_idx], glitches_snr[idx][newsnr_idx],
                    glitches_newsnr[idx][newsnr_idx], glitches_redchisq[idx][newsnr_idx] ] )
    else:
        blips.append( [ glitches[idx][newsnr_idx], glitches_snr[idx][newsnr_idx],
                    glitches_newsnr[idx][newsnr_idx], glitches_redchisq[idx][newsnr_idx] ] )

# Remove blips with SNR > 100, SNR < 7.5 and newSNR < 6
# The newSNR cutoff can be set in pycbc_inspiral, but I keep it here
# just in case
# Save the rejected glitches in a separate list
del_blips = []
for blip_idx in range( len( blips ) ):
    if blips[blip_idx][1] > 150:
        del_blips.append( blip_idx )
    elif blips[blip_idx][1] < 7.5:
        del_blips.append( blip_idx )
    elif blips[blip_idx][2] < 6:
        del_blips.append( blip_idx )
    # Remove also blips with red_chisq > 200
    elif blips[blip_idx][3] > 200:
        del_blips.append( blip_idx )

[ rejected.append( blips[del_idx] ) for del_idx in del_blips ]
del_blips.reverse()
for del_idx in del_blips:
    del blips[del_idx]

blips.sort()
rejected.sort()
# Write blip glitches in output_blips (time, snr, newsnr, red_chisq)
# and rejected glitches in output_rejected for checking purposes
with open(blips_list,'a') as output:
    [ output.write( '%.6f %f %f %f \n' % ( blips[idx][0], blips[idx][1], blips[idx][2], blips[idx][3] ) ) for idx in range( len( blips ) ) ]
with open(rejected_list,'a') as output:
    [ output.write( '%.6f %f %f %f \n' % ( rejected[idx][0], rejected[idx][1], rejected[idx][2], rejected[idx][3] ) ) for idx in range( len( rejected ) ) ]
