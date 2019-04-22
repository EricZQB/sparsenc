from pysnc import *
import math
import os
import sys
import struct


def print_usage():
    print("Usage: ./ProgramName snum gfpower")

if len(sys.argv) != 3:
    print_usage()
    sys.exit(1)

code_t = BAND_SNC
dec_t = CBD_DECODER
snum = int(sys.argv[1])
GFpower = int(sys.argv[2])

# Generate random data for the test
pktsize = 105
data = os.urandom(snum*pktsize)
# print(data)
buf = create_string_buffer(data)        # Create a buffer contains random source data
buf_p = cast(buf, POINTER(c_ubyte))

cnum = 0     # if precoding, use cnum = LDPC_check_num(snum, 0.01)
binary_precode = 1
systematic = 0
sp = snc_parameters(snum*pktsize, pktsize, cnum, snum, snum, code_t, 
                        binary_precode, GFpower, systematic, -1)   # RLNC

sc = snc.snc_create_enc_context(buf_p, byref(sp))
sp.seed = (snc.snc_get_parameters(sc))[0].seed
relay = snc.snc_create_buffer(byref(sp), 128)
decoder = snc.snc_create_decoder(byref(sp), dec_t)  # Create decoder
while not snc.snc_decoder_finished(decoder):
    pkt_p = snc.snc_generate_packet(sc)
    # Emulate (de)-serialization
    pktstr = pkt_p.contents.serialize(sp.size_g, sp.size_p, sp.gfpower)
    snc.snc_free_packet(pkt_p)
    pkt_p2 = snc.snc_alloc_empty_packet(snc.snc_get_parameters(sc))
    pkt_p2.contents.deserialize(pktstr, sp.size_g, sp.size_p, sp.gfpower)
    snc.snc_process_packet(decoder, pkt_p2)
    snc.snc_free_packet(pkt_p2)
    # snc.snc_process_packet(decoder, pkt_p)
    # snc.snc_free_packet(pkt_p)

snc.print_code_summary(snc.snc_get_enc_context(decoder), 0, 0)
buf_r = snc.snc_recover_data(snc.snc_get_enc_context(decoder))      # Recover the decoded packets a buffer
data_r = bytes(cast(buf_r, POINTER(c_ubyte * pktsize * snum))[0])   # Covert the buffer data to python types
if (data == data_r):
    print("SUCCESS: decoded data is IDENTICAL to the source.")
else:
    print("ERROR: decoded data is NOT identical to the source!")
snc.snc_free_decoder(decoder)
snc.snc_free_enc_context(sc)
