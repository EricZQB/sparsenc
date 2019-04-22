from pysnc import *
import math
import os
import sys
import struct


def print_usage():
    print("Usage: ./ProgramName c_type d_type filename")
    print("        c_type   - code type, available options: RAND BAND")
    print("        d_type   - decoder type, available options: GG OA BD CBD")
    print("        filename - filename to be encoded and decoded.")

if len(sys.argv) != 4:
    print_usage()
    sys.exit(1)

c_type = 0
if sys.argv[1] == "RAND":
    c_type = RAND_SNC
elif sys.argv[1] == "BAND":
    c_type = BAND_SNC
else:
    print_usage()
    sys.exit(1)

d_type = 0
if sys.argv[2] == "GG":
    d_type = GG_DECODER
elif sys.argv[2] == "OA":
    d_type = OA_DECODER
elif sys.argv[2] == "BD":
    d_type = BD_DECODER
elif sys.argv[2] == "CBD":
    d_type = CBD_DECODER
else:
    print_usage()
    sys.exit(1)

filename = sys.argv[3]
if not os.path.isfile(filename):
    print("Filename is invalid.")
    sys.exit(1)

copyname = filename + '.dec.copy'
if os.path.isfile(copyname):
    print("%s is already exist, delete it first" % (copyname))
    os.remove(copyname)

#filename = filename.encode('UTF-8') # byte string literal for compatibility in Python 3.x
statinfo = os.stat(filename)
filesize = statinfo.st_size
datasize = 5120

chunks = math.ceil(float(filesize) / datasize)
offset = 0
remaining = filesize
while chunks > 0:
    print ("%d segments remain " % (chunks))
    if remaining > datasize:
        toEncode = datasize
    else:
        toEncode = remaining
    cnum = LDPC_check_num(source_num(toEncode, 1280), 0.01)
    sp = snc_parameters(toEncode, 1280, cnum, 32, 64, c_type, 1, 1, 0, -1)
    sc = snc.snc_create_enc_context(None, byref(sp))
    snc.snc_load_file_to_context(c_char_p(filename.encode()),
                                 offset, sc)  # Load file to snc_context
    sp.seed = (snc.snc_get_parameters(sc))[0].seed
    buf = snc.snc_create_buffer(byref(sp), 128)
    decoder = snc.snc_create_decoder(byref(sp), d_type)  # Create decoder
    while not snc.snc_decoder_finished(decoder):
        pkt_p = snc.snc_generate_packet(sc)
        # Emulate (de)-serialization
        pktstr = pkt_p.contents.serialize(sp.size_g, sp.size_p, sp.bnc)
        snc.snc_free_packet(pkt_p)
        pkt_p2 = snc.snc_alloc_empty_packet(snc.snc_get_parameters(sc))
        pkt_p2.contents.deserialize(pktstr, sp.size_g, sp.size_p, sp.bnc)
        snc.snc_process_packet(decoder, pkt_p2)
        snc.snc_free_packet(pkt_p2)
        # snc.snc_process_packet(decoder, pkt_p)
        # snc.snc_free_packet(pkt_p)

    snc.snc_recover_to_file(c_char_p(copyname.encode()),
                            snc.snc_get_enc_context(decoder))
    snc.snc_free_decoder(decoder)
    snc.snc_free_enc_context(sc)
    offset += toEncode
    remaining -= toEncode
    chunks -= 1

