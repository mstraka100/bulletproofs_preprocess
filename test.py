import struct

with open("./testfilepython", 'bw') as f:
	f.write(struct.pack('<L', 52))