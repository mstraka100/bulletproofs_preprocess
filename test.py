import struct

with open("./testfilepython", 'bw') as f:
	f.write(struct.pack('<LLQ', 1, 0, 1363))