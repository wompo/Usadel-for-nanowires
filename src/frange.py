import sys

def frange(start, end=None, inc=None):
    "A range function, that does accept floats"

    if end == None:
        end = start + 0.0
        start = 0.0
    else: start += 0.0 # force it to be a float

    if inc == None:
        inc = 1.0

    count = int((end - start) / inc)
    if start + count * inc != end + inc:
        count += 1

    L = [None,] * count
    for i in xrange(count):
        L[i] = start + i * inc

    return L

start	= float(sys.argv[1])
end		= float(sys.argv[2])
inc		= float(sys.argv[3])

for item in frange(start, end, inc):
	print("%.2f" % item)
