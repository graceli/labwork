# Yield demo
def get():
	"""docstring for get"""
	print "call started"
	array = []
	for i in range(0,100):
		if i%10 == 0:
			yield array
			print "after yield"
			print "length: ", len(array)
			array = []
		else:
			array.append(i)
	
	print "call ended"

for a in get():
	print a
