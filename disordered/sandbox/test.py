import pickle

def main():
	"""docstring for main"""
	

if __name__ == '__main__':
	main()

d = { "abc" : [1, 2, 3], "qwerty" : [4,5,6] }
afile = open(r'playground', 'wb')
pickle.dump(d, afile)
afile.close()

#reload object from file
file2 = open(r'playground', 'rb')
new_d = pickle.load(file2)
file2.close()

#print dictionary object loaded from file
print new_d