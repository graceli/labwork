import result
import xvgfile
	
# This class appends new data into underlying tables

class Loader
	_xvgfile = xvgfile.XVGFile()
	_result = result.Result()

	def load(self, analysisName, TableStructure, analysisRoot):
		dirList = os.listdir(analysisRoot)
		for dir in dirList:
			filesList = os.listdir(dir)
			for file in filesList:
				data = _xvgfile.parse(file)
				table = _result.addToTable(data, {'group':'Protein', 'tableName':analysisName, 'tableStruct':TableStructure})
		return table
		

