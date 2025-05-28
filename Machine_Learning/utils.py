#!usr/bin/env/python

'''
Script created by Gemma I. Martinez-Redondo based on a function found in https://stackoverflow.com/questions/7156539/how-do-i-transpose-pivot-a-csv-file-with-python-without-loading-the-whole-file to transpose a dataframe
To use this function, import it in another script like:
	from utils import transpose_file
'''


'''
def remove_columns_from_list(input_list,columns_indexes):
        cols_to_remove = [i if i >= 0 else len(input_list) + i for i in columns_indexes] #Negative values may not be valid
        return [input_list[i] for i in range(len(input_list)) if i not in cols_to_remove]
'''

def transpose_file(file_path, output_file_path='transposed.csv', delimiter=','): #, cols_to_remove=None):
	import csv

	transposed_iterator = zip(*csv.reader(open(file_path),delimiter=delimiter))
	with open(output_file_path, 'w') as out:
		for row in transposed_iterator:
			out.write(delimiter.join(row) + '\n')
			'''
			if cols_to_remove is None:
				out.write(delimiter.join(row) + '\n')
			else:
				out.write(delimiter.join(remove_columns_from_list(list(row),cols_to_remove)) + '\n')
			'''
