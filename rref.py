import random
import traceback


class Utils:
	def __init__(self):
		pass

	@staticmethod
	def pretty_printer(matrix):
		from pprint import PrettyPrinter

		# Create a PrettyPrinter with custom settings
		# custom_printer = PrettyPrinter(indent=4, width=30, depth=3, compact=False)
		custom_printer = PrettyPrinter(compact=False)

		# Use the custom printer to pretty-print the nested list
		custom_printer.pprint(matrix)

	def generate_random_matrix(self, r, c, max_values=50, zero_matrix=False):
		"""Generate a random matrix

		Args:
			r (int): Number of rows
			c (int): Number of columns
			max_values (int): Maximum number of values
			zero_matrix (boolen): Whether to generate zero values
		
		Returns:
			matrix (matrix): randomly generated matrix
		"""
		matrix = list()
		for _ in range(r):
			row = list()
			for _ in range(c):
				if zero_matrix:
					row.append(0)
				else:
					row.append(random.randint(0, max_values))
			matrix.append(row)
		return matrix

	def augment(self, a, b, axis=1):
		"""Matrix Augmentation by concatenate two matrices into one matrix on the right.

		Args:
			a (list): Matrix A
			b (list): Matrix B
			axis (int, optional): _description_. Defaults to 1.

		Returns:
			matrix: Augmented Matrix.
		"""

		assert type(a) == list
		assert type(b) == list
		assert len(a) == len(b), 'A and b must be the same length'

		matrix = list()
		for r1, r2 in zip(a,b):
			r1.extend(r2)
			matrix.append(r1)
		return matrix

	def insertion_sort_by_zeros(self, matrix):
		"""Perform insertion sort by preceding zeros

		Args:
			matrix (list): augmented matrix
		Returns:
			matrix (list): Sorted list
		"""
		def count_preceding_zeros(row):
			# Count the number of preceding zeros in a row
			count = 0
			for element in row:
				if element == 0:
					count += 1
				else:
					break
			return count

		for i in range(1, len(matrix)):
			current_row = matrix[i]
			current_zeros = count_preceding_zeros(current_row)
			j = i - 1
			while j >= 0 and count_preceding_zeros(matrix[j]) > current_zeros:
				matrix[j + 1] = matrix[j]
				j -= 1
			matrix[j + 1] = current_row
		return matrix

	def row_echelon_form(self, matrix, float_precision=2):
		"""Perform row echelon transformation

		Args:
			matrix (list): augmented matrix
			float_precision (int, optional): limit float precision. Defaults to 2.

		Returns:
			matrix: row echelon form augmented matrix
		"""
		rows = len(matrix)
		cols = len(matrix[0]) 

		for i in range(min(rows, cols)):
			# Make the diagonal element in the current column 1
			diagonal_element = matrix[i][i]
			if diagonal_element != 0:
				for j in range(cols):
					matrix[i][j] /= diagonal_element
					matrix[i][j] = round(matrix[i][j], float_precision)
					if matrix[i][j] == -0.0: matrix[i][j] = 0.0
			

			# Make all other elements in the current column 0
			for k in range(i + 1, rows):
				factor = matrix[k][i]
				for j in range(cols):
					matrix[k][j] -= factor * matrix[i][j]
		return matrix

	def reduced_row_echelon_form(self, matrix, float_precision=2):
		"""Perform reduced row echelon transformation

		Args:
			matrix (list): row echelon augmented matrix
			float_precision (int, optional): limit float precision. Defaults to 2.

		Returns:
			matrix: reduced row echelon form augmented matrix
		"""
		rows = len(matrix)
		cols = len(matrix[0]) - 1   #TODO: test me

		for i in range(min(rows, cols) - 1, 0, -1):
			# Make the leading entry in the current row 1
			try:
				leading_entry = matrix[i][i]
				if leading_entry != 0:
					for j in range(cols):
						matrix[i][j] /= leading_entry
						matrix[i][j] = round(matrix[i][j], float_precision)
						# if matrix[i][j] == -0.0: matrix[i][j] = 0.0

					# Eliminate other entries in the current column
					for k in range(i - 1, -1, -1):
						factor = matrix[k][i]
						for j in range(cols):
							matrix[k][j] -= factor * matrix[i][j]
			except:
				traceback.print_exc()
				print(i)
				exit()
		return matrix

	def find_pivot_nonpivot_columns(self, matrix):
		"""Find pivot and nonpivot columns in the rref given matrix

		Args:
			matrix (list): rref matrix

		Returns:
			pivot_columns (list): pivot columns indexes
			nonpivot_column (list): nonpivot columns indexes
		"""
		rows = len(matrix)
		cols = len(matrix[0])

		pivot_columns = []
		nonpivot_columns = []

		for j in range(cols):
			# Check if the column contains a leading entry (pivot)
			has_pivot = False
			pivot_row = -1
			for i in range(j ,rows): #TODO: verify me later verified by ronak
				if matrix[i][j] != 0:
					has_pivot = True
					pivot_row = i
					break

			if has_pivot:
				pivot_columns.append(j)
			else:
				nonpivot_columns.append(j)

		return pivot_columns, nonpivot_columns

	def split_rref(self, rref_matrix):
		"""Remove b from rref_matrix

		Args:
			rref_matrix (list): rref matrix

		Returns:
			mat (list): b removed matrix from rref_matrix
		"""
		Arr = list()
		Brr = list()
		for r in rref_matrix:
			Arr.append(r[:-1])
			Brr.append(r[-1])
		return Arr, Brr

	def extract_particular_solution(self, rref_matrix):
		"""Find particular souln form given rref matrix

		Args:
			rref_matrix (list): rref matrix

		Returns:
			particular_solution: particular solution of rref matrix
		"""
		rows = len(rref_matrix)
		cols = len(rref_matrix[0])

		pivot_columns, nonpivot_columns = self.find_pivot_nonpivot_columns(rref_matrix)

		particular_solution = list()
		for col in range(cols):
			if col in pivot_columns:
				particular_solution.append(rref_matrix[col][-1])
			else:
				particular_solution.append(0)


		return particular_solution

	def mul(self, a, b):
		mat = a
		for ridx, (ra, rb) in enumerate(zip(a,b)):
			for cidx, (ca, cb) in enumerate(zip(ra, rb)):
				mat[ridx][cidx] = ca * cb
		return mat

	def get_A0_soln(self, Arr):
		pivot_columns, nonpivot_columns = self.find_pivot_nonpivot_columns(Arr)

		rows = len(Arr)
		cols = len(Arr[0])

		soln_matrix = [[0] * cols for _ in nonpivot_columns]
		for idx, nnp_idx in enumerate(nonpivot_columns):
			for ridx, r in enumerate(Arr):
				soln_matrix[idx][ridx] = r[nnp_idx]
			soln_matrix[idx][nnp_idx] = -1
		return soln_matrix



if __name__ == "__main__":

	m, n = 5, 7
	utils = Utils()
	pprint = utils.pretty_printer

	A = utils.generate_random_matrix(m, n)
	# A = [[1, 0, 8, -4],
	#   	 [0 ,1, 2, 12]
	# 	]
	A = [[48, 19, 32, 46, 3, 32, 29],
		[8, 50, 17, 48, 18, 11, 28],
		[50, 28, 1, 20, 12, 43, 36],
		[29, 16, 28, 6, 49, 30, 4],
		[8, 38, 43, 25, 46, 39, 20]]

	# B = utils.generate_random_matrix(m, 1, zero_matrix=True)
	B = utils.generate_random_matrix(m, 1)
	# B = [[42], [8]]
	B = [[25], [18], [44], [41], [12]]

	print('---- A ---- ')
	pprint(A)
	print('---- B ---- ')
	pprint(B)

	augmented_matrix = utils.augment(A, B)
	print('---- augmented_matrix ---- ')
	pprint(augmented_matrix)

	print('---- insertion_sort_by_zeros ---- ')
	sorted_augmented_matrix = utils.insertion_sort_by_zeros(augmented_matrix)
	pprint(sorted_augmented_matrix)

	print('---- row_echelon_form ---- ')
	ref_matrix = utils.row_echelon_form(sorted_augmented_matrix)
	pprint(ref_matrix)

	
	print('---- reduced_row_echelon_form ---- ')
	rref_matrix = utils.reduced_row_echelon_form(ref_matrix)
	pprint(rref_matrix)

	print('---- find_pivot_nonpivot_columns ---- ')
	pivot_columns, nonpivot_columns = utils.find_pivot_nonpivot_columns(rref_matrix)
	pprint(pivot_columns)
	pprint(nonpivot_columns)

	print('---- extract_particular_solution ---- ')
	Arr, Brr = utils.split_rref(rref_matrix) # TODO: change method name later
	perticular_solution = utils.extract_particular_solution(rref_matrix)
	pprint(perticular_solution)

	print('---- find_special_solutions ---- ')
	# B0 = utils.generate_random_matrix(m, 1, zero_matrix=True)
	special_solutions = utils.get_A0_soln(Arr)
	pprint(special_solutions)

	print('---- general_solutions ---- ')
	ps = list()
	for idx, solution in enumerate(special_solutions):
		ps.append(f'lambda{idx+1}*({solution})')
	print(f"x = {perticular_solution} + {' + '.join(ps)}")

